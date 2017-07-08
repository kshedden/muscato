// Copyright 2017, Kerby Shedden and the Muscato contributors.

// muscato_confirm takes the pre-screening results and determines
// which read x gene pairs match sufficiently well over the entire
// read.
//
// For each k-mer window used in the pre-screening step, and for each
// high-entropy k-mer sequence that appears in the reads at that
// position, all reads containing the subsequence (at a fixed offset)
// are matched against all target genes containing the subsequence (at
// any position where extension to the full read is possible).  This
// "all pairs" matching is done for each k-mer sequence, and the
// results that match sufficiently well, as determined by the PMatch
// parameter, are retained for further processing.

package main

import (
	"bufio"
	"bytes"
	"fmt"
	"log"
	"os"
	"path"
	"strconv"
	"strings"

	"github.com/golang/snappy"
	"github.com/kshedden/seqmatch/utils"
	"github.com/pkg/profile"
)

const (
	// Number of simultaneous goroutines, maybe should scale with
	// number of available cores.
	concurrency = 100

	doProfile = false

	// Maintain a pool of byte arrays of length bufsize
	poolsize = 10000
)

var (
	logger *log.Logger

	config *utils.Config

	tmpdir string

	// Pool of reusable byte slices
	pool chan []byte

	win int // The window to process, win=0,1,...

	// Pass results to driver then write to disk
	rsltChan chan []byte

	bufsize int

	alldone chan bool
)

type rec struct {
	buf    []byte
	fields [][]byte
}

func (r *rec) Print() {
	fmt.Printf("len(buf)=%d\n", len(r.buf))
	for k, f := range r.fields {
		fmt.Printf("%d %s\n", k, string(f))
	}
}

func (r *rec) release() {
	if r.buf == nil {
		logger.Print("nothing to release")
		panic("nothing to release")
	}
	putbuf(r.buf)
	r.buf = nil
	r.fields = nil
}

func (r *rec) init() {
	if r.buf != nil {
		logger.Print("cannot init non-nil rec")
		panic("cannot init non-nil rec")
	}
	r.buf = getbuf()
	r.buf = r.buf[0:0]
}

func (r *rec) setfields() {
	r.fields = bytes.Split(r.buf, []byte("\t"))
}

// breader iterates through a set of sequences, combining blocks of
// contiguous records with the same window sequence.  A breader can be
// used to iterate through either the match or the raw read data.  The
// input sequence windows must be sorted.
type breader struct {

	// The input sequences
	scanner *bufio.Scanner

	// The caller can access the block data through this field
	recs []*rec

	// If we read past the end of a block, put it here so it can
	// be included in the next iteration.
	stash *rec

	// True if all sequences have been read.  At this point, the
	// recs field will continue to hold the final block of
	// sequences.
	done bool

	// The current line number in the input file
	lnum int

	// The name of the source of sequences (either "match" or
	// "source").
	name string

	// Used to confirm that file is sorted
	last *rec
}

// Next advances a breader to the next block.
func (b *breader) Next() bool {

	if b.done {
		return false
	}

	for _, b := range b.recs {
		b.release()
	}
	b.recs = b.recs[0:0]

	if b.stash != nil {
		b.recs = append(b.recs, b.stash)
		b.stash = nil
	}

	for ii := 0; b.scanner.Scan(); ii++ {

		// Process a line
		bb := b.scanner.Bytes()
		if len(bb) > bufsize {
			logger.Print("line too long")
			panic("line too long")
		}
		rx := new(rec)
		rx.init()
		rx.buf = rx.buf[0:len(bb)]
		copy(rx.buf, bb)
		rx.setfields()

		b.lnum++
		if b.lnum%100000 == 0 {
			logger.Printf("%s: %d\n", b.name, b.lnum)
		}

		if (len(b.recs) > 0) && !bytes.Equal(b.recs[0].fields[0], rx.fields[0]) {
			b.stash = rx
			return true
		}
		// Check sorting (harder to check in other branch of the if).
		if ii > 0 {
			if bytes.Compare(b.last.fields[0], rx.fields[0]) > 0 {
				logger.Print("file is not sorted")
				panic("file is not sorted")
			}
		}
		b.last = rx
		b.recs = append(b.recs, rx)
	}

	if err := b.scanner.Err(); err != nil {
		logger.Print(err)
		panic(err)
	}

	b.done = true
	logger.Printf("%s done", b.name)
	return true
}

// cdiff returns the number of unequal values in two byte sequences
func cdiff(x, y []byte) int {
	var c int
	for i, v := range x {
		if v != y[i] {
			c++
		}
	}
	return c
}

func putbuf(buf []byte) {
	select {
	case pool <- buf[0:0]:
	default:
		// pool is full, buffer goes to garbage
	}
}

func getbuf() []byte {
	var buf []byte
	select {
	case buf = <-pool:
		buf = buf[0:0]
	default:
		buf = make([]byte, 0, bufsize)
	}
	return buf
}

type qrect struct {
	mismatch int
	gob      []byte
}

func searchpairs(source, match []*rec, limit chan bool) {

	defer func() { <-limit }()
	if len(match)*len(source) > 100000 {
		logger.Printf("searching %d %d ...", len(match), len(source))
	}

	var qvals []*qrect

	first := config.MatchMode == "first"

	var stag []byte
	for _, mrec := range match {

		mtag := mrec.fields[0]
		mlft := mrec.fields[1]
		mrgt := mrec.fields[2]
		mgene := mrec.fields[3]
		mpos := mrec.fields[4]

		for _, srec := range source {

			stag = srec.fields[0] // must equal mtag
			slft := srec.fields[1]
			srgt := srec.fields[2]

			// Allowed number of mismatches
			nmiss := int((1 - config.PMatch) * float64(len(stag)+len(slft)+len(srgt)))

			// Gene ends before read would end, can't match.
			if len(srgt) > len(mrgt) {
				continue
			}

			// Count differences
			mk := len(srgt)
			nx := cdiff(mlft, slft)
			nx += cdiff(mrgt[0:mk], srgt)
			if nx > nmiss {
				continue
			}

			// unavoidable []byte to string copy
			mposi, err := strconv.Atoi(strings.TrimRight(string(mpos), " "))
			if err != nil {
				logger.Print(err)
				panic(err)
			}

			// Found a match, pass to output
			buf := getbuf()
			bbuf := bytes.NewBuffer(buf)
			bbuf.Write(slft)
			bbuf.Write(stag)
			bbuf.Write(srgt)
			bbuf.Write([]byte("\t"))
			bbuf.Write(mlft)
			bbuf.Write(mtag)
			bbuf.Write(mrgt[0:mk])
			x := fmt.Sprintf("\t%d\t%d\t%s\n", mposi-len(mlft), nx, mgene)
			bbuf.Write([]byte(x))

			qq := &qrect{mismatch: nx, gob: bbuf.Bytes()}
			if first {
				// Make no attempt to rank matches, just keep first ones.
				qvals = append(qvals, qq)
				if len(qvals) > config.MaxMatches {
					goto E
				}
			} else {
				// A priority queue of top matches.
				qvals = qinsert(qvals, qq)
			}
		}
	}

E:
	for _, v := range qvals {
		rsltChan <- v.gob
	}

	for _, x := range source {
		x.release()
	}
	for _, x := range match {
		x.release()
	}
}

func setupLog(win int) {
	logname := path.Join(tmpdir, fmt.Sprintf("mergebloom_%d.log", win))
	fid, err := os.Create(logname)
	if err != nil {
		panic(err)
	}
	logger = log.New(fid, "", log.Ltime)
}

// rcpy deeply copies its argument.
func rcpy(r []*rec) []*rec {
	x := make([]*rec, len(r))
	for j := range x {
		x[j] = new(rec)
		x[j].init()
		x[j].buf = x[j].buf[0:len(r[j].buf)]
		copy(x[j].buf, r[j].buf)
		x[j].setfields()
	}
	return x
}

func main() {

	if len(os.Args) != 4 {
		panic("wrong number of arguments")
	}

	config = utils.ReadConfig(os.Args[1])

	if config.TempDir == "" {
		tmpdir = os.Args[3]
	} else {
		tmpdir = config.TempDir
	}

	bufsize = 2*config.MaxReadLength + 50

	var err error
	win, err = strconv.Atoi(os.Args[2])
	if err != nil {
		panic(err)
	}
	setupLog(win)

	if doProfile && win == 0 {
		p := profile.Start(profile.ProfilePath("."))
		defer p.Stop()
	}

	f := fmt.Sprintf("win_%d_sorted.txt.sz", win)
	sourcefile := path.Join(tmpdir, f)
	logger.Printf("sourcefile: %s", sourcefile)

	f = fmt.Sprintf("smatch_%d.txt.sz", win)
	matchfile := path.Join(tmpdir, f)
	logger.Printf("matchfile: %s", matchfile)

	f = fmt.Sprintf("rmatch_%d.txt.sz", win)
	outfile := path.Join(tmpdir, f)
	logger.Printf("outfile: %s", outfile)

	pool = make(chan []byte, poolsize)

	// Read source sequences
	fid, err := os.Open(sourcefile)
	if err != nil {
		logger.Print(err)
		panic(err)
	}
	defer fid.Close()
	szr := snappy.NewReader(fid)
	scanner := bufio.NewScanner(szr)
	source := &breader{scanner: scanner, name: "source"}

	// Read candidate match sequences
	gid, err := os.Open(matchfile)
	if err != nil {
		logger.Print(err)
		panic(err)
	}
	defer gid.Close()
	szq := snappy.NewReader(gid)
	scanner = bufio.NewScanner(szq)
	match := &breader{scanner: scanner, name: "match"}

	// Place to write results
	fi, err := os.Create(outfile)
	if err != nil {
		logger.Print(err)
		panic(err)
	}
	defer fi.Close()
	out := snappy.NewBufferedWriter(fi)
	defer out.Close()

	source.Next()
	match.Next()

	rsltChan = make(chan []byte, 5*concurrency)
	limit := make(chan bool, concurrency)
	alldone = make(chan bool)

	// Harvest the results
	go func() {
		for r := range rsltChan {
			_, err := out.Write(r)
			if err != nil {
				panic(err)
			}
			putbuf(r)
		}
		alldone <- true
	}()

lp:
	for ii := 0; ; ii++ {

		if ii%100000 == 0 {
			logger.Printf("%d", ii)
		}

		s := source.recs[0].fields[0]
		m := match.recs[0].fields[0]
		c := bytes.Compare(s, m)

		ms := true
		mb := true

		switch {
		case c == 0:
			// Window sequences match, check if it is a real match.
			limit <- true
			go searchpairs(rcpy(source.recs), rcpy(match.recs), limit)
			ms = source.Next()
			mb = match.Next()
			if !(ms || mb) {
				break lp
			}
		case c < 0:
			// The source sequence is behind, move it up.
			ms = source.Next()
			if !ms {
				break lp
			}
		case c > 0:
			// The match sequence is behind, move it up.
			mb = match.Next()
			if !mb {
				break lp
			}
		}
		if !(ms && mb) {
			// One of the files is done
			logger.Printf("ms=%v, mb=%v\n", ms, mb)
		}
	}

	logger.Print("clearing channel")
	for k := 0; k < cap(limit); k++ {
		limit <- true
	}

	close(rsltChan)
	<-alldone

	logger.Print("done")
}

// qinsert inserts a into the array q, maintaining a heap structure on
// q in which the heap is ordered by decreasing mismatch values.  The
// length of the heap is limited to MaxMatches.
func qinsert(q []*qrect, a *qrect) []*qrect {

	q = append(q, a)
	ii := len(q) - 1 // Position of just-inserted node

	for ii > 0 {
		// Position of parent
		jj := (ii - 1) / 2

		if q[jj].mismatch > q[ii].mismatch {
			q[jj], q[ii] = q[ii], q[jj]
			ii = jj
		} else {
			break
		}
	}

	// Not guaranteed to retain the best matches, but approximate
	// and fast.
	if len(q) > config.MaxMatches {
		q = q[0:config.MaxMatches]
	}

	return q
}
