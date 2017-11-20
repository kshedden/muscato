// Copyright 2017, Kerby Shedden and the Muscato contributors.

// muscato_screen is an initial screening step used by Muscato to
// identify candidate matches of a set of reads into a set of target
// gene sequences.  The results of the screen may contain false
// positives, but will not contain any false negatives.
//
// The approach is to use a Bloom filter to sketch the reads based on
// the subsequences that appear at defined offsets within the reads.
// For example, if position 10 is an offset and we are looking at
// subequences of width 15, then the read subsequences from position
// 10 through position 25 are entered into a Bloom filter.  Then, we
// scan through every target gene looking for matches to the Bloom
// filter.  When a match occurs, the match position (in the target)
// and flanking sequences are saved for subequent checking against the
// full read sequence.
//
// A simple entropy check is used to avoid considering subsequences
// that could match large numbers of reads or genes (and hence would
// be uninformative).  Currently, this check is based on the number of
// distinct dinucleotide subsequences in the window (e.g. in the
// 15-mer in the example above).
//
// The results are saved in files named bmatch*.txt.sz, where * is the
// window number.
//
// The format of the bmatch files is:
//
// (window sequence) (left tail) (right tail) (gene id) (position)

package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"log"
	"math/rand"
	"os"
	"path"
	"runtime/pprof"
	"strings"
	"sync"

	"github.com/chmduquesne/rollinghash"
	"github.com/chmduquesne/rollinghash/buzhash32"
	"github.com/golang-collections/go-datastructures/bitarray"
	"github.com/golang/snappy"
	"github.com/kshedden/muscato/utils"
)

const (
	// Number of goroutines, around 5-10x the typical number of
	// cores seems to work well.
	concurrency int = 200
)

var (
	// A log
	logger *log.Logger

	// Configuration information
	config *utils.Config

	// All working files are stored here
	tmpdir string

	// Bitarrays that back the Bloom filters
	smp []bitarray.BitArray

	// Tables to produce independent running hashes
	tables [][256]uint32

	// Communicate results back to driver
	hitchan chan rec

	// Semaphore for limiting goroutines
	limit chan bool

	// Line length for output
	bufsize int
)

// genTables generates base hash functions for a collection of rolling hashes.
func genTables() {
	tables = make([][256]uint32, config.NumHash)
	for j := 0; j < config.NumHash; j++ {
		mp := make(map[uint32]bool)
		for i := 0; i < 256; i++ {
			for {
				x := uint32(rand.Int63())
				if !mp[x] {
					tables[j][i] = x
					mp[x] = true
					break
				}
			}
		}
	}
}

// buildBloom constructs bloom filters for each window
func buildBloom() error {

	logger.Printf("Building Bloom sketch of read collection...")

	hashes := make([]rollinghash.Hash32, config.NumHash)
	for j := range hashes {
		hashes[j] = buzhash32.NewFromUint32Array(tables[j])
	}

	fname := path.Join(tmpdir, "reads_sorted.txt.sz")
	fid, err := os.Open(fname)
	if err != nil {
		return err
	}
	defer fid.Close()
	snr := snappy.NewReader(fid)
	scanner := bufio.NewScanner(snr)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	// Workspace for sequence diversity checker
	wk := make([]int, 25)

	var j int
	for ; scanner.Scan(); j++ {

		if j%1000000 == 0 {
			logger.Printf("%d\n", j)
		}

		line := scanner.Bytes()
		seq := bytes.Fields(line)[0]

		for k := 0; k < len(config.Windows); k++ {
			q1 := config.Windows[k]
			q2 := q1 + config.WindowWidth
			if q2 > len(seq) {
				continue
			}
			seqw := seq[q1:q2]

			// Check entropy
			if utils.CountDinuc(seqw, wk) < config.MinDinuc {
				continue
			}

			// Update the Bloom filter for this sequence.
			// This could probably be made concurrent, but
			// there may be too much contention for a big
			// payoff.
			for _, ha := range hashes {
				ha.Reset()
				_, err = ha.Write(seqw)
				if err != nil {
					return err
				}
				x := uint64(ha.Sum32()) % config.BloomSize
				err := smp[k].SetBit(x)
				if err != nil {
					return err
				}
			}
		}
	}

	if err := scanner.Err(); err != nil {
		msg := fmt.Sprintf("Problem reading reads_sorted.txt.sz on line %d\n", j)
		os.Stderr.WriteString(msg)
		return err
	}

	logger.Printf("Done constructing Bloom filters")
	return nil
}

type rec struct {
	mseq  string
	left  string
	right string
	win   int
	tnum  int
	pos   uint32
}

// checkWin returns the indices of the Bloom filters that match the
// current state of the hashes.  iw is workspace and hashses contains
// the hashes that define the Bloom filters.
func checkWin(ix []int, iw []uint64, hashes []rollinghash.Hash32) ([]int, error) {

	// Get the hash states
	for j, ha := range hashes {
		iw[j] = uint64(ha.Sum32()) % config.BloomSize
	}

	ix = ix[0:0]

	// Loop over Bloom filters
	for k, ba := range smp {

		// Determine if the Bloom filter matches
		g := true
		for j := range hashes {
			f, err := ba.GetBit(iw[j])
			if err != nil {
				return nil, err
			}
			if !f {
				g = false
				break
			}
		}
		if g {
			ix = append(ix, k)
		}
	}

	return ix, nil
}

// process one target sequence, runs concurrently with main loop.
func processseq(seq []byte, genenum int, errc chan error) {

	defer func() { <-limit }()

	hashes := make([]rollinghash.Hash32, config.NumHash)
	for j := range hashes {
		hashes[j] = buzhash32.NewFromUint32Array(tables[j])
	}

	// Initialize the hashes with the first window.
	hlen := config.WindowWidth
	if len(seq) < hlen {
		// Not long enough to fit even one window.
		return
	}
	for j := range hashes {
		_, err := hashes[j].Write(seq[0:hlen])
		if err != nil {
			errc <- err
			return
		}
	}

	// Will contain the indices of the matching windows
	ix := make([]int, len(smp))

	// Workspace
	iw := make([]uint64, config.NumHash)

	// Check if the initial window is a match
	var err error
	ix, err = checkWin(ix, iw, hashes)
	if err != nil {
		errc <- err
		return
	}

	for _, i := range ix {

		q1 := config.Windows[i]
		if q1 != 0 {
			// The only way the full read can match at the
			// beginning of the target is if the first
			// window starts at the beginning of the read.
			continue
		}
		q2 := q1 + config.WindowWidth

		jz := 100 - q2
		if jz > len(seq) {
			jz = len(seq)
		}
		hitchan <- rec{
			mseq:  string(seq[0:hlen]),
			left:  "",
			right: string(seq[hlen:jz]),
			tnum:  genenum,
			win:   i,
			pos:   0,
		}
	}

	// Check the rest of the windows
	for j := hlen; j < len(seq); j++ {

		for _, ha := range hashes {
			ha.Roll(seq[j])
		}
		ix, err = checkWin(ix, iw, hashes)
		if err != nil {
			errc <- err
			return
		}

		// Process a match
		for _, i := range ix {

			q1 := config.Windows[i]
			q2 := q1 + config.WindowWidth
			if j < q2-1 {
				// The read would not fit
				continue
			}

			// Matching sequence is jx:jy
			jx := j - hlen + 1
			jy := j + 1

			// Left tail is jw:jx
			jw := jx - q1

			// Right tail is jy:jz
			jz := jy + config.MaxReadLength - q2
			if jz > len(seq) {
				// May not be long enough to fit, but
				// we don't know until we merge.
				jz = len(seq)
			}

			if jw >= 0 {
				hitchan <- rec{
					mseq:  string(seq[jx:jy]),
					left:  string(seq[jw:jx]),
					right: string(seq[jy:jz]),
					tnum:  genenum,
					win:   i,
					pos:   uint32(j - hlen + 1),
				}
			}
		}
	}
}

// Retrieve the results and write to disk
func harvest(wg *sync.WaitGroup) {

	var warn int

	var wtrs []io.Writer
	var allwtrs []io.Closer
	for k := 0; k < len(config.Windows); k++ {
		f := fmt.Sprintf("bmatch_%d.txt.sz", k)
		outname := path.Join(tmpdir, f)
		out, err := os.Create(outname)
		if err != nil {
			logger.Print(err)
			panic(err)
		}
		wtr := snappy.NewBufferedWriter(out)
		wtrs = append(wtrs, wtr)
		allwtrs = append(allwtrs, wtr, out)
	}

	bb := bytes.Repeat([]byte(" "), bufsize)
	bb[bufsize-1] = byte('\n')

	for r := range hitchan {

		if len(hitchan) > cap(hitchan)/2 {
			if warn%100 == 0 {
				logger.Print("hitchan more than half full")
			}
			warn++
		}

		wtr := wtrs[r.win]

		n1, err1 := wtr.Write([]byte(fmt.Sprintf("%s\t", r.mseq)))
		n2, err2 := wtr.Write([]byte(fmt.Sprintf("%s\t", r.left)))
		n3, err3 := wtr.Write([]byte(fmt.Sprintf("%s\t", r.right)))
		n4, err4 := wtr.Write([]byte(fmt.Sprintf("%011d\t", r.tnum)))
		n5, err5 := wtr.Write([]byte(fmt.Sprintf("%d", r.pos)))

		for _, err := range []error{err1, err2, err3, err4, err5} {
			if err != nil {
				logger.Print(err)
				panic("writing error")
			}
		}

		n := n1 + n2 + n3 + n4 + n5
		if n > bufsize {
			panic("output line is too long")
		}

		// The rest of the line is spaces, then newline.
		_, err := wtr.Write(bb[n:bufsize])
		if err != nil {
			logger.Print(err)
			panic(err)
		}
	}

	for _, wtr := range allwtrs {
		wtr.Close()
	}
	wg.Done()
	logger.Printf("Exiting harvest")
}

// search loops through the target sequences, checking each window
// within each target gene for possible matches to the read
// collection.
func search() error {

	logger.Printf("Checking target sequences for matches...")

	fid, err := os.Open(config.GeneFileName)
	if err != nil {
		return err
	}
	defer fid.Close()
	snr := snappy.NewReader(fid)

	// Target file contains some very long lines
	scanner := bufio.NewScanner(snr)
	sbuf := make([]byte, 1024*1024)
	scanner.Buffer(sbuf, 1024*1024)

	hitchan = make(chan rec, 10000)
	limit = make(chan bool, concurrency)
	errc := make(chan error, concurrency)

	var wg sync.WaitGroup
	wg.Add(1)
	go harvest(&wg)

	var i int
	for ; scanner.Scan(); i++ {

		if i%1000000 == 0 {
			logger.Printf("%dM\n", i/1000000)
		}

		line := scanner.Text() // need a copy here

		toks := strings.Split(line, "\t")
		seq := toks[0] // The sequence

		limit <- true
		go processseq([]byte(seq), i, errc)
	}

	if err := scanner.Err(); err != nil {
		msg := fmt.Sprintf("Problem reading %s on line %d\n", config.GeneFileName, i)
		os.Stderr.WriteString(msg)
		logger.Print(err)
		return err
	}

	for k := 0; k < concurrency; k++ {
		limit <- true
	}

	// Get an error if one was generated
	select {
	case e := <-errc:
		log.Fatal(e)
	default:
	}

	close(hitchan)
	wg.Wait()
	logger.Printf("Done checking target sequences for matches")
	return nil
}

func setupLogger() error {
	logname := path.Join(config.LogDir, "muscato_screen.log")
	logfid, err := os.Create(logname)
	if err != nil {
		return err
	}
	logger = log.New(logfid, "", log.Ltime)
	return nil
}

func estimateFullness() error {

	n := 1000
	logger.Printf("Bloom filter fill rates:\n")

	for j, ba := range smp {
		c := 0
		for k := 0; k < n; k++ {
			i := uint64(rand.Int63()) % config.BloomSize
			f, err := ba.GetBit(i)
			if err != nil {
				return err
			}
			if f {
				c++
			}
		}
		logger.Printf("%3d %.3f\n", j, float64(c)/float64(n))
	}

	return nil
}

func main() {

	if len(os.Args) != 2 && len(os.Args) != 3 {
		os.Stderr.WriteString(fmt.Sprintf("%s: wrong number of arguments", os.Args[0]))
		os.Exit(1)
	}

	config = utils.ReadConfig(os.Args[1])

	if config.TempDir == "" {
		tmpdir = os.Args[2]
	} else {
		tmpdir = config.TempDir
	}

	if config.CPUProfile {
		f, err := os.Create(path.Join(config.LogDir, "muscato_screen_cpu.prof"))
		if err != nil {
			panic(err)
		}
		defer f.Close()
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	bufsize = config.MaxReadLength + 50

	err := setupLogger()
	if err != nil {
		log.Fatal(err)
	}

	genTables()

	smp = make([]bitarray.BitArray, len(config.Windows))
	for k := range smp {
		smp[k] = bitarray.NewBitArray(config.BloomSize)
	}

	err = buildBloom()
	if err != nil {
		log.Fatal(err)
	}

	err = estimateFullness()
	if err != nil {
		log.Fatal(err)
	}

	err = search()
	if err != nil {
		log.Fatal(err)
	}
}
