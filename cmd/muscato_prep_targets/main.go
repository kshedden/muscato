// Copyright 2017, Kerby Shedden and the Muscato contributors.

// muscato_prep_targets converts a gene sequence file to a simple text
// format used internally by Muscato.  The ids and sequences are
// placed into newline-delimited text files, with one id or sequence
// per row.
//
// The input can be either a fasta file, or a text format with each
// line containing an id followed by a tab followed by a sequence.
// Letters other than A/T/G/C are replaced with X.

package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"strings"

	"github.com/golang/snappy"
)

const (
	// Maximum sequence length.  If there are sequences longer
	// than this, the program will exit with an error.
	maxline int = 1024 * 1024
)

var (
	// If true, data are fasta format, else they follow a format
	// with one line per sequence, having format id<tab>sequence.
	fasta bool

	logger *log.Logger
)

// revcomp reverse complements its argument.
func revcomp(seq []byte) []byte {
	m := len(seq) - 1
	b := make([]byte, len(seq))
	for i, x := range seq {
		switch x {
		case 'A':
			b[m-i] = 'T'
		case 'T':
			b[m-i] = 'A'
		case 'G':
			b[m-i] = 'C'
		case 'C':
			b[m-i] = 'G'
		case 'X':
			b[m-i] = 'X'
		}
	}
	return b
}

// subx replaces non A/T/G/C with X
func subx(seq []byte) {
	for i, c := range seq {
		switch c {
		case 'A':
		case 'T':
		case 'C':
		case 'G':
		default:
			seq[i] = 'X'
		}
	}
}

func processText(scanner *bufio.Scanner, idout, seqout io.Writer, rev bool) {

	logger.Print("Processing text format file...")

	var lnum int
	for scanner.Scan() {

		if lnum%1000000 == 0 {
			logger.Printf("%d\n", lnum)
		}

		line := scanner.Bytes()
		if len(line) == 0 {
			break
		}

		toks := bytes.Split(line, []byte("\t"))
		if len(toks) != 2 {
			logger.Printf("Text format gene file should have two tab-delimited tokens per row.  Line %d has %d tokens.\n",
				lnum+1, len(toks))
			os.Exit(0)
		}

		nam := toks[0]
		seq := toks[1]

		subx(seq)

		// Write the sequence
		_, err := seqout.Write(append(seq, '\n'))
		if err != nil {
			panic(err)
		}
		if rev {
			_, err := seqout.Write(append(revcomp(seq), '\n'))
			if err != nil {
				panic(err)
			}
		}

		// Write the gene id
		_, err = idout.Write([]byte(fmt.Sprintf("%011d\t%s\t%d\n", lnum, nam, len(seq))))
		if err != nil {
			panic(err)
		}
		lnum++
		if rev {
			_, err = idout.Write([]byte(fmt.Sprintf("%011d\t%s_r\t%d\n", lnum, nam, len(seq))))
			if err != nil {
				panic(err)
			}
			lnum++
		}
	}

	if err := scanner.Err(); err != nil {
		logger.Printf("Failed on line %d", lnum)
		panic(err)
	}
}

func processFasta(scanner *bufio.Scanner, idout, seqout io.Writer, rev bool) {

	logger.Print("Processing FASTA format file...")

	var seqname string
	var seq []byte
	var lnum int

	flush := func(r bool) {

		// Write the sequence
		_, err := seqout.Write(append(seq, '\n'))
		if err != nil {
			panic(err)
		}

		// Write the gene id
		x := ""
		if r {
			x = "_r"
		}

		_, err = idout.Write([]byte(fmt.Sprintf("%011d\t%s%s\t%d\n", lnum, seqname, x, len(seq))))
		if err != nil {
			panic(err)
		}
	}

	for scanner.Scan() {

		if lnum%1000000 == 0 {
			logger.Printf("%d\n", lnum)
		}

		line := scanner.Bytes()

		if line[0] == '>' {
			if len(seq) > 0 {
				subx(seq)
				flush(false)
				lnum++
				if rev {
					seq = revcomp(seq)
					flush(true)
					lnum++
				}
			}
			seqname = string(line)
			seq = seq[0:0]
			continue
		}

		seq = append(seq, line...)
	}

	if err := scanner.Err(); err != nil {
		logger.Printf("Failed on line %d", lnum)
		logger.Printf("Final sequence name: %s", seqname)
		panic(err)
	}

	if len(seq) > 0 {
		flush(false)
		lnum++
		if rev {
			seq = revcomp(seq)
			flush(true)
			lnum++
		}
	}
}

func targets(genefile string, rev bool) {

	var rdr io.ReadCloser

	// Setup for reading the input file
	var err error
	rdr, err = os.Open(genefile)
	if err != nil {
		panic(err)
	}
	defer func(r io.Closer) { r.Close() }(rdr)

	// The input file is gzipped
	ext := filepath.Ext(genefile)
	if strings.ToLower(ext) == ".gz" {
		logger.Printf("Reading gzipped gene sequence file")
		rdr, err = gzip.NewReader(rdr)
		if err != nil {
			panic(err)
		}
		defer func(r io.Closer) { r.Close() }(rdr)
		genefile = strings.Replace(genefile, ext, "", -1)
		ext = filepath.Ext(genefile)
	}

	// Setup for writing the sequence output
	geneoutfile := strings.Replace(genefile, ext, ".txt.sz", 1)
	gid, err := os.Create(geneoutfile)
	if err != nil {
		panic(err)
	}
	defer gid.Close()
	seqout := snappy.NewBufferedWriter(gid)
	defer seqout.Close()

	// Setup for writing the identifier output
	geneidfile := strings.Replace(genefile, ext, "_ids.txt.sz", 1)
	idwtr, err := os.Create(geneidfile)
	if err != nil {
		panic(err)
	}
	defer idwtr.Close()
	idout := snappy.NewBufferedWriter(idwtr)
	defer idout.Close()

	// Setup a scanner to read long lines
	scanner := bufio.NewScanner(rdr)
	sbuf := make([]byte, 64*1024)
	scanner.Buffer(sbuf, maxline)

	if fasta {
		processFasta(scanner, idout, seqout, rev)
	} else {
		processText(scanner, idout, seqout, rev)
	}

	logger.Printf("Done processing targets")
}

func setupLog() {
	fid, err := os.Create("muscato_prep_targets.log")
	if err != nil {
		panic(err)
	}
	logger = log.New(fid, "", log.Ltime)
}

func main() {

	rev := flag.Bool("rev", false, "Include reverse complement sequences")
	flag.Parse()
	args := flag.Args()

	if len(args) != 1 {
		os.Stderr.WriteString("muscato_prep_targets: usage\n")
		os.Stderr.WriteString("  muscato_prep_targets [-rev] genefile\n\n")
		os.Exit(1)
	}

	genefile := args[0]

	gl := strings.ToLower(genefile)
	if strings.HasSuffix(gl, "fasta") {
		fasta = true
	}

	setupLog()
	if *rev {
		logger.Printf("Including reverse complements")
	} else {
		logger.Printf("Not including reverse complements")
	}

	targets(genefile, *rev)
	logger.Printf("Done")
}
