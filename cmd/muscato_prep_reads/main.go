// Copyright 2017, Kerby Shedden and the Muscato contributors.

// muscato_prep_reads converts a source file of sequencing reads from
// fastq format to a simple format with one sequence per row, used
// internally by Muscato.

package main

import (
	"bytes"
	"fmt"
	"log"
	"os"
	"path"

	"github.com/kshedden/muscato/utils"
)

const (
	// The maximum length of a read identifier
	maxNameLen = 1000
)

var (
	config *utils.Config

	tmpdir string

	logger *log.Logger
)

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

func source() {

	ris := utils.NewReadInSeq(config.ReadFileName, "")

	var bbuf bytes.Buffer

	nskip := 0

	var lnum int
	for lnum = 0; ris.Next(); lnum++ {

		bbuf.Reset()

		if len(ris.Seq) < config.MinReadLength {
			nskip++
			continue
		}

		xseq := []byte(ris.Seq)
		subx(xseq)

		if len(xseq) > config.MaxReadLength {
			xseq = xseq[0:config.MaxReadLength]
		}

		_, err := bbuf.Write(append(xseq, '\t'))
		if err != nil {
			panic(err)
		}

		rn := ris.Name
		if len(rn) > maxNameLen {
			rn = rn[0:(maxNameLen-5)] + "..."
		}
		bbuf.Write([]byte(rn))

		bbuf.Write([]byte("\n"))

		_, err = os.Stdout.Write(bbuf.Bytes())
		if err != nil {
			panic(err)
		}
	}

	logger.Printf("Skipped %d reads for being too short", nskip)
}

func setupLog() {
	logname := path.Join(config.LogDir, "muscato_prep_reads.log")
	fid, err := os.Create(logname)
	if err != nil {
		panic(err)
	}
	logger = log.New(fid, "", log.Ltime)
}

func main() {
	if len(os.Args) != 2 && len(os.Args) != 3 {
		os.Stderr.WriteString(fmt.Sprintf("%s: wrong number of arguments\n", os.Args[0]))
		os.Exit(1)
	}

	config = utils.ReadConfig(os.Args[1])

	if config.TempDir == "" {
		tmpdir = os.Args[2]
	} else {
		tmpdir = config.TempDir
	}

	setupLog()
	logger.Printf("Starting prep_reads")
	source()
	logger.Printf("Done")
}
