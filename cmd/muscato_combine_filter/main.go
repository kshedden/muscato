// Copyright 2017, Kerby Shedden and the Muscato contributors.
//
// muscato_combine_filter reads newline-delimited lines of text
// from multiple Snappy-compressed input files, and prints
// non-duplicated lines to stdio.
//
// Usage:
//
// > muscato_combine_filter n f mode file1 file2...
//
// where n is the approximate number of lines in all files combined,
// f is the desired false positive rate, and mode is either 'check'
// or 'run'.  If mode is 'check', the bit field size and number of
// hashes required to meet the given false positive rate are computed
// and returned.  If mode is 'run', the files are read and filtered.

package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"

	"github.com/golang/snappy"
	"github.com/willf/bloom"
)

// makeReaders creates scanners for reading the source files.  These are
// needed throughout the execution so there is no need to close the underlying
// files.
func makeReaders(files []string) []*bufio.Scanner {

	var scanners []*bufio.Scanner
	for _, f := range files {
		r, err := os.Open(f)
		if err != nil {
			panic(err)
		}

		s := snappy.NewReader(r)

		scanner := bufio.NewScanner(s)
		scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

		scanners = append(scanners, scanner)
	}

	return scanners
}

func main() {

	if len(os.Args) < 5 {
		msg := fmt.Sprintf("Usage: %s num_objects fpr mode file1...\n", os.Args[0])
		os.Stderr.WriteString(msg)
		os.Exit(1)
	}

	mode := os.Args[3]
	if mode != "run" && mode != "check" {
		msg := "The 'mode' argument must be equal to 'run' or 'check'.\n"
		os.Stderr.WriteString(msg)
		os.Exit(1)
	}

	files := os.Args[4:len(os.Args)]
	scanners := makeReaders(files)

	nlines, err := strconv.Atoi(os.Args[1])
	if err != nil {
		panic(err)
	}

	fpr, err := strconv.ParseFloat(os.Args[2], 64)
	if err != nil {
		panic(err)
	}

	m, k := bloom.EstimateParameters(uint(nlines), fpr)
	if mode == "check" {
		fmt.Printf("n=%d\nk=%d\n", m, k)
		os.Exit(0)
	}

	filter := bloom.New(m, k)

	var ix []int
	for j := range scanners {
		ix = append(ix, j)
	}

	var n int
	for len(ix) > 0 {

		for _, i := range ix {

			if scanners[i].Scan() {
				n++
				line := scanners[i].Bytes()
				if !filter.Test(line) {
					fmt.Println(string(line))
					filter.Add(line)
				}
			} else {
				var ixnew []int
				for _, j := range ix {
					if j != i {
						ixnew = append(ixnew, j)
					}
				}
				ix = ixnew
				break // The loop object has changed
			}
		}
	}

	// Check for errors
	for _, scn := range scanners {
		if err := scn.Err(); err != nil {
			panic(err)
		}
	}
}
