// Copyright 2017, Kerby Shedden and the Muscato contributors.

/*
Generate simple data sets for testing.

In the first half of the genes, gene i contains an exact copy of read
i % 10, starting at position i % 10. The remainder of these gene sequences
are random.

The second half of the gene sequences are random and should contain
few or no matches.
*/

package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io"
	"math/rand"
	"os"
	"path"
)

var (
	numRead int
	readLen int
	numGene int
	geneLen int
	dir     string

	reads []string
)

func generateReads() {

	fmt.Printf("Writing %d reads\n", numRead)

	fname := path.Join(dir, "reads.fastq")
	fid, err := os.Create(fname)
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	w := bufio.NewWriter(fid)
	defer w.Flush()

	buf := new(bytes.Buffer)
	seq := new(bytes.Buffer)

	for i := 0; i < numRead; i++ {

		buf.Reset()
		seq.Reset()

		io.WriteString(buf, fmt.Sprintf("read_%d\n", i))

		for j := 0; j < readLen; j++ {
			x := rand.Float64()
			switch {
			case x < 0.25:
				io.WriteString(seq, "A")
			case x < 0.5:
				io.WriteString(seq, "T")
			case x < 0.75:
				io.WriteString(seq, "G")
			default:
				io.WriteString(seq, "C")
			}
		}
		buf.Write(seq.Bytes())

		io.WriteString(buf, "\n+\n")
		for j := 0; j < readLen; j++ {
			io.WriteString(buf, "!")
		}
		io.WriteString(buf, "\n")

		_, err := w.Write(buf.Bytes())
		if err != nil {
			panic(err)
		}

		if i < 10 {
			reads = append(reads, string(seq.Bytes()))
		}
	}
}

func writeRand(w io.Writer, n int) {

	seq := make([]byte, n)

	for j := 0; j < n; j++ {
		x := rand.Float64()
		switch {
		case x < 0.25:
			seq[j] = 'A'
		case x < 0.5:
			seq[j] = 'T'
		case x < 0.75:
			seq[j] = 'G'
		default:
			seq[j] = 'C'
		}
	}

	_, err := w.Write(seq)
	if err != nil {
		panic(err)
	}
}

func generateGenes() {

	fname := path.Join(dir, "genes.txt")
	fid, err := os.Create(fname)
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	w := bufio.NewWriter(fid)
	defer w.Flush()

	fmt.Printf("Writing %d genes\n", numGene)
	for i := 0; i < numGene; i++ {

		_, err := io.WriteString(w, fmt.Sprintf("gene_%d\t", i))
		if err != nil {
			panic(err)
		}

		if i < numGene/2 {

			j := i % 10
			writeRand(w, j)

			_, err := io.WriteString(w, reads[j])
			if err != nil {
				panic(err)
			}

			writeRand(w, geneLen-(readLen+j))
		} else {
			writeRand(w, geneLen)
		}

		_, err = io.WriteString(w, "\n")
		if err != nil {
			panic(err)
		}
	}
}

func main() {

	flag.IntVar(&numRead, "NumRead", 10000, "Number of reads")
	flag.IntVar(&readLen, "ReadLen", 100, "Read length")
	flag.IntVar(&numGene, "NumGene", 10000, "Number of genes")
	flag.IntVar(&geneLen, "GeneLen", 1000, "Gene length")
	flag.StringVar(&dir, "Dir", ".", "Directory")

	flag.Parse()

	if numRead < 10 {
		panic("numRead must be at least 10")
	}

	generateReads()
	generateGenes()
}
