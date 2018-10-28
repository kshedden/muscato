// Copyright 2017, Kerby Shedden and the Muscato contributors.

/*
Generate simple data sets for testing.  Gene i contains an exact
copy of read i, starting at position i % 10.  If there are more
genes than reads, then the additional genes do not have any
matches beyond those occuring by chance.
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
)

var (
	numRead int
	readLen int
	numGene int
	geneLen int

	reads []string
)

func generateReads() {

	fid, err := os.Create("reads.fastq")
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	w := bufio.NewWriter(fid)
	defer w.Flush()

	var buf bytes.Buffer

	for i := 0; i < numRead; i++ {
		w.WriteString(fmt.Sprintf("read_%d\n", i))

		buf.Reset()
		for j := 0; j < readLen; j++ {
			x := rand.Float64()
			switch {
			case x < 0.25:
				buf.Write([]byte("A"))
			case x < 0.5:
				buf.Write([]byte("T"))
			case x < 0.75:
				buf.Write([]byte("G"))
			default:
				buf.Write([]byte("C"))
			}
		}
		w.Write(buf.Bytes())

		if i < 10 {
			reads = append(reads, string(buf.Bytes()))
		}

		w.WriteString("\n+\n")
		for j := 0; j < readLen; j++ {
			w.WriteString("!")
		}
		w.WriteString("\n")
	}
}

func writeRand(w io.Writer, n int) {

	for j := 0; j < n; j++ {
		x := rand.Float64()
		switch {
		case x < 0.25:
			w.Write([]byte("A"))
		case x < 0.5:
			w.Write([]byte("T"))
		case x < 0.75:
			w.Write([]byte("G"))
		default:
			w.Write([]byte("C"))
		}
	}
}

func generateGenes() {

	fid, err := os.Create("genes.txt")
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	w := bufio.NewWriter(fid)
	defer w.Flush()

	for i := 0; i < numGene; i++ {
		w.WriteString(fmt.Sprintf("gene_%d\t", i))

		if i < len(reads) {
			j := i % 10
			writeRand(w, j)
			w.WriteString(reads[i])
			writeRand(w, geneLen-(readLen+j))
		} else {
			writeRand(w, geneLen)
		}
		w.WriteString("\n")
	}
}

func main() {

	flag.IntVar(&numRead, "NumRead", 10000, "Number of reads")
	flag.IntVar(&readLen, "ReadLen", 100, "Read length")
	flag.IntVar(&numGene, "NumGene", 10000, "Number of genes")
	flag.IntVar(&geneLen, "GeneLen", 1000, "Gene length")

	flag.Parse()

	generateReads()
	generateGenes()
}
