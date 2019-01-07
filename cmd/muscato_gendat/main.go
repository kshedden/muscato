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

	"github.com/golang/snappy"
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
	seq := make([]byte, readLen+geneLen)

	for i := 0; i < numRead; i++ {

		buf.Reset()

		io.WriteString(buf, fmt.Sprintf("read_%d\n", i))

		seq = genRand(readLen, seq)
		buf.Write(seq)

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
			reads = append(reads, string(seq))
		}
	}
}

func genRand(n int, seq []byte) []byte {

	bases := []byte{'A', 'T', 'G', 'C'}

	if cap(seq) < n {
		seq = make([]byte, n)
	}
	seq = seq[0:n]

	for j := 0; j < n; j++ {
		x := rand.Float64()
		k := int(4 * x)
		seq[j] = bases[k]
	}

	return seq
}

func generateGenes() {

	seq := make([]byte, geneLen+readLen)

	fname := path.Join(dir, "genes.txt.sz")
	fid, err := os.Create(fname)
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	w := snappy.NewBufferedWriter(fid)
	defer w.Close()

	fmt.Printf("Writing %d genes\n", numGene)
	for i := 0; i < numGene; i++ {

		_, err := io.WriteString(w, fmt.Sprintf("gene_%d\t", i))
		if err != nil {
			panic(err)
		}

		seq = genRand(geneLen, seq)

		if i < numGene/2 {
			j := i % 10
			copy(seq[j:len(seq)], reads[j])
		}

		if _, err := w.Write(seq); err != nil {
			panic(err)
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
