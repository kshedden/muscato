package main

import (
	"flag"
	"fmt"
	"math/rand"
	"os"
)

var (
	numRead int
	readLen int
	numGene int
	geneLen int
)

func generateReads() {

	fid, err := os.Create("reads.fastq")
	if err != nil {
		panic(err)
	}
	defer fid.Close()

	for i := 0; i < numRead; i++ {
		fid.WriteString(fmt.Sprintf("read_%d\n", i))
		for j := 0; j < readLen; j++ {
			x := rand.Float64()
			switch {
			case x < 0.25:
				fid.WriteString("A")
			case x < 0.5:
				fid.WriteString("T")
			case x < 0.75:
				fid.WriteString("G")
			default:
				fid.WriteString("C")
			}
		}
		fid.WriteString("\n+\n")
		for j := 0; j < readLen; j++ {
			fid.WriteString("!")
		}
		fid.WriteString("\n")
	}
}

func generateGenes() {

	fid, err := os.Create("genes.txt")
	if err != nil {
		panic(err)
	}
	defer fid.Close()

	for i := 0; i < numGene; i++ {
		fid.WriteString(fmt.Sprintf("gene_%d\t", i))
		for j := 0; j < geneLen; j++ {
			x := rand.Float64()
			switch {
			case x < 0.25:
				fid.WriteString("A")
			case x < 0.5:
				fid.WriteString("T")
			case x < 0.75:
				fid.WriteString("G")
			default:
				fid.WriteString("C")
			}
		}
		fid.WriteString("\n")
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
