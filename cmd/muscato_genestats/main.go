package main

import (
	"bufio"
	"bytes"
	"fmt"
	"os"
)

func main() {

	fid, err := os.Open(os.Args[1])
	if err != nil {
		panic(err)
	}

	scanner := bufio.NewScanner(fid)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	var oldgene, gene []byte
	var first bool = true
	var n int

	writeout := func(gene []byte) {
		fmt.Printf("%s\t%d\t\n", gene, n)
	}

	for scanner.Scan() {
		fields := bytes.Fields(scanner.Bytes())
		gene = fields[4]

		if first {
			oldgene = gene
			first = false
		}

		if bytes.Compare(gene, oldgene) != 0 {
			writeout(oldgene)
			oldgene = []byte(string(gene))
			n = 0
		}

		n++
	}

	writeout(gene)

	if err := scanner.Err(); err != nil {
		panic(err)
	}
}
