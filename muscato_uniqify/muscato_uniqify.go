// Copyright 2017, Kerby Shedden and the Muscato contributors.

// muscato_uniqify is a simple stream processor...

package main

import (
	"bufio"
	"bytes"
	"fmt"
	"os"
	"strings"

	"github.com/golang/snappy"
)

func main() {

	if len(os.Args) != 2 {
		msg := fmt.Sprintf("%s: wrong number of arguments", os.Args[0])
		os.Stderr.WriteString(msg)
		os.Exit(2)
	}

	fid, err := os.Open(os.Args[1])
	if err != nil {
		panic(err)
	}
	defer fid.Close()

	rdr := bufio.NewReader(fid)
	scanner := bufio.NewScanner(rdr)
	buf := make([]byte, 1024*1024)
	scanner.Buffer(buf, 1024*1024)

	wtr := snappy.NewWriter(os.Stdout)

	if !scanner.Scan() {
		if err := scanner.Err(); err != nil {
			panic(err)
		}
		os.Stderr.WriteString(fmt.Sprintf("%s: no input", os.Args[0]))
		os.Exit(1)
	}

	// Current read sequence
	var seq []byte

	// All names matching the current read sequence
	var names []string

	line := scanner.Bytes()
	toks := bytes.Split(line, []byte("\t"))

	seq = append(seq, toks[0]...)
	names = append(names, string(toks[1]))

	printrow := func(seq []byte, names []string) {
		na := strings.Join(names, ";")
		if len(na) > 1000 {
			na = na[0:996] + "..."
		}
		wtr.Write(seq)
		wtr.Write([]byte(fmt.Sprintf("\t%d\t", len(names))))
		wtr.Write([]byte(na))
		wtr.Write([]byte("\n"))
	}

	for scanner.Scan() {

		line = scanner.Bytes()
		toks := bytes.Split(line, []byte("\t"))

		if bytes.Compare(toks[0], seq) != 0 {
			printrow(seq, names)
			seq = seq[0:0]
			names = names[0:0]
			seq = append(seq, toks[0]...)
		}
		names = append(names, string(toks[1]))
	}

	if err := scanner.Err(); err != nil {
		panic(err)
	}

	printrow(seq, names)
}
