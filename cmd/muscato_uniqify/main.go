// Copyright 2017, Kerby Shedden and the Muscato contributors.

// muscato_uniqify is a simple stream processor...

package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"log"
	"os"
	"path"
	"strings"

	"github.com/golang/snappy"
	"github.com/kshedden/muscato/utils"
)

var (
	logger *log.Logger

	config *utils.Config
)

func setupLog() {
	logname := path.Join(config.LogDir, "muscato_uniqify.log")
	fid, err := os.Create(logname)
	if err != nil {
		panic(err)
	}
	logger = log.New(fid, "", log.Ltime)
}

func main() {

	if len(os.Args) != 3 {
		msg := fmt.Sprintf("%s: wrong number of arguments", os.Args[0])
		os.Stderr.WriteString(msg)
		os.Exit(1)
	}

	config = utils.ReadConfig(os.Args[1])

	setupLog()

	var fid io.ReadCloser
	if os.Args[2] == "-" {
		fid = os.Stdin
	} else {
		fid, err := os.Open(os.Args[2])
		if err != nil {
			log.Fatal(err)
		}
		defer fid.Close()
	}

	rdr := bufio.NewReader(fid)
	scanner := bufio.NewScanner(rdr)
	buf := make([]byte, 1024*1024)
	scanner.Buffer(buf, 1024*1024)

	wtr := snappy.NewBufferedWriter(os.Stdout)
	defer wtr.Close()

	// Try to read one line to prime the pipeline.
	if !scanner.Scan() {
		// Can't read even one line
		if err := scanner.Err(); err != nil {
			log.Fatal(err)
		}
		log.Fatal(fmt.Errorf("%s: no input from %s", os.Args[0], os.Args[2]))
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

		_, err := wtr.Write(seq)
		if err != nil {
			panic(err)
		}
		_, err = wtr.Write([]byte(fmt.Sprintf("\t%d\t", len(names))))
		if err != nil {
			panic(err)
		}
		_, err = wtr.Write([]byte(na))
		if err != nil {
			panic(err)
		}
		_, err = wtr.Write([]byte("\n"))
		if err != nil {
			panic(err)
		}
	}

	var nseq, nunq int
	for scanner.Scan() {

		line = scanner.Bytes()
		toks := bytes.Split(line, []byte("\t"))
		nseq++

		if bytes.Compare(toks[0], seq) != 0 {
			printrow(seq, names)
			nunq++
			seq = seq[0:0]
			names = names[0:0]
			seq = append(seq, toks[0]...)
		}
		names = append(names, string(toks[1]))
	}

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}

	printrow(seq, names)
	nunq++
	nseq++

	os.Stderr.WriteString(fmt.Sprintf("Found %d total sequences\n", nseq))
	os.Stderr.WriteString(fmt.Sprintf("Found %d unique sequences\n", nunq))

	writeSeqInfo(nseq, nunq)
}

func writeSeqInfo(nseq, nunq int) {

	seqinfo := struct {
		NumUnique int
		NumTotal  int
	}{
		NumUnique: nunq,
		NumTotal:  nseq,
	}

	fid, err := os.Create(path.Join(config.LogDir, "seqinfo.json"))
	if err != nil {
		log.Fatal(err)
	}
	enc := json.NewEncoder(fid)
	enc.Encode(seqinfo)
	fid.Close()
}
