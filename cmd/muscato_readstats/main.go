// Copyright 2017, Kerby Shedden and the Muscato contributors.
//
// readStats calculates statistics for each read, using a results
// datafile that is sorted by read.

package main

import (
	"bufio"
	"bytes"
	"fmt"
	"log"
	"os"
	"path"

	"github.com/kshedden/muscato/utils"
)

var (
	config *utils.Config

	tmpdir string
)

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

	fid, err := os.Open(config.ResultsFileName)
	if err != nil {
		if os.IsNotExist(err) {
			msg := fmt.Sprintf("Cannot open results file %s, see log files for details.\n", config.ResultsFileName)
			os.Stderr.WriteString(msg)
		}
		log.Fatal(err)
	}
	defer fid.Close()

	var outfile string
	ext := path.Ext(config.ResultsFileName)
	if ext != "" {
		m := len(config.ResultsFileName)
		outfile = config.ResultsFileName[0:m-len(ext)] + "_readstats" + ext
	} else {
		outfile = config.ResultsFileName + "_readstats"
	}
	out, err := os.Create(outfile)
	if err != nil {
		msg := fmt.Sprintf("Cannot create %s, see log files for details.\n", outfile)
		os.Stderr.WriteString(msg)
		log.Fatal(err)
	}
	defer out.Close()

	scanner := bufio.NewScanner(fid)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	var oldread, read []byte
	var first bool = true
	var n int
	genes := make(map[string]bool)

	writeout := func(read []byte) error {
		var buf bytes.Buffer
		for g, _ := range genes {
			buf.Write([]byte(g))
			buf.Write([]byte(";"))
		}
		_, err := out.WriteString(fmt.Sprintf("%s\t%s\n", read, buf.String()))
		if err != nil {
			return err
		}
		return nil
	}

	for scanner.Scan() {
		fields := bytes.Fields(scanner.Bytes())
		read = fields[7]

		if first {
			oldread = read
			first = false
		}

		if bytes.Compare(read, oldread) != 0 {
			err := writeout(oldread)
			if err != nil {
				os.Stderr.WriteString("Error in readStats, see log files for details.\n")
				log.Fatal(err)
			}
			oldread = []byte(string(read))
			n = 0
			genes = make(map[string]bool)
		}

		n++
		genes[string(fields[4])] = true
	}

	err = writeout(read)
	if err != nil {
		os.Stderr.WriteString("Error in readStats, see log files for details.\n")
		log.Fatal(err)
	}

	if err := scanner.Err(); err != nil {
		os.Stderr.WriteString("Error in readStats, see log files for details.\n")
		log.Fatal(err)
	}
}
