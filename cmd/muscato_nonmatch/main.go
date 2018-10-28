package main

import (
	"bufio"
	"bytes"
	"fmt"
	"log"
	"os"
	"path"
	"strings"

	"github.com/golang/snappy"
	"github.com/kshedden/muscato/utils"
	"github.com/willf/bloom"
)

var (
	config *utils.Config

	tmpdir string

	logger *log.Logger
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

	// Reader for the match file
	inf, err := os.Open(config.ResultsFileName)
	if err != nil {
		if os.IsNotExist(err) {
			msg := fmt.Sprintf("Cannot open file %s.", config.ResultsFileName)
			os.Stderr.WriteString(msg)
			os.Exit(1)
		}
		log.Fatal(err)
	}
	defer inf.Close()

	// Build a bloom filter based on the matched sequences
	billion := uint(1000 * 1000 * 1000)
	bf := bloom.New(4*billion, 5)
	scanner := bufio.NewScanner(inf)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)
	for scanner.Scan() {
		f := bytes.Fields(scanner.Bytes())
		bf.Add(f[0])
	}
	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}

	// Open the nonmatch output file
	a, b := path.Split(config.ResultsFileName)
	c := strings.Split(b, ".")
	d := c[len(c)-1]
	c[len(c)-1] = "nonmatch"
	c = append(c, d+".fastq")
	outname := path.Join(a, strings.Join(c, "."))
	out, err := os.Create(outname)
	if err != nil {
		msg := fmt.Sprintf("Cannot create file %s.", outname)
		if os.IsNotExist(err) {
			os.Stderr.WriteString(msg)
			os.Exit(1)
		}
		log.Fatal(msg)
	}
	defer out.Close()
	wtr := bufio.NewWriter(out)
	defer wtr.Flush()

	// Check each read to see if it was matched.
	rfname := path.Join(config.TempDir, "reads_sorted.txt.sz")
	inf, err = os.Open(rfname)
	if err != nil {
		log.Fatal(err)
	}
	defer inf.Close()
	rdr := snappy.NewReader(inf)
	scanner = bufio.NewScanner(rdr)
	var buf bytes.Buffer
	for scanner.Scan() {
		f := bytes.Fields(scanner.Bytes())
		if !bf.Test(f[0]) {
			buf.Reset()
			buf.Write(f[2])
			buf.WriteString("#")
			buf.Write(f[1])
			buf.WriteString("\n")
			buf.Write(f[0])
			buf.WriteString("\n+\n")
			for k := 0; k < len(f[0]); k++ {
				buf.WriteString("!")
			}
			buf.WriteString("\n")
			_, err = wtr.Write(buf.Bytes())
			if err != nil {
				log.Fatal(err)
			}
		}
	}
}
