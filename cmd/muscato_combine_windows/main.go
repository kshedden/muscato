// Copyright 2017, Kerby Shedden and the Muscato contributors.
//
// muscato_combine_windows takes all matches for the same read, then
// retains only those with nmiss equal to at most one greater than
// the lowest nmiss.

package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"path"
	"strconv"
	"strings"

	"github.com/kshedden/muscato/utils"
)

var (
	config *utils.Config

	sortpar     string
	sortmem     string
	sortTmpFlag string
	tmpdir      string

	logger *log.Logger
)

// writebest accepts a set of lines (lines), which have also been
// broken into fields (bfr).  Every line represents a candidate match.
// The matches with at most mmtol more matches than the best match are
// printed out.  ibuf is provided workspace.
func writebest(lines []string, bfr [][]string, ibuf []int, mmtol int) ([]int, error) {

	// Find the best fit, determine the number of mismatches for each sequence.
	ibuf = ibuf[0:0]
	best := -1
	for _, x := range bfr {
		y, err := strconv.Atoi(x[3]) // 3 is position of nmiss
		if err != nil {
			return nil, err
		}
		if best == -1 || y < best {
			best = y
		}
		ibuf = append(ibuf, y)
	}

	// Output the sequences with acceptable number of mismatches.
	for i, x := range lines {
		if ibuf[i] <= best+mmtol {
			fmt.Println(x)
		}
	}

	return ibuf, nil
}

func setupLog() {

	logname := path.Join(config.LogDir, "muscato_combine_windows.log")

	fid, err := os.Create(logname)
	if err != nil {
		panic(err)
	}

	logger = log.New(fid, "", log.Ltime)
}

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

	setupLog()
	logger.Print("starting combineWindows")

	mmtol := config.MMTol

	scanner := bufio.NewScanner(os.Stdin)
	var lines []string
	var fields [][]string
	var ibuf []int
	var current string
	var err error
	for scanner.Scan() {

		line := scanner.Text()
		field := strings.Fields(line)

		// Add to the current block.
		if current == "" || field[0] == current {
			lines = append(lines, line)
			fields = append(fields, field)
			current = field[0]
			continue
		}

		// Process a block
		ibuf, err = writebest(lines, fields, ibuf, mmtol)
		if err != nil {
			msg := "Error in combineWindows, see log file for details.\n"
			os.Stderr.WriteString(msg)
			log.Fatal(err)
		}
		lines = lines[0:0]
		lines = append(lines, line)
		fields = fields[0:0]
		fields = append(fields, field)
		current = field[0]
	}

	if err := scanner.Err(); err == nil {
		// Process the final block if possible
		_, err := writebest(lines, fields, ibuf, mmtol)
		if err != nil {
			msg := "Error in combineWindows, see log file for details.\n"
			os.Stderr.WriteString(msg)
			log.Fatal(err)
		}
	} else {
		// Should never get here, but just in case log
		// the error but don't try to process the
		// remaining lines which may be corrupted.
		logger.Printf("%v", err)
	}

	logger.Print("combineWindows done")
}
