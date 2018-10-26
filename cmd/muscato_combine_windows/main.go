// Copyright 2017, Kerby Shedden and the Muscato contributors.

// muscato_combine_windows...

package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"path"
	"strconv"
	"strings"

	"github.com/golang/snappy"
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
// written to the io writer (wtr).  ibuf is provided workspace.
func writebest(lines []string, bfr [][]string, wtr io.Writer, ibuf []int, mmtol int) ([]int, error) {

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
			_, err := wtr.Write([]byte(x))
			if err != nil {
				return nil, err
			}
			_, err = wtr.Write([]byte("\n"))
			if err != nil {
				return nil, err
			}
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

	sortpar = fmt.Sprintf("--parallel=%d", config.SortPar)
	sortmem = fmt.Sprintf("-S %s", config.SortMem)

	sortTmpFlag = fmt.Sprintf("--temporary-directory=%s", config.SortTemp)

	mmtol := config.MMTol

	// Pipe everything into one sort/unique
	var c0 *exec.Cmd
	if sortTmpFlag != "" {
		c0 = exec.Command("sort", sortmem, sortpar, sortTmpFlag, "-u", "-")
	} else {
		c0 = exec.Command("sort", sortmem, sortpar, "-u", "-")
	}
	c0.Env = os.Environ()
	c0.Stderr = os.Stderr
	cmds := []*exec.Cmd{c0}

	// The sorted results go to disk
	outname := path.Join(config.TempDir, "matches.txt.sz")
	out, err := os.Create(outname)
	if err != nil {
		msg := "Error in combineWindows, see log files for details.\n"
		os.Stderr.WriteString(msg)
		log.Fatal(err)
	}
	defer out.Close()
	wtr := snappy.NewBufferedWriter(out)
	defer wtr.Close()

	// TODO: Add Bloom filter here to screen out duplicates
	var fd []io.Reader
	for j := 0; j < len(config.Windows); j++ {
		f := fmt.Sprintf("rmatch_%d.txt.sz", j)
		fname := path.Join(config.TempDir, f)
		c := exec.Command("sztool", "-d", fname)
		c.Env = os.Environ()
		c.Stderr = os.Stderr
		cmds = append(cmds, c)
		p, err := c.StdoutPipe()
		if err != nil {
			msg := "Error in combineWindows, see log files for details.\n"
			os.Stderr.WriteString(msg)
			log.Fatal(err)
		}
		fd = append(fd, p)
	}
	c0.Stdin = io.MultiReader(fd...)
	da, err := c0.StdoutPipe()
	if err != nil {
		msg := "Error in combineWindows, see log files for details.\n"
		os.Stderr.WriteString(msg)
		log.Fatal(err)
	}

	for _, c := range cmds {
		err := c.Start()
		if err != nil {
			msg := "Error in combineWindows, see log files for details.\n"
			os.Stderr.WriteString(msg)
			log.Fatal(err)
		}
	}

	// Taking all matches for the same read, retain only those
	// with nmiss equal to at most one greater than the lowest
	// nmiss.
	sem := make(chan bool, 1)
	sem <- true
	// DEBUG used to be go func()
	func() {
		scanner := bufio.NewScanner(da)
		var lines []string
		var fields [][]string
		var ibuf []int
		var current string
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
			ibuf, err = writebest(lines, fields, wtr, ibuf, mmtol)
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
			_, err := writebest(lines, fields, wtr, ibuf, mmtol)
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

		<-sem
	}()

	// OK to call Wait, done reading.
	for _, c := range cmds {
		err := c.Wait()
		if err != nil {
			msg := "Error in combineWindows, see log file for details.\n"
			os.Stderr.WriteString(msg)
			log.Fatal(err)
		}
	}
	sem <- true

	logger.Print("combineWindows done")
}
