// Copyright 2017, Kerby Shedden and the Muscato contributors.

// Muscato (Multi-Genome Scalable Alignment Tool) is a software tool
// for matching large collections of reads into large collections of
// target sequence (e.g. transcript sequences).
//
// Muscato uses a two-stage approach.  First, high-entropy
// subsequences of the reads are used to produce Bloom filter sketches
// of the read collection.  These sketches are used to identify
// candidate matches.  For example, if three offsets are chosen at
// positions 0, 20, and 40 of the reads, then three Bloom filter
// sketches are constructed.  Then, every window in the target
// sequence collection is queried against these sketches, identifying
// a set of candidate matches.  In the next step, for every
// subsequence appearing at each read offset, all reads and all genes
// containing the subsequence are assessed for pairwise similarity,
// and read/target pairs showing sufficiently high similarity are
// retained.
//
// This script is the entry point for the Muscato tool.  Normally,
// this is the only script that will be run directly.  It calls the
// other Muscato scripts in turn.
//
// Muscato can be invoked either using a configuration file in JSON
// format, or using command-line flags.  A typical invocation using
// flags is:
//
// muscato --ResultsFileName=results.txt --ReadFileName=reads.fastq --GeneFileName=genes.txt.sz --GeneIdFileName=1genes_ids.txt.sz
//    --Windows=0,20,40,60,80 --WindowWidth=15 --BloomSize=4000000000 --NumHash=20 --PMatch=0.96 --MinDinuc=5 --MinReadLength=50
//    --MaxMatches=1000000 --MaxMergeProcs=5 --MaxReadLength=300 --MatchMode=best --MMTol=2
//
// To use a JSON config file, create a file with the flag information in JSON format, e.g.
//
//    {ResultsFileName: "results.txt", MaxReadLength: 300, ...}
//
// Then provide the configuration file path when invoking Muscato, e.g.
//
// muscato --ConfigFileName=config.json
//
// Note that before running muscato, it is necessary to produce a
// processed version of the target sequence data.  This can be done
// using the muscato_prep_targets tool, invoked as follows.
//
// muscato_prep_targets genes.fasta
//
// See utils/Config.go for the full set of configuration parameters.
//
// Muscato generates a number of intermediate files and logs that by
// default are placed into the directory tmp/#####, where ##### is a
// generated number.  This temporary directory can be deleted after a
// successful run if desired.  The log files in the tmp directory may
// contain useful information for troubleshooting.
//
// Since Muscato uses Unix-style FIFOs for interprocess communication,
// it can only be run on Unix-like systems at present.  For the same
// reason, Muscato may not be runnable from AFS or NFS implementations
// that do not support FIFOs.

package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"math/rand"
	"os"
	"os/exec"
	"path"
	"strconv"
	"strings"
	"time"

	"github.com/golang/snappy"
	"github.com/kshedden/seqmatch/utils"
	"github.com/scipipe/scipipe"
	"github.com/willf/bloom"
	"golang.org/x/sys/unix"
)

var (
	jsonfile    string
	startpoint  int
	tmpjsonfile string
	config      *utils.Config
	basename    string
	tmpdir      string
	pipedir     string
	logger      *log.Logger
)

const (
	sortbuf string = "-S 2G"
	sortpar string = "--parallel=8"
)

func pipename() string {
	f := fmt.Sprintf("%09d", rand.Int63()%1e9)
	return path.Join(pipedir, f)
}

// pipefromsz creates a fifo and starts decompressing the given snappy
// file into it.
func pipefromsz(fname string) string {

	rand.Seed(int64(time.Now().UnixNano() + int64(os.Getpid())))

	for k := 0; k < 10; k++ {
		name := pipename()
		err := unix.Mkfifo(name, 0755)
		if err == nil {
			go func() {
				cmd := exec.Command("sztool", "-d", fname, name)
				cmd.Env = os.Environ()
				cmd.Stderr = os.Stderr
				err := cmd.Run()
				if err != nil {
					panic(err)
				}
			}()
			return name
		}
		print(fmt.Sprintf("%v\n", err))
	}

	panic("unable to create pipe")
}

func sortSource() {

	logger.Printf("starting sortSource")

	logger.Printf("Running prep_reads %s %s", tmpjsonfile, tmpdir)
	cmd0 := exec.Command("muscato_prep_reads", tmpjsonfile, tmpdir)
	cmd0.Env = os.Environ()
	cmd0.Stderr = os.Stderr

	cmd1 := exec.Command("sort", sortbuf, sortpar)
	cmd1.Env = os.Environ()
	cmd1.Stderr = os.Stderr
	var err error
	cmd1.Stdin, err = cmd0.StdoutPipe()
	if err != nil {
		panic(err)
	}
	pip, err := cmd1.StdoutPipe()
	if err != nil {
		panic(err)
	}

	cmds := []*exec.Cmd{cmd0, cmd1}

	for _, cmd := range cmds {
		err = cmd.Start()
		if err != nil {
			panic(err)
		}
	}

	scanner := bufio.NewScanner(pip)
	buf := make([]byte, 1024*1024)
	scanner.Buffer(buf, len(buf))

	// File for sequences
	outname := path.Join(tmpdir, "reads_sorted.txt.sz")
	logger.Printf("Writing sequences to %s", outname)
	fid, err := os.Create(outname)
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	wtr := snappy.NewBufferedWriter(fid)
	defer wtr.Close()

	// Get the first line
	if !scanner.Scan() {
		logger.Printf("no input")
		panic("no input (is the read file empty?)")
	}
	if err := scanner.Err(); err != nil {
		panic(err)
	}
	fields := strings.Fields(scanner.Text())
	seq := fields[0]
	name := []string{fields[1]}
	n := 1
	nseq := 0

	dowrite := func(seq string, name []string, n int) {
		xn := strings.Join(name, ";")
		if len(xn) > 1000 {
			xn = xn[0:995]
			xn += "..."
		}
		nseq++
		_, err = wtr.Write([]byte(seq))
		if err != nil {
			panic(err)
		}
		_, err = wtr.Write([]byte("\t"))
		if err != nil {
			panic(err)
		}
		s := fmt.Sprintf("%d\t%s\n", n, xn)
		_, err = wtr.Write([]byte(s))
		if err != nil {
			panic(err)
		}
	}

	for scanner.Scan() {
		line := scanner.Text()
		fields1 := strings.Fields(line)
		seq1 := fields1[0]
		name1 := fields1[1]

		if strings.Compare(seq, seq1) == 0 {
			n++
			name = append(name, name1)
			continue
		}

		dowrite(seq, name, n)
		seq = seq1
		name = name[0:1]
		name[0] = name1
		n = 1
	}

	if err := scanner.Err(); err != nil {
		panic(err)
	}

	// Read to EOF before calling wait.
	dowrite(seq, name, n)

	for _, cmd := range cmds {
		if err := cmd.Wait(); err != nil {
			log.Fatal(err)
		}
	}

	logger.Printf(fmt.Sprintf("Wrote %d read sequences", nseq))
	logger.Printf("sortSource done")
}

func windowReads() {
	logger.Printf("starting windowReads")

	cmd := exec.Command("muscato_window_reads", tmpjsonfile, tmpdir)
	cmd.Env = os.Environ()
	cmd.Stderr = os.Stderr
	err := cmd.Run()
	if err != nil {
		panic(err)
	}

	logger.Printf("windowReads done")
}

func sortWindows() {

	logger.Printf("starting sortWindows")

	for k := 0; k < len(config.Windows); k++ {
		f := fmt.Sprintf("win_%d.txt.sz", k)
		fname := path.Join(tmpdir, f)
		pname1 := pipefromsz(fname)

		cmd1 := exec.Command("sort", sortbuf, sortpar, "-k1", pname1)
		cmd1.Env = os.Environ()
		cmd1.Stderr = os.Stderr

		fname = strings.Replace(fname, ".txt.sz", "_sorted.txt.sz", 1)
		cmd2 := exec.Command("sztool", "-c", "-", fname)
		cmd2.Env = os.Environ()
		cmd2.Stderr = os.Stderr
		var err error
		cmd2.Stdin, err = cmd1.StdoutPipe()
		if err != nil {
			panic(err)
		}

		cmds := []*exec.Cmd{cmd2, cmd1}

		for _, cmd := range cmds {
			err := cmd.Start()
			if err != nil {
				panic(err)
			}
		}

		// Order is important here, need to wait on cmd2
		// before waiting on cmd1.
		for _, cmd := range cmds {
			err := cmd.Wait()
			if err != nil {
				panic(err)
			}
		}
	}

	logger.Printf("sortWindows done")
}

func screen() {
	logger.Printf("starting screening")

	cmd := exec.Command("muscato_screen", tmpjsonfile, tmpdir)
	cmd.Env = os.Environ()
	cmd.Stderr = os.Stderr
	err := cmd.Run()
	if err != nil {
		panic(err)
	}
	logger.Printf("screening done")
}

func sortBloom() {

	logger.Printf("starting sortBloom")

	for k := range config.Windows {
		f := fmt.Sprintf("bmatch_%d.txt.sz", k)
		fname := path.Join(tmpdir, f)
		pname1 := pipefromsz(fname)

		cmd1 := exec.Command("sort", sortbuf, sortpar, "-k1", pname1)
		cmd1.Env = os.Environ()
		cmd1.Stderr = os.Stderr

		f = fmt.Sprintf("smatch_%d.txt.sz", k)
		fname = path.Join(tmpdir, f)
		cmd2 := exec.Command("sztool", "-c", "-", fname)
		cmd2.Env = os.Environ()
		cmd2.Stderr = os.Stderr
		var err error
		cmd2.Stdin, err = cmd1.StdoutPipe()
		if err != nil {
			panic(err)
		}

		// Order matters here
		cmds := []*exec.Cmd{cmd2, cmd1}

		for _, cmd := range cmds {
			err := cmd.Start()
			if err != nil {
				panic(err)
			}
		}

		// Order is important, must wait on cmd2 first
		for _, cmd := range cmds {
			err := cmd.Wait()
			if err != nil {
				panic(err)
			}
		}
	}

	logger.Printf("sortBloom done")
}

func confirm() {

	logger.Printf("starting match confirmation")
	fp := 0
	for {
		nproc := config.MaxMergeProcs
		if nproc > len(config.Windows)-fp {
			nproc = len(config.Windows) - fp
		}
		if nproc == 0 {
			break
		}

		var cmds []*exec.Cmd
		for k := fp; k < fp+nproc; k++ {
			logger.Printf("Starting a round of confirmation processes")
			cmd := exec.Command("muscato_confirm", tmpjsonfile, fmt.Sprintf("%d", k), tmpdir)
			cmd.Env = os.Environ()
			cmd.Stderr = os.Stderr
			err := cmd.Start()
			if err != nil {
				panic(err)
			}
			cmds = append(cmds, cmd)
		}

		for _, cmd := range cmds {
			err := cmd.Wait()
			if err != nil {
				panic(err)
			}
		}
		fp += nproc
	}

	logger.Printf("match confirmation done")
}

// writebest accepts a set of lines (lines), which have also been
// broken into fields (bfr).  Every line represents a candidate match.
// The matches with at most mmtol more matches than the best match are
// written to the io writer (wtr).  ibuf is provided workspace.
func writebest(lines []string, bfr [][]string, wtr io.Writer, ibuf []int, mmtol int) []int {

	// Find the best fit, determine the number of mismatches for each sequence.
	ibuf = ibuf[0:0]
	best := -1
	for _, x := range bfr {
		y, err := strconv.Atoi(x[3]) // 3 is position of nmiss
		if err != nil {
			panic(err)
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
				panic(err)
			}
			_, err = wtr.Write([]byte("\n"))
			if err != nil {
				panic(err)
			}
		}
	}

	return ibuf
}

func combineWindows() {

	logger.Printf("starting combineWindows")

	mmtol := config.MMTol

	// Pipe everything into one sort/unique
	c0 := exec.Command("sort", sortbuf, sortpar, "-u", "-")
	c0.Env = os.Environ()
	c0.Stderr = os.Stderr
	cmds := []*exec.Cmd{c0}

	// The sorted results go to disk
	outname := path.Join(tmpdir, "matches.txt.sz")
	out, err := os.Create(outname)
	if err != nil {
		panic(err)
	}
	wtr := snappy.NewBufferedWriter(out)

	// TODO: Add Bloom filter here to screen out duplicates
	var fd []io.Reader
	for j := 0; j < len(config.Windows); j++ {
		f := fmt.Sprintf("rmatch_%d.txt.sz", j)
		fname := path.Join(tmpdir, f)
		c := exec.Command("sztool", "-d", fname)
		c.Env = os.Environ()
		c.Stderr = os.Stderr
		cmds = append(cmds, c)
		p, err := c.StdoutPipe()
		if err != nil {
			panic(err)
		}
		fd = append(fd, p)
	}
	c0.Stdin = io.MultiReader(fd...)
	da, err := c0.StdoutPipe()
	if err != nil {
		panic(err)
	}

	for _, c := range cmds {
		err := c.Start()
		if err != nil {
			panic(err)
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
			ibuf = writebest(lines, fields, wtr, ibuf, mmtol)
			lines = lines[0:0]
			lines = append(lines, line)
			fields = fields[0:0]
			fields = append(fields, field)
			current = field[0]
		}

		if err := scanner.Err(); err == nil {
			// Process the final block if possible
			writebest(lines, fields, wtr, ibuf, mmtol)
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
			panic(err)
		}
	}
	sem <- true

	wtr.Close()
	out.Close()

	logger.Printf("combineWindows done")
}

func sortByGeneId() {

	logger.Printf("starting sortByGeneid")
	inname := path.Join(tmpdir, "matches.txt.sz")
	outname := path.Join(tmpdir, "matches_sg.txt.sz")

	// Sort by gene number
	cmd1 := exec.Command("sztool", "-d", inname)
	cmd1.Env = os.Environ()
	cmd1.Stderr = os.Stderr
	// k5 is position of gene id
	cmd2 := exec.Command("sort", sortbuf, sortpar, "-k5", "-")
	cmd2.Env = os.Environ()
	cmd2.Stderr = os.Stderr
	var err error
	cmd2.Stdin, err = cmd1.StdoutPipe()
	if err != nil {
		panic(err)
	}
	cmd3 := exec.Command("sztool", "-c", "-", outname)
	cmd3.Env = os.Environ()
	cmd3.Stderr = os.Stderr
	cmd3.Stdin, err = cmd2.StdoutPipe()
	if err != nil {
		panic(err)
	}

	// Order matters
	cmds := []*exec.Cmd{cmd3, cmd2, cmd1}
	for _, c := range cmds {
		err := c.Start()
		if err != nil {
			panic(err)
		}
	}

	// Call Wait from end to beginning of chained commands
	for _, c := range cmds {
		err := c.Wait()
		if err != nil {
			panic(err)
		}
	}

	logger.Printf("sortbyGeneId done")
}

func joinGeneNames() {

	logger.Printf("starting joinGeneNames")

	// Decompress matches
	ma := scipipe.NewProc("ma", fmt.Sprintf("sztool -d %s > {os:ma}", path.Join(tmpdir, "matches_sg.txt.sz")))
	ma.SetPathStatic("ma", path.Join(pipedir, "jgn_ma.txt"))

	// Decompress gene ids
	gn := scipipe.NewProc("gn", fmt.Sprintf("sztool -d %s > {os:gn}", config.GeneIdFileName))
	gn.SetPathStatic("gn", path.Join(pipedir, "jgn_gn.txt"))

	// Join genes and matches
	jo := scipipe.NewProc("jo", "join -1 5 -2 1 -t'\t' {i:mx} {i:gx} > {os:jx}")
	jo.SetPathStatic("jx", path.Join(pipedir, "jgn_joined.txt"))

	// Cut out unwanted column
	ct := scipipe.NewProc("ct", "cut -d'\t' -f 1 --complement {i:jy} > {os:co}")
	ct.SetPathStatic("co", path.Join(pipedir, "jgn_cut.txt"))

	// Compress the result
	sz := scipipe.NewProc("sz", fmt.Sprintf("sztool -c {i:zi} %s", path.Join(tmpdir, "matches_sn.txt.sz")))

	jo.In("mx").Connect(ma.Out("ma"))
	jo.In("gx").Connect(gn.Out("gn"))
	ct.In("jy").Connect(jo.Out("jx"))
	sz.In("zi").Connect(ct.Out("co"))

	wf := scipipe.NewWorkflow("jgn")
	wf.AddProcs(ma, gn, jo, ct, sz)
	wf.SetDriver(sz)
	wf.Run()

	logger.Printf("joinGeneNames done")
}

func joinReadNames() {

	logger.Printf("starting joinReadNames")

	// The workflow hangs if the results file already exists, so
	// remove it.
	_, err := os.Stat(config.ResultsFileName)
	if err == nil {
		err := os.Remove(config.ResultsFileName)
		if err != nil {
			panic(err)
		}
	} else if os.IsNotExist(err) {
		// do nothing
	} else {
		panic(err)
	}

	// Decompress matches
	ma := scipipe.NewProc("ma", fmt.Sprintf("sztool -d %s > {os:ma}",
		path.Join(tmpdir, "matches_sn.txt.sz")))
	ma.SetPathStatic("ma", path.Join(pipedir, "jrn_ma.txt"))

	// Decompress sorted reads
	rd := scipipe.NewProc("rd", fmt.Sprintf("sztool -d %s > {os:rd}",
		path.Join(tmpdir, "reads_sorted.txt.sz")))
	rd.SetPathStatic("rd", path.Join(pipedir, "jrn_rd.txt"))

	// Sort the matches
	sm := scipipe.NewProc("sm", fmt.Sprintf("sort %s %s -k1 {i:in} > {os:sort}", sortbuf, sortpar))
	sm.SetPathStatic("sort", path.Join(pipedir, "jrn_sort.txt"))

	// Join the sorted matches with the reads
	jo := scipipe.NewProc("jo", "join -1 1 -2 1 -t'\t' {i:srx} {i:rdx} > {o:out}")
	jo.SetPathStatic("out", config.ResultsFileName)

	snk := scipipe.NewSink("snk")

	// Connect the network
	sm.In("in").Connect(ma.Out("ma"))
	jo.In("srx").Connect(sm.Out("sort"))
	jo.In("rdx").Connect(rd.Out("rd"))
	snk.Connect(jo.Out("out"))

	wf := scipipe.NewWorkflow("jrn")
	wf.AddProcs(ma, rd, sm, jo)
	wf.SetDriver(snk)
	wf.Run()

	logger.Printf("joinReadNames done")
}

func setupLog() {
	logname := path.Join(tmpdir, "muscato.log")
	fid, err := os.Create(logname)
	if err != nil {
		panic(err)
	}
	logger = log.New(fid, "", log.Ltime)
}

func copyconfig(config *utils.Config, tmpdir string) {

	fid, err := os.Create(path.Join(tmpdir, "config.json"))
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	enc := json.NewEncoder(fid)
	err = enc.Encode(config)
	if err != nil {
		panic(err)
	}
	tmpjsonfile = path.Join(tmpdir, "config.json")
}

func handleArgs() {

	ConfigFileName := flag.String("ConfigFileName", "", "JSON file containing configuration parameters")
	ReadFileName := flag.String("ReadFileName", "", "Sequencing read file (fastq format)")
	GeneFileName := flag.String("GeneFileName", "", "Gene file name (processed form)")
	GeneIdFileName := flag.String("GeneIdFileName", "", "Gene ID file name (processed form)")
	ResultsFileName := flag.String("ResultsFileName", "", "File name for results")
	WindowsRaw := flag.String("Windows", "", "Starting position of each window")
	WindowWidth := flag.Int("WindowWidth", 0, "Width of each window")
	BloomSize := flag.Int("BloomSize", 0, "Size of Bloom filter, in bits")
	NumHash := flag.Int("NumHash", 0, "Number of hashses")
	PMatch := flag.Float64("PMatch", 0, "Required proportion of matching positions")
	MinDinuc := flag.Int("MinDinuc", 0, "Minimum number of dinucleotides to check for match")
	TempDir := flag.String("TempDir", "", "Workspace for temporary files")
	MinReadLength := flag.Int("MinReadLength", 0, "Reads shorter than this length are skipped")
	MaxReadLength := flag.Int("MaxReadLength", 0, "Reads longer than this length are truncated")
	MaxMatches := flag.Int("MaxMatches", 0, "Return no more than this number of matches per window")
	MaxMergeProcs := flag.Int("MaxMergeProcs", 0, "Run this number of merge processes concurrently")
	MMTol := flag.Int("MMTol", 0, "Number of mismatches allowed above best fit")
	StartPoint := flag.Int("StartPoint", 0, "Restart at a given point in the procedure")
	MatchMode := flag.String("MatchMode", "", "'first' (retain first matches meeting criteria) or 'best' (returns best matches meeting criteria)")

	flag.Parse()

	if *ConfigFileName != "" {
		jsonfile = *ConfigFileName
		config = utils.ReadConfig(jsonfile)
	} else {
		config = new(utils.Config)
	}

	if *ReadFileName != "" {
		config.ReadFileName = *ReadFileName
	}
	if *GeneFileName != "" {
		config.GeneFileName = *GeneFileName
	}
	if *GeneIdFileName != "" {
		config.GeneIdFileName = *GeneIdFileName
	}
	if *WindowWidth != 0 {
		config.WindowWidth = *WindowWidth
	}
	if *BloomSize != 0 {
		config.BloomSize = uint64(*BloomSize)
	}
	if *NumHash != 0 {
		config.NumHash = *NumHash
	}
	if *PMatch != 0 {
		config.PMatch = *PMatch
	}
	if *MinDinuc != 0 {
		config.MinDinuc = *MinDinuc
	}
	if *TempDir != "" {
		config.TempDir = *TempDir
	}
	if *MinReadLength != 0 {
		config.MinReadLength = *MinReadLength
	}
	if *MaxReadLength != 0 {
		config.MaxReadLength = *MaxReadLength
	}
	if *MaxMatches != 0 {
		config.MaxMatches = *MaxMatches
	}
	if *MaxMergeProcs != 0 {
		config.MaxMergeProcs = *MaxMergeProcs
	}
	if *MatchMode != "" {
		config.MatchMode = *MatchMode
	}
	if *MMTol != 0 {
		config.MMTol = *MMTol
	}
	if *ResultsFileName != "" {
		config.ResultsFileName = *ResultsFileName
	}

	if config.ResultsFileName == "" {
		print("ResultsFileName must be specified.  Run 'muscato --help' for more information.\n\n")
		os.Exit(1)
	}

	startpoint = *StartPoint

	if *WindowsRaw != "" {
		toks := strings.Split(*WindowsRaw, ",")
		var itoks []int
		for _, x := range toks {
			y, err := strconv.Atoi(x)
			if err != nil {
				panic(err)
			}
			itoks = append(itoks, y)
		}
		config.Windows = itoks

		// TODO: this could be relaxed
		if config.Windows[0] != 0 {
			print("Windows[0] must be set to 0.\n\n")
			os.Exit(1)
		}
	} else {
		print("Windows must be specified.  Run `muscato --help` for more information.\n\n")
		os.Exit(1)
	}
}

func checkArgs() {

	if config.ReadFileName == "" {
		os.Stderr.WriteString("ReadFileName not provided, run 'muscato --help for more information.\n\n")
		os.Exit(1)
	}
	if config.GeneFileName == "" {
		os.Stderr.WriteString("GeneFileName not provided, run 'muscato --help for more information.\n\n")
		os.Exit(1)
	}
	if config.GeneIdFileName == "" {
		os.Stderr.WriteString("GeneIdFileName not provided, run 'muscato --help for more information.\n\n")
		os.Exit(1)
	}
	if config.ResultsFileName == "" {
		config.ResultsFileName = "results.txt"
		os.Stderr.WriteString("ResultsFileName not provided, defaulting to 'results.txt'\n\n")
	}
	if len(config.Windows) == 0 {
		os.Stderr.WriteString("Windows not provided, run 'muscato --help for more information.\n\n")
		os.Exit(1)
	}
	if config.WindowWidth == 0 {
		os.Stderr.WriteString("WindowWidth not provided, run 'muscato --help for more information.\n\n")
		os.Exit(1)
	}
	if config.BloomSize == 0 {
		os.Stderr.WriteString("BloomSize not provided, defaulting to 4 billion.\n\n")
		config.BloomSize = 4 * 1000 * 1000 * 1000
	}
	if config.NumHash == 0 {
		os.Stderr.WriteString("NumHash not provided, defaulting to 20.\n\n")
		config.NumHash = 20
	}
	if config.PMatch == 0 {
		os.Stderr.WriteString("PMatch not provided, defaulting to 1.\n\n")
		config.PMatch = 1
	}
	if config.MaxReadLength == 0 {
		os.Stderr.WriteString("MaxReadLength not provided, run 'muscato --help for more information.\n\n")
		os.Exit(1)
	}
	if config.MaxMatches == 0 {
		os.Stderr.WriteString("MaxMatches not provided, defaulting to 1 million\n\n")
		config.MaxMatches = 1000 * 1000
	}
	if config.MaxMergeProcs == 0 {
		os.Stderr.WriteString("MaxMergeProcs not provided, defaulting to 3\n\n")
		config.MaxMergeProcs = 3
	}
	if !strings.HasSuffix(config.ReadFileName, ".fastq") {
		msg := fmt.Sprintf("Warning: %s may not be a fastq file, continuing anyway\n\n",
			config.ReadFileName)
		os.Stderr.WriteString(msg)
	}
	if config.MatchMode == "" {
		os.Stderr.WriteString("MatchMode not provided, defaulting to 'best'\n")
		config.MatchMode = "best"
	}
}

func setupEnvs() {
	err := os.Setenv("LC_ALL", "C")
	if err != nil {
		panic(err)
	}
	home := os.Getenv("HOME")
	gopath := path.Join(home, "go")
	err = os.Setenv("GOPATH", gopath)
	if err != nil {
		panic(err)
	}
	err = os.Setenv("PATH", os.Getenv("PATH")+":"+home+"/go/bin")
	if err != nil {
		panic(err)
	}
}

// Create the directory for all temporary files, if needed
func makeTemp() {
	var d string
	var err error
	d, basename = path.Split(config.ReadFileName)
	if config.TempDir == "" {
		d = path.Join(d, "tmp")
		err = os.MkdirAll(d, 0755)
		if err != nil {
			panic(err)
		}
		tmpdir, err = ioutil.TempDir(d, "")
		if err != nil {
			panic(err)
		}
	} else {
		tmpdir = config.TempDir
		err = os.MkdirAll(tmpdir, 0755)
		if err != nil {
			panic(err)
		}
	}

	pipedir = path.Join(tmpdir, "pipes")
	err = os.MkdirAll(pipedir, 0755)
	if err != nil {
		panic(err)
	}
}

func writeNonMatch() {

	logger.Print("Starting writeNonMatch")

	// Reader for the match file
	inf, err := os.Open(config.ResultsFileName)
	if err != nil {
		panic(err)
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
		panic(err)
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
		panic(err)
	}
	defer out.Close()
	wtr := bufio.NewWriter(out)
	defer wtr.Flush()

	// Check each read to see if it was matched.
	rfname := path.Join(tmpdir, "reads_sorted.txt.sz")
	inf, err = os.Open(rfname)
	if err != nil {
		panic(err)
	}
	defer inf.Close()
	rdr := snappy.NewReader(inf)
	scanner = bufio.NewScanner(rdr)
	for scanner.Scan() {
		f := bytes.Fields(scanner.Bytes())
		if !bf.Test(f[0]) {
			_, err := wtr.Write(f[2])
			if err != nil {
				panic(err)
			}
			_, err = wtr.Write([]byte("\n"))
			if err != nil {
				panic(err)
			}
			_, err = wtr.Write(f[0])
			if err != nil {
				panic(err)
			}
			_, err = wtr.Write([]byte("\n+\n"))
			if err != nil {
				panic(err)
			}
			_, err = wtr.Write(bytes.Repeat([]byte{'!'}, len(f[0])))
			if err != nil {
				panic(err)
			}
			_, err = wtr.Write([]byte("\n"))
			if err != nil {
				panic(err)
			}
		}
	}

	logger.Printf("writeNonMatch done")
}

func run() {
	if startpoint <= 0 {
		sortSource()
	}

	if startpoint <= 1 {
		windowReads()
	}

	if startpoint <= 2 {
		sortWindows()
	}

	if startpoint <= 3 {
		screen()
	}

	if startpoint <= 4 {
		sortBloom()
	}

	if startpoint <= 5 {
		confirm()
	}

	if startpoint <= 6 {
		combineWindows()
	}

	if startpoint <= 7 {
		sortByGeneId()
	}

	if startpoint <= 8 {
		joinGeneNames()
	}

	if startpoint <= 9 {
		joinReadNames()
	}

	if startpoint <= 10 {
		writeNonMatch()
	}
}

func main() {

	handleArgs()
	checkArgs()
	setupEnvs()
	makeTemp()
	copyconfig(config, tmpdir)
	setupLog()

	logger.Printf("Storing temporary files in %s", tmpdir)

	run()

	logger.Printf("All done, exiting")
}
