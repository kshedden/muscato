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
//    --MaxMatches=1000000 --MaxConfirmProcs=5 --MaxReadLength=300 --MatchMode=best --MMTol=2
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
	"io/ioutil"
	"log"
	"math"
	"os"
	"os/exec"
	"path"
	"strconv"
	"strings"

	"github.com/golang/snappy"
	"github.com/google/uuid"
	"github.com/kshedden/muscato/utils"
	"github.com/scipipe/scipipe"
	"github.com/willf/bloom"
)

var (
	configFilePath string

	wf *scipipe.Workflow

	config   *utils.Config
	basename string
	pipedir  string
	logger   *log.Logger

	// Flag for setting the tmp file location for sorting.
	sortTmpFlag string

	sortpar string
	sortmem string
)

// geneStats
func geneStats() {

	c := fmt.Sprintf("sort %s %s %s %s -k 5 > {os:outsort}", sortmem, sortpar, sortTmpFlag, config.ResultsFileName)
	logger.Print(c)
	gsrt := wf.NewProc("gsrt", c)
	gsrt.SetOut("outsort", path.Join(pipedir, "gs_outsort"))

	var outfile string
	ext := path.Ext(config.ResultsFileName)
	if ext != "" {
		m := len(config.ResultsFileName)
		outfile = config.ResultsFileName[0:m-len(ext)] + "_genestats" + ext
	} else {
		outfile = config.ResultsFileName + "_genestats"
	}
	os.Remove(outfile)

	gsta := wf.NewProc("gst", "muscato_genestats {i:insort} > {o:out}")
	gsta.SetOut("out", outfile)
	gsta.In("insort").From(gsrt.Out("outsort"))
}

func prepReads() {

	logger.Print("Starting prepReads")

	// Run muscato_prep_reads
	dc := wf.NewProc("mpr", fmt.Sprintf("muscato_prep_reads %s > {os:mpr_out}", configFilePath))
	dc.SetOut("mpr_out", path.Join(pipedir, "pr_mpr"))

	// Sort the output of muscato_prep_reads
	c := fmt.Sprintf("sort %s %s %s {i:insort} > {os:outsort}", sortmem, sortpar, sortTmpFlag)
	logger.Print(c)
	sr := wf.NewProc("sr", c)
	sr.SetOut("outsort", path.Join(pipedir, "pr_outsort"))

	// Uniqify and count duplicates
	c = fmt.Sprintf("muscato_uniqify %s {i:inuniq} > {o:outfinal}", configFilePath)
	logger.Print(c)
	mu := wf.NewProc("mu", c)
	outname := path.Join(config.TempDir, "reads_sorted.txt.sz")
	mu.SetOut("outfinal", outname)

	// Connect the network
	sr.In("insort").From(dc.Out("mpr_out"))
	mu.In("inuniq").From(sr.Out("outsort"))

	logger.Print("prepReads done")
}

func windowReads() {
	logger.Print("starting windowReads")
	logger.Printf("Running command: 'muscato_window_reads %s'\n", configFilePath)
	wf.NewProc("mwr", fmt.Sprintf("muscato_window_reads %s", configFilePath))
	logger.Print("windowReads done")
}

func sortWindows() {

	logger.Print("starting sortWindows")

	for k := 0; k < len(config.Windows); k++ {

		logger.Printf("sortWindows %d...", k)

		// Decompress matches
		fn := path.Join(config.TempDir, fmt.Sprintf("win_%d.txt.sz", k))
		dw := fmt.Sprintf("dec_win_%d", k)
		dc := wf.NewProc(dw, fmt.Sprintf("sztool -d %s > {os:dx}", fn))
		dc.SetOut("dx", path.Join(pipedir, fmt.Sprintf("sw_dc_%d", k)))

		// Sort the matches
		sc := fmt.Sprintf("sort %s %s -k1 %s {i:in} > {os:sort}", sortmem, sortpar, sortTmpFlag)
		smn := fmt.Sprintf("swin_%d", k)
		sm := wf.NewProc(smn, sc)
		logger.Print(sc)
		sm.SetOut("sort", path.Join(pipedir, fmt.Sprintf("sw_sort_%d", k)))

		// Compress results
		fn = strings.Replace(fn, ".txt.sz", "_sorted.txt.sz", 1)
		rwn := fmt.Sprintf("rec_win_%d", k)
		rc := wf.NewProc(rwn, fmt.Sprintf("sztool -c {i:ins} %s", fn))

		// Connect the network
		sm.In("in").From(dc.Out("dx"))
		rc.In("ins").From(sm.Out("sort"))

		logger.Print("done\n")
	}

	logger.Print("sortWindows done")
}

func screen() {
	logger.Print("Starting screening")
	logger.Printf("Running command: 'muscato_screen %s'\n", configFilePath)
	wf.NewProc("mscr", fmt.Sprintf("muscato_screen %s", configFilePath))
	logger.Print("Screening done")
}

func sortBloom() {

	logger.Print("Starting sortBloom")

	for k := range config.Windows {

		logger.Printf("sortBloom %d...", k)

		// Decompress matches
		fn := path.Join(config.TempDir, fmt.Sprintf("bmatch_%d.txt.sz", k))
		dcn := fmt.Sprintf("dcb_%d", k)
		dc := wf.NewProc(dcn, fmt.Sprintf("sztool -d %s > {os:dx}", fn))
		dc.SetOut("dx", path.Join(pipedir, fmt.Sprintf("sb_dc_%d", k)))

		// Sort the matches
		c := fmt.Sprintf("sort %s %s -k1 %s {i:in} > {os:sort}", sortmem, sortpar, sortTmpFlag)
		logger.Print(c)
		smn := fmt.Sprintf("sb_%d", k)
		sm := wf.NewProc(smn, c)
		sm.SetOut("sort", path.Join(pipedir, fmt.Sprintf("sb_sort_%d", k)))

		// Compress results
		fn = path.Join(config.TempDir, fmt.Sprintf("smatch_%d.txt.sz", k))
		rcn := fmt.Sprintf("rc_%d", k)
		rc := wf.NewProc(rcn, fmt.Sprintf("sztool -c {i:ins} %s", fn))

		// Connect the network
		sm.In("in").From(dc.Out("dx"))
		rc.In("ins").From(sm.Out("sort"))

		logger.Print("done")
	}

	logger.Print("sortBloom done")
}

func confirm() {

	logger.Print("Starting match confirmation")
	for k := 0; k < len(config.Windows); k++ {
		logger.Printf("Running command: 'muscato_confirm %s %d'\n", configFilePath, k)
		crn := fmt.Sprintf("mc_%d", k)
		wf.NewProc(crn, fmt.Sprintf("muscato_confirm %s %d", configFilePath, k))
	}

	logger.Print("Match confirmation done")
}

func combineWindows() {

	logger.Print("Starting combineWindows")
	logger.Printf("Running command: 'muscato_combine_windows %s'\n", configFilePath)
	wf.NewProc("cwn", fmt.Sprintf("muscato_combine_windows %s", configFilePath))
	logger.Print("combineWindows done")
}

func sortByGeneId() {

	logger.Print("starting sortByGeneid")
	inname := path.Join(config.TempDir, "matches.txt.sz")
	outname := path.Join(config.TempDir, "matches_sg.txt.sz")

	// Sort by gene number
	cmd1 := exec.Command("sztool", "-d", inname)
	cmd1.Env = os.Environ()
	cmd1.Stderr = os.Stderr
	// k5 is position of gene id
	var cmd2 *exec.Cmd
	if sortTmpFlag != "" {
		cmd2 = exec.Command("sort", sortmem, sortpar, sortTmpFlag, "-k5", "-")
	} else {
		cmd2 = exec.Command("sort", sortmem, sortpar, "-k5", "-")
	}
	cmd2.Env = os.Environ()
	cmd2.Stderr = os.Stderr
	var err error
	cmd2.Stdin, err = cmd1.StdoutPipe()
	if err != nil {
		msg := "Error in sortByGeneId, see log files for details.\n"
		os.Stderr.WriteString(msg)
		log.Fatal(err)
	}
	cmd3 := exec.Command("sztool", "-c", "-", outname)
	cmd3.Env = os.Environ()
	cmd3.Stderr = os.Stderr
	cmd3.Stdin, err = cmd2.StdoutPipe()
	if err != nil {
		msg := "Error in sortByGeneId, see log files for details.\n"
		os.Stderr.WriteString(msg)
		log.Fatal(err)
	}

	// Order matters
	cmds := []*exec.Cmd{cmd3, cmd2, cmd1}
	for _, c := range cmds {
		err := c.Start()
		if err != nil {
			msg := "Error in sortByGeneId, see log files for details.\n"
			os.Stderr.WriteString(msg)
			log.Fatal(err)
		}
	}

	// Call Wait from end to beginning of chained commands
	for _, c := range cmds {
		err := c.Wait()
		if err != nil {
			msg := "Error in sortByGeneId, see log files for details.\n"
			os.Stderr.WriteString(msg)
			log.Fatal(err)
		}
	}

	logger.Print("sortbyGeneId done")
}

func joinGeneNames() {

	logger.Print("starting joinGeneNames")

	// Decompress matches
	ma := wf.NewProc("ma", fmt.Sprintf("sztool -d %s > {os:ma}", path.Join(config.TempDir, "matches_sg.txt.sz")))
	ma.SetOut("ma", path.Join(pipedir, "jgn_ma.txt"))

	// Decompress gene ids
	gn := wf.NewProc("gn", fmt.Sprintf("sztool -d %s > {os:gn}", config.GeneIdFileName))
	gn.SetOut("gn", path.Join(pipedir, "jgn_gn.txt"))

	// Join genes and matches
	jo := wf.NewProc("jo", "join -1 5 -2 1 -t'\t' {i:mx} {i:gx} > {os:jx}")
	jo.SetOut("jx", path.Join(pipedir, "jgn_joined.txt"))

	// Cut out unwanted column
	ct := wf.NewProc("ct", "cut -d'\t' -f 1 --complement {i:jy} > {os:co}")
	ct.SetOut("co", path.Join(pipedir, "jgn_cut.txt"))

	// Compress the result
	sz := wf.NewProc("sz", fmt.Sprintf("sztool -c {i:zi} %s", path.Join(config.TempDir, "matches_sn.txt.sz")))

	jo.In("mx").From(ma.Out("ma"))
	jo.In("gx").From(gn.Out("gn"))
	ct.In("jy").From(jo.Out("jx"))
	sz.In("zi").From(ct.Out("co"))

	logger.Print("joinGeneNames done")
}

func joinReadNames() {

	logger.Print("starting joinReadNames")

	// The workflow hangs if the results file already exists, so
	// remove it.
	_, err := os.Stat(config.ResultsFileName)
	if err == nil {
		err := os.Remove(config.ResultsFileName)
		if err != nil {
			msg := "Error in joinReadNames, see log files for details.\n"
			os.Stderr.WriteString(msg)
			log.Fatal(err)
		}
	} else if os.IsNotExist(err) {
		// do nothing
	} else {
		msg := "Error in joinReadNames, see log files for details.\n"
		os.Stderr.WriteString(msg)
		log.Fatal(err)
	}

	// Decompress matches
	ma := wf.NewProc("ma", fmt.Sprintf("sztool -d %s > {os:ma}",
		path.Join(config.TempDir, "matches_sn.txt.sz")))
	ma.SetOut("ma", path.Join(pipedir, "jrn_ma.txt"))

	// Decompress sorted reads
	rd := wf.NewProc("rd", fmt.Sprintf("sztool -d %s > {os:rd}",
		path.Join(config.TempDir, "reads_sorted.txt.sz")))
	rd.SetOut("rd", path.Join(pipedir, "jrn_rd.txt"))

	// Sort the matches
	sm := wf.NewProc("sm", fmt.Sprintf("sort %s %s -k1 %s {i:in} > {os:sort}", sortmem, sortpar, sortTmpFlag))
	sm.SetOut("sort", path.Join(pipedir, "jrn_sort.txt"))

	// Join the sorted matches with the reads
	jo := wf.NewProc("jo", "join -1 1 -2 1 -t'\t' {i:srx} {i:rdx} > {o:out}")
	jo.SetOut("out", config.ResultsFileName)

	// Connect the network
	sm.In("in").From(ma.Out("ma"))
	jo.In("srx").From(sm.Out("sort"))
	jo.In("rdx").From(rd.Out("rd"))

	logger.Print("joinReadNames done")
}

func setupLog() {
	logname := path.Join(config.LogDir, "muscato.log")
	fid, err := os.Create(logname)
	if err != nil {
		msg := fmt.Sprintf("Error creating %s, see log files for details.\n", logname)
		os.Stderr.WriteString(msg)
		log.Fatal(err)
	}
	logger = log.New(fid, "", log.Ltime)
}

// saveConfig saves the configuration file in json format into the log
// directory.
func saveConfig(config *utils.Config) {

	fid, err := os.Create(path.Join(config.LogDir, "config.json"))
	if err != nil {
		msg := "Error in saveConfig, see log files for details."
		os.Stderr.WriteString(msg)
		log.Fatal(err)
	}
	defer fid.Close()
	enc := json.NewEncoder(fid)
	err = enc.Encode(config)
	if err != nil {
		msg := "Error in saveConfig, see log files for details."
		os.Stderr.WriteString(msg)
		log.Fatal(err)
	}
	configFilePath = path.Join(config.LogDir, "config.json")
}

func handleArgs() {

	// Can't use logger in here because it has not been created yet, write errors/warnings to Stderr.

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
	MaxConfirmProcs := flag.Int("MaxConfirmProcs", 0, "Run this number of match confirmation processes concurrently")
	MMTol := flag.Int("MMTol", 0, "Number of mismatches allowed above best fit")
	MatchMode := flag.String("MatchMode", "", "'first' or 'best' (retain first/best 'MaxMatches' matches meeting criteria)")
	NoCleanTemp := flag.Bool("NoCleanTemp", false, "Do not delete temporary files from TempDir")
	SortPar := flag.Int("SortPar", 0, "Number of parallel sort processes")
	SortTemp := flag.String("SortTemp", "", "Directory to use for sort temp files")
	SortMem := flag.String("SortMem", "", "Gnu sort -S parameter")
	CPUProfile := flag.Bool("CPUProfile", false, "Capture CPU profile data")

	flag.Parse()

	if *ConfigFileName != "" {
		config = utils.ReadConfig(*ConfigFileName)
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
	if *MaxConfirmProcs != 0 {
		config.MaxConfirmProcs = *MaxConfirmProcs
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
	if *NoCleanTemp {
		config.NoCleanTemp = true
	}
	if *CPUProfile {
		config.CPUProfile = true
	}
	if *SortPar != 0 {
		config.SortPar = *SortPar
	}
	if *SortMem != "" {
		config.SortMem = *SortMem
	}

	// Configure the temporary directory for sort.
	if *SortTemp != "" {
		config.SortTemp = *SortTemp
		os.MkdirAll(config.SortTemp, 0770)
	}
	if config.SortTemp != "" {
		sortTmpFlag = fmt.Sprintf("--temporary-directory=%s", config.SortTemp)
	}

	if config.ResultsFileName == "" {
		config.ResultsFileName = "results.txt"
		os.Stderr.WriteString("ResultsFileName not specified, defaulting to 'results.txt'\n")
	}

	if *WindowsRaw != "" {
		toks := strings.Split(*WindowsRaw, ",")
		var itoks []int
		for _, x := range toks {
			y, err := strconv.Atoi(x)
			if err != nil {
				msg := "Error in handleArgs, see log files for details.\n"
				os.Stderr.WriteString(msg)
				log.Fatal(err)
			}
			itoks = append(itoks, y)
		}
		config.Windows = itoks
	}
}

func checkArgs() {

	if config.ReadFileName == "" {
		os.Stderr.WriteString("\nReadFileName not provided, run 'muscato --help for more information.\n\n")
		os.Exit(1)
	}
	if config.GeneFileName == "" {
		os.Stderr.WriteString("\nGeneFileName not provided, run 'muscato --help for more information.\n\n")
		os.Exit(1)
	}
	if config.GeneIdFileName == "" {
		os.Stderr.WriteString("\nGeneIdFileName not provided, run 'muscato --help for more information.\n\n")
		os.Exit(1)
	}
	if config.ResultsFileName == "" {
		config.ResultsFileName = "results.txt"
		os.Stderr.WriteString("ResultsFileName not provided, defaulting to 'results.txt'\n")
	}
	if len(config.Windows) == 0 {
		os.Stderr.WriteString("\nWindows not provided, run 'muscato --help for more information.\n\n")
		os.Exit(1)
	}
	if config.WindowWidth == 0 {
		os.Stderr.WriteString("\nWindowWidth not provided, run 'muscato --help for more information.\n\n")
		os.Exit(1)
	}
	if config.BloomSize == 0 {
		os.Stderr.WriteString("BloomSize not provided, defaulting to 4 billion\n")
		config.BloomSize = 4 * 1000 * 1000 * 1000
	}
	if config.NumHash == 0 {
		os.Stderr.WriteString("NumHash not provided, defaulting to 20\n")
		config.NumHash = 20
	}
	if config.PMatch == 0 {
		os.Stderr.WriteString("PMatch not provided, defaulting to 1\n")
		config.PMatch = 1
	}
	if config.MaxReadLength == 0 {
		os.Stderr.WriteString("MaxReadLength not provided, run 'muscato --help for more information.\n\n")
		os.Exit(1)
	}
	if config.MaxMatches == 0 {
		os.Stderr.WriteString("MaxMatches not provided, defaulting to 1 million\n")
		config.MaxMatches = 1000 * 1000
	}
	if config.MaxConfirmProcs == 0 {
		os.Stderr.WriteString("MaxConfirmProcs not provided, defaulting to 3\n")
		config.MaxConfirmProcs = 3
	}
	if !strings.HasSuffix(config.ReadFileName, ".fastq") {
		msg := fmt.Sprintf("Warning: %s may not be a fastq file, continuing anyway\n",
			config.ReadFileName)
		os.Stderr.WriteString(msg)
	}
	if config.MatchMode == "" {
		os.Stderr.WriteString("MatchMode not provided, defaulting to 'best'\n")
		config.MatchMode = "best"
	}

	if config.SortPar == 0 {
		// warning not needed
		config.SortPar = 8
	}
	sortpar = fmt.Sprintf("--parallel=%d", config.SortPar)

	if config.SortMem == "" {
		os.Stderr.WriteString("SortMem not provided, defaulting to 50%\n")
		config.SortMem = "50%"
	}
	sortmem = fmt.Sprintf("-S %s", config.SortMem)
}

func setupEnvs() {
	err := os.Setenv("LC_ALL", "C")
	if err != nil {
		msg := "Error in setupEnvs, see log files for details.\n"
		os.Stderr.WriteString(msg)
		log.Fatal(err)
	}
	home := os.Getenv("HOME")
	gopath := path.Join(home, "go")
	err = os.Setenv("GOPATH", gopath)
	if err != nil {
		msg := "Error in setupEnvs, see log files for details.\n"
		os.Stderr.WriteString(msg)
		log.Fatal(err)
	}
	err = os.Setenv("PATH", os.Getenv("PATH")+":"+home+"/go/bin")
	if err != nil {
		msg := "Error in setupEnvs, see log files for details.\n"
		os.Stderr.WriteString(msg)
		log.Fatal(err)
	}
}

// Create the directory for all temporary files, if needed
func makeTemp() {

	// temp files, log files, etc. are stored in directories defined by this unique id.
	xuid, err := uuid.NewUUID()
	if err != nil {
		msg := "Error in makeTemp, see log files for details.\n"
		os.Stderr.WriteString(msg)
		log.Fatal(err)
	}
	uid := xuid.String()

	if config.TempDir == "" {
		config.TempDir = path.Join("muscato_tmp", uid)
	} else {
		// Overwrite the provided TempDir with a subdirectory.
		config.TempDir = path.Join(config.TempDir, uid)
	}
	err = os.MkdirAll(config.TempDir, 0770)
	if err != nil {
		if os.IsNotExist(err) {
			msg := fmt.Sprintf("Directory %s does not exist and cannot be created.", config.TempDir)
			os.Stderr.WriteString(msg)
			os.Exit(1)
		}
		log.Fatal(err)
	}

	// The directory where all pipes are written, needs to be in a
	// filesystem that supports pipes.
	pipedir, err = ioutil.TempDir("/tmp", "muscato-pipes-")
	if err != nil {
		if os.IsNotExist(err) {
			msg := "Cannot create temporary directory in /tmp for pipes."
			os.Stderr.WriteString(msg)
			os.Exit(1)
		}
		log.Fatal(err)
	}

	// Setup the directory for logging.
	if config.LogDir == "" {
		config.LogDir = "muscato_logs"
	}
	config.LogDir = path.Join(config.LogDir, uid)
	err = os.MkdirAll(config.LogDir, 0770)
	if err != nil {
		if os.IsNotExist(err) {
			msg := fmt.Sprintf("Cannot create directory %s for log files.", config.LogDir)
			os.Stderr.WriteString(msg)
			os.Exit(1)
		}
		log.Fatal(err)
	}

}

func writeNonMatch() {

	logger.Print("Starting writeNonMatch")

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

	logger.Print("writeNonMatch done")
}

// readStats calculates statistics for each read, using a results
// datafile that is sorted by read.
func readStats() {

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

func run() {
	prepReads()
	windowReads()
	sortWindows()
	screen()
	sortBloom()
	confirm()
	combineWindows()
	sortByGeneId()
	joinGeneNames()
	joinReadNames()
}

func cleanPipes() {
	logger.Printf("Removing pipes from %s", pipedir)
	err := os.RemoveAll(pipedir)
	if err != nil {
		logger.Print("Can't remove pipes:")
		logger.Print(err)
		logger.Print("Continuing anyway...\n")
	}
}

func cleanTmp() {
	if !config.NoCleanTemp {
		logger.Printf("Removing temporary files from %s", config.TempDir)
		err := os.RemoveAll(config.TempDir)
		if err != nil {
			logger.Print("Can't remove temporary files:")
			logger.Print(err)
			logger.Print("Continuing anyway...\n")
		}
	}
}

func perfInfo() {

	// Number of possible k-mers
	wf := math.Pow(4, float64(config.WindowWidth))

	var inf struct {
		NumUnique int
		NumTotal  int
	}

	fid, err := os.Open(path.Join(config.LogDir, "seqinfo.json"))
	if err != nil {
		msg := "Error in perfInfo, see log files for details.\n"
		os.Stderr.WriteString(msg)
		log.Fatal(err)
	}

	dec := json.NewDecoder(fid)
	err = dec.Decode(&inf)
	if err != nil {
		msg := "Error in perfInfo, see log files for details.\n"
		os.Stderr.WriteString(msg)
		log.Fatal(err)
	}

	// Fill rate of the k-mer set inclusion function.
	ff := float64(inf.NumUnique) / wf

	logger.Printf("k-mer sketch fill rate: %.5f", ff)
}

func main() {

	wf = scipipe.NewWorkflow("muscato", 4)

	defer cleanPipes()
	defer cleanTmp()

	handleArgs()
	checkArgs()
	setupEnvs()
	makeTemp()
	saveConfig(config)
	setupLog()

	logger.Printf("Storing temporary files in %s", config.TempDir)
	logger.Printf("Storing pipes in %s", pipedir)
	logger.Printf("Storing log files in %s", config.LogDir)

	run()
	wf.Run()

	perfInfo()
	writeNonMatch()
	readStats()
	geneStats()

	logger.Print("All done, exit after cleanup")
}