// Copyright 2017, Kerby Shedden and the Muscato contributors.
//
// Muscato (Multi-Genome Scalable Alignment Tool) is a software tool
// for matching large collections of sequencing reads into large
// collections of target sequence (e.g. exon or gene sequences).
//
// Muscato uses a two-stage approach:
//
// 1. High entropy subsequences of the reads are used to produce Bloom
// filter sketches of the read collection.  Separate sketches are
// obtained for each of several read position offsets, e.g. 20-mers at
// positions 20, 40, and 60 in the read. Every window in the target
// sequence collection is then queried against these sketches,
// identifying a set of candidate matches.
//
// 2. The reads and candidate matching target sequences are sorted for
// each offset position, allowing read/target pairs showing
// sufficiently high similarity to be retained.
//
// This script is the entry point for the Muscato tool.  Normally,
// this is the only script that will be run directly.  It calls the
// other Muscato scripts.  All muscato scripts begin with `muscato_`
// and are installed to the GOBIN directory.
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
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"path"
	"strconv"
	"strings"

	"github.com/google/uuid"
	"github.com/kshedden/muscato/utils"
	"golang.org/x/sys/unix"
)

var (
	configFilePath string

	config   *utils.Config
	basename string

	// Flag for setting the tmp file location for sorting.
	sortTmpFlag string

	logger *log.Logger

	sortpar string
	sortmem string
)

// geneStats
func geneStats() {

	io.WriteString(os.Stderr, "Generating gene statistics...\n")

	pr1, pw1, err := os.Pipe()
	if err != nil {
		panic(err)
	}

	args := []string{sortmem, sortpar, "-k5"}
	if sortTmpFlag != "" {
		args = append(args, sortTmpFlag)
	}
	args = append(args, config.ResultsFileName)
	cmd1 := exec.Command("sort", args...)
	cmd1.Stdout = pw1

	var outfile string
	ext := path.Ext(config.ResultsFileName)
	if ext != "" {
		m := len(config.ResultsFileName)
		outfile = config.ResultsFileName[0:m-len(ext)] + "_genestats" + ext
	} else {
		outfile = config.ResultsFileName + "_genestats"
	}

	cmd2 := exec.Command("muscato_genestats", "-")
	cmd2.Stdin = pr1
	fid, err := os.Create(outfile)
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	cmd2.Stdout = fid

	for _, c := range []*exec.Cmd{cmd1, cmd2} {
		c.Stderr = os.Stderr
		if err := c.Start(); err != nil {
			panic(err)
		}
	}

	if err := cmd1.Wait(); err != nil {
		panic(err)
	}

	pw1.Close()

	if err := cmd2.Wait(); err != nil {
		panic(err)
	}
}

func mkfifo(pa string) *os.File {

	err := unix.Mkfifo(pa, 0600)
	if err != nil {
		panic(err)
	}

	file, err := os.OpenFile(pa, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0777)
	if err != nil {
		panic(err)
	}

	return file
}

func prepReads() {

	io.WriteString(os.Stderr, "Preparing reads...\n")

	pr1, pw1, err := os.Pipe()
	if err != nil {
		panic(err)
	}

	pr2, pw2, err := os.Pipe()
	if err != nil {
		panic(err)
	}

	// Run muscato_prep_reads
	cmd1 := exec.Command("muscato_prep_reads", configFilePath)
	cmd1.Stdout = pw1

	// Sort the output of muscato_prep_reads
	args := []string{sortmem, sortpar}
	if sortTmpFlag != "" {
		args = append(args, sortTmpFlag)
	}
	cmd2 := exec.Command("sort", args...)
	cmd2.Stdin = pr1
	cmd2.Stdout = pw2

	// Uniqify and count duplicates
	outfinal := path.Join(config.TempDir, "reads_sorted.txt.sz")
	cmd3 := exec.Command("muscato_uniqify", configFilePath, "-")
	cmd3.Stdin = pr2
	fid, err := os.Create(outfinal)
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	cmd3.Stdout = fid

	for _, cmd := range []*exec.Cmd{cmd1, cmd2, cmd3} {
		cmd.Stderr = os.Stderr
		if err := cmd.Start(); err != nil {
			panic(err)
		}
	}

	if err := cmd1.Wait(); err != nil {
		panic(err)
	}

	pw1.Close()

	if err := cmd2.Wait(); err != nil {
		panic(err)
	}

	pw2.Close()

	if err := cmd3.Wait(); err != nil {
		panic(err)
	}
}

func windowReads() {

	io.WriteString(os.Stderr, "Windowing reads...\n")

	// Run muscato_prep_reads
	cmd := exec.Command("muscato_window_reads", configFilePath)
	cmd.Stderr = os.Stderr

	if err := cmd.Run(); err != nil {
		panic(err)
	}
}

func sortWindows() {

	for k := 0; k < len(config.Windows); k++ {

		io.WriteString(os.Stderr, fmt.Sprintf("Sorting windows %d...\n", k))

		pr1, pw1, err := os.Pipe()
		if err != nil {
			panic(err)
		}

		pr2, pw2, err := os.Pipe()
		if err != nil {
			panic(err)
		}

		// Decompress matches
		fn := path.Join(config.TempDir, fmt.Sprintf("win_%d.txt.sz", k))
		cmd1 := exec.Command("sztool", "-d", fn)
		cmd1.Stdout = pw1

		// Sort the matches
		args := []string{sortmem, sortpar, "-k1"}
		if sortTmpFlag != "" {
			args = append(args, sortTmpFlag)
		}
		args = append(args, "-")
		cmd2 := exec.Command("sort", args...)
		cmd2.Stdin = pr1
		cmd2.Stdout = pw2

		// Compress results
		fn = strings.Replace(fn, ".txt.sz", "_sorted.txt.sz", 1)
		cmd3 := exec.Command("sztool", "-c", "-", fn)
		cmd3.Stdin = pr2

		for _, cmd := range []*exec.Cmd{cmd1, cmd2, cmd3} {
			cmd.Stderr = os.Stderr
			if err := cmd.Start(); err != nil {
				panic(err)
			}
		}

		if err := cmd1.Wait(); err != nil {
			panic(err)
		}

		pw1.Close()

		if err := cmd2.Wait(); err != nil {
			panic(err)
		}

		pw2.Close()

		if err := cmd3.Wait(); err != nil {
			panic(err)
		}
	}
}

func screen() {

	io.WriteString(os.Stderr, "Screening...\n")

	cmd := exec.Command("muscato_screen", configFilePath)
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		panic(err)
	}
}

func sortBloom() {

	for k := range config.Windows {

		pr1, pw1, err := os.Pipe()
		if err != nil {
			panic(err)
		}

		pr2, pw2, err := os.Pipe()
		if err != nil {
			panic(err)
		}

		io.WriteString(os.Stderr, fmt.Sprintf("Sorting Bloom %d...\n", k))

		// Decompress matches
		fn := path.Join(config.TempDir, fmt.Sprintf("bmatch_%d.txt.sz", k))
		cmd1 := exec.Command("sztool", "-d", fn)
		cmd1.Stdout = pw1

		// Sort the matches
		args := []string{sortmem, sortpar, "-k1"}
		if sortTmpFlag != "" {
			args = append(args, sortTmpFlag)
		}
		args = append(args, "-")
		cmd2 := exec.Command("sort", args...)
		cmd2.Stdin = pr1
		cmd2.Stdout = pw2

		// Compress results
		fn = path.Join(config.TempDir, fmt.Sprintf("smatch_%d.txt.sz", k))
		cmd3 := exec.Command("sztool", "-c", "-", fn)
		cmd3.Stdin = pr2

		for _, cmd := range []*exec.Cmd{cmd1, cmd2, cmd3} {
			cmd.Stderr = os.Stderr
			if err := cmd.Start(); err != nil {
				panic(err)
			}
		}

		if err := cmd1.Wait(); err != nil {
			panic(err)
		}

		pw1.Close()

		if err := cmd2.Wait(); err != nil {
			panic(err)
		}

		pw2.Close()

		if err := cmd3.Wait(); err != nil {
			panic(err)
		}
	}
}

func confirm() {

	io.WriteString(os.Stderr, "Confirming...\n")

	for k := 0; k < len(config.Windows); k++ {
		cmd := exec.Command("muscato_confirm", configFilePath, fmt.Sprintf("%d", k))
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			panic(err)
		}
	}
}

func combineWindows() {

	io.WriteString(os.Stderr, "Combining windows...\n")

	cmd := exec.Command("muscato_combine_windows", configFilePath)
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		panic(err)
	}
}

func sortByGeneId() {

	io.WriteString(os.Stderr, "Sorting by gene id...\n")

	inname := path.Join(config.TempDir, "matches.txt.sz")
	outname := path.Join(config.TempDir, "matches_sg.txt.sz")

	pr1, pw1, err := os.Pipe()
	if err != nil {
		panic(err)
	}

	pr2, pw2, err := os.Pipe()
	if err != nil {
		panic(err)
	}

	// Sort by gene number
	cmd1 := exec.Command("sztool", "-d", inname)
	cmd1.Stdout = pw1

	// k5 is position of gene id
	args := []string{sortmem, sortpar, "-k5"}
	if sortTmpFlag != "" {
		args = append(args, sortTmpFlag)
	}
	args = append(args, "-")
	cmd2 := exec.Command("sort", args...)
	cmd2.Stdin = pr1
	cmd2.Stdout = pw2

	// Compress the results
	cmd3 := exec.Command("sztool", "-c", "-", outname)
	cmd3.Stdin = pr2

	for _, cmd := range []*exec.Cmd{cmd1, cmd2, cmd3} {
		cmd.Stderr = os.Stderr
		if err := cmd.Start(); err != nil {
			panic(err)
		}
	}

	if err := cmd1.Wait(); err != nil {
		panic(err)
	}

	pw1.Close()

	if err := cmd2.Wait(); err != nil {
		panic(err)
	}

	pw2.Close()

	if err := cmd3.Wait(); err != nil {
		panic(err)
	}
}

func joinGeneNames() {

	io.WriteString(os.Stderr, "Joining gene names...\n")

	pr1, pw1, err := os.Pipe()
	if err != nil {
		panic(err)
	}

	pr2, pw2, err := os.Pipe()
	if err != nil {
		panic(err)
	}

	// Join genes and matches
	fn := path.Join(config.TempDir, "matches_sg.txt.sz")
	bs := fmt.Sprintf("join -1 5 -2 1 -t $'\t' <(sztool -d %s) <(sztool -d %s)\n", fn, config.GeneIdFileName)
	fid, err := os.Create("bs.sh")
	io.WriteString(fid, bs)
	fid.Close()
	cmd1 := exec.Command("/bin/bash", "bs.sh")
	cmd1.Stdout = pw1

	// Cut out unwanted column
	// The first argument after cur is -d(tab)
	cmd2 := exec.Command("cut", "-d	", "-f1", "--complement", "-")
	cmd2.Stdin = pr1
	cmd2.Stdout = pw2

	// Compress the result
	cmd3 := exec.Command("sztool", "-c", "-", path.Join(config.TempDir, "matches_sn.txt.sz"))
	cmd3.Stdin = pr2

	for _, cmd := range []*exec.Cmd{cmd1, cmd2, cmd3} {
		cmd.Stderr = os.Stderr
		if err := cmd.Start(); err != nil {
			panic(err)
		}
	}

	if err := cmd1.Wait(); err != nil {
		panic(err)
	}

	pw1.Close()

	if err := cmd2.Wait(); err != nil {
		panic(err)
	}

	pw2.Close()

	if err := cmd3.Wait(); err != nil {
		panic(err)
	}
}

func joinReadNames() {

	io.WriteString(os.Stderr, "Joining read names...\n")

	fn := path.Join(config.TempDir, "reads_sorted.txt.sz")
	gn := path.Join(config.TempDir, "matches_sn.txt.sz")

	if _, err := os.Stat(fn); os.IsNotExist(err) {
		err := fmt.Errorf("reads_sorted.txt.sz does not exist")
		panic(err)
	}

	if _, err := os.Stat(gn); os.IsNotExist(err) {
		err := fmt.Errorf("matches_sn.txt.sz does not exist")
		panic(err)
	}

	c1 := fmt.Sprintf("<(sort -k1 %s %s %s <(sztool -d %s))", sortmem, sortpar, sortTmpFlag, gn)
	c2 := fmt.Sprintf("<(sztool -d %s)", fn)
	bs := fmt.Sprintf("join -1 1 -2 1 -t'\t' %s %s > %s", c1, c2, config.ResultsFileName)
	fid, err := os.Create("bs.sh")
	if err != nil {
		panic(err)
	}
	_, err = io.WriteString(fid, bs)
	if err != nil {
		panic(err)
	}
	fid.Close()

	cmd := exec.Command("/bin/bash", "bs.sh")
	if err := cmd.Run(); err != nil {
		panic(err)
	}
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

func setupLog() {
	logname := path.Join(config.LogDir, "muscato.log")
	fid, err := os.Create(logname)
	if err != nil {
		panic(err)
	}
	logger = log.New(fid, "", log.Ltime)
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
		os.MkdirAll(config.SortTemp, os.ModePerm)
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
	err = os.MkdirAll(config.TempDir, os.ModePerm)
	if err != nil {
		if os.IsNotExist(err) {
			msg := fmt.Sprintf("Directory %s does not exist and cannot be created.", config.TempDir)
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

	err = os.MkdirAll(config.LogDir, os.ModePerm)
	if err != nil {
		panic(err)
	}
}

func cleanTmp() {

	if config.NoCleanTemp {
		return
	}

	err := os.RemoveAll(config.TempDir)
	if err != nil {
		panic(err)
	}
}

func genReadStats() {

	io.WriteString(os.Stderr, "Generating read statistics...\n")

	cmd := exec.Command("muscato_readstats", configFilePath)
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		panic(err)
	}
}

func writeNonMatch() {

	io.WriteString(os.Stderr, "Writing non-matching sequences...\n")

	cmd := exec.Command("muscato_nonmatch", configFilePath)
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		panic(err)
	}
}

func main() {

	defer cleanTmp()

	handleArgs()
	checkArgs()
	setupEnvs()
	makeTemp()

	// The logger is not available until after makeTemp runs.
	setupLog()

	logger.Printf("Starting saveConfig...\n")
	saveConfig(config)

	logger.Printf("Starting prepReads...\n")
	prepReads()

	logger.Printf("Starting windowReads...\n")
	windowReads()

	logger.Printf("Starting sortWindows...\n")
	sortWindows()

	logger.Printf("Starting screen...\n")
	screen()

	logger.Printf("Starting sortBloom...\n")
	sortBloom()

	logger.Printf("Starting confirm...\n")
	confirm()

	logger.Printf("Starting combineWindows...\n")
	combineWindows()

	logger.Printf("Starting sortByGeneId...\n")
	sortByGeneId()

	logger.Printf("Starting joinGeneNames...\n")
	joinGeneNames()

	logger.Printf("Starting joinReadNames...\n")
	joinReadNames()

	logger.Printf("Starting writeNoneMatch...\n")
	writeNonMatch()

	logger.Printf("Starting genReadStats...\n")
	genReadStats()

	logger.Printf("Starting geneStats...\n")
	geneStats()
}
