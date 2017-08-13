// Copyright 2017, Kerby Shedden and the Muscato contributors.

package utils

import (
	"encoding/json"
	"os"
)

type Config struct {

	// The name of the fastq file containing the reads.
	ReadFileName string

	// The name of the fasta or plain text file containing the
	// target sequences (genes).
	GeneFileName string

	// The name of the file containing the target sequence (gene)
	// identifiers.
	GeneIdFileName string

	// The file path where the results are written.
	ResultsFileName string

	// The left end point of each window with a read.
	Windows []int

	// The width of each window.
	WindowWidth int

	// The size of the Bloom filter in bits.
	BloomSize uint64

	// The number of hash functions to use in the Bloom filter.
	NumHash int

	// The minimum allowed proportion of matching bases.
	PMatch float64

	// The exact-match subsequence must have this many distinct
	// dinucleotide subsequences.
	MinDinuc int

	// Use this location to place temporary files.  If blank or
	// missing, a temporary directory is generated of the form
	// tmp/######## in the local directory.
	TempDir string

	// The directory where log files are written.  By default the
	// logs are placed into muscato_logs/###### in the local
	// directory, where the number matches the default prefix of
	// the temporary directory.
	LogDir string

	// Skip all reads shorter than this length.
	MinReadLength int

	// Truncate all reads at this length.
	MaxReadLength int

	// The confirmatory matching step returns at most this many
	// matches for each k-mer seqeunces.  Since a k-mer sequence
	// may match many reads and many genes, setting MaxMatches to
	// a low value may lead to some reads not being mapped, or not
	// multi-mapping as well as possible.
	MaxMatches int

	// The maximum number of merge processes that are run
	// simultaneously.
	MaxMergeProcs int

	// Number of additional mismatches beyond the best possible
	// number of mismatches that are allowed when retaining the
	// target sequence matches to each read.
	MMTol int

	// Either "first" (default) or "best".  If first, returns the
	// first MaxMatches matches for each window.  If best, returns
	// the MaxMatches matches for each window with the fewest
	// mismatched values.
	MatchMode string

	// If true, temporary files are not removed upon program
	// completion.  If false, which is the default, the temporary
	// files are removed.
	NoCleanTmp bool
}

func ReadConfig(filename string) *Config {
	fid, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	dec := json.NewDecoder(fid)
	config := new(Config)
	err = dec.Decode(config)
	if err != nil {
		panic(err)
	}

	return config
}
