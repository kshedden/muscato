```
Usage of muscato:
  -BloomSize int
    	Size of Bloom filter, in bits
  -ConfigFileName string
    	JSON file containing configuration parameters
  -GeneFileName string
    	Gene file name (processed form)
  -GeneIdFileName string
    	Gene ID file name (processed form)
  -MMTol int
    	Number of mismatches allowed above best fit
  -MatchMode string
    	'first' or 'best' (retain first/best 'MaxMatches' matches meeting criteria)
  -MaxConfirmProcs int
    	Run this number of match confirmation processes concurrently
  -MaxMatches int
    	Return no more than this number of matches per window
  -MaxReadLength int
    	Reads longer than this length are truncated
  -MinDinuc int
    	Minimum number of dinucleotides to check for match
  -MinReadLength int
    	Reads shorter than this length are skipped
  -NoCleanTemp
    	Do not delete temporary files from TempDir
  -NumHash int
    	Number of hashses
  -PMatch float
    	Required proportion of matching positions
  -ReadFileName string
    	Sequencing read file (fastq format)
  -ResultsFileName string
    	File name for results
  -SortPar int
    	Number of parallel sort processes (default 8)
  -SortTemp string
    	Directory to use for sort temp files
  -TempDir string
    	Workspace for temporary files
  -WindowWidth int
    	Width of each window
  -Windows string
    	Starting position of each window
```