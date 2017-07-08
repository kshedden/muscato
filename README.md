__Muscato (Multi-Genome Scalable Alignment Tool)__

Muscato is a software tool for matching a collection of sequence reads
into a collection of target sequences (e.g. gene sequences).  The
approach effectively scales to hundreds of millions of reads and
target sequences.  A major goal of Muscato is to perform exhaustive
multi-mapping, meaning that each read is mapped to as many gene
sequences as possible, subject to specified match quality constraints.

__Installation__

Muscato is written in [Go](https://golang.org) and uses several of the
[Gnu core
utilities](http://www.gnu.org/software/coreutils/coreutils.html).  It
should run on any Unix-like system on which the [Go
tool](https://golang.org/dl) and Gnu utilities are available.

In most cases installation of Muscato should only require invoking the
following commands in the shell:

```
go get github.com/kshedden/muscato/...

go get github.com/kshedden/sztool/...
```

The executables for muscato and its auxiliary scripts should appear in
your GOBIN directory (usually ${HOME}/go/bin if installed in a user
account).  You will need to add GOBIN to your PATH environment
variable when using the tool.

__Basic usage__

Before running Muscato, you should prepare a version of your target
sequence file using the `muscato_prep_targets` program.  If your
targets are in a fasta format file, you can simply run:

```
muscato_prep_targets genes.fasta
```

Instead of using a fasta input file, it is also possible to use a
plain text file with the format `id<tab>sequence<newline>` for each
target sequence.  The sequence should consist of the upper-case
characters A, T, G, and C.  Any other letters are replaced with 'X'.

The `muscato_prep_targets` script accepts a `-rev` flag in which
reverse complement target sequences are added to the database along
with the original sequences.

Now you can run muscato.  A basic invocation is:

```
muscato --ReadFileName=reads.fastq --GeneFileName=genes.fasta.sz --GeneIdFileName=genes_ids.sz\
        --Windows=0,20 --WindowWidth=15
```

Note that the target files `genes.fasta.sz` and `genes_ids.sz` were
produced by the `muscato_prep_targets` script, to be run as shown
above.

Many other command-line flags are available, run `muscato --help` for
more information.

The results by default are written to a file named `results.txt`, a
tab delimited file with the following columns:

1. Read sequence

2. Matching subsequence of a target sequence

3. Position within the target where the read matches (counting from 0)

4. Number of mismatches

5. Target sequence identifier

6. Target sequence length

7. Number of copies of the read in the read pool

8. Read identifier

The tool also generates a fastq file containing all non-matching reads.

Muscato uses a temporary directory for intermediate and logging files,
by default named `tmp/######`, where ###### is a integer id generated
by the tool.  This directory can be safely deleted after the results
are obtained.

__Logging__

Several log files are written to the temporary directory discussed
above.  High-level logging messages are written to 'muscato.log'.
More detailed logging information written to logs specific to each
component of the tool, e.g. 'muscato_screen.log'.

__Testing__

There is currently a small collection of unit tests in the `tests`
directory.  To run the tests, enter the test directory and type:

```
go run test.go
```

Any errors will be printed to the terminal.  Detailed results of the
tests are written to the file `test.log`.

__Dependencies__

Muscato has the following dependencies, which should be automatically
installed by `go get` when installing muscato, as described above
(sztool must be installed directly).

[github.com/kshedden/sztool](http://github.com/kshedden/sztool)

[github.com/chmduquesne/rollinghash](http://github.com/chmduquesne/rollinghash)

[github.com/golang-collections/go-datastructures/bitarray](http://github.com/golang-collections/go-datastructures/bitarray)

[github.com/golang/snappy](http://github.com/golang/snappy)

[github.com/willf/bloom](http://github.com/willf/bloom)

[github.com/scipipe/scipipe](http://github.com/scipipe/scipipe)