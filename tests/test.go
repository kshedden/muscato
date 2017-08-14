// Copyright 2017, Kerby Shedden and the Muscato contributors.

// test is a script that runs a series of unit tests on the Muscato
// code base.
//
// To run the tests, use:
//
// go run test.go

package main

import (
	"bufio"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path"
	"strings"

	"github.com/BurntSushi/toml"
	"github.com/golang/snappy"
)

var (
	logger *log.Logger
)

type Test struct {
	Name    string
	Base    string
	Command string
	Opts    []string
	Args    []string
	Files   [][2]string
}

func getTests() []Test {

	fid, err := os.Open("tests.toml")
	if err != nil {
		panic(err)
	}
	s, err := ioutil.ReadAll(fid)
	if err != nil {
		panic(err)
	}
	fid.Close()

	type vd struct {
		Test []Test
	}

	var v vd
	_, err = toml.Decode(string(s), &v)
	if err != nil {
		panic(err)
	}

	logger.Printf("Found %d tests\n", len(v.Test))

	return v.Test
}

// getScanner returns a scanner for reading the contents of a file.
// Snappy compression is handled automatically.  An array of values
// that should be closed when the scanner is no longer needed is also
// returned.
func getScanner(f string) (*bufio.Scanner, []io.Closer) {

	var toclose []io.Closer
	var g io.Reader

	h, err := os.Open(f)
	if err != nil {
		panic(err)
	}
	toclose = append(toclose, h)
	g = h

	if strings.HasSuffix(f, ".sz") {
		g = snappy.NewReader(g)
	}

	s := bufio.NewScanner(g)
	return s, toclose
}

// compare returns true if and only if the contents of the files named
// by the arguments f1 and f2 are identical.  Snappy compression is
// handled automatically.
func compare(f1, f2 string) bool {

	s1, tc1 := getScanner(f1)
	s2, tc2 := getScanner(f2)

	for {
		q1 := s1.Scan()
		q2 := s2.Scan()

		if q1 != q2 {
			msg := fmt.Sprintf("files %s and %s have different numbers of lines\n", f1, f2)
			panic(msg)
		}
		if !q1 {
			break
		}

		v1 := s1.Text()
		v2 := s2.Text()
		if v1 != v2 {
			msg := fmt.Sprintf("%s\nin file %s\ndiffers from\n%v\nin file %s\n", v1, f1, v2, f2)
			panic(msg)
		}
	}

	if err := s1.Err(); err != nil {
		panic(err)
	}
	if err := s2.Err(); err != nil {
		panic(err)
	}

	for _, x := range tc1 {
		x.Close()
	}
	for _, x := range tc2 {
		x.Close()
	}

	return true
}

func run(tests []Test) {

	for _, t := range tests {

		c := []string{path.Join(t.Command)}
		for _, o := range t.Opts {
			c = append(c, o)
		}
		for _, f := range t.Args {
			c = append(c, path.Join(t.Base, f))
		}
		logger.Printf("%s\n", t.Name)
		logger.Printf("Running command %s\n", c[0])
		logger.Printf("with arguments: %v\n", c[1:])
		cmd := exec.Command(c[0], c[1:len(c)]...)
		cmd.Stderr = os.Stderr
		err := cmd.Run()
		if err != nil {
			panic(err)
		}
		for _, fp := range t.Files {
			compare(path.Join(t.Base, fp[0]), path.Join(t.Base, fp[1]))
		}

		logger.Printf("done\n\n")
	}
}

func setupLog() {
	fid, err := os.Create("test.log")
	if err != nil {
		panic(err)
	}
	logger = log.New(fid, "", log.Ltime)
}

func main() {
	setupLog()
	tests := getTests()
	run(tests)
}
