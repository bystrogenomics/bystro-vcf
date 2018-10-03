package main

import (
	"bufio"
	"os"
	"testing"
)

func benchmarkReader(b *testing.B) {
	config := &Config{}

	outFh, err := os.OpenFile("/dev/null", os.O_WRONLY|os.O_CREATE, 0644)
	if err != nil {
		panic(err)
	}
	writer := bufio.NewWriter(outFh)
	// make sure it gets closed
	defer outFh.Close()

	for n := 0; n < b.N; n++ {
		inFh, err := os.Open("examples/test.query.vcf")
		if err != nil {
			panic(err)
		}

		reader := bufio.NewReader(inFh)

		readVcf(config, reader, writer)
	}

	writer.Flush()
	outFh.Close()
}

func BenchmarkLazy(b *testing.B) { benchmarkReader(b) }
