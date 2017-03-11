// estab exports elasticsearch fields as tab separated values
package main

import (
	"bufio"
	"bytes"
	// "encoding/csv"
	// "encoding/json"
	"flag"
	"fmt"
	"io"
	// "io/ioutil"
	"log"
	// "math"
	"os"
	// "sort"
	"strconv"
	"strings"

	// "github.com/davecgh/go-spew/spew"
	// "math/big"
	"runtime/pprof"
	// "sync"
	// "errors"
)

// NOTE: For now this only supports \n end of line characters
// If we want to support CLRF or whatever, use either csv package, or set a different delimiter
func main() {
	inputFilePath := flag.String("inPath", "", "The input file path (optional: default is stdin)")
	errPath := flag.String("errPath", "", "The output path for the JSON output (optional)")
	emptyField := flag.String("emptyField", "!", "The output path for the JSON output (optional)")
	// chrPrefix := flag.Bool("ucscChr", "", "Whether or not to use UCSC style chromosome designations, i.e chrX")

	cpuprofile := flag.String("cpuProfile", "", "write cpu profile to file")
	flag.Parse()

	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	inFh := (*os.File)(nil)

	if *inputFilePath != "" {
		var err error
		inFh, err = os.Open(*inputFilePath)

		if err != nil {
			log.Fatal(err)
		}
	} else {
		inFh = os.Stdin
	}

	// make sure it gets closed
	defer inFh.Close()

	if *errPath != "" {
		var err error
		os.Stderr, err = os.Open(*errPath)

		if err != nil {
			log.Fatal(err)
		}
	}

	reader := bufio.NewReader(inFh)

	foundHeader := false
	// checkedChrType := false

	//Predeclar sampleNames to be a large item
	header := make([]string, 0, 10000)

	var lastIndex int

	c := make(chan string)

	for {
		// http://stackoverflow.com/questions/8757389/reading-file-line-by-line-in-go
		// http://www.jeffduckett.com/blog/551119d6c6b86364cef12da7/golang---read-a-file-line-by-line.html
		// Scanner doesn't work well, has buffer restrictions that we need to manually get around
		// and we don't expect any newline characters in a Seqant output body
		row, err := reader.ReadString('\n') // 0x0A separator = newline

		if err == io.EOF {
			// do something here
			break
		} else if err != nil {
			log.Fatal(err)
		}

		// // remove the trailing \n
		// // equivalent of chomp https://groups.google.com/forum/#!topic/golang-nuts/smFU8TytFr4
		record := strings.Split(row[:len(row)-1], "\t")

		if foundHeader == false {
			if record[0] == "#CHROM" {

				lastIndex = len(record) - 1

				if lastIndex < 9 {
					log.Fatal("Expected to find at least 1 sample, 0 found")
				}

				header = record

				// spew.Dump(header)

				// homs = make([]string, 0, len(header)-9)
				// homs = make([]string, 0, len(header)-9)

				foundHeader = true
				break
			}
		}
	}

	// spew.Dump(header)

	go func() {
		for {
			row, err := reader.ReadString('\n') // 0x0A separator = newline

			if err == io.EOF {
				// do something here
				close(c)
				break
			} else if err != nil {
				log.Fatal(err)
			}

			// // remove the trailing \n
			// // equivalent of chomp https://groups.google.com/forum/#!topic/golang-nuts/smFU8TytFr4
			record := strings.Split(row[:len(row)-1], "\t")

			// spew.Dump(record)
			if record[4] == "." || record[3] == record[4] || (record[3][0:1] != "A" && record[3][0:1] != "C" && record[3][0:1] != "T" && record[3][0:1] != "G") {
				continue
			}

			// homs = homs[:0]
			// hets = homs[:0]

			if strings.Contains(record[4], ",") {
				go processMultiLine(record, header, lastIndex, emptyField, c)
			} else {
				go processLine(record, header, lastIndex, emptyField, c)
			}
		}
	}()

	// Write all the data
	// Somewhat surprisingly this is faster than building up array and writing in builk
	for data := range c {
		fmt.Print(data)
	}

}

func processMultiLine(record []string, header []string, lastIndex int, emptyField *string, results chan<- string) {
	var output bytes.Buffer

	if record[0] != "c" && record[1] != "h" {
		output.WriteString("chr")
		output.WriteString(record[0])
		output.WriteString("\t")
	}
	// altAlleles := strings.Split(record[4], ",")

	// spew.Dump(altAlleles)

	// records := make([]string, 0, len(altAlleles))

	for idx, allele := range strings.Split(record[4], ",") {
		if allele[0:1] != "A" && allele[0:1] != "C" && allele[0:1] != "T" && allele[0:1] != "G" {
			continue
		}

		siteType, pos, ref, alt, err := updateFieldsWithAlt(record[3], allele, record[1])

		if err != nil {
			log.Fatal(err)
		}

		if pos == "" {
			continue
		}

		// Attempt at reducing malloc
		var homs []string
		var hets []string

		// fmt.Printf("%d %s", idx+1, strconv.Itoa(idx+1))
		err = makeHetHomozygotes(record, header, lastIndex, &homs, &hets, strconv.Itoa(idx+1))

		if err != nil {
			log.Fatal(err)
		}

		if len(homs) == 0 && len(hets) == 0 {
			continue
		}

		output.WriteString(pos)
		output.WriteString("\t")
		output.WriteString(siteType)
		output.WriteString("\t")
		output.WriteString(ref)
		output.WriteString("\t")
		output.WriteString(alt)

		output.WriteString("\t")

		if len(hets) == 0 {
			output.WriteString(*emptyField)
		} else {
			output.WriteString(strings.Join(hets, ","))
		}

		output.WriteString("\t")

		if len(homs) == 0 {
			output.WriteString(*emptyField)
		} else {
			output.WriteString(strings.Join(homs, ","))
		}

		output.WriteString("\n")

		results <- output.String()
	}
}

func processLine(record []string, header []string, lastIndex int, emptyField *string, results chan<- string) {
	var homs []string
	var hets []string

	var output bytes.Buffer

	if record[0] != "c" && record[1] != "h" {
		output.WriteString("chr")
		output.WriteString(record[0])
		output.WriteString("\t")
	}

	siteType, pos, ref, alt, err := updateFieldsWithAlt(record[3], record[4], record[1])

	if err != nil {
		log.Fatal(err)
	}

	if pos == "" {
		fmt.Fprint(os.Stderr, "NO POSS!!!", record)
	}

	// spew.Dump(hets)
	err = makeHetHomozygotes(record, header, lastIndex, &homs, &hets, "1")

	if err != nil {
		log.Fatal(err)
	}

	if len(homs) == 0 && len(hets) == 0 {
		return
	}

	output.WriteString(pos)
	output.WriteString("\t")
	output.WriteString(siteType)
	output.WriteString("\t")
	output.WriteString(ref)
	output.WriteString("\t")
	output.WriteString(alt)

	output.WriteString("\t")

	if len(hets) == 0 {
		output.WriteString(*emptyField)
	} else {
		output.WriteString(strings.Join(hets, ","))
	}

	output.WriteString("\t")

	if len(homs) == 0 {
		output.WriteString(*emptyField)
	} else {
		output.WriteString(strings.Join(homs, ","))
	}

	output.WriteString("\n")

	results <- output.String()
}

func updateFieldsWithAlt(ref string, alt string, pos string) (string, string, string, string, error) {
	var siteType string

	if len(alt) == len(ref) {
		if len(ref) > 1 {
			// We don't support MNPs
			return "", "", "", "", nil
		}

		return "SNP", pos, ref, alt, nil
	}

	if len(ref) > len(alt) {
		siteType = "DEL"

		alt = strconv.Itoa(len(alt) - len(ref))
		intPos, err := strconv.Atoi(pos)

		if err != nil {
			return "", "", "", "", err
		}

		pos = strconv.Itoa(intPos + len(alt))
		ref = ref[1:2]
	} else {
		siteType = "INS"
		var insBuffer bytes.Buffer
		insBuffer.WriteString("+")
		insBuffer.WriteString(alt[len(ref):])

		alt = insBuffer.String()

		// Most cases are simple, handle complex ones as well
		if len(ref) > 1 {
			insIndex := strings.Index(alt, ref)

			if insIndex == -1 {
				return "", "", "", "", nil
			}

			if insIndex != 0 {
				log.Println("Don't support complex insertion alleles")
			}

			intPos, err := strconv.Atoi(pos)
			if err != nil {
				log.Fatal("Failed to convert position to integer")
			}

			intPos = intPos + len(ref) - 1
			pos = strconv.Itoa(intPos)
			ref = ref[len(ref)-1:]
		}
	}

	return siteType, pos, ref, alt, nil
}

func makeHetHomozygotes(fields []string, header []string, lastIdx int, homsArr *[]string, hetsArr *[]string, alleleIdx string) error {
	simpleGT := fields[8] == "GT"

	gt := make([]string, 0, 2)
	gtCount := 0
	altCount := 0

	for i := 9; i <= lastIdx; i++ {
		if strings.Contains(fields[i], "|") {
			gt = strings.Split(fields[i], "|")
		} else {
			gt = strings.Split(fields[i], "/")
		}

		altCount = 0
		gtCount = 0
		if simpleGT {
			for _, val := range gt {
				if val == alleleIdx {
					altCount++
				}

				gtCount++
			}
		} else {

			for _, val := range gt {
				// val[0] doesn't work...
				// fmt.Printf("Genotype %s, %s", gt, val[0:1])

				if val[0:1] == alleleIdx {
					altCount++
				}

				gtCount++
			}
		}

		if altCount == 0 {
			continue
		}

		if altCount == gtCount {
			*homsArr = append(*homsArr, header[i])
		} else {
			*hetsArr = append(*hetsArr, header[i])
		}
	}

	// spew.Dump(hetsArr)
	// spew.Dump(homsArr)

	return nil
}

// func coercePositionsAlleles(emptyField string, primaryDelim string) func(string) bool {
