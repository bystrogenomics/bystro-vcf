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
	fieldDelimiter := flag.String("fieldDelimiter", ";", "The output path for the JSON output (optional)")
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
	var header []string

	var lastIndex int

	c := make(chan string)

	for {
		// http://stackoverflow.com/questions/8757389/reading-file-line-by-line-in-go
		// http://www.jeffduckett.com/blog/551119d6c6b86364cef12da7/golang---read-a-file-line-by-line.html
		// Scanner doesn't work well, has buffer restrictions that we need to manually get around
		// and we don't expect any newline characters in a Seqant output body
		row, err := reader.ReadString('\n') // 0x0A separator = newline

		if err == io.EOF {
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

				foundHeader = true
				break
			}
		}
	}

	// TODO: will we be more efficient if we pre-make these and clear them each round?
	// homSingle := make([]string, 0, lastIndx-8)
	// hetsSingle := make([]string, 0, lastIndx-8)
	// homMulti := make([]string, 0, lastIndx-8)
	// hetsMulti := make([]string, 0, lastIndx-8)

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

			if len(record) < len(header) {
				log.Println("Truncated line. Skipping: ", record[0], record[1])
				continue
			}

			if strings.Contains(record[4], ",") {
				go processMultiLine(record, header, lastIndex, emptyField, fieldDelimiter, c)
			} else {
				go processLine(record, header, lastIndex, emptyField, fieldDelimiter, c)
			}
		}
	}()

	// Write all the data
	// Somewhat surprisingly this is faster than building up array and writing in builk
	for data := range c {
		fmt.Print(data)
	}
}

func lineIsValid(alt string) bool {
	if len(alt) == 1 {
		if alt != "A" && alt != "C" && alt != "T" && alt != "G" {
			return false
		}
	} else {
		// Most common case is <CNV|DUP|DEL>
		if alt[0] != 'A' && alt[0] != 'C' && alt[0] != 'T' && alt[0] != 'G' {
			return false
		}

		for _, val := range alt[1:] {
			if val != 'A' && val != 'C' && val != 'T' && val != 'G' {
				return false
			}
		}
	}

	return true
}

func processMultiLine(record []string, header []string, lastIndex int, emptyField *string,
	fieldDelimiter *string, results chan<- string) {
	var chr string

	if len(record[0]) < 4 || record[0][0:2] != "ch" {
		var buff bytes.Buffer
		buff.WriteString("chr")
		buff.WriteString(record[0])

		chr = buff.String()
	}

	for idx, allele := range strings.Split(record[4], ",") {
		if lineIsValid(allele) == false {
			log.Printf("Non-ACTG Alt #%d, skipping: %s %s", idx+1, record[0], record[1])
			continue
		}

		siteType, pos, ref, alt, err := updateFieldsWithAlt(record[3], allele, record[1])

		if err != nil {
			log.Fatal(err)
		}

		if pos == "" {
			log.Printf("Invalid Alt #%d, skipping: %s %s", idx+1, record[0], record[1])
			continue
		}

		// Attempt at reducing malloc
		var homs []string
		var hets []string

		err = makeHetHomozygotes(record, header, &homs, &hets, strconv.Itoa(idx+1))

		if err != nil {
			log.Fatal(err)
		}

		if len(homs) == 0 && len(hets) == 0 {
			continue
		}

		var output bytes.Buffer

		if chr != "" {
			output.WriteString(chr)
		} else {
			output.WriteString(record[0])
		}

		output.WriteString("\t")
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
			output.WriteString(strings.Join(hets, *fieldDelimiter))
		}

		output.WriteString("\t")

		if len(homs) == 0 {
			output.WriteString(*emptyField)
		} else {
			output.WriteString(strings.Join(homs, *fieldDelimiter))
		}

		output.WriteString("\n")

		results <- output.String()
	}
}

func processLine(record []string, header []string, lastIndex int,
	emptyField *string, fieldDelimiter *string, results chan<- string) {
	if lineIsValid(record[4]) == false {
		log.Println("Non-ACTG Alt, skipping: ", record[0], record[1])
		return
	}

	var homs []string
	var hets []string

	var output bytes.Buffer

	if len(record[0]) < 4 || record[0][0:2] != "ch" {
		output.WriteString("chr")
		output.WriteString(record[0])
		output.WriteString("\t")
	}

	siteType, pos, ref, alt, err := updateFieldsWithAlt(record[3], record[4], record[1])

	if err != nil {
		log.Fatal(err)
	}

	if pos == "" {
		log.Println("Invalid Alt, skipping: ", record[0], record[1])
		return
	}

	err = makeHetHomozygotes(record, header, &homs, &hets, "1")

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
		output.WriteString(strings.Join(hets, *fieldDelimiter))
	}

	output.WriteString("\t")

	if len(homs) == 0 {
		output.WriteString(*emptyField)
	} else {
		output.WriteString(strings.Join(homs, *fieldDelimiter))
	}

	output.WriteString("\n")

	results <- output.String()
}

func updateFieldsWithAlt(ref string, alt string, pos string) (string, string, string, string, error) {
	var siteType string

	if len(alt) == len(ref) {
		if alt == ref {
			// No point in returning ref sites
			return "", "", "", "", nil
		}

		if len(ref) > 1 {
			count := 0
			var diffIdx int
			for index, _ := range ref {
				if ref[index] != alt[index] {
					if count != 0 {
						// MNPs are not supported
						return "", "", "", "", nil
					}

					count += 1
					diffIdx = index
				}
			}

			intPos, err := strconv.Atoi(pos)

			if err != nil {
				return "", "", "", "", err
			}

			return "SNP", strconv.Itoa(intPos + diffIdx), ref[diffIdx : diffIdx+1], alt[diffIdx : diffIdx+1], nil
		}

		return "SNP", pos, ref, alt, nil
	}

	if len(ref) > len(alt) {
		siteType = "DEL"

		intPos, err := strconv.Atoi(pos)

		if err != nil {
			return "", "", "", "", err
		}

		if len(alt) == 1 {
			if ref[0] != alt[0] {
				return "", "", "", "", nil
			}

			alt = strconv.Itoa(1 - len(ref))
			pos = strconv.Itoa(intPos + 1)
			ref = ref[1:2]
		} else {
			altLen := len(alt)

			if ref[0:altLen] != alt {
				return "", "", "", "", nil
			}

			pos = strconv.Itoa(intPos + altLen)
			alt = strconv.Itoa(altLen - len(ref))
			ref = ref[altLen : altLen+1]
		}
	} else {
		siteType = "INS"

		// Most cases are simple, handle complex ones as well
		if len(ref) > 1 {
			insIndex := strings.Index(alt, ref)

			if insIndex != 0 {
				return "", "", "", "", nil
			}

			intPos, err := strconv.Atoi(pos)

			if err != nil {
				log.Fatal("Failed to convert position to integer")
			}

			intPos = intPos + len(ref) - 1
			pos = strconv.Itoa(intPos)

			var insBuffer bytes.Buffer
			insBuffer.WriteString("+")
			insBuffer.WriteString(alt[len(ref):])

			alt = insBuffer.String()

			ref = ref[len(ref)-1:]
		} else {
			if ref[0] != alt[0] {
				return "", "", "", "", nil
			}

			var insBuffer bytes.Buffer
			insBuffer.WriteString("+")
			insBuffer.WriteString(alt[1:])

			alt = insBuffer.String()
		}
	}

	return siteType, pos, ref, alt, nil
}

func makeHetHomozygotes(fields []string, header []string, homsArr *[]string, hetsArr *[]string, alleleIdx string) error {
	simpleGT := fields[8] == "GT"

	gt := make([]string, 0, 2)
	gtCount := 0
	altCount := 0

SAMPLES:
	for i := 9; i < len(header); i++ {
		if strings.Contains(fields[i], "|") {
			gt = strings.Split(fields[i], "|")
		} else {
			gt = strings.Split(fields[i], "/")
		}

		altCount = 0
		gtCount = 0
		if simpleGT {
			for _, val := range gt {
				if val == "." {
					continue SAMPLES
				}

				if val == alleleIdx {
					altCount++
				}

				gtCount++
			}
		} else {

			for _, val := range gt {
				if val == "." {
					continue SAMPLES
				}

				// val[0] doesn't work...because byte array?
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

	return nil
}

// func coercePositionsAlleles(emptyField string, primaryDelim string) func(string) bool {
