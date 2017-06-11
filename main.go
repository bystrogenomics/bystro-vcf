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
	"sync"

	// "github.com/davecgh/go-spew/spew"
	// "math/big"
	"runtime/pprof"
	// "sync"
	// "errors"
)

const chromIdx int = 0
const posIdx int = 1
const idIdx int = 2
const refIdx int = 3
const altIdx int = 4
const qualIdx int = 5
const filterIdx int = 6
const infoIdx int = 7
const formatIdx int = 8

// NOTE: For now this only supports \n end of line characters
// If we want to support CLRF or whatever, use either csv package, or set a different delimiter
func main() {
	inputFilePath := flag.String("inPath", "", "The input file path (optional: default is stdin)")
	errPath := flag.String("errPath", "", "The output path for the JSON output (optional)")
	emptyFieldOpt := flag.String("emptyField", "!", "The output path for the JSON output (optional)")
	fieldDelimiterOpt := flag.String("fieldDelimiter", ";", "The output path for the JSON output (optional)")
	retainIdOpt := flag.Bool("retainID", false, "Retain the ID field in the output (1 additional output field, before the reatinInfo output fields, should those be present")
	retainInfoOpt := flag.Bool("retainInfo", false, "Should we retain INFO field data (if so, will output 2 additional fields, the index of the allele (to handle multiallelic segregation of INFO data properly), and the INFO field")
	// chrPrefix := flag.Bool("ucscChr", "", "Whether or not to use UCSC style chromosome designations, i.e chrX")

	cpuprofile := flag.String("cpuProfile", "", "write cpu profile to file")
	flag.Parse()

	log.SetFlags(0)
	//Dereference variables that will be passed into loop, supposedly cheaper to pass by value
	//http://goinbigdata.com/golang-pass-by-pointer-vs-pass-by-value/
	//https://stackoverflow.com/questions/24452323/go-performance-whats-the-difference-between-pointer-and-value-in-struct
	emptyField := *emptyFieldOpt
	fieldDelimiter := *fieldDelimiterOpt
	retainInfo := *retainInfoOpt
	retainId := *retainIdOpt

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

	c := make(chan string)
	// I think we need a wait group, not sure.
	wg := new(sync.WaitGroup)

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
			if record[chromIdx] == "#CHROM" {
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
	normalizeSampleNames(header)

	go func() {
		for {
			row, err := reader.ReadString('\n') // 0x0A separator = newline

			if err == io.EOF {
				break
			} else if err != nil {
				log.Fatal(err)
			}

			// // remove the trailing \n
			// // equivalent of chomp https://groups.google.com/forum/#!topic/golang-nuts/smFU8TytFr4
			record := strings.Split(row[:len(row)-1], "\t")

			if linePasses(record, header) == false {
				continue
			}

			wg.Add(1)
			go processLine(record, header, emptyField, fieldDelimiter, retainId, retainInfo, c, wg)
		}

		wg.Wait()
		close(c)
	}()

	// Write all the data
	// Somewhat surprisingly this is faster than building up array and writing in builk
	for data := range c {
		fmt.Print(data)
	}
}

func chrToUCSC(chr string) string {
	if len(chr) < 4 || chr[0:2] != "ch" {
		var buff bytes.Buffer
		buff.WriteString("chr")
		buff.WriteString(chr)

		return buff.String()
	}

	return chr
}

func linePasses(record []string, header []string) bool {
	return len(record) == len(header) && (record[filterIdx] == "." || record[filterIdx] == "PASS")
}

func altIsValid(alt string) bool {
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

// Without this (calling each separately) real real	1m0.753s, with:	1m0.478s
func processLine(record []string, header []string, emptyField string,
	fieldDelimiter string, retainId bool, retainInfo bool, c chan<- string, wg *sync.WaitGroup) {

	if strings.Contains(record[altIdx], ",") {
		processMultiLine(record, header, emptyField, fieldDelimiter, retainId, retainInfo, c, wg)
	} else {
		processSingleLine(record, header, emptyField, fieldDelimiter, retainId, retainInfo, c, wg)
	}
}

func processMultiLine(record []string, header []string, emptyField string,
	fieldDelimiter string, retainId bool, retainInfo bool, results chan<- string, wg *sync.WaitGroup) {

	defer wg.Done()

	chr := chrToUCSC(record[chromIdx])

	var homs []string
	var hets []string
	for idx, allele := range strings.Split(record[altIdx], ",") {
		if altIsValid(allele) == false {
			log.Printf("%s:%s Skip ALT #%d (not ACTG)", record[chromIdx], record[posIdx], idx+1)
			continue
		}

		siteType, pos, ref, alt, err := updateFieldsWithAlt(record[refIdx], allele, record[posIdx])

		if err != nil {
			log.Fatal(err)
		}

		if pos == "" {
			log.Printf("%s:%s Skip ALT #%d (complex)", record[chromIdx], record[posIdx], idx+1)
			continue
		}

		// If no sampels are provided, annotate what we can, skipping hets and homs
		if len(header) > 9 {
			homs = homs[:0]
			hets = hets[:0]

			err = makeHetHomozygotes(record, header, &homs, &hets, strconv.Itoa(idx+1))

			if err != nil {
				log.Fatal(err)
			}

			if len(homs) == 0 && len(hets) == 0 {
				continue
			}
		}

		var output bytes.Buffer

		output.WriteString(chr)
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
			output.WriteString(emptyField)
		} else {
			output.WriteString(strings.Join(hets, fieldDelimiter))
		}

		output.WriteString("\t")

		if len(homs) == 0 {
			output.WriteString(emptyField)
		} else {
			output.WriteString(strings.Join(homs, fieldDelimiter))
		}

		if retainId == true {
			output.WriteString("\t")
			output.WriteString(record[idIdx])
		}

		if retainInfo == true {
			// Write the index of the allele, to allow users to segregate data in the INFO field
			output.WriteString("\t")
			output.WriteString(strconv.Itoa(idx))
			output.WriteString("\t")
			// INFO index is 7
			output.WriteString(record[infoIdx])
		}

		output.WriteString("\n")

		results <- output.String()
	}
}

func processSingleLine(record []string, header []string,
	emptyField string, fieldDelimiter string, retainId bool, retainInfo bool, results chan<- string, wg *sync.WaitGroup) {
	defer wg.Done()

	if altIsValid(record[altIdx]) == false {
		log.Printf("%s:%s Skip ALT (not ACTG)", record[chromIdx], record[posIdx])
		return
	}

	var homs []string
	var hets []string

	siteType, pos, ref, alt, err := updateFieldsWithAlt(record[refIdx], record[altIdx], record[posIdx])

	if err != nil {
		log.Fatal(err)
	}

	if pos == "" {
		log.Printf("%s:%s Skip ALT (complex)", record[chromIdx], record[posIdx])
		return
	}

	// If no sampels are provided, annotate what we can, skipping hets and homs
	if len(header) > 9 {
		err = makeHetHomozygotes(record, header, &homs, &hets, "1")

		if err != nil {
			log.Fatal(err)
		}

		if len(homs) == 0 && len(hets) == 0 {
			return
		}
	}

	var output bytes.Buffer

	output.WriteString(chrToUCSC(record[chromIdx]))
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
		output.WriteString(emptyField)
	} else {
		output.WriteString(strings.Join(hets, fieldDelimiter))
	}

	output.WriteString("\t")

	if len(homs) == 0 {
		output.WriteString(emptyField)
	} else {
		output.WriteString(strings.Join(homs, fieldDelimiter))
	}

	if retainId == true {
		output.WriteString("\t")
		output.WriteString(record[idIdx])
	}

	if retainInfo == true {
		// Write the index of the allele, to allow users to segregate data in the INFO field
		// Of course in singl allele case, index is 0 (index is relative to alt alleles, not ref + alt)
		output.WriteString("\t")
		output.WriteString("0")
		output.WriteString("\t")
		// INFO index is 7
		output.WriteString(record[infoIdx])
		output.WriteString("\t")
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
	simpleGT := fields[formatIdx] == "GT"

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

func normalizeSampleNames(header []string) {
	for i := 9; i < len(header); i++ {
		header[i] = strings.Replace(header[i], ".", "_", -1)
	}
}

// func coercePositionsAlleles(emptyField string, primaryDelim string) func(string) bool {
