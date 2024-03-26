package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"regexp"
	"runtime"
	"runtime/pprof"
	"strconv"
	"strings"
	"sync"

	"github.com/akotlar/bystro-utils/parse"
	"github.com/apache/arrow/go/v14/arrow"
	"github.com/apache/arrow/go/v14/arrow/ipc"
	bystroArrow "github.com/bystrogenomics/bystro-vcf/arrow"
)

var fileMutex sync.Mutex

var concurrency = runtime.NumCPU()

const (
	chromIdx  int = 0
	posIdx    int = 1
	idIdx     int = 2
	refIdx    int = 3
	altIdx    int = 4
	qualIdx   int = 5
	filterIdx int = 6
	infoIdx   int = 7
	formatIdx int = 8
	sampleIdx int = 9
)

const (
	insError1      string = "1st base ALT != REF"
	delError1      string = "1st base REF != ALT"
	badPrefixError string = "No shared suffix && prefix doesn't match"
	sameError      string = "REF == ALT"
	missError      string = "ALT == '.'"
	posError       string = "Invalid POS"
	badAltError    string = "ALT not ACTG"
	mixedError     string = "Mixed indel/snp sites not supported"
	complexNotice  string = "Complex site"

	errorLvl string = "Error: "
)

const tabByte = byte('\t')
const clByte = byte('\n')
const chrByte = byte('c')
const zeroByte = byte('0')

// Decimal places to round floats to
const precision = 3

type Config struct {
	inPath              string
	outPath             string
	noOut               bool
	dosageMatrixOutPath string
	sampleListPath      string
	famPath             string
	errPath             string
	emptyField          string
	fieldDelimiter      string
	keepID              bool
	keepInfo            bool
	keepQual            bool
	keepPos             bool
	cpuProfile          string
	allowedFilters      map[string]bool
	excludedFilters     map[string]bool
}

func setup(args []string) *Config {
	config := &Config{}
	flag.StringVar(&config.inPath, "in", "", "The input file path (optional: default stdin)")
	flag.StringVar(&config.famPath, "fam", "", "The fam file path (optional)")
	flag.StringVar(&config.errPath, "err", "", "The log path (optional: default stderr)")
	flag.StringVar(&config.outPath, "out", "", "The output path (optional: default stdout)")
	flag.BoolVar(&config.noOut, "noOut", false, "Skip writing output (useful in conjunction with dosageOutput)")
	flag.StringVar(&config.dosageMatrixOutPath, "dosageOutput", "", "The output path for the dosage matrix (optional). If not provided, dosage matrix will not be output.")
	flag.StringVar(&config.sampleListPath, "sample", "", "The output path of the sample list (optional: default stdout)")
	flag.StringVar(&config.emptyField, "emptyField", "!", "The output path for the JSON output (optional)")
	flag.StringVar(&config.fieldDelimiter, "fieldDelimiter", ";", "The output path for the JSON output (optional)")
	flag.BoolVar(&config.keepID, "keepId", false, "Retain the ID field in output")
	flag.BoolVar(&config.keepQual, "keepQual", false, "Retain the QUAL field in output")
	flag.BoolVar(&config.keepPos, "keepPos", false, "Retain the original VCF position in output")
	flag.BoolVar(&config.keepInfo, "keepInfo", false, "Retain INFO field in output (2 appended output fields: allele index and the INFO field. Will appear after id field if --keepId flag set.")
	flag.StringVar(&config.cpuProfile, "cpuProfile", "", "Write cpu profile to file at this path")
	filteredVals := flag.String("allowFilter", "PASS,.", "Allow rows that have this FILTER value (comma separated)")
	excludeFilterVals := flag.String("excludeFilter", "", "Exclude rows that have this FILTER value (comma separated)")
	// allows args to be mocked https://github.com/nwjlyons/email/blob/master/inputs.go
	// can only run 1 such test, else, redefined flags error
	a := os.Args[1:]
	if args != nil {
		a = args
	}
	flag.CommandLine.Parse(a)

	if *filteredVals != "" && *filteredVals != "*" {
		config.allowedFilters = make(map[string]bool)

		for _, val := range strings.Split(*filteredVals, ",") {
			config.allowedFilters[strings.TrimSpace(val)] = true
		}
	}

	// We don't allow exclude all, that would be nonsensical
	if *excludeFilterVals != "" {
		config.excludedFilters = make(map[string]bool)

		for _, val := range strings.Split(*excludeFilterVals, ",") {
			config.excludedFilters[strings.TrimSpace(val)] = true
		}
	}

	return config
}

func init() {
	log.SetFlags(0)
}

// NOTE: For now this only supports \n end of line characters
// If we want to support CLRF or whatever, use either csv package, or set a different delimiter
func main() {
	config := setup(nil)

	inFh := (*os.File)(nil)
	if config.inPath != "" {
		var err error

		inFh, err = os.Open(config.inPath)
		if err != nil {
			log.Fatal(err)
		}
	} else {
		inFh = os.Stdin
	}

	defer inFh.Close()
	if config.errPath != "" {
		var err error
		os.Stderr, err = os.Open(config.errPath)
		if err != nil {
			log.Fatal(err)
		}
	}

	outFh := (*os.File)(nil)

	if config.noOut && config.outPath != "" {
		log.Fatal("Cannot specify --noOut and --out")
	}

	if config.noOut && config.dosageMatrixOutPath == "" {
		log.Fatal("When specifying --noOut, must specify --dosageOutput")
	}

	if !config.noOut {
		if config.outPath != "" {
			var err error

			outFh, err = os.OpenFile(config.outPath, os.O_WRONLY|os.O_CREATE, 0644)
			if err != nil {
				log.Fatal(err)
			}
		} else {
			outFh = os.Stdout
		}

		defer outFh.Close()
	}

	if config.cpuProfile != "" {
		f, err := os.Create(config.cpuProfile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	reader := bufio.NewReaderSize(inFh, 48*1024*1024)

	var writer *bufio.Writer

	if !config.noOut {
		writer = bufio.NewWriterSize(outFh, 48*1024*1024)

		fmt.Fprintln(writer, stringHeader(config))
	}

	readVcf(config, reader, writer)

	if !config.noOut {
		err := writer.Flush()

		if err != nil {
			log.Fatal(err)
		}

		err = outFh.Close()

		if err != nil {
			log.Print(err)
		}
	}
}

func stringHeader(config *Config) string {
	return strings.Join(header(config), string(tabByte))
}

func header(config *Config) []string {
	header := parse.Header

	if config.keepPos {
		header = append(header, "vcfPos")
	}

	if config.keepID {
		header = append(header, "id")
	}

	if config.keepInfo {
		header = append(header, "alleleIdx", "info")
	}

	return header
}

func readVcf(config *Config, reader *bufio.Reader, writer *bufio.Writer) {
	foundHeader := false

	var header []string

	// Read buffer
	workQueue := make(chan [][]byte, 16)
	complete := make(chan bool)

	endOfLineByte, numChars, versionLine, err := parse.FindEndOfLine(reader, "")

	if err != nil {
		log.Fatal(err)
	}

	vcfMatch, err := regexp.MatchString("##fileformat=VCFv4", versionLine)

	if err != nil {
		log.Fatal(err)
	}

	if !vcfMatch {
		log.Fatal("Not a VCF file")
	}

	for {
		// http://stackoverflow.com/questions/8757389/reading-file-line-by-line-in-go
		// http://www.jeffduckett.com/blog/551119d6c6b86364cef12da7/golang---read-a-file-line-by-line.html
		row, err := reader.ReadString(endOfLineByte) // 0x0A separator = newline

		if err == io.EOF {
			break
		} else if err != nil {
			log.Fatal(err)
		} else if row == "" {
			// This shouldn't occur, however, in case
			continue
		}

		// Chomp equivalent: https://groups.google.com/forum/#!topic/golang-nuts/smFU8TytFr4
		record := strings.Split(row[:len(row)-numChars], "\t")

		if foundHeader == false {
			if record[chromIdx] == "#CHROM" {
				header = record
				foundHeader = true
				break
			}
		}
	}

	if !foundHeader {
		log.Fatal("No header found")
	}

	parse.NormalizeHeader(header)

	if !config.noOut {
		err = writeSampleListIfWanted(config, header)

		if err != nil {
			log.Fatal("Couldn't write sample list file")
		}
	}

	var arrowWriter *bystroArrow.ArrowWriter
	if config.dosageMatrixOutPath != "" {
		if len(header) <= sampleIdx {
			log.Print("No samples found in VCF file; writing empty dosage matrix file")
			// Write empty file
			file, err := os.Create(config.dosageMatrixOutPath)
			if err != nil {
				log.Fatal(err)
			}

			file.Close()

			config.dosageMatrixOutPath = ""
		} else {
			sampleNames := header[sampleIdx:]

			fieldNames := append([]string{"locus"}, sampleNames...)

			fieldTypes := make([]arrow.DataType, len(fieldNames))
			fieldTypes[0] = arrow.BinaryTypes.String
			for i := 1; i < len(fieldNames); i++ {
				fieldTypes[i] = arrow.PrimitiveTypes.Uint16
			}

			file, err := os.Create(config.dosageMatrixOutPath)
			if err != nil {
				log.Fatal(err)
			}
			defer file.Close()

			arrowWriter, err = bystroArrow.NewArrowIPCFileWriter(file, fieldNames, fieldTypes, ipc.WithZstd())
			if err != nil {
				log.Fatal(err)
			}
			defer arrowWriter.Close()
		}
	}

	// Spawn threads
	for i := 0; i < concurrency; i++ {
		go processLines(header, numChars, config, workQueue, writer, complete, arrowWriter)
	}

	maxCapacity := 64

	// Fill the work queue.
	buff := make([][]byte, 0, maxCapacity)
	for {
		row, err := reader.ReadBytes(endOfLineByte) // 0x0A separator = newline

		if err == io.EOF {
			break
		} else if err != nil {
			log.Fatal(err)
		} else if len(row) == 0 {
			// We may have not closed the pipe, but not have any more information to send
			continue
		}

		if len(buff) >= maxCapacity {
			workQueue <- buff

			// if we re-assign it it will data race
			// i.e don't do buff = buff[:0]
			// buff = nil also works, but will set capacity to 0
			buff = make([][]byte, 0, maxCapacity)
		}

		buff = append(buff, row)
	}

	if len(buff) > 0 {
		workQueue <- buff
		buff = nil
	}

	// Indicate to all processing threads that no more work remains
	close(workQueue)

	// Wait for everyone to finish.
	for i := 0; i < concurrency; i++ {
		<-complete
	}

	if arrowWriter != nil {
		err = arrowWriter.Close()
		if err != nil {
			log.Fatal(err)
		}
	}
}

func writeSampleListIfWanted(config *Config, header []string) error {
	if config.sampleListPath == "" {
		return nil
	}

	outFh, err := os.OpenFile(config.sampleListPath, os.O_WRONLY|os.O_CREATE, 0644)

	if err != nil {
		return err
	}

	sList := makeSampleList(header)

	_, err = outFh.WriteString(sList.String())

	if err != nil {
		return err
	}

	err = outFh.Sync()

	if err != nil {
		return err
	}

	err = outFh.Close()

	if err != nil {
		return err
	}

	return nil
}

func makeSampleList(header []string) bytes.Buffer {
	var buf bytes.Buffer

	if len(header) < 10 {
		return buf
	}

	for i := sampleIdx; i < len(header); i++ {
		buf.WriteString(header[i])
		buf.WriteByte(clByte)
	}

	return buf
}

func linePasses(record []string, header []string, allowedFilters map[string]bool,
	excludedFilters map[string]bool) bool {
	return len(record) == len(header) &&
		// whitelist: if true it's present in the map
		(allowedFilters == nil || allowedFilters[record[filterIdx]] == true) &&
		// blacklist: if false it's not present in the map, and we allow it
		(excludedFilters == nil || excludedFilters[record[filterIdx]] == false)
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

		for i := 1; i < len(alt); i++ {
			if alt[i] != 'A' && alt[i] != 'C' && alt[i] != 'T' && alt[i] != 'G' {
				return false
			}
		}
	}
	return true
}

func processLines(header []string, numChars int, config *Config, queue chan [][]byte,
	writer *bufio.Writer, complete chan bool, arrowWriter *bystroArrow.ArrowWriter) {
	var multiallelic bool

	// Declare sample-related variables outside loop, in case this helps us
	// reduce allocations
	// Safe, because the property of having samples is invariant across lines within
	// a single file
	var numSamples float64
	var homs []string
	var hets []string
	var missing []string
	var dosages []any
	var effectiveSamples float64
	var ac int
	var an int
	var chrom string

	emptyField := config.emptyField
	fieldDelim := config.fieldDelimiter
	keepID := config.keepID
	keepInfo := config.keepInfo
	allowedFilters := config.allowedFilters
	excludedFilters := config.excludedFilters
	keepPos := config.keepPos

	needsLabels := !config.noOut
	needsDosages := config.dosageMatrixOutPath != ""

	if len(header) > sampleIdx {
		numSamples = float64(len(header) - sampleIdx)
	} else if len(header) == sampleIdx {
		log.Printf("Found 9 header fields. When genotypes present, we expect 1+ samples after FORMAT (10 fields minimum)")
	}

	var output bytes.Buffer
	var record []string

	var arrowBuilder *bystroArrow.ArrowRowBuilder
	var err error
	if arrowWriter != nil {
		arrowBuilder, err = bystroArrow.NewArrowRowBuilder(arrowWriter, 5e3)
		if err != nil {
			log.Fatal(err)
		}
	}

	for lines := range queue {
		if !config.noOut && output.Len() >= 2e6 {
			fileMutex.Lock()

			writer.Write(output.Bytes())

			fileMutex.Unlock()

			output.Reset()
		}

		for _, row := range lines {
			record = strings.Split(string(row[:len(row)-numChars]), "\t")

			if !linePasses(record, header, allowedFilters, excludedFilters) {
				continue
			}

			siteType, positions, refs, alts, altIndices := getAlleles(record[chromIdx], record[posIdx], record[refIdx], record[altIdx])

			if len(altIndices) == 0 {
				continue
			}

			multiallelic = siteType == parse.Multi

			for i := range alts {
				var arrowRow []any

				strAlt := strconv.Itoa(altIndices[i] + 1)
				// If no samples are provided, annotate what we can, skipping hets and homs
				// If samples are provided, but only missing genotypes, skip the site altogether
				if numSamples > 0 {
					homs, hets, missing, dosages, ac, an = makeHetHomozygotes(record, header, strAlt, needsLabels, needsDosages)

					if ac == 0 {
						continue
					}

					// homozygosity and heterozygosity should be relative to complete genotypes
					effectiveSamples = numSamples - float64(len(missing))
				}

				// output is [chr, pos, type, ref, alt, trTv, het, heterozygosity, hom, homozygosity, missing, missingness, sampleMaf]
				// if keepID append id
				// if keepInfo append [alleleIndex, info]

				if len(record[chromIdx]) < 4 || record[chromIdx][0] != chrByte {
					chrom = "chr" + record[chromIdx]
				} else {
					chrom = record[chromIdx]
				}

				if arrowBuilder != nil {
					arrowRow = append(arrowRow, fmt.Sprintf("%s:%s:%s:%s", chrom, positions[i], string(refs[i]), alts[i]))

					if numSamples > 0 {
						arrowRow = append(arrowRow, dosages...)
					}

					arrowBuilder.WriteRow(arrowRow)
				}

				if needsLabels {
					output.WriteString(chrom)
					output.WriteByte(tabByte)

					output.WriteString(positions[i])
					output.WriteByte(tabByte)

					output.WriteString(siteType)
					output.WriteByte(tabByte)

					output.WriteByte(refs[i])
					output.WriteByte(tabByte)

					output.WriteString(alts[i])
					output.WriteByte(tabByte)

					if multiallelic {
						output.WriteString(parse.NotTrTv)
					} else {
						output.WriteString(parse.GetTrTv(string(refs[i]), alts[i]))
					}

					output.WriteByte(tabByte)

					// Write missing samples
					// heterozygotes \t heterozygosity
					if len(hets) == 0 {
						output.WriteString(emptyField)
						output.WriteByte(tabByte)
						output.WriteByte(zeroByte)
					} else {
						output.WriteString(strings.Join(hets, fieldDelim))
						output.WriteByte(tabByte)

						// This gives plenty precision; we are mostly interested in
						// the first or maybe 2-3 significant digits
						// https://play.golang.org/p/Ux-QmClaJG
						// Also, gnomAD seems to use 6 bits of precision
						// the bitSize == 64 allows us to round properly past 6 s.f
						// Note: 'G' requires these numbers to be < 0 for proper precision
						// (elase only 6 s.f total, rather than after decimal)
						output.WriteString(strconv.FormatFloat(float64(len(hets))/effectiveSamples, 'G', precision, 64))
					}

					output.WriteByte(tabByte)

					// Write missing samples
					// homozygotes \t homozygosity
					if len(homs) == 0 {
						output.WriteString(emptyField)
						output.WriteByte(tabByte)
						output.WriteByte(zeroByte)
					} else {
						output.WriteString(strings.Join(homs, fieldDelim))
						output.WriteByte(tabByte)
						output.WriteString(strconv.FormatFloat(float64(len(homs))/effectiveSamples, 'G', precision, 64))
					}

					output.WriteByte(tabByte)

					// Write missing samples
					// missingGenos \t missingness
					if len(missing) == 0 {
						output.WriteString(emptyField)
						output.WriteByte(tabByte)
						output.WriteByte(zeroByte)
					} else {
						output.WriteString(strings.Join(missing, fieldDelim))
						output.WriteByte(tabByte)
						output.WriteString(strconv.FormatFloat(float64(len(missing))/numSamples, 'G', precision, 64))
					}

					// Write the sample minor allele frequency
					output.WriteByte(tabByte)

					output.WriteString(strconv.Itoa(ac))
					output.WriteByte(tabByte)
					output.WriteString(strconv.Itoa(an))
					output.WriteByte(tabByte)

					// TODO: can ac == 0 && (len(het) > 0 || len(hom) > 0) occur?
					if ac == 0 {
						output.WriteByte(zeroByte)
					} else {
						output.WriteString(strconv.FormatFloat(float64(ac)/float64(an), 'G', precision, 64))
					}

					/******************* Optional Fields ***********************/
					if keepPos == true {
						output.WriteByte(tabByte)
						output.WriteString(record[posIdx])
					}

					if keepID == true {
						output.WriteByte(tabByte)
						output.WriteString(record[idIdx])
					}

					if keepInfo == true {
						// Write the index of the allele, to allow users to segregate data in the INFO field
						output.WriteByte(tabByte)
						output.WriteString(strconv.Itoa(altIndices[i]))

						// Write info for all indices
						output.WriteByte(tabByte)
						output.WriteString(record[infoIdx])
					}

					output.WriteByte(clByte)
				}

			}
		}

		if arrowBuilder != nil {
			arrowBuilder.WriteRow(nil)
		}
	}

	if !config.noOut && output.Len() > 0 {
		fileMutex.Lock()

		writer.Write(output.Bytes())

		fileMutex.Unlock()
	}

	if arrowBuilder != nil {
		err = arrowBuilder.Release()
		if err != nil {
			log.Fatal(err)
		}
	}

	complete <- true
}

func getAlleles(chrom string, pos string, ref string, alt string) (string, []string, []byte, []string, []int) {
	// Indel format:
	// Deletion: -N : "-" followed by # of deleted bases (inclusive of ref)
	// Insertion: \+[ATCG]+ : "+" followed by the inserted bases, which occur after the ref
	// ref is always 1 base,  which follows from not requiring deletions to have a ACTG base

	if alt == ref {
		log.Printf("%s:%s : %s\n", chrom, pos, sameError)
		return "", nil, nil, nil, nil
	}

	// optimize for the cases where no "," could be present, i.e len(alt) == 1
	if len(alt) == 1 {
		if alt != "A" && alt != "C" && alt != "G" && alt != "T" {
			log.Printf("%s:%s ALT #1 %s\n", chrom, pos, badAltError)

			return "", nil, nil, nil, nil
		}

		if len(ref) == 1 {
			return parse.Snp, []string{pos}, []byte{ref[0]}, []string{alt}, []int{0}
		}

		// simple deletion must have 1 base padding match
		if alt[0] != ref[0] {
			log.Printf("%s:%s ALT #1 %s", chrom, pos, delError1)
			return "", nil, nil, nil, nil
		}

		intPos, err := strconv.Atoi(pos)

		if err != nil {
			log.Printf("%s:%s ALT #1 %s", chrom, pos, posError)
			return "", nil, nil, nil, nil
		}

		// pos is the next base over (first deleted base)
		// ref is also the first deleted base, since alt is of 1 padding, that's idx 1 (2nd ref base)
		// alt == len(alt) - len(ref) for len(alt) < len(ref)
		// example: alt = A (len == 1), ref = AAATCC (len == 6)
		// 1 - 6 = -5 (then conver to string)
		return parse.Del, []string{strconv.Itoa(intPos + 1)}, []byte{ref[1]}, []string{strconv.Itoa(1 - len(ref))}, []int{0}
	}

	var intPos int

	var alleles []string
	var references []byte
	var positions []string
	var indexes []int
	var multi bool
	for altIdx, tAlt := range strings.Split(alt, ",") {
		// It can be a MULTIALLELIC and have errors that
		// reduce output allele count to 1, so record here
		if !multi && altIdx > 0 {
			multi = true
		}

		if altIsValid(tAlt) == false {
			log.Printf("%s:%s ALT #%d %s\n", chrom, pos, altIdx+1, badAltError)
			continue
		}

		if len(ref) == 1 {
			if len(tAlt) == 1 {
				//tAlt isn't modified
				positions = append(positions, pos)
				references = append(references, ref[0])
				alleles = append(alleles, tAlt)
				indexes = append(indexes, altIdx)
				continue
			}

			// Simple insertion : tAlt > 1 base, and ref == 1 base
			if tAlt[0] != ref[0] {
				log.Printf("%s:%s ALT #%d %s", chrom, pos, altIdx+1, insError1)
				continue
			}

			// we take the allele from the 2nd base, since insertion occurs after the reference
			var buffer bytes.Buffer
			buffer.WriteString("+")
			buffer.WriteString(tAlt[1:])

			// simple insertion; our annotations also use 1 base padding for insertions, so keep pos same
			positions = append(positions, pos)
			// since pos same, ref is just the first (only) base
			references = append(references, ref[0])
			alleles = append(alleles, buffer.String())
			indexes = append(indexes, altIdx)

			continue
		}

		// len(ref) > 1

		// If given 0-based file, this will be re-generated potentially
		// Notice we exit the loop here; position is invalid, should leave
		// We convert Atoi here to avoid wasting performance
		if intPos == 0 {
			var err error
			intPos, err = strconv.Atoi(pos)

			if err != nil {
				log.Printf("%s:%s %s", chrom, pos, posError)
				break
			}
		}

		if len(tAlt) == 1 {
			// Simple deletion, padding of 1 base, padding must match
			if tAlt[0] != ref[0] {
				log.Printf("%s:%s ALT#%d %s", chrom, pos, altIdx+1, delError1)
				continue
			}

			// 1 base deletion; we use 0 padding for deletions, showing 1st deleted base
			// as ref; so shift pos, ref by 1, return len(ref) - 1 for alt
			positions = append(positions, strconv.Itoa(intPos+1))
			references = append(references, ref[1])
			alleles = append(alleles, strconv.Itoa(1-len(ref)))
			indexes = append(indexes, altIdx)

			continue
		}

		// If we're here, ref and alt are both > 1 base long
		// could be a weird SNP (multiple bases are SNPS, len(ref) == len(alt))
		// could be a weird deletion/insertion
		// could be a completely normal multiallelic (due to padding, shifted)

		//1st check for MNPs and extra-padding SNPs
		if len(ref) == len(tAlt) {
			// Let's check each base; if there is more than 1 change, this is an MNP
			// whether a sparse MNP or not.
			// We'll report each modified base
			for i := 0; i < len(ref); i++ {
				if ref[i] != tAlt[i] {
					positions = append(positions, strconv.Itoa(intPos+i))
					references = append(references, ref[i])
					alleles = append(alleles, string(tAlt[i]))

					// Here we append the index of the VCF ALT, not the MNP idx
					// all bases in an MNP will have the same index
					// else we won't be able to determine homozygous, heterozygous, reference status
					indexes = append(indexes, altIdx)
				}
			}

			continue
		}

		// Find the allele representation that minimizes padding, while still checking
		// that the site isn't a mixed type (indel + snp) and checking for intercolation
		// Essentially, Occam's Razor for padding: minimize the number of steps away
		// from left edge to explan the allele
		// EX:
		// If ref == AATCG
		// If alt == AG
		// One interpretation of this site is mixed A->G -3 (-TCG)
		// Another is -3 (-ATC) between the A (0-index) and G (4-index) in ref
		// We prefer the latter approach
		// Ex2: ref: TT alt: TCGATT
		// We prefer +CGAT

		// Like http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/
		// we will use a simple heuristic:
		// 1) For insertions, figure out the shared right edge, from 1 base downstream of first ref base
		// Then, check if the remaining ref bases match the left edge of the alt
		// If they don't, skip that site
		// 2) For deletions, the same, except replace the role of ref with the tAlt

		// Our method should be substantially faster, since we don't need to calculate
		// the min(len(ref), len(tAlt))
		// and because we don't create a new slice for every shared ref/alt at right edges and left

		if len(tAlt) > len(ref) {
			rIdx := 0
			// we won't allow the entire ref to match as a suffix, only up to 1 base downstream from beginning of ref
			for len(tAlt)+rIdx > 0 && len(ref)+rIdx > 1 && tAlt[len(tAlt)+rIdx-1] == ref[len(ref)+rIdx-1] {
				rIdx--
			}

			// Then, we require an exact match from left edge, for the difference between the
			// length of the ref, and the shared suffix
			// Ex: alt: TAGCTT ref: TAT
			// We shared 1 base at right edge, so expect that len(ref) - 1, or 3 - 1 = 2 bases of ref
			// match the left edge of alt
			// Here that is TA, for an insertion of +GCT
			// Ex2: alt: TAGCAT ref: TAT
			// Here the AT of the ref matches the last 2 bases of alt
			// So we expect len(ref) - 2 == 1 base of ref to match left edge of the alt (T), for +AGC
			// Ex3: alt TAGTAT ref: TAT
			// Since our loop doesn't check the last base of ref, as in ex2, +AGC
			// This mean we always prefer a 1-base padding, when possible
			// Ex4: alt TAGTAT ref: TGG
			// In this case, we require len(ref) - 0 bases in the ref to match left edge of alt
			// Since they don't (TAG != TGG), we call this complex and move on

			// Insertion
			// If pos is 100 and ref is AATCG
			// and alt is AAAAATCG (len == 7)
			// we expect lIdx to be 2
			// and rIdx to be -3
			// alt[2] is the first non-ref base
			// and alt[len(alt) - 3] == alt[4] is the last non-ref base
			// The position is intPos + lIdx or 100 + 2 - 1 == 101 (100, 101 are padding bases,
			// and we want to keep the last reference base
			// The ref is ref[2 - 1] or ref[1]
			offset := len(ref) + rIdx
			if ref[:offset] != tAlt[:offset] {
				log.Printf("%s:%s ALT#%d %s", chrom, pos, altIdx+1, mixedError)
				continue
			}

			// position is offset by len(ref) + 1 - rIdx
			// ex1: alt: TAGCTT ref: TAT
			// here we match the first base, so -1
			// we require remainder of left edge to be present,
			// or len(ref) - 1 == 2
			// so intPos + 2 - 1 for last padding base (the A in TA) (intPos + 2 is first unique base)
			positions = append(positions, strconv.Itoa(intPos+offset-1))
			references = append(references, ref[offset-1])

			// Similarly, the alt allele starts from len(ref) + rIdx, and ends at len(tAlt) + rIdx
			// from ex: TAGCTT ref: TAT :
			// rIdx == -1 , real alt == tAlt[len(ref) - 1:len(tAlt) - 1] == tALt[2:5]
			var insBuffer bytes.Buffer
			insBuffer.WriteString("+")
			insBuffer.WriteString(tAlt[offset : len(tAlt)+rIdx])

			alleles = append(alleles, insBuffer.String())
			indexes = append(indexes, altIdx)

			continue
		}

		// Deletion
		// If pos is 100 and alt is AATCG
		// and ref is AAAAATCG (len == 7)
		// we expect lIdx to be 2
		// and rIdx to be -3
		// and alt is -3 or len(ref) + rIdx - lIdx == 8 + -3 - 2
		// position is the first deleted base, or intPos + lIdx == 100 + 2 == 102
		// where (100, 101) are the two padding bases
		// ref is the first deleted base or ref[lIdx] == ref[2]

		// Just like insertion, but try to match all bases from 1 base downstream of tAlt to ref
		rIdx := 0
		// we won't allow the enitre ALT to match as a suffix
		// this is the inverse of the insertion case
		// can have an intercolated deletion, where some bases in middle of deletion are missing w.r.t ref
		// or deletions where teh entire left edge of the alt match the ref
		for len(tAlt)+rIdx > 1 && len(ref)+rIdx > 0 && tAlt[len(tAlt)+rIdx-1] == ref[len(ref)+rIdx-1] {
			rIdx--
		}

		//ex: if post: 100 alt: AATCG and ref: AAAAATCG, we expect deletion to
		//occurs on base 101, alt: -AAAA == -4 == -(len(ref) - 3 - (len(tAlt) - 3)) == -8 - 3 - 2 = -3
		//rIdx is 3 since matches 3 bases
		//position gets shifted by len(tAlt) + rIdx, since we don't want any padding in our output
		offset := len(tAlt) + rIdx
		if ref[:offset] != tAlt[:offset] {
			log.Printf("%s:%s ALT#%d %s", chrom, pos, altIdx+1, mixedError)
			continue
		}

		positions = append(positions, strconv.Itoa(intPos+offset))
		// we want the base after the last shared
		references = append(references, ref[offset])

		// the allele if -(len(ref) + rIdx - (len(tAlt) + rIdx))
		alleles = append(alleles, strconv.Itoa(-(len(ref) + rIdx - offset)))
		indexes = append(indexes, altIdx)

		continue
	}

	// Any of the logged errors generated in loop
	// This method saves us having to build up an "errors" slice
	// While still allowing us to report/log bases
	if len(alleles) == 0 {
		return "", nil, nil, nil, nil
	}

	// MULTI is basically the presence of a comma
	// We don't just look at the length of the resulting alleles or references or indexes
	// slices, because these can include MNPs, which expand the array size, but are
	// somewhat distinct in their meaning
	if multi {
		return parse.Multi, positions, references, alleles, indexes
	}

	// Anything below here is guaranteed to be a single allele from the VCF file
	// whether that is an MNP, DEL, INS, or SNP
	if len(alleles[0]) > 1 {
		if alleles[0][0] == '-' {
			return parse.Del, positions, references, alleles, indexes
		}

		return parse.Ins, positions, references, alleles, indexes
	}

	// MNPs result in multiple alleles
	// They may be sparse or complete, so we empirically check for their presence
	// If the MNP is really just a snp, there is only 1 allele, and reduces to snp
	// > 1, these are labeled differently to allow people to jointly consider the effects
	// of the array of SNPs, since we at the moment consider their effects only independently
	// (which has advantages for CADD, phyloP, phastCons, clinvar, etc reporting)
	if len(alleles) > 1 {
		return parse.Mnp, positions, references, alleles, indexes
	}

	// MNPs and SNPs are both labeled SNP
	return parse.Snp, positions, references, alleles, indexes
}

// makeHetHomozygotes process all sample genotype fields, and for a single alleleNum, which is the allele index (1 based)
// returns the homozygotes, heterozygotes, missing samples, total alt counts, genotype counts, and missing counts
func makeHetHomozygotes(fields []string, header []string, alleleNum string, needsLabels bool, needsDosages bool) ([]string, []string, []string, []any, int, int) {
	var homs []string
	var hets []string
	var missing []string
	var dosages []any

	var gtCount int
	var altCount int

	var totalAltCount int
	var totalGtCount int

SAMPLES:
	// NOTE: If any errors encountered, all genotypes in row will be skipped and logged, since
	// this represents a likely corruption of data
	for i := sampleIdx; i < len(header); i++ {
		sampleGenotypeField := fields[i]

		// We want to speed up the common case, where the genotype is bi-allelic
		// e.g. we allow 0|0, 0/0, 0|1, 1|0, 1/0, 0/1, 1|1, 1/1
		// or those genotypes with other information, e.g. 0|0:DP:AD:GQ:PL
		if (len(sampleGenotypeField) == 3 || (len(sampleGenotypeField) > 3 && sampleGenotypeField[3] == ':')) && // 0|0, 0/0, 1|1, 1/1 or with format, e.g. 0|1:DP
			(sampleGenotypeField[1] == '|' || sampleGenotypeField[1] == '/') { // Not a haploid site with 100+ alleles, e.g. {0-9}|{0-9} or {0-9}/{0-9}
			// Reference is the most common case, e.g. 0|0, 0/0
			if sampleGenotypeField[0] == '0' && sampleGenotypeField[2] == '0' {
				totalGtCount += 2

				if needsDosages {
					dosages = append(dosages, uint16(0))
				}

				continue SAMPLES
			}

			// In this function, we only care about the alleleNum allele
			// Any diploid genotype with an allele that is longer than 1 number will be longer than 3 characters
			// E.g., if alleleNum is 10, then 0|10, 10|0, 10|10 will be the shortest possible genotype, 4 characters
			if len(alleleNum) == 1 {
				if (sampleGenotypeField[0] == '0' && sampleGenotypeField[2] == alleleNum[0]) || (sampleGenotypeField[0] == alleleNum[0] && sampleGenotypeField[2] == '0') {
					totalGtCount += 2
					totalAltCount += 1

					if needsLabels {
						hets = append(hets, header[i])
					}

					if needsDosages {
						dosages = append(dosages, uint16(1))
					}

					continue SAMPLES
				}

				// Homozygote
				if sampleGenotypeField[0] == alleleNum[0] && sampleGenotypeField[2] == alleleNum[0] {
					totalGtCount += 2
					totalAltCount += 2

					if needsLabels {
						homs = append(homs, header[i])
					}

					if needsDosages {
						dosages = append(dosages, uint16(2))
					}

					continue SAMPLES
				}
			}

			// N|., .|N, .|., N/., ./N, ./. are all considered missing samples, because if one site is missing, the other is likely unreliable
			if sampleGenotypeField[0] == '.' || sampleGenotypeField[2] == '.' {
				if needsLabels {
					missing = append(missing, header[i])
				}

				if needsDosages {
					dosages = append(dosages, nil)
				}

				continue SAMPLES
			}
		}

		// Split the field on the colon to separate alleles from additional information
		parts := strings.SplitN(sampleGenotypeField, ":", 2)
		alleleField := parts[0]

		var sep string
		if strings.Contains(alleleField, "|") {
			sep = "|"
		} else if strings.Contains(alleleField, "/") {
			sep = "/"
		} else {
			sep = ""
		}

		var alleles []string
		if sep != "" {
			alleles = strings.Split(alleleField, sep)
		} else {
			alleles = []string{alleleField}
		}

		altCount = 0
		gtCount = 0

		for _, allele := range alleles {
			if allele == "." {
				if needsLabels {
					missing = append(missing, header[i])
				}

				if needsDosages {
					dosages = append(dosages, nil)
				}

				continue SAMPLES
			}

			if allele == alleleNum {
				altCount++
			}

			gtCount++
		}

		totalGtCount += gtCount
		totalAltCount += altCount

		if needsDosages {
			dosages = append(dosages, uint16(altCount))
		}

		if altCount == 0 {
			continue
		}

		if needsLabels {
			if int(altCount) == gtCount {
				homs = append(homs, header[i])
			} else {
				hets = append(hets, header[i])
			}
		}
	}

	return homs, hets, missing, dosages, totalAltCount, totalGtCount
}
