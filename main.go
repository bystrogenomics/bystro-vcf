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
	inPath          string
	outPath         string
	sampleListPath  string
	famPath         string
	errPath         string
	emptyField      string
	fieldDelimiter  string
	keepID          bool
	keepInfo        bool
	keepQual        bool
	keepPos         bool
	cpuProfile      string
	allowedFilters  map[string]bool
	excludedFilters map[string]bool
}

func setup(args []string) *Config {
	config := &Config{}
	flag.StringVar(&config.inPath, "in", "", "The input file path (optional: default stdin)")
	flag.StringVar(&config.famPath, "fam", "", "The fam file path (optional)")
	flag.StringVar(&config.errPath, "err", "", "The log path (optional: default stderr)")
	flag.StringVar(&config.outPath, "out", "", "The output path (optional: default stdout")
	flag.StringVar(&config.sampleListPath, "sample", "", "The output path of the sample list (optional: default stdout")
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

	// make sure it gets closed
	defer inFh.Close()
	if config.errPath != "" {
		var err error
		os.Stderr, err = os.Open(config.errPath)
		if err != nil {
			log.Fatal(err)
		}
	}

	outFh := (*os.File)(nil)
	if config.outPath != "" {
		var err error

		outFh, err = os.OpenFile(config.outPath, os.O_WRONLY|os.O_CREATE, 0644)
		if err != nil {
			log.Fatal(err)
		}
	} else {
		outFh = os.Stdout
	}
	// make sure it gets closed
	defer outFh.Close()

	if config.cpuProfile != "" {
		f, err := os.Create(config.cpuProfile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	reader := bufio.NewReaderSize(inFh, 48*1024*1024)

	writer := bufio.NewWriterSize(outFh, 48*1024*1024)

	fmt.Fprintln(writer, stringHeader(config))

	// if config.famPath != "" {
	//   fillSampleIdx(config.famPath)
	// }

	readVcf(config, reader, writer)
	writer.Flush()
	outFh.Close()
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

// func fillSampleIdx (config *Config) {

// }

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

	// Get the header
	for {
		// http://stackoverflow.com/questions/8757389/reading-file-line-by-line-in-go
		// http://www.jeffduckett.com/blog/551119d6c6b86364cef12da7/golang---read-a-file-line-by-line.html
		// Scanner doesn't work well, has buffer restrictions that we need to manually get around
		// and we don't expect any newline characters in a Seqant output body
		row, err := reader.ReadString(endOfLineByte) // 0x0A separator = newline

		if err == io.EOF {
			break
		} else if err != nil {
			log.Fatal(err)
		} else if row == "" {
			// This shouldn't occur, however, in case
			continue
		}

		// remove the trailing \n or \r
		// equivalent of chomp https://groups.google.com/forum/#!topic/golang-nuts/smFU8TytFr4
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

	// Remove periods from sample names
	parse.NormalizeHeader(header)

	err = writeSampleListIfWanted(config, header)

	if err != nil {
		log.Fatal("Couldn't write sample list file")
	}

	// Now read them all off, concurrently.
	for i := 0; i < concurrency; i++ {
		go processLines(header, numChars, config, workQueue, writer, complete)
	}

	maxCapacity := 64
	// Read the lines into the work queue.
	// idx := 0
	buff := make([][]byte, 0, maxCapacity)
	for {
		row, err := reader.ReadBytes(endOfLineByte) // 0x0A separator = newline

		if err == io.EOF {
			break
		} else if err != nil {
			log.Fatal(err)
		} else if len(row) == 0 {
			// We may have not closed the pipe, but not have any more information to send
			// Wait for EOF
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

	// Close the channel so everyone reading from it knows we're done.
	close(workQueue)

	// Wait for everyone to finish.
	for i := 0; i < concurrency; i++ {
		<-complete
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

	// make sure it gets closed
	defer outFh.Close()
	writer := bufio.NewWriter(outFh)

	writeSampleList(writer, header)

	return nil
}

func writeSampleList(writer *bufio.Writer, header []string) {
	if len(header) < 10 {
		return
	}

	for i := 9; i < len(header); i++ {
		writer.WriteString(header[i])
		writer.WriteByte(clByte)
	}

	writer.Flush()
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

func processLines(header []string, numChars int, config *Config, queue chan [][]byte, writer *bufio.Writer, complete chan bool) {
	var multiallelic bool

	// Declare sample-related variables outside loop, in case this helps us
	// reduce allocations
	// Safe, because the property of having samples is invariant across lines within
	// a single file
	var numSamples float64
	var homs []string
	var hets []string
	var missing []string
	var effectiveSamples float64
	var ac int
	var an int

	emptyField := config.emptyField
	fieldDelim := config.fieldDelimiter
	keepID := config.keepID
	keepInfo := config.keepInfo
	allowedFilters := config.allowedFilters
	excludedFilters := config.excludedFilters
	keepPos := config.keepPos

	if len(header) > 9 {
		numSamples = float64(len(header) - 9)
	} else if len(header) == 9 {
		log.Printf("Found 9 header fields. When genotypes present, we expect 1+ samples after FORMAT (10 fields minimum)")
	}

	iLookup := []rune{'1', '2', '3', '4', '5', '6', '7', '8', '9'}

	var output bytes.Buffer
	var record []string

	for lines := range queue {
		if output.Len() >= 2e6 {
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

			// if last index > 9 then we can't accept the site, since won't be able
			// to identify het/hom status
			if altIndices[len(altIndices)-1] > 9 {
				log.Printf("%s %s:%s: We currently don't support sites with > 9 minor alleles, found %d", record[chromIdx], record[posIdx], errorLvl, len(altIndices))
				continue
			}

			for i := range alts {
				// If no samples are provided, annotate what we can, skipping hets and homs
				// If samples are provided, but only missing genotypes, skip the site altogether
				if numSamples > 0 {
					homs, hets, missing, ac, an = makeHetHomozygotes(record, header, iLookup[altIndices[i]])

					if len(homs) == 0 && len(hets) == 0 {
						continue
					}

					// homozygosity and heterozygosity should be relative to complete genotypes
					effectiveSamples = numSamples - float64(len(missing))
				}

				// output is [chr, pos, type, ref, alt, trTv, het, heterozygosity, hom, homozygosity, missing, missingness, sampleMaf]
				// if keepID append id
				// if keepInfo append [alleleIndex, info]

				if len(record[chromIdx]) < 4 || record[chromIdx][0] != chrByte {
					output.WriteString("chr")
				}

				output.WriteString(record[chromIdx])
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
				// This can be 0 in one of wo situations
				// First, if we have only missing genotypes at this site
				// However, in this case, we don't reach this code, because of line
				// 302 (if len(homs) == 0 && len(hets) == 0)
				// Else if there are truly no minor allele
				output.WriteByte(tabByte)

				output.WriteString(strconv.Itoa(ac))
				output.WriteByte(tabByte)
				output.WriteString(strconv.Itoa(an))
				output.WriteByte(tabByte)

				if ac == 0 {
					output.WriteByte(zeroByte)
				} else {
					output.WriteString(strconv.FormatFloat(float64(ac)/float64(an), 'G', precision, 64))
				}

				/******************* Optional Fields ***********************/
				if keepPos == true {
					// Write the input VCF position; we normalize this above
					// and here you can validate that transformation
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

	if output.Len() > 0 {
		fileMutex.Lock()

		writer.Write(output.Bytes())

		fileMutex.Unlock()
	}

	// log.Println("Worker hit, missed this many times: ", hitCount, missCount)
	// Let the main process know we're done.
	complete <- true
}

// @returns type, pos []string, ref []byte, alt []string, altIndices []int
// ref is always 1 base
// could return error object, but logging works as well, and seems cheaper
func getAlleles(chrom string, pos string, ref string, alt string) (string, []string, []byte, []string, []int) {
	// First check the simple cases, for performance reasons
	if alt == ref {
		// No point in returning ref sites
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
		// ref is also the first deleted base, since alt is of 1 padding, that's idx 1/ 2nd ref base
		// alt == len(alt) - len(ref) for len(alt) < len(ref)
		// example: alt = A (len == 1), ref = AAATCC (len == 6)
		// 1 - 6 = -5 (then conver to string)
		return parse.Del, []string{strconv.Itoa(intPos + 1)}, []byte{ref[1]}, []string{strconv.Itoa(1 - len(ref))}, []int{0}
	}

	// len(ref) > 1
	// this means position errors only generated if complex
	// consuming script will probaly also need position, and will generate its own
	// errors
	// we don't pad our deletions; the first deleted base is -1
	// and the reference is that first deleted base
	// don't modify this value, reuse
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

		// if a site doesn't have a valid form, like ACTG, skip it
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

		// I we're here, ref and alt are both > 1 base long
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

// Current limitations: Does not support alleleIdx > 9, or sites which have 2 digits allele numbers
// This allows us to improve performance, decrease code verbosity
// Most multialellics with > 10 alleles will be false positives
// because a true site would require a mutation rate of >> 1e-8 (say in chromosomal instability) and 10k samples
// or an effective population size of billions (vs 10k expected) ; else .001^10 == 10^-30 == never happens

// TODO: decide whether to be more strict about missing genotypes
// Currently some garbage like .... would be considered "missing"
func makeHetHomozygotes(fields []string, header []string, alleleNum rune) ([]string, []string, []string, int, int) {
	simpleGt := !strings.Contains(fields[formatIdx], ":")

	var homs []string
	var hets []string
	var missing []string

	var gt []string

	var gtCount int
	var altCount int

	var totalAltCount int
	var totalGtCount int

	// Unfortunately there is no guarantee that genotypes will be consistently phased or unphased
	/*  From https://samtools.github.io/hts-specs/VCFv4.1.pdf
	#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
	20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
	20 17330 . T A 3 q10 NS=3;DP=11;AF=0.017 GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3 0/0:41:3
	20 1110696 rs6040355 A G,T 67 PASS NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2 2/2:35:4
	20 1230237 . T . 47 PASS NS=3;DP=13;AA=T GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
	20 1234567 microsat1 GTC G,GTCT 50 PASS NS=3;DP=9;AA=G GT:GQ:DP 0/1:35:4 0/2:17:2 1/1:40:3
	*/
SAMPLES:
	// NOTE: If any errors encountered, all genotypes in row will be skipped and logged, since
	// this represents a likely corruption of data
	for i := 9; i < len(header); i++ {
		// If any 1 allele missing the genotype is, by definition missing
		// Applies to haploid as well as diploid+ sites
		if fields[i][0] == '.' {
			missing = append(missing, header[i])
			continue
		}

		// haploid
		if len(fields[i]) == 1 || fields[i][1] == ':' {
			if fields[i][0] == '0' {
				totalGtCount++
				continue
			}

			// We don't support haploid genotypes very well; I will count such sites
			// homozygous, because Dave Cutler says that is what people would mostly expect
			// Note: downstream tools will typicaly consider homozygotes to have 2 copies of the alt
			// and hets to have 1 copy of the alt, so special consideration must be made for haploids
			// in such tools
			if rune(fields[i][0]) == alleleNum {
				totalAltCount++
				totalGtCount++
				homs = append(homs, header[i])
				continue
			}

			continue
		}

		if len(fields[i]) < 3 {
			log.Printf("%s:%s : Skipping. Couldn't decode genotype %s (expected at least 3 characters)", fields[chromIdx], fields[posIdx], fields[i])
			return nil, nil, nil, 0, 0
		}

		// If for some reason the first allele call isn't . but 2nd is,
		// send that to missing too
		// Important to check here, because even if !simpleGt, GATK
		// will not output genotype QC data for missing genotypes
		// i.e, even with a format string (0/0:1,2:3:4:etc), we will see './.'
		if fields[i][2] == '.' {
			missing = append(missing, header[i])
			continue
		}

		// Allow for some rare cases where > 10 alleles (including reference)
		if fields[i][1] == '|' || fields[i][2] == '|' {
			// Speed up the most common cases
			if simpleGt {
				if fields[i] == "0|0" {
					totalGtCount += 2
					continue
				}

				// No longer strictly needed because of line 863
				// if fields[i] == ".|." {
				//   missing = append(missing, header[i])
				//   continue
				// }

				if alleleNum == '1' {
					if fields[i] == "0|1" || fields[i] == "1|0" {
						totalGtCount += 2
						totalAltCount++
						hets = append(hets, header[i])
						continue
					}

					if fields[i] == "1|1" {
						totalGtCount += 2
						totalAltCount += 2
						homs = append(homs, header[i])
						continue
					}
				}
			} else {
				if fields[i][0:4] == "0|0:" {
					totalGtCount += 2
					continue
				}

				if alleleNum == '1' {
					if fields[i][0:4] == "0|1:" || fields[i][0:4] == "1|0:" {
						totalGtCount += 2
						totalAltCount++
						hets = append(hets, header[i])
						continue
					}

					if fields[i][0:4] == "1|1:" {
						totalGtCount += 2
						totalAltCount += 2
						homs = append(homs, header[i])
						continue
					}
				}
			}

			gt = strings.Split(fields[i], "|")
		} else {
			// alleles separated by /, or some very malformed file
			if simpleGt {
				if fields[i] == "0/0" {
					totalGtCount += 2
					continue
				}

				if alleleNum == '1' {
					if fields[i] == "0/1" || fields[i] == "1/0" {
						totalGtCount += 2
						totalAltCount++
						hets = append(hets, header[i])
						continue
					}

					if fields[i] == "1/1" {
						totalGtCount += 2
						totalAltCount += 2
						homs = append(homs, header[i])
						continue
					}
				}
			} else {
				if len(fields[i]) < 4 {
					log.Printf("%s:%s : Skipping. Couldn't decode genotype %s (expected FORMAT data)", fields[chromIdx], fields[posIdx], fields[i])
					return nil, nil, nil, 0, 0
				}

				if fields[i][0:4] == "0/0:" {
					totalGtCount += 2
					continue
				}

				if alleleNum == '1' {
					if fields[i][0:4] == "0/1:" || fields[i][0:4] == "1/0:" {
						totalGtCount += 2
						totalAltCount++
						hets = append(hets, header[i])
						continue
					}

					if fields[i] == "1/1:" {
						totalGtCount += 2
						totalAltCount += 2
						homs = append(homs, header[i])
						continue
					}
				}
			}

			gt = strings.Split(fields[i], "/")
		}

		//https://play.golang.org/p/zjUf2rhBHn
		if len(gt) == 1 {
			log.Printf("%s:%s : Skipping. Couldn't decode genotype %s", fields[chromIdx], fields[posIdx], fields[i])
			return nil, nil, nil, 0, 0
		}

		// We should only get here for triploid+ and multiallelics
		altCount = 0
		gtCount = 0
		// log.Print(gt)
		for _, val := range gt {
			if val[0] == '.' {
				missing = append(missing, header[i])
				continue SAMPLES
			}

			if rune(val[0]) == alleleNum {
				altCount++
			}

			gtCount++
		}

		totalGtCount += gtCount
		totalAltCount += altCount

		if altCount == 0 {
			continue
		}

		if altCount == gtCount {
			homs = append(homs, header[i])
		} else {
			hets = append(hets, header[i])
		}
	}

	return homs, hets, missing, totalAltCount, totalGtCount
}
