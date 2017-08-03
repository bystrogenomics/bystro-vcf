// estab exports elasticsearch fields as tab separated values
package main
import (
  "bufio"
  "bytes"
  "flag"
  "fmt"
  "io"
  "log"
  "os"
  "strconv"
  "strings"
  "sync"
  "regexp"
  "runtime/pprof"
  "github.com/akotlar/sequtils/parse"
)

const concurrency int = 8
const chromIdx int = 0
const posIdx int = 1
const idIdx int = 2
const refIdx int = 3
const altIdx int = 4
const qualIdx int = 5
const filterIdx int = 6
const infoIdx int = 7
const formatIdx int = 8

type Config struct {
  inPath string
  errPath string
  emptyField string
  fieldDelimiter string
  keepId bool
  keepInfo bool
  cpuProfile string
  keepFiltered map[string]bool
}

func setup(args []string) *Config {
  config := &Config{}
  flag.StringVar(&config.inPath, "inPath", "", "The input file path (optional: default is stdin)")
  flag.StringVar(&config.errPath, "errPath", "", "The output path for the JSON output (optional)")
  flag.StringVar(&config.emptyField, "emptyField", "!", "The output path for the JSON output (optional)")
  flag.StringVar(&config.fieldDelimiter, "fieldDelimiter", ";", "The output path for the JSON output (optional)")
  flag.BoolVar(&config.keepId, "keepId", false, "Retain the ID field in output")
  flag.BoolVar(&config.keepInfo, "keepInfo", false, "Retain INFO field in output (2 appended output fields: allele index and the INFO field. Will appear after id field if --keepId flag set.")
  flag.StringVar(&config.cpuProfile, "cpuProfile", "", "Write cpu profile to file at this path")
  filteredVals := flag.String("keepFilter", "PASS,.", "Allow rows that have this FILTER value (comma separated)")
  excludeFilterVals := flag.String("excludeFilter", "", "Exclude rows that have this FILTER value (comma separated)")
  // allows args to be mocked https://github.com/nwjlyons/email/blob/master/inputs.go
  // can only run 1 such test, else, redefined flags error
  a := os.Args[1:]
  if args != nil {
    a = args
  }
  flag.CommandLine.Parse(a)
  config.keepFiltered = map[string]bool{"PASS": true, ".": true}
  if *filteredVals != "" {
    for _, val := range strings.Split(*filteredVals, ",") {
      config.keepFiltered[val] = true
    }
  }
  if *excludeFilterVals != "" {
    for _, val := range strings.Split(*excludeFilterVals, ",") {
      config.keepFiltered[val] = false
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

  if config.cpuProfile != "" {
    f, err := os.Create(config.cpuProfile)
    if err != nil {
        log.Fatal(err)
    }
    pprof.StartCPUProfile(f)
    defer pprof.StopCPUProfile()
  }

  reader := bufio.NewReader(inFh)

  readVcf(config, reader, func(row string) {fmt.Print(row)})
}

func readVcf (config *Config, reader *bufio.Reader, resultFunc func(row string)) {
  foundHeader := false

  var header []string

  // Read buffer
  workQueue := make(chan string, 100)
  complete := make(chan bool)
  // Write buffer
  results := make(chan string, 100)
  var wg sync.WaitGroup

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
    record := strings.Split(row[:len(row) - numChars], "\t")

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

  // Read the lines into the work queue.
  go func() {
    for {
      row, err := reader.ReadString(endOfLineByte) // 0x0A separator = newline

      if err == io.EOF {
        break
      } else if err != nil {
        log.Fatal(err)
      } else if row == "" {
        // We may have not closed the pipe, but not have any more information to send
        // Wait for EOF
        continue
      }

      workQueue <- row[:len(row) - numChars];
    }

    // Close the channel so everyone reading from it knows we're done.
    close(workQueue)
  }()

  wg.Add(1)
  go func() {
    defer wg.Done()
    for line := range results {
      resultFunc(line)
    }
  }()

  // Now read them all off, concurrently.
  for i := 0; i < concurrency; i++ {
    go processLines(header, config.emptyField, config.fieldDelimiter, config.keepId,
      config.keepInfo, config.keepFiltered, workQueue, results, complete)
  }

  // Wait for everyone to finish.
  for i := 0; i < concurrency; i++ {
    <-complete
  }

  close(results)

  wg.Wait()
}

func linePasses(record []string, header []string, filterKeys map[string]bool) bool {
  return len(record) == len(header) && len(filterKeys) == 0 || filterKeys[record[filterIdx]] == true
}

// func gatherAlt(ref byte, alt string, pos string, multiallelic bool, altCache map[byte]map[string][]string) []string {
//  if altCache[ref] == nil {
//    alt[ref] = make(map[string][]string
//  } else if altCache[ref][alt] != nil {
//    return 
//  }

// }

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

func processLines(header []string, emptyField string, fieldDelimiter string, keepId bool, keepInfo bool,
keepFiltered map[string]bool, queue chan string, results chan string, complete chan bool) {
  var alleles []string
  var multiallelic bool

  for line := range queue {
    record := strings.Split(line, "\t")
    
    if !linePasses(record, header, keepFiltered) {
      continue
    }

    var homs []string
    var hets []string
    var missing []string

    alleles = strings.Split(record[altIdx], ",")
    multiallelic = len(alleles) > 1

    for idx, allele := range alleles {
      if altIsValid(allele) == false {
        log.Printf("%s:%s Skip ALT #%d (not ACTG)", record[chromIdx], record[posIdx], idx+1)
        continue
      }

      siteType, pos, ref, alt, err := updateFieldsWithAlt(record[refIdx], allele, record[posIdx], multiallelic)
      if err != nil {
        log.Fatal(err)
      }

      if pos == "" {
        log.Printf("%s:%s Skip ALT #%d (complex)", record[chromIdx], record[posIdx], idx+1)
        continue
      }

      // If no sampels are provided, annotate what we can, skipping hets and homs
      if len(header) > 9 {
        homs, hets, missing = makeHetHomozygotes(record, header, strconv.Itoa(idx+1))

        if len(homs) == 0 && len(hets) == 0 {
          continue
        }
      }

      var output bytes.Buffer
      if len(record[chromIdx]) < 4 || record[chromIdx][0:2] != "ch" {
        output.WriteString("chr")
      }
      output.WriteString(record[chromIdx])
      output.WriteString("\t")
      output.WriteString(pos)
      output.WriteString("\t")
      if multiallelic {
        output.WriteString("MULTIALLELIC")
      } else {
        output.WriteString(siteType)
      }

      output.WriteString("\t")
      output.WriteString(ref)
      output.WriteString("\t")
      output.WriteString(alt)
      output.WriteString("\t")

      if(multiallelic) {
        output.WriteString(parse.NotTrTv)
      } else {
        output.WriteString(parse.GetTrTv(ref, alt))
      }

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

      output.WriteString("\t")

      if len(missing) == 0 {
        output.WriteString(emptyField)
      } else {
        output.WriteString(strings.Join(missing, fieldDelimiter))
      }

      if keepId == true {
        output.WriteString("\t")
        output.WriteString(record[idIdx])
      }

      if keepInfo == true {
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

  // log.Println("Worker hit, missed this many times: ", hitCount, missCount)
  // Let the main process know we're done.
  complete <- true
}

func updateFieldsWithAlt(ref string, alt string, pos string, multiallelic bool) (string, string, string, string, error) {
  /*********************** SNPs *********************/
  if len(alt) == len(ref) {
    if alt == ref {
      // No point in returning ref sites
      return "", "", "", "", nil
    }

    if len(ref) > 1 {
      // SNPs that are multiallelic with indels, can be longer than 1 base long
      // Confusingly enough, so can MNPs
      // So, we will check for both
      if ref[0] != alt[0] && ref[1] != alt[1]{
        // This is most likely an MNP
        // As MNPs are contiguous
        // Currently we haven't enabled MNP parsing in the caller
        return "", "", "", "", nil
      }

      diffIdx := -1
      // Let's check each base; if there is more than 1 change, that is an error
      for i := 0; i < len(ref); i++ {
        if ref[i] != alt[i] {
          // Major red flag, there should never be a len(ref) == len(alt) allele that isn't an MNP or SNP
          // TODO: should we relax this? If we allow MNP, may as well allow sparse MNPs
          if diffIdx > -1 {
            return "", "", "", "", nil
          }

          diffIdx = i
        }
      }

      // Most cases are diffIdx == 0, allow us to skip 1 strconv.Atoi, 1 assignment, 1 strconv.Itoa, and 1 addition
      if diffIdx == 0 {
        return "SNP", pos, string(ref[diffIdx]), string(alt[diffIdx]), nil
      }

      intPos, _ := strconv.Atoi(pos)
      return "SNP", strconv.Itoa(intPos + diffIdx), string(ref[diffIdx]), string(alt[diffIdx]), nil
    }

    return "SNP", pos, ref, alt, nil
  }

  /*********************** INSERTIONS AND DELETIONS *********************/
  // TODO: Handle case where first base of contig is deleted, and padded as the first unmodified base downstream
  //First base is always padding
  if ref[0] != alt[0] {
    return "", "", "", "", nil
  }

  /*************************** DELETIONS FIRST **************************/
  if len(ref) > len(alt) {
    intPos, err := strconv.Atoi(pos)
    if err != nil {
      return "", "", "", "", err
    }

    // Simple insertions delete the entire reference, sans the padding base to the left
    // TODO: handle 1st base deleted in contig, padded to right
    if len(alt) == 1 {
      return "DEL", strconv.Itoa(intPos + 1), string(ref[1]), strconv.Itoa(1 - len(ref)), nil
    }

    // log.Println("Complex del", pos, ref, alt, multiallelic)

    // Complex deletions, inside of a reference
    // Ex: Ref: TCT Alt: T, TT (the TT is a 1 base C deletion)
    //this typically only comes up with multiallelics that have a 2nd deletion, that covers all bases (excepting 1 padding base) in the reference
    //In other cases, just skip for now, mostly seems like an error
    if multiallelic == false {
      return "", "", "", "", nil
    }

    // Our deletion should happen within the reference, so the non-padded
    // portion of the reference is what we'll check
    if strings.Contains(ref, alt[1: ]) == false {
      return "", "", "", "", nil
    }
    // TODO: More precise checking; for instance we can check if the alt is contained within the end of the ref (sans the 1 base deletion)
    return "DEL", strconv.Itoa(intPos + 1), string(ref[1]), strconv.Itoa(len(alt) - len(ref)), nil
  }

  /*********************** INSERTIONS *********************/
  // len(ref) > 1 should always indicate that this is a multiallelic that contains a deletion
  // therefore requiring 1 base of padding to the left
  // there may be cases where VCF variants are unnecessarily padded, but lets skip these
  if len(ref) > 1 {
    if multiallelic == false {
      return "", "", "", "", nil
    }

    // log.Println("Complex ins", pos, ref, alt, multiallelic)

    // Our insertion should happen within the reference, so the non-padded
    // portion of the reference is what we'll check
    if strings.Contains(alt, ref[1: ]) == false {
      return "", "", "", "", nil
    }
    // TODO: More precise checking; for instance we can check if the alt is contained within the end of the ref (sans the 1 base deletion)
    var insBuffer bytes.Buffer
    insBuffer.WriteString("+")
    insBuffer.WriteString(alt[1:len(alt) - len(ref) + 1])
    return "INS", pos, string(ref[0]), insBuffer.String(), nil
  }

  var insBuffer bytes.Buffer
  insBuffer.WriteString("+")
  insBuffer.WriteString(alt[1:])
  return "INS", pos, ref, insBuffer.String(), nil
}

func makeHetHomozygotes(fields []string, header []string, alleleIdx string) ([]string, []string, []string) {
  simpleGt := strings.Contains(fields[formatIdx], ":")

  gt := make([]string, 0, 2)
  gtCount := 0
  altCount := 0

  var homs []string
  var hets []string
  var missing []string

  SAMPLES:
    for i := 9; i < len(header); i++ {
      if strings.Contains(fields[i], "|") {
        if (simpleGt && strings.Contains(fields[i], "0|0:")) || fields[i] == "0|0" {
          continue
        }

        if strings.Contains(fields[i], ".|.") {
          missing = append(missing, header[i])
          continue SAMPLES
        }

        gt = strings.Split(fields[i], "|")
      } else {
        if (simpleGt && strings.Contains(fields[i], "0/0:")) || fields[i] == "0/0" {
          continue
        }

        if strings.Contains(fields[i], "./.") {
          missing = append(missing, header[i])
          continue SAMPLES
        }

        gt = strings.Split(fields[i], "/")
      }

      altCount = 0
      gtCount = 0
      for _, val := range gt {
        if val == "." {
          missing = append(missing, header[i])
          continue SAMPLES
        }

        // val[0] doesn't work...because byte array?
        if string(val[0]) == alleleIdx {
          altCount++
        }

        gtCount++
      }

      if altCount == 0 {
        continue
      }

      if altCount == gtCount {
        homs = append(homs, header[i])
      } else {
        hets = append(hets, header[i])
      }
    }

  return homs, hets, missing
}