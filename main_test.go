package main

import (
	"strings"
	"testing"
	"strconv"
	"fmt"
	"bufio"
	"github.com/akotlar/sequtils/parse"
)

func TestKeepFlagsTrue(t *testing.T) {
	args := []string{
		"--keepInfo",
		"--keepId",
		"--inPath", "/path/to/file",
		"--errPath", "/path/to/err",
		"--cpuProfile", "/path/to/profile",
		"--emptyField", ".",
		"--fieldDelimiter", "&",
	}

	config := setup(args)

  if config.keepInfo != true || config.keepId != true {
  	t.Error("NOT OK: parse keepInfo and keepId args")
  }

  if config.inPath != "/path/to/file" || config.errPath != "/path/to/err" ||
  config.cpuProfile != "/path/to/profile" {
  	t.Error("NOT OK: parse inPath and errPath args")
  }

  if config.emptyField != "." || config.fieldDelimiter != "&" {
  	t.Error("NOT OK: parse emptyField and fieldDelimiter args")
  }
}

func TestUpdateFieldsWithAlt(t *testing.T) {
	expType, exPos, expRef, expAlt := "SNP", "100", "T", "C"

	sType, pos, ref, alt, err := updateFieldsWithAlt("T", "C", "100", false)

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed")
	} else {
		t.Log("OK: SNP")
	}

	expType, exPos, expRef, expAlt = "SNP", "103", "T", "A"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TCCT", "TCCA", "100", false)

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: SNPs that are longer than 1 base are suported when SNP at end")
	}

	expType, exPos, expRef, expAlt = "SNP", "102", "C", "A"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TGCT", "TGAT", "100", false)

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: SNPs that are longer than 1 base are suported when SNP in middle")
	}

	expType, exPos, expRef, expAlt = "SNP", "100", "T", "A"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TGCT", "AGCT", "100", false)

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: SNPs that are longer than 1 base are suported when SNP at beginning")
	}

	expType, exPos, expRef, expAlt = "", "", "", ""

	sType, pos, ref, alt, err = updateFieldsWithAlt("TCCT", "GTAA", "100", false)

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: MNPs are not supported")
	}

	expType, exPos, expRef, expAlt = "DEL", "101", "C", "-1"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TC", "T", "100", false)

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: 1-based deletions ")
	}

	expType, exPos, expRef, expAlt = "DEL", "101", "A", "-5"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TAGCGT", "T", "100", false)

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Deletions with references longer than 2 bases")
	}

	// Test multiallelic intercolated deletion
	expType, exPos, expRef, expAlt = "DEL", "101", "A", "-4"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TAGCTT", "TA", "100", true)

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Deletions longer than 1 base")
	}

	// If not multiallelic, the same deletion should be skipped, as it is odd,
	// and the position doesn't make sense
	expType, exPos, expRef, expAlt = "", "", "", ""

	sType, pos, ref, alt, err = updateFieldsWithAlt("TAGCTT", "TA", "100", false)

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Deletions longer than 1 base")
	}

	expType, exPos, expRef, expAlt = "", "", "", ""

	sType, pos, ref, alt, err = updateFieldsWithAlt("TAGCTT", "TAC", "100", false)

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Left edge of deletion must match exactly")
	}

	expType, exPos, expRef, expAlt = "INS", "100", "T", "+AGCTT"

	sType, pos, ref, alt, err = updateFieldsWithAlt("T", "TAGCTT", "100", false)

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Insertions where reference is 1 base long")
	}

	// Test multiallelic
	// say came from TT T, TAGCTT
	// I believe the answer is, this should be a +AGCT in between then two T's
	expType, exPos, expRef, expAlt = "INS", "100", "T", "+AGCT"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TT", "TAGCTT", "100", true)

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Multiallelics insertion where reference is 2 bases long")
	}

	expType, exPos, expRef, expAlt = "", "", "", ""

	sType, pos, ref, alt, err = updateFieldsWithAlt("TT", "TAGCTT", "100", false)

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Non-multiallelic insertions, whose ref are 2 bases long are skipped")
	}
}

func TestPassesLine(t *testing.T) {
	expect := true

	header := []string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"}
	//Test example 5.2.4 https://samtools.github.io/hts-specs/VCFv4.2.pdf
	record := []string{"20", "4", ".", "GCG", "G,GCGCG", ".", "PASS", "DP=100"}

	keepFiltered := map[string]bool{ "PASS": true, ".": true}
	actual := linePasses(record, header, keepFiltered)

	if actual == expect {
		t.Log("OK: PASS lines pass")
	} else {
		t.Error("NOT OK: PASS lines should pass")
	}

	record = []string{"20", "4", ".", "GCG", "G,GCGCG", ".", ".", "DP=100"}

	actual = linePasses(record, header, keepFiltered)

	if actual == expect {
		t.Log("OK: lines with missing (.) values under FILTER pass")
	} else {
		t.Error("NOT OK: lines with missing (.) values under FILTER pass")
	}

	expect = false
	actual = altIsValid(".")

	if expect != actual {
		t.Error("NOT OK: Can't handle missing Alt alleles")
	} else {
		t.Log("OK: Handles missing Alt alleles")
	}

	expect = false
	actual = altIsValid("]13 : 123456]T")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("NOT OK: Can't handle single breakends")
	} else {
		t.Log("OK: Handles single breakend ']13 : 123456]T'")
	}

	expect = false
	actual = altIsValid("C[2 : 321682[")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("NOT OK: Can't handle single breakends")
	} else {
		t.Log("OK: Handles single breakend 'C[2 : 321682['")
	}

	expect = false
	actual = altIsValid(".A")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("NOT OK: Can't handle single breakends")
	} else {
		t.Log("OK: Handles single breakend '.A'")
	}

	expect = false
	actual = altIsValid("G.")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("NOT OK: Can't handle single breakends")
	} else {
		t.Log("OK: Handles single breakend 'G.'")
	}

	expect = false
	actual = altIsValid("<DUP>")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("NOT OK: Can't handle complex tags")
	} else {
		t.Log("OK: Handles complex Alt tags '<DUP>'")
	}

	expect = false
	actual = altIsValid("A,C")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("NOT OK: Allows multiallelics", actual)
	} else {
		t.Log("OK: multiallelics are not supported. altIsValid() requires multiallelics to be split")
	}
}

func TestAltIsValid(t *testing.T) {
	expect := true

	actual := altIsValid("ACTG")

	if expect != actual {
		t.Error()
	} else {
		t.Log("NOT OK: Support ACTG-containing alleles")
	}

	expect = false
	actual = altIsValid(".")

	if expect != actual {
		t.Error("NOT OK: Can't handle missing Alt alleles")
	} else {
		t.Log("OK: Handles missing Alt alleles")
	}

	expect = false
	actual = altIsValid("]13 : 123456]T")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("NOT OK: Can't handle single breakends")
	} else {
		t.Log("OK: Handles single breakend ']13 : 123456]T'")
	}

	expect = false
	actual = altIsValid("C[2 : 321682[")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("NOT OK: Can't handle single breakends")
	} else {
		t.Log("OK: Handles single breakend 'C[2 : 321682['")
	}

	expect = false
	actual = altIsValid(".A")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("NOT OK: Can't handle single breakends")
	} else {
		t.Log("OK: Handles single breakend '.A'")
	}

	expect = false
	actual = altIsValid("G.")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("NOT OK: Can't handle single breakends")
	} else {
		t.Log("OK: Handles single breakend 'G.'")
	}

	expect = false
	actual = altIsValid("<DUP>")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("NOT OK: Can't handle complex tags")
	} else {
		t.Log("OK: Handles complex Alt tags '<DUP>'")
	}

	expect = false
	actual = altIsValid("A,C")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("NOT OK: Should require single alleles in altIsValid", actual)
	} else {
		t.Log("OK: multiallelics are not supported. altIsValid() requires multiallelics to be split")
	}
}

func TestMakeHetHomozygotes(t *testing.T) {
	header := []string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1", "S2", "S3", "S4"}

	sharedFieldsGT := []string{"10", "1000", "rs#", "C", "T", "100", "PASS", "AC=1", "GT"}

	// GT:DS:GL is what 1000 genomes phase 1 provides
	sharedFieldsGTcomplex := []string{"10", "1000", "rs#", "C", "T", "100", "PASS", "AC=1", "GT:DS:GL"}

	fields := append(sharedFieldsGT, "0|0", "0|0", "0|0", "0|0")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing := makeHetHomozygotes(fields, header, "1")

	if len(actualHoms) == 0 && len(actualHets) == 0 && len(missing) == 0 {
		t.Log("OK: Homozygous reference samples are skipped")
	} else {
		t.Error("NOT OK: 0 alleles give unexpected results", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, "0|1", "0|1", "0|1", "0|1")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing = makeHetHomozygotes(fields, header, "1")

	if len(actualHoms) == 0 && len(actualHets) == 4 && len(missing) == 0 {
		t.Log("OK: handles hets")
	} else {
		t.Error("NOT OK: 0 alleles give unexpected results", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, ".|1", "0|1", "0|1", "0|1")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing = makeHetHomozygotes(fields, header, "1")

	if len(actualHoms) == 0 && len(actualHets) == 3 && len(missing) == 1 {
		t.Log("OK: GT's containing missing data are entirely uncertain, therefore skipped. Missing genotypes are called if any of the calls are missing")
	} else {
		t.Error("NOT OK: Fails to handle missing data", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, "1|.", "0|1", "0|1", "0|1")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing = makeHetHomozygotes(fields, header, "1")

	if len(actualHoms) == 0 && len(actualHets) == 3 && len(missing) == 1 {
		t.Log("OK: GT's containing missing data are entirely uncertain, therefore skipped. Missing genotypes are called if any of the calls are missing")
	} else {
		t.Error("NOT OK: Fails to handle missing data", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, "1|1", "1|1", "0|1", "0|1")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing = makeHetHomozygotes(fields, header, "1")

	if len(actualHoms) == 2 && len(actualHets) == 2 {
		t.Log("OK: handles homs and hets")
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, "1|2", "1|1", "0|1", "0|1")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing = makeHetHomozygotes(fields, header, "1")

	if len(actualHoms) == 1 && len(actualHets) == 3 {
		t.Log("OK: a sample heterozygous for a wanted allele is heterozygous for that allele even if its other allele is unwanted (for multiallelic phasing)")
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, "1|2", "1|1", "0|1", "0|1")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing = makeHetHomozygotes(fields, header, "2")

	if len(actualHoms) == 0 && len(actualHets) == 1 {
		t.Log("OK: Het / homozygous status is based purely on the wanted allele, rather than total non-ref count")
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGTcomplex, "1|2:-0.03,-1.12,-5.00", "1|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing = makeHetHomozygotes(fields, header, "2")

	if len(actualHoms) == 0 && len(actualHets) == 1 {
		t.Log("OK: handles complicated GTs, with non-1 alleles", fields)
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGTcomplex, "1|2|1:-0.03,-1.12,-5.00", "1|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing = makeHetHomozygotes(fields, header, "2")

	if len(actualHoms) == 0 && len(actualHets) == 1 {
		t.Log("OK: Complicated GT: Triploids are considered het if only 1 present desired allele", fields)
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, "1|2|1", "1|1", "0|1", "0|1")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing = makeHetHomozygotes(fields, header, "2")

	if len(actualHoms) == 0 && len(actualHets) == 1 {
		t.Log("OK: Triploids are considered het if only 1 present desired allele")
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGTcomplex, "2|2|2:-0.03,-1.12,-5.00", "1|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing = makeHetHomozygotes(fields, header, "2")

	if len(actualHoms) == 1 && len(actualHets) == 0 {
		t.Log("OK: Complicated GT: Triploids are considered hom only if all alleles present are desired", fields)
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, "2|2|2", "1|1", "0|1", "0|1")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing = makeHetHomozygotes(fields, header, "2")

	if len(actualHoms) == 1 && len(actualHets) == 0 {
		t.Log("OK: Triploids are considered hom only if all alleles present are desired")
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}
}

func TestNomralizeSampleNames(t *testing.T) {
	header := []string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1.HAHAHAH", "S2.TRYINGTO.MESSYOUUP", "S3", "S-4"}

	parse.NormalizeHeader(header)

	for i := 9; i < len(header); i++ {
		if strings.Contains(header[i], ".") {
			t.Error("NOT OK: Couldn't replace period")
		} else {
			t.Log("OK: no periods found in", header[i])
		}
	}

	if header[9] == "S1_HAHAHAH" {
		t.Log("OK: replaced period in S1.HAHAHAH", header[9])
	} else {
		t.Error("NOT OK: Couldn't replace period in S1.HAHAHAH", header[9])
	}

	if header[10] == "S2_TRYINGTO_MESSYOUUP" {
		t.Log("OK: replaced two periods in S2.TRYINGTO.MESSYOUUP", header[10])
	} else {
		t.Error("NOT OK: Couldn't replace periods in S2.TRYINGTO.MESSYOUUP", header[10])
	}

	if header[11] == "S3" {
		t.Log("OK: didn't mess up name S3", header[11])
	} else {
		t.Error("NOT OK: Messed up name S3", header[11])
	}

	if header[12] == "S-4" {
		t.Log("OK:  didn't mess up name without a period", header[12])
	} else {
		t.Error("NOT OK: Messed up name S-4", header[12])
	}
}

func TestOutputsInfo(t *testing.T) {
	versionLine := "##fileformat=VCFv4.x"
	header := strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"}, "\t")
	record := strings.Join([]string{"10", "1000", "rs#", "C", "T", "100", "PASS", "AC=1", "GT"}, "\t")

	lines := versionLine + "\n" + header + "\n" + record	+ "\n"
	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", keepId: false, keepInfo: true}

	readVcf(&config, reader, func(row string) {
  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr10" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}
  	
  	if len(resultRow) == 11 {
  		t.Log("OK: With keepInfo flag set, but not keepId, should output 11 fields")
  	} else {
  		t.Error("NOT OK: With keepInfo flag set, but not keepId, should output 11 fields", resultRow)
  	}
  	
		if resultRow[len(resultRow) - 2] == "0" && resultRow[len(resultRow) - 1] == "AC=1" {
			t.Log("OK: add INFO field correctly for single field")
		} else {
			t.Error("NOT OK: Couldn't add INFO field", resultRow)
		}
	})

	record = strings.Join([]string{"10", "1000", "rs#", "C", "T,G", "100", "PASS", "AC=1", "GT"}, "\t")

	lines = versionLine + "\n" + header + "\n" + record	+ "\n"
	reader = bufio.NewReader(strings.NewReader(lines))

	index := -1;
  readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")
  	
  	if resultRow[0] != "chr10" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if len(resultRow) == 11 {
  		t.Log("OK: With keepInfo flag set, but not keepId, should output 11 fields")
  	} else {
  		t.Error("NOT OK: With keepInfo flag set, but not keepId, should output 11 fields", resultRow)
  	}

  	altIdx, err := strconv.Atoi(resultRow[len(resultRow) - 2])

  	if err != nil {
  		t.Error("NOT OK: The 8th column should be numeric")
  	}

  	if altIdx == index && resultRow[len(resultRow) - 1] == "AC=1" {
			t.Log("OK: add INFO field correctly for multiple field, index", altIdx)
		} else {
			t.Error("NOT OK: Couldn't add INFO field", resultRow)
		}
	})

	if index != 1 {
		t.Error("NOT OK: Expected to parse 2 alleles, parsed fewer")
	}
}

func TestOutputsId(t *testing.T) {
	versionLine := "##fileformat=VCFv4.x"
	header := strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"}, "\t")
	record := strings.Join([]string{"10", "1000", "rs123", "C", "T", "100", "PASS", "AC=1", "GT"}, "\t")

	lines := versionLine + "\n" + header + "\n" + record	+ "\n"

	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", keepId: true, keepInfo: false}

  readVcf(&config, reader, func(row string) {
  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr10" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if len(resultRow) == 10 {
  		t.Log("OK: With keepId flag set, but not keepInfo, should output 10 fields")
  	} else {
  		t.Error("NOT OK: With keepId flag set, but not keepInfo, should output 10 fields", resultRow)
  	}

		if resultRow[len(resultRow) - 1] == "rs123" {
			t.Log("OK: add ID field correctly for single field")
		} else {
			t.Error("NOT OK: Couldn't add ID field", resultRow)
		}
	})

	record = strings.Join([]string{"10", "1000", "rs456", "C", "T,G", "100", "PASS", "AC=1", "GT"}, "\t")

	lines = versionLine + "\n" + header + "\n" + record	+ "\n"
	reader = bufio.NewReader(strings.NewReader(lines))

	index := -1
  readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr10" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if len(resultRow) == 10 {
  		t.Log("OK: With keepId flag set, but not keepInfo, should output 10 fields")
  	} else {
  		t.Error("NOT OK: With keepId flag set, but not keepInfo, should output 10 fields", resultRow)
  	}

  	if resultRow[len(resultRow) - 1] == "rs456" {
			t.Log("OK: add ID field correctly for multiple field, index", altIdx)
		} else {
			t.Error("NOT OK: Couldn't add ID field", resultRow)
		}
	})

	if index != 1 {
		t.Error("NOT OK: Expected to parse 2 alleles, parsed fewer")
	}
}

func TestOutputsSamplesIdAndInfo(t *testing.T) {
	versionLine := "##fileformat=VCFv4.x"
	header := strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
		"FILTER", "INFO", "FORMAT", "Sample1", "Sample2", "Sample3", "Sample4"}, "\t")
	record := strings.Join([]string{"10", "1000", "rs123", "C", "T", "100", "PASS",
		"AC=1", "GT", "0/0", "0/1", "1/1", "./."}, "\t")

	lines := versionLine + "\n" + header + "\n" + record	+ "\n"

	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", keepId: true, keepInfo: true}

  readVcf(&config, reader, func(row string) {
  	resultRow := strings.Split(row[:len(row)-1], "\t")
  	
  	if resultRow[0] != "chr10" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if len(resultRow) == 12 {
  		t.Log("OK: With both keepId and retainInfo flags set, should output 12 fields")
  	} else {
  		t.Error("NOT OK: With both keepId and retainInfo flags set, should output 12 fields", resultRow)
  	}

		if resultRow[len(resultRow) - 3] == "rs123" {
			t.Log("OK: add ID field correctly for single field")
		} else {
			t.Error("NOT OK: Couldn't add ID field", resultRow)
		}

		if resultRow[len(resultRow) - 2] == "0" && resultRow[len(resultRow) - 1] == "AC=1" {
			t.Log("OK: add INFO field correctly for single field")
		} else {
			t.Error("NOT OK: Couldn't add INFO field", resultRow)
		}

		if resultRow[6] == "Sample2" {
			t.Log("OK: Recapitualte the het", resultRow)
		} else {
			t.Error("NOT OK: Couldn't recapitualte the het", resultRow)
		}

		if resultRow[7] == "Sample3" {
			t.Log("OK: Recapitualte the homozygote", resultRow)
		} else {
			t.Error("NOT OK: Couldn't recapitualte the homozygote", resultRow)
		}

		if resultRow[8] == "Sample4" {
			t.Log("OK: Recapitualte the missing sample", resultRow)
		} else {
			t.Error("NOT OK: Couldn't recapitualte the missing sample", resultRow)
		}
	})

	record = strings.Join([]string{"10", "1000", "rs456", "C", "T,G", "100", "PASS",
		"AC=1", "GT", "1|1", "0|0", "0|2", ".|."}, "\t")

	lines = versionLine + "\n" + header + "\n" + record	+ "\n"
	reader = bufio.NewReader(strings.NewReader(lines))

	index := -1
  readVcf(&config, reader, func(row string) {
  	index++
 
  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr10" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if len(resultRow) == 12 {
  		t.Log("OK: With both keepId and retainInfo flags set, should output 12 fields")
  	} else {
  		t.Error("With both keepId and retainInfo flags set, should output 12 fields", resultRow)
  	}

  	if resultRow[len(resultRow) - 3] == "rs456" {
			t.Log("OK: add ID field correctly for multiple field, index", altIdx)
		} else {
			t.Error("NOT OK: Couldn't add ID field", altIdx, resultRow)
		}

		altIdx, err := strconv.Atoi(resultRow[len(resultRow) - 2])

  	if err != nil {
  		t.Error("NOT OK: The 9th column should be numeric")
  	}

  	if altIdx == index {
  		t.Log("OK: Multiallelic index is in 10th column when keepInfo is true")
		} else {
			t.Error("NOT OK: Multiallelic index isn't in 10th column when keepInfo is true", resultRow)
		}

		if resultRow[len(resultRow) - 1] == "AC=1" {
			t.Log("OK: add INFO field correctly for multiallelic field in column 10")
		} else {
			t.Error("NOT OK: Couldn't add INFO field", resultRow)
		}
	})

	if index != 1 {
		t.Error("NOT OK: Expected to parse 2 alleles, parsed fewer")
	}
}

func TestOutputMultiallelic(t *testing.T) {
	versionLine := "##fileformat=VCFv4.x"
	header := strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
	"Format", "Sample1", "Sample2", "Sample3", "Sample4"}, "\t")
	//Test example 5.2.4 https://samtools.github.io/hts-specs/VCFv4.2.pdf
	//Except use A as last base, to properly distinguish first and last base of reference
	record := strings.Join([]string{"20", "4", ".", "GCACG", "G,GTCACACG", ".", "PASS", "DP=100",
	"GT", "0|0", "0|1", "2|2", ".|."}, "\t")

	keepFiltered := map[string]bool{ "PASS": true, ".": true}

	lines := versionLine + "\n" + header + "\n" + record	+ "\n"
	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", keepFiltered: keepFiltered}

	index := -1
  readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr20" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if index == 0 {
  		if resultRow[altIdx] == "-4" && resultRow[posIdx] == "5" && resultRow[refIdx] == "C" {
  			t.Log("OK: 2 base deletion recapitulated", resultRow)
  		} else {
  			t.Error("NOT OK: 2 base deletion not recapitulated", resultRow)
  		}
			
			if resultRow[6] == "Sample2" {
				t.Log("OK: Recapitualte 1st allele het", resultRow)
  		} else {
  			t.Error("NOT OK: Couldn't recapitualte 1st allele het", resultRow)
  		}

  		if resultRow[7] == config.emptyField {
				t.Log("OK: Recapitualte 1st allele het", resultRow)
  		} else {
  			t.Error("NOT OK: Couldn't recapitualte 1st allele het", resultRow)
  		}

  		if resultRow[8] == "Sample4" {
				t.Log("OK: Recapitualte 1st allele het", resultRow)
  		} else {
  			t.Error("NOT OK: Couldn't recapitualte 1st allele het", resultRow)
  		}
		}

		if index == 1 {
  		if resultRow[altIdx] == "+TCA" && resultRow[posIdx] == "4" && resultRow[refIdx] == "G" {
  			t.Log("OK: 2 base deletion recapitulated", resultRow)
  		} else {
  			t.Error("NOT OK: 2 base deletion not recapitulated", resultRow)
  		}

  		if resultRow[6] == config.emptyField {
				t.Log("OK: Recapitualte 2nd allele het", resultRow)
  		} else {
  			t.Error("NOT OK: Couldn't recapitualte 2nd allele het", resultRow)
  		}

  		if resultRow[7] == "Sample3" {
				t.Log("OK: Recapitualte 2nd allele het", resultRow)
  		} else {
  			t.Error("NOT OK: Couldn't recapitualte 2nd allele hom", resultRow)
  		}

  		if resultRow[8] == "Sample4" {
				t.Log("OK: Recapitualte 2nd allele het", resultRow)
  		} else {
  			t.Error("NOT OK: Couldn't recapitualte 2nd allele het", resultRow)
  		}
		}
	})

	if index == 1 {
		t.Log("OK: parsed 2 alleles")
	} else {
		t.Error("NOT OK: Expected to parse 2 alleles")
	}
}

func TestOutputComplexMultiDel(t *testing.T) {
	versionLine := "##fileformat=VCFv4.x"
	header := strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
		"FILTER", "INFO"}, "\t")
	//Test example 5.2.4 https://samtools.github.io/hts-specs/VCFv4.2.pdf
	//Test from gnomad genomes
	record := strings.Join([]string{"16", "84034434", "rs141446650", "GAGGGAGACAGAGGGAAGT",
		"G,GGGGAGACAGAGGGAAGT", ".", "PASS", "DP=100"}, "\t")

	keepFiltered := map[string]bool{ "PASS": true, ".": true}

	lines := versionLine + "\n" + header + "\n" + record	+ "\n"
	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", keepFiltered: keepFiltered}

	index := -1;

	readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr16" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if index == 0 {
  		if resultRow[refIdx] == "A" && resultRow[posIdx] == strconv.Itoa(84034434 + 1) && resultRow[altIdx] == "-18" {
  			t.Log("OK: 18 base deletion recapitulated", resultRow)
  		} else {
  			t.Error("NOT OK: 18 base deletion not recapitulated", resultRow)
  		}
		}

		if index == 1 {
  		if resultRow[refIdx] == "A" && resultRow[posIdx] == strconv.Itoa(84034434 + 1) && resultRow[altIdx] == "-1" {
  			t.Log("OK: Complex 1 base deletion recapitulated", resultRow)
  		} else {
  			t.Error("NOT OK: 1 base deletion not recapitulated", resultRow)
  		}
		}
	})

	if index == 1 {
		t.Log("OK: parsed 2 alleles")
	} else {
		t.Error("NOT OK: Expected to parse 2 alleles")
	}
}

func TestOutputComplexDel(t *testing.T) {
	versionLine := "##fileformat=VCFv4.x"
	header := strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
		"FILTER", "INFO"}, "\t")
	//Test example 5.2.4 https://samtools.github.io/hts-specs/VCFv4.2.pdf
	//Test from gnomad genomes
	record := strings.Join([]string{"1", "874816", "rs200996316",
		"CCCCCTCATCACCTCCCCAGCCACGGTGAGGACCCACCCTGGCATGATCT",
		"CCCCCTCATCACCTCCCCAGCCACGGTGAGGACCCACCCTGGCATGATCTCCCCTCATCACCTCCCCAGCCACGGTGAGGACCCACCCTGGCATGATCT,GCCCCTCATCACCTCCCCAGCCACGGTGAGGACCCACCCTGGCATGATCT,C,CTCCCCTCATCACCTCCCCAGCCACGGTGAGGACCCACCCTGGCATGATCT",
		".", "PASS", "DP=100"}, "\t")

	keepFiltered := map[string]bool{ "PASS": true, ".": true}

	lines := versionLine + "\n" + header + "\n" + record	+ "\n"
	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", keepFiltered: keepFiltered}

	index := -1;
  readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr1" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	fmt.Println(index, resultRow)
  	if index == 0 {
  		if resultRow[refIdx] == "C" && resultRow[posIdx] == "874816" && resultRow[altIdx] == "+CCCCTCATCACCTCCCCAGCCACGGTGAGGACCCACCCTGGCATGATCT" {
  			t.Log("OK: Intercolated insertion +CCCCTCATCACCTCCCCAGCCACGGTGAGGACCCACCCTGGCATGATCT in complex, multiallelic microsatellite recapitulated", resultRow)
  		} else {
  			t.Error("NOT OK: Intercolated insertion +CCCCTCATCACCTCCCCAGCCACGGTGAGGACCCACCCTGGCATGATCT in complex, multiallelic microsatellite not recapitulated", resultRow)
  		}
		}

		if index == 1 {
  		if resultRow[refIdx] == "C" && resultRow[posIdx] == "874816" && resultRow[altIdx] == "G" {
  			t.Log("OK: SNP in complex ultiallelic microsatellite recapitulated", resultRow)
  		} else {
  			t.Error("NOT OK: SNP in complex multiallelic microsatellite recapitulated", resultRow)
  		}
		}

		if index == 2 {
  		if resultRow[refIdx] == "C" && resultRow[posIdx] == "874817" && resultRow[altIdx] == "-49" {
  			t.Log("OK: 49bp DEL in complex microsatellite multiallelic recapitulated", resultRow)
  		} else {
  			t.Error("NOT OK: 49bp DEL in complex microsatellite multiallelic recapitulated", resultRow)
  		}
		}

  	if index == 3 {
  		if resultRow[refIdx] == "C" && resultRow[posIdx] == "874816" && resultRow[altIdx] == "+T" {
  			t.Log("OK: Intercolated insertion +T recapitulated", resultRow)
  		} else {
  			t.Error("NOT OK: Intercolated insertion +T not recapitulated", resultRow)
  		}
		}
	})

	if index == 3 {
		t.Log("Ok, recapitulated 4 alleles")
	} else {
		t.Error("expected 4 alleles")
	}
}