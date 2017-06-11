package main

import (
	"strings"
	"testing"
	"sync"
	"strconv"
)

func TestUpdateFieldsWithAlt(t *testing.T) {
	expType, exPos, expRef, expAlt := "SNP", "100", "T", "C"

	sType, pos, ref, alt, err := updateFieldsWithAlt("T", "C", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed")
	} else {
		t.Log("OK: SNP")
	}

	expType, exPos, expRef, expAlt = "SNP", "103", "T", "A"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TCCT", "TCCA", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: SNPs that are longer than 1 base are suported when SNP at end")
	}

	expType, exPos, expRef, expAlt = "SNP", "102", "C", "A"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TGCT", "TGAT", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: SNPs that are longer than 1 base are suported when SNP in middle")
	}

	expType, exPos, expRef, expAlt = "SNP", "100", "T", "A"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TGCT", "AGCT", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: SNPs that are longer than 1 base are suported when SNP at beginning")
	}

	expType, exPos, expRef, expAlt = "", "", "", ""

	sType, pos, ref, alt, err = updateFieldsWithAlt("TCCT", "GTAA", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: MNPs are not supported")
	}

	expType, exPos, expRef, expAlt = "DEL", "101", "C", "-1"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TC", "T", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: 1-based deletions ")
	}

	expType, exPos, expRef, expAlt = "DEL", "101", "A", "-5"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TAGCGT", "T", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Deletions with references longer than 2 bases")
	}

	expType, exPos, expRef, expAlt = "DEL", "102", "G", "-4"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TAGCTT", "TA", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Deletions longer than 1 base")
	}

	expType, exPos, expRef, expAlt = "", "", "", ""

	sType, pos, ref, alt, err = updateFieldsWithAlt("TAGCTT", "TAC", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Left edge of deletion must match exactly")
	}

	expType, exPos, expRef, expAlt = "INS", "100", "T", "+AGCTT"

	sType, pos, ref, alt, err = updateFieldsWithAlt("T", "TAGCTT", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Insertions where reference is 1 base long")
	}

	expType, exPos, expRef, expAlt = "INS", "101", "A", "+GCTT"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TA", "TAGCTT", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Insertions where reference is 2 bases long")
	}

	expType, exPos, expRef, expAlt = "", "", "", ""

	sType, pos, ref, alt, err = updateFieldsWithAlt("TT", "TAGCTT", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Insertions where reference and alt don't share a left edge are skipped")
	}
}

func TestPassesLine(t *testing.T) {
	expect := true

	header := []string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"}
	//Test example 5.2.4 https://samtools.github.io/hts-specs/VCFv4.2.pdf
	record := []string{"20", "4", ".", "GCG", "G,GCGCG", ".", "PASS", "DP=100"}

	actual := linePasses(record, header)

	if actual == expect {
		t.Log("OK: PASS lines pass")
	} else {
		t.Error("NOT OK: PASS lines should pass")
	}

	record = []string{"20", "4", ".", "GCG", "G,GCGCG", ".", ".", "DP=100"}

	actual = linePasses(record, header)

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
	var actualHoms []string
	var actualHets []string

	header := []string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1", "S2", "S3", "S4"}

	sharedFieldsGT := []string{"10", "1000", "rs#", "C", "T", "100", "PASS", "AC=1", "GT"}

	// GT:DS:GL is what 1000 genomes phase 1 provides
	sharedFieldsGTcomplex := []string{"10", "1000", "rs#", "C", "T", "100", "PASS", "AC=1", "GT:DS:GL"}

	fields := append(sharedFieldsGT, "0|0", "0|0", "0|0", "0|0")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	makeHetHomozygotes(fields, header, &actualHoms, &actualHets, "1")

	if len(actualHoms) == 0 && len(actualHets) == 0 {
		t.Log("OK: Homozygous reference samples are skipped")
	} else {
		t.Error("NOT OK: 0 alleles give unexpected results", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, "0|1", "0|1", "0|1", "0|1")

	actualHoms = actualHoms[:0]
	actualHets = actualHets[:0]
	// The allele index we want to test is always 1...unless it's a multiallelic site
	makeHetHomozygotes(fields, header, &actualHoms, &actualHets, "1")

	if len(actualHoms) == 0 && len(actualHets) == 4 {
		t.Log("OK: handles hets")
	} else {
		t.Error("NOT OK: 0 alleles give unexpected results", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, ".|1", "0|1", "0|1", "0|1")

	actualHoms = actualHoms[:0]
	actualHets = actualHets[:0]
	// The allele index we want to test is always 1...unless it's a multiallelic site
	makeHetHomozygotes(fields, header, &actualHoms, &actualHets, "1")

	if len(actualHoms) == 0 && len(actualHets) == 3 {
		t.Log("OK: GT's containing missing data are entirely uncertain, therefore skipped")
	} else {
		t.Error("NOT OK: Fails to handle missing data", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, "1|1", "1|1", "0|1", "0|1")

	actualHoms = actualHoms[:0]
	actualHets = actualHets[:0]
	// The allele index we want to test is always 1...unless it's a multiallelic site
	makeHetHomozygotes(fields, header, &actualHoms, &actualHets, "1")

	if len(actualHoms) == 2 && len(actualHets) == 2 {
		t.Log("OK: handles homs and hets")
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, "1|2", "1|1", "0|1", "0|1")

	actualHoms = actualHoms[:0]
	actualHets = actualHets[:0]
	// The allele index we want to test is always 1...unless it's a multiallelic site
	makeHetHomozygotes(fields, header, &actualHoms, &actualHets, "1")

	if len(actualHoms) == 1 && len(actualHets) == 3 {
		t.Log("OK: a sample heterozygous for a wanted allele is heterozygous for that allele even if its other allele is unwanted (for multiallelic phasing)")
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, "1|2", "1|1", "0|1", "0|1")

	actualHoms = actualHoms[:0]
	actualHets = actualHets[:0]
	// The allele index we want to test is always 1...unless it's a multiallelic site
	makeHetHomozygotes(fields, header, &actualHoms, &actualHets, "2")

	if len(actualHoms) == 0 && len(actualHets) == 1 {
		t.Log("OK: Het / homozygous status is based purely on the wanted allele, rather than total non-ref count")
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGTcomplex, "1|2:-0.03,-1.12,-5.00", "1|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00")

	actualHoms = actualHoms[:0]
	actualHets = actualHets[:0]
	// The allele index we want to test is always 1...unless it's a multiallelic site
	makeHetHomozygotes(fields, header, &actualHoms, &actualHets, "2")

	if len(actualHoms) == 0 && len(actualHets) == 1 {
		t.Log("OK: handles complicated GTs, with non-1 alleles", fields)
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGTcomplex, "1|2|1:-0.03,-1.12,-5.00", "1|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00")

	actualHoms = actualHoms[:0]
	actualHets = actualHets[:0]
	// The allele index we want to test is always 1...unless it's a multiallelic site
	makeHetHomozygotes(fields, header, &actualHoms, &actualHets, "2")

	if len(actualHoms) == 0 && len(actualHets) == 1 {
		t.Log("OK: Complicated GT: Triploids are considered het if only 1 present desired allele", fields)
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, "1|2|1", "1|1", "0|1", "0|1")

	actualHoms = actualHoms[:0]
	actualHets = actualHets[:0]
	// The allele index we want to test is always 1...unless it's a multiallelic site
	makeHetHomozygotes(fields, header, &actualHoms, &actualHets, "2")

	if len(actualHoms) == 0 && len(actualHets) == 1 {
		t.Log("OK: Triploids are considered het if only 1 present desired allele")
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGTcomplex, "2|2|2:-0.03,-1.12,-5.00", "1|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00")

	actualHoms = actualHoms[:0]
	actualHets = actualHets[:0]
	// The allele index we want to test is always 1...unless it's a multiallelic site
	makeHetHomozygotes(fields, header, &actualHoms, &actualHets, "2")

	if len(actualHoms) == 1 && len(actualHets) == 0 {
		t.Log("OK: Complicated GT: Triploids are considered hom only if all alleles present are desired", fields)
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, "2|2|2", "1|1", "0|1", "0|1")

	actualHoms = actualHoms[:0]
	actualHets = actualHets[:0]
	// The allele index we want to test is always 1...unless it's a multiallelic site
	makeHetHomozygotes(fields, header, &actualHoms, &actualHets, "2")

	if len(actualHoms) == 1 && len(actualHets) == 0 {
		t.Log("OK: Triploids are considered hom only if all alleles present are desired")
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}
}

func TestNomralizeSampleNames(t *testing.T) {
	header := []string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1.HAHAHAH", "S2.TRYINGTO.MESSYOUUP", "S3", "S-4"}

	normalizeSampleNames(header)

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
	header := []string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"}
	record := []string{"10", "1000", "rs#", "C", "T", "100", "PASS", "AC=1", "GT"}

	emptyField := "!"
	fieldDelimiter := ";"
	retainId := false;
	retainInfo := true;

	c := make(chan string)
	// I think we need a wait group, not sure.
	wg := new(sync.WaitGroup)

	go func(){
		wg.Add(1)
		processSingleLine(record, header, emptyField, fieldDelimiter, retainId, retainInfo, c, wg)
		wg.Wait();
		close(c)
	}()

  for row := range c {
  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if len(resultRow) == 9 {
  		t.Log("OK: With retainInfo flag set, but not retainID, should output 8 fields")
  	} else {
  		t.Error("NOT OK: With retainInfo flag set, but not retainID, should output 8 fields")
  	}
  	
		if resultRow[7] == "0" && resultRow[8] == "AC=1" {
			t.Log("OK: add INFO field correctly for single field")
		} else {
			t.Error("NOT OK: Couldn't add INFO field", resultRow)
		}
	}

	record = []string{"10", "1000", "rs#", "C", "T,G", "100", "PASS", "AC=1", "GT"}

	c = make(chan string)

	go func(){
		wg.Add(1)
		processMultiLine(record, header, emptyField, fieldDelimiter, retainId, retainInfo, c, wg)
		wg.Wait();
		close(c)
	}()

	index := -1;
  for row := range c {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if len(resultRow) == 9 {
  		t.Log("OK: With retainInfo flag set, but not retainID, should output 8 fields")
  	} else {
  		t.Error("NOT OK: With retainInfo flag set, but not retainID, should output 8 fields")
  	}

  	altIdx, err := strconv.Atoi(resultRow[7])

  	if err != nil {
  		t.Error("NOT OK: The 8th column should be numeric")
  	}

  	if altIdx == index && resultRow[8] == "AC=1" {
			t.Log("OK: add INFO field correctly for multiple field, index", altIdx)
		} else {
			t.Error("NOT OK: Couldn't add INFO field", resultRow)
		}
	}

	if index != 1 {
		t.Error("NOT OK: Expected to parse 2 alleles, parsed fewer")
	}
}

func TestOutputsId(t *testing.T) {
	header := []string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"}
	record := []string{"10", "1000", "rs123", "C", "T", "100", "PASS", "AC=1", "GT"}

	emptyField := "!"
	fieldDelimiter := ";"
	retainId := true;
	retainInfo := false;

	c := make(chan string)
	// I think we need a wait group, not sure.
	wg := new(sync.WaitGroup)

	go func(){
		wg.Add(1)
		processSingleLine(record, header, emptyField, fieldDelimiter, retainId, retainInfo, c, wg)
		wg.Wait();
		close(c)
	}()

  for row := range c {
  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if len(resultRow) == 8 {
  		t.Log("OK: With retainID flag set, but not retainInfo, should output 8 fields")
  	} else {
  		t.Error("NOT OK: With retainID flag set, but not retainInfo, should output 8 fields")
  	}

		if resultRow[7] == "rs123" {
			t.Log("OK: add ID field correctly for single field")
		} else {
			t.Error("NOT OK: Couldn't add ID field", resultRow)
		}
	}

	record = []string{"10", "1000", "rs456", "C", "T,G", "100", "PASS", "AC=1", "GT"}

	c = make(chan string)

	go func(){
		wg.Add(1)
		processMultiLine(record, header, emptyField, fieldDelimiter, retainId, retainInfo, c, wg)
		wg.Wait();
		close(c)
	}()

	index := -1
  for row := range c {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if len(resultRow) == 8 {
  		t.Log("OK: With retainID flag set, but not retainInfo, should output 8 fields")
  	} else {
  		t.Error("NOT OK: With retainID flag set, but not retainInfo, should output 8 fields")
  	}

  	if resultRow[7] == "rs456" {
			t.Log("OK: add ID field correctly for multiple field, index", altIdx)
		} else {
			t.Error("NOT OK: Couldn't add ID field", resultRow)
		}
	}

	if index != 1 {
		t.Error("NOT OK: Expected to parse 2 alleles, parsed fewer")
	}
}

func TestOutputsIdAndInfo(t *testing.T) {
	header := []string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"}
	record := []string{"10", "1000", "rs123", "C", "T", "100", "PASS", "AC=1", "GT"}

	emptyField := "!"
	fieldDelimiter := ";"
	retainId := true;
	retainInfo := true;

	c := make(chan string)
	// I think we need a wait group, not sure.
	wg := new(sync.WaitGroup)

	go func(){
		wg.Add(1)
		processSingleLine(record, header, emptyField, fieldDelimiter, retainId, retainInfo, c, wg)
		wg.Wait();
		close(c)
	}()

  for row := range c {
  	resultRow := strings.Split(row[:len(row)-1], "\t")
  	
  	if len(resultRow) == 10 {
  		t.Log("OK: With both retainID and ratinInfo flags set, should output 10 fields")
  	} else {
  		t.Error("NOT OK: With both retainID and ratinInfo flags set, should output 10 fields")
  	}

		if resultRow[7] == "rs123" {
			t.Log("OK: add ID field correctly for single field")
		} else {
			t.Error("NOT OK: Couldn't add ID field", resultRow)
		}

		if resultRow[8] == "0" && resultRow[9] == "AC=1" {
			t.Log("OK: add INFO field correctly for single field")
		} else {
			t.Error("NOT OK: Couldn't add INFO field", resultRow)
		}
	}

	record = []string{"10", "1000", "rs456", "C", "T,G", "100", "PASS", "AC=1", "GT"}

	c = make(chan string)

	go func(){
		wg.Add(1)
		processMultiLine(record, header, emptyField, fieldDelimiter, retainId, retainInfo, c, wg)
		wg.Wait();
		close(c)
	}()

	index := -1
  for row := range c {
  	index++
 
  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if len(resultRow) == 10 {
  		t.Log("OK: With both retainID and ratinInfo flags set, should output 10 fields")
  	} else {
  		t.Error("With both retainID and ratinInfo flags set, should output 10 fields")
  	}

  	if resultRow[7] == "rs456" {
			t.Log("OK: add ID field correctly for multiple field, index", altIdx)
		} else {
			t.Error("NOT OK: Couldn't add ID field", altIdx, record[8])
		}

		altIdx, err := strconv.Atoi(resultRow[8])

  	if err != nil {
  		t.Error("NOT OK: The 9th column should be numeric")
  	}

  	if altIdx == index {
  		t.Log("OK: Multiallelic index is in 9th column when retainInfo is true")
		} else {
			t.Error("NOT OK: Multiallelic index isn't in 9th column when retainInfo is true", resultRow)
		}

		if resultRow[9] == "AC=1" {
			t.Log("OK: add INFO field correctly for multiallelic field in column 10")
		} else {
			t.Error("NOT OK: Couldn't add INFO field", resultRow)
		}
	}

	if index != 1 {
		t.Error("NOT OK: Expected to parse 2 alleles, parsed fewer")
	}
}

func TestOutputMultiallelic(t *testing.T) {
	header := []string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"}
	//Test example 5.2.4 https://samtools.github.io/hts-specs/VCFv4.2.pdf
	//Except use A as last base, to properly distinguish first and last base of reference
	record := []string{"20", "4", ".", "GCA", "G,GCACG", ".", "PASS", "DP=100"}

	emptyField := "!"
	fieldDelimiter := ";"

	c := make(chan string)
	// I think we need a wait group, not sure.
	wg := new(sync.WaitGroup)

	go func(){
		if linePasses(record, header) == false {
			t.Error("NOT OK: Line should pass", record, header)
		}

		wg.Add(1)
		go processLine(record, header, emptyField, fieldDelimiter, false, false, c, wg)
		wg.Wait();
		close(c)
	}()

	index := -1;
  for row := range c {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if index == 0 {
  		if resultRow[altIdx] == "-2" {
  			t.Log("OK: 2 base deletion recapitulated", resultRow)
  		} else {
  			t.Error("NOT OK: 2 base deletion not recapitulated", resultRow)
  		}
			
			if resultRow[posIdx] == "5" {
				t.Log("OK: 2 base deletion should be shifted by the length of the deleted allele", resultRow)
			} else {
				// This means by len(G) in this case
				t.Error("NOT OK: 2 base deletion should be shifted by the length of the deleted allele", resultRow)
			}

			if resultRow[refIdx] == "C" {
				t.Log("OK: In deletion, reference base is the first deleted base", resultRow)
			} else {
				t.Error("NOT OK: In deletion, reference base is the first deleted base", resultRow)
			}
		}

		if index == 1 {
  		if resultRow[altIdx] == "+CG" {
  			t.Log("OK: 2 base deletion recapitulated", resultRow)
  		} else {
  			t.Error("NOT OK: 2 base deletion not recapitulated", resultRow)
  		}
			
			if resultRow[posIdx] == "6" {
				t.Log("OK: Insertion should be shifted by padded base to last reference base", resultRow)
			} else {
				// This means by len(G) in this case
				t.Error("NOT OK: Insertion should be shifted by padded base to last reference base", resultRow)
			}

			if resultRow[refIdx] == "A" {
				t.Log("OK: In insertion, reference base should be last reference in REF string", resultRow)
			} else {
				t.Error("NOT OK: In insertion, reference base should be last reference in REF string", resultRow)
			}
		}
	}

	if index == 1 {
		t.Log("OK: parsed 2 alleles")
	} else {
		t.Error("NOT OK: Expected to parse 2 alleles")
	}
}