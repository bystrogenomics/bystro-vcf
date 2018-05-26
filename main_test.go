package main

import (
	"strings"
	"testing"
	"strconv"
	"fmt"
	"bufio"
)

func TestKeepFlagsTrue(t *testing.T) {
	args := []string{
		"--keepInfo",
		"--keepId",
		"--keepPos",
		"--inPath", "/path/to/file",
		"--errPath", "/path/to/err",
		"--cpuProfile", "/path/to/profile",
		"--emptyField", ".",
		"--fieldDelimiter", "&",
	}

	config := setup(args)

  if !(config.keepInfo == true && config.keepId == true && config.keepPos == true)  {
  	t.Error("NOT OK: parse keepInfo, keepId, keepPos args")
  }

  if config.inPath != "/path/to/file" || config.errPath != "/path/to/err" ||
  config.cpuProfile != "/path/to/profile" {
  	t.Error("NOT OK: parse inPath and errPath args")
  }

  if config.emptyField != "." || config.fieldDelimiter != "&" {
  	t.Error("NOT OK: parse emptyField and fieldDelimiter args")
  }
}

func TestHeader(t *testing.T) {
	config := Config{keepId: false, keepInfo: false}

	header := stringHeader(&config)

	expected := strings.Join([]string{"chrom", "pos", "type", "ref", "alt", "trTv", "heterozygotes",
    "heterozygosity", "homozygotes", "homozygosity", "missingGenos", "missingness", "sampleMaf"}, "\t")

	if header == expected + "\n" {
		t.Log("OK: print header", header, expected)
	} else {
		t.Error("NOT OK: print header", header, expected)
	}

	config = Config{keepPos: true, keepId: false, keepInfo: false}

	header = stringHeader(&config)

	expected = strings.Join([]string{"chrom", "pos", "type", "ref", "alt", "trTv", "heterozygotes",
    "heterozygosity", "homozygotes", "homozygosity", "missingGenos", "missingness", "sampleMaf", "vcfPos"}, "\t")

	if header == expected + "\n" {
		t.Log("OK: print header with --keepPos true", header, expected)
	} else {
		t.Error("NOT OK: print header with --keepPostrue true", header, expected)
	}

	config = Config{keepPos: true, keepId: true, keepInfo: false}

	header = stringHeader(&config)

	expected = strings.Join([]string{"chrom", "pos", "type", "ref", "alt", "trTv", "heterozygotes",
    "heterozygosity", "homozygotes", "homozygosity", "missingGenos", "missingness", "sampleMaf", "vcfPos", "id"}, "\t")

	if header == expected + "\n" {
		t.Log("OK: print header with --keepId true", header, expected)
	} else {
		t.Error("NOT OK: print header with --keepId true", header, expected)
	}

	config = Config{keepPos: false, keepId: false, keepInfo: true}

	header = stringHeader(&config)

	expected = strings.Join([]string{"chrom", "pos", "type", "ref", "alt", "trTv", "heterozygotes",
    "heterozygosity", "homozygotes", "homozygosity", "missingGenos", "missingness",
    "sampleMaf", "alleleIdx", "info"}, "\t")

	if header == expected + "\n" {
		t.Log("OK: print header with --keepInfo true", header, expected)
	} else {
		t.Error("NOT OK: print header with --keepInfo true", header, expected)
	}

	config = Config{keepPos: false, keepId: true, keepInfo: true}

	header = stringHeader(&config)

	expected = strings.Join([]string{"chrom", "pos", "type", "ref", "alt", "trTv", "heterozygotes",
    "heterozygosity", "homozygotes", "homozygosity", "missingGenos", "missingness",
    "sampleMaf", "id", "alleleIdx", "info"}, "\t")

	if header == expected + "\n" {
		t.Log("OK: print header with -keepId true --keepInfo true", header)
	} else {
		t.Error("NOT OK: print header with --keepId true --keepInfo true", header)
	}

	config = Config{keepPos: true, keepId: true, keepInfo: true}

	header = stringHeader(&config)

	expected = strings.Join([]string{"chrom", "pos", "type", "ref", "alt", "trTv", "heterozygotes",
    "heterozygosity", "homozygotes", "homozygosity", "missingGenos", "missingness",
    "sampleMaf", "vcfPos", "id", "alleleIdx", "info"}, "\t")

	if header == expected + "\n" {
		t.Log("OK: print header with --keepPos true --keepId true --keepInfo true", header)
	} else {
		t.Error("NOT OK: print header with --keepId true --keepInfo true", header)
	}

	config = Config{keepPos: true, keepId: false, keepInfo: true}

	header = stringHeader(&config)

	expected = strings.Join([]string{"chrom", "pos", "type", "ref", "alt", "trTv", "heterozygotes",
    "heterozygosity", "homozygotes", "homozygosity", "missingGenos", "missingness",
    "sampleMaf", "vcfPos", "alleleIdx", "info"}, "\t")

	if header == expected + "\n" {
		t.Log("OK: print header with --keepPos true --keepId false --keepInfo true", header)
	} else {
		t.Error("NOT OK: print header with --keepId true --keepInfo true", header)
	}
}

func TestUpdateFieldsWithAlt(t *testing.T) {
	expType, exPos, expRef, expAlt := "SNP", "100", byte('T'), "C"

	sType, pos, refs, alts, altIndices := getAlleles("chr1", "100", "T", "C")
	test := "Has alt indices"

	if len(altIndices) != 1 && altIndices[1] != 0 {
		t.Errorf("NOT OK: %s", test)
	} else {
		t.Logf("OK: %s", test)
	}

	if sType != expType || pos[0] != exPos || refs[0] != expRef || alts[0] != expAlt {
		t.Error("NOT OK: Test failed")
	} else {
		t.Log("OK: SNP")
	}

	expType, exPos, expRef, expAlt = "SNP", "103", byte('T'), "A"

	sType, pos, refs, alts, altIndices = getAlleles("chr1", "100", "TCCT", "TCCA")

	if sType != expType || pos[0] != exPos || refs[0] != expRef || alts[0] != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, refs, alts)
	} else {
		t.Log("OK: SNPs that are longer than 1 base are suported when SNP at end")
	}

	expType, exPos, expRef, expAlt = "SNP", "102", byte('C'), "A"

	sType, pos, refs, alts, altIndices = getAlleles("chr1", "100", "TGCT", "TGAT")

	if sType != expType || pos[0] != exPos || refs[0] != expRef || alts[0] != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, refs, alts)
	} else {
		t.Log("OK: SNPs that are longer than 1 base are suported when SNP in middle")
	}

	expType, exPos, expRef, expAlt = "SNP", "100", byte('T'), "A"

	sType, pos, refs, alts, altIndices = getAlleles("chr1", "100", "TGCT", "AGCT")

	if sType != expType || pos[0] != exPos || refs[0] != expRef || alts[0] != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, refs, alts)
	} else {
		t.Log("OK: SNPs that are longer than 1 base are suported when SNP at beginning")
	}

	/******************** Test contiguous MNPs ***************/
	expType = "MNP"
	expRefs := []byte{'T', 'C', 'G', 'T'}
	expAlts := []string{"G", "T", "A", "A"}
	expPositions := []string{"100", "101", "102", "103"}
	// indices are relative to the VCF order; 1 MNP is 1 allele
	expIndices := []int{0, 0, 0, 0}

	sType, pos, refs, alts, altIndices = getAlleles("chr1", "100", "TCGT", "GTAA")
	test = "Supports MNPs"
	// TODO: Should we call MNPs SNP or MNP?
	if sType != expType {
		t.Errorf("NOT OK: %s : MNPs should be labeled 'MNP' if they have > 1 allele", test, sType, expType)
	} else {
		t.Logf("OK: %s : MNPs should be labeled 'MNP' if they have > 1 allele", test)
	}

	if len(altIndices) != 4 {
		t.Errorf("NOT OK: %s : All MNP alleles should be reported", test)
	} else {
		t.Logf("OK: %s : All MNP alleles should be reported", test)
	}

	for i := 0; i < len(altIndices); i++ {
		if pos[i] != expPositions[i] || refs[i] != expRefs[i] || alts[i] != expAlts[i] || altIndices[i] != expIndices[i] {
			t.Errorf("NOT OK: %s : mangled MNP annotation", test, sType, pos, refs, alts)
		} else {
			t.Logf("OK: %s : correct MNP annotation", test)
		}
	}

	/******************** Test sparse MNPs ***************/
	// not certain if these are possible, but will be annotated correctly
	expType = "MNP"
	expRefs = []byte{'C', 'T'}
	expAlts = []string{"A", "C"}
	expPositions = []string{"101", "103"}
	// indices are relative to the VCF order; 1 MNP is 1 allele
	expIndices = []int{0, 0}

	sType, pos, refs, alts, altIndices = getAlleles("chr1", "100", "TCGT", "TAGC")
	test = "Supports sparse MNPs"

	// TODO: Should we call MNPs SNP or MNP?
	if sType != expType {
		t.Errorf("NOT OK: %s : MNPs should be labeled 'MNP' if they have > 1 allele", test, sType, expType)
	} else {
		t.Logf("OK: %s : MNPs should be labeled 'MNP' if they have > 1 allele", test)
	}

	if len(altIndices) != 2 {
		t.Errorf("NOT OK: %s : All MNP alleles should be reported", test)
	} else {
		t.Logf("OK: %s : All MNP alleles should be reported", test)
	}

	for i := 0; i < len(altIndices); i++ {
		if pos[i] != expPositions[i] || refs[i] != expRefs[i] || alts[i] != expAlts[i] || altIndices[i] != expIndices[i] {
			t.Errorf("NOT OK: %s : mangled MNP annotation", test, sType, pos, refs, alts)
		} else {
			t.Logf("OK: %s : correct MNP annotation", test)
		}
	}

	/******************** Test MNP-like SNP***************/
	// not certain if these are possible, but will be annotated correctly
	expType = "SNP"
	expRefs = []byte{'T'}
	expAlts = []string{"C"}
	expPositions = []string{"103"}
	// indices are relative to the VCF order; 1 MNP is 1 allele
	expIndices = []int{0}

	sType, pos, refs, alts, altIndices = getAlleles("chr1", "100", "TCGT", "TCGC")
	test = "Supports SNPs that are reported as an N-long block"

	if sType != expType {
		t.Errorf("NOT OK: %s : sites should be labeled 'SNP' if they have only 1 change, that is a single base", test, sType, expType)
	} else {
		t.Logf("OK: %s : sites should be labeled 'SNP' if they have only 1 change, that is a single base", test)
	}

	if len(altIndices) != 1 {
		t.Errorf("NOT OK: %s : SNPs that are reported as N-block should have only 1 allele", test)
	} else {
		t.Logf("OK: %s : SNPs that are reported as N-block should have only 1 allele", test)
	}

	for i := 0; i < len(altIndices); i++ {
		if pos[i] != expPositions[i] || refs[i] != expRefs[i] || alts[i] != expAlts[i] || altIndices[i] != expIndices[i] {
			t.Errorf("NOT OK: %s : mangled long-SNP annotation", test, sType, pos, refs, alts)
		} else {
			t.Logf("OK: %s : correct long-SNP annotation", test)
		}
	}

	/******************** Test simple 1 base deletion ***************/
	expType, exPos, expRef, expAlt = "DEL", "101", byte('C'), "-1"

	sType, pos, refs, alts, altIndices = getAlleles("chr1", "100", "TC", "T")

	if sType != expType || pos[0] != exPos || refs[0] != expRef || alts[0] != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, refs, alts)
	} else {
		t.Log("OK: 1-based deletions ")
	}

	expType, exPos, expRef, expAlt = "DEL", "101", byte('A'), "-5"

	sType, pos, refs, alts, altIndices = getAlleles("chr1", "100", "TAGCGT", "T")

	if sType != expType || pos[0] != exPos || refs[0] != expRef || alts[0] != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, refs, alts)
	} else {
		t.Log("OK: Deletions with references longer than 2 bases")
	}

	// Test multiallelic intercolated deletion
	expType, exPos, expRef, expAlt = "DEL", "102", byte('G'), "-4"

	sType, pos, refs, alts, altIndices = getAlleles("chr1", "100", "TAGCTT", "TA")

	if sType != expType || pos[0] != exPos || refs[0] != expRef || alts[0] != expAlt {
		t.Error("NOT OK: Test failed", sType, pos[0], string(refs[0]), alts[0])
	} else {
		t.Log("OK: Deletions longer than 1 base")
	}

	//Check a malformed allele
	expType, exPos, expRef, expAlt = "", "", 0, ""

	sType, pos, refs, alts, altIndices = getAlleles("chr1", "100", "TAGCTT", "TAC")

	if sType != expType || len(pos) != 0 || len(refs) != 0 || len(alts) != 0 {
		t.Error("NOT OK: expect ref:TAGCTT alt:TAC to return 0 alleles", sType, pos, refs, alts)
	} else {
		t.Log("OK: expected ref:TAGCTT alt:TAT to return 0 allele")
	}

	//Check an intercolated deletion
	expType, exPos, expRef, expAlt = "DEL", "102", byte('G'), "-3"

	sType, pos, refs, alts, altIndices = getAlleles("chr1", "100", "TAGCTT", "TAT")

	if len(pos) != 1 || len(refs) != 1 || len(alts) != 1 || len(altIndices) != 1 {
		t.Error("NOT OK: expected ref:TAGCTT alt:TAT to return 1 intercolated allele", sType, pos, refs, alts)
	} else {
		t.Log("OK: expected ref:TAGCTT alt:TAT to return 1 intercolated allele")
	}

	if sType != expType || pos[0] != exPos || refs[0] != expRef || alts[0] != expAlt {
		t.Error("NOT OK: expect ref:TAGCTT alt:TAC to be pos:102, ref:G, alt:-3", sType, pos, refs, alts)
	} else {
		t.Log("OK: expect ref:TAGCTT alt:TAC to be pos:102, ref:G, alt:-3")
	}

	expType, exPos, expRef, expAlt = "INS", "100", byte('T'), "+AGCTT"

	sType, pos, refs, alts, altIndices = getAlleles("chr1", "100", "T", "TAGCTT")

	test = "Insertions where reference is 1 base long"

	if sType != expType || pos[0] != exPos || refs[0] != expRef || alts[0] != expAlt {
		t.Errorf("NOT OK: %s", test, sType, pos, refs, alts)
	} else {
		t.Logf("OK: %s", test)
	}

	// Test intercolated insertions
	expType, exPos, expRef, expAlt = "INS", "100", byte('T'), "+AGCT"
	test = "Intercolated insertions"

	sType, pos, refs, alts, altIndices = getAlleles("chr1", "100", "TT", "TAGCTT")

	if sType != expType || pos[0] != exPos || refs[0] != expRef || alts[0] != expAlt {
		t.Error("NOT OK: Test failed", sType, pos, string(refs[0]), alts)
	} else {
		t.Log("OK: Multiallelics insertion where reference is 2 bases long")
	}
}

func TestPassesLine(t *testing.T) {
	expect := true

	header := []string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"}
	//Test example 5.2.4 https://samtools.github.io/hts-specs/VCFv4.2.pdf
	record := []string{"20", "4", ".", "GCG", "G,GCGCG", ".", "PASS", "DP=100"}

	allowedFilters := map[string]bool{ "PASS": true, ".": true}
	actual := linePasses(record, header, allowedFilters)

	if actual == expect {
		t.Log("OK: PASS lines pass")
	} else {
		t.Error("NOT OK: PASS lines should pass")
	}

	record = []string{"20", "4", ".", "GCG", "G,GCGCG", ".", ".", "DP=100"}

	actual = linePasses(record, header, allowedFilters)

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
	actualHoms, actualHets, missing, sampleMaf := makeHetHomozygotes(fields, header, '1')

	if len(actualHoms) == 0 && len(actualHets) == 0 && len(missing) == 0 {
		t.Log("OK: Homozygous reference samples are skipped")
	} else {
		t.Error("NOT OK: 0 alleles give unexpected results", actualHoms, actualHets, missing)
	}

	if sampleMaf == 0 {
		t.Log("OK: sampleMaf should be 0 for all reference sites", strconv.FormatFloat(sampleMaf, 'G', 3, 64))
	} else {
		t.Error("NOT OK: sampleMaf should be 0 for all reference sites", strconv.FormatFloat(sampleMaf, 'G', 3, 64))
	}

	fields = append(sharedFieldsGT, "0|1", "0|1", "0|1", "0|1")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing, sampleMaf = makeHetHomozygotes(fields, header, '1')

	if len(actualHoms) == 0 && len(actualHets) == 4 && len(missing) == 0 {
		t.Log("OK: handles hets")
	} else {
		t.Error("NOT OK: 0 alleles give unexpected results", actualHoms, actualHets, missing)
	}

	if sampleMaf == float64(4)/float64(8) {
		t.Log("OK: sampleMaf should be .5 for 0|1, 0|1, 0|1, 0|1", strconv.FormatFloat(sampleMaf, 'G', 3, 64))
	} else {
		t.Error("NOT OK: sampleMaf miscounted for hets")
	}

	fields = append(sharedFieldsGT, ".|1", "0|1", "0|1", "0|1", strconv.FormatFloat(sampleMaf, 'G', 3, 64))

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing, sampleMaf = makeHetHomozygotes(fields, header, '1')

	if len(actualHoms) == 0 && len(actualHets) == 3 && len(missing) == 1 {
		t.Log("OK: GT's containing missing data are entirely uncertain, therefore skipped. Missing genotypes are called if any of the calls are missing")
	} else {
		t.Error("NOT OK: Fails to handle missing data", actualHoms, actualHets, missing)
	}

	if sampleMaf == float64(3)/float64(6) {
		t.Log("OK: sampleMaf should be .5 for .|1, 0|1, 0|1, 0|1", strconv.FormatFloat(sampleMaf, 'G', 3, 64))
	} else {
		t.Error("NOT OK: sampleMaf miscounted for hets in presence of missing data")
	}

	fields = append(sharedFieldsGT, "1|.", "0|1", "0|1", "0|1", strconv.FormatFloat(sampleMaf, 'G', 3, 64))

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing, sampleMaf = makeHetHomozygotes(fields, header, '1')

	if len(actualHoms) == 0 && len(actualHets) == 3 && len(missing) == 1 {
		t.Log("OK: GT's containing missing data are entirely uncertain, therefore skipped. Missing genotypes are called if any of the calls are missing")
	} else {
		t.Error("NOT OK: Fails to handle missing data", actualHoms, actualHets)
	}

	if sampleMaf == float64(3)/float64(6) {
		t.Log("OK: sampleMaf should be .5 for 1|., 0|1, 0|1, 0|1", strconv.FormatFloat(sampleMaf, 'G', 3, 64))
	} else {
		t.Error("NOT OK: sampleMaf miscounted for hets in presence of missing data")
	}

	fields = append(sharedFieldsGT, "1|1", "1|1", "0|1", "0|1", strconv.FormatFloat(sampleMaf, 'G', 3, 64))

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing, sampleMaf = makeHetHomozygotes(fields, header, '1')

	if len(actualHoms) == 2 && len(actualHets) == 2 {
		t.Log("OK: handles homs and hets")
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	if sampleMaf == float64(6)/float64(8) {
		t.Log("OK: sampleMaf should be 6/8 for 1|1, 1|1, 0|1, 0|1", strconv.FormatFloat(sampleMaf, 'G', 3, 64))
	} else {
		t.Error("NOT OK: sampleMaf miscounted for hets + homs admixture", strconv.FormatFloat(sampleMaf, 'G', 3, 64))
	}

	fields = append(sharedFieldsGT, "1|2", "1|1", "0|1", "0|1")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing, sampleMaf = makeHetHomozygotes(fields, header, '1')

	if len(actualHoms) == 1 && len(actualHets) == 3 {
		t.Log("OK: a sample heterozygous for a wanted allele is heterozygous for that allele even if its other allele is unwanted (for multiallelic phasing)")
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	if sampleMaf == float64(5)/float64(8) {
		t.Log("OK: sampleMaf should be 5/8 for 1|2, 1|1, 0|1, 0|1 and allele 1, because we consider only the allele being counted", strconv.FormatFloat(sampleMaf, 'G', 3, 64))
	} else {
		t.Error("NOT OK: sampleMaf miscounted for hets + homs admixture in multiallelic case", strconv.FormatFloat(sampleMaf, 'G', 3, 64))
	}

	fields = append(sharedFieldsGT, "1|2", "1|1", "0|1", "0|1")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing, sampleMaf = makeHetHomozygotes(fields, header, '2')

	if len(actualHoms) == 0 && len(actualHets) == 1 {
		t.Log("OK: Het / homozygous status is based purely on the wanted allele, rather than total non-ref count")
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	if sampleMaf == float64(1)/float64(8) {
		t.Log("OK: sampleMaf should be 1/8 for 1|2, 1|1, 0|1, 0|1 and allele 2", strconv.FormatFloat(sampleMaf, 'G', 3, 64))
	} else {
		t.Error("NOT OK: sampleMaf miscounted for hets + homs admixture in multiallelic case", strconv.FormatFloat(sampleMaf, 'G', 3, 64))
	}

	fields = append(sharedFieldsGTcomplex, "1|2:-0.03,-1.12,-5.00", "1|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing, sampleMaf = makeHetHomozygotes(fields, header, '2')

	if len(actualHoms) == 0 && len(actualHets) == 1 {
		t.Log("OK: handles complicated GTs, with non-1 alleles", fields)
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGTcomplex, "1|2|1:-0.03,-1.12,-5.00", "1|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing, sampleMaf = makeHetHomozygotes(fields, header, '2')

	if len(actualHoms) == 0 && len(actualHets) == 1 {
		t.Log("OK: Complicated GT: Triploids are considered het if only 1 present desired allele", fields)
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, "1|2|1", "1|1", "0|1", "0|1")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing, sampleMaf = makeHetHomozygotes(fields, header, '2')

	if len(actualHoms) == 0 && len(actualHets) == 1 {
		t.Log("OK: Triploids are considered het if only 1 present desired allele")
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGTcomplex, "2|2|2:-0.03,-1.12,-5.00", "1|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00", "0|1:-0.03,-1.12,-5.00")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing, sampleMaf = makeHetHomozygotes(fields, header, '2')

	if len(actualHoms) == 1 && len(actualHets) == 0 {
		t.Log("OK: Complicated GT: Triploids are considered hom only if all alleles present are desired", fields)
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}

	fields = append(sharedFieldsGT, "2|2|2", "1|1", "0|1", "0|1")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing, sampleMaf = makeHetHomozygotes(fields, header, '2')

	if len(actualHoms) == 1 && len(actualHets) == 0 {
		t.Log("OK: Triploids are considered hom only if all alleles present are desired")
	} else {
		t.Error("NOT OK: Fails to handle homs and hets", actualHoms, actualHets)
	}
}

func TestMakeHetHomozygotesHaploid(t *testing.T) {
	header := []string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1", "S2", "S3", "S4"}

	sharedFieldsGT := []string{"10", "1000", "rs#", "C", "T", "100", "PASS", "AC=1", "GT"}

	// GT:DS:GL is what 1000 genomes phase 1 provides
	sharedFieldsGTcomplex := []string{"10", "1000", "rs#", "C", "T", "100", "PASS", "AC=1", "GT:DS:GL"}

	fields := append(sharedFieldsGT, "0", ".", "1", "0")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing, sampleMaf := makeHetHomozygotes(fields, header, '1')

	if len(actualHoms) == 1 && len(actualHets) == 0 && len(missing) == 1 {
		t.Log("OK: Haploid non-reference sites called homozygous")
	} else {
		t.Error("NOT OK: Haploid non-reference sites called homozygous", actualHoms, actualHets, missing)
	}

	if sampleMaf == float64(1)/float64(3) {
		t.Log("OK: Calculate sample maf as 1/3 when 1 non-ref allele, 1 misisng, in haploid", sampleMaf)
	} else {
		t.Error("NOT OK: Calculate sample maf as 1/3 when 1 non-ref allele, 1 misisng, in haploid", sampleMaf)
	}

	fields = append(sharedFieldsGTcomplex, "0:1", ".:1", "1:1", "0:1")

	// The allele index we want to test is always 1...unless it's a multiallelic site
	actualHoms, actualHets, missing, sampleMaf = makeHetHomozygotes(fields, header, '1')

	if len(actualHoms) == 1 && len(actualHets) == 0 && len(missing) == 1 {
		t.Log("OK: Haploid non-reference sites called heterozygous in complex case", sampleMaf)
	} else {
		t.Error("NOT OK: Haploid non-reference sites called heterozygous in complex case", sampleMaf)
	}

	if sampleMaf == float64(1)/float64(3) {
		t.Log("OK: Calculate sample maf as 1/3 when 1 non-ref allele, 1 misisng, in haploid complex GT")
	} else {
		t.Error("NOT OK: Calculate sample maf as 1/3 when 1 non-ref allele, 1 misisng, in haploid complex GT")
	}
}

// Headers are expected to be:
// chrom: 0, pos: 1, siteType: 2, ref: 3, alt: 4, trTv: 5, heterozygotes: 6, heterozygosity: 7
// homozyogtes: 8, homozygosity: 9, missingGenos: 10, missingness: 11, sampleMaf: 12,
// id: 13 (if keepId), alleleIdx: 14 (if keepId and keepInfo), info: 15 (if keepId and keepInfo)
// if keepInfo only: alleleIdx: 13, info: 14
func TestOutputsInfo(t *testing.T) {
	versionLine := "##fileformat=VCFv4.x"
	header := strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"}, "\t")
	record := strings.Join([]string{"10", "1000", "rs#", "C", "T", "100", "PASS", "AC=1"}, "\t")

	lines := versionLine + "\n" + header + "\n" + record	+ "\n"
	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", keepId: false, keepInfo: true}

	readVcf(&config, reader, func(row string) {
  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr10" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if resultRow[4] != "T" {
  		t.Error("Couldn't parse T allele", resultRow)
  	}

  	if resultRow[5] != "1" {
  		t.Error("Transitions should have value 1", resultRow)
  	}

  	if len(resultRow) == 15 {
  		t.Log("OK: With keepInfo flag set, but not keepId, should output 15 fields", resultRow)
  	} else {
  		t.Error("NOT OK: With keepInfo flag set, but not keepId, should output 15 fields", resultRow)
  	}

		if resultRow[len(resultRow) - 2] == "0" && resultRow[len(resultRow) - 1] == "AC=1" {
			t.Log("OK: add INFO field correctly for single field")
		} else {
			t.Error("NOT OK: Couldn't add INFO field", resultRow)
		}
	})

	record = strings.Join([]string{"10", "1000", "rs#", "C", "T,G", "100", "PASS", "AC=1"}, "\t")

	lines = versionLine + "\n" + header + "\n" + record	+ "\n"
	reader = bufio.NewReader(strings.NewReader(lines))

	index := -1;
  readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr10" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if resultRow[5] != "0" {
  		t.Error("Multiallelics should be called neither transitions nor transversions with value 0", resultRow)
  	}

  	if len(resultRow) == 15 {
  		t.Log("OK: With keepInfo flag set, but not keepId, should output 15 fields", resultRow)
  	} else {
  		t.Error("NOT OK: With keepInfo flag set, but not keepId, should output 15 fields", resultRow)
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
	header := strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"}, "\t")
	record := strings.Join([]string{"10", "1000", "rs123", "C", "T", "100", "PASS", "AC=1"}, "\t")

	lines := versionLine + "\n" + header + "\n" + record	+ "\n"

	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", keepId: true, keepInfo: false}

  readVcf(&config, reader, func(row string) {
  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr10" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if resultRow[5] != "1" {
  		t.Error("Transitions should have value 1", resultRow)
  	} else {
  		t.Log("Parsed transition with value 1", resultRow)
  	}

  	if len(resultRow) == 14 {
  		t.Log("OK: With keepId flag set, but not keepInfo, should output 14 fields")
  	} else {
  		t.Error("NOT OK: With keepId flag set, but not keepInfo, should output 14 fields", resultRow)
  	}

		if resultRow[len(resultRow) - 1] == "rs123" {
			t.Log("OK: add ID field correctly for single field")
		} else {
			t.Error("NOT OK: Couldn't add ID field", resultRow)
		}
	})

	record = strings.Join([]string{"10", "1000", "rs456", "C", "T,G", "100", "PASS", "AC=1"}, "\t")

	lines = versionLine + "\n" + header + "\n" + record	+ "\n"
	reader = bufio.NewReader(strings.NewReader(lines))

	index := -1
  readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr10" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if resultRow[5] != "0" {
  		t.Error("Multiallelics should be called neither transitions nor transversions with value 0", resultRow)
  	}

  	if len(resultRow) == 14 {
  		t.Log("OK: With keepId flag set, but not keepInfo, should output 14 fields")
  	} else {
  		t.Error("NOT OK: With keepId flag set, but not keepInfo, should output 14 fields", resultRow)
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

func TestOutputsVcfPos(t *testing.T) {
	versionLine := "##fileformat=VCFv4.x"
	header := strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"}, "\t")

	// Define a interstital insertion, between C & T
	record := strings.Join([]string{"10", "1000", "rs#", "CTT", "CT", "100", "PASS", "AC=1"}, "\t")

	lines := versionLine + "\n" + header + "\n" + record	+ "\n"
	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", keepPos: true}

	readVcf(&config, reader, func(row string) {
  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if len(resultRow) != 14 {
  		t.Error("With keepPos only set, expect 14 fields", resultRow)
  	}

  	vcfPos, err := strconv.Atoi(resultRow[len(resultRow) - 1])

		if err != nil || vcfPos != 1000 {
			t.Error("Deletions should get the vcf input position with --keepPos, in the vcfPos column", resultRow)
		}
	})

	record = strings.Join([]string{"10", "1003", "rs#", "C", "T,G", "100", "PASS", "AC=1"}, "\t")

	lines = versionLine + "\n" + header + "\n" + record	+ "\n"
	reader = bufio.NewReader(strings.NewReader(lines))

	index := -1;
  readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

		if len(resultRow) != 14 {
  		t.Error("With keepPos flag set only, each allele in multiallelics should have 14 fields in its row", resultRow)
  	}

  	vcfPos, err := strconv.Atoi(resultRow[len(resultRow) - 1])

  	if err != nil || vcfPos != 1003 {
			t.Error("Each allele in a multiallelic should get the vcf input position with --keepPos, in the vcfPos column", resultRow)
		}
	})

	if index != 1 {
		t.Error("NOT OK: Expected to parse 2 alleles, parsed fewer")
	}
}

func TestOutputsVcfPosIdAndInfo(t *testing.T) {
	versionLine := "##fileformat=VCFv4.x"
	header := strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"}, "\t")

	// Define a interstital insertion, between C & T
	record := strings.Join([]string{"10", "1000", "rs1", "CTT", "CT", "100", "PASS", "AC=1"}, "\t")

	lines := versionLine + "\n" + header + "\n" + record	+ "\n"
	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", keepPos: true, keepId: true, keepInfo: true}

	readVcf(&config, reader, func(row string) {
  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if len(resultRow) != 17 {
  		t.Error("With keepId, keepInfo, and keepPos flags set, expected 18 fields", resultRow)
  	}

  	vcfPos, err := strconv.Atoi(resultRow[len(resultRow) - 4])

		if err != nil || vcfPos != 1000 {
			t.Error("vcfPos should be the first optional field")
		}

		if resultRow[len(resultRow) - 3] != "rs1" {
			t.Error("id should come after vcfPos")
		}

		altIdx, err := strconv.Atoi(resultRow[len(resultRow) - 2])

		if err != nil || altIdx != 0 {
			t.Error("altIdx should come after vcfPos and after id")
		}

		if resultRow[len(resultRow) - 1] != "AC=1" {
			t.Error("INFO should come after vcfPos, id, and altIdx")
		}
	});
}

func TestOutputsSamplesVcfPosIdAndInfo(t *testing.T) {
	versionLine := "##fileformat=VCFv4.x"
	header := strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
		"FILTER", "INFO", "FORMAT", "Sample1", "Sample2", "Sample3", "Sample4"}, "\t")
	record := strings.Join([]string{"10", "1000", "rs123", "A", "T", "100", "PASS",
		"AC=1", "GT", "0/0", "0/1", "1/1", "./."}, "\t")

	lines := versionLine + "\n" + header + "\n" + record	+ "\n"

	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", keepPos: true, keepId: true, keepInfo: true}

  readVcf(&config, reader, func(row string) {
  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr10" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if resultRow[5] != "2" {
  		t.Error("Transversions should have value 2", resultRow)
  	}

  	if len(resultRow) != 17 {
  		t.Error("With both keepId and retainInfo flags set, should output 16 fields")
  	}

  	if resultRow[len(resultRow) - 4] != "1000" {
			t.Error("with keepPos, should get original input position as first optional field")
		}

		if resultRow[len(resultRow) - 3] != "rs123" {
			t.Error("Expect ID field after vcfPos")
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

		// heterozygosity
		// the denominator excludes the missing samples
		if resultRow[7] == strconv.FormatFloat(float64(1)/float64(3), 'G', 3, 64) {
			t.Log("OK: Recapitualte the heterozygosity", resultRow)
		} else {
			t.Error("NOT OK: Couldn't recapitualte the heterozygosity", resultRow)
		}

		if resultRow[8] == "Sample3" {
			t.Log("OK: Recapitualte the homozygote", resultRow)
		} else {
			t.Error("NOT OK: Couldn't recapitualte the homozygote", resultRow)
		}

		// homozygosity
		// the denominator excludes the missing samples
		if resultRow[9] == strconv.FormatFloat(float64(1)/float64(3), 'G', 3, 64) {
			t.Log("OK: Recapitualte the homozygosity", resultRow)
		} else {
			t.Error("NOT OK: Couldn't recapitualte the homozygosity", resultRow)
		}

		if resultRow[10] == "Sample4" {
			t.Log("OK: Recapitualte the missing sample", resultRow)
		} else {
			t.Error("NOT OK: Couldn't recapitualte the missing sample", resultRow)
		}

		// missingness
		// missingness denominator does not exclude missing samples
		if resultRow[11] == strconv.FormatFloat(float64(1)/float64(4), 'G', 3, 64) {
			t.Log("OK: Recapitualte the missingness", resultRow)
		} else {
			t.Error("NOT OK: Couldn't recapitualte the missingness", resultRow)
		}

		// sampleMaf ; 2 alleles for 1|1, 1 allele for 0|1 and 0 in denominator for .|.
		if resultRow[12] == strconv.FormatFloat(float64(3)/float64(6), 'G', 3, 64) {
			t.Log("OK: sampleMaf will count homozygotes, heterozygotes, and will exclude missing alleles from denominator", resultRow)
		} else {
			t.Error("NOT OK: sampleMaf will count homozygotes, heterozygotes, and will exclude missing alleles from denominator", resultRow)
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

  	if len(resultRow) != 17 {
  		t.Error("With keepPos, keepId, and keepInfo flags set, should output 17 fields", resultRow)
  	}

  	if resultRow[len(resultRow) - 4] != "1000" {
			t.Error("with keepPos, should get original input position as first optional field, for each allele in multialelic")
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

		// missingness
		// denominator does not exclude missingness
		if resultRow[11] == strconv.FormatFloat(float64(1)/float64(4), 'G', 3, 64) {
			t.Log("OK: Missingness is identical for each allele in multiallelic output", resultRow)
		} else {
			t.Error("NOT OK: Missingness is identical for each allele in multiallelic output", resultRow)
		}

		if index == 0 {
			// heterozygosity
			if resultRow[7] == "0" {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case for 1st allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case for 1st allele", resultRow)
			}

			// homozygosity
			// the denominator excludes the missing samples
			if resultRow[9] == strconv.FormatFloat(float64(1)/float64(3), 'G', 3, 64) {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case for 1st allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case for 1st allele", resultRow)
			}

			// sampleMaf ; 2 alleles for 1|1 and 0 in denominator for .|.
			if resultRow[12] == strconv.FormatFloat(float64(2)/float64(6), 'G', 3, 64) {
				t.Log("OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			}
		} else if index == 1 {
			// heterozygosity
			// the denominator excludes the missing samples
			if resultRow[7] == strconv.FormatFloat(float64(1)/float64(3), 'G', 3, 64) {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// homozygosity
			if resultRow[9] == "0" {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// sampleMaf ; 2 alleles for 0|2 and 0 in denominator for .|.
			if resultRow[12] == strconv.FormatFloat(float64(1)/float64(6), 'G', 3, 64) {
				t.Log("OK: sampleMaf in multialellic case for 2nd allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 2nd allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			}
		}
	})

	if index != 1 {
		t.Error("NOT OK: Expected to parse 2 alleles, parsed fewer")
	}

	// repeat but with "/"
	record = strings.Join([]string{"10", "1000", "rs456", "C", "T,G", "100", "PASS",
		"AC=1", "GT", "1/1", "0/0", "0/2", "./."}, "\t")

	lines = versionLine + "\n" + header + "\n" + record	+ "\n"
	reader = bufio.NewReader(strings.NewReader(lines))

	index = -1
  readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr10" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if len(resultRow) != 17 {
  		t.Error("With keepPos, keepId, and keepInfo flags set, should output 17 fields", resultRow)
  	}

  	if resultRow[len(resultRow) - 4] != "1000" {
			t.Error("with keepPos, should get original input position as first optional field")
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

		// missingness
		// denominator does not exclude missingness
		if resultRow[11] == strconv.FormatFloat(float64(1)/float64(4), 'G', 3, 64) {
			t.Log("OK: Missingness is identical for each allele in multiallelic output", resultRow)
		} else {
			t.Error("NOT OK: Missingness is identical for each allele in multiallelic output", resultRow)
		}

		if index == 0 {
			// heterozygosity
			if resultRow[7] == "0" {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case for 1st allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case for 1st allele", resultRow)
			}

			// homozygosity
			// the denominator excludes the missing samples
			if resultRow[9] == strconv.FormatFloat(float64(1)/float64(3), 'G', 3, 64) {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case for 1st allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case for 1st allele", resultRow)
			}

			// sampleMaf ; 2 alleles for 1|1 and 0 in denominator for .|.
			if resultRow[12] == strconv.FormatFloat(float64(2)/float64(6), 'G', 3, 64) {
				t.Log("OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			}
		} else if index == 1 {
			// heterozygosity
			// the denominator excludes the missing samples
			if resultRow[7] == strconv.FormatFloat(float64(1)/float64(3), 'G', 3, 64) {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// homozygosity
			if resultRow[9] == "0" {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// sampleMaf ; 2 alleles for 0|2 and 0 in denominator for .|.
			if resultRow[12] == strconv.FormatFloat(float64(1)/float64(6), 'G', 3, 64) {
				t.Log("OK: sampleMaf in multialellic case for 2nd allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 2nd allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			}
		}
	})

	if index != 1 {
		t.Error("NOT OK: Expected to parse 2 alleles, parsed fewer")
	}

	// repeat but with "|" and GT
	record = strings.Join([]string{"10", "1000", "rs456", "C", "T,G", "100", "PASS",
		"AC=1", "GT:GQ", "1|1:1,2,3", "0|0:4,5,6", "0|2:1,3,5", ".|.:0,0,0"}, "\t")

	lines = versionLine + "\n" + header + "\n" + record	+ "\n"
	reader = bufio.NewReader(strings.NewReader(lines))

	index = -1
  readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr10" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if len(resultRow) != 17 {
  		t.Error("With keepPos, keepId, and keepInfo flags set, should output 17 fields", resultRow)
  	}

  	if resultRow[len(resultRow) - 4] != "1000" {
			t.Error("with keepPos, should get original input position as first optional field")
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

		// missingness
		// denominator does not exclude missingness
		if resultRow[11] == strconv.FormatFloat(float64(1)/float64(4), 'G', 3, 64) {
			t.Log("OK: Missingness is identical for each allele in multiallelic output", resultRow)
		} else {
			t.Error("NOT OK: Missingness is identical for each allele in multiallelic output", resultRow)
		}

		if index == 0 {
			// heterozygosity
			if resultRow[7] == "0" {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case for 1st allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case for 1st allele", resultRow)
			}

			// homozygosity
			// the denominator excludes the missing samples
			if resultRow[9] == strconv.FormatFloat(float64(1)/float64(3), 'G', 3, 64) {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case for 1st allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case for 1st allele", resultRow)
			}

			// sampleMaf ; 2 alleles for 1|1 and 0 in denominator for .|.
			if resultRow[12] == strconv.FormatFloat(float64(2)/float64(6), 'G', 3, 64) {
				t.Log("OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			}
		} else if index == 1 {
			// heterozygosity
			// the denominator excludes the missing samples
			if resultRow[7] == strconv.FormatFloat(float64(1)/float64(3), 'G', 3, 64) {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// homozygosity
			if resultRow[9] == "0" {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// sampleMaf ; 2 alleles for 0|2 and 0 in denominator for .|.
			if resultRow[12] == strconv.FormatFloat(float64(1)/float64(6), 'G', 3, 64) {
				t.Log("OK: sampleMaf in multialellic case for 2nd allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 2nd allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			}
		}
	})

	if index != 1 {
		t.Error("NOT OK: Expected to parse 2 alleles, parsed fewer")
	}

	// repeat but with "/" and GT
	// Note, that we no longer expect missing records to have genotype quality strings
	record = strings.Join([]string{"10", "1000", "rs456", "C", "T,G", "100", "PASS",
		"AC=1", "GT:GQ", "1/1:1,2,3", "0/0:4,5,6", "0/2:1,3,5", "./."}, "\t")

	lines = versionLine + "\n" + header + "\n" + record	+ "\n"
	reader = bufio.NewReader(strings.NewReader(lines))

	index = -1
  readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr10" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if len(resultRow) != 17 {
  		t.Error("With keepPos, keepId, and keepInfo flags set, should output 17 fields", resultRow)
  	}

  	if resultRow[len(resultRow) - 4] != "1000" {
			t.Error("with keepPos, should get original input position as first optional field")
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

		// missingness
		// denominator does not exclude missingness
		if resultRow[11] == strconv.FormatFloat(float64(1)/float64(4), 'G', 3, 64) {
			t.Log("OK: Missingness is identical for each allele in multiallelic output", resultRow)
		} else {
			t.Error("NOT OK: Missingness is identical for each allele in multiallelic output", resultRow)
		}

		if index == 0 {
			// heterozygosity
			if resultRow[7] == "0" {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case for 1st allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case for 1st allele", resultRow)
			}

			// homozygosity
			// the denominator excludes the missing samples
			if resultRow[9] == strconv.FormatFloat(float64(1)/float64(3), 'G', 3, 64) {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case for 1st allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case for 1st allele", resultRow)
			}

			// sampleMaf ; 2 alleles for 1|1 and 0 in denominator for .|.
			if resultRow[12] == strconv.FormatFloat(float64(2)/float64(6), 'G', 3, 64) {
				t.Log("OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			}
		} else if index == 1 {
			// heterozygosity
			// the denominator excludes the missing samples
			if resultRow[7] == strconv.FormatFloat(float64(1)/float64(3), 'G', 3, 64) {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// homozygosity
			if resultRow[9] == "0" {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// missingness
			if resultRow[11] == strconv.FormatFloat(float64(1)/float64(4), 'G', 3, 64) {
				t.Log("OK: could recapitulate missingness even when the misssing genotype doesn't have the expected FORMAT data", resultRow)
			} else {
				t.Error("NOT OK: couldn't recapitulate missingness even when the misssing genotype doesn't have the expected FORMAT data", resultRow)
			}

			// sampleMaf ; 2 alleles for 0|2 and 0 in denominator for .|.
			if resultRow[12] == strconv.FormatFloat(float64(1)/float64(6), 'G', 3, 64) {
				t.Log("OK: sampleMaf in multialellic case for 2nd allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 2nd allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			}
		}
	})

	if index != 1 {
		t.Error("NOT OK: Expected to parse 2 alleles, parsed fewer")
	}

	// Test with no missing alleles
	header = strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
		"FILTER", "INFO", "FORMAT", "Sample1", "Sample2", "Sample3", "Sample4", "Sample5"}, "\t")
	record = strings.Join([]string{"15", "1001", "rs457", "C", "T,G", "100", "PASS",
		"AC=1", "GT", "0|1", "2|0", "2|2", "0|0", "1|0"}, "\t")

	lines = versionLine + "\n" + header + "\n" + record	+ "\n"

	reader = bufio.NewReader(strings.NewReader(lines))

	index = -1
  readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr15" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if len(resultRow) != 17 {
  		t.Error("With keepPos, keepId, and keepInfo flags set, should output 17 fields", resultRow)
  	}

  	if resultRow[len(resultRow) - 4] != "1001" {
			t.Error("with keepPos, should get original input position as first optional field")
		}

  	if resultRow[len(resultRow) - 3] == "rs457" {
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

		// missingness
		// denominator does not exclude missingness
		if resultRow[11] == "0" {
			t.Log("OK: Missingness is identical for each allele in multiallelic output", resultRow)
		} else {
			t.Error("NOT OK: Missingness is identical for each allele in multiallelic output", resultRow)
		}

		if index == 0 {
			// heterozygosity
			if resultRow[7] == strconv.FormatFloat(float64(2)/float64(5), 'G', 3, 64) {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case for 1st allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case for 1st allele", resultRow)
			}

			// homozygosity
			// the denominator excludes the missing samples
			if resultRow[9] == "0" {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case for 1st allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case for 1st allele", resultRow)
			}

			// sampleMaf ; 2 alleles for 1|1 and 0 in denominator for .|.
			if resultRow[12] == strconv.FormatFloat(float64(2)/float64(10), 'G', 3, 64) {
				t.Log("OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			}
		} else if index == 1 {
			// heterozygosity: 1 het for 2|0
			// the denominator excludes the missing samples
			if resultRow[7] == strconv.FormatFloat(float64(1)/float64(5), 'G', 3, 64) {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// homozygosity: 1 homozygote (2|2)
			if resultRow[9] == strconv.FormatFloat(float64(1)/float64(5), 'G', 3, 64) {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// sampleMaf ; 2 alleles for 2|2
			if resultRow[12] == strconv.FormatFloat(float64(3)/float64(10), 'G', 3, 64){
				t.Log("OK: sampleMaf in multialellic case for 2nd allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 2nd allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			}
		}
	})

	if index != 1 {
		t.Error("NOT OK: Expected to parse 2 alleles, parsed fewer")
	}

	// Test with no missing alleles, and "/" delimiter
	header = strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
		"FILTER", "INFO", "FORMAT", "Sample1", "Sample2", "Sample3", "Sample4", "Sample5"}, "\t")
	record = strings.Join([]string{"15", "1002", "rs457", "C", "T,G", "100", "PASS",
		"AC=2", "GT", "0/1", "2/0", "2/2", "0/0", "1/0"}, "\t")

	lines = versionLine + "\n" + header + "\n" + record	+ "\n"

	reader = bufio.NewReader(strings.NewReader(lines))

	index = -1
  readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr15" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if len(resultRow) != 17 {
  		t.Error("With keepPos, keepId, and keepInfo flags set, should output 17 fields", resultRow)
  	}

  	if resultRow[len(resultRow) - 4] != "1002" {
			t.Error("with keepPos, should get original input position as first optional field")
		}

  	if resultRow[len(resultRow) - 3] == "rs457" {
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

		if resultRow[len(resultRow) - 1] == "AC=2" {
			t.Log("OK: add INFO field correctly for multiallelic field in column 10")
		} else {
			t.Error("NOT OK: Couldn't add INFO field", resultRow)
		}

		// missingness
		// denominator does not exclude missingness
		if resultRow[11] == "0" {
			t.Log("OK: Missingness is identical for each allele in multiallelic output", resultRow)
		} else {
			t.Error("NOT OK: Missingness is identical for each allele in multiallelic output", resultRow)
		}

		if index == 0 {
			// heterozygosity
			if resultRow[7] == strconv.FormatFloat(float64(2)/float64(5), 'G', 3, 64) {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case for 1st allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case for 1st allele", resultRow)
			}

			// homozygosity
			// the denominator excludes the missing samples
			if resultRow[9] == "0" {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case for 1st allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case for 1st allele", resultRow)
			}

			// sampleMaf ; 2 alleles for 1|1 and 0 in denominator for .|.
			if resultRow[12] == strconv.FormatFloat(float64(2)/float64(10), 'G', 3, 64) {
				t.Log("OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			}
		} else if index == 1 {
			// heterozygosity: 1 het for 2|0
			// the denominator excludes the missing samples
			if resultRow[7] == strconv.FormatFloat(float64(1)/float64(5), 'G', 3, 64) {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// homozygosity: 1 homozygote (2|2)
			if resultRow[9] == strconv.FormatFloat(float64(1)/float64(5), 'G', 3, 64) {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// sampleMaf ; 2 alleles for 2|2
			if resultRow[12] == strconv.FormatFloat(float64(3)/float64(10), 'G', 3, 64){
				t.Log("OK: sampleMaf in multialellic case for 2nd allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 2nd allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			}
		}
	})

	if index != 1 {
		t.Error("NOT OK: Expected to parse 2 alleles, parsed fewer")
	}

	// Test with no missing alleles, "/" delimiter and stuff after the genotype
	header = strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
		"FILTER", "INFO", "FORMAT", "Sample1", "Sample2", "Sample3", "Sample4", "Sample5"}, "\t")
	record = strings.Join([]string{"15", "1001", "rs457", "C", "T,G", "100", "PASS",
		"AC=2", "GT:GL", "0/1:4,5,6", "2/0:7,8,9", "2/2:1,2,3", "0/0:.,.,.", "1/0:1,2,5"}, "\t")

	lines = versionLine + "\n" + header + "\n" + record	+ "\n"

	reader = bufio.NewReader(strings.NewReader(lines))

	index = -1
  readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr15" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if len(resultRow) != 17 {
  		t.Error("With keepPos, keepId, and keepInfo flags set, should output 17 fields", resultRow)
  	}

  	if resultRow[len(resultRow) - 4] != "1001" {
			t.Error("with keepPos, should get original input position as first optional field")
		}

  	if resultRow[len(resultRow) - 3] == "rs457" {
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

		if resultRow[len(resultRow) - 1] == "AC=2" {
			t.Log("OK: add INFO field correctly for multiallelic field in column 10")
		} else {
			t.Error("NOT OK: Couldn't add INFO field", resultRow)
		}

		// missingness
		// denominator does not exclude missingness
		if resultRow[11] == "0" {
			t.Log("OK: Missingness is identical for each allele in multiallelic output", resultRow)
		} else {
			t.Error("NOT OK: Missingness is identical for each allele in multiallelic output", resultRow)
		}

		if index == 0 {
			// heterozygosity
			if resultRow[7] == strconv.FormatFloat(float64(2)/float64(5), 'G', 3, 64) {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case for 1st allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case for 1st allele", resultRow)
			}

			// homozygosity
			// the denominator excludes the missing samples
			if resultRow[9] == "0" {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case for 1st allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case for 1st allele", resultRow)
			}

			// sampleMaf ; 2 alleles for 1|1 and 0 in denominator for .|.
			if resultRow[12] == strconv.FormatFloat(float64(2)/float64(10), 'G', 3, 64) {
				t.Log("OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			}
		} else if index == 1 {
			// heterozygosity: 1 het for 2|0
			// the denominator excludes the missing samples
			if resultRow[7] == strconv.FormatFloat(float64(1)/float64(5), 'G', 3, 64) {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// homozygosity: 1 homozygote (2|2)
			if resultRow[9] == strconv.FormatFloat(float64(1)/float64(5), 'G', 3, 64) {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// sampleMaf ; 2 alleles for 2|2
			if resultRow[12] == strconv.FormatFloat(float64(3)/float64(10), 'G', 3, 64){
				t.Log("OK: sampleMaf in multialellic case for 2nd allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 2nd allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			}
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

	allowedFilters := map[string]bool{ "PASS": true, ".": true}

	lines := versionLine + "\n" + header + "\n" + record	+ "\n"
	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", allowedFilters: allowedFilters}

	index := -1
  readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr20" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	if resultRow[5] != "0" {
  		t.Error("Multiallelics should be called neither transitions nor transversions with value 0", resultRow)
  	}

  	if index == 0 {
  		if resultRow[altIdx] == "-4" && resultRow[posIdx] == "5" && resultRow[refIdx] == "C" {
  			t.Log("OK: 4 base deletion recapitulated", resultRow)
  		} else {
  			t.Error("NOT OK: 4 base deletion not recapitulated", resultRow)
  		}

			if resultRow[6] == "Sample2" {
				t.Log("OK: Recapitualte 1st allele het", resultRow)
  		} else {
  			t.Error("NOT OK: Couldn't recapitualte 1st allele het", resultRow)
  		}

  		if resultRow[8] == config.emptyField {
				t.Log("OK: Recapitualte 1st allele has no homozygotes", resultRow)
  		} else {
  			t.Error("NOT OK: Couldn't recapitualte 1st allele has no homozygotes", resultRow)
  		}

  		if resultRow[10] == "Sample4" {
				t.Log("OK: Recapitualte 1st allele missing for Sample4", resultRow)
  		} else {
  			t.Error("NOT OK: Couldn't recapitualte 1st allele missing for Sample4", resultRow)
  		}

  		// heterozygosity
  		// the denominator excludes the missing samples
			if resultRow[7] == strconv.FormatFloat(float64(1)/float64(3), 'G', 3, 64) {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case", resultRow)
			}

			// homozygosity
			if resultRow[9] == "0" {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case", resultRow)
			}

			// missingness
			// the denominator of missingness does not exclude the missing samples
			if resultRow[11] == strconv.FormatFloat(float64(1)/float64(4), 'G', 3, 64) {
				t.Log("OK: Recapitualte the missingness in multialellic case", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the missingness in multiallelic case", resultRow)
			}

			// sampleMaf ; 2 alleles for 0|1 and 0 in denominator for .|.
			if resultRow[12] == strconv.FormatFloat(float64(1)/float64(6), 'G', 3, 64) {
				t.Log("OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 1st allele will count in denominator all multiallelic alleles, but will exclude missing alleles", resultRow)
			}
		}

		if index == 1 {
  		if resultRow[altIdx] == "+TCA" && resultRow[posIdx] == "4" && resultRow[refIdx] == "G" {
  			t.Log("OK: 3 base insertion recapitulated", resultRow)
  		} else {
  			t.Error("NOT OK: 3 base insertion not recapitulated", resultRow)
  		}

  		if resultRow[6] == config.emptyField {
				t.Log("OK: Recapitualte 2nd allele het", resultRow)
  		} else {
  			t.Error("NOT OK: Couldn't recapitualte 2nd allele het", resultRow)
  		}

  		if resultRow[8] == "Sample3" {
				t.Log("OK: Recapitualte 2nd allele has Sample3 as homozygous", resultRow)
  		} else {
  			t.Error("NOT OK: Couldn't recapitualte 2nd allele has Sample3 as homozygous", resultRow)
  		}

  		if resultRow[10] == "Sample4" {
				t.Log("OK: Recapitualte 2nd allele is missing for Sample4", resultRow)
  		} else {
  			t.Error("NOT OK: Couldn't recapitualte 2nd allele missing for Sample4", resultRow)
  		}

  		// heterozygosity
			if resultRow[7] == "0" {
				t.Log("OK: Recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the heterozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// homozygosity
			// the denominator excludes the missing samples
			if resultRow[9] == strconv.FormatFloat(float64(1)/float64(3), 'G', 3, 64) {
				t.Log("OK: Recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the homozygosity in multiallelic case for 2nd allele", resultRow)
			}

			// missingness
			// the denominator of missingness does not exclude the missing samples
			if resultRow[11] == strconv.FormatFloat(float64(1)/float64(4), 'G', 3, 64) {
				t.Log("OK: Recapitualte the missingness in multialellic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: Couldn't recapitualte the missingness in multiallelic case for 2nd allele", resultRow)
			}

			// sampleMaf ; 2 alleles for 2|2 and 0 in denominator for .|.
			if resultRow[12] == strconv.FormatFloat(float64(2)/float64(6), 'G', 3, 64) {
				t.Log("OK: sampleMaf in multialellic case for 2nd allele", resultRow)
			} else {
				t.Error("NOT OK: sampleMaf in multialellic case for 2nd allele", resultRow)
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

	allowedFilters := map[string]bool{ "PASS": true, ".": true}

	lines := versionLine + "\n" + header + "\n" + record	+ "\n"
	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", allowedFilters: allowedFilters}

	index := -1;

	readVcf(&config, reader, func(row string) {
  	index++

  	resultRow := strings.Split(row[:len(row)-1], "\t")

  	if resultRow[0] != "chr16" {
  		t.Error("chromosome should have chr appended", resultRow)
  	}

  	// heterozygosity
		if resultRow[9] == "0" {
			t.Log("OK: When no samples provided, heterozygosity is 0", resultRow)
		} else {
			t.Error("NOT OK: When no samples provided, heterozygosity is not 0", resultRow)
		}

		// homozygosity
		if resultRow[9] == "0" {
			t.Log("OK: When no samples provided, homozygosity is 0", resultRow)
		} else {
			t.Error("NOT OK: When no samples provided, homozygosity is not 0", resultRow)
		}

		// missingness
		if resultRow[11] == "0" {
			t.Log("OK: When no samples provided, missingness is 0", resultRow)
		} else {
			t.Error("NOT OK: When no samples provided, missingness is not 0", resultRow)
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

	allowedFilters := map[string]bool{ "PASS": true, ".": true}

	lines := versionLine + "\n" + header + "\n" + record	+ "\n"
	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", allowedFilters: allowedFilters}

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

func TestOutputMultiallelicSnp(t *testing.T) {
  versionLine := "##fileformat=VCFv4.x"
  header := strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
    "FILTER", "INFO"}, "\t")
  //Test example 5.2.4 https://samtools.github.io/hts-specs/VCFv4.2.pdf
  //Test from 1000 genomes
  record := strings.Join([]string{"1", "1265061", "rs138351882;rs563042459",
    "CGT", "TGT,C", ".", "PASS", "DP=100"}, "\t")

  allowedFilters := map[string]bool{ "PASS": true, ".": true}

  lines := versionLine + "\n" + header + "\n" + record  + "\n"
  reader := bufio.NewReader(strings.NewReader(lines))

  config := Config{emptyField: "!", fieldDelimiter: ";", allowedFilters: allowedFilters}

  index := -1;
  readVcf(&config, reader, func(row string) {
    index++

    resultRow := strings.Split(row[:len(row)-1], "\t")

    if resultRow[0] != "chr1" {
      t.Error("chromosome should have chr appended", resultRow)
    }

    fmt.Println(index, resultRow)
    if index == 0 {
      if resultRow[refIdx] == "C" && resultRow[posIdx] == "1265061" && resultRow[altIdx] == "T" {
        t.Log("OK: Complex SNP in multiallelic SNP/DEL annotated correctly", resultRow)
      } else {
        t.Error("NOT OK: Complex SNP in multiallelic SNP/DEL annotated correctly", resultRow)
      }
    }

    if index == 1 {
      if resultRow[refIdx] == "G" && resultRow[posIdx] == "1265062" && resultRow[altIdx] == "-2" {
          t.Log("OK: DEL in multiallelic SNP/DEL annotated correctly", resultRow)
      } else {
        t.Error("NOT OK: DEL in multiallelic SNP/DEL not annotated correctly", resultRow)
      }
    }
  })

  if index == 1 {
    t.Log("Ok, recapitulated 4 alleles")
  } else {
    t.Error("expected 4 alleles")
  }
}

func TestComplexSnp(t *testing.T) {
  versionLine := "##fileformat=VCFv4.x"
  header := strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
    "FILTER", "INFO"}, "\t")
  //Test example 5.2.4 https://samtools.github.io/hts-specs/VCFv4.2.pdf
  //Hypothetical, where a silly site lists the T->C 2 bases before it occurs
  record := strings.Join([]string{"1", "1265062", "rs138351882;rs563042459",
    "CGT", "CGA", ".", "PASS", "DP=100"}, "\t")

  allowedFilters := map[string]bool{ "PASS": true, ".": true}

  lines := versionLine + "\n" + header + "\n" + record  + "\n"
  reader := bufio.NewReader(strings.NewReader(lines))

  config := Config{emptyField: "!", fieldDelimiter: ";", allowedFilters: allowedFilters}

  index := -1;
  readVcf(&config, reader, func(row string) {
    index++

    resultRow := strings.Split(row[:len(row)-1], "\t")

    if resultRow[0] != "chr1" {
      t.Error("chromosome should have chr appended", resultRow)
    }

    if index == 0 {
      if resultRow[refIdx] == "T" && resultRow[posIdx] == "1265064" && resultRow[altIdx] == "A" {
        t.Log("OK: Complex SNP annotated correctly", resultRow)
      } else {
        t.Error("NOT OK: Complex SNP not annotated correctly", resultRow)
      }
    }
  })

  if index == 0 {
    t.Log("Ok, recapitulated 1 alleles")
  } else {
    t.Error("expected 1 alleles")
  }
}

func TestMNP(t *testing.T) {
  versionLine := "##fileformat=VCFv4.x"
  header := strings.Join([]string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
    "FILTER", "INFO"}, "\t")
  //Hypothetical 
  record := strings.Join([]string{"1", "1000", "rs138351882;rs563042459",
    "ACGT", "GATC", ".", "PASS", "DP=100"}, "\t")

  test := "Finds MNPs"

  allowedFilters := map[string]bool{ "PASS": true, ".": true}

  lines := versionLine + "\n" + header + "\n" + record  + "\n"
  reader := bufio.NewReader(strings.NewReader(lines))

  config := Config{emptyField: "!", fieldDelimiter: ";", allowedFilters: allowedFilters}

  index := -1;
  readVcf(&config, reader, func(row string) {
    index++

    resultRow := strings.Split(row[:len(row)-1], "\t")

    if resultRow[0] != "chr1" {
      t.Error("chromosome should have chr appended", resultRow)
    }

    if index == 0 {
      if resultRow[refIdx] == "A" && resultRow[posIdx] == "1000" && resultRow[altIdx] == "G" {
        t.Logf("OK: %s", test)
      } else {
        t.Errorf("NOT OK: %s", test, resultRow)
      }
    }

    if index == 1 {
      if resultRow[refIdx] == "C" && resultRow[posIdx] == "1001" && resultRow[altIdx] == "A" {
        t.Logf("OK: %s", test)
      } else {
        t.Errorf("NOT OK: %s", test, resultRow)
      }
    }

    if index == 2 {
      if resultRow[refIdx] == "G" && resultRow[posIdx] == "1002" && resultRow[altIdx] == "T" {
        t.Logf("OK: %s", test)
      } else {
        t.Errorf("NOT OK: %s", test, resultRow)
      }
    }

    if index == 3 {
      if resultRow[refIdx] == "T" && resultRow[posIdx] == "1003" && resultRow[altIdx] == "C" {
        t.Logf("OK: %s", test)
      } else {
        t.Errorf("NOT OK: %s", test, resultRow)
      }
    }
  })

  if index == 3 {
    t.Logf("OK: %s", test)
  } else {
    t.Errorf("NOT OK: %s", test)
  }
}