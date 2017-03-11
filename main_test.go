package main

import "testing"

func TestUpdateFieldsWithAlt(t *testing.T) {
	expType, exPos, expRef, expAlt := "SNP", "100", "T", "C"

	sType, pos, ref, alt, err := updateFieldsWithAlt("T", "C", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("Test failed")
	} else {
		t.Log("OK: SNP")
	}

	expType, exPos, expRef, expAlt = "SNP", "103", "T", "A"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TCCT", "TCCA", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: SNPs that are longer than 1 base are suported when SNP at end")
	}

	expType, exPos, expRef, expAlt = "SNP", "102", "C", "A"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TGCT", "TGAT", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: SNPs that are longer than 1 base are suported when SNP in middle")
	}

	expType, exPos, expRef, expAlt = "SNP", "100", "T", "A"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TGCT", "AGCT", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: SNPs that are longer than 1 base are suported when SNP at beginning")
	}

	expType, exPos, expRef, expAlt = "", "", "", ""

	sType, pos, ref, alt, err = updateFieldsWithAlt("TCCT", "GTAA", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: MNPs are not supported")
	}

	expType, exPos, expRef, expAlt = "DEL", "101", "C", "-1"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TC", "T", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: 1-based deletions ")
	}

	expType, exPos, expRef, expAlt = "DEL", "101", "A", "-5"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TAGCGT", "T", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Deletions with references longer than 2 bases")
	}

	expType, exPos, expRef, expAlt = "DEL", "102", "G", "-4"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TAGCTT", "TA", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Deletions longer than 1 base")
	}

	expType, exPos, expRef, expAlt = "", "", "", ""

	sType, pos, ref, alt, err = updateFieldsWithAlt("TAGCTT", "TAC", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Left edge of deletion must match exactly")
	}

	expType, exPos, expRef, expAlt = "INS", "100", "T", "+AGCTT"

	sType, pos, ref, alt, err = updateFieldsWithAlt("T", "TAGCTT", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Insertions where reference is 1 base long")
	}

	expType, exPos, expRef, expAlt = "INS", "101", "A", "+GCTT"

	sType, pos, ref, alt, err = updateFieldsWithAlt("TA", "TAGCTT", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Insertions where reference is 2 bases long")
	}

	expType, exPos, expRef, expAlt = "", "", "", ""

	sType, pos, ref, alt, err = updateFieldsWithAlt("TT", "TAGCTT", "100")

	if err != nil || sType != expType || pos != exPos || ref != expRef || alt != expAlt {
		t.Error("Test failed", sType, pos, ref, alt)
	} else {
		t.Log("OK: Insertions where reference and alt don't share a left edge are skipped")
	}
}

func TestLineIsValid(t *testing.T) {
	expect := true

	actual := lineIsValid("ACTG")

	if expect != actual {
		t.Error()
	} else {
		t.Log("Support ACTG-containing alleles")
	}

	expect = false
	actual = lineIsValid(".")

	if expect != actual {
		t.Error("Can't handle missing Alt alleles")
	} else {
		t.Log("OK: Handles missing Alt alleles")
	}

	expect = false
	actual = lineIsValid("]13 : 123456]T")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("Can't handle single breakends")
	} else {
		t.Log("OK: Handles single breakend ']13 : 123456]T'")
	}

	expect = false
	actual = lineIsValid("C[2 : 321682[")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("Can't handle single breakends")
	} else {
		t.Log("OK: Handles single breakend 'C[2 : 321682['")
	}

	expect = false
	actual = lineIsValid(".A")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("Can't handle single breakends")
	} else {
		t.Log("OK: Handles single breakend '.A'")
	}

	expect = false
	actual = lineIsValid("G.")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("Can't handle single breakends")
	} else {
		t.Log("OK: Handles single breakend 'G.'")
	}

	expect = false
	actual = lineIsValid("<DUP>")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("Can't handle complex tags")
	} else {
		t.Log("OK: Handles complex Alt tags '<DUP>'")
	}

	expect = false
	actual = lineIsValid("A,C")

	// https://samtools.github.io/hts-specs/VCFv4.1.pdf
	if expect != actual {
		t.Error("Allows multiallelics", actual)
	} else {
		t.Log("OK: multiallelics are not supported. lineIsValid() requires multiallelics to be split")
	}
}
