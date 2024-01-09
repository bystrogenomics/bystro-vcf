package main

import (
	"testing"
)

// Benchmark for comparing a rune element with a rune variable
func BenchmarkRuneEquality(b *testing.B) {
	// Sample rune slice and a rune variable to compare with
	data := []rune{'a', 'b', 'c'}
	compareRune := 'a'

	b.ResetTimer() // Reset the timer to exclude setup time
	for i := 0; i < b.N; i++ {
		for _, r := range data {
			if r == compareRune {
				// Perform some action if equal
			}
		}
	}
}

// Benchmark for comparing a rune element cast to a string with a string variable
func BenchmarkStringEquality(b *testing.B) {
	// Sample rune slice and a string variable to compare with
	data := []rune{'a', 'b', 'c'}
	compareString := "a"

	b.ResetTimer() // Reset the timer to exclude setup time
	for i := 0; i < b.N; i++ {
		for _, r := range data {
			if string(r) == compareString {
				// Perform some action if equal
			}
		}
	}
}

func BenchmarkStringToRuneEquality(b *testing.B) {
	// Sample rune slice and a string variable to compare with
	data := []string{"a", "b", "c"}
	compareRune := 'a'

	b.ResetTimer() // Reset the timer to exclude setup time
	for i := 0; i < b.N; i++ {
		for _, r := range data {
			if rune(r[0]) == compareRune {
				// Perform some action if equal
			}
		}
	}
}

func BenchmarkStringToStringEquality(b *testing.B) {
	// Sample rune slice and a string variable to compare with
	data := []string{"a", "b", "c"}
	compareString := "a"

	b.ResetTimer() // Reset the timer to exclude setup time
	for i := 0; i < b.N; i++ {
		for _, r := range data {
			if string(r[0]) == compareString {
				// Perform some action if equal
			}
		}
	}
}

func BenchmarkGenotypeToStringEquality(b *testing.B) {
	// Sample rune slice and a string variable to compare with
	data := []string{"0|1", "1|0", "0|1"}
	compareString := "0"

	b.ResetTimer() // Reset the timer to exclude setup time
	for i := 0; i < b.N; i++ {
		for _, r := range data {
			if string(r[0]) == compareString {
				// Perform some action if equal
			}
		}
	}
}

func BenchmarkGenotypeToRuneEquality(b *testing.B) {
	// Sample rune slice and a string variable to compare with
	data := []string{"0|1", "1|0", "0|1"}
	compareString := "0"

	b.ResetTimer() // Reset the timer to exclude setup time
	for i := 0; i < b.N; i++ {
		for _, r := range data {
			if rune(r[0]) == rune(compareString[0]) {
				// Perform some action if equal
			}
		}
	}
}

func BenchmarkGenotypeToRuneEqualityWithLengthTest(b *testing.B) {
	data := []string{"0|1", "1|0", "0|1"}
	compareString := "0"

	b.ResetTimer() // Reset the timer to exclude setup time
	for i := 0; i < b.N; i++ {
		if len(compareString) == 1 {
			for _, r := range data {
				if rune(r[0]) == rune(compareString[0]) {
					// Perform some action if equal
				}
			}
		}
	}
}
