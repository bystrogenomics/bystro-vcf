# bystro-vcf [![Build Status](https://travis-ci.org/akotlar/bystro-vcf.svg?branch=master)](https://travis-ci.org/akotlar/bystro-vcf)

<br>

## TL:DR

A really fast, simple VCF pre-processor and annotator.

Performs several important functions:
1) Splits multiallelics and MNP alleles, keeping track of each allele's index with respect to the original alleles for downstream INFO property segregation
2) Performs QC on variants: checks whether allele contains ACTG, that padding bases match reference, and more
3) Allows filtering of variants by any number of FILTER properties (by default allows PASS/. variants)
4) Normalizes indel representations by removing padding, left shifting alleles to their parsimonious representations
5) Calculates whether site is transition, transversion, or neither
6) Processes all available samples
    - calculates homozygosity, heterozygosity, missingness
    - labels samples as homozygous, heterozygous, or missing

<br>

## Publication

bystro-vcf is used to pre-proces VCF files for [Bystro](https://bystro.io) ([github](https://github.com/akotlar/bystro))

If you use bystro-vcf please cite https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1387-3 

<br>

## Performance
Millions of variants/rows per minute. Performance is dependent on the # of samples.

Amazon i3.2xlarge (4 core), 1K Genomes Phase 3 (2,504 samples): chromosome 1 (6.2M variants) in 5m30s.

<br>

## Installation

```shell
go get github.com/akotlar/bystro-vcf && go install $_;

pigz -d -c /some/vcf.gz | bystro-vcf --keepId --keepInfo --allowFilter "PASS,." 
```

<br>

## Use

Piping:
```shell
pigz -d -c /some/vcf.gz | bystro-vcf --keepId --keepInfo --allowFilter "PASS,." | pigz -c - > /path/to/out.gz
```

Input argument
```shell
bystro-vcf --inPath /path/to/vcf --keepId --keepInfo --allowFilter "PASS,." > /path/to/out
```

<br>

## Output
```tsv
chrom <String>   pos <Int>   type <String[SNP|DEL|INS|MULTIALLELIC]>    ref <String>    alt <String>    trTv <Int[0|1|2]>     heterozygotes <String>     heterozygosity <Float64>    homozygotes <String>     homozygosity <Float64>     missingGenos <String>    missingness <Float64>    sampleMaf <Float64>    id <String?>    alleleIndex <Int?>   info <String?>
```

<br>

## Optional arguments

```shell
--keepId <Bool>
```

Retain the "ID" field in the output.

<br>

```shell
--keepInfo <Bool>
```

Retain the "INFO" field in the output. 
  - Since we decompose multiallelics, an "alleleIdx" field is added to the output. It contains the 0-based index of that allele in the multiallelic
  - This is necessary for downstream programs to decompose the INFO field per-allele


Results in 2 output fields, following `missingGenos` or `id` should `--keepId` be set
  1. `alleleIdx` will contain the index of allele in a split multiallelic. 0 by default.
  2. `info` will contain the entire `INFO` string

<br>

```shell
--allowFilter <String>
```

Which `FILTER` values to keep. Comma separated. Defaults to `"PASS,."`
- Similar to [https://samtools.github.io/bcftools/bcftools.html](bcftools) `-f, --apply-filters LIST`

<br>

```shell
--excludeFilter <String>
```

Which `FILTER` values to exclude. Comma separated. Defaults to `""`
- Opposite of [https://samtools.github.io/bcftools/bcftools.html](bcftools) `-f, --apply-filters LIST`

<br>

```shell
--inPath /path/to/uncompressedFile.vcf
```

An input file path, to an uncompressed VCF file. Defaults to `stdin`

<br>

```shell
--errPath /path/to/log.txt
```

Where to store log messages. Defaults to `stderr`

<br>

```shell
--emptyField "!"
```

Which value to assign to missing data. Defaults to `!`

<br>

```shell
--fieldDelimiter ";"
```

Which delimiter to use when joining multiple values. Defaults to `;`
