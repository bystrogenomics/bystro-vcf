# bystro-vcf [![Build Status](https://travis-ci.org/akotlar/bystro-vcf.svg?branch=master)](https://travis-ci.org/akotlar/bystro-vcf)



## TL:DR

A really fast, simple VCF pre-processor and annotator. Takes a VCF file, splits multiallelics, normalizes indel representations by removing padding. 

```shell
go get github.com/akotlar/bystro-vcf && go install $_;

pigz -d -c /some/vcf.gz | bystro-vcf --keepId --keepInfo --allowFilter "PASS,."
```

## Output
```tsv
chrom <String>   pos <Int>   type <String[SNP|DEL|INS|MULTIALLELIC]>    ref <String>    alt <String>    trTv <Int[0|1|2]>     heterozygotes <String>     heterozygosity <Float64>    homozygotes <String>     homozygosity <Float64>     missingGenos <String>    missingness <Float64>    sampleMaf <Float64>    id <String?>    alleleIndex <Int?>   info <String?>
```

## Options

```shell
--keepId <Bool>
```
**Optional**. Keep delimiter to use when joining multiple values in the output of one field. Defaults to `;`

Will add the `id` after missingGenos

<br/>

```shell
--keepInfo <Bool>
```

**Optional**. What delimiter to use when joining multiple values in the output of one field. Defaults to `;`

Will add two fields: after either `missingGenos`, or `id` should `--keepId` be set
  1. `alleleIdx` will contain the index of allele in a split multiallelic. 0 by default.
  2. `info` will contain the entire `INFO` string

<br/>


```shell
--allowFilter <String>
```

**Optional**. Which `FILTER` values to keep. Comma separated. Defaults to `"PASS,."`
- Similar to [https://samtools.github.io/bcftools/bcftools.html](bcftools) `-f, --apply-filters LIST`

<br/>

```shell
--excludeFilter <String>
```

**Optional**. Which `FILTER` values to exclude. Comma separated. Defaults to `""`
- Opposite of [https://samtools.github.io/bcftools/bcftools.html](bcftools) `-f, --apply-filters LIST`

<br/>

```shell
--inPath /path/to/uncompressedFile.vcf
```

**Optional**. An input file path, to an uncompressed VCF file. Defaults to `stdin`

<br/>

```shell
--errPath /path/to/log.txt
```

**Optional**. Where to store log messages. Defaults to `stderr`

<br/>

```shell
--emptyField "!"
```
**Optional**. What value to give to missing data. Defaults to `!`

<br/>

```shell
--fieldDelimiter ";"
```
**Optional**. What delimiter to use when joining multiple values in the output of one field. Defaults to `;`
