# 8_17_18 commit c44b19278bc2d913f4b488fb9d5658a87cc7dc79

## Test (sort multiallelics consistently, using -k5,5)

```sh
pigz -d -c ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.20klines.vcf.gz | bystro-vcf | pigz -c - >  out_check_new_8_17_18.vcf.gz

pigz -d -c out_check_new_5_26_18.vcf.gz | sort -k1,1 -k2,2n -k5,5 | pigz -c > out_check_new_5_26_18.vcf.sorted.gz

pigz -d -c out_check_new_8_17_18.vcf.gz | sort -k1,1 -k2,2n -k5,5 | pigz -c > out_check_new_8_17_18.vcf.sorted.gz

diff <(pigz -d -c out_check_new_8_17_18.vcf.sorted.gz) <(pigz -d -c out_check_new_5_26_18.vcf.sorted.gz)

# Nothing

pigz -d -c out_check_old_5_26_18.vcf.gz | sort -k1,1 -k2,2n -k5,5 | pigz -c - > out_check_old_5_26_18.vcf.sorted.gz

diff <(pigz -d -c out_check_new_8_17_18.vcf.sorted.gz) <(pigz -d -c out_check_old_5_26_18.vcf.sorted.gz)
# Nothing
```
