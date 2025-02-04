# VCF File Comparison

This R file can be used to compare two vcf files to understand the difference in terms of number of SNPs and Indels.

The script creates a txt file that contains the differences and similarities in the two vcf files. Then, following logic is used to compute the error (i.e.discrepencies)
 
```NP discrepancies = Total SNPs in combined dataset−(SNPs in Sample 1+SNPs in Sample 2)```
 
```Indel discrepancies = Total Indels in combined dataset−(Indels in Sample 1+Indels in Sample 2)```
 
```Total discrepancies = sum of SNP and indel discrepancies```

To run the code, you can use ```Rscript vcf_file_comparison.R``` in linux terminal. The file names and the ```common_data_path``` should be adjusted to input two vcf files to get final results.
