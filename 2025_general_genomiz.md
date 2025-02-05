# General Genomics

After hardfiltering, I removed positions with missing genotypes:
```
vcftools --gzvcf all_162_maqs_chrY.vcf_malez_only.vcf.gz --max-missing-count 0 --minQ 30 --recode --recode-INFO-all --out all_162_maqs_chrY.vcf_malez_maxmissingcount_0_genoqual30.vcf
```

Then I converted these files to geno format (on info2020):
```
python3 /home/ben/2025_genomics_general/genomics_general/VCF_processing/parseVCF.py -i all_162_maqs_chrX_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=5 max=100 -o all_162_maqs_chrX_maxmissingcount_0_genoqual30.geno.gz
```
then I removed positions with "N"s like this:
```
zcat all_162_maqs_chr1_maxmissingcount_0_genoqual30.geno.gz | sed '/N\/N/d' > all_162_maqs_chr1_maxmissingcount_0_genoqual30_a.geno

cat all_162_maqs_chr1_maxmissingcount_0_genoqual30_a.geno | sed '/N|N/d' > all_162_maqs_chr1_maxmissingcount_0_genoqual30_b.geno

```
It will be interesting to compare 
* Fst of autosomes and the X (expect higher for the X)
* Fst of the X of females and the X of males (expect higher for female X)
* Fst of the X of males and the Y of males (expect higher for X)
