# General Genomics

Convert the hard filtered vcf files to geno format (on info2020):
```
python3 /home/ben/2025_genomics_general/genomics_general/VCF_processing/parseVCF.py -i all_162_maqs_chr1.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=5 max=100 -o all_162_maqs_chr1.geno.gz
```

It will be interesting to compare 
* Fst of autosomes and the X (expect higher for the X)
* Fst of the X of females and the X of males (expect higher for female X)
* Fst of the X of males and the Y of males (expect higher for X)
