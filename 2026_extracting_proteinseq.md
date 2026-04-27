# Extract protein seq from vcf

This worked the best because the variants are annotated by Alan using Ensembl VEP:
```
zgrep '|ENSMMUG00000030317' mac_chr5_concat.vcf.gz_filtered.vcf.gz_all160_snpsonly.vcf.gz | grep 'missense_variant' > chr5_MRPL1_missense_positions.txt
```

# Below did not work as well:
```
zgrep -E "ID=gene:ENSMMUG00000030317|Parent=gene:ENSMMUG00000030317|transcript:ENSMMUT00000042226" ../../../rheMac10/Macaca_mulatta.Mmul_10.115.gff3.gz > MRPL1.gff3
```


# Below could be useful but never worked:
```
module load apptainer/1.4.5    StdEnv/2023
singularity run agat_1.4.2--pl5321hdfd78af_0.sif
```
# Convert gtf to gff3 (not used)
```
agat_convert_sp_gxf2gxf.pl --gtf ../2021_rhemac_v10/rheMac10.refGene.gtf -o rheMac10.refGene.gff3
```
# Extract gff3 from one gene from a genomewide gff3 file:
```
agat_sp_extract_attributes.pl --gff ../../../rheMac10/Macaca_mulatta.Mmul_10.115.gff3.gz -t gene --attribute ENSMMUG00000030317 -o MRPL1.gff3
```
