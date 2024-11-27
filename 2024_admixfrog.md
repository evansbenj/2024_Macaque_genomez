# Admixfrog

Admixfrog is a program written by Ben Peter that uses an HMM to estimate locations of introgression blocks in genomic data.

I installed the latest version (admixfrog 0.7.2.post1.dev0+1c27c43) like this:
```
module load StdEnv/2023 python/3.12.4
pip install cython scipy --upgrade --user
pip install git+https://github.com/benjaminpeter/admixfrog --user
```

and tested like this:

```
/home/ben/.local/bin/admixfrog --help
```

# Filtering
I first removed hard filtered positions using bcftools:
```
#!/bin/sh
#SBATCH --job-name=bcftools_filter
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=2gb
#SBATCH --output=bcftools_filter.%J.out
#SBATCH --error=bcftools_filter.%J.err
#SBATCH --account=rrg-ben

# execute like this: sbatch 2024_bcftools_concat.sh output.filtered.snps.1.*.AB.vep.vcf.gz_snpsonly.vcf.gz
# load these modules before running:
module load StdEnv/2023  gcc/12.3 bcftools/1.19 gatk/4.4.0.0 java/17.0.6 tabix

bcftools view -f PASS ${1} -O z -o ${1}_filtered.vcf.gz
bcftools index ${1}_filtered.vcf.gz
gatk --java-options -Xmx10G IndexFeatureFile -I ${1}_filtered.vcf.gz
```
Then I filtered to include only sites with the number of missing genotypes being less than or equal to 10 and setting the minimum genotype quality at 30 as follows:
```
module load StdEnv/2023 vcftools
```
```
vcftools --gzvcf output.filtered.snps.1.all.AB.vep.vcf.gz_snpsonly.vcf.gz_concat.vcf.gz_Sulasnpsonly.vcf.gz_filtered.vcf.gz --max-missing-count 10 --minQ 30 --recode --recode-INFO-all --stdout | gzip -c >  ./SulaSNPs.chr1_maxmissingcount_10_genoqual30.vcf.gz
```
