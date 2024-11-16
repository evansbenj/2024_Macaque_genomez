# Combining GVCFs

Alan and Ravi generated gvcfs for all the samples, which is amazing!  They called all callable SNPs, which will allow me to estimate pi per site.  Here is the information Alan provided:
```
I generated the GATK intervals by splitting the rhesus assembly on assembly gaps. That means there are gaps between intervals, but those are assembly gaps consisting on Ns where SNVs cannot be called.

Below are UCSC views of the two gaps you pointed out below showing the assembly gaps: 

https://genome.ucsc.edu/s/Rharris1/rheMac10.gap1

https://genome.ucsc.edu/s/Rharris1/rheMac10.gap2

I completed the all-sites GenotypeGVCFs SNV calling on the Sulawesi samples. I used GATK VariantAnnotator for AlleleBalance annotations and Ensembl VEP for gene annotations.

We apply the GATK hard filters described here:

https://software.broadinstitute.org/gatk/documentation/article?id=11097

Normally we remove SNVs that fail the hard filters which are marked with "hardFilter.snp" in the FILTER field. However, since you want all sites vcfs I left those in the vcfs.

```

I also need to work with only SNPs, so I am generating a concatenated SNPs only files for each chromosome like this.

Working directory:
```
/home/ben/projects/rrg-ben/ben/2024_macaques/concatenated_allsites_vcf/
```
# First extract SNPs from allsites files:
```
#!/bin/sh
#SBATCH --job-name=bcftools
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=2gb
#SBATCH --output=bcftools.%J.out
#SBATCH --error=bcftools.%J.err
#SBATCH --account=rrg-ben

# execute like this: ./2024_bcftools_extract_snps.sh input
# load these modules before running:
module load StdEnv/2023  gcc/12.3 bcftools/1.19
module load StdEnv/2023 gatk/4.4.0.0 java/21.0.1
module load tabix

bcftools view --include 'TYPE="snp"' ${1} > ${1}_snpsonly.vcf
bgzip ${1}_snpsonly.vcf
bcftools index ${1}_snpsonly.vcf.gz
gatk --java-options -Xmx10G IndexFeatureFile -I ${1}_snpsonly.vcf.gz

# I think this option retains indels; it generates a larger file than the above
# bcftools view -e 'TYPE="ref"' ${1} > ${1}_noref.vcf
# bgzip ${1}_noref.vcf
# bcftools index ${1}_noref.vcf.gz
# gatk --java-options -Xmx10G IndexFeatureFile -I ${1}_noref.vcf.gz
```
# Now concatenate these for each chr:
```
#!/bin/sh
#SBATCH --job-name=bcftools
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=2gb
#SBATCH --output=bcftools.%J.out
#SBATCH --error=bcftools.%J.err
#SBATCH --account=rrg-ben

# execute like this: sbatch 2024_bcftools_concat.sh output.filtered.snps.1.*.AB.vep.vcf.gz_snpsonly.vcf.gz
# load these modules before running:
module load StdEnv/2023  gcc/12.3 bcftools/1.19 gatk/4.4.0.0 java/17.0.6 tabix

command="bcftools concat --allow-overlaps --rm-dups all "

for file in ./output.filtered.snps.19.*.AB.vep.vcf.gz_snpsonly.vcf.gz; do         # Use ./* ... NEVER bare *    
    if [ -e "${file}" ] ; then   # Check whether file exists.
	command="${command} ${file} "
    fi
done

command="${command} -O z -o ${1}_concat.vcf.gz"
echo ${command}
${command}

#bgzip -c ${1}_concat.vcf > ${1}_concat.vcf.gz
# make csi and tbi indexes for the concatenated and compressed vcf.gz
bcftools index ${1}_concat.vcf.gz
gatk --java-options -Xmx10G IndexFeatureFile -I ${1}_concat.vcf.gz
```
# Now extract only Sulawesi for PCA:
```
#!/bin/sh
#SBATCH --job-name=bcftools_sula
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=2gb
#SBATCH --output=bcftools_sula.%J.out
#SBATCH --error=bcftools_sula.%J.err
#SBATCH --account=rrg-ben

# execute like this: sbatch 2024_bcftools_extractSulaSamples.sh output.filtered.snps.1.all.AB.vep.vcf.gz_snpsonly.vcf.gz_concat.vcf.gz
# load these modules before running:
module load StdEnv/2023  gcc/12.3 bcftools/1.19 gatk/4.4.0.0 java/17.0.6 tabix

bcftools view -S Sula_Samples.txt ${1} > ${1}_Sulasnpsonly.vcf
bgzip ${1}_Sulasnpsonly.vcf
bcftools index ${1}_Sulasnpsonly.vcf.gz
gatk --java-options -Xmx10G IndexFeatureFile -I ${1}_Sulasnpsonly.vcf.gz
```
