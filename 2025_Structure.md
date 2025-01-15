# Structure analysis

I'm going to use our vcf files as input to do Structure analysis using AdmixPipe3 (https://github.com/stevemussmann/admixturePipeline). I'd like to summarize multiple runs using CLUMPAK (https://clumpak.tau.ac.il/download.html).

First step is to combine the autosomal chromosomes and then remove extraneous info from the vcf file:

```
module load StdEnv/2020  gcc/9.3.0 bcftools/1.9
bcftools concat ../SulaSNPs_only_hardfiltered_and_thinned/thinned_vcfz/SulaSNPs.Chr{1..20}_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf -Ov -o  ../SulaSNPs_only_hardfiltered_and_thinned/thinned_vcfz/SulaSNPs.AllChr_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf
bcftools annotate --remove FORMAT ../SulaSNPs_only_hardfiltered_and_thinned/thinned_vcfz/SulaSNPs.Chr1_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf -Ov -o temp.vcf
```

To install the AdmixPipe3 pipeline (python3-based) on computecanada do this:
```
module load apptainer
apptainer build --sandbox bb.dir docker://mussmann/admixpipe:3.2
```

To run the analysis separately for each k value try this:
```
#!/bin/sh
#SBATCH --job-name=doAdmixPipe3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=120:00:00
#SBATCH --mem=12gb
#SBATCH --output=doAdmixPipe3.%J.out
#SBATCH --error=doAdmixPipe3.%J.err
#SBATCH --account=def-ben

module load StdEnv/2020 scipy-stack/2022a python/3.10.2 vcftools/0.1.16 plink/1.9b_6.21-x86_64 plink/1.9b_6.21-x86_64 admixture/1.3.0

#app/scripts/python/admixturePipeline/admixturePipeline.py -m popmap_Sulaonly.txt -v ../SulaSNPs_only_hardfiltered_and_thinned/thinned_vcfz/SulaSNPs.AllChr_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode_simplifiedsimplified.vcf --rep 2 -k 1 -K 2 -n 16 -S 1 -a 0.05 # -S 1 excludes all sites with missing data, -a 0.05 requires minimum alllele freq of 0.05

app/scripts/python/admixturePipeline/admixturePipeline.py -m popmap_Sulaonly.txt -v tempp.vcf --rep 10 -k ${1} -K ${1} -n 1 -S 1 -a 0.05 # you can also do parallel jobs by adding a flad -n 8
```

In addition to plotting a bar chart, maybe do pie charts on a map:
```
https://github.com/Tom-Jenkins/mapmixture
```
