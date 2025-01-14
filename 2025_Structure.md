# Structure analysis

I'm going to use our vcf files as input to do Structure analysis using AdmixPipe3 (https://github.com/stevemussmann/admixturePipeline). I'd like to summarize multiple runs using CLUMPAK (https://clumpak.tau.ac.il/download.html).

First step is to remove extraneous info from the vcf file:
```
module load StdEnv/2020  gcc/9.3.0 bcftools/1.9
bcftools annotate --remove FORMAT ../SulaSNPs_only_hardfiltered_and_thinned/thinned_vcfz/SulaSNPs.Chr1_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf -Ov -o temp.vcf
```

To install the AdmixPipe3 pipeline (python3-based) on computecanada do this:
```
module load apptainer
apptainer build --sandbox bb.dir docker://mussmann/admixpipe:3.2
```
