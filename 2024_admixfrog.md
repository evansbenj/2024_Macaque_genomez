# Admixfrog

Directory:
```
/home/ben/projects/rrg-ben/ben/2024_macaques/concatenated_vcfs/SulaSNPs_only_hardfiltered_and_thinned
```

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
Then I filtered to include only sites with the number of missing genotypes being less than or equal to 0 and setting the minimum genotype quality at 30 as follows:
```
#!/bin/sh
#SBATCH --job-name=vcftools_filter
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --mem=2gb
#SBATCH --output=vcftools_filter.%J.out
#SBATCH --error=vcftools_filter.%J.err
#SBATCH --account=rrg-ben

# execute like this: sbatch 2024_bcftools_concat.sh output.filtered.snps.1.*.AB.vep.vcf.gz_snpsonly.vcf.gz
# load these modules before running:
module load StdEnv/2023  gcc/12.3 bcftools/1.19 gatk/4.4.0.0 java/17.0.6 tabix vcftools
#vcftools --gzvcf ${1} --max-missing 1 --recode --recode-INFO-all --out ${1}_nomissing.vcf
vcftools --gzvcf ${1} --max-missing-count 0 --minQ 30 --recode --recode-INFO-all --out SulaSNPs.${2}_maxmissingcount_0_genoqual30.vcf

#bgzip -c SulaSNPs.${2}_maxmissingcount_0_genoqual30.vcf
bgzip -c SulaSNPs.${2}_maxmissingcount_0_genoqual30.vcf > SulaSNPs.${2}_maxmissingcount_0_genoqual30.vcf.gz
bcftools index SulaSNPs.${2}_maxmissingcount_0_genoqual30.vcf.gz
gatk --java-options -Xmx10G IndexFeatureFile -I SulaSNPs.${2}_maxmissingcount_0_genoqual30.vcf.gz
```
Then I thinned, saving only positions that are at least 500 bp apart.
```
#!/bin/sh
#SBATCH --job-name=vcftools_thin
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --mem=2gb
#SBATCH --output=vcftools_thin.%J.out
#SBATCH --error=vcftools_thin.%J.err
#SBATCH --account=rrg-ben

# execute like this: sbatch 2024_bcftools_concat.sh output.filtered.snps.1.*.AB.vep.vcf.gz_snpsonly.vcf.gz
# load these modules before running:
module load StdEnv/2023  gcc/12.3 bcftools/1.19 gatk/4.4.0.0 java/17.0.6 tabix vcftools
#vcftools --gzvcf ${1} --max-missing 1 --recode --recode-INFO-all --out ${1}_nomissing.vcf
vcftools --gzvcf ${1} --max-missing-count 0 --minQ 30 --recode --recode-INFO-all --out SulaSNPs.${2}_maxmissingcount_0_genoqual30.vcf
vcftools --gzvcf ${1} --thin 500 --recode --recode-INFO-all --out ${1}_thinned
bgzip -c ${1}_thinned.recode.vcf > ${1}_thinned.recode.vcf.gz
tabix -p vcf ${1}_thinned.recode.vcf.gz

#bgzip -c SulaSNPs.${2}_maxmissingcount_0_genoqual30.vcf
#bgzip -c SulaSNPs.${2}_maxmissingcount_0_genoqual30.vcf > SulaSNPs.${2}_maxmissingcount_0_genoqual30.vcf.gz
#bcftools index SulaSNPs.${2}_maxmissingcount_0_genoqual30.vcf.gz
#gatk --java-options -Xmx10G IndexFeatureFile -I SulaSNPs.${2}_maxmissingcount_0_genoqual30.vcf.gz
```

# Make reference

I had to change the sample names to include an alphabetic componenet because admixfrog doesn't handle sample names that are numeric:

First turn off tab complete in terminal
```
bind '"\t":self-insert'
```
now search replace each vcf file
```
sed -i 's/103200	103202	103203	103204	103206	103209	103210	103211	103212	103214	103215	103219	103222	103224	103226	103227	103228	103230	103233	103234	103235	103238	103243	103246	103247	103249	103251	103252	103255	103257	103258	103259	103260	103264	103265	103266	103267	103270	103271	103272	103273	103274	103275	103278	103279	103280	103281	103282	103283	103286	103288	103289	103291	103294	103302	103303	103306	103312	103313	103314	103315	103318	103319	103320	103326	103327	103334	103335	103346	103371	103374	103376	103391	105178	105179	105180	105181	105183	105184	105185	105186	105187	105188	105190	105191	105192	105193	105194	105195	105196	105197	105198	105199	105200	105201	105202	105203	105204	105205	105206	105207	105208	105209	105210	105211	105212	105213	105214	105215	105222	105223	105224	105225	105226	SAMN07503429	SAMN07503430	SAMN07508135	SAMN07508136	SAMN07508138	SAMN07508142	SAMN07508143	SAMN07508144	SAMN07508145	SAMN07508146	SAMN07508155	SAMN07508156	SAMN07508157	SAMN07508162	SAMN18570966	SAMN18570967	SAMN18570968	SAMN18571119	SAMN18571130	SAMN18571131	SAMN18571132	SAMN18571133/s103200	s103202	s103203	s103204	s103206	s103209	s103210	s103211	s103212	s103214	s103215	s103219	s103222	s103224	s103226	s103227	s103228	s103230	s103233	s103234	s103235	s103238	s103243	s103246	s103247	s103249	s103251	s103252	s103255	s103257	s103258	s103259	s103260	s103264	s103265	s103266	s103267	s103270	s103271	s103272	s103273	s103274	s103275	s103278	s103279	s103280	s103281	s103282	s103283	s103286	s103288	s103289	s103291	s103294	s103302	s103303	s103306	s103312	s103313	s103314	s103315	s103318	s103319	s103320	s103326	s103327	s103334	s103335	s103346	s103371	s103374	s103376	s103391	s105178	s105179	s105180	s105181	s105183	s105184	s105185	s105186	s105187	s105188	s105190	s105191	s105192	s105193	s105194	s105195	s105196	s105197	s105198	s105199	s105200	s105201	s105202	s105203	s105204	s105205	s105206	s105207	s105208	s105209	s105210	s105211	s105212	s105213	s105214	s105215	s105222	s105223	s105224	s105225	s105226	SAMN07503429	SAMN07503430	SAMN07508135	SAMN07508136	SAMN07508138	SAMN07508142	SAMN07508143	SAMN07508144	SAMN07508145	SAMN07508146	SAMN07508155	SAMN07508156	SAMN07508157	SAMN07508162	SAMN18570966	SAMN18570967	SAMN18570968	SAMN18571119	SAMN18571130	SAMN18571131	SAMN18571132	SAMN18571133/g' SulaSNPs.Chr18_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf
```
now turn tab complete on again
```
bind '"\C-i":complete'
```
Now zip and index each one:
```
#!/bin/sh
#SBATCH --job-name=bgzip
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --mem=2gb
#SBATCH --output=bgzip.%J.out
#SBATCH --error=bgzip.%J.err
#SBATCH --account=rrg-ben

# execute like this: sbatch 2024_bgzip_and_tabix.sh vcffile.vcf
# load these modules before running:
module load StdEnv/2023  gcc/12.3 bcftools/1.19 gatk/4.4.0.0 java/17.0.6 tabix vcftools

bgzip -c ${1} > ${1}.gz
tabix -p vcf ${1}.gz
```
OK now make the refs
```
#!/bin/sh
#SBATCH --job-name=AF_ref
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=2gb
#SBATCH --output=AF_ref.%J.out
#SBATCH --error=AF_ref.%J.err
#SBATCH --account=rrg-ben

# BRU HEC MAU NEM SUM NGA NGE TOG TON

# execute like this:
# sbatch 2024_admixfrog_make_ref.sh vcf chr
module load StdEnv/2023 python/3.12.4
# run the analyses for each chr
/home/ben/.local/bin/admixfrog-ref --out Sula_only_MAU_TON_HEC_${2}.ref.xz --vcf-ref ${1} --state-file spops.yaml --states MAU TON HEC --chroms 1-20,X,Y
```
where the file `spops.yaml` is this (with only spaces, no tabs):
```
MAU:                 
        - s103200
        - s103206
        - s103214
        - s103215
        - s103222
        - s103230
        - s103238
        - s103246
        - s103247
        - s103255
        - s103270
        - s103271
        - s103275
        - s103278
        - s103279
        - s103283
        - s103286
        - s103291
        - s103302
        - s105212
        - s105213
        - SAMN07508142
        - SAMN07508143
        - SAMN07508144
        - SAMN07508145
        - SAMN07508146
TON:             
        - s103235
        - s103243
        - s103251
        - s103267
        - s103303
        - s103312
        - s103318
        - s103319
        - s103320
        - s103327
        - s103335
        - s105178
        - s105184
        - s105201
        - s105206
        - s105214
        - s105222
        - SAMN07503429
        - SAMN18571130
        - SAMN18571131
        - SAMN18571132
HEC:                 
        - s103306
        - s103314
        - s103346
        - s105179
        - s105180
        - SAMN07508136
        - SAMN07508138
        - SAMN18570967
        - SAMN18570968
TOG:                 
        - s103326
        - s103334
        - s103374
        - s105186
        - s105187
        - s105188
        - s105190
        - s105191
        - s105192
        - s105193
        - s105194
        - s105195
        - s105196
        - s105197
        - s105198
        - SAMN07508162
THH:                 
        - s105225
GTH:
        - s105199
OCH:                 
        - s103376
        - s105202
        - s105207
        - s105209
MTH:                 
        - s103202
        - s103209
        - s103212
        - s103228
        - s103233
        - s103258
        - s103264
        - s103272
        - s103280
        - s103288
        - s103294
        - s105211
        - s105215
TMH:                 
        - s103203
        - s103204
        - s103210
        - s103211
        - s103219
        - s103224
        - s103226
        - s103227
        - s103234
        - s103249
        - s103252
        - s103257
        - s103259
        - s103260
        - s103265
        - s103266
        - s103273
        - s103274
        - s103281
        - s103282
        - s103289
        - s103313
        - s105210
        - s105223
        - SAMN18571133
NGE:                 
        - s103315
        - SAMN07508157
        - SAMN18570966
NGA:                 
        - s103391
        - s103371
        - SAMN07503430
        - SAMN07508155
        - SAMN07508156
BRU:                 
        - SAMN07508135     
HTH:                 
        - s105181
        - s105224
        - s105204
        - SAMN18571119
OTH:                 
        - s105203
        - s105205
        - s105208
TGH:
        - s105183
        - s105185
        - s105200
        - s105226
```

OK now make targets for each sample and for each chromosome:
```
#!/bin/sh
#SBATCH --job-name=admixfrog_target
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=2gb
#SBATCH --output=admixfrog_target.%J.out
#SBATCH --error=admixfrog_target.%J.err
#SBATCH --account=rrg-ben

# BRU HEC MAU NEM SUM NGA NGE TOG TON

# execute like this:
# sbatch admixfrog_make_target_input.sh 

module load StdEnv/2023 python/3.12.4
# run the analyses for each chr
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr1_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr1.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr2_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr2.ref.xz --out ${1}_Chr2.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr3_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr3.ref.xz --out ${1}_Chr3.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr4_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr4.ref.xz --out ${1}_Chr4.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr5_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr5.ref.xz --out ${1}_Chr5.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr6_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr6.ref.xz --out ${1}_Chr6.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr7_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr7.ref.xz --out ${1}_Chr7.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr8_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr8.ref.xz --out ${1}_Chr8.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr9_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr9.ref.xz --out ${1}_Chr9.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr10_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr10.ref.xz --out ${1}_Chr10.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr11_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr11.ref.xz --out ${1}_Chr11.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr12_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr12.ref.xz --out ${1}_Chr12.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr13_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr13.ref.xz --out ${1}_Chr13.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr14_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr14.ref.xz --out ${1}_Chr14.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr15_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr15.ref.xz --out ${1}_Chr15.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr16_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr16.ref.xz --out ${1}_Chr16.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr17_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr17.ref.xz --out ${1}_Chr17.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr18_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr18.ref.xz --out ${1}_Chr18.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr19_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr19.ref.xz --out ${1}_Chr19.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.Chr20_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_Chr20.ref.xz --out ${1}_Chr20.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.ChrX_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_ChrX.ref.xz --out ${1}_ChrX.in.xz --chroms 1-20,X,Y
/home/ben/.local/bin/admixfrog-bam --vcfgt SulaSNPs.ChrY_maxmissingcount_0_genoqual30.vcf.recode.vcf.gz_thinned.recode.vcf.gz --target ${1} --ref Sula_only_MAU_TON_HEC_ChrY.ref.xz --out ${1}_ChrY.in.xz --chroms 1-20,X,Y
```
Ok now do the analysis for each sample and for each chromosome; note that I had to use "--female" for the sex chromosomes to force them to be diploid for everyone because the "--male" option is broken (according to the help menu) and it would not work for the X without this.
```
#!/bin/sh
#SBATCH --job-name=admixfrog_analysis
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=6:00:00
#SBATCH --mem=2gb
#SBATCH --output=admixfrog_analysis.%J.out
#SBATCH --error=admixfrog_analysis.%J.err
#SBATCH --account=def-ben

# BRU HEC MAU NEM SUM NGA NGE TOG TON

# execute like this:
# sbatch admixfrog_do_analysis.sh samplename

module load StdEnv/2023 python/3.12.4
# run the analyses for each chr
/home/ben/.local/bin/admixfrog --infile ${1}_Chr1.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr1_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr2.in.xz --ref Sula_only_MAU_TON_HEC_Chr2.ref.xz --out ${1}_Chr2_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr3.in.xz --ref Sula_only_MAU_TON_HEC_Chr3.ref.xz --out ${1}_Chr3_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr4.in.xz --ref Sula_only_MAU_TON_HEC_Chr4.ref.xz --out ${1}_Chr4_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr5.in.xz --ref Sula_only_MAU_TON_HEC_Chr5.ref.xz --out ${1}_Chr5_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr6.in.xz --ref Sula_only_MAU_TON_HEC_Chr6.ref.xz --out ${1}_Chr6_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr7.in.xz --ref Sula_only_MAU_TON_HEC_Chr7.ref.xz --out ${1}_Chr7_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr8.in.xz --ref Sula_only_MAU_TON_HEC_Chr8.ref.xz --out ${1}_Chr8_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr9.in.xz --ref Sula_only_MAU_TON_HEC_Chr9.ref.xz --out ${1}_Chr9_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr10.in.xz --ref Sula_only_MAU_TON_HEC_Chr10.ref.xz --out ${1}_Chr10_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr11.in.xz --ref Sula_only_MAU_TON_HEC_Chr11.ref.xz --out ${1}_Chr11_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr12.in.xz --ref Sula_only_MAU_TON_HEC_Chr12.ref.xz --out ${1}_Chr12_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr13.in.xz --ref Sula_only_MAU_TON_HEC_Chr13.ref.xz --out ${1}_Chr13_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr14.in.xz --ref Sula_only_MAU_TON_HEC_Chr14.ref.xz --out ${1}_Chr14_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr15.in.xz --ref Sula_only_MAU_TON_HEC_Chr15.ref.xz --out ${1}_Chr15_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr16.in.xz --ref Sula_only_MAU_TON_HEC_Chr16.ref.xz --out ${1}_Chr16_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr17.in.xz --ref Sula_only_MAU_TON_HEC_Chr17.ref.xz --out ${1}_Chr17_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr18.in.xz --ref Sula_only_MAU_TON_HEC_Chr18.ref.xz --out ${1}_Chr18_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr19.in.xz --ref Sula_only_MAU_TON_HEC_Chr19.ref.xz --out ${1}_Chr19_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr20.in.xz --ref Sula_only_MAU_TON_HEC_Chr20.ref.xz --out ${1}_Chr20_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_ChrX.in.xz --ref Sula_only_MAU_TON_HEC_ChrX.ref.xz --out ${1}_ChrX_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination --female

/home/ben/.local/bin/admixfrog --infile ${1}_ChrY.in.xz --ref Sula_only_MAU_TON_HEC_ChrY.ref.xz --out ${1}_ChrY_MAU_TON_HEC.oout -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination --female
```

# Square plot

```
## Working directory
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2023_macaque_genomes/2024_macaques/admixfrog")
library(tidyverse)
library(ggplot2)
library(circlize)
library(dplyr)
library(stringr)
library(data.table)
# https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html
# https://cran.r-project.org/web/packages/circlize/circlize.pdf
# initilize
rm(list=ls()) # removes all variables


#MAU
#sample_vector <- c("s103291","s103275","s103283","s103302","s105212","s105213","SAMN07508142","SAMN07508143","SAMN07508144","SAMN07508145","SAMN07508146","s103270","s103278","s103247","s103206","s103214","s103222","s103230","s103238","s103246","s103286","s103215","s103255","s103271","s103279","s103200")
#pop <- "MAU"

#TON
sample_vector <- c("s105223","s105222","s105225","s105226","s103318","s105206","s105178","s105183","s105184","s105185","s105199","s105200","s105201","s103303","s103319","s103327","s103335","s103312","s103320","SAMN07503429","SAMN18571130","SAMN18571131","SAMN18571132","SAMN18571133","s103235","s103243","s103251","s103267")
pop <- "TON"

#HEC
# sample_vector <- c("s103306","s103314","SAMN07508136","s103346","s105180","SAMN07508138","SAMN18570967","SAMN18570968")
# pop <- "HEC"

#MTH
# sample_vector <- c("s103294","s103260","s105211","s105215","s103212","s103228","s103264","s103280","s103288","s103272","s103209","s103233","s103202","s103258")
# pop <- "MTH"

#TMH
# sample_vector <- c("s103204","s105210","s103252","s103313","s105214","s103224","s103274","s103282","s103219","s103226","s103266","s103273","s103281","s103289","s103227","s103249","s103257","s103265","s103210","s103234","s103203","s103211","s103259")
# pop <- "TMH"

#HTH
# sample_vector <- c("s105224","s105181","s105179","SAMN18571119")
# pop <- "HTH"

#TOG
# sample_vector <- c("s105187","s103334","s105188","s105190","s105191","s105194","s105198","s103326","s105186","s105192","s105193","s105195","s105196","s103374","s105197","SAMN07508162")
# pop <- "TOG"

#OCH
# sample_vector <- c("s105202","s105203","s105205","s105207","s105204","s105208","s105209")
# pop <- "OCH"

#NGE
# sample_vector <- c("s103315","SAMN07508157","SAMN18570966")
# pop <- "NGE"

#NGA
# sample_vector <- c("s103391","s103371","SAMN07503430","SAMN07508155","SAMN07508156")
# pop <- "NGA"

# BRU
# sample_vector <- c("SAMN07508135")
# pop <- "BRU"

analysis <-"_MAU_TON_HEC"
chrs <- factor(c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9",
         "Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19","Chr20","ChrX"))


# this is a list of chr lengths
chr_length_list = list("1" = 223616942,
                       "2" = 196197964,
                       "3" = 185288947,
                       "4" = 169963040,
                       "5" = 187317192,
                       "6" = 179085566,
                       "7" = 169868564,
                       "8" = 145679320,
                       "9" = 134124166,
                       "10" = 99517758,
                       "11" = 133066086,
                       "12" = 130043856,
                       "13" = 108737130,
                       "14" = 128056306,
                       "15" = 113283604,
                       "16" = 79627064,
                       "17" = 95433459,
                       "18" = 74474043,
                       "19" = 58315233,
                       "20" = 77137495,
                       "X" = 153388924
                       # Y = 11753682
)

# this is a vector of chr lenths
begins <- rep(1,21)
ends = c(223616942,
         196197964,
         185288947,
         169963040,
         187317192,
         179085566,
         169868564,
         145679320,
         134124166,
         99517758,
         133066086,
         130043856,
         108737130,
         128056306,
         113283604,
         79627064,
         95433459,
         74474043,
         58315233,
         77137495,
         153388924)

# this is needed to initalize the graph below
chr_lengths <- cbind(begins,ends)

# chr names
chr_names <- c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10",
               "Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19",
               "Chr20","ChrX")
chr_names_ordered <- factor(chr_names, ordered = TRUE, 
                            levels = c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10",
                                       "Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19",
                                       "Chr20","ChrX"))

chr_names_ordered_simple <- factor(chr_names, ordered = TRUE, 
                                   levels = c("1","2","3","4","5","6","7","8","9","10",
                                              "11","12","13","14","15","16","17","18","19","20",
                                              "X"))


# loop through each sample
for (sample in sample_vector){
  # sample <- 's105223'
  # loop through chrs
  for(i in levels(chrs)){
    print(eval(i))
    a <- read_csv(paste(eval(sample), "_",i,eval(analysis),".oout.bin.xz", sep=""))
    # a <- read_csv("s105224_Chr1_MAU_TON_HEC.oout.bin.xz")
    assign(i,a)
  }  

  # rename chromosome column
  Chr1$chrom <- "1"
  Chr2$chrom <- "2"
  Chr3$chrom <- "3"
  Chr4$chrom <- "4"
  Chr5$chrom <- "5"
  Chr6$chrom <- "6"
  Chr7$chrom <- "7"
  Chr8$chrom <- "8"
  Chr9$chrom <- "9"
  Chr10$chrom <- "10"
  Chr11$chrom <- "11"
  Chr12$chrom <- "12"
  Chr13$chrom <- "13"
  Chr14$chrom <- "14"
  Chr15$chrom <- "15"
  Chr16$chrom <- "16"
  Chr17$chrom <- "17"
  Chr18$chrom <- "18"
  Chr19$chrom <- "19"
  Chr20$chrom <- "20"
  ChrX$chrom <- "X"
  
  Allchrs <- rbind(Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,
                   Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,
                   Chr20,ChrX)
  
  # make the chr a factor so we can use it for faceting
  Allchrs$chrom <- as.factor(Allchrs$chrom)
  
  # make a dataframe for circular plotting
  Allchr_circular <- as.data.frame(Allchrs[,c(1,3,8,9,10,11,12,13)])
  # create end coordinates based on next start site and have NA for last start site
  Allchr_circular$end <- c(Allchr_circular$pos[-1]-1, NA)
  #View(Allchr_circular)
  # check for changes in chr at the last window
  temp <- ifelse(Allchr_circular$end != Allchr_circular$pos + 99999,
                                Allchr_circular$pos+99999,
                                Allchr_circular$end)
  # add this to the dataframe
  Allchr_circular$end <- temp
  # make an entry for last end positioin 
  Allchr_circular$end[nrow(Allchr_circular)]<-Allchr_circular$pos[nrow(Allchr_circular)]+99999
  #View(Allchr_circular)
  # Now fix the first entry of each chr
  Allchr_circular$end <- ifelse(Allchr_circular$pos < 99999,
                                99999,
                                Allchr_circular$end)
  
  # reorder the columns
  Allchr_circular <- Allchr_circular[, c(1,2,9,3,4,5,6,7,8)]
  names(Allchr_circular)[names(Allchr_circular) == "pos"] <- "start"
  # ok looks good
  Allchr_circular$sample <- sample
  if(sample == sample_vector[1]){
    # make a big df
    big_monkey_df <- data.frame(matrix(ncol = 10, nrow = 0))
    colnames(big_monkey_df) <- colnames(Allchr_circular)
  }
  # now bind this to the big_monkey_df
  big_monkey_df <- rbind(big_monkey_df,Allchr_circular)

long_big_monkey_df <- melt(setDT(big_monkey_df[,-3]), id.vars = c("chrom","start","sample"), variable.name = "ancestry")

long_big_monkey_df$color[long_big_monkey_df$ancestry == "MAU" ] <- "yellow"
long_big_monkey_df$color[long_big_monkey_df$ancestry == "TON" ] <- "red"
long_big_monkey_df$color[long_big_monkey_df$ancestry == "HEC" ] <- "blue"
long_big_monkey_df$color[long_big_monkey_df$ancestry == "MAUTON" ] <- "orange"
long_big_monkey_df$color[long_big_monkey_df$ancestry == "MAUHEC" ] <- "black"
long_big_monkey_df$color[long_big_monkey_df$ancestry == "TONHEC" ] <- "purple"

long_big_monkey_df$color <- factor(long_big_monkey_df$color, ordered = TRUE, 
                                   levels = c("yellow","orange","red","purple","blue","black"))
long_big_monkey_df$chrom <- factor(long_big_monkey_df$chrom, ordered = TRUE, 
                                   levels = c("1","2","3","4","5","6","7","8","9","10",
                                              "11","12","13","14","15","16","17","18","19","20",
                                              "X"))

}

# get rid of rows with zero probability that do not need to be plotted
long_big_monkey_df_smaller <- long_big_monkey_df[long_big_monkey_df$value != 0.0e+00 ]

# temp <- long_big_monkey_df_smaller[long_big_monkey_df_smaller$sample == "s103291" ]

# fiddle with x axis tics
# limits_fun <- function(x) {
#  if (max(x) < 150) {  # This will identify the current "Wind" panel
#    seq(0, 100, by = 100)  # New limits for the "Wind" panel
#  } else if (max(x) < 250) {  # This will identify the current "Temp" panel
#    seq(0, 200, by = 100)  # New limits for the "Temp" panel
#  } 
# }

 

png(paste(eval(pop),eval(analysis),"_stacked.png",sep=""),
    width = 1500, height = eval(length(sample_vector))*35, units='mm', res = 300) 
    # adjust the height depending on the number of samples
   #  ggplot(long_big_monkey_df, aes(x = start, y = value, fill = ancestry, color = color)) + 
  # ggplot(temp %>% arrange(sample,chrom),   
  ggplot(long_big_monkey_df_smaller %>% arrange(sample,chrom), 
            aes(x = start/1000000, y = value, fill = ancestry)) + 
        geom_bar(position='fill', stat='identity') +
    #scale_fill_brewer(type = "seq", palette = 6) +
    scale_fill_manual(name="Ancestry", values = c("MAU"="yellow",
                                                    "MAUTON"="orange",
                                                    "TON"="red",
                                                    "TONHEC"="purple",
                                                    "HEC"="blue",
                                                    "MAUHEC"="black"),
      breaks=c("MAU","MAUTON","TON","TONHEC","HEC","MAUHEC"),
      labels=c("MAU","MAUTON","TON","TONHEC","HEC","MAUHEC"))+
      # facet_wrap(sample ~ chrom, ncol=21, scale="free_x")+
      facet_grid(sample ~ chrom, scales = "free", space="free") +
      # scale_x_continuous(breaks = limits_fun) +
      scale_x_continuous(breaks = c(0,100,200)) +
      scale_y_continuous(breaks = c(0,1)) +
      labs(x = "Chromosome and Coordinates (100Mb)", y = "Probability") +
      theme_classic(base_size = 38) + 
      # guides(color = FALSE) +
      theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
      theme(strip.text.y.right = element_text(angle = 0)) +
      theme(strip.background = element_blank()) +
      guides(fill=guide_legend(title="Ancestry")) +
      ggtitle(eval(pop))
      # +
      # Change horizontal spacing between facets
      # theme(panel.spacing.x = unit(0.5, "lines")) +
      # Change vertical spacing between facets
      # theme(panel.spacing.y = unit(0.5, "lines"))
dev.off()

```
