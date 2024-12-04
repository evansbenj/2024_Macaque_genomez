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
        - s103291
        - s103275
        - s103283
        - s103302
        - s105212
        - s105213
        - SAMN07508142
        - SAMN07508143
        - SAMN07508144
        - SAMN07508145
        - SAMN07508146
        - s103270
        - s103278
        - s103247
        - s103206
        - s103214
        - s103222
        - s103230
        - s103238
        - s103246
        - s103286
        - s103215
        - s103255
        - s103271
        - s103279
        - s103200
TON:             
        - s105223
        - s105222
        - s105225
        - s105226
        - s103318
        - s105206
        - s105178
        - s105183
        - s105184
        - s105185
        - s105199
        - s105200
        - s105201
        - s103303
        - s103319
        - s103327
        - s103335
        - s103312
        - s103320
        - s103376
        - SAMN07503429
        - SAMN18571130
        - SAMN18571131
        - SAMN18571132
        - SAMN18571133
        - s103235
        - s103243
        - s103251
        - s103267
HEC:                 
        - s103346
        - s105180
        - s103306
        - s103314
        - SAMN07508136
        - SAMN07508137
        - SAMN07508138
        - SAMN18570967
        - SAMN18570968
TOG:                 
        - s105187
        - s103334
        - s105188
        - s105190
        - s105191
        - s105194
        - s105198
        - s103326
        - s105186
        - s105192
        - s105193
        - s105195
        - s105196
        - s103374
        - s105197
        - SAMN07508162
OCH:                 
        - s105202
        - s105203
        - s105205
        - s105207
        - s105204
        - s105208
        - s105209
MTH:                 
        - s103294
        - s103260
        - s105211
        - s105215
        - s103212
        - s103228
        - s103264
        - s103280
        - s103288
        - s103272
        - s103209
        - s103233
        - s103202
        - s103258
TMH:                 
        - s103204
        - s105210
        - s103252
        - s103313
        - s105214
        - s103224
        - s103274
        - s103282
        - s103219
        - s103226
        - s103266
        - s103273
        - s103281
        - s103289
        - s103227
        - s103249
        - s103257
        - s103265
        - s103210
        - s103234
        - s103203
        - s103211
        - s103259
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
        - s105224
        - s105181
        - s105179
        - SAMN18571119
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
Ok now do the analysis for each sample and for each chromosome:
```
#!/bin/sh
#SBATCH --job-name=admixfrog_analysis
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=6:00:00
#SBATCH --mem=2gb
#SBATCH --output=admixfrog_analysis.%J.out
#SBATCH --error=admixfrog_analysis.%J.err
#SBATCH --account=rrg-ben

# BRU HEC MAU NEM SUM NGA NGE TOG TON

# execute like this:
# sbatch admixfrog_do_analysis.sh 

module load StdEnv/2023 python/3.12.4
# run the analyses for each chr
/home/ben/.local/bin/admixfrog --infile ${1}_Chr1.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr1_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr2.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr2_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr3.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr3_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr4.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr4_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr5.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr5_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr6.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr6_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr7.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr7_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr8.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr8_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr9.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr9_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr10.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr10_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr11.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr11_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr12.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr12_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr13.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr13_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr14.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr14_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr15.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr15_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr16.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr16_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr17.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr17_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr18.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr18_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr19.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr19_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_Chr20.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_Chr20_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_ChrX.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_ChrX_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination

/home/ben/.local/bin/admixfrog --infile ${1}_ChrY.in.xz --ref Sula_only_MAU_TON_HEC_Chr1.ref.xz --out ${1}_ChrY_MAU_TON_HEC.out -b 10000 --states MAU TON HEC --c0 0 --dont-est-contamination
```
