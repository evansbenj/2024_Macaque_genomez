# Admixture analysis

I'm going to use our vcf files as input to do Structure analysis using AdmixPipe3 (https://github.com/stevemussmann/admixturePipeline). I'd like to summarize multiple runs using CLUMPAK (https://clumpak.tau.ac.il/download.html).

I'm working with thinned files that have the '0/0' genotypes with no data removed in this directory
```
/home/ben/projects/rrg-ben/ben/2024_macaques/concatenated_vcfs/all_hardfiltered/160_vcfz
```
The data are thinned to include every 5000 bp as detailed for the admix frog analysis

on info I then concatenated the autosomal chrs:
```
bcftools concat all160.chr{1..20}_maxmissingcount_0_genoqual30.vcf.gz_snpsonly.vcf.gz_5000_thinned.recode.vcf.gz -Ov -o all_160_maqs_allautsomal_chrs_maxmissingcount_0_genoqual30_thin_5000.recode.vcf
```
and then I used plink to make the input files:
```
plink --vcf all_160_maqs_allautsomal_chrs_maxmissingcount_0_genoqual30_thin_5000.recode.vcf --make-bed --out all_160_maqs_allautsomal_chrs_maxmissingcount_0_genoqual30_thin_5000 --allow-extra-chr

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just change the first column from 0 to 1
# this did not work:
# awk '{$1="0";print $0}' all_160_maqs_allautsomal_chrs_maxmissingcount_0_genoqual30_thin_5000.bim > all_160_maqs_allautsomal_chrs_maxmissingcount_0_genoqual30_thin_5000.bim.tmp
but this did:
cp all_160_maqs_allautsomal_chrs_maxmissingcount_0_genoqual30_thin_5000.bim tmp
sed -i 's/0 . 0/1 . 1/g' tmp
mv tmp all_160_maqs_allautsomal_chrs_maxmissingcount_0_genoqual30_thin_5000.bim
```
I then opened 20 screens and did independent analyses (with different seeds) in separate folders:
```
for i in {2..12}
do
/home/ben/2024_macaques/admixture/releases/admixture_linux-1.3.0/admixture --cv --seed $((1 + $RANDOM % 1000)) ../all_160_maqs_allautsomal_chrs_maxmissingcount_0_genoqual30_thin_5000.bed $i > log${i}.out
done
```
I did this also for the X chromosome.

# Rename folders
```
#!/usr/bin/env perl
use strict;
use warnings;

# This script will rename directories for admixXtureplotter

# Execute it in a directory that has 2 folders in it called
# admixX1...admixX20


my $y= 12; # this is the number of ancestral populations
my $x;
my $z;
my $command;


#my @directorynamez = ("admixX1/","admixX2/","admixX3/","admixX4/","admixX5/","admixX6/",
#						"admixX7/","admixX8/","admixX9/","admixX10/",
#						"admixX11/","admixX12/","admixX13/","admixX14/","admixX15/","admixX16/",
#						"admixX17/","admixX18/","admixX19/","admixX20/",);

# make directories
for ($x = 2 ; $x <= $y ; $x++ ) {
	$command = "mkdir ".$x;
	system($command);
	$command = "mkdir ".$x."/Logs";
	system($command);
	for ($z = 1 ; $z <= 20 ; $z++ ) {
		$command = "mkdir ".$x."/".$z;
		system($command);
	}	
}


for ($z = 1 ; $z <= 20 ; $z++ ) {
	for ($x = 2 ; $x <= $y ; $x++ ) {
		$command = "mv admixX".$z."/*.".$x.".P ".$x."/".$z."/.";
		#print $command,"\n";
		system($command);
		$command = "mv admixX".$z."/*.".$x.".Q ".$x."/".$z."/.";
		#print $command,"\n";
		system($command);
		$command = "mv admixX".$z."/log".$x.".out ".$x."/Logs/".$x."_".$z.".log";
		#print $command,"\n";
		system($command);
	}	
}		

```

# AdmixturePlotter

I used scripts in the AdmixturePlotter pipeline (https://github.com/TCLamnidis/AdmixturePlotter) to plot the admixture results.

To compile the results after the files, filenames, and directories have been set up as detailed in the github page, you can run the bash script: 
```
./CompileData_Automes.sh 
```

Then to plot, run the R script `AdmixturePlotter.R`:
```
Rscript AdmixturePlotter.R -i ../Plotting/compound.labelled.QperK.txt -c ../colorlist.txt -p ../poporder.txt
```



# Below is for Admix pipeline

Now combine the autosomal chromosomes and then remove extraneous info from the vcf file:

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
