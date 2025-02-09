# General Genomics

After hardfiltering, I removed positions with missing genotypes:
```
vcftools --gzvcf all_162_maqs_chr1.vcf.gz --max-missing-count 0 --min-alleles 2 --max-alleles 2 --minQ 30 --recode --recode-INFO-all --stdout | gzip -c > all_162_maqs_chr1_maxmissingcount_0_biallelic_genoqual30.vcf.gz
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

Now calculate Fst in windows:
```
python3 /home/ben/2025_genomics_general/genomics_general/popgenWindows.py -w 5000000 -m 100 -g all_162_maqs_chr14_maxmissingcount_0_biallelic_genoqual30_b.geno.gz -o all_162_maqs_chr14_maxmissingcount_0_biallelic_genoqual30_b.geno_diversity_SUM_BOR_PAG.csv.gz -f phased -T 5 -p SUM -p BOR -p PAG --popsFile pops_all.txt
```

The gz output for each chr can be downloaded. Then concated them using this script:
```
#!/usr/bin/env perl
use strict;
use warnings;


# This program globs a bunch of chromosome output files from the script called Performs_ABBA_BABA_on_populations.pl
# and concatenates them for analysis by "Does_block_jackknife.pl"

# Run like this "Makes_inputfile_for_jackknife.pl inputprefix"

# The output prefix will be: inputprefix.concat

my $inputfile = $ARGV[0];


my @files = glob("'*${inputfile}*'");

print "hi @files\n";

# open an outputfile
unless (open(OUTFILE, ">$inputfile.concat"))  {
	print "I can\'t write to $inputfile.concat\n";
	exit;
}
print "Creating output file: $inputfile.concat\n";


unless (open DATAINPUT, "gunzip -c $files[0] |") {
	print "Can not find the input file.\n";
	exit;
}

while ( my $line = <DATAINPUT>) {
	print OUTFILE $line;
}		
close DATAINPUT;

my $counter=0;

foreach my $infile (1..$#files){
	unless (open DATAINPUT, "gunzip -c $files[$infile] |") {
		print "Can not find the input file.\n";
		exit;
	}
	$counter=0;
	while ( my $line = <DATAINPUT>) {
		if($counter>0){
			print OUTFILE $line;
		}
		$counter+=1;
	}		
	close DATAINPUT;
}	

close OUTFILE
```

It will be interesting to compare 
* Fst of autosomes and the X (expect higher for the X)
* Fst of the X of females and the X of males (expect higher for female X)
* Fst of the X of males and the Y of males (expect higher for X)
