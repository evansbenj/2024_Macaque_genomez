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

Now calculate Fst in 30Kb windows (on info2020 only):
```
python3 /home/ben/2025_genomics_general/genomics_general/popgenWindows.py -w 30000 -m 50 -g all_162_maqs_chr14_maxmissingcount_0_biallelic_genoqual30_b.geno.gz -o all_162_maqs_chr14_maxmissingcount_0_biallelic_genoqual30_b.geno_diversity_SUM_BOR_PAG.csv.gz -f phased -T 5 -p SUM -p BOR -p PAG --popsFile pops_all.txt --writeFailedWindows
```
30kb windows ended up having many windows with less than 50 variable positions. This may be due to the agressive filtering, which removed positions with any missing genotype. Instead let's try 100kb windows and also relax the minmum number of sites to be 30 instead of 50 (on info2020 only):
```
python3 /home/ben/2025_genomics_general/genomics_general/popgenWindows.py -w 100000 -m 30 -g all_162_maqs_chr14_maxmissingcount_0_biallelic_genoqual30_b.geno.gz -o all_162_maqs_chr14_maxmissingcount_0_biallelic_genoqual30_b.geno_diversity_SUM_BOR_PAG.csv.gz -f phased -T 5 -p SUM -p BOR -p PAG --popsFile pops_all.txt --writeFailedWindows

```
where the pops_all.txt file is this:
```
SAMN07503410	BOR
SAMN07508148	BOR
SAMN07508153	BOR
SAMN07508154	BOR
SAMN18570964	BOR
SAMN07508135	BRU
105199		GTH
103306		HEC
103314		HEC
103346		HEC
105180		HEC
SAMN07508136	HEC
SAMN07508138	HEC
SAMN18570967	HEC
SAMN18570968	HEC
105179		HTH
105181		HTH
105204		HTH
105224		HTH
SAMN18571119	HTH
103200		MAU
103206		MAU
103214		MAU
103215		MAU
103222		MAU
103228		MAU
103230		MAU
103238		MAU
103246		MAU
103247		MAU
103255		MAU
103270		MAU
103271		MAU
103275		MAU
103278		MAU
103279		MAU
103283		MAU
103286		MAU
103291		MAU
103302		MAU
105212		MAU
105213		MAU
SAMN07508142	MAU
SAMN07508143	MAU
SAMN07508144	MAU
SAMN07508145	MAU
SAMN07508146	MAU
103202		MTH
103209		MTH
103212		MTH
103233		MTH
103258		MTH
103264		MTH
103272		MTH
103280		MTH
103288		MTH
103294		MTH
105211		MTH
105215		MTH
103371		NGA
103391		NGA
SAMN07503430	NGA
SAMN07508155	NGA
SAMN07508156	NGA
103315		NGS
SAMN07508157	NGS
SAMN18570966	NGS
103376		OCH
105202		OCH
105207		OCH
105209		OCH
105203		OTH
105205		OTH
105208		OTH
103404		PAG
103412		PAG
103467		PAG
105227		SIN
105230		SIN
105169		SUM
105170		SUM
105171		SUM
105172		SUM
105173		SUM
105174		SUM
105175		SUM
105176		SUM
105177		SUM
105216		SUM
105217		SUM
105218		SUM
105219		SUM
105220		SUM
105221		SUM
SAMN07508151	SUM
105225		THH
103203		TMH
103204		TMH
103210		TMH
103211		TMH
103219		TMH
103224		TMH
103226		TMH
103227		TMH
103234		TMH
103249		TMH
103252		TMH
103257		TMH
103259		TMH
103260		TMH
103265		TMH
103266		TMH
103273		TMH
103274		TMH
103281		TMH
103282		TMH
103289		TMH
103313		TMH
105210		TMH
105223		TMH
SAMN18571133	TMH
103326		TOG
103334		TOG
103374		TOG
105186		TOG
105187		TOG
105188		TOG
105190		TOG
105191		TOG
105192		TOG
105193		TOG
105194		TOG
105195		TOG
105196		TOG
105197		TOG
105198		TOG
SAMN07508162	TOG
103235		TON
103243		TON
103251		TON
103267		TON
103303		TON
103312		TON
103318		TON
103319		TON
103320		TON
103327		TON
103335		TON
105178		TON
105183		TON
105184		TON
105185		TON
105200		TON
105201		TON
105206		TON
105214		TON
105222		TON
105226		TON
SAMN07503429	TON
SAMN18571130	TON
SAMN18571131	TON
SAMN18571132	TON
```

The gz output for each chr can be downloaded. Then concated them using this script:
```
#!/usr/bin/env perl
use strict;
use warnings;


# This program globs a bunch of chromosome output files 
# and concatenates them

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

# Making density plots and doing permutations
This perl script will do permutations to test whether the difference between the mean Fst of Ninteract and nonNinteract windows is likely to occur by chance. It also makes and input file for plotting density plots of Fst using a R script that follows.

```perl
#!/usr/bin/env perl
use strict;
use warnings;


# This program will read in two files. The first contains the coordinates of
# all N-mt interact genes, their acronyms, and whether or not (1 or 0) they interact
# directly with mt genes.

# The other file is a file with Fst (or pi) in windows with coordinates.  First
# the mean Fst of N-mt interacting (1) and non-interacting (0) genes will be calculated
# then permutations will be performed where the difference between these categories is 
# recalculated after the interaction is permuted n times. This will allow a p value of the
# Fst value to be estimated.

# this program deliverately ignores all genes in chrX.  A separate script will be generated
# that is only for chrX

# execute like this:
# ./2025_All_N_mt_allinteract_fst_permutation.pl rheMac10.ncbiRefSeq.CDS.marked.corrected.orientation.oout diversity.concat 11

my $inputfile1 = $ARGV[0];
my $inputfile2 = $ARGV[1];
my $column = $ARGV[2]; # column with Fst of interest (in 0 based numbering so 13th column is a 12)

my @windowsites;
my @Fst_values;
my $sumsites=0;
my $counter=0;
my @temp;
my $y;
my $x;
my %OXPHOS;
my %windows;

# first open up the OXPHOS gene info (# rheMac10.ncbiRefSeq.CDS.marked.corrected.new.oout)
unless (open DATAINPUT, $inputfile1) {
	print "Can not find the input file.\n";
	exit;
}

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if(($temp[0] ne 'gene')&&($temp[2] ne 'chrX')){ # deliberately ignores chrX
		# deliberately ignoring orientation because we are considering any overlap
		$OXPHOS{$temp[2]."_".$temp[3]."_".$temp[4]}{"gene"} = $temp[0];
		$OXPHOS{$temp[2]."_".$temp[3]."_".$temp[4]}{"complex"} = $temp[1];
		$OXPHOS{$temp[2]."_".$temp[3]."_".$temp[4]}{"mt_interact"} = $temp[5];
	}	
}		
close DATAINPUT;

# this will print out a file that I can use for plotting later

my $outputfile = $inputfile2."_FST__density.txt"; # the name of the output file is from the commandline
unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}

# now open up the Fst data
unless (open DATAINPUT2, $inputfile2) {
	print "Can not find the input file.\n";
	exit;
}

my @temp1;
my $N_interact_window=0;
my $gene_containing_window=0;
my $number_of_genes_in_this_window=0;
my $number_of_gene_beginnings_in_this_window=0;
my $number_of_Ninteract_gene_beginnings_in_this_window=0;
my $number_of_Ninteract_genes_in_this_window=0;
my $Ninteract_acronym="-";
my $number_of_Ninteract_genes_spanning_a_window=0;
my $Fst_no_genez=0; # no genes
my $n_Fst_no_genez=0; # this includes all non associated genes, including those anywhere 

print OUTFILE "chr\tpos\tFst\tcontainsgenes\tcontainsNinteractgenez\tnumber_of_genes\tnumber_of_Ninteractgenez\tNinteract_acronym\n";

while ( my $line = <DATAINPUT2>) {
	chomp($line);
	@temp=split(',',$line); # temp[0] = window_chr, temp[1] = window_start, temp[2] = window_end
	if(($temp[0] ne 'scaffold')&&($temp[$column] ne 'nan')){ # this line is not the beginning and 
															  # also the stat is defined
			$N_interact_window=0;
			$gene_containing_window=0;
			$number_of_genes_in_this_window=0;
			$number_of_gene_beginnings_in_this_window=0;
			$number_of_Ninteract_genes_in_this_window=0;
			$Ninteract_acronym="-";
			# cycle through each gene
			foreach my $key (keys %OXPHOS){
				@temp1=split('_',$key); # temp1[0] = gene_chr, temp1[1] = gene_start, temp1[2] = gene_end
				# check if this window contains one or more N_mt genes
				if(
				($temp1[0] eq $temp[0])&&($temp1[1] >= $temp[1])&&($temp1[1] <= $temp[2])|| 
				#the window includes the the lower coordinates of the gene
				($temp1[0] eq $temp[0])&&($temp1[2] >= $temp[1])&&($temp1[2] <= $temp[2] )|| 
				# the window includes the upper coordinates of the gene 
				# (not necessarily the end depending on orientation)
				($temp1[0] eq $temp[0])&&($temp1[1] <= $temp[1])&&($temp1[2] >= $temp[2])
				# the window has the middle of the gene but not the lower or upper coordinates
				){
						$gene_containing_window=1; 
						$number_of_genes_in_this_window+=1; # will include any portion of a gene
						if(($temp1[0] eq $temp[0])&&($temp1[1] >= $temp[1])&&($temp1[1] <= $temp[2])){
							$number_of_gene_beginnings_in_this_window+=1; # includes only beginning so genes get counted only once across all windows
						}
						# I am not using the OXPHOS hash to record Fst now because instead we are doing everything
						# from the perspective of windows and one gene may span multiple windows and one window
						# may have multiple genes
						#$OXPHOS{$key}{"start_fst"} = $temp[$column];
						#$OXPHOS{$key}{"start_fst_sites"} = $temp[4];
						if($OXPHOS{$key}{"mt_interact"} == 1){
							$N_interact_window=1;
							$number_of_Ninteract_genes_in_this_window+=1; # count any portion of genes in windows
							if(($temp1[0] eq $temp[0])&&($temp1[1] >= $temp[1])&&($temp1[1] <= $temp[2])){
								# includes only beginning so Ninteract genes get counted only once
								$number_of_Ninteract_gene_beginnings_in_this_window+=1;
							}	
							if($number_of_Ninteract_genes_in_this_window == 1){
								$Ninteract_acronym=$OXPHOS{$key}{"gene"};
							}	
							else{
								$Ninteract_acronym=$Ninteract_acronym.",".$OXPHOS{$key}{"gene"}; # just keep track of all the Ninteract genes in all windows (no need to record names for each window)
							}	
						}	
				} # check for any genez
				if($gene_containing_window == 0){
					# this window has no genez
					$Fst_no_genez+=$temp[$column];
					$n_Fst_no_genez+=1;
				}
			}
		}
		if($temp[0] ne 'scaffold'){
			# print information to density file for plotting
			print OUTFILE $temp[0],"\t",$temp[1],"\t",$temp[$column],"\t",$gene_containing_window,"\t",$N_interact_window,
			"\t",$number_of_genes_in_this_window,"\t",$number_of_Ninteract_genes_in_this_window,"\t",$Ninteract_acronym,"\n";
			# record information in a hash for Fst calculations
			$windows{$temp[0]."_".$temp[1]."_".$temp[2]}{"contains_genes"}=$gene_containing_window;
			$windows{$temp[0]."_".$temp[1]."_".$temp[2]}{"contains_mt_interact"}=$N_interact_window;
			$windows{$temp[0]."_".$temp[1]."_".$temp[2]}{"fst"}=$temp[$column];
		}	
}

close OUTFILE;
close DATAINPUT2;

# Now the %windows hash has Fst and coordinates of all windows 
my @fst_for_perms; # this has only the mtinteractors and all genes
my $Fst_associated=0; # all OXPHOS,MRP, ARP2 genes
my $Fst_non_associated=0;
my $n_Fst_associated=0; # this also works for only_N_mt
my $n_Fst_non_associated=0; # this includes all non associated genes, including those anywhere 


my $N_interact_counter=0;

# now calculate the average fst for associated and non-associated OXPHOS genes
foreach my $key (keys %windows){
	if((exists($windows{$key}{"fst"})) &&
		($windows{$key}{"fst"} ne 'nan')){

		if(($windows{$key}{"contains_genes"} != 0)&&($windows{$key}{"contains_mt_interact"} == 1)){
			$N_interact_counter+=1;
			$Fst_associated += $windows{$key}{"fst"};
			$n_Fst_associated += 1;	
			print "hi ",$key," ",$windows{$key}{"fst"},"\n";
		}
		elsif(($windows{$key}{"contains_genes"} != 0)&&($windows{$key}{"contains_mt_interact"} == 0)){
			$Fst_non_associated += $windows{$key}{"fst"};
			$n_Fst_non_associated += 1;
		}
		push(@fst_for_perms,$windows{$key}{"fst"});
	}
}


# now report values
print "Number of windows that have at least a portion of an N_interact gene: ",$N_interact_counter,"\n";
print "Number of windows that have an N_interact beginning: ",$number_of_Ninteract_gene_beginnings_in_this_window,"\n";
print "Mean Fst all N_interact windows with at least a portion of at least one N_interact gene: ",$Fst_associated/$n_Fst_associated,"\n";
print "Mean Fst of windows with only non-N_interact genes: ",$Fst_non_associated/$n_Fst_non_associated,"\n";
print "Mean Fst of windows with no genes: ",$Fst_no_genez/$n_Fst_no_genez,"\n";
#print "Mean Fst non-associated only N_mt ",$Fst_non_associated_only_N_mt/$n_Fst_non_associated_only_N_mt,"\n";


#################
# ALL COMPLEXES perms including all genes
#################

# calculate test statistic for all genes
# this is the difference between the mean Fst for Ninteract windows and the mean for non_Ninteract windows
# permutations will tell us whether this difference is likely to have arisen by chance
my $test_stat = ($Fst_associated/$n_Fst_associated) - ($Fst_non_associated/$n_Fst_non_associated);

# do permutations for all_complex_mt_interact vs all_other_genes
# first make an array that will be shuffled with the same number of 1s and 0s as the $OXPHOS{$key}{"complex"} variable
my @associated_or_not_array = (('1') x $n_Fst_associated, ('0') x $n_Fst_non_associated);

my $perms=2;
my $Fst_associated_perm=0;
my $Fst_not_associated_perm=0;
my $n_Fst_associated_perm=0;
my $n_Fst_not_associated_perm=0;

# first check if the length of the permutation array is the same as the fst array
# for this analysis
if($#associated_or_not_array ne $#fst_for_perms){
	print "Problem: length associated_or_not_array ",$#associated_or_not_array,
	" length fst array ",$#fst_for_perms,"\n";
}

my @perm_diffs;
@perm_diffs={};
for ($y = 0 ; $y < $perms; $y++ ) {
	fisher_yates_shuffle( \@associated_or_not_array );    # permutes @array in place
	$Fst_associated_perm=0;
	$Fst_not_associated_perm=0;
	$n_Fst_associated_perm=0;
	$n_Fst_not_associated_perm=0;
	for ($x = 0 ; $x <= $#associated_or_not_array; $x++ ) {
		if($associated_or_not_array[$x] == 1){
			$Fst_associated_perm+= $fst_for_perms[$x];
			$n_Fst_associated_perm +=1;
		}
		elsif($associated_or_not_array[$x] == 0){	
			$Fst_not_associated_perm+= $fst_for_perms[$x];
			$n_Fst_not_associated_perm +=1;
		}
	}
	push(@perm_diffs,($Fst_associated_perm/$n_Fst_associated_perm) - ($Fst_not_associated_perm/$n_Fst_not_associated_perm));	
}

my @perm_diffs_sorted = sort { $a <=> $b } @perm_diffs;
my $switch=0;
my $pval=0;
# now figure out where the test stat is
for ($y = 0 ; $y <= $#perm_diffs_sorted; $y++ ) {
	if(($test_stat <= $perm_diffs_sorted[$y])&&($switch==0)){
		$pval=$counter;
		$switch = 1;
	}
	$counter+=1;
}	

#print "@perm_diffs_sorted\n";
print "Test stat for test including all genes (if negative, then not significant):",$test_stat,"\n";
print "This is the difference between the mean Fst for Ninteract windows and the mean for non_Ninteract windows\n";
print "P = ",1-($pval/$perms),"\n";




# fisher_yates_shuffle( \@array ) : 
    # generate a random permutation of @array in place
    sub fisher_yates_shuffle {
        my $array = shift;
        my $i;
        for ($i = @$array; --$i; ) {
            my $j = int rand ($i+1);
            next if $i == $j;
            @$array[$i,$j] = @$array[$j,$i];
        }
    }



```

# Density plot
We can use the output from the perl script above to make density plots:
```R
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2023_macaque_genomes/2024_macaques/general_genomix_Fst_pi')

# I made the permutation perl scripts print out TajD and FW_H values of genomic windows in one 
# column and whether or not the window contains an N_interact gene in the other column

# this script will make faceted density plots for each of the 5 populations for all N_interact

# read in the TajD data 
bor_sum <- read.table("./diversity.concat_FST__density.txt", header = T, sep="\t")
bor_ton <- read.table("./nem_ton_windowstats.concat_FST__density.txt", header = T, sep="\t")
bor_hec <- read.table("./nem_hec_windowstats.concat_FST__density.txt", header = T, sep="\t")
bor_nig <- read.table("./nem_nig_windowstats.concat_FST__density.txt", header = T, sep="\t")
mau_ton <- read.table("./mau_ton_windowstats.concat_FST__density.txt", header = T, sep="\t")
mau_hec <- read.table("./mau_hec_windowstats.concat_FST__density.txt", header = T, sep="\t")
mau_nig <- read.table("./mau_nig_windowstats.concat_FST__density.txt", header = T, sep="\t")
ton_hec <- read.table("./ton_hec_windowstats.concat_FST__density.txt", header = T, sep="\t")
ton_nig <- read.table("./ton_nig_windowstats.concat_FST__density.txt", header = T, sep="\t")
hec_nig <- read.table("./hec_nig_windowstats.concat_FST__density.txt", header = T, sep="\t")

# define a pairwise variable
bor_sum$pairwise <- "bor_sum"
bor_mau$pairwise <- "bor_mau"
bor_ton$pairwise <- "bor_ton"
bor_hec$pairwise <- "bor_hec"
bor_nig$pairwise <- "bor_nig"
mau_ton$pairwise <- "mau_ton"
mau_hec$pairwise <- "mau_hec"
mau_nig$pairwise <- "mau_nig"
ton_hec$pairwise <- "ton_hec"
ton_nig$pairwise <- "ton_nig"
hec_nig$pairwise <- "hec_nig"

# some quick statistics 
# independent 2-group Mann-Whitney U Test
wilcox.test(bor_sum$Fst~bor_sum$containsgenes)
wilcox.test(bor_mau$Fst~bor_mau$containsgenes)
wilcox.test(bor_ton$Fst~bor_ton$containsgenes)
wilcox.test(bor_hec$Fst~bor_hec$containsgenes)
wilcox.test(bor_nig$Fst~bor_nig$containsgenes)
wilcox.test(mau_ton$Fst~mau_ton$containsgenes)
wilcox.test(mau_hec$Fst~mau_hec$containsgenes)
wilcox.test(mau_nig$Fst~mau_nig$containsgenes)
wilcox.test(ton_hec$Fst~ton_hec$containsgenes)
wilcox.test(ton_nig$Fst~ton_nig$containsgenes)
wilcox.test(hec_nig$Fst~hec_nig$containsgenes)
# where y is numeric and A is A binary factor

# all are significant.
# sanity check to make sure that Fst is higher for genic compared to non-genic windows
only_nongenez <- bor_sum[which(bor_sum$containsgenes == 0),] ; mean(only_nongenez$Fst, na.rm = T)
only_genez <- bor_sum[which(bor_sum$containsgenes == 1),]; mean(only_genez$Fst, na.rm = T) 


only_nongenez <- hec_nig[which(hec_nig$containsgenes == 0),] ; mean(only_nongenez$Fst, na.rm = T)
only_genez <- hec_nig[which(hec_nig$containsgenes == 1),]; mean(only_genez$Fst, na.rm = T) 

only_nongenez <- bor_nig[which(bor_nig$containsgenes == 0),] ; mean(only_nongenez$Fst, na.rm = T)
only_genez <- bor_nig[which(bor_nig$containsgenes == 1),]; mean(only_genez$Fst, na.rm = T) 

my_data <- rbind(bor_sum)
# select only windows with genes
only_genez <- my_data[which(my_data$containsgenes != 0),] 
# independent 2-group Mann-Whitney U Test
wilcox.test(only_genez$TajD~only_genez$containsNinteractgenez)
wilcox.test(only_genez$FW_H~only_genez$containsNinteractgenez)
# where y is numeric and A is A binary factor



# now rbind them together
all_my_data <- rbind(bor_sum)
all_my_data <- rbind(bor_mau,bor_ton,bor_hec,bor_nig,mau_ton,mau_hec,mau_nig,ton_hec,ton_nig,hec_nig)

# define a grouping variable
all_my_data$group <- "NA"
# numbering the categories like this makes the one sided p-value sensible
all_my_data$group[which((all_my_data$containsgenes == 0)&(all_my_data$containsNinteractgenez == 0))] <- "No genes"
all_my_data$group[which((all_my_data$containsgenes == 1)&(all_my_data$containsNinteractgenez == 0))] <- "Other genes"
all_my_data$group[which((all_my_data$containsgenes == 1)&(all_my_data$containsNinteractgenez == 1))] <- "Ninteract genes"
table(all_my_data$group)

# re-order the species
# re-order the species
all_my_data$pairwise_f = factor(all_my_data$pairwise, levels=c('bor_sum'), ordered = T)
all_my_data$pairwise_f = factor(all_my_data$pairwise, levels=c('bor_mau','bor_ton','bor_hec','bor_nig',
                                                             'mau_ton','mau_hec','mau_nig','ton_hec',
                                                             'ton_nig','hec_nig'), ordered = T)
# reorder the group
all_my_data$group_f = factor(all_my_data$group, levels=c('No genes','Other genes','Ninteract genes'), ordered = T)


#Plot.
library(ggplot2)
png(paste("Fst_density_allNinteract.png",sep=""),
    width = 300, height = 50, units='mm', res = 100)
      
      ggplot(all_my_data) + 
      geom_density(aes(x=Fst, colour=group_f),show_guide=FALSE)+
      stat_density(aes(x=Fst, colour=group_f),
                   geom="line",position="identity")+
      scale_color_manual(values = c("black", "blue", "red")) +
      #facet_wrap(vars(pairwise_f), nrow=2, ncol = 5, scales = "free_x") +
      theme_bw() +
      # remove legend key border color & background
      theme(legend.key=element_blank()) +
      theme(legend.box.background = element_blank())+
      #theme(legend.key = element_rect(colour = NA, fill = NA)) +
      theme(legend.title = element_blank()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(strip.background = element_blank()) +
      theme(axis.title.x = element_blank()) +
      theme(
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 20),
          strip.text = element_text(size = 20),
          legend.text = element_text(size = 20)) +
          scale_x_continuous(limits = c(0.05, .12) #,
       # labels = scales::number_format(accuracy = 0.1)
       ) 
dev.off()

```
