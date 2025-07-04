# Permutations

This script quantifies how many Ninteract genes are on admixture blocks. It also calculates how many homozygous or heterozygous admixture blocks have a portion of at least one Ninteract gene. Below this I have another script that quantifies only homozygous blocks.

It does a permutation that uses the proportion of Ninteract genes out of the total number of genes on admixture blocks as a test statistic.

Another way (that probably is quite similar) would be to use the number of admixture blocks with a proportion of an Ninteract gene out of the total number of admixture blocks with a portion of any gene as the test statistic.

# Script for homoz and heteroz Ninteract genes


``` perl
#!/usr/bin/env perl
use strict;
use warnings;


# This program will read in two files. The first contains the coordinates of
# OXPHOS and non OXPHOS genes and their acronyms, and all other genes

# The other file is a file admix introgression blocks concatenated for all chrs.  

# First identify how many interacting OXPHOS genes are within introgression blocks
# then scramble them and check how many are expected by chance.

# run like this:
# perl 2024_Introgression_block_permutation_homoz_or_heteroz_introgression.pl rheMac10.ncbiRefSeq.CDS.marked.corrected.new.out MAU s105224_PM500_concatforperms.txt

# where the second argument (MAU in the example above) is the species from which the mtDNA is derived

# make the concatenated file like this:
# xzcat s105224_Chr1_MAU_TON_HEC.oout.bin.xz s105224_Chr2_MAU_TON_HEC.oout.bin.xz s105224_Chr3_MAU_TON_HEC.oout.bin.xz s105224_Chr4_MAU_TON_HEC.oout.bin.xz s105224_Chr5_MAU_TON_HEC.oout.bin.xz s105224_Chr6_MAU_TON_HEC.oout.bin.xz s105224_Chr7_MAU_TON_HEC.oout.bin.xz s105224_Chr8_MAU_TON_HEC.oout.bin.xz s105224_Chr9_MAU_TON_HEC.oout.bin.xz s105224_Chr10_MAU_TON_HEC.oout.bin.xz s105224_Chr11_MAU_TON_HEC.oout.bin.xz s105224_Chr12_MAU_TON_HEC.oout.bin.xz s105224_Chr13_MAU_TON_HEC.oout.bin.xz s105224_Chr14_MAU_TON_HEC.oout.bin.xz s105224_Chr15_MAU_TON_HEC.oout.bin.xz s105224_Chr16_MAU_TON_HEC.oout.bin.xz s105224_Chr17_MAU_TON_HEC.oout.bin.xz s105224_Chr18_MAU_TON_HEC.oout.bin.xz s105224_Chr19_MAU_TON_HEC.oout.bin.xz s105224_Chr20_MAU_TON_HEC.oout.bin.xz > s105224_PM500_concatforperms.txt



my $inputfile1 = $ARGV[0];
my $inputfile2 = $ARGV[1];
my $inputfile3 = $ARGV[2];
my $outputfile = $inputfile3."_hethomoz_".$inputfile2."_perms.oout";

print $outputfile,"\n";

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";

my @windowsites;
my @Fst_values;
my $sumsites=0;
my @temp;
my $y;
my $x;
my %N_interact_hash;
my @interact_perm;

# first open up the N_interact_hash gene info 
unless (open DATAINPUT, $inputfile1) {
	print "Can not find the input file1.\n";
	exit;
}

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] ne 'gene'){ 
		$N_interact_hash{$temp[2]."_".$temp[3]."_".$temp[4]}{"gene"} = $temp[0]; # key is chr_start_stop; value is gene_acronym
		$N_interact_hash{$temp[2]."_".$temp[3]."_".$temp[4]}{"mt_interact"} = $temp[5]; # key is chr_start_stop; value is Ninteract_or_not
		push(@interact_perm,$temp[5]); # this will be used for the permutations later
	}	
}		
close DATAINPUT;
#fisher_yates_shuffle( \@interact_perm );
# print "hello ",@interact_perm,"\n";

# now open up the introgression data
# consider an introgression window as
# any window with the homoz TON probability <0.5
unless (open DATAINPUT3, $inputfile3) {
	print "Can not find the input file3.\n";
	exit;
}

my @temp1;
my $n_introgression_blocks_with_interacting_genez=0;
my $n_introgression_blocks_with_other_genez=0;
my $n_introgression_blocks_without_genez=0;
my $n_Ninteract_genes_on_introgression_blockz=0;
my $n_genes_on_introgression_blockz=0;
my $n_genes_on_non_introgression_blockz=0;
my $admixfrog_block_size=30000; # the 2024 admixfrog blocks are 30000 bp
my %introgression_blocks;
my @Ninteract_genez_on_one_or_more_admixture_block;
my %perm_blockz; 

while ( my $line = <DATAINPUT3>) {
	chomp($line);
	@temp=split(',',$line);
	# ignore first line
	if($temp[0] ne 'chrom'){
		# check if this is an introgression block
		if($temp[5] =~ /$inputfile2/){ # this means N_interact blocks must be homoz or heteroz 
			$perm_blockz{$temp[0]."_".($temp[1]*1000000)}{"introgression"}=1; # this is an introgression block
			# this is an introgression block; 
			# block size is $admixfrog_block_size bp
			# initially assign the block to have no genes 
			$introgression_blocks{$temp[0]."_".($temp[1]*1000000)}{"genes"} = 0; # key has the lower limit of window
																	   # from admix frog.  the upper limit
																	   # is this plus $admixfrog_block_size
																	   # minus 1
			# also assume that it does not have any interacting genes
			$introgression_blocks{$temp[0]."_".($temp[1]*1000000)}{"interacting"} = 0;

			# Now cycle through all the genes to see if the beginning of any gene is in this block
			foreach my $key (keys %N_interact_hash){
				@temp1=split('_',$key);
				# print $temp[1]*1000000,"\n";
				# now check if this block contains any genes
				# print $temp1[0]," ",$temp1[1]," ",$temp1[2]," ",$temp[0]," ",($temp[1]*1000000)," ",($temp[1]*1000000+$admixfrog_block_size-1),"\n";
				if(
					($temp1[0] eq "chr".$temp[0])&&($temp1[1] >= ($temp[1]*1000000))&&($temp1[1] <= (($temp[1]*1000000+$admixfrog_block_size-1)))
					# gene beginning is in this block
					||
					($temp1[0] eq "chr".$temp[0])&&($temp1[2] >= ($temp[1]*1000000))&&($temp1[2] <= (($temp[1]*1000000+$admixfrog_block_size-1)))
					# gene end is in this block	
					||
					($temp1[0] eq "chr".$temp[0])&&($temp1[1] <= ($temp[1]*1000000))&&($temp1[2] >= (($temp[1]*1000000+$admixfrog_block_size-1)))
					# gene middle is in this block
					)
					{
						# this block has the beginning or end or middle of a gene						
						$introgression_blocks{$temp[0]."_".($temp[1]*1000000)}{"genes"} = 1;
						if(($temp1[0] eq "chr".$temp[0])&&($temp1[1] >= ($temp[1]*1000000))&&($temp1[1] <= ($temp[1]*1000000+$admixfrog_block_size-1))){	
							# for permutations count only the admixture blocks that contain the beginning of an Ninteract gene
							$n_genes_on_introgression_blockz+=1;
						}	# end if
						# check if it is an interacting gene
						if($N_interact_hash{$key}{"mt_interact"} == 1){
							$introgression_blocks{$temp[0]."_".($temp[1]*1000000)}{"interacting"} = 1; # this records overlap with any portion of an Ninteract genes 
							# key is chr and end; value is Ninteract or not
							print $N_interact_hash{$key}{"gene"}," is in introgression block ",$temp[0]."_".($temp[1]*1000000),"\n";
							print OUTFILE $N_interact_hash{$key}{"gene"}," is in introgression block ",$temp[0]."_".($temp[1]*1000000),"\n";
							$n_Ninteract_genes_on_introgression_blockz+=1; 
							# note that this will count one gene multiple times if it spans multiple admixture blocks
							push (@Ninteract_genez_on_one_or_more_admixture_block,$N_interact_hash{$key}{"gene"});
							# getting unique values from this array will provide the total unique Ninteract genes in one or more admixture blocks
						} # end if
						# no else statement because we don't want to erase that assignment if we had an interacting
						# gene in this block already
				} # end if
				else{
					$perm_blockz{$temp[0]."_".($temp[1]*1000000)}{"introgression"}=0; # this is not an introgression block
				} # end else
			} # end foreach
		}	# end test for introgression block
	} # end test for chrom	
} # endwhile

close DATAINPUT3;

# ok now I have a hash (%introgression_blocks) that has information on whether or not a block has any genes ($introgression_blocks{$key}{"genes"})
# and whether or not any of these genes have any N_interact genes ($introgression_blocks{$key}{"interacting"}).
# print these numbers

foreach my $key (keys %introgression_blocks){
	if($introgression_blocks{$key}{"genes"} == 1){ # there is a portion of a gene on this admixture block
		if($introgression_blocks{$key}{"interacting"} == 1){ # there is a portion of an Ninteract gene on on this admixture block
			$n_introgression_blocks_with_interacting_genez+=1;
			# print "Introgression_with_interacting ",$key,"\n";
		}
		else{
			$n_introgression_blocks_with_other_genez+=1;
			# print "Introgression_with_NON_interacting ",$key,"\n";
		}	
	}
	else{
		$n_introgression_blocks_without_genez+=1;
	}	
}

# get the unique Ninteract genes that are in one or more admixture blocks
print "This is the number of introgression blocks that overlap with at least a portion of an Ninteract genes ",$#Ninteract_genez_on_one_or_more_admixture_block+1,"\n"; 
print OUTFILE "This is the number of introgression blocks that overlap with at least a portion of an Ninteract genes ",$#Ninteract_genez_on_one_or_more_admixture_block+1,"\n"; 
my %seen;
@Ninteract_genez_on_one_or_more_admixture_block = grep { ! $seen{ $_ }++ } @Ninteract_genez_on_one_or_more_admixture_block;
print "This is the number of Ninteract genes that are on at least one admixture block ",$#Ninteract_genez_on_one_or_more_admixture_block+1,"\n"; 
print OUTFILE "This is the number of Ninteract genes that are on at least one admixture block ",$#Ninteract_genez_on_one_or_more_admixture_block+1,"\n"; 
print "This can be more than the number of admixture block with one or more N_interact genes because some admixture blocks can have more than one Ninteract gene\n"; 
print OUTFILE "This can be more than the number of admixture block with one or more N_interact genes because some admixture blocks can have more than one Ninteract gene\n"; 

print "These are the Ninteract genes that are on at least one admixture block  @Ninteract_genez_on_one_or_more_admixture_block\n"; 	
print OUTFILE "These are the Ninteract genes that are on at least one admixture block  @Ninteract_genez_on_one_or_more_admixture_block\n"; 	


print "Number of introgression blocks with one or more N_interact genes: ",$n_introgression_blocks_with_interacting_genez,"\n";
print OUTFILE"Number of introgression blocks with one or more N_interact genes: ",$n_introgression_blocks_with_interacting_genez,"\n";

print "Number of Ninteract genes on introgression blocks ",$n_Ninteract_genes_on_introgression_blockz,"\n";
print OUTFILE "Number of Ninteract genes on introgression blocks ",$n_Ninteract_genes_on_introgression_blockz,"\n";

print "Number of introgression blocks with other genes: ",$n_introgression_blocks_with_other_genez,"\n";
print "Number of introgression blocks without genes: ",$n_introgression_blocks_without_genez,"\n";
print OUTFILE "Number of introgression blocks with other genes: ",$n_introgression_blocks_with_other_genez,"\n";
print OUTFILE "Number of introgression blocks without genes: ",$n_introgression_blocks_without_genez,"\n";

print "Number of genes on introgression blocks: ",$n_genes_on_introgression_blockz,"\n";
print OUTFILE "Number of genes on introgression blocks: ",$n_genes_on_introgression_blockz,"\n";

print "Number of Ninteract genes with at least a portion on at least one admixture block: ",$#Ninteract_genez_on_one_or_more_admixture_block+1,"\n";
print OUTFILE "Number of Ninteract genes with at least a portion on at least one admixture block: ",$#Ninteract_genez_on_one_or_more_admixture_block+1,"\n";


if(($n_introgression_blocks_with_interacting_genez+$n_introgression_blocks_with_other_genez)>0){
	print "Proportion of introgression blocks with genes that have N_interact genes ",$n_introgression_blocks_with_interacting_genez/
	($n_introgression_blocks_with_interacting_genez+$n_introgression_blocks_with_other_genez),"\n";
	print "Proportion of genes on introgression blocks that are N_interact genes 	",$n_Ninteract_genes_on_introgression_blockz/$n_genes_on_introgression_blockz,"\n";
	print OUTFILE "Proportion of introgression blocks with genes that have N_interact genes ",$n_introgression_blocks_with_interacting_genez/
	($n_introgression_blocks_with_interacting_genez+$n_introgression_blocks_with_other_genez),"\n";
	print OUTFILE "Proportion of genes on introgression blocks that are N_interact genes 	",$n_Ninteract_genes_on_introgression_blockz/$n_genes_on_introgression_blockz,"\n";

#my $size = keys %N_interact_hash;
# print "hello ",$#interact_perm,"\n";

######################
# Permutations using only the beginning of Ninteract genes
######################

	my $perms=1000;
	my $counter=0;
	my $introg_interacter=0;
	my $introg_withgenez=0;
	my $introg_without_genez=0;
	my @permed_interactorz;
	my @permed_interactorz_gene_proportions;
	my $introg_withgenez_switch=0;
	my $introg_interacter_switch=0;
	my $genes_on_introgression_blocks=0;
	my $Ninteract_genes_on_introgression_blocks=0;


	my $Ngenez_on_blockz=0;
	my @quick_perm_genez;
	#print "yo ",$n_genes_on_introgression_blockz,"\n";
	#print "hey ",@interact_perm,"\n";
	for ($y = 0 ; $y < $perms; $y++ ) {
		$Ngenez_on_blockz=0;
		fisher_yates_shuffle( \@interact_perm );    # permutes the N_interact assignment for each gene
	# the quickest way is to just assume the first $n_genes_on_introgression_blockz of these genes are the 
	# ones on introgressin blocks
		for ($x = 0 ; $x < $n_genes_on_introgression_blockz; $x++ ) {
			if($interact_perm[$x] == 1){
				$Ngenez_on_blockz +=1;
			}
		}
		push(@quick_perm_genez,($Ngenez_on_blockz/$n_genes_on_introgression_blockz));	
	}	

	if($#quick_perm_genez != $perms-1){
		print "Hey, something wrong with perms\n";
	}
	my $test_stat2 = ($#Ninteract_genez_on_one_or_more_admixture_block+1)/$n_genes_on_introgression_blockz;

	my @quick_perm_genez_sorted = sort { $a <=> $b } @quick_perm_genez;
	my $switch=0;
	my $pval=0;
	$counter=0;

	#print "@quick_perm_genez_sorted\n";
	# now figure out where the test stat is
	for ($y = 0 ; $y <= $#quick_perm_genez_sorted; $y++ ) {
		if(($test_stat2 <= $quick_perm_genez_sorted[$y])&&($switch==0)){
			print $counter," out of ",$perms," permutations have a smaller number of Ninteract genes in introgression blocks.\n";
			print OUTFILE $counter," out of ",$perms," permutations have a smaller number of Ninteract genes in introgression blocks.\n";
			$pval=$counter;
			$switch = 1;
		}
		$counter+=1;
	}
	if($switch==0){ # this means that all of the perms were less than the test stat
		print $counter," out of ",$perms," permutations have a smaller number of Ninteract genes in introgression blocks.\n";
		print OUTFILE $counter," out of ",$perms," permutations have a smaller number of Ninteract genes in introgression blocks.\n";
		$pval = $counter;
	}

	print "Test stat:",$test_stat2,"\n";
	print "QuickP = ",1-($pval/$perms),"\n";
	print OUTFILE "Test stat:",$test_stat2,"\n";
	print OUTFILE "QuickP = ",1-($pval/$perms),"\n";


}
else{
	print "Permutations not performed because there are no introgression blocks that have genes.\n";
	print OUTFILE "Permutations not performed because there are no introgression blocks that have genes.\n";
}	



######################
# Another option would be to do the permutations
# using the proportion of admixture blocks that have 
# part of an Ninteract gene as the test statistic
######################



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


close OUTFILE;

```


# Script for only homoz Ninteract genes

``` perl
#!/usr/bin/env perl
use strict;
use warnings;


# This program will read in two files. The first contains the coordinates of
# OXPHOS and non OXPHOS genes and their acronyms, and all other genes

# The other file is a file admix introgression blocks concatenated for all chrs.  

# First identify how many interacting OXPHOS genes are within introgression blocks
# then scramble them and check how many are expected by chance.

# run like this:
# perl 2024_Introgression_block_permutation_homoz_only_introgression_new.pl rheMac10.ncbiRefSeq.CDS.marked.corrected.new.out MAU s105224_PM500_concatforperms.txt

# where the second argument (MAU in the example above) is the species from which the mtDNA is derived

# make the concatenated file like this:
# xzcat s105224_Chr1_MAU_TON_HEC.oout.bin.xz s105224_Chr2_MAU_TON_HEC.oout.bin.xz s105224_Chr3_MAU_TON_HEC.oout.bin.xz s105224_Chr4_MAU_TON_HEC.oout.bin.xz s105224_Chr5_MAU_TON_HEC.oout.bin.xz s105224_Chr6_MAU_TON_HEC.oout.bin.xz s105224_Chr7_MAU_TON_HEC.oout.bin.xz s105224_Chr8_MAU_TON_HEC.oout.bin.xz s105224_Chr9_MAU_TON_HEC.oout.bin.xz s105224_Chr10_MAU_TON_HEC.oout.bin.xz s105224_Chr11_MAU_TON_HEC.oout.bin.xz s105224_Chr12_MAU_TON_HEC.oout.bin.xz s105224_Chr13_MAU_TON_HEC.oout.bin.xz s105224_Chr14_MAU_TON_HEC.oout.bin.xz s105224_Chr15_MAU_TON_HEC.oout.bin.xz s105224_Chr16_MAU_TON_HEC.oout.bin.xz s105224_Chr17_MAU_TON_HEC.oout.bin.xz s105224_Chr18_MAU_TON_HEC.oout.bin.xz s105224_Chr19_MAU_TON_HEC.oout.bin.xz s105224_Chr20_MAU_TON_HEC.oout.bin.xz > s105224_PM500_concatforperms.txt



my $inputfile1 = $ARGV[0];
my $inputfile2 = $ARGV[1];
my $inputfile3 = $ARGV[2];
my $outputfile = $inputfile3."_homoz_only_".$inputfile2."_perms.oout";

print $outputfile,"\n";

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";

my @windowsites;
my @Fst_values;
my $sumsites=0;
my @temp;
my $y;
my $x;
my %N_interact_hash;
my @interact_perm;

# first open up the N_interact_hash gene info 
unless (open DATAINPUT, $inputfile1) {
	print "Can not find the input file1.\n";
	exit;
}

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] ne 'gene'){ 
		$N_interact_hash{$temp[2]."_".$temp[3]."_".$temp[4]}{"gene"} = $temp[0]; # key is chr_start_stop; value is gene_acronym
		$N_interact_hash{$temp[2]."_".$temp[3]."_".$temp[4]}{"mt_interact"} = $temp[5]; # key is chr_start_stop; value is Ninteract_or_not
		push(@interact_perm,$temp[5]); # this will be used for the permutations later
	}	
}		
close DATAINPUT;
#fisher_yates_shuffle( \@interact_perm );
# print "hello ",@interact_perm,"\n";

# now open up the introgression data
# consider an introgression window as
# any window with the homoz TON probability <0.5
unless (open DATAINPUT3, $inputfile3) {
	print "Can not find the input file3.\n";
	exit;
}

my @temp1;
my $n_introgression_blocks_with_interacting_genez=0;
my $n_introgression_blocks_with_other_genez=0;
my $n_introgression_blocks_without_genez=0;
my $n_Ninteract_genes_on_introgression_blockz=0;
my $n_genes_on_introgression_blockz=0;
my $n_genes_on_non_introgression_blockz=0;
my $admixfrog_block_size=30000; # the 2024 admixfrog blocks are 30000 bp
my %introgression_blocks;
my @Ninteract_genez_on_one_or_more_admixture_block;
my %perm_blockz; 

while ( my $line = <DATAINPUT3>) {
	chomp($line);
	@temp=split(',',$line);
	# ignore first line
	if($temp[0] ne 'chrom'){
		# check if this is an introgression block
		if($temp[5] eq $inputfile2){ # this means N_interact blocks must be homoz 
			$perm_blockz{$temp[0]."_".($temp[1]*1000000)}{"introgression"}=1; # this is an introgression block
			# this is an introgression block; 
			# block size is $admixfrog_block_size bp
			# initially assign the block to have no genes 
			$introgression_blocks{$temp[0]."_".($temp[1]*1000000)}{"genes"} = 0; # key has the lower limit of window
																	   # from admix frog.  the upper limit
																	   # is this plus $admixfrog_block_size
																	   # minus 1
			# also assume that it does not have any interacting genes
			$introgression_blocks{$temp[0]."_".($temp[1]*1000000)}{"interacting"} = 0;

			# Now cycle through all the genes to see if the beginning of any gene is in this block
			foreach my $key (keys %N_interact_hash){
				@temp1=split('_',$key);
				# print $temp[1]*1000000,"\n";
				# now check if this block contains any genes
				# print $temp1[0]," ",$temp1[1]," ",$temp1[2]," ",$temp[0]," ",($temp[1]*1000000)," ",($temp[1]*1000000+$admixfrog_block_size-1),"\n";
				if(
					($temp1[0] eq "chr".$temp[0])&&($temp1[1] >= ($temp[1]*1000000))&&($temp1[1] <= (($temp[1]*1000000+$admixfrog_block_size-1)))
					# gene beginning is in this block
					||
					($temp1[0] eq "chr".$temp[0])&&($temp1[2] >= ($temp[1]*1000000))&&($temp1[2] <= (($temp[1]*1000000+$admixfrog_block_size-1)))
					# gene end is in this block	
					||
					($temp1[0] eq "chr".$temp[0])&&($temp1[1] <= ($temp[1]*1000000))&&($temp1[2] >= (($temp[1]*1000000+$admixfrog_block_size-1)))
					# gene middle is in this block
					)
					{
						# this block has the beginning or end or middle of a gene						
						$introgression_blocks{$temp[0]."_".($temp[1]*1000000)}{"genes"} = 1;
						if(($temp1[0] eq "chr".$temp[0])&&($temp1[1] >= ($temp[1]*1000000))&&($temp1[1] <= ($temp[1]*1000000+$admixfrog_block_size-1))){	
							# for permutations count only the admixture blocks that contain the beginning of an Ninteract gene
							$n_genes_on_introgression_blockz+=1;
						}	# end if
						# check if it is an interacting gene
						if($N_interact_hash{$key}{"mt_interact"} == 1){
							$introgression_blocks{$temp[0]."_".($temp[1]*1000000)}{"interacting"} = 1; # this records overlap with any portion of an Ninteract genes 
							# key is chr and end; value is Ninteract or not
							print $N_interact_hash{$key}{"gene"}," is in introgression block ",$temp[0]."_".($temp[1]*1000000),"\n";
							print OUTFILE $N_interact_hash{$key}{"gene"}," is in introgression block ",$temp[0]."_".($temp[1]*1000000),"\n";
							$n_Ninteract_genes_on_introgression_blockz+=1; 
							# note that this will count one gene multiple times if it spans multiple admixture blocks
							push (@Ninteract_genez_on_one_or_more_admixture_block,$N_interact_hash{$key}{"gene"});
							# getting unique values from this array will provide the total unique Ninteract genes in one or more admixture blocks
						} # end if
						# no else statement because we don't want to erase that assignment if we had an interacting
						# gene in this block already
				} # end if
				else{
					$perm_blockz{$temp[0]."_".($temp[1]*1000000)}{"introgression"}=0; # this is not an introgression block
				} # end else
			} # end foreach
		}	# end test for introgression block
	} # end test for chrom	
} # endwhile

close DATAINPUT3;

# ok now I have a hash (%introgression_blocks) that has information on whether or not a block has any genes ($introgression_blocks{$key}{"genes"})
# and whether or not any of these genes have any N_interact genes ($introgression_blocks{$key}{"interacting"}).
# print these numbers

foreach my $key (keys %introgression_blocks){
	if($introgression_blocks{$key}{"genes"} == 1){ # there is a portion of a gene on this admixture block
		if($introgression_blocks{$key}{"interacting"} == 1){ # there is a portion of an Ninteract gene on on this admixture block
			$n_introgression_blocks_with_interacting_genez+=1;
			# print "Introgression_with_interacting ",$key,"\n";
		}
		else{
			$n_introgression_blocks_with_other_genez+=1;
			# print "Introgression_with_NON_interacting ",$key,"\n";
		}	
	}
	else{
		$n_introgression_blocks_without_genez+=1;
	}	
}

# get the unique Ninteract genes that are in one or more admixture blocks
print "This is the number of introgression blocks that overlap with at least a portion of an Ninteract genes ",$#Ninteract_genez_on_one_or_more_admixture_block+1,"\n"; 
print OUTFILE "This is the number of introgression blocks that overlap with at least a portion of an Ninteract genes ",$#Ninteract_genez_on_one_or_more_admixture_block+1,"\n"; 
my %seen;
@Ninteract_genez_on_one_or_more_admixture_block = grep { ! $seen{ $_ }++ } @Ninteract_genez_on_one_or_more_admixture_block;
print "This is the number of Ninteract genes that are on at least one admixture block ",$#Ninteract_genez_on_one_or_more_admixture_block+1,"\n"; 
print OUTFILE "This is the number of Ninteract genes that are on at least one admixture block ",$#Ninteract_genez_on_one_or_more_admixture_block+1,"\n"; 
print "This can be more than the number of admixture block with one or more N_interact genes because some admixture blocks can have more than one Ninteract gene\n"; 
print OUTFILE "This can be more than the number of admixture block with one or more N_interact genes because some admixture blocks can have more than one Ninteract gene\n"; 

print "These are the Ninteract genes that are on at least one admixture block  @Ninteract_genez_on_one_or_more_admixture_block\n"; 	
print OUTFILE "These are the Ninteract genes that are on at least one admixture block  @Ninteract_genez_on_one_or_more_admixture_block\n"; 	


print "Number of introgression blocks with one or more N_interact genes: ",$n_introgression_blocks_with_interacting_genez,"\n";
print OUTFILE"Number of introgression blocks with one or more N_interact genes: ",$n_introgression_blocks_with_interacting_genez,"\n";

print "Number of Ninteract genes on introgression blocks ",$n_Ninteract_genes_on_introgression_blockz,"\n";
print OUTFILE "Number of Ninteract genes on introgression blocks ",$n_Ninteract_genes_on_introgression_blockz,"\n";

print "Number of introgression blocks with other genes: ",$n_introgression_blocks_with_other_genez,"\n";
print "Number of introgression blocks without genes: ",$n_introgression_blocks_without_genez,"\n";
print OUTFILE "Number of introgression blocks with other genes: ",$n_introgression_blocks_with_other_genez,"\n";
print OUTFILE "Number of introgression blocks without genes: ",$n_introgression_blocks_without_genez,"\n";

print "Number of genes on introgression blocks: ",$n_genes_on_introgression_blockz,"\n";
print OUTFILE "Number of genes on introgression blocks: ",$n_genes_on_introgression_blockz,"\n";

print "Number of Ninteract genes with at least a portion on at least one admixture block: ",$#Ninteract_genez_on_one_or_more_admixture_block+1,"\n";
print OUTFILE "Number of Ninteract genes with at least a portion on at least one admixture block: ",$#Ninteract_genez_on_one_or_more_admixture_block+1,"\n";


if(($n_introgression_blocks_with_interacting_genez+$n_introgression_blocks_with_other_genez)>0){
	print "Proportion of introgression blocks with genes that have N_interact genes ",$n_introgression_blocks_with_interacting_genez/
	($n_introgression_blocks_with_interacting_genez+$n_introgression_blocks_with_other_genez),"\n";
	print "Proportion of genes on introgression blocks that are N_interact genes 	",$n_Ninteract_genes_on_introgression_blockz/$n_genes_on_introgression_blockz,"\n";
	print OUTFILE "Proportion of introgression blocks with genes that have N_interact genes ",$n_introgression_blocks_with_interacting_genez/
	($n_introgression_blocks_with_interacting_genez+$n_introgression_blocks_with_other_genez),"\n";
	print OUTFILE "Proportion of genes on introgression blocks that are N_interact genes 	",$n_Ninteract_genes_on_introgression_blockz/$n_genes_on_introgression_blockz,"\n";

#my $size = keys %N_interact_hash;
# print "hello ",$#interact_perm,"\n";

######################
# Permutations using only the beginning of Ninteract genes
######################

	my $perms=1000;
	my $counter=0;
	my $introg_interacter=0;
	my $introg_withgenez=0;
	my $introg_without_genez=0;
	my @permed_interactorz;
	my @permed_interactorz_gene_proportions;
	my $introg_withgenez_switch=0;
	my $introg_interacter_switch=0;
	my $genes_on_introgression_blocks=0;
	my $Ninteract_genes_on_introgression_blocks=0;


	my $Ngenez_on_blockz=0;
	my @quick_perm_genez;
	#print "yo ",$n_genes_on_introgression_blockz,"\n";
	#print "hey ",@interact_perm,"\n";
	for ($y = 0 ; $y < $perms; $y++ ) {
		$Ngenez_on_blockz=0;
		fisher_yates_shuffle( \@interact_perm );    # permutes the N_interact assignment for each gene
	# the quickest way is to just assume the first $n_genes_on_introgression_blockz of these genes are the 
	# ones on introgressin blocks
		for ($x = 0 ; $x < $n_genes_on_introgression_blockz; $x++ ) {
			if($interact_perm[$x] == 1){
				$Ngenez_on_blockz +=1;
			}
		}
		push(@quick_perm_genez,($Ngenez_on_blockz/$n_genes_on_introgression_blockz));	
	}	

	if($#quick_perm_genez != $perms-1){
		print "Hey, something wrong with perms\n";
	}
	my $test_stat2 = ($#Ninteract_genez_on_one_or_more_admixture_block+1)/$n_genes_on_introgression_blockz;

	my @quick_perm_genez_sorted = sort { $a <=> $b } @quick_perm_genez;
	my $switch=0;
	my $pval=0;
	$counter=0;

	#print "@quick_perm_genez_sorted\n";
	# now figure out where the test stat is
	for ($y = 0 ; $y <= $#quick_perm_genez_sorted; $y++ ) {
		if(($test_stat2 <= $quick_perm_genez_sorted[$y])&&($switch==0)){
			print $counter," out of ",$perms," permutations have a smaller number of Ninteract genes in introgression blocks.\n";
			print OUTFILE $counter," out of ",$perms," permutations have a smaller number of Ninteract genes in introgression blocks.\n";
			$pval=$counter;
			$switch = 1;
		}
		$counter+=1;
	}
	if($switch==0){ # this means that all of the perms were less than the test stat
		print $counter," out of ",$perms," permutations have a smaller number of Ninteract genes in introgression blocks.\n";
		print OUTFILE $counter," out of ",$perms," permutations have a smaller number of Ninteract genes in introgression blocks.\n";
		$pval = $counter;
	}

	print "Test stat:",$test_stat2,"\n";
	print "QuickP = ",1-($pval/$perms),"\n";
	print OUTFILE "Test stat:",$test_stat2,"\n";
	print OUTFILE "QuickP = ",1-($pval/$perms),"\n";


}
else{
	print "Permutations not performed because there are no introgression blocks that have genes.\n";
	print OUTFILE "Permutations not performed because there are no introgression blocks that have genes.\n";
}	



######################
# Another option would be to do the permutations
# using the proportion of admixture blocks that have 
# part of an Ninteract gene as the test statistic
######################



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


close OUTFILE;

```

# Tabulating Ninteract genes

Of interest is whether there are some Ninteract genes that tend to be in introgression blocks more than others. I used grep to quantify the presence of Ninteract genes in introgression blocks by searching through the output files generated from the scripts above:
```
grep -c 'TFB2M'  het.Ninteract.txt > het_counts.txt
grep -c 'COX20'  het.Ninteract.txt >> het_counts.txt
grep -c 'DARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS14'  het.Ninteract.txt >> het_counts.txt
grep -c 'PDC'  het.Ninteract.txt >> het_counts.txt
grep -c 'IARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL55'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL24'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL9'  het.Ninteract.txt >> het_counts.txt
grep -c 'TARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS21'  het.Ninteract.txt >> het_counts.txt
grep -c 'WARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5PB'  het.Ninteract.txt >> het_counts.txt
grep -c 'PARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL37'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATPAF1'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFS5'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS15'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5IF1'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL20'  het.Ninteract.txt >> het_counts.txt
grep -c 'LOC716161'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL40'  het.Ninteract.txt >> het_counts.txt
grep -c 'PET117'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS26'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5F1E'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL51'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS35'  het.Ninteract.txt >> het_counts.txt
grep -c 'YARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'COX14'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5F1B'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL42'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFA12'  het.Ninteract.txt >> het_counts.txt
grep -c 'TMEM177'  het.Ninteract.txt >> het_counts.txt
grep -c 'MARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFB3'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFS1'  het.Ninteract.txt >> het_counts.txt
grep -c 'PNKD'  het.Ninteract.txt >> het_counts.txt
grep -c 'BCS1L'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL44'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFA10'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS5'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL30'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS9'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL35'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL19'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL53'  het.Ninteract.txt >> het_counts.txt
grep -c 'COX7C'  het.Ninteract.txt >> het_counts.txt
grep -c 'LRPPRC'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL33'  het.Ninteract.txt >> het_counts.txt
grep -c 'ACP1'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL23'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL21'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFV1'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL11'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL49'  het.Ninteract.txt >> het_counts.txt
grep -c 'COX8A'  het.Ninteract.txt >> het_counts.txt
grep -c 'UQCC3'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL16'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL17'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL48'  het.Ninteract.txt >> het_counts.txt
grep -c 'KCTD14'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFC2'  het.Ninteract.txt >> het_counts.txt
grep -c 'NARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'TMEM126B'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5MG'  het.Ninteract.txt >> het_counts.txt
grep -c 'FOXRED1'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL41'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'LOC722212'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFA8'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL50'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFB6'  het.Ninteract.txt >> het_counts.txt
grep -c 'DMAC1'  het.Ninteract.txt >> het_counts.txt
grep -c 'TMEM220'  het.Ninteract.txt >> het_counts.txt
grep -c 'COX10'  het.Ninteract.txt >> het_counts.txt
grep -c 'TTC19'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATPAF2'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS23'  het.Ninteract.txt >> het_counts.txt
grep -c 'COX11'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL27'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5MC1'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL10'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL45'  het.Ninteract.txt >> het_counts.txt
grep -c 'CCDC56'  het.Ninteract.txt >> het_counts.txt
grep -c 'TACO1'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL58'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5PD'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS7'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL38'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL12'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL57'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS31'  het.Ninteract.txt >> het_counts.txt
grep -c 'CARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5F1A'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFV2'  het.Ninteract.txt >> het_counts.txt
grep -c 'POLRMT'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5F1D'  het.Ninteract.txt >> het_counts.txt
grep -c 'UQCR'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL54'  het.Ninteract.txt >> het_counts.txt
grep -c 'PET100'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFA7'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL4'  het.Ninteract.txt >> het_counts.txt
grep -c 'ECSIT'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFB7'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL34'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFA13'  het.Ninteract.txt >> het_counts.txt
grep -c 'COX6B1'  het.Ninteract.txt >> het_counts.txt
grep -c 'SARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS12'  het.Ninteract.txt >> het_counts.txt
grep -c 'DMAC2'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFA3'  het.Ninteract.txt >> het_counts.txt
grep -c 'CMC1'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL3'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS22'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL47'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFB5'  het.Ninteract.txt >> het_counts.txt
grep -c 'HIGD1A'  het.Ninteract.txt >> het_counts.txt
grep -c 'LARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFAF3'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS25'  het.Ninteract.txt >> het_counts.txt
grep -c 'ACAD9'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFB4'  het.Ninteract.txt >> het_counts.txt
grep -c 'COX17'  het.Ninteract.txt >> het_counts.txt
grep -c 'TIMMDC1'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL28'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS34'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFB10'  het.Ninteract.txt >> het_counts.txt
grep -c 'TMEM186'  het.Ninteract.txt >> het_counts.txt
grep -c 'EARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFAB1'  het.Ninteract.txt >> het_counts.txt
grep -c 'COX4I1'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFV3'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS6'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5PO'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5PF'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL39'  het.Ninteract.txt >> het_counts.txt
grep -c 'COX19'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5MF'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS17'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS24'  het.Ninteract.txt >> het_counts.txt
grep -c 'COA1'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL32'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFA4'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFB2'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS33'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL18'  het.Ninteract.txt >> het_counts.txt
grep -c 'TFB1M'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFAF4'  het.Ninteract.txt >> het_counts.txt
grep -c 'RARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'AARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL14'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS18A'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL2'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS10'  het.Ninteract.txt >> het_counts.txt
grep -c 'C6orf125'  het.Ninteract.txt >> het_counts.txt
grep -c 'VARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS18B'  het.Ninteract.txt >> het_counts.txt
grep -c 'FARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5ME'  het.Ninteract.txt >> het_counts.txt
grep -c 'SMIM20'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS18C'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL1'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFC1'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL36'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFS6'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS30'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFS4'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFAF2'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS36'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS27'  het.Ninteract.txt >> het_counts.txt
grep -c 'LYRM7'  het.Ninteract.txt >> het_counts.txt
grep -c 'UQCRQ'  het.Ninteract.txt >> het_counts.txt
grep -c 'CD14'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFA2'  het.Ninteract.txt >> het_counts.txt
grep -c 'HARS2'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL22'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFAF1'  het.Ninteract.txt >> het_counts.txt
grep -c 'COX5A'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL46'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS11'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL52'  het.Ninteract.txt >> het_counts.txt
grep -c 'DMAC2L'  het.Ninteract.txt >> het_counts.txt
grep -c 'COX16'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFB1'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5MPL'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL15'  het.Ninteract.txt >> het_counts.txt
grep -c 'TMEM70'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS28'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFAF6'  het.Ninteract.txt >> het_counts.txt
grep -c 'COX6C'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL13'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFB9'  het.Ninteract.txt >> het_counts.txt
grep -c 'CYC1'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5F1C'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPS16'  het.Ninteract.txt >> het_counts.txt
grep -c 'TFAM'  het.Ninteract.txt >> het_counts.txt
grep -c 'NDUFB8'  het.Ninteract.txt >> het_counts.txt
grep -c 'MRPL43'  het.Ninteract.txt >> het_counts.txt
grep -c 'ATP5MD'  het.Ninteract.txt >> het_counts.txt
```
I should make some density plots or histograms to see what these distributions look like...
