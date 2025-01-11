# Permutations

This script quantifies how many Ninteract genes are on admixture blocks. It also calculates how many admixture blocks have a portion of at least one Ninteract gene.

It does a permutation that uses the proportion of Ninteract genes out of the total number of genes on admixture blocks as a test statistic.

Another way (that probably is quite similar) would be to use the number of admixture blocks with a proportion of an Ninteract gene out of the total number of admixture blocks with a portion of any gene as the test statistic.

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
# perl 2024_Introgression_block_permutation_homoz_introgression.pl rheMac10.ncbiRefSeq.CDS.marked.corrected.new.out 3 s105224_PM500_concatforperms.txt

# where the second argument (3 in the example above) is the species column that is the same as the mtDNA

# make the concatenated file like this:
# xzcat s105224_Chr1_MAU_TON_HEC.oout.bin.xz s105224_Chr2_MAU_TON_HEC.oout.bin.xz s105224_Chr3_MAU_TON_HEC.oout.bin.xz s105224_Chr4_MAU_TON_HEC.oout.bin.xz s105224_Chr5_MAU_TON_HEC.oout.bin.xz s105224_Chr6_MAU_TON_HEC.oout.bin.xz s105224_Chr7_MAU_TON_HEC.oout.bin.xz s105224_Chr8_MAU_TON_HEC.oout.bin.xz s105224_Chr9_MAU_TON_HEC.oout.bin.xz s105224_Chr10_MAU_TON_HEC.oout.bin.xz s105224_Chr11_MAU_TON_HEC.oout.bin.xz s105224_Chr12_MAU_TON_HEC.oout.bin.xz s105224_Chr13_MAU_TON_HEC.oout.bin.xz s105224_Chr14_MAU_TON_HEC.oout.bin.xz s105224_Chr15_MAU_TON_HEC.oout.bin.xz s105224_Chr16_MAU_TON_HEC.oout.bin.xz s105224_Chr17_MAU_TON_HEC.oout.bin.xz s105224_Chr18_MAU_TON_HEC.oout.bin.xz s105224_Chr19_MAU_TON_HEC.oout.bin.xz s105224_Chr20_MAU_TON_HEC.oout.bin.xz > s105224_PM500_concatforperms.txt



my $inputfile1 = $ARGV[0];
my $inputfile2 = $ARGV[1];
my $inputfile3 = $ARGV[2];
my $outputfile = $inputfile3."_perms.oout";

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
my $admixfrog_block_size=30000; # the 2024 admixfrog blocks are 10000 bp
my %introgression_blocks;
my @Ninteract_genez_on_one_or_more_admixture_block;
my %perm_blockz; 

while ( my $line = <DATAINPUT3>) {
	chomp($line);
	@temp=split(',',$line);
	# ignore first line
	if($temp[0] ne 'chrom'){
		# check if this is an introgression block
		# if($temp[8] > 0.5){ # this means introgression blocks must be homoz hecki (use for PF511)
		# if($temp[9] > 0.5){ # this means introgression blocks must be homoz hecki (use for PM500)	
		# if($temp[7] < 0.5){ # this means introgression blocks must not be homoz arc
		if($temp[6 + $inputfile2] > 0.5){ # this means N_interact blocks must be homoz hecki (use for PM500 which has hecki mtDNA)
			# print "hi ",$temp[6 + $inputfile2],"\n";
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
	my $test_stat2 = $n_Ninteract_genes_on_introgression_blockz/$n_genes_on_introgression_blockz;

	my @quick_perm_genez_sorted = sort { $a <=> $b } @quick_perm_genez;
	my $switch=0;
	my $pval=0;
	$counter=0;
	# print "@quick_perm_genez_sorted\n";
	# now figure out where the test stat is
	for ($y = 0 ; $y <= $#quick_perm_genez_sorted; $y++ ) {
		if(($test_stat2 <= $quick_perm_genez_sorted[$y])&&($switch==0)){ # use for PF626; test stat is zero
			print $counter," out of ",$perms," permutations have a smaller number of Ninteract genes in introgression blocks.\n";
			print OUTFILE $counter," out of ",$perms," permutations have a smaller number of Ninteract genes in introgression blocks.\n";
			$pval=$counter;
			$switch = 1;
		}
		$counter+=1;
	}
	if($switch==0){
		print $counter," out of ",$perms," permutations have a smaller number of Ninteract genes in introgression blocks.\n";
		print OUTFILE $counter," out of ",$perms," permutations have a smaller number of Ninteract genes in introgression blocks.\n";
	}
	# if all the perms are less than the test stat, then
	# we still need to set pval to be equal to the number of
	# perms
	if($counter == $perms){
	    $pval = $counter;
	}
		
	print "Test stat:",$test_stat2,"\n";
	print "QuickP = ",$pval/$perms,"\n";
	print OUTFILE "Test stat:",$test_stat2,"\n";
	print OUTFILE "QuickP = ",1-$pval/$perms,"\n";


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
