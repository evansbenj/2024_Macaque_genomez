# First get longest transcript 

(from support):
```
This software can be used with Singularity (which is Apptainer on our clusters):
https://github.com/NBISweden/AGAT#using-singularity
If you run (on a login node):
module load apptainer
apptainer build agat-1.7.0.sif docker://quay.io/biocontainers/agat:1.7.0--pl5321hdfd78af_0
If that is the version you want to use.
apptainer run agat-1.7.0.sif
And once in the container you can run commands. That said, its instructions are interactive. To submit jobs, you will need to run the needed commands in a non-interactive way, e.g., in your job script:
module load apptainer
apptainer run agat-1.7.0.sif <COMMAND and ARGS...>
e.g., apptainer run agat-1.7.0.sif agat_convert_sp_gxf2gxf.pl --help
where <COMMAND and ARGS...> is whatever command you need to run.
```
Here is the command I used:

```
apptainer run agat-1.7.0.sif agat_sp_keep_longest_isoform.pl \
  --gff /home/ben/projects/rrg-ben/ben/2024_macaques/rheMac10/GCF_003339765.1_Mmul_10_genomic.gff \
  -o /home/ben/projects/rrg-ben/ben/2024_macaques/rheMac10/GCF_003339765.1_Mmul_10_genomic_longest_isoforms.gff
```


# Extract CDS, concatenate, and translate each gene; then align and calculate Fst

```
All linez
### set the gene name
gff=AARS2
#### get information about a gene (need to check if there are multiple transcripts - if yes, then remove some of the CDS)
grep ${gff} /home/ben/projects/rrg-ben/ben/2024_macaques/rheMac10/GCF_003339765.1_Mmul_10_genomic_longest_isoforms.gff | grep 'CDS' > ${gff}.gff
#### get chromosome:
silly_chr=$(head -n1 ${gff}.gff | cut -f1)
better_chr=$(grep "${silly_chr}" /home/ben/projects/rrg-ben/ben/2024_macaques/rheMac10/GCF_003339765.1_Mmul_10_genomic.gff_codez.txt | cut -f1)
#### make bed out of gff
source make_bed_out_of_gff.sh ${gff}
### extract the CDS coordinates
source get_span_from_gff.sh ${gff}.gff
### Get all sample names
# --> bcftools query -l mac_chr${better_chr}_concat.vcf.gz > samples.txt
### Change the $REGION environmental variable to match chr name (need to change chromosome number!)
### Transfer the chromosome names
REGION="${better_chr}:${REGION#*:}"
### Generate haplotype consensus sequences (need to change chromosome number!)
source get_the_haplotypez.sh mac_chr${better_chr}_concat.vcf.gz ${gff}
### Convert CDS coordinates into region-relative coordinates
source get_relative_bed.sh ${gff}
### Change the name of the chromosome
sed -i "s/${silly_chr}/${REGION}/g" ${gff}.relative.bed
### Extract CDS from every sample
source extract_CDZ_using_relative_bed.sh ${gff}
### Concatenate CDS exons
source concat_cdz_exonz.sh ${gff}
### Translate CDS to protein
source translate_dat_seq.sh ${gff}
### Combine proteins into a single FASTA
cat *.${gff}.hap[1,2].region.pep.fa > ${gff}_all_proteins.fa
### Remove stop codon
sed -i 's/*//g' ${gff}_all_proteins.fa
### align
module load StdEnv/2023 mafft/7.526
mafft ${gff}_all_proteins.fa > ${gff}_all_proteins_aligned.fa
### make a directory and move everything
mkdir ${gff}; mv *${gff}* ${gff}
### calculate pi
cd ${gff}
perl ../aa_pi_per_site.pl --fasta ${gff}_all_proteins_aligned.fa --subset "$( < ../MAU_binary.txt )"
perl ../aa_pi_per_site.pl --fasta ${gff}_all_proteins_aligned.fa --subset "$( < ../TON_binary.txt )"
perl ../aa_pi_per_site.pl --fasta ${gff}_all_proteins_aligned.fa --subset "$( < ../MAUTON_binary.txt )" 
```

# Scripts:
