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
