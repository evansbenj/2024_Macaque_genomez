# ChrY and ChrM

Directory on graham (ChrY):
```
/home/ben/projects/rrg-ben/ben/2024_macaques/concatenated_vcfs/all_hardfiltered/ChrY
```

Directory on graham (ChrM):
```
/home/ben/projects/rrg-ben/ben/2024_macaques/concatenated_vcfs/all_hardfiltered/ChrM
```


I'm going to make a ChrY phylogeny from the male individual maqz. I'm using this python script (https://github.com/edgardomortiz/vcf2phylip/tree/master) to generate a nexus file:
```
python3 vcf2phylip/vcf2phylip.py -i all_162_maqs_chrY.vcf_malez_only.vcf -n --output-folder ChrY --output-prefix malz_
```
The ChrY nexus file is on info2020 here:
```
/home/ben/2024_macaques/all_hardfiltered/ChrY
```

And then I'll use iqtree to do model selection and ML analysis...
