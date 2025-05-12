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

# New ChrY with rhesus
I want to use the rhesus ref as an outgroup. I needed to filter biallelic sites and also save only the males. So I used bcftools in two steps like this:
```
bcftools view -S males_list.txt -M2 -v snps ../mac_chrY_concat.vcf.gz_filtered.vcf.gz > ../mac_chrY_concat.vcf.gz_filtered_onlybiallelic4real_onlymalez.vcf
```


Then I compressed this file and ran vcf2phylip like this:
```
python3 ../vcf2phylip/vcf2phylip.py -i ../mac_chrY_concat.vcf.gz_filtered_onlybiallelic4real_onlymalez.vcf.gz -m 0 --output-prefix theY_min0_ --write-used-sites --nexus
```

save the ref seq:
```
zcat ../mac_chrY_concat.vcf.gz_filtered_onlybiallelic4real_onlymalez.vcf.gz | cut -f4 > rhesus_seq.txt
````

find out where the REF is:
```
grep -n 'REF' rhesus_seq.txt
````

remove the first 2992 lines
```
sed '1,2992 d' rhesus_seq.txt > rhesus_seq_only.txt
```
get rid of the carriage returns:
```
sed ':a;N;$!ba;s/\n//g' rhesus_seq_only.txt > rhesus_seq_only_oneline.txt
```


And then I'll use iqtree to do model selection and ML analysis...
