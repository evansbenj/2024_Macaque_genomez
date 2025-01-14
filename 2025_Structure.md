# Structure analysis

I'm going to use our vcf files as input to do Structure analysis using AdmixPipe3 (https://github.com/stevemussmann/admixturePipeline). I'd like to summarize multiple runs using CLUMPAK (https://clumpak.tau.ac.il/download.html).

To install the AdmixPipe3 pipeline (python3-based) on computecanada do this:
```
module load apptainer
apptainer build --sandbox bb.dir docker://mussmann/admixpipe:3.2
```
