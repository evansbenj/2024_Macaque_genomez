# Admixfrog

Admixfrog is a program written by Ben Peter that uses an HMM to estimate locations of introgression blocks in genomic data.

I installed the latest version (admixfrog 0.7.2.post1.dev0+1c27c43) like this:
```
module load StdEnv/2023 python/3.12.4
pip install cython scipy --upgrade --user
pip install git+https://github.com/benjaminpeter/admixfrog --user
```

and tested like this:

```
/home/ben/.local/bin/admixfrog --help
```
