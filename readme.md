<!--pandoc
t: html
toc:
s:
self-contained:
highlight-style:tango
-->



# Overview of processing scripts

SuRE data is processed using 4 bash scripts. The scripts implement the entire data processing pipeline for both iPCR data and cDNA data from plasmid library samples and RNA expression samples. The 4 scripts roughly perform the following functions:

- **cDNA-plDNA-count-BC.bash:** parse cDNA fastq files into count tables

- **iPCR-map-BC.bash:** parse iPCR fastq to bed-like files per sample

- **iPCR-merge-bedpe-Filter-BC-multi-pos.bash:** merge multiple iPCR bed-like files into 1

- **merge-iPCR-cDNA-plDNA.bash:** merge cDNA count data and iPCR bed-like file

## Required software
The scripts use various software tools, some standard on most linux
environments but most should be installed prior to using these scripts. A
complete list is shown in the following table:

| name          | version                |
|---------------|------------------------|
| gawk          | 4.0.1                  |
| cutadapt      | 1.9.1                  |
| gzip          | -                      |
| parallel (GNU)| -                      |
| bowtie2       | 2.1.0                  |
| samtools      | 1.2                    |
| python        | 2.7.6                  |


## cDNA-plDNA-count-BC.bash

```
DESCRIPTION:
   This is a bash script (mostly gawk and cutadapt) to process raw fastq files containing
   data from cDNA/plDNA samples. The barcodes are extracted from the reads and
   counted. For exatrcting the barcodes the adapter sequence is aligned with
   the read (using cutadapt) and the preceding part of the read is defined as
   the barcode.\\
   Barcodes with length != 20, or which contain N's are discarded. The
   filtered barcodes and counts are written to stdout, sorted (alphabetically)
   on barcode.
USAGE/OPTIONS:
     cDNA-plDNA-count-BC.bash [options] fastqfiles
   required:
     -o: output directory
     -a: adapter sequence
   optional:
     -l: log-filename [stdout]
     -n: number of cores used in parallel processes [10]
INPUT:
   cDNA/plDNA fastq files
OUTPUT:
   tabular txt files (one for each input file) with count and barcode-sequence
```





# Notes

## Filename pattern

Some scripts take fastq filenames as input. Different samples are named using
these fastq filenames. The filenames are assumed to contain certain string
patterns. If your filenames are not compatible you can either adjust the
scripts or your filenames. This dependency on string patterns is in the
scripts:

- **cDNA-plDNA-count-BC.bash:** lines 144-149
-  **iPCR-map-BC.bash:** lines 179-187

