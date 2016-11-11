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

## Additional software
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

