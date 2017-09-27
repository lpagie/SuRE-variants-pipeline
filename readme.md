<!--pandoc
t: html
toc:
s:
self-contained:
highlight-style:tango
-->

# SUBMODULE ALERT!!!!!

This pipeline uses WASP (https://github.com/bmvdgeijn/WASP) for some of the SNP
related data processing. I have a forked/adapted version of WASP
(https://github.com/lpagie/WASP) which I added to this repos as a submodule.

Most importantly: in order to clone this repos use the following command:
```
git clone --recursive https://github.com/lpagie/SuRE-pipeline.git
```

Also note that you will not be notified of changes in the submodule. ( I don't
expect changes, so this is just a disclaimer ;-)


# Overview of processing steps

SuRE data consists of potentially 3 types of NGS data:
- plasmid library cDNA data
- RNA cDNA data
- iPCR data

## cDNA data

The cDNA data generated from plasmid libraries or RNA samples are processed
similarly. The reads are parsed to extract the sequences representing barcodes,
the barcodes are counted; tables of the barcode counts are generated.

## iPCR data

The iPCR data is processed in multiple steps:

- Parsing raw reads
  The read sequences are parsed to extract the barcode sequences and the gDNA sequences
- Aligning gDNA sequences to genomic reference
- Assciating barcodes and genomic positions
- Filtering barcode/position pairs
  Used criteria are, e.g, (relative) frequency of occurrence of a
  barcode/position pair, relative to other positions with the same barcode, or
  other barcodes at the same position.
- Annotation of (raw) gDNA sequences with base identites at SNP positions.
- Generate tables of 'SuRE fragments, annotated by:
  - barcode sequence
  - genomic position
  - observed base identities at SNP positions

## Merging cDNA and iPCR data

The data are merged in  table with count columns for all iPCR and cDNA samples


# Pipeline

The pipeline is build using
[snakemake](https://snakemake.readthedocs.io/en/stable/), managing a set of
bash and python scripts. The scripts use additional, external software. All
software is managed using [(bio-)conda](https://bioconda.github.io/). A [conda
environment
file](https://github.com/lpagie/SuRE-pipeline/blob/snakemake/code/conda-wasp-environment.yml)
specifies the full set of software required.

## Setup of pipeline

The git repository contains a
[Makefile](https://github.com/lpagie/SuRE-pipeline/blob/snakemake/Makefile)
which (supposedly) installs the pipeline with all required additional software.
Before running the makefile bioconda should be installed first.

Snakemake itself is not included in the conda environment file (not sure why
not), so it should be installed separately as well.

# Running the pipeline

Before running the pipeline activate the conda environment *SuRE-pipeline*, and
additionally make sure *snakemake* is in your path.

The pipeline is run as follows:
```
CORES=10 # the maximum (absolute) number of cores that can be used concurrently
RAM=100  # the (approximate) maximum amount of RAM (in Gb) to be used by the pipeline
CONFIGFILE=SuRE-config.yaml # snakemake config file

snakemake  --cores ${CORES} --resources ram=${RAM} -Tprs path-to-pipelineRootDir/SuRE-snakemake --configfile ${CONFIGFILE} all
```

This will run the pipeline with a specified number of cores, a specified max
amount of RAM to be used, and a particular config file. Finally, the rule
*'all'* is triggered.


## Notes

The amount of memory consumed can be limited only to a proximate extend. Expect
that in fact up to 150% of the configured maximum may be used.

# Configuration of the snakemake run

See config file in repos.

