#!/bin/bash

# AUTHOR / DATE
#   Ludo Pagie; April 08, 2016; counts2coverage.bash

# INTRO / BACKGROUND

# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#     -o name of output file
#     -i name of input file
#   optional:
#     -l a logfilename (default is to print to stderr only)
#     -n number of cores to use (default 10)
# INPUT:
#   tabular data file with read counts for all SuRE-fragments, for all iPCR, plDNA, cDNA samples
# OUTPUT:
#   big wig files for various combinations of samples:
#   - for each bio replicate of all samples
#   - for sum of normalised bio reps for all samples
#   - normalized SuRE: ratio of summed cDNA, and summed (plDNA+1) samples

# VERSIONS:
#   -160408: initial version, VERSION set to 0.0.1



# VERSIONS

# TODO

SCRIPTNAME=counts2coverage.bash
VERSION=0.0.1 # 160408

# EXTERNAL SOFTWARE

# GLOBAL VARIABLES
#DATETAG=`date +"%y%m%d"`
#TIMETAG=`date +"%H%M%S"`

# PARSE OPTIONS
OPTIND=1         # Reset in case getopts has been used previously in the shell.
USAGE=
usage() {
  echo >&2 "usage: ${SCRIPTNAME} -lch"
  echo >&2 "OPTIONS:"
  echo >&2 "  -l: set name of logfile (in all cases log is printed to stderr)"
  echo >&2 "  -c: number of cores used for Bowtie2 [32]"
#  echo >&2 "  -k: keep intermediate files [false]"
#  echo >&2 "  -f: forward_adapter_sequence [$ADPTR_FORW_SEQ]"
#  echo >&2 "  -r: reverse_adapter_sequence [$ADPTR_REV_SEQ]"
#  echo >&2 "  -d: restriction enzyme site [$RESTRICT_SITE]"
#  echo >&2 "  -o: name of output directory [pipeline_output_SIGTAG]"
#  echo >&2 "  -g: name of BOWTIE2 index file for genome reference sequence [hg19_ch1-22_XYM]"
  echo >&2 ""
  exit 1;
}

while getopts "h?kf:r:c:o:g:d:l:" opt; do
  case $opt in
    l)
      LOG=$OPTARG;
      ;;
    c)
      NCORES=$OPTARG;
      ;;
#    k)
#      CLEAN=false;
#      ;;
#    f)
#      ADPTR_FORW_SEQ=$OPTARG;
#      ;;
#    r)
#      ADPTR_REV_SEQ=$OPTARG;
#      ;;
#    d)
#      RESTRICT_SITE=$OPTARG;
#      ;;
#    o)
#      OUTDIR=$OPTARG;
#      ;;
#    g)
#      BOWTIE2_REFSEQ=$OPTARG;
#      ;;
    \?)
      usage
      ;;
  esac
done
shift $(( OPTIND - 1 ))

# CHECK ALL REQUIRED OPTIONS ARE SET BY USER
if [ -z ${LOG+x} ]; then echo "option log (-l) not set, aborting"; exit 1; fi

# define function log which writes (status lines) to stderr and (if logfile is given) to LOG
#  if [ -z ${LOG+x} ]; then LOGREDIRECT=""; else LOGREDIRECT=" | tee ${LOG}"; fi
#  function log {
#    msg=$1
#    >&2 eval "echo " ${msg} ${LOGREDIRECT}
#  }# define function log which writes (status lines) to stderr and (if logfile is given) to LOG
if [ ! -z ${LOG+x} ]; then
   exec 1>>${LOG}
fi


# print values of variables and CLI args for log
# print header for log
######################
LINE="running "${SCRIPTNAME}" (version: "$VERSION")"
SEPARATOR=$(head -c ${#LINE} </dev/zero | tr '\0' '=')
echo $SEPARATOR; ECHO $LINE; echo $SEPARATOR
echo "User set variables:"
echo "==================="
echo "script context"
echo "=============="
starttime=$(date +%c)
echo "starting date/time = "${starttime}
echo "User set variables:"
echo "==================="
echo "LOG=${LOG}"
echo "NCORES=${NCORES}"
#echo "CLEAN=${CLEAN}"
#echo "ADPTR_FORW_SEQ=${ADPTR_FORW_SEQ}"
#echo "ADPTR_REV_SEQ=${ADPTR_REV_SEQ}"
#echo "RESTRICT_SITE=${RESTRICT_SITE}"
#echo "OUTDIR=${OUTDIR}"
#echo "BOWTIE REFSEQ INDEX=${BOWTIE2_REFSEQ}"
#echo ""
#echo "FASTQ input files:"
#printf -- '%s\n' "${INPUTFILES[@]}"
echo ""
# print some software version info
echo "Used software:"
echo "=============="
echo "unix/host"
uname -a
echo ""
echo "bash:"
bash --version 2>&1 head -3
echo ""
#echo "bowtie2:"
#echo "executable used: ${BOWTIE2}"
#${BOWTIE2} --version
#echo ""
#echo "bedtools:"
#echo "executable used: ${BEDTOOLS}"
#${BEDTOOLS} -version
#echo ""
#echo "custom perl script:"
#echo "executable used: ${REPAIR_PAIRS_PERL}"
#echo ""
#echo "cutadapt:"
#echo "executable used: ${CUTADAPT}"
#${CUTADAPT} --version
#echo ""
#echo "python:"
#echo "executable used: python"
#((python --version) 2>&1)
#echo ""
#echo "awk:"
#echo "executable used: ${AWK}"
#${AWK} --version
#echo ""
echo "=============="
echo ""



# read sample file; this describes for each samples: tech repl #, bio repl #, type (iPCR, cDNA, plDNA), sample short name

# create bigwig coverage files for all samples

# exp 42
EXP="K562_42"
BGOUT="${EXP}.bedgraph"
WOUT="${EXP}.wig"
BWOUT="${EXP}.bw"
zcat SuRE_pipeline_OUTPUT_160413_154017/SuRE-counts.txt_K562.gz | \
  head -10000 | \
gawk ' 
BEGIN {
  FS="\t"; 
  OFS="\t"
} 
{ 
  tot = $6+$7+$8+$9+$10+$19
  for (i=1; i<=tot; i++)
    print $1, $2, $3
} ' | \
bedtools genomecov -d -i stdin -g hg19.chrom.sizes | sort -S50% -k1,1 -k2,2n > ${WOUT}
bedGraphToBigWig ${BGOUT} hg19.chrom.sizes ${BWOUT}

# exp 45.B1
EXP="K562_45.B1"
BGOUT="${EXP}.bedgraph"
BWOUT="${EXP}.bw"
zcat SuRE_pipeline_OUTPUT_160413_154017/SuRE-counts.txt_K562.gz | \
gawk ' 
BEGIN {
  FS="\t"; 
  OFS="\t"
} 
{ 
  tot = $11+$12
  for (i=1; i<=tot; i++)
    print $1, $2, $3
} ' | \
bedtools genomecov -bga -i stdin -g hg19.chrom.sizes > ${BGOUT}
bedGraphToBigWig ${BGOUT} hg19.chrom.sizes ${BWOUT}

# exp 45.B2
EXP="K562_45.B2"
BGOUT="${EXP}.bedgraph"
BWOUT="${EXP}.bw"
zcat SuRE_pipeline_OUTPUT_160413_154017/SuRE-counts.txt_K562.gz | \
gawk ' 
BEGIN {
  FS="\t"; 
  OFS="\t"
} 
{ 
  tot = $13+$14
  for (i=1; i<=tot; i++)
    print $1, $2, $3
} ' | \
bedtools genomecov -bga -i stdin -g hg19.chrom.sizes > ${BGOUT}
bedGraphToBigWig ${BGOUT} hg19.chrom.sizes ${BWOUT}




########## DONE ##############
LINE="finished "${SCRIPTNAME}
SEPARATOR=$(head -c ${#LINE} </dev/zero | tr '\0' '=')
echo $SEPARATOR; ECHO $LINE; echo $SEPARATOR
endtime=$(date +%c)
echo "end date/time = "${endtime}
echo "==================================="
echo ""
