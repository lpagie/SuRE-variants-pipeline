#!/bin/bash

# AUTHOR / DATE
#   Ludo Pagie; March 22, 2016; iPCR-map.bash

# INTRO / BACKGROUND
#   bash script (awk, bowtie2, samtools and cutadapt) to process raw fastq
#   files containing data from iPCR samples. The barcodes and gDNA are
#   extracted from the reads and the gDNA is aligned to the reference genome. 
#   Barcodes with length != 20, or which contain N's are discarded, as are
#   reads with a MAPQ less then 20. The aligned and filtered paired reads are
#   written to stdout in bedpelike format including the barcode and sorted
#   (alphabetically) on barcode.
#
# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#   1 pair of fastq files for input are given as final arguments 
#   -o: output directory
#   optional:
#   -l: write to logfile instead of stdout
#   -t: statsfile [default: based on basename, in output directory]
#   -m: max insert length
#   -b: basename [based on input file name]
#   -c: if set; do not clean intermediate files
#   -n: number of cores used in parallel processes (10)
#   -s: name of bowtie2 index file
# INPUT:
#   trimmed fastq files
# OUTPUT:
#   bam files with aligned reads

# VERSIONS:
#   -170206: initial version, VERSION set to 0.0.1

# TODO
#   - parameterize filter criteria (min MAPQ score, BC length, etc)

SCRIPTNAME=iPCR-map.bash
VERSION=0.0.1 # YYMMDD

# EXTERNAL SOFTWARE
GAWK=/usr/bin/gawk
BOWTIE2=bowtie2
SAMTOOLS=$HOME/vanSteensel/bin/samtools

# GLOBAL VARIABLES
NCORES=10
MAX_INSERT_LENGTH=1000
export BOWTIE2_INDEXES=$HOME/data/bowtie2-indexes/
BOWTIE2_REFSEQ="hg19_ch1-22_XYM"
CLEAN=true;
LOG="false"

# PARSE OPTIONS
OPTIND=1         # Reset in case getopts has been used previously in the shell.
USAGE=
usage() {
  echo >&2 "usage: ${SCRIPTNAME} -o[lbncsh] forw-reads.fastq[.gz/bz2] rev-reads.fastq[.gz/.bz2]"
  echo >&2 "OPTIONS:"
  echo >&2 "  -o: directory for generated count-table files  [required]"
  echo >&2 "  -l: write messages to logfile (OUTDIR/BASENAME.log) instead of stdout"
  echo >&2 "  -b: sets basename used for all output files [default: based on input filename]"
  echo >&2 "  -n: number of cores used where possible [default: 10]"
  echo >&2 "  -c: do not clean up intermediate files [default: clean]"
  echo >&2 "  -m: max insert length for aligned read pair to be considered 'concordant' by bowtie [default: 1000]"
  echo >&2 "  -s: name of bowtie2 index file with genome reference sequence [default: hg19_ch1-22_XYM]"
  echo >&2 "  -h: print this message"
  echo >&2 ""
  exit 1;
}

while getopts "h?f:o:lb:n:m:s:t:c" opt; do
  case $opt in
    l)
      LOG="true";
      ;;
    t)
      STATS=$OPTARG;
      ;;
    n)
      NCORES=$OPTARG;
      ;;
    o)
      OUTDIR=$OPTARG;
      ;;
    b)
      BASENAME=$OPTARG;
      ;;
    c)
      CLEAN=false;
      ;;
    m)
      MAX_INSERT_LENGTH=$OPTARG;
      ;;
    s)
      BOWTIE2_REFSEQ=$OPTARG;
      ;;
    h)
      usage;
      ;;
    \?)
      echo "option not recognized: "$opt
      usage
      ;;
  esac
done
shift $(( OPTIND - 1 ))

# the remaining CLI arguments should be a pair of filenames which are the
# trimmed forward and reverse reads fastq files
# check we have exactly 2 remaining arguments
if [ ! $# -eq 2 ]; then
  echo -e "\nerror: too few, or too many, arguments left after options are parsed (should be 1 pair of fastq filenames).\nThe remaining args are:"
  while test $# -gt 0; do
    echo $1
    shift
  done
  echo -e "Aborting\n\n"
  usage
fi

# retrieve input fastq files from command line
declare -a FASTQ_FNAMES=( "$@" );

# check; file exists, whether compressed (set CAT)
# file exists
abort_flag="false"
for f in ${FASTQ_FNAMES[@]}; do
  if [ ! -f ${f} ]; then
    echo -e "error; fastq file (${f}) doesn't exist.\n" 
    abort_flag=true
  fi
done
if [ $abort_flag == 'true' ]; then
  echo -e "Aborting\n\n"
  usage
  exit 1
fi
unset abort_flag
# determine file extension of the first file (assuming both files are
# compressed identically)
f=${FASTQ_FNAMES[0]}
extension="${f##*.}"
# determine CAT
case ${extension} in
  gz)
    CAT="gzip -cd ";
    ;;
  bz2)
    CAT="bzip2 -cd ";
    ;;
  *)
    CAT="cat ";
    ;;
esac
unset extension
unset f

# Make name of fastq file absolute
for (( i=0; i<2; i++ )); do
  D=`dirname "${FASTQ_FNAMES[$i]}"`
  B=`basename "${FASTQ_FNAMES[$i]}"`
  DD="`cd $D 2>/dev/null && pwd || echo $D`"
  FASTQ_FNAMES[$i]="$DD/$B"
done

# check all required options are set
if [ -z ${OUTDIR+x} ]; then echo "option -o not set (directory for output files)"; usage; exit 1; fi
# check required subdirectories exist
if [ ! -d ${OUTDIR} ]; then mkdir -p ${OUTDIR}; echo "making directory \"${OUTDIR}\" for output"; echo ""; fi
# make path to OUTDIR absolute
OUTDIR="`cd \"$OUTDIR\" 2>/dev/null && pwd || echo \"$OUTDIR\"`"
# BASENAME
if [ -z ${BASENAME+x} ]; then
  # create BASENAME based on 1st input fastq filename remove ".fastq.*" (or ".fq.*") from filename
  BASENAME=$(basename ${FASTQ_FNAMES[0]} | sed -e 's/.[fF]\(ast\|AST\)\?[qQ].*//')
fi
# STATS FILE
if [ -z ${STATS} ]; then
  STATS="${OUTDIR}/${BASENAME}.stats"
fi

######################################
# write stdout to stdout or a log file
######################################
if [ ${LOG} == "true" ]; then 
  LOG="${OUTDIR}/${BASENAME}.log"
  exec 1>>${LOG}
fi

# print values of variables and CLI args for log
# print header for log
######################
LINE="running "${SCRIPTNAME}" (version: "$VERSION")"
SEPARATOR=$(head -c ${#LINE} </dev/zero | tr '\0' '=')
echo $SEPARATOR; 
echo $SEPARATOR
echo $LINE; 
echo $SEPARATOR
echo $SEPARATOR
echo "script context"
echo "=============="
starttime=$(date +%c)
echo "starting date/time = "${starttime}
echo "User set variables:"
echo "==================="
echo "directory for output files=${OUTDIR}"
echo "basename for output files=${BASENAME}"
echo "LOG=${LOG}"
echo "NCORES=${NCORES}"
echo "MAX_INSERT_LENGTH=${MAX_INSERT_LENGTH}"
echo "BOWTIE2_INDEXES=${BOWTIE2_INDEXES}"
echo "BOWTIE2_REFSEQ=${BOWTIE2_REFSEQ}"
echo "CLEAN=${CLEAN}"
echo ""
echo "fastq files for input:"
echo "================================="
for f in ${FASTQ_FNAMES[@]}; do echo $f; done
echo ""
# print some software version info
echo "Used software:"
echo "=============="
echo "unix/host"; uname -a; echo "---------------";
echo "bash:"; bash --version 2>&1 head -3; echo "---------------";
echo "gawk:"; echo "executable used: ${GAWK}"; ${GAWK} --version; echo "---------------";
echo "bowtie2:"; echo "executable used: ${BOWTIE2}"; ${BOWTIE2} --version; echo "---------------";
echo "samtools:"; echo "executable used: ${SAMTOOLS}"; ${SAMTOOLS} 2>&1 | head -3; echo "---------------";
echo "=============="
echo ""

echo -e "finished prepping for processing"
echo -e "================================\n\n"

echo "===================================="
echo "===================================="
echo "MAIN: starting to process fastq file" 
echo "===================================="
echo "===================================="
echo ""

#################################
#######  MAIN  ##################
#################################

# setwd processing directory
cd ${OUTDIR}

### ALIGNMENT OF PAIRED_END READS, plus filtering on concordant reads and sorting on readID ######
##################################################################################################
echo "starting alignment"
FORW="${FASTQ_FNAMES[0]}"
REV="${FASTQ_FNAMES[1]}"
BAM="${OUTDIR}/${BASENAME}.bam"

CMD="(${BOWTIE2} -p ${NCORES} -x ${BOWTIE2_REFSEQ} -1 $FORW -2 $REV -X ${MAX_INSERT_LENGTH} | \
  ${SAMTOOLS} view -b -f2 -u - -o - | \
  ${SAMTOOLS} sort -n - -o ${BAM} -@ ${NCORES} -m 2G -T ${BAM%.bam}_srt ) 2>> ${STATS}"

echo "command to run bowtie = ${CMD}"
eval "${CMD}"
# record some stats in file STATS
nreads=$(cat ${STATS} | gawk '/^[[:digit:]]+.*reads; of these:/{print $1}')
echo -e "\nalignedReadCount\t${nreads}\n\n" >> ${STATS}
echo -e "alignment done\n"

##############################
########## DONE ##############
##############################
LINE="finished "${SCRIPTNAME}
SEPARATOR=$(head -c ${#LINE} </dev/zero | tr '\0' '=')
echo $SEPARATOR; 
echo $SEPARATOR; 
echo $LINE; 
echo $SEPARATOR
echo $SEPARATOR; 
endtime=$(date +%c)
echo "end date/time = "${endtime}
echo $SEPARATOR; 
echo $SEPARATOR; 
