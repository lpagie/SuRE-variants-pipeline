#!/bin/bash

# AUTHOR / DATE
#   Ludo Pagie; March 22, 2016; cDNA-plDNA-count-BC_160322.bash

# INTRO / BACKGROUND
#   bash script (mostly awk and cutadapt) to process raw fastq files containing
#   data from cDNA/plDNA samples. The barcodes are extracted from the reads and
#   counted. For exatrcting the barcodes the adapter sequence is aligned with
#   the read (using cutadapt) and the preceding part of the read is defined as
#   the barcode.
#   Barcodes with length != 20, or which contain N's are discarded. The
#   filtered barcodes and counts are written to stdout, sorted (alphabetically)
#   on barcode.
#
# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#   -o: output directory
#   -a: adapter sequence
#   optional:
#   -b: basename [based on input file name]
#   -l: write to logfile instead of stdout
#   fastq files are expected as last arguments on the commandline
# INPUT:
#   cDNA/plDNA fastq files
# OUTPUT:
#   tabular txt file with count and barcode-sequence

# VERSIONS:
#   -160322: initial version, VERSION set to 0.0.1
#   -161129: VERSION set to 0.0.2; changed from multi fastq input to single fastq input


SCRIPTNAME=cDNA-plDNA-count-BC_160322.bash
VERSION=0.0.2 # YYMMDD

# EXTERNAL SOFTWARE
# GAWK=/usr/bin/gawk
GAWK=gawk
# CUTADAPT=/home/NFS/users/l.pagie/vanSteensel/src/cutadapt-stable_1.2.1/bin//cutadapt
# CUTADAPT=cutadapt
CUTADAPT=/home/NFS/users/l.pagie/python_virt_env_cutadapt/bin/cutadapt

# GLOBAL VARIABLES
LOG="false"

# PARSE OPTIONS
OPTIND=1         # Reset in case getopts has been used previously in the shell.
USAGE=
usage() {
  echo >&2 "usage: ${SCRIPTNAME} -oa[lbh] fastq-filename"
  echo >&2 "OPTIONS:"
  echo >&2 "  -o: directory for generated count-table files  [required]"
  echo >&2 "  -a: adapter sequence to parse read sequences [required]"
  echo >&2 "  -b: sets basename used for all output files [default: based on input filename]"
  echo >&2 "  -l: write messages to logfile (OUTDIR/BASENAME.log) instead of stdout"
  echo >&2 "  fastq-filename: final argument; single fastq filename (ending with fastq|fq, and optionally compressed with gz/bz2)"
  echo >&2 "  -h: print this message"
  echo >&2 ""
  exit 1;
}
while getopts "h?b:o:a:ln:" opt; do
  case $opt in
    l)
      LOG="true";
      ;;
    o)
      OUTDIR=$OPTARG;
      ;;
    a)
      ADPTR_SEQ=$OPTARG;
      ;;
    b)
      BASENAME=$OPTARG;
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

# CONFIGURE INPUT FASTQ FILES
# check we have exactly 1 remaining argument
if [ ! $# -eq 1 ]; then
  echo -e "\nerror: no, or too many, arguments left after options are parsed (should be a single fastq filename):"
  while test $# -gt 0; do
    echo $1
    shift
  done
  echo -e "Aborting\n\n"
  usage
fi

# retrieve input fastq files from command line
FASTQ_FNAME=( "$1" );
# check; file exists, whether compressed (set CAT)
# file exists
if [ ! -f ${FASTQ_FNAME} ]; then
  echo -e "error; fastq file (${FASTQ_FNAME}) doesn't exist.\nAborting\n\n" 
  exit 1
fi
# determine file extension
extension="${FASTQ_FNAME##*.}"
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

# Make name of fastq file absolute
D=`dirname "${FASTQ_FNAME}"`
B=`basename "${FASTQ_FNAME}"`
FASTQ_FNAME="`cd \"$D\" 2>/dev/null && pwd || echo \"$D\"`/$B"

# CHECK VARIABLES/OPTIONS
# check all required options are set
if [ -z ${OUTDIR+x} ]; then echo "option -o not set (directory for output files)"; usage; exit 1; fi
if [ -z ${ADPTR_SEQ+x} ]; then echo "option -a not set (adapter sequence for parsing reads)"; usage; exit 1; fi
# check required subdirectories exist
if [ ! -d ${OUTDIR} ]; then mkdir -p ${OUTDIR}; echo "making directory \"${OUTDIR}\" for output"; fi
# make path to OUTDIR absolute
OUTDIR="`cd \"$OUTDIR\" 2>/dev/null && pwd || echo \"$OUTDIR\"`"
# check, or create if necessary, BASENAME
if [ -z ${BASENAME+x} ]; then 
  # create BASENAME based on input fastq filename remove ".fastq.*" (or ".fq.*") from filename
  BASENAME=$(basename $FASTQ_FNAME | sed -e 's/.[fF]\(ast\|AST\)\?[qQ].*//')
fi

# LOG FUNCTION
# define function log which writes (status lines) to stderr and (if logfile is given) to LOG
if [ ${LOG} == "true" ]; then 
  LOG="${OUTDIR}/${BASENAME}.log"
  exec 1>>${LOG}
fi

# PRINT HEADER FOR LOG
# print values of variables and CLI args for log
LINE="running "${SCRIPTNAME}" (version: "$VERSION")"
SEPARATOR=$(head -c ${#LINE} </dev/zero | tr '\0' '=')
echo $SEPARATOR; echo $LINE; echo $SEPARATOR
echo "script context"
echo "=============="
starttime=$(date +%c)
echo "starting date/time = "${starttime}
echo "starting work directory = "`pwd`
echo "User set variables:"
echo "==================="
echo "directory for output files=${OUTDIR}"
echo "basename for output files=${BASENAME}"
echo "adapter sequence=${ADPTR_SEQ}"
echo "LOG=${LOG}"
echo ""
echo "fastq file for input:"
echo "================================="
echo ${FASTQ_FNAME}
echo ""
# print some software version info
echo "Used software:"
echo "=============="
echo "unix/host"; uname -a; echo "";
echo "bash:"; bash --version 2>&1 head -3; echo "";
echo "gawk:"; echo "executable used: ${GAWK}"; ${GAWK} --version; echo "";
echo "cutadapt:"; echo "executable used: ${CUTADAPT}"; ${CUTADAPT} --version; echo "";
echo "=============="
echo ""


##################################################
#########  MAIN  #################################
##################################################
# set some filenames for intermediate output
INFO=${OUTDIR}/${BASENAME}_info
OUT=${OUTDIR}/${BASENAME}_out
STATS=${OUTDIR}/${BASENAME}_trimmed.stats
TABLE=${OUTDIR}/${BASENAME}_trimmed_table.txt.gz
BCLENGTHFNAME=${OUTDIR}/${BASENAME}_BClengths.txt
NNNFNAME=${OUTDIR}/${BASENAME}_NNN.txt
# parse the reads
${CAT} ${FASTQ_FNAME} | "${CUTADAPT}" -g "${ADPTR_SEQ}" -o ${OUT} --info-file=${INFO} - > ${STATS} 

# count and filter the parsing results for barcode sequences of incorrect
# length or containing Ns. The barcode sequence is retrieved from the info
# files generated by cutadapt. The barcode is the 3rd element in the read, ie
# the 5th column in the info file.

cat ${INFO} | cut -f 5 | sort -S 50% | uniq -c | \
  ${GAWK} -v NNNfname=${NNNFNAME} -v BClengthfname=${BCLENGTHFNAME} '
BEGIN { 
OFS = "\t"; 
  } 
  { 
    l = length($2) # length of the barcode sequence
    N = $2 ~ /N/   # boolean specifying whether barcode sequence contains Ns
    if(l==20 && N==0) print $1, $2 # only print correct barcodes
      len[l]++  # collect barcode length statistics 
      NNN[N]++  # collect barcode containing Ns statistics
    }
  END {
  # print the collected statistics of the barcodes to file
  PROCINFO["sorted_in"] = "@ind_num_asc"
  for (i in len) 
    printf("%d\t%d\n", i, len[i]) > BClengthfname
    for (i in NNN)
      printf("%d\t%d\n", i, NNN[i]) > NNNfname
    } ' | \
      gzip -c > ${TABLE} 

rm -f ${OUT}
gzip ${INFO}

echo ""

########## DONE ##############
LINE="finished "${SCRIPTNAME}
SEPARATOR=$(head -c ${#LINE} </dev/zero | tr '\0' '=')
echo $SEPARATOR; echo $LINE; echo $SEPARATOR
endtime=$(date +%c)
echo "end date/time = "${endtime}
echo "==================================="
echo ""
