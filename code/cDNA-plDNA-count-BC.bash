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
#   -l: log-filename [stdout]
#   -n: number of cores used in parallel processes [10]
#   fastq files are expected as last arguments on the commandline
# INPUT:
#   cDNA/plDNA fastq files
# OUTPUT:
#   tabular txt file with count and barcode-sequence

# VERSIONS:
#   -160322: initial version, VERSION set to 0.0.1


SCRIPTNAME=cDNA-plDNA-count-BC_160322.bash
VERSION=0.0.1 # YYMMDD

# EXTERNAL SOFTWARE
# GAWK=/usr/bin/gawk
GAWK=gawk
# CUTADAPT=/home/NFS/users/l.pagie/vanSteensel/src/cutadapt-stable_1.2.1/bin//cutadapt
CUTADAPT=cutadapt
CUTADAPT=/home/NFS/users/l.pagie/python_virt_env_cutadapt/bin/cutadapt

# GLOBAL VARIABLES
NCORES=10

# PARSE OPTIONS
OPTIND=1         # Reset in case getopts has been used previously in the shell.
USAGE=
usage() {
  echo >&2 "usage: ${SCRIPTNAME} -foaln"
  echo >&2 "OPTIONS:"
  echo >&2 "  -o: directory for generated count-table files  [required]"
  echo >&2 "  -a: adapter sequence to parse read sequences [required]"
  echo >&2 "  -b: sets basename used for all output files [default: based on input filename]"
  echo >&2 "  -l: set name of logfile [default: stdout]"
  echo >&2 "  -n: number of cores used where possible [default: 10]"
  echo >&2 ""
  exit 1;
}
while getopts "h?f:o:a:l:n:" opt; do
  case $opt in
    l)
      LOG=$OPTARG;
      ;;
    n)
      NCORES=$OPTARG;
      ;;
    o)
      OUTDIR=$OPTARG;
      # make path to OUTDIR absolute
      OUTDIR="`cd \"$OUTDIR\" 2>/dev/null && pwd || echo \"$OUTDIR\"`"
      ;;
    a)
      ADPTR_SEQ=$OPTARG;
      ;;
    \?)
      echo "option not recognized: "$opt
      usage
      ;;
  esac
done
shift $(( OPTIND - 1 ))

# CONFIGURE INPUT FASTQ FILES
# retrieve input fastq files from command line
declare -a FASTQ_FNAMES=( "$@" );
# Make names of fastq files absolute
for ((i=0; i<${#FASTQ_FNAMES[@]}; i++));
do
  D=`dirname "${FASTQ_FNAMES[i]}"`
  B=`basename "${FASTQ_FNAMES[i]}"`
  FASTQ_FNAMES[i]="`cd \"$D\" 2>/dev/null && pwd || echo \"$D\"`/$B"
done

# CHECK VARIABLES/OPTIONS
# check all required options are set
if [ -z ${OUTDIR+x} ]; then echo "option -o not set (directory for output files)"; usage; exit 1; fi
if [ -z ${ADPTR_SEQ+x} ]; then echo "option -a not set (adapter sequence for parsing reads)"; usage; exit 1; fi
# check required subdirectories exist
if [ ! -d ${OUTDIR} ]; then mkdir -p ${OUTDIR}; echo "making directory \"${OUTDIR}\" for output"; fi

# LOG FUNCTION
# define function log which writes (status lines) to stderr and (if logfile is given) to LOG
if [ ! -z ${LOG+x} ]; then 
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
echo "adapter sequence=${ADPTR_SEQ}"
echo "LOG=${LOG}"
echo "NCORES=${NCORES}"
echo ""
echo "fastq files for input:"
echo "================================="
for f in ${FASTQ_FNAMES[@]}; do echo $f; done
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

# FUNCTION DEFINITIONS
function proc_func() {
  # function which takes a gzipped fastq file and parses it using cutadapt.
  # Parsed reads are then filtered for barcode sequences with incorrect length,
  # or which contain N's

  f=$(basename "$1")
  # set some filenames for intermediate output
  INFO=${f%.fastq*}_info
  OUT=${f%.fastq*}_out
  STATS=${f%.fastq*}_trimmed.stats
  TABLE=${f%.fastq*}_trimmed_table.txt.gz
  BCLENGTHFNAME=${f%.fastq*}_BClengths.txt
  NNNFNAME=${f%.fastq*}_NNN.txt
  # parse the reads
  gzip -d -c $1 | "${CUTADAPT}" -g "${ADPTR_SEQ}" -o ${OUT} --info-file=${INFO} - > ${STATS} 

  # count and filter the parsing results for barcode sequences of incorrect
  # length or containing Ns. The barcode sequence is retrieved from the info
  # files generated by cutadapt. The barcode is the 3rd element in the read, ie
  # the 5th column in the info file.
  cat ${INFO} | cut -f 5| sort -S 50% --parallel=${NCORES} | uniq -c | \
  ${GAWK} -v NNNfname=${NNNFNAME}   -v BClengthfname=${BCLENGTHFNAME} '
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
  parallel -j ${NCORES} gzip ::: ${INFO}
}

# export necessary variables for parallel processing
export CUTADAPT
export NCORES
export ADPTR_SEQ
export -f proc_func

# MAIN: PARALLEL PROCESSING OF READS
OWD=$PWD
cd ${OUTDIR}
echo "starting parallel processing"
echo "----------------------------"
parallel --joblog parallelJoblog -j ${NCORES} proc_func ::: ${FASTQ_FNAMES[@]}
cat parallelJoblog # echo the joblog to stdout
rm -f parallelJoblog
echo ""
echo "finished parallel processing"
echo "----------------------------"
echo ""
cd $OWD

########## DONE ##############
LINE="finished "${SCRIPTNAME}
SEPARATOR=$(head -c ${#LINE} </dev/zero | tr '\0' '=')
echo $SEPARATOR; echo $LINE; echo $SEPARATOR
endtime=$(date +%c)
echo "end date/time = "${endtime}
echo "==================================="
echo ""
