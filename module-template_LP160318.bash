#!/bin/bash

# AUTHOR / DATE

# INTRO / BACKGROUND

# USAGE / INPUT / ARGUMENTS / OUTPUT



# VERSIONS

# TODO

SCRIPTNAME=cDNA-plDNA-count-BC_160322.bash
VERSION=0.0.1 # YYMMDD

# EXTERNAL SOFTWARE

# GLOBAL VARIABLES
DATETAG=`date +"%y%m%d"`
TIMETAG=`date +"%H%M%S"`

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


## howto's ##
## run command and test exit status:
# command && echo "command succeeded" || echo "command failed"

########## DONE ##############
LINE="finished "${SCRIPTNAME}
SEPARATOR=$(head -c ${#LINE} </dev/zero | tr '\0' '=')
echo $SEPARATOR; ECHO $LINE; echo $SEPARATOR
endtime=$(date +%c)
echo "end date/time = "${endtime}
echo "==================================="
echo ""
