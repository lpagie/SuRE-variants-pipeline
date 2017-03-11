#!/bin/bash

# AUTHOR / DATE
#   Ludo Pagie; March 22, 2016; iPCR-trim.bash

# INTRO / BACKGROUND
#   bash script (awk, bowtie2, samtools and cutadapt) to process raw fastq
#   files containing data from iPCR samples. The barcodes and gDNA are
#   extracted from the reads Barcodes with length != 20, or which contain N's
#   are discarded
#
# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#   1 pair of fastq files for input are given as final arguments 
#   -o: output directory
#   optional:
#   -f: forward adapter sequence [CCTAGCTAACTATAACGGTCCTAAGGTAGCGAACCAGTGAT]
#   -r: reverse adapter sequence [CCAGTCGT]
#   -d: digestion site [CATG]
#   -l: write to logfile instead of stdout
#   -b: basename [based on input file name]
#   -c: if set; do not clean intermediate files
#   -n: number of cores used in parallel processes (10)
# INPUT:
#   iPCR fastq files
# OUTPUT:
#   fastq files

# VERSIONS:
#   -170203: initial version, VERSION set to 0.0.1

# TODO
#   - parameterize filter criteria (min MAPQ score, BC length, etc)

SCRIPTNAME=iPCR-trim.bash
VERSION=0.0.1 # YYMMDD

# EXTERNAL SOFTWARE
GAWK=/usr/bin/gawk
CUTADAPT=/home/NFS/users/l.pagie/python_virt_env_cutadapt/bin/cutadapt

# GLOBAL VARIABLES
NCORES=1
MIN_READ_LENGTH=5
MAX_INSERT_LENGTH=4000
ADPTR_FORW_SEQ="CCTAGCTAACTATAACGGTCCTAAGGTAGCGAACCAGTGAT"
ADPTR_REV_SEQ="CCAGTCGT"
RESTRICT_SITE="CATG"
CLEAN=true;
LOG="false"

# PARSE OPTIONS
OPTIND=1         # Reset in case getopts has been used previously in the shell.
USAGE=
usage() {
  echo >&2 "usage: ${SCRIPTNAME} -o[frdlbnch] forw-reads.fastq[.gz/bz2] rev-reads.fastq[.gz/.bz2]"
  echo >&2 "OPTIONS:"
  echo >&2 "  -o: directory for generated output files  [required]"
  echo >&2 "  -f: forward read adapter sequence [CCTAGCTAACTATAACGGTCCTAAGGTAGCGAACCAGTGAT]"
  echo >&2 "  -r: reverse read adapter sequence [CCAGTCGT]"
  echo >&2 "  -d: digestion site used for minimzing iPCR circle [CATG]"
  echo >&2 "  -l: write messages to logfile (OUTDIR/BASENAME.log) instead of stdout"
  echo >&2 "  -b: sets basename used for all output files [default: based on input filename]"
  echo >&2 "  -n: number of cores used where possible [default: 1]"
  echo >&2 "  -c: do not clean up intermediate files [default: clean]"
  echo >&2 "  -h: print this message"
  echo >&2 ""
  exit 1;
}

while getopts "h?f:o:r:d:lb:n:c" opt; do
  case $opt in
    l)
      LOG="true";
      ;;
    n)
      NCORES=$OPTARG;
      ;;
    f)
      ADPTR_FORW_SEQ=$OPTARG;
      ;;
    r)
      ADPTR_REV_SEQ=$OPTARG;
      ;;
    d)
      RESTRICT_SITE=$OPTARG;
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
# forward and reverse reads fastq files
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
if [ -z ${BASENAME+x} ]; then
  # create BASENAME based on 1st input fastq filename remove ".fastq.*" (or ".fq.*") from filename
  BASENAME=$(basename ${FASTQ_FNAMES[0]} | sed -e 's/.[fF]\(ast\|AST\)\?[qQ].*//')
fi

# check required subdirectories exist
if [ ! -d ${OUTDIR} ]; then mkdir -p ${OUTDIR}; echo "making directory \"${OUTDIR}\" for output"; echo ""; fi
# make path to OUTDIR absolute
OUTDIR="`cd \"$OUTDIR\" 2>/dev/null && pwd || echo \"$OUTDIR\"`"

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
echo "adapter sequence=${ADPTR_SEQ}"
echo "LOG=${LOG}"
echo "NCORES=${NCORES}"
echo "MIN_READ_LENGTH=${MIN_READ_LENGTH}"
echo "ADPTR_FORW_SEQ=${ADPTR_FORW_SEQ}"
echo "ADPTR_REV_SEQ=${ADPTR_REV_SEQ}"
echo "RESTRICT_SITE=${RESTRICT_SITE}"
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
echo "cutadapt:"; echo "executable used: ${CUTADAPT}"; ${CUTADAPT} --version; echo "---------------";
echo "=============="
echo ""

# setwd processing directory
cd ${OUTDIR}

echo -e "finished prepping for processing"
echo -e "================================\n"

echo "===================================="
echo "===================================="
echo "MAIN: starting to process fastq file" 
echo "===================================="
echo "===================================="
echo ""

#################################
#######  MAIN  ##################
#################################

STATS="${OUTDIR}/${BASENAME}.stats"

### READS TRIMMING ######
#########################
# construct cutadapt command:
createCMD () {
#CMD="${CUTADAPT} -g ${ADPTR} -o ${OUTDIR}/${BASENAME}_${DIR}.fastq --discard-untrimmed \
#       --info-file=${OUTDIR}/${BASENAME}_${DIR}.info -O4 ${FASTQ} >> ${OUTDIR}/${BASENAME}_${DIR}.stats"
CMD="${CUTADAPT} -g ${ADPTR} -o ${OUTDIR}/${BASENAME}_${DIR}.fastq \
       --info-file=${OUTDIR}/${BASENAME}_${DIR}.info -O4 ${FASTQ} >> ${OUTDIR}/${BASENAME}_${DIR}.stats"
if [[ ! -z ${RESTRICT_SITE} ]]
then
  echo "extending cutadapt CMD because restrict-site is not empty (${RESTRICT_SITE})"
  CMD="${CMD}; \
         mv ${OUTDIR}/${BASENAME}_${DIR}.fastq ${OUTDIR}/tmp.${BASENAME}_${DIR}.fastq; \
         ${CUTADAPT} -a ${RESTRICT_SITE} -o ${OUTDIR}/${BASENAME}_${DIR}.fastq -O4 \
         ${OUTDIR}/tmp.${BASENAME}_${DIR}.fastq >> ${OUTDIR}/${BASENAME}_${DIR}.stats; \
	 rm -f ${OUTDIR}/tmp.${BASENAME}_${DIR}.fastq"
fi
}

# trim forw read #
##################
echo "forward reads file = ${FASTQ_FNAMES[0]}"
echo "starting to trim adapter in forward reads; trim adapter from 5'" 
echo "and trim all after digest restriction site (${RESTRICT_SITE}) site on 3'"
ADPTR=${ADPTR_FORW_SEQ}
FASTQ=${FASTQ_FNAMES[0]}
DIR=forw
createCMD;
echo "cutadapt command = ${CMD}"
eval "${CMD}"
# record some stats in file STATS
cat "${OUTDIR}/${BASENAME}_${DIR}.stats" | \
  ${GAWK} '
    BEGIN {OFS="\t"}
    /Total reads processed/{
      printf ("totalReadCount\t%d\n", 
              gensub(/,/,"","g", gensub(/^.*:\s*(.*)$/,"\\1", "g"))); 
      exit}' >> ${STATS}
cat "${OUTDIR}/${BASENAME}_${DIR}.stats" | \
  ${GAWK} '
    BEGIN {OFS="\t"}
    /Reads with adapters/{
      printf ("trimmedForwReadCount\t%d\n", 
              gensub(/,/,"","g", gensub(/^.*:\s*(.*)$/,"\\1", "g"))); 
      exit}' >> ${STATS}
echo -e "finished trimming adapter in forward reads\n\n"

# trim reverse read #
#####################
echo "reverse reads file = ${FASTQ_FNAMES[1]}"
echo "starting to trim adapter in reverse reads; trim adapter from 5'" 
echo "and trim all after digest restriction site (${RESTRICT_SITE}) site on 3'"
ADPTR=${ADPTR_REV_SEQ}
FASTQ=${FASTQ_FNAMES[1]}
DIR=rev
createCMD;
echo "cutadapt command = ${CMD}"
eval "${CMD}"
# record some stats in file STATS
cat "${OUTDIR}/${BASENAME}_${DIR}.stats" | \
  ${GAWK} '
    BEGIN {OFS="\t"}
    /Reads with adapters/{
      printf ("trimmedRevReadCount\t%d\n", 
              gensub(/,/,"","g", gensub(/^.*:\s*(.*)$/,"\\1", "g"))); 
      exit}' >> ${STATS}
echo -e "finished trimming adapter in reverse reads\n\n"

##  # the trimmed fastq files may not contain the same set of reads
##  # unify the two fastq files:
##  # read forward reads, store with readID in array
##  # read reverse reads; for every read also in forward array print both reads
##  mv "${OUTDIR}/${BASENAME}_forw.fastq" "${OUTDIR}/tmp.${BASENAME}_forw.fastq"
##  mv "${OUTDIR}/${BASENAME}_rev.fastq" "${OUTDIR}/tmp.${BASENAME}_rev.fastq"
##  
##  ${GAWK} -v forw="${OUTDIR}/${BASENAME}_forw.fastq" -v rev="${OUTDIR}/${BASENAME}_rev.fastq" -v statsfile="${STATS}" ' 
##  BEGIN { OFS="\t" }
##  ## process forward reads first
##  (NR==FNR) && (NR%4==1) {
##    # set readID for current read
##    readID=$1
##  }
##  (NR==FNR) {
##    # store 4 fastq lines in array with current readID
##    # weird indexing to get the 4 fastq lines in proper order (1,2,3,4) into the array
##    reads[readID][((NR-1)%4)+1]=$0
##    next
##  }
##  ## process reverse reads
##  {
##    if ($1 in reads) {
##      # readID of reverse reads is in array with forward reads print both forward and reverse read sets
##      for (i=1; i<=4; i++) {
##        print reads[$1][i] >> forw
##      }
##      # print current reverse read to rev file
##      print $0 >> rev
##      toread=3
##      while (toread-- > 0) {
##        getline
##        print $0 >> rev
##      }
##      totalReads++
##    }
##  }
##  END {
##    # print total readcount to stats file
##    printf ("trimmedReadCount\t%d\n", totalReads) >> statsfile
##  }
##  ' ${OUTDIR}/tmp.${BASENAME}_forw.fastq ${OUTDIR}/tmp.${BASENAME}_rev.fastq
##  # delete tmp fastq files
##  # rm -f ${OUTDIR}/tmp.${BASENAME}_forw.fastq ${OUTDIR}/tmp.${BASENAME}_rev.fastq

# remove all reads which are $MIN_READ_LENGTH basepairs or shorter
##################################################################
echo "starting to filtered reads too short for aligning to genome"
# filter_read_length
FORW="${OUTDIR}/${BASENAME}_forw.fastq.tmp"
REV="${OUTDIR}/${BASENAME}_rev.fastq.tmp"
FORW_FLTR="${OUTDIR}/${BASENAME}_forw.fastq"
REV_FLTR="${OUTDIR}/${BASENAME}_rev.fastq"
INFO_FORW="${OUTDIR}/${BASENAME}_forw.info"

mv $FORW_FLTR $FORW
mv $REV_FLTR $REV

${GAWK} -v file1=${FORW}  -v file2=${REV} -v out1=${FORW_FLTR} -v out2=${REV_FLTR} -v min_length=${MIN_READ_LENGTH} '
  BEGIN {
    OFS="\n"
    FS="\t"
    incl=0;
    excl=0;
    while((getline id1 < file1)>0) { 
      getline seq1 < file1; getline p1 < file1; getline qual1 < file1; 
      getline id2 <file2; getline seq2 < file2; getline p2 < file2; getline qual2 < file2;
      if ( length(seq1) > min_length && length(seq2) > min_length ) {
        print id1,seq1,p1,qual1 > out1;
        print id2,seq2,p2,qual2 > out2;
        incl++;
      } else {
        excl++;
      }
    }
    print "removed "excl" paired-end reads for which either forward or reverse read were "min_length"bp or shorter\n";
    print incl" reads passed this filter\n\n";
    printf ("lengthFilteredReadCount\t%d\n\n", incl);
  }' >> ${STATS}
# delete temporary fastq files
rm -f *fastq.tmp
echo "finished filtered read on length"
echo ""

echo "finished processing ${BASENAME} files"

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
