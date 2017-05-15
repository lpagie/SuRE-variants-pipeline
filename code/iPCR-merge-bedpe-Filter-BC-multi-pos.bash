#!/bin/bash

# AUTHOR / DATE
#   Ludo Pagie; March 19, 2016; iPCR-merge-bedpe-Filter-BC-multi-pos.bash

# INTRO / BACKGROUND
#   bash script (mostly awk) to merge a (set of) raw iPCR bedpe file and
#   remove BC-position pairs if the BC is associated with multiple positions.
#   The BC-position pair with highest iPCR count is retained and remaining
#   pairs are removed.
#
# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#     -o: output directory
#   optional:
#     -b: bedpe output file [default: OUTDIR/iPCR-combined-bedpe.txt.gz]
#     -l: write to logfile instead of stdout
#     -n number of cores to use (default 10)
#   The iPCR bedpe filenames are given as last argument to the script
# INPUT:
#   raw iPCR bedpe file(s)
# OUTPUT:
#   tabular txt file with position, strand, alignment info, gDNA sequences, barcode.
#   This output file is filtered for barcodes which occur at multiple
#   positions. But positions which associate with multiple barcodes are still
#   included.

# VERSIONS:
#   -160319: initial version, VERSION set to 0.0.1
#   -161213: adapted for use with snakemake, version 0.0.2

# TODO
# - check supplied bedpe files exist etc

VERSION=0.0.2 # YYMMDD
SCRIPTNAME=iPCR-merge-bedpe-Filter-BC-multi-pos.bash

# EXTERNAL SOFTWARE
GAWK=gawk

# GLOBAL VARIABLES
NCORES=10
LOG="false"

###############
# PARSE OPTIONS
###############
OPTIND=1         # Reset in case getopts has been used previously in the shell.
USAGE=
usage() {
  echo >&2 "usage: ${SCRIPTNAME} -lbno bedpe1 bedpe2 bedpe3 ..."
  echo >&2 "OPTIONS:"
  echo >&2 "  -o: name of output directory [required]"
  echo >2& "  -b: bedpe output file [default: OUTDIR/iPCR-combined-bedpe.txt.gz]"
  echo >&2 "  -l: write messages to logfile (OUTDIR/iPCR-merge.log) instead of stdout"
  echo >&2 "  -n: number of cores used where possible [10]"
  echo >&2 ""
  exit 1;
}

while getopts "h?ln:o:b:" opt; do
  case $opt in
    l)
      LOG="true";
      ;;
    n)
      NCORES=$OPTARG;
      ;;
    o)
      OUTDIR=$OPTARG;
      ;;
    b)
      OUTPUT=$OPTARG;
      ;;
    \?)
      usage
      ;;
  esac
done
shift $(( OPTIND - 1 ))

############################################
# set remaining arguments as BEDPE filenames
############################################
BEDPE_FNAME=$*
# check; file exists, whether compressed (set CAT)
# file exists
abort_flag="false"
for f in ${BEDPE_FNAME}; do
  if [ ! -f ${f} ]; then
    echo -e "error; bedpe file (${f}) doesn't exist.\n" 
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
f=${BEDPE_FNAME[0]}
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
# Make name of bedpe file absolute
for (( i=0; i<$#; i++ )); do
  D=`dirname "${BEDPE_FNAME[$i]}"`
  B=`basename "${BEDPE_FNAME[$i]}"`
  DD="`cd $D 2>/dev/null && pwd || echo $D`"
  BEDPE_FNAME[$i]="$DD/$B"
done

####################################
# check all required options are set
####################################
if [ -z ${OUTDIR+x} ]; 
then 
  echo "option -o not set (directory for output files)"; usage; 
  exit 1; 
fi
# check required subdirectories exist
if [ ! -d ${OUTDIR} ]; then mkdir -p ${OUTDIR}; echo "making directory \"${OUTDIR}\" for output"; fi
# make path to OUTDIR absolute
OUTDIR="`cd \"$OUTDIR\" 2>/dev/null && pwd || echo \"$OUTDIR\"`"
# set OUTPUT file
if [ -z ${OUTDIR+x} ]; then
  OUTPUT="${OUTDIR}/iPCR-combined-bedpe.txt.gz"
fi
# make path to OUTPUT absolute
D=`dirname "${OUTPUT}"`
B=`basename "${OUTPUT}"`
DD="`cd $D 2>/dev/null && pwd || echo $D`"
OUTPUT="$DD/$B"

######################################
# write stdout to stdout or a log file
######################################
if [ ${LOG} == "true" ]; then 
  LOG="${OUTDIR}/iPCR-merge.log"
  exec 1>>${LOG}
fi

################################################
# print values of variables and CLI args for log
################################################
LINE="running "${SCRIPTNAME}" (version: "$VERSION")"
SEPARATOR=$(head -c ${#LINE} </dev/zero | tr '\0' '=')
echo $SEPARATOR; echo $LINE; echo $SEPARATOR
echo "script context"
echo "=============="
starttime=$(date +%c)
echo "starting date/time = "${starttime}
echo "User set variables:"
echo "==================="
echo "LOG=${LOG}"
echo "NCORES=${NCORES}"
echo "directory for output files=${OUTDIR}"
echo "main output file=${OUTPUT}"
echo ""
echo "iPCR bedpe input files";
echo "======================";
for f in $BEDPE_FNAME; do echo $f; done
echo ""
# print some software version info
echo "Used software:"
echo "=============="
echo "unix/host"; echo "---------"; uname -a; echo "";
echo "bash"; echo "----"; bash --version 2>&1 head -3; echo "";
echo "gawk"; echo "----"; echo "executable used: ${GAWK}"; ${GAWK} --version; echo "";
echo "=============="
echo ""


##################################################
#########  MAIN  #################################
##################################################
# pipe all bedpe files into ...
${CAT} $BEDPE_FNAME | \
  # sort bedpe data on barcode (dictionary order)
  sort -S 50% --parallel=${NCORES} -k 6,6d | \
  ${GAWK} -v BC_MULTI_POS_FNAME="${OUTDIR}/BC_multi_pos_${DATETAG}_${TIMETAG}.txt" '
  # this piece of awk collapses iPCR bedpe data "points": 
  # while consecutive lines contain identical barcodes the lines are collected in a array "PREVLINE"
  # once a line with a different barcode is read the "PREVLINE" array is
  #   processed before resetting the array with the new line as only element
  # processing of the array means that all bedpe data "points" are grouped on
  #   genomic position; the position with max readcount is chosen as "true"
  #   barcode; this barcode is printed with corrected values for read-count and
  #   mapq
  BEGIN {
    FS="\t"
    OFS="\t"

    # set some column numbers
    COLLENGTH=4
    COLMAPQ=10
    COLBC=6
    COLCHR=1
    COLS = 2
    COLE = 3
    COLSTR = 5
    COLCNT = 7
    COLSEQ1 = 15
    COLSEQ2 = 16
  }
  NR == 1 {
    PREVBC      = $COLBC
    PREVLINE[1] = $0
    next
  }
  {
    # here the fragment is considered good quality (ie properly alligned and such).
    # the fragments are sorted on barcode. If consecutive barcodes are
    # identical all but the most frequent one is discarded
    if ($COLBC == PREVBC) {
      # we read an identical barcode; add this line to the line stack
      PREVLINE[length(PREVLINE)+1] = $0
    } else {
      # we read a novel barcode; 
      # process the previous barcode
      if (length(PREVLINE) == 1) {
	      # if previous barcode is present at a single position just print the corresponding line
    	  print PREVLINE[1]
      } else {
      	# if the barcode exists at potentially multiple positions; clean the
        # set (BC2POS) and then print the corresponding line
      	delete BC2POS
      	for (line in PREVLINE) {
      	  # loop over all previous line, extract position, compare, select best, and print selection
      	  n = split(PREVLINE[line], records)
      	  POS = records[COLCHR]"\t"records[COLS]"\t"records[COLE]"\t"records[COLLENGTH]"\t"records[COLSTR]
      	  if (POS in BC2POS)
      	    BC2POS[POS] += records[COLCNT]
      	  else
      	    BC2POS[POS] = records[COLCNT]
      	  # print each case to a log file for statistics
      	  print PREVLINE[line] > BC_MULTI_POS_FNAME
      	}
      	# select position with highest count
      	max = 0
      	for (pos in BC2POS) {
      	  if (BC2POS[pos] > max) {
      	    # update best element
      	    max = BC2POS[pos]
      	    maxpos = pos
      	  }
      	}
      	# print position maxpos, with sum of counts of that position, and max mapq of that position
      	mapq = 0
      	cnt = 0
      	readseq = ""
      	for (line in PREVLINE) {
      	  if (match(PREVLINE[line], maxpos)) {
      	    n = split(PREVLINE[line], w)
      	    cnt += w[COLCNT]
      	    if (w[COLMAPQ] > mapq) {
      	      mapq = w[COLMAPQ]
      	      readseq = w[COLSEQ1]"\t"w[COLSEQ2]
      	    }
      	  }
      	}
        # update last line at position maxpos with new values for count/snp/etc
      	w[COLCNT] = cnt
      	w[COLMAPQ] = mapq
      	split(readseq, s, "\t")
      	w[COLSEQ1] = s[1]
      	w[COLSEQ2] = s[2]
      	# print a line composed of w
      	l=w[1]
      	for (i=2; i<=length(w); i++)
      	  l=l "\t" w[i]
      	print l
        # print final line to logfile as well, prepended with "#" to indicate which line is selected
      	print "#"l > BC_MULTI_POS_FNAME
      }
      # previous barcode is printed. start a new line stack
      delete PREVLINE
      PREVLINE[1] = $0
      PREVBC = $COLBC
    }
  } ' | \
gzip -c > ${OUTPUT}


LINE="finished "${SCRIPTNAME}
SEPARATOR=$(head -c ${#LINE} </dev/zero | tr '\0' '=')
echo $SEPARATOR; echo $LINE; echo $SEPARATOR
endtime=$(date +%c)
echo "end date/time = "${endtime}
echo "==================================="
echo ""
