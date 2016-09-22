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
#     -o name of output file
#   optional:
#     -l a logfilename (default is to print to stderr only)
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

# TODO
# - check supplied bedpe files exist etc

VERSION=0.0.1 # YYMMDD
SCRIPTNAME=iPCR-merge-bedpe-Filter-BC-multi-pos.bash

# EXTERNAL SOFTWARE
GAWK=/usr/bin/gawk

# GLOBAL VARIABLES
NCORES=10

# PARSE OPTIONS
OPTIND=1         # Reset in case getopts has been used previously in the shell.
USAGE=
usage() {
  echo >&2 "usage: ${SCRIPTNAME} -lno bedpe-file-names"
  echo >&2 "OPTIONS:"
  echo >&2 "  -o: name of outputfile [required]"
  echo >&2 "  -l: set name of logfile (in all cases log is printed to stderr)"
  echo >&2 "  -n: number of cores used where possible [10]"
  echo >&2 ""
  exit 1;
}

while getopts "h?l:n:o:" opt; do
  case $opt in
    l)
      LOG=$OPTARG;
      ;;
    n)
      NCORES=$OPTARG;
      ;;
    o)
      OUTPUT=$OPTARG;
      ;;
    \?)
      usage
      ;;
  esac
done
shift $(( OPTIND - 1 ))

BEDPE_FNAME=$*

# check all required options are set
if [ -z ${OUTPUT+x} ]; 
then 
  echo "option -o not set (output filename)"; usage;
  exit 1; 
fi

# define function log which writes (status lines) to stderr and (if logfile is given) to LOG
if [ ! -z ${LOG+x} ]; then 
  exec 1>>${LOG}
fi


# print values of variables and CLI args for log
# print header for log
######################
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
echo "OUTPUT=${OUTPUT}"
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


# check required subdirectories exist
# if [ ! -d iPCR ]; then mkdir -p "iPCR"; echo "making directory \"iPCR\" for output"; fi
# if [ ! -d output ]; then mkdir -p "output"; echo "making directory \"output\" for general output"; fi
# check/make directory which will contain output file
OUTDIR=`dirname ${OUTPUT}`
echo "OUTDIR = ${OUTDIR}"
if [ ! -d ${OUTDIR} ]; 
then
  mkdir -p ${OUTDIR}
  echo "made OUTDIR ${OUTDIR}"
fi

gunzip -c $BEDPE_FNAME | \
#   # tee ${TMP}/interOutput_0.txt | \
  sort -S 50% --parallel=${NCORES} -k 6,6d | \
#   ${TMP}/interOutput_0.txt
#  tee ${TMP}/interOutput_1.txt | \
# 
# exit 0;
#cat ${TMP}/interOutput_0.txt | \
  ${GAWK} -v BC_MULTI_POS_FNAME="${OUTDIR}/BC_multi_pos_${DATETAG}_${TIMETAG}.txt" '
  BEGIN {
    FS="\t"
    OFS="\t"

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

    # BC_MULTI_POS_FNAME = "iPCR/BC_multi_pos.txt"
  }
  NR == 1 {
    PREVBC      = $COLBC
    PREVLINE[1] = $0
  }
  NR > 1 {
    # here the fragment is considered good quality.
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
	# set and then print the corresponding line
	delete BC2POS
	for (line in PREVLINE) {
	  # loop over all previous line, extract position, compare, select best, and print selection
	  n = split(PREVLINE[line], records)
	  POS = records[1]"\t"records[2]"\t"records[3]"\t"records[4]"\t"records[5]
	  if (POS in BC2POS)
	    BC2POS[POS] += records[COLCNT]
	  else
	    BC2POS[POS] = records[COLCNT]
	  # print each case for statistics
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
