#!/bin/bash

# AUTHOR / DATE
#   Ludo Pagie; March 02, 2017; iPCR-single-QA.bash

# INTRO / BACKGROUND
#   bash script (awk, R) to collect various simple statistics of the iPCR data
#   processing for a single iPCR sample. The collected statistics are:
#   - total raw readcount
#   - readcount after; trimming, aligning, conversion to bedpe
#   - length distribution of extracted barcodes
#   - distribution of occurences of N's in barcode
#   - count of uniq bedpe (ie SuRE) fragments
#   The statistics will be collected/visualized in:
#   - a tabular text file (all read-/fragment-counts)
#   - barplot (read counts at various stages))
#   - barplots (distribution barcode length and N's)
#
# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#   -o: output directory
#   -b: basename [based on input file name]
#   optional:
#   -l: write to logfile instead of stdout
# INPUT:
#   basename; used for selecting input filenames and generating output file names
# OUTPUT:
#   tabular text file and 3 png's

# VERSIONS:
#   -170302: initial version, VERSION set to 0.0.1

# TODO

SCRIPTNAME=iPCR-single-QA.bash
VERSION=0.0.1 # YYMMDD

# EXTERNAL SOFTWARE
GAWK=/usr/bin/gawk
R=/usr/bin/R

# GLOBAL VARIABLES
NCORES=1
LOG="false"

# PARSE OPTIONS
OPTIND=1         # Reset in case getopts has been used previously in the shell.
USAGE=
usage() {
  echo >&2 "usage: ${SCRIPTNAME} -ob[l]"
  echo >&2 "OPTIONS:"
  echo >&2 "  -o: directory for generated output files  [required]"
  echo >&2 "  -b: sets basename used for all output files [required]"
  echo >&2 "  -l: write messages to logfile (OUTDIR/BASENAME.log) instead of stdout"
  echo >&2 "  -h: print this message"
  echo >&2 ""
  exit 1;
}

while getopts "h?o:b:l" opt; do
  case $opt in
    l)
      LOG="true";
      ;;
    o)
      OUTDIR=$OPTARG;
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

# check stats file containing statistics exists
STATS="${OUTDIR}/${BASENAME}.stats"
if [ ! -f ${STATS} ] ; then
  echo "File with statistics (${STATS}) not found. Aborting!"
  exit 1;
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
echo ""
# print some software version info
echo "Used software:"
echo "=============="
echo "unix/host"; uname -a; echo "---------------";
echo "bash:"; bash --version 2>&1 head -3; echo "---------------";
echo "gawk:"; echo "executable used: ${GAWK}"; ${GAWK} --version; echo "---------------";
echo "R:"; echo "executable used: ${R}"; ${R} --version; echo "---------------";
echo "=============="
echo ""

# setwd processing directory
cd ${OUTDIR}

echo -e "finished prepping for processing"
echo -e "================================\n"

echo "===================================="
echo "===================================="
echo "MAIN: starting to process iPCR stats" 
echo "===================================="
echo "===================================="
echo ""

#################################
#######  MAIN  ##################
#################################

R --no-save -q << EOR
basename <- "${BASENAME}"; 
odir <- "${OUTDIR}"
statsfile <- paste0(basename,".stats")

statlines <- readLines(statsfile)

# statlables <- c("totalReadCount", "trimmedForwReadCount", "trimmedRevReadCount", "trimmedReadCount", "lengthFilteredReadCount", "alignedReadCount", "bedpeFragmentCount", "bedpeFragmentUniqCount")
statlables <- c("totalReadCount", "trimmedForwReadCount", "trimmedRevReadCount", "lengthFilteredReadCount", "alignedReadCount", "bedpeFragmentCount", "bedpeFragmentUniqCount")

stats <- vector('integer', length=length(statlables))
names(stats) <- statlables
for (lab in statlables) 
  stats[[lab]] <- as.integer(strsplit(split="\t", grep(paste0('^', lab), statlines, value=TRUE))[[1]][[2]])

for (lab in statlables)
  cat(lab, " \t ", stats[[lab]], "\n")

ofile <- paste0(basename, "_counts.png")
bitmap(file=ofile, type='png16m', res=144, taa=4)
opar=par(mar=c(7,4,4,1)+.1)
x <- barplot(stats/1e6, beside=TRUE, names.arg=names(stats),srt=30, xaxt='n', width=1, space=0.2, ylab="count (x1e6)", main=paste(basename,": iPCR read/SuREfragment counts"))
text(cex=1, x=x+.4, y=-max(stats/1e6)*0.025, sub("Count$","",names(stats)),xpd=NA, srt=45, pos=2)
par(opar)
dev.off()

ofile <- paste0(basename, "_BClengths.png")
bitmap(file=ofile, type='png16m', res=144, taa=4)
BClengths <- data.frame(length=as.integer(strsplit(split="\t",grep('^iPCR_BC_lengths[[:space:]]+',statlines, value=T))[[1]][-1]),
			count=as.integer(strsplit(split="\t",grep('^iPCR_BC_lengths_counts[[:space:]]+',statlines, value=T))[[1]][-1]))
BClengths <- as.matrix(BClengths[order(BClengths[[1]]),])
plot(BClengths[,1],BClengths[,2]/1e6, type='h',lwd=3, ylab='frequency (x1e6)', xlab='barcode length', col='black', 
     main=paste0(basename, ": iPCR BC lengths"))
dev.off()

ofile <- paste0(basename, "_BCNNNs.png")
bitmap(file=ofile, type='png16m', res=144, taa=4)
Ncounts <- data.frame("N#"=as.integer(strsplit(split="\t",grep('iPCR_BC_NNN[[:space:]]+',statlines, value=T))[[1]][-1]),
		      count=as.integer(strsplit(split="\t",grep('iPCR_BC_NNN_counts[[:space:]]+',statlines, value=T))[[1]][-1]))
Ncounts <- as.matrix(Ncounts[order(Ncounts[[1]]),])
plot(Ncounts[,1], Ncounts[,2]/1e6, type='h',lwd=3, ylab='frequency (x1e6)', xlab='#N\'s in BC', col='black', 
     main=paste0(basename, ": iPCR N\'s in BC"), xaxt='n')
axis(1, at=Ncounts[,1])
dev.off()

EOR

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
