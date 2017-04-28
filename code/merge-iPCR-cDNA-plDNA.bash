#!/bin/bash

# AUTHOR / DATE
#   Ludo Pagie; March 18, 2016; merge-iPCR-cDNA-plDNA.bash

# INTRO / BACKGROUND
#   bash script (mostly awk) to merge data of iPCR, cDNA, and plDNA data from a
#   SuRE experiment into a single tabular text file
#   All input files are sorted on barcode, the iPCR input is still redundant,
#   i.e. multiple barcodes may be associated with a singe position
#   (SuRE-fragment). But the iPCR data is cleaned from barcodes associated with
#   multiple positions.
#   The output is a single tabular txt file/pipe with only positional
#   information from the iPCR bedpe file and columns for each cDNA/plDNA file
#   with corresponding counts.
#   The ouput file is ordered along genomic position.
#
# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   merge-iPCR-cDNA-plDNA.bash -i iPCR-input -o output-filename -[ln] cDNA-input-files plDNA-input-files
#   required:
#     -i: name of iPCR bedpe file
#     -o: output directory
#   optional:
#     -l: log-filename [stdout]
#     -n: number of cores used in parallel processes (10)
#     -h: print usage
# INPUT:
#   iPCR bedpe file; redundant, sorted on barcode
#   cDNA/plDNA count-tables, sorted on barcode
# OUTPUT:
#   tabular txt file with position, strand, sample-counts, ordered on position

# VERSIONS:
#   -160318: initial version, VERSION set to 0.0.1
#   -161214: adapted to be used with snakemake, VERSION set to 0.0.2

# TODO
# merge counts if sample names for groups of identical sample names (ie PL)

VERSION=0.0.2 # YYMMDD
SCRIPTNAME=merge-iPCR-cDNA-plDNA.bash

# EXTERNAL SOFTWARE
GAWK=/usr/bin/gawk

# GLOBAL VARIABLES
NCORES=10
LOG="false"
# PLDNA_DIR=plDNA/count-tables/
# CDNA_DIR=cDNA/count-tables/

# PARSE OPTIONS
OPTIND=1         # Reset in case getopts has been used previously in the shell.
USAGE=
usage() {
  echo >&2 "usage: SCRIPTNAME -io[lnh] cDNA-input-files plDNA-input-files"
  echo >&2 "OPTIONS:"
  echo >&2 "  -o: directory for generated count-table files  [required]"
  echo >&2 "  -i: iPCR bedpe filename [required]"
  echo >&2 "  -l: set name of logfile [default: stdout]"
  echo >&2 "  -n: number of cores used where possible [default: 10]"
  echo >&2 "  -h: print this message"
  echo >&2 ""
  exit 1;
}

while getopts "h?s:o:i:ln:" opt; do
  case $opt in
    l)
      LOG=true;
      ;;
    n)
      NCORES=$OPTARG;
      ;;
    i)
      IPCR_FNAME=$OPTARG;
      ;;
    o)
      OUTDIR=$OPTARG;
      ;;
    h)
      # execute next code block
      ;&
    \?)
      usage
      ;;
  esac
done
shift $(( OPTIND - 1 ))

# remaining args are cDNA/plDNA inputfiles
# loop through remaining args, from filename paths extract data type (cDNA/plDNA) and sample name
declare -A FNAMES
declare -a SAMPLES
for f in $*; do
  # echo "arg=${f}"
  sample=$(basename $(dirname "${f}"))
  type=$(basename $(dirname $(dirname "${f}")))
  # echo "sample=${sample}"
  # echo "type=${type}"
  # make filename path absolute
  D=`dirname "${f}"`
  B=`basename "${f}"`
  DD="`cd $D 2>/dev/null && pwd || echo $D`"
  f="$DD/$B"
  SAMPLES+=("${sample}")
  FNAMES["${sample}"]="${f}"
done

# check all required options are set
if [ -z ${IPCR_FNAME+x} ]; then echo "option -i not set (iPCR bedpe filename)"; usage; exit 1; fi
if [ -z ${OUTDIR+x} ]; then echo "option -o not set (directory for output files)"; usage; exit 1; fi
# check required subdirectories exist
if [ ! -d ${OUTDIR} ]; then mkdir -p ${OUTDIR}; echo "making directory \"${OUTDIR}\" for output"; echo ""; fi
# make path to OUTDIR absolute
OUTDIR="`cd \"$OUTDIR\" 2>/dev/null && pwd || echo \"$OUTDIR\"`"
OUTPUT="${OUTDIR}/SuRE-counts.txt"
OUTPUT_BC="${OUTDIR}/SuRE-counts_BC.txt"

######################################
# write stdout to stdout or a log file
######################################
if [ ${LOG} == "true" ]; then 
  LOG="${OUTDIR}/mergeAll.log"
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
echo "IPCR_FNAME=${IPCR_FNAME}"
echo "directory for output files=${OUTDIR}"
echo ""
echo "cDNA and plDNA count-table files for input (sample_name):"
echo "==============================================="
for k in "${!FNAMES[@]}"; 
do 
  echo -e "${FNAMES[$k]}\t($k)"; 
done
echo ""
# print some software version info
echo "Used software:"
echo "=============="
echo "unix/host"; uname -a; echo "";
echo "bash:"; bash --version 2>&1 head -3; echo "";
echo "gawk:"; echo "executable used: ${GAWK}"; ${GAWK} --version; echo "";
echo "=============="
echo ""





##################################################
#########  MAIN  #################################
##################################################

{ echo -e "iPCR\t${IPCR_FNAME}";
  for (( i=0; i<${#SAMPLES[@]}; i++ )); do
    echo -e "${SAMPLES[$i]}\t${FNAMES["${SAMPLES[$i]}"]}"; 
  done } | \
${GAWK} ' 
BEGIN {
  # global variables
  FS="\t"

  split("", pipes)
  split("", samples)

}
NR==1 {
  # iPCR filename
  ipcr_pipe  = "gzip -dc "$2
  next
}
{
  # cDNA/plDNA filename
  samples[length(samples)+1] = $1
  pipes[$1]                  = "gzip -dc "$2
}
END {
  # print header to stdout 
  printf ("chr\tstart\tend\tstrand\tSNPrelpos\tSNPbase\tSNPvar\tSNPabspos\tSNPidx\tiPCR")
  for (i=1; i<=length(samples); i++)
    printf("\t%s", samples[i])
  printf("\n")

  # read input from iPCR datafile; for every BC read from iPCR iterate over all
  # other input pipes and for each pipe extract all lines untill the BC read
  # from that input pipe is (alphabetically) larger than the iPCR barcode
  cnt=0
  while ((ipcr_pipe | getline) > 0) {
    # initialize/empty an array to collect fields for output
    split("", lineout)
    BCipcr = $6
    lineout["iPCR"] = $7
    lineout["CHR"]  = $1
    lineout["START"] = $2
    lineout["END"] = $3
    lineout["STRAND"] = $5
    lineout["SNPrelpos"] = $19
    lineout["SNPbase"] = $20
    lineout["SNPabspos"] = $21
    lineout["SNPvar"] = $22
    lineout["SNPidx"] = $23

    # iterate over all open pipes
    for (i in samples) {
      # initialize count for current sample to zero
      lineout[samples[i]] = 0
      # check whether pipe to current sample file is closed; skip to next sample in that case
      if ( !(samples[i] in pipes) )
        continue

# THE ABOVE CHECK SHOULD MAKE SURE THE PROBLEM DOES NOT OCCUR!!!!!

      # for current iPCR barcode check whether the previously read barcode for current sample is identical or not
      if (BCprev[i] == BCipcr)
        # identical barcodes; set count for current iPCR barcode
        lineout[samples[i]] = CNTprev[i]

      if (BCprev[i] > BCipcr) {
        # previous sample barcode is "larger" than iPCR barcode, therefor
        # current sample does not contain current iPCR barcode; continue to
        # next sample
        continue
      }

      # at this point the previous sample barcode is <= current iPCR barcode;
      # we need to read from the current sample pipe untill the sample barcode
      # is > iPCR barcode. If a newly read sample barcode is identical to iPCR
      # barcode update the count for the current sample
      while(1) {
        # read from current sample pipe
        # print "in while loop: i="i", sample="samples[i]", pipe="pipes[samples[i]] > "/dev/stderr"
        status = (pipes[samples[i]]) | getline line
        # checkfile read status; if EOF close this pipe and delete the pipe from array _pipes_
        if (status == 0) {
          # print "in while loop, status==0: i="i", sample="samples[i]", pipe="pipes[samples[i]] > "/dev/stderr"
          close (pipes[samples[i]])
          delete pipes[samples[i]]
          # break from _while(1)_ loop which reads from current sample pipe, to read next sample pipe
          break
        }

        # process record read from sample pipe; compare sample barcode and iPCR barcode
        split(line, fields)
        if (fields[2] == BCipcr) {
          # sample barcode equal to iPCR barcode; update count for current sample
          lineout[samples[i]] = lineout[samples[i]] + fields[1]
        }
        if (fields[2] > BCipcr) {
          # newly read sample barcode > iPCR barcode, therefor store sample barcode for checking with subsequent iPCR barcodes
          BCprev[i] = fields[2]
          CNTprev[i] = fields[1]
          # break from _while(1)_ loop which reads from current sample pipe, to read next sample pipe
          break
        }
      } # end loop _while(1)_
    } # end loop _for (i in pipes)_

    # all sample pipes for current iPCR barcode have been processed
    # print record to stdout
    printf("%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%d", 
      BCipcr, lineout["CHR"], lineout["START"], lineout["END"], lineout["STRAND"], 
      lineout["SNPrelpos"], lineout["SNPbase"], lineout["SNPvar"], lineout["SNPabspos"],
      lineout["SNPidx"], lineout["iPCR"])
    for (i=1; i<=length(samples); i++)
      printf("\t%d", lineout[samples[i]])
    printf("\n")
  } # end loop _while ((ipcrpipe | getline) > 0)_
  close (ipcr file)
  # close any pipe left open
  for (p in pipes) {
    # print "pipe "p" still open" > "/dev/stderr"
    # (pipes[p]) | getline > "/dev/stderr"
    close (pipes[p])
    delete pipes[p]
    # print
  }
}' | \
  tee >(bzip2 -c > ${OUTPUT_BC}".bz2") | \
  # remove barcode from intermediate output in column 1 (but leave first line intact)
  # cut -f 2- - | \
  ${GAWK} ' BEGIN {OFS="\t"} 
    NR==1 {print; next} 
    { sub(/^\w+\s+/, "")} 1 ' | \
  # sort on chromosomal position but leave the header:
  ${GAWK} ' 
    NR==1 {
      print $0; 
      next
    }
    { 
      print $0 | "sort -S 50%  -k1.4,1V -k2,2g -k3,3g"
    }' | \
${GAWK} -v POS_MULTI_BC_FNAME="${OUTDIR}/pos_multi_BC.txt" '
# awk script to merge counts for duplicated positions
# the input is position sorted: CHR/START/END/STRAND/IPCR/SAMPLES.....
# for every input line a position-label is created as a concatenation of the fields chr/start/end/strand
# if position-labels of subsequent lines are identical the lines are merged;
#   ie. the counts for all samples in all lines are summed

function PROC_PREV_POS(     maxcnt, maxline, BCs, BCmax) {
  # function to process the previous collected lines which have the same position-label
  if (length(PREVLINE) == 1) {
    print PREVLINE[1]
  }
  else {
    split(PREVLINE[1], outline)
    for (i=10; i<=length(outline); i++)
      outline[i] = 0
    for (line in PREVLINE) {
      print PREVLINE[line] > POS_MULTI_BC_FNAME
      split(PREVLINE[line], w)
      for (i=10; i<=length(outline); i++)
	outline[i] = outline[i] + w[i]
    }
    # print the merged output line
    printf("%s", outline[1])
    for (i=2; i<=length(outline); i++)
      printf("\t%s", outline[i])
    printf("\n")
  }
}

BEGIN {
  FS="\t"
  OFS="\t"
}

NR == 1 { # read and print header from input
  print
}

NR == 2 { # initialize PREVPOS and PREVLINE with values from 1st input line
  PREVPOS     = $1"\t"$2"\t"$3"\t"$4
  PREVLINE[1] = $0
}

NR > 2 {
  # the fragment is good quality and duplicate barcodes are filtered
  # the fragments are filtered on position. If consecutive positions are
  # identical (ie, either barcode, mapq, or snp differs) determine which entry
  # has the maximum count and discard all others
  POS = $1"\t"$2"\t"$3"\t"$4
  if (POS == PREVPOS)
    PREVLINE[length(PREVLINE)+1] = $0
  else {
    PROC_PREV_POS()
    delete PREVLINE
    PREVLINE[1] = $0
    PREVPOS     = POS
  }
}

END {
  PROC_PREV_POS()
} ' | \
gzip -c > ${OUTPUT}".gz"

    
    if [ -f "${OUTDIR}/pos_multi_BC.txt" ]; then
  bzip2 "${OUTDIR}/pos_multi_BC.txt"
fi

LINE="finished "${SCRIPTNAME}
SEPARATOR=$(head -c ${#LINE} </dev/zero | tr '\0' '=')
echo $SEPARATOR; echo $LINE; echo $SEPARATOR
endtime=$(date +%c)
echo "end date/time = "${endtime}
echo "==================================="
echo ""


