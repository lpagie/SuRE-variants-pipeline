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
#   required:
#   -s: sample meta file which has columns with fastq file names for all
#       samples (iPCR, cDNA, and plDNA) and a column with (short) sample names. 
#   -i: name of iPCR bedpe file
#   -o name of output file
#   optional:
#   -c: directory containing cDNA count-table file(s) [cDNA/count-tables]
#   -p: directory containing plDNA count-table file(s) [plDNA/count-tables]
#   -l: log-filename [stdout]
#   -n: number of cores used in parallel processes (10)
# INPUT:
#   iPCR bedpe file; redundant, sorted on barcode
#   cDNA/plDNA count-tables, sorted on barcode
# OUTPUT:
#   tabular txt file with position, strand, sample-counts, ordered on position

# VERSIONS:
#   -160318: initial version, VERSION set to 0.0.1

# TODO
# merge counts if sample names for groups of identical sample names (ie PL)

VERSION=0.0.1 # YYMMDD
SCRIPTNAME=merge-iPCR-cDNA-plDNA.bash

# EXTERNAL SOFTWARE
GAWK=/usr/bin/gawk

# GLOBAL VARIABLES
NCORES=10
PLDNA_DIR=plDNA/count-tables/
CDNA_DIR=cDNA/count-tables/

# PARSE OPTIONS
OPTIND=1         # Reset in case getopts has been used previously in the shell.
USAGE=
usage() {
  echo >&2 "usage: SCRIPTNAME -sicplno"
  echo >&2 "OPTIONS:"
  echo >&2 "  -s: sample meta file [required]"
  echo >&2 "  -o: output file [required]"
  echo >&2 "  -i: iPCR bedpe filename [required]"
  echo >&2 "  -p: directory with plDNA count-table file(s) [default: plDNA/count-tables/]"
  echo >&2 "  -c: directory with cDNA count-table file(s) [default: cDNA/count-tables/]"
  echo >&2 "  -l: set name of logfile [default: stdout]"
  echo >&2 "  -n: number of cores used where possible [default: 10]"
  echo >&2 ""
  exit 1;
}

while getopts "h?s:o:i:p:c:l:n:" opt; do
  case $opt in
    l)
      LOG=$OPTARG;
      ;;
    n)
      NCORES=$OPTARG;
      ;;
    s)
      SAMPLES=$OPTARG;
      ;;
    i)
      IPCR_FNAME=$OPTARG;
      ;;
    p)
      PLDNA_DIR=$OPTARG;
      ;;
    c)
      CDNA_DIR=$OPTARG;
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

# check all required options are set
if [ -z ${SAMPLES+x} ]; then echo "option -s not set (sample meta file)"; usage; exit 1; fi
if [ -z ${IPCR_FNAME+x} ]; then echo "option -i not set (iPCR bedpe filename)"; usage; exit 1; fi
if [ -z ${OUTPUT+x} ]; then echo "option -o not set (output filename)"; usage; exit 1; fi

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
echo "SAMPLES=${SAMPLES}"
echo "IPCR_FNAME=${IPCR_FNAME}"
echo "CDNA_DIR=${CDNA_DIR}"
echo "PLDNA_DIR=${PLDNA_DIR}"
echo "output file=${OUTPUT}"
echo ""
echo "cDNA count-table files for input:"
echo "================================="
for f in $CDNA_DIR/*; do echo $f; done
echo ""
echo "plDNA count-table files for input:"
echo "=================================="
for f in $PLDNA_DIR/*; do echo $f; done
echo ""
# print some software version info
echo "Used software:"
echo "=============="
echo "unix/host"; uname -a; echo "";
echo "bash:"; bash --version 2>&1 head -3; echo "";
echo "gawk:"; echo "executable used: ${GAWK}"; ${GAWK} --version; echo "";
echo "=============="
echo ""

# check required subdirectories exist
if [ ! -d iPCR ]; then mkdir -p "iPCR"; echo "making directory \"iPCR\" for output"; fi

echo "starting loop over samples in sample-file"

${GAWK} -v samplesfile=${SAMPLES} -v ipcrfname=${IPCR_FNAME} \
  -v cdnadir=${CDNA_DIR} -v pldnadir=${PLDNA_DIR} -v bc2bcfname=${BC2BC_FNAME} ' 
function basename(file) {
  sub(".*/", "", file)
  return file
}
function ripExtension(file,    ext) {
  if (ext == "") {
    ext=file
    sub(".*\\.", ".", ext)
  }
  split(file, parts, ext)
  return parts[1]
}
function read_file_into_array(file, array     ,status, record, count ) {
  # define a function to read the entire config file into an array
  # from http://www.unix.com/shell-programming-scripting/135373-put-lines-file-array-awk.html
  # this function an entire text file into an array. It returns the number of read lines.
  count  = 0;
  while (1) {
    status = getline record < file
    if (status == -1) {
      print "Failed to read file " file;
      exit 1;
    }
    if (status == 0) break;
    array[++count] = record;
  }
  close(file);
  return count
}
BEGIN {
# global variables
FS="\t"

# read sample meta file to extract sample- and file-names
n = read_file_into_array(samplesfile, LINES)
print "read samplefile with n="n" lines" > "/dev/stderr"
# initialize/empty arrays for samplenames and filenames, extracted from sample meta file
split("", pipes)
split("", fnames)
split("", samplenames)
# iterate over records (ie. samples), read from sample metafile, store sample names and file names
for (i=2; i <= n; i++) {
  print "ith line = "LINES[i] > "/dev/stderr"
  # for each iteration split the line in chunks into FIELDS
  split(LINES[i], FIELDS)
  # skip if sample file 5th column says so
  if ( FIELDS[5] == "no" ) {
    print "skipping this sample" > "/dev/stderr"
    continue
  }
  print "fields[4] = "FIELDS[4] > "/dev/stderr"
  switch (FIELDS[4]) { # field 4 has the type of data for current record
  case "cDNA":
    # count-tabels file 
    fname = cdnadir"/"ripExtension(basename(FIELDS[2]), ".fastq.gz")"_trimmed_table.txt.gz"
    if ( system("test -f " fname)!=0 )
      break
    fnames[length(fnames)+1] = fname
    pipes[fname] = "gzip -dc "fname
    samplenames[length(samplenames)+1] = FIELDS[1]
    break
  case "plDNA":
    # count-tabels file 
    fname = pldnadir"/"ripExtension(basename(FIELDS[2]), ".fastq.gz")"_trimmed_table.txt.gz"
    if ( system("test -f " fname)!=0 ) {
      print "this file doesnt exist: "fname > "/dev/stderr"
      break
    }
  else {
      print "this file exist: "fname > "/dev/stderr"
    }
    fnames[length(fnames)+1] = fname
    pipes[fname] = "gzip -dc "fname
    samplenames[length(samplenames)+1] = FIELDS[1]
    break
  case "iPCR":
    # skip as iPCR datafile is supplied on commandline
    break
  case "?":
  default:
    printf ("default (case =#%s#)", FIELDS[4], "\n") > "/dev/stderr"
    usage()
    break
  }
}

  # open ipcr datafile as pipe
  ipcrpipe = "gzip -dc "ipcrfname

  # print header to stdout 
  printf ("chr\tstart\tend\tstrand\tiPCR")
  for (i in samplenames)
    printf("\t%s", samplenames[i])
  printf("\n")

  # read input from iPCR datafile; for every BC read from iPCR iterate over all
  # other input pipes and for each pipe extract all lines untill the BC read
  # from that input pipe is (alphabetically) larger than the iPCR barcode
  cnt=0
  while ((ipcrpipe | getline) > 0) {
    # initialize/empty an array to collect fields for output
    split("", lineout)
    BCipcr = $6
    lineout["iPCR"] = $7
    lineout["CHR"]  = $1
    lineout["START"] = $2
    lineout["END"] = $3
    lineout["STRAND"] = $5

    # iterate over all open pipes
    for (i in fnames) {
      # initialize count for current sample to zero
      lineout[samplenames[i]] = 0
      # check whether pipe to current sample file is closed; skip to next sample in that case
      if (! fnames[i] in pipes)
	continue

      # for current iPCR barcode check whether the previously read barcode for current sample is identical or not
      if (BCprev[i] == BCipcr)
	# identical barcodes; set count for current iPCR barcode
	lineout[samplenames[i]] = CNTprev[i]

      if (BCprev[i] > BCipcr)
        # previous sample barcode is "larger" than iPCR barcode, therefor
        # current sample does not contain current iPCR barcode; continue to
        # next sample
        continue

      # at this point the previous sample barcode is <= current iPCR barcode;
      # we need to read from the current sample pipe untill the sample barcode
      # is > iPCR barcode. If a newly read sample barcode is identical to iPCR
      # barcode update the count for the current sample
      while(1) {
	# read from current sample pipe
	# print "in while loop: i="i", fname="fnames[i]", pipe="pipes[fnames[i]] > "/dev/stderr"
        status = (pipes[fnames[i]]) | getline line
	# checkfile read status; if EOF close this pipe and delete the pipe from array _pipes_
        if (status == 0) {
          close (pipes[fnames[i]])
          delete pipes[fnames[i]]
	  delete fnames[i]
	  # break from _while(1)_ loop which reads from current sample pipe, to read next sample pipe
          break
	}

	# process record read from sample pipe; compare sample barcode and iPCR barcode
	split(line, fields)
        if (fields[2] == BCipcr) {
	  # sample barcode equal to iPCR barcode; update count for current sample
          lineout[samplenames[i]] = lineout[samplenames[i]] + fields[1]
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
    printf("%s\t%s\t%d\t%d\t%s\t%d", 
      BCipcr, lineout["CHR"], lineout["START"], lineout["END"], lineout["STRAND"], lineout["iPCR"])
    for (i in samplenames) 
      printf("\t%d", lineout[samplenames[i]])
    printf("\n")
  } # end loop _while ((ipcrpipe | getline) > 0)_
  close (ipcr file)
}' | \
  tee ${OUTPUT}".plusBC" | \
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
      print $0 | "sort -S 50% --parallel=${NCORES}  -k1.4,1V -k2,2g -k3,3g"
    }' | \
${GAWK} '
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
    for (i=5; i<length(outline); i++)
      outline[i] = 0
    for (line in PREVLINE) {
      print PREVLINE[line] > POS_MULTI_BC_FNAME
      split(PREVLINE[line], w)
      for (i=5; i<length(outline); i++)
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

  COLLENGTH=4
  COLMAPQ=10
  COLBC=6
  COLCHR=1
  COLS = 2
  COLE = 3
  COLSTR = 5
  COLCNT = 7

  POS_MULTI_BC_FNAME = "iPCR/pos_multi_BC.txt"
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
gzip -c > ${OUTPUT}


LINE="finished "${SCRIPTNAME}
SEPARATOR=$(head -c ${#LINE} </dev/zero | tr '\0' '=')
echo $SEPARATOR; echo $LINE; echo $SEPARATOR
endtime=$(date +%c)
echo "end date/time = "${endtime}
echo "==================================="
echo ""


