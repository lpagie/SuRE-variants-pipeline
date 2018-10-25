#!/bin/bash

# AUTHOR / DATE
#   Ludo Pagie; March 22, 2016; iPCR-bam2bedpe.bash

# INTRO / BACKGROUND
#   bash script (awk) to process bam files containing data from iPCR samples.
#   The aligned and filtered paired reads are written to stdout in bedpelike
#   format including the barcode and sorted (alphabetically) on barcode.
#
# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#   1 bam file + 1 info file for input as final arguments 
#   -o: output directory
#   optional:
#   -l: write to logfile instead of stdout
#   -b: basename [based on input file name]
#   -c: if set; do not clean intermediate files
# INPUT:
#   bam file
# OUTPUT:
#   tabular txt file in bedpe-like format

# VERSIONS:
#   -170206: initial version, VERSION set to 0.0.1, based on iPCR-map-BC.bash

# TODO
#   - parameterize filter criteria (min MAPQ score, BC length, etc)

SCRIPTNAME=iPCR-bam2bedpe.bash
VERSION=0.0.1 # YYMMDD

# EXTERNAL SOFTWARE
GAWK=gawk
# SAMTOOLS=$HOME/vanSteensel/bin/samtools
SAMTOOLS=samtools

# GLOBAL VARIABLES
CLEAN=true;
LOG="false"

# PARSE OPTIONS
OPTIND=1         # Reset in case getopts has been used previously in the shell.
USAGE=
usage() {
  echo >&2 "usage: ${SCRIPTNAME} -o[lbch] file.bam file.info[.gz/.bz2]"
  echo >&2 "OPTIONS:"
  echo >&2 "  -o: directory for generated count-table files  [required]"
  echo >&2 "  -l: write messages to logfile (OUTDIR/BASENAME.log) instead of stdout"
  echo >&2 "  -b: sets basename used for all output files [default: based on input filename]"
  echo >&2 "  -c: do not clean up intermediate files [default: clean]"
  echo >&2 "  -h: print this message"
  echo >&2 ""
  exit 1;
}

while getopts "h?o:lb:c" opt; do
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

# the remaining CLI argument should be two filenames which are the bam file
# containing aligned reads and the info file with readname-barcode association
# check we have exactly 2 remaining arguments
if [ ! $# -eq 2 ]; then
  echo -e "\nerror: too few, or too many, arguments left after options are parsed (should be 1 bamfile + 1 info file).\nThe remaining args are:"
  while test $# -gt 0; do
    echo $1
    shift
  done
  echo -e "Aborting\n\n"
  usage
fi

# retrieve input bam file from command line
BAM=$1;
# retrieve info filename
INFO=$2;

# check; files exists
if [ ! -f ${BAM} ]; then
  echo -e "error; bam file (${BAM}) doesn't exist.\n" 
  echo -e "Aborting\n\n"
  usage
  exit 1
fi
if [ ! -f ${INFO} ]; then
  echo -e "error; bam file (${INFO}) doesn't exist.\n" 
  echo -e "Aborting\n\n"
  usage
  exit 1
fi

# Make name of files absolute
D=`dirname "${BAM}"`
B=`basename "${BAM}"`
DD="`cd $D 2>/dev/null && pwd || echo $D`"
BAM="$DD/$B"
D=`dirname "${INFO}"`
B=`basename "${INFO}"`
DD="`cd $D 2>/dev/null && pwd || echo $D`"
INFO="$DD/$B"
extension="${INFO##*.}"
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

# check all required options are set
if [ -z ${OUTDIR+x} ]; then echo "option -o not set (directory for output files)"; usage; exit 1; fi
# check required subdirectories exist
if [ ! -d ${OUTDIR} ]; then mkdir -p ${OUTDIR}; echo "making directory \"${OUTDIR}\" for output"; echo ""; fi
# make path to OUTDIR absolute
OUTDIR="`cd \"$OUTDIR\" 2>/dev/null && pwd || echo \"$OUTDIR\"`"

if [ -z ${BASENAME+x} ]; then
  # create BASENAME based on 1st input fastq filename remove ".fastq.*" (or ".fq.*") from filename
  BASENAME=$(basename ${BAM} | sed -e 's/.bam//')
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
echo "CLEAN=${CLEAN}"
echo ""
echo "bam file for input: ${BAM}"
echo ""
# print some software version info
echo "Used software:"
echo "=============="
echo "unix/host"; uname -a; echo "---------------";
echo "bash:"; bash --version 2>&1 head -3; echo "---------------";
echo "gawk:"; echo "executable used: ${GAWK}"; ${GAWK} --version; echo "---------------";
echo "samtools:"; echo "executable used: ${SAMTOOLS}"; ${SAMTOOLS} 2>&1 | head -3; echo "---------------";
echo "=============="
echo ""

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

# setwd processing directory
cd ${OUTDIR}

STATS="${OUTDIR}/${BASENAME}.stats"

### CONVERT BAM FILE INTO BEDPE FILE ###########
################################################
echo "starting conversion of bam file to bedpe file"
## convert bam to bed file (in bedpe format)
BEDPE=${BAM%.bam}.bedpe

${SAMTOOLS} view ${BAM} | \
  ${GAWK} '
function RevComp( theBases ) {
  # from http://www.blossomassociates.net/molbio/revcomp.awk
  answer = "";
  l = length( theBases );
  for ( i = l; 0 < i; i-- ) {
    b = substr( theBases, i, 1 );
  
    if ( "c" == b ) b = "g";
    else if ( "g" == b ) b = "c";
    else if ( "a" == b ) b = "t";
    else if ( "t" == b ) b = "a";
    else if ( "u" == b ) b = "a";
  
    else if ( "C" == b ) b = "G";
    else if ( "G" == b ) b = "C";
    else if ( "A" == b ) b = "T";
    else if ( "T" == b ) b = "A";
    else if ( "U" == b ) b = "A";
  
    else if ( "m" == b ) b = "k";
    else if ( "r" == b ) b = "y";
    else if ( "y" == b ) b = "r";
    else if ( "k" == b ) b = "m";
    else if ( "v" == b ) b = "b";
    else if ( "h" == b ) b = "d";
    else if ( "d" == b ) b = "h";
    else if ( "b" == b ) b = "v";
    else if ( "n" == b ) b = "x";
  
    else if ( "M" == b ) b = "K";
    else if ( "R" == b ) b = "Y";
    else if ( "Y" == b ) b = "R";
    else if ( "K" == b ) b = "M";
    else if ( "V" == b ) b = "B";
    else if ( "H" == b ) b = "D";
    else if ( "D" == b ) b = "H";
    else if ( "B" == b ) b = "V";
    else if ( "N" == b ) b = "N";

    answer = answer b;
  }
  return answer;
}

function CIGAR2length( cigar ) {
  # compute from a CIGAR string the length of the aligned genomic region:
  # first; CIGAR operators which appear only once do not have a counter in
  # front of them. First I will place a '1' in front of those so afterwards I
  # can easily split the CIGAR string and have counts and operators
  # separately.
  # if first character is not a count but a CIGAR operator (MDINSHP=X), put a 1 in front of it
  mod=gensub("^([MDINSHP=X])", "1\\1",1,cigar);
  # for all consecutive pairs of CIGAR operators insert a 1 as count for the 2nd operator
  mod=gensub("([MDINSHP=X])([MDINSHP=X])", "\\11\\2",1,cigar);
  while (mod != cigar) {
    cigar=mod;
    mod=gensub("([MDINSHP=X])([MDINSHP=X])", "\\11\\2",1,cigar);
  }
  
  # split CIGAR string into an array of the counts and an array of the CIGAR operators
  split(cigar, counts, /[MDINSHP=X]+/);
  split(cigar, ops, /[[:digit:]]+/);
  # map the operators to 0,1; depending on whether the operators corresponds to
  # an increase in length of the aligned genomic region
  for (i in ops) {
    gsub("[MDNPX=]", "1", ops[i]) # these ops increase genomic alignment
    gsub("[ISH]", "0", ops[i]) # these ops do not increase genomic alignment
  }
  # compute the total length by multiplying the counts and the lengths of the operators
  len=0;
  for (i=1; i<length(counts); i++) {
    len = len + (counts[i]*ops[i+1])
  }
  return len;
}

BEGIN { 
  FS  = "\t";
  OFS = "\t";
  flag_properPair = lshift(1, 1); # bit-4 flag
  flag_firstInPair = lshift(1, 6);
  flag_readReverseStrand = lshift(1, 4);
};

{
  flag_1=$2;
  # This read should always be the first in pair
  if (! and (flag_1, flag_firstInPair))
    print("The first read is not first in pair!!!!\tFLAG = "flag_1", ("flag_firstInPair")");

    # check this read pair is a concordant pair, otherwise skip the entire pair
    if (! and(flag_1, flag_properPair) ) {
      # if we read a discordant read-pair skip the entire pair (ie also read next read and then go the next pair)
      getline;
      next;
    }
    # store alignment data of first read in pair
    read[1, "MAPQ"]  = $5;
    read[1, "SEQ"]   = $10;
    read[1, "START"] = $4;
    read[1, "FLAG"]  = $2;
    read[1, "CIGAR"] = $6;
    read[1, "RNAME"] = $3;
    read[1, "ID"]    = $1;
    # get MD field from optional fields
    read[1, "MD"]    = "";
    for (i=12; i<NF; i++) {
      if ($i ~ /MD:Z:/) {
	read[1, "MD"]=$i;
	break;
      }
    }
    # get optional XS:i: field from optional fields
    read[1, "XS"]    = "F";
    for (i=12; i<NF; i++) {
      if ($i ~ /XS:i:/) {
	read[1, "XS"]="T";
	break;
      }
    }
    # compute end position on genome
    read[1, "END"] = read[1,"START"] - 1 + CIGAR2length(read[1, "CIGAR"])

    # read second read in pair and parse data
    getline;
    read[2, "MAPQ"]  = $5;
    read[2, "SEQ"]   = $10;
    read[2, "START"] = $4;
    read[2, "FLAG"]  = $2;
    read[2, "CIGAR"] = $6;
    # get MD field from optional fields
    read[2, "MD"]    = "";
    for (i=12; i<NF; i++) {
      if ($i ~ /MD:Z:/) {
	read[2, "MD"]=$i;
	break;
      }
    }
    # get optional XS:i: field from optional fields
    read[2, "XS"]    = "F";
    for (i=12; i<NF; i++) {
      if ($i ~ /XS:i:/) {
	read[2, "XS"]="T";
	break;
      }
    }
    # compute end position on genome
    read[2, "END"] = read[2,"START"] - 1 + CIGAR2length(read[2, "CIGAR"])

    # is the fragment on the forward or on the reverse strand
    if ( and(read[1, "FLAG"], flag_readReverseStrand) ) {
      strand="-";
      # read on reverse strand; reverse-complement the sequence
      # read[1, "SEQ"]=RevComp(read[1, "SEQ"]);
    }
  else {
    strand="+";
    # read on reverse strand; reverse-complement the sequence
    # read[2, "SEQ"]=RevComp(read[2, "SEQ"]);
  }

  # which read is first on forw strand
  if (strand == "+") {
    first_read=1;
  } else {
  first_read=2;
}
last_read=3-first_read;

# print the bedpe output, including the extra data, in the following format:
# readID seqname start end strand end.2 start.2 MAPQ.1 MAPQ.2 MD.1 MD.2 SEQ.1 SEQ.2
print( read[1, "ID"],
read[1, "RNAME"],
read[first_read, "START"],
read[last_read, "END"],
strand,
read[first_read, "END"],
read[last_read, "START"],
read[first_read, "MAPQ"],
read[last_read, "MAPQ"],
read[first_read, "MD"],
read[last_read, "MD"],
read[first_read, "XS"],
read[last_read, "XS"],
read[first_read, "SEQ"],
read[last_read, "SEQ"],
read[first_read, "CIGAR"],
read[last_read, "CIGAR"])
}
' > ${BEDPE}

## add barcode sequence from INFO file to BEDPE file
# $INFO is the info-file created while trimming the forward read by
# cutadapt (ie very first step)
mv ${BEDPE} ${BEDPE}.tmp
# merge barcodes from $INFO into $BEDPE using awk
# (store readID from column-5 from $BEDPE as key in array 'a' with
# entire line, except readID, as value. read $INFO, look for key and add barcode
# sequence $5) to array element
# INFO="${OUTDIR}/${BASENAME}_forw_trimmed.info"
# info file may contain reads for which no adapter sequence was found: discard
# those reads from the info file prior to adding the barcodes into the bedpe
# file

${GAWK} -F '\t' -v statsfile=${STATS} ' 
BEGIN {
  incl=0;
  excl=0;
}
FNR==NR{
# print FILENAME, ARGV[1], ARGV[2] >> statsfile
  if (FILENAME != ARGV[1]) {exit} # if 1st input file is empty, abort
  a[$1]=substr($0, index($0,$2)); 
  next
}
{
  # LP140424; trim readIDs differently for NKI formatted readIDs or BGI
  # formatted readIDs
  sub(/\s.*$/,"",$1); # trim readID for NKI format
  sub(/\/1$/,"",$1); # trim readID for bgi format
  # end LP140424
}
$2==-1 { 
  # here, the info file indicates the read did not contain an identifiable
  # barcode. Thus this read should be skipped and removed from the array
  
  excl++
  delete a[$1]
  next; 
}
{ 
  if ($1 in a) {
    BC = $5
    len  = length(BC)
    BClen[len]++
    hasN = BC~/N/
    NNN[hasN]++
    if (len==20 && !hasN ) {
      incl++
      a[$1]=a[$1] FS BC;
    }
    else {
      excl++
      delete a[$1]
    }
  }
}
END {
# print FNR, NR, FILENAME >> statsfile
  if (FNR == NR) {print "2nd file empty" >> statsfile; exit} # the 2nd input file is empty; abort
  for (key in a) {cnt++; print a[key];}
  print "while filtering for proper barcodes: included = "incl", discarded = "excl >> statsfile
  
  printf "\n\n" >> statsfile
  
  printf ("bedpeFragmentCount\t%d\n", cnt) >> statsfile
  
  printf "iPCR_BC_lengths" >> statsfile
  for (k in BClen) printf("\t%d", k) >> statsfile
  printf "\n" >> statsfile
  
  printf "iPCR_BC_lengths_counts" >> statsfile
  for (k in BClen) printf("\t%d", BClen[k]) >> statsfile
  printf "\n" >> statsfile
  
  printf "iPCR_BC_NNN" >> statsfile
  for (k in NNN) printf("\t%d", k) >> statsfile
  printf "\n" >> statsfile
  
  printf "iPCR_BC_NNN_counts" >> statsfile
  for (k in NNN) printf("\t%d", NNN[k]) >> statsfile
  printf "\n" >> statsfile
  
  printf "\n\n" >> statsfile
} ' ${BEDPE}.tmp <(${CAT} ${INFO}) | \
# sort resulting bedpe file on seqname and start (latter numeric
# sort)
sort -k1,1 -k2,2n | \
# remove duplicates and add a count column
uniq -c |\
# reorder columns: move count to 7th column
# and write to $BEDPE the following columns:
# chr start end length strand barcode count internal-end internal-start MAPQ MD1 MD2 XS1 XS2 SEQ1 SEQ2 CIGAR1 CIGAR2
${GAWK} -v statsfile=${STATS} ' 
BEGIN{ 
  OFS="\t";
  print("seqname","start","end","length","strand","barcode","count","end.intr","start.intr","MAPQ","MD.1","MD.2","alt.1","alt.2","seq.1","seq.2","cigar.1","cigar.2");
} 
{ 
  print($2,$3,$4,$4-$3+1,$5,$18,$1,$6,$7,$8,$10,$11,$12,$13,$14,$15,$16,$17)
} 
END {
  printf("bedpeFragmentUniqCount\t%d\n\n", NR) >> statsfile
}' > ${BEDPE}
echo "conversion bam to bedpe done"


# clean or compress text files
##############################
if $CLEAN; then
  rm -f *.tmp
else
  gzip *.tmp
fi
# compress bedpe file
gzip ${BEDPE}

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
