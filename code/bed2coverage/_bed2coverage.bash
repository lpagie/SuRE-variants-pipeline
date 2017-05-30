#!/bin/bash

# BED2COVERAGE="/home/NFS/users/l.pagie/usr/local/bin/bed2coverage"
WIG2BIGWIG="wigToBigWig" # should be in $PATH

INPUT=$1
OUTALL=$2
OUTPLUS=$3
OUTMINUS=$4
COLUMN=$5
CHROMSIZES=$6

tmpfileall=$(mktemp wig.all.XXXXXX)
tmpfileplus=$(mktemp wig.plus.XXXXXX)
tmpfileminus=$(mktemp wig.minus.XXXXXX)

# trick to enforce waiting for subprocesses created by tee
# from https://unix.stackexchange.com/questions/351780/wait-for-bash-subshells
# and https://unix.stackexchange.com/questions/29851/shell-script-mktemp-whats-the-best-method-to-create-temporary-named-pipe
tempdir=$(mktemp -d waitfifo.tempdir.XXXXXX)
trap 'rm -rf ${tempdir}' EXIT INT TERM HUP
waitfifo="${tempdir}/pipe"
mkfifo "${waitfifo}"

zcat ${INPUT} |\
  tee >( { 
         awk -v col="${COLUMN}" '
           BEGIN { OFS="\t"; FS="\t"}
           NR==1 { for (i=1; i<=NF; i++) { if($i==col) { col=i; break; } }; next }
	   { print $1, $2, $3, $col ; }' |\
	 tee >( awk 'BEGIN{OFS="\t"; FS="\t"} $4>0{$4=1} 1' | ${BED2COVERAGE_EXE} > "flat.${tmpfileall}" ; : >${waitfifo} ) |\
         ${BED2COVERAGE_EXE} > $tmpfileall; : >${waitfifo}
       } ) \
      >( { awk -v col="${COLUMN}" '
           BEGIN { OFS="\t"; FS="\t"}
           NR==1 { for (i=1; i<=NF; i++) { if($i==col) { col=i; break; } }; next }
                 { print $1, $2, $3, ($4=="-"?$col:0) };' |\
	 tee >( awk 'BEGIN{OFS="\t"; FS="\t"} $4>0{$4=1} 1' | ${BED2COVERAGE_EXE} > "flat.${tmpfileminus}" ; : >${waitfifo} ) |\
	 ${BED2COVERAGE_EXE} > $tmpfileminus; : >${waitfifo}
       } ) \
      >( { awk -v col="${COLUMN}" '
	   BEGIN { OFS="\t"; FS="\t"}
	   NR==1 { for (i=1; i<=NF; i++) { if($i==col) { col=i; break; } }; next }
	         { print $1, $2, $3, ($4=="+"?$col:0) };' |\
	 tee >( awk 'BEGIN{OFS="\t"; FS="\t"} $4>0{$4=1} 1' | ${BED2COVERAGE_EXE} > "flat.${tmpfileplus}" ; : >${waitfifo} ) |\
	 ${BED2COVERAGE_EXE} > $tmpfileplus; : >${waitfifo}
       } ) > /dev/null

# wait for all subprocesses to finish
for (( i=0;i<6;i++ )); do read <${waitfifo}; done

echo "$tmpfileall ${CHROMSIZES} $OUTALL
$tmpfileminus ${CHROMSIZES} $OUTMINUS
$tmpfileplus ${CHROMSIZES} $OUTPLUS" |\
  parallel --colsep ' ' ${WIG2BIGWIG}

echo "flat.$tmpfileall ${CHROMSIZES} ${OUTALL/.cov./.covflat.}
flat.$tmpfileminus ${CHROMSIZES} ${OUTMINUS/.cov./.covflat.}
flat.$tmpfileplus ${CHROMSIZES} ${OUTPLUS/.cov./.covflat.}" |\
  parallel --colsep ' ' ${WIG2BIGWIG}

rm -f $tmpfileall $tmpfileminus $tmpfileplus
rm -f flat.$tmpfileall flat.$tmpfileminus flat.$tmpfileplus
