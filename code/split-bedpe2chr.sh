#!/bin/bash

fname=$1
echo "fname = ${fname}"
# out_fname_base=$(basename ${fname%.txt.gz})
out_fname_base=${fname%.txt.gz}
echo "out_fname_base=$out_fname_base"

zcat $fname | gawk -v base=${out_fname_base} '
BEGIN {
# open pipes
  p["chr1"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr1.txt.gz"
  p["chr2"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr2.txt.gz"
  p["chr3"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr3.txt.gz"
  p["chr4"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr4.txt.gz"
  p["chr5"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr5.txt.gz"
  p["chr6"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr6.txt.gz"
  p["chr7"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr7.txt.gz"
  p["chr8"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr8.txt.gz"
  p["chr9"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr9.txt.gz"
  p["chr10"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr10.txt.gz"
  p["chr11"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr11.txt.gz"
  p["chr12"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr12.txt.gz"
  p["chr13"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr13.txt.gz"
  p["chr14"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr14.txt.gz"
  p["chr15"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr15.txt.gz"
  p["chr16"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr16.txt.gz"
  p["chr17"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr17.txt.gz"
  p["chr18"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr18.txt.gz"
  p["chr19"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr19.txt.gz"
  p["chr20"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr20.txt.gz"
  p["chr21"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr21.txt.gz"
  p["chr22"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chr22.txt.gz"
  p["chrX"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chrX.txt.gz"
  p["chrY"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chrY.txt.gz"
  p["chrM"]="sort -S 15G  -k1.4,1V -k2,2g -k3,3g | gzip -c > " base"_chrM.txt.gz"
}
{
  print $0 | p[$1]
}'
