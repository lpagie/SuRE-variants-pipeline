#!/bin/bash

# input:
#   - SuRE-counts
#   - VCF file
#   - SAMPLE ID

GAWK=gawk

# check number of arguments
if [[ $# -ne 3 ]]; then
    (>&2 echo "expected 3 arguments, got $# instead. Aborting")
    exit 0
fi

VCF=$1
SURE=$2
SAMPLE=$3

${GAWK} -v smpl=${SAMPLE} '
# following is from https://stackoverflow.com/questions/42875915/how-do-you-convert-an-array-to-a-string-in-awk
# Usage:
#    arr2str(flds[,seps,[sortOrder]])
#
# flds:
#    This function converts the mandatory "flds" array argument into a string.
#
# seps:
#    If "seps" is not present then the "flds" values will simply be concatenated
#    in the returned string.
#
#    If "seps" is present and is a string then that "seps" value will be inserted
#    between each "flds" value in the returned string.
#
#    If "seps" is present and is an array then each "seps" value with the same index
#    as a "flds" index will be inserted in the returned string before or after
#    (sort order dependent) the corresponding "flds" value with that same index.
#    - All "seps" values that do not have an index in "flds" will be inserted in
#      the returned string before or after all of the "flds" and other "seps" values.
#      This ensures that a "seps" array that, for example, starts at zero as a result
#      of a previous split(str,flds,re,seps) will have its zeroth entry included.
#
# sortOrder:
#    If "sortOrder" is present then it will be used as the order the "flds" values
#    are visited in, otherwise it uses PROCINFO["sorted_in"] if set, otherwise it
#    uses ascending numeric indices.
#    - If the sort order is descending (ends in "desc") and "seps" is an array then
#      the "seps" values are inserted before each "flds" value, otherwise after them.
#
# Example:
#    $ cat tst.awk
#    BEGIN {
#        orig = ",a+b:c-d="
#        split(orig,flds,/[^[:alpha:]]/,seps)
#
#        printf "orig: <%s>\n", orig
#        printf "asc:  <%s>\n", arr2str(flds,seps)
#        printf "desc: <%s>\n", arr2str(flds,seps,"@ind_num_desc")
#    }
#    $ awk -f arr2str.awk -f tst.awk
#    orig: <,a+b:c-d=>
#    asc:  <,a+b:c-d=>
#    desc: <=d-c:b+a,>

function arr2str(flds, seps, sortOrder,      sortedInPresent, sortedInValue, currIdx, prevIdx, idxCnt, outStr) {

if ( "sorted_in" in PROCINFO ) {
  sortedInPresent = 1
  sortedInValue = PROCINFO["sorted_in"]
}

if ( sortOrder == "" ) {
  sortOrder = (sortedInPresent ? sortedInValue : "@ind_num_asc")
}
PROCINFO["sorted_in"] = sortOrder

if ( isarray(seps) ) {
  # An array of separators.
  if ( sortOrder ~ /desc$/ ) {
    for (currIdx in flds) {
      outStr = outStr (currIdx in seps ? seps[currIdx] : "") flds[currIdx]
    }
  }

  for (currIdx in seps) {
    if ( !(currIdx in flds) ) {
      outStr = outStr seps[currIdx]
    }
  }

  if ( sortOrder !~ /desc$/ ) {
    for (currIdx in flds) {
      outStr = outStr flds[currIdx] (currIdx in seps ? seps[currIdx] : "")
    }
  }
}
else {
  # Fixed scalar separator.
  # We would use this if we could distinguish an unset variable arg from a missing arg:
  #    seps = (magic_argument_present_test == true ? seps : OFS)
  # but we cant so just use whatever value was passed in.
#  print "new"
  for (currIdx in flds) {
#    print currIdx, flds[currIdx], outStr
    outStr = outStr (idxCnt++ ? seps : "") flds[currIdx]
  }
}

if ( sortedInPresent ) {
  PROCINFO["sorted_in"] = sortedInValue
}
else {
  delete PROCINFO["sorted_in"]
}

return outStr
}
BEGIN {
  FS="\t"
  OFS="\t"
  split("", zygosity)
  split("", bases)
  PROCINFO["sorted_in"] = "@ind_num_asc"
}
NR==FNR {
  if (/^##/) {
    next
  }
  if (/#CHROM/) {
    # select correct column based on SAMPLE
    for (i=1; i<=NF; i++) {
      if ($i == smpl) {
	GT=i
	continue
      }
    }

    # GT=11
    # print "GTcol="GT, $GT > "/dev/stderr"
    next
  }
  split($GT, alleles, "|")
  idx = length(zygosity)
  if (alleles[1] == alleles[2]) {
    # homozygous case; store GT (0,1) and corresponding base identity
    zygosity[idx][0] = alleles[1]
    zygosity[idx][1] = $(4+alleles[1]) # value of alleles[i] is 0 or 1, corresponding ref/alt allele column is 4 or 5
    zygosity[idx][2] = $3
  } else {
    zygosity[idx][0] = 3
    zygosity[idx][2] = $3
  }
  next
}
# reading SuRE-counts
FNR==1 {
# read header of SuRE-counts; print and add 2 column names
    print $0, "SNPbaseInf","SNPvarInf","SNP_ID"
    next
  }
$7 == "" {
  # no SNPs in this fragment; simply add empty entries for the 3 additional columns
  print $0, "", "", ""
  next
}
{
  # all fragments need SNP_IDs added
  split($7, snpvar, ",")
  split($9, snpidx, ",")
  split("", ID) # will hold SNP_IDs for all SNPs
  split("", unread) # used to record all SNPs in fragment which are not read
  baseout = $6 # string with base identities, possibly inferred
  varout = $7 # string with new vars, possibly inferred
  for (k in snpvar) {
    # loop pver 
    ID[k] = zygosity[snpidx[k]][2]
    if (snpvar[k] == 3) 
      unread[length(unread)]=k
  }
  IDout = arr2str(ID, ",")

  # we looped over all SNPs and marked those which are not read
  inferred=0
  # if any SNP is unread check whether we can infer the base identity using zygosyti of the alleles
  if (length(unread) > 0 ) {
    split($6, snpbase, ",")
    split($7, newvar, ",")
    split($6, newbase, ",")
    for (k in unread) {
      # loop over unread SNPs
      if ( zygosity[snpidx[unread[k]]][0] != 3 ) {
	# SNP in this genome is not heterozygous (ie 3) assign new base identities and SNPvar
	newvar[unread[k]] = zygosity[snpidx[unread[k]]][0]
	newbase[unread[k]] = zygosity[snpidx[unread[k]]][1]
	inferred = 1
      } else {
        # otherwise simply copy the old (ambiguous) base identity and SNPvar
        newvar[unread[k]] = snpvar[unread[k]]
	newbase[unread[k]] = snpbase[unread[k]]
      }
    }
    # create new output if any of the SNPs are now inferred
    if (inferred) {
      varout = arr2str(newvar, ",")
      baseout = arr2str(newbase, ",")
    }
  }
  print $0, baseout, varout, IDout
}
' <( gzip -cd ${VCF} ) <( gzip -cd ${SURE} )
