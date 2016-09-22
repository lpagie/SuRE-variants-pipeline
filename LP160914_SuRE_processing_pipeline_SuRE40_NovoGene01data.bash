#!/bin/bash




# HARD-CODED DIRECTORY PATHS!!!!!
#################################
# SCRATCH="/home/NKI/l.pagie/scratch/" # scratch directory on local harddisk for optimal file access time
# export TMP="${SCRATCH}/tmp" # tmp directory on loacl scratch, again for optimal file access times
# directory with bash scripts:
CODEDIR="$PWD/code"

# GLOBAL VARIABLES
##################
TIMETAG=`date +"%H%M%S"`
DATETAG=`date +"%y%m%d"`
# adapter sequences 
ADPT_SEQ="CCTAGCTAACTATAACGGTCCTAAGGTAGCGAA"
ADPTR_IPCR_FORW_SEQ="CCTAGCTAACTATAACGGTCCTAAGGTAGCGAACCAGTGAT"
ADPTR_IPCR_REV_SEQ="CCAGTCGT"
RESTRICT_SITE="NO_RESTRICT"
BOWTIE2_REFSEQ="$HOME/projects/LP140430_SureSeq_JvArensbergen/data/LP160714_CHO_refseq/bowtie2-index/criGri1"

# running/OS VARIABLES
RUN_PARSE_DNA="false"
RUN_PARSE_IPCR="true"
RUN_MERGE_IPCR="true"
RUN_MERGE_IPCR_DNA="true"
export NCORES=30
export LC_ALL=C # to get the gnu sort in correct order on all machines

# FILENAMES
###########
OUTDIR=SuRE_pipeline_OUTPUT_SuRE40_NovoGene01data_"${DATETAG}"_"${TIMETAG}"
# OUTDIR="SuRE_pipeline_OUTPUT_160401_111919"
PAR_FNAME="pipeline_runtime_parameters_SuRE40_NovoGene01data_${DATETAG}_${TIMETAG}.bash"
# scripts:
# PARSE_DNA_SCRIPT="${CODEDIR}/cDNA-plDNA-count-BC_160322.bash"
PARSE_IPCR_SCRIPT="${CODEDIR}/iPCR-map-BC_160322.bash"
MERGE_IPCR_SCRIPT="${CODEDIR}/iPCR-merge-bedpe-Filter-BC-multi-pos.bash"
MERGE_IPCR_DNA_SCRIPT="${CODEDIR}/merge-iPCR-cDNA-plDNA.bash"
# data:
SAMPLES="./SuRE-sample-listing_SuRE40_NovoGene01data_LP160912.txt"
IPCRFASTQ="/home/NFS/users/l.pagie/projects/LP140430_SureSeq_JvArensbergen/data/LP160909_LP160714_CHO_iPCR_SuRE40_NHHW160047-03/links2fastq/*fq.gz"
# CDNAFASTQ="data/cDNA/*fastq.gz"
# PLDNAFASTQ="data/PL/*fastq.gz"


# create various directories for output
#######################################
mkdir ${OUTDIR}
# check/create scratch and tmp directory
##  if [ ! -d  ${SCRATCH} ]; then
##    mkdir -p ${SCRATCH}
##  fi
##  if [ ! -d  ${TMP} ]; then
##    mkdir -p ${TMP}
##  fi


##########
## MAIN ##
##########


# PARSE CDNA, PLDNA FASTQ FILES
###############################
if [ $RUN_PARSE_DNA == "true" ]; then
  # run script for cDNA files
  /usr/bin/time -v nice -19 bash ${PARSE_DNA_SCRIPT} -o ${OUTDIR}/cDNA -a ${ADPT_SEQ} -n ${NCORES} ${CDNAFASTQ}
  # run same script for plDNA files
  /usr/bin/time -v nice -19 bash ${PARSE_DNA_SCRIPT} -o ${OUTDIR}/plDNA -a ${ADPT_SEQ} -n ${NCORES} ${PLDNAFASTQ}
fi

# PARSE IPCR FASTQ
##################
if [ $RUN_PARSE_IPCR == "true" ]; then
  # run the script for fastq files which are generated using restriction enzyme XYZ for shortening circular products in iPCR
  /usr/bin/time -v nice -19 bash ${PARSE_IPCR_SCRIPT} -o ${OUTDIR}/iPCR/bedpe/ -f ${ADPTR_IPCR_FORW_SEQ} -r ${ADPTR_IPCR_REV_SEQ} -d ${RESTRICT_SITE} -n ${NCORES} -i ${BOWTIE2_REFSEQ} ${IPCRFASTQ}
fi


# MERGE IPCR BEDPE FILES INTO ONE BEDPE
#######################################
if [ $RUN_MERGE_IPCR == "true" ]; then
  /usr/bin/time -v nice -19 bash ${MERGE_IPCR_SCRIPT} -o ${OUTDIR}/iPCR/bedpe/SuRE-combined-bedpe-${DATETAG}_${TIMETAG}.txt.gz -n ${NCORES}  ${OUTDIR}/iPCR/bedpe/*_bedpe.gz
fi


# MERGE iPCR/CDNA/PLDNA TO SURE-COUNTS
######################################
if [ $RUN_MERGE_IPCR_DNA == "true" ]; then
  #IPCR=/home/NFS/users/l.pagie/projects/LP140430_SureSeq_JvArensbergen/analyses/LP160621_SuRE34_pipeline_results/SuRE_pipeline_OUTPUT_160621_135013/iPCR/bedpe/SuRE-combined-bedpe-160621_135013.txt
  #IPCR=/home/NFS/users/l.pagie/projects/LP140430_SureSeq_JvArensbergen/analyses/LP160715_ChinHamster_SuRE40/SuRE_pipeline_OUTPUT_160801_163909/iPCR/bedpe/SuRE-combined-bedpe-160801_163909.txt.gz
  IPCR=${OUTDIR}/iPCR/bedpe/SuRE-combined-bedpe-${DATETAG}_${TIMETAG}.txt.gz
  /usr/bin/time -v nice -19 bash ${MERGE_IPCR_DNA_SCRIPT}  -s ${SAMPLES} -i $IPCR -c ${OUTDIR}/cDNA -p ${OUTDIR}/plDNA -o "${OUTDIR}/SuRE-counts_${DATETAG}_${TIMETAG}.txt.gz" -n ${NCORES}
fi

