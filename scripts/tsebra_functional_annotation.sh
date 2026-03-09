#!/bin/bash

## This script will read:
##      polypeptide fasta file from augustus (via BRAKER)
##      gtf file
##      reference fasta file
##      output directory path
## And write out:
##      Annotated GFF3 file
##      Genome FASTA file
##
## This script will run the following annotation programs:
##      RNAMMER
##      TRNASCAN
##      TMHMM
##      HMMSCAN
##      RAPSEARCH2
##      ATTRIBUTOR
##
## Tools required:
##      GFFTOOLS
##      python v2.7
##      python v3.8
##      Perl
##


#--------------------------------------------------------
# Set Usage Parametres
#--------------------------------------------------------

programmename=$0
function usage {

    echo "$(basename $programmename) [-a <AA_FASTA>] [-s <GTF_FILE>] [-o <output_directory>] [-f <reference_fasta_file>] [-g <gff3_file_name>] [-h <help>]
    
    This script will run the automated IGS functional annotation pipeline via qsub commands
    
    -a          This is the amino acid fasta file from the structural GTF file

    -s          The structural annotation file

    -o          Specify the output directory. Five Subdirectories will be created within this output directory (attributor, hmm_scan, rapsearch, rnammer, tmhmm, trnascan)

    -f          The reference FASTA file for the organism of interest
    
    -g          Full name of final annotated gff3 file

    -h          Print the help menu
                
                
                *** FULL PATHS MUST BE SPECIFICED FOR ALL ARGUMENTS ***


                "
    
}

if [ "$#" -eq 0  ]; then
    usage
    exit
fi
if [ "$#" -eq 1 ] && [ "$1" == "-h" ]; then
    usage
    exit
fi

while getopts a:s:o:f:g:h flag
 do
     case "${flag}" in
         a) AA_FASTA=${OPTARG};;
         s) GTF=${OPTARG};;
         o) OUTPUT_DIRECTORY=${OPTARG};;
         f) REFERENCE_FASTA=${OPTARG};;
         g) FINAL_GFF3_FILE=${OPTARG};;
     esac
 done

#--------------------------------------------------------
# Set Variables Needed
#--------------------------------------------------------


## Main Directory Paths
PACKAGE_DIR=
USER_PACKAGE_DIR=
DBS_DIR=
SCRIPT_DIR=


## Programmes Needed
RNAMMER_DIR=${USER_PACKAGE_DIR}/rnammer-1.2
TRNASCAN_BIN_DIR=${USER_PACKAGE_DIR}/trnascan-se-2.0.3/bin
TMHMM_BIN_DIR=${PACKAGE_DIR}/tmhmm-2.0c/bin
HMM_BIN_DIR=${PACKAGE_DIR}/hmmer-3.3/bin
RAPSEARCH2_BIN_DIR=${PACKAGE_DIR}/rapsearch2-2.24/bin
ATTRIBUTOR_DIR=${USER_PACKAGE_DIR}/Attributor/
GFFTOOLS_DIR=${USER_PACKAGE_DIR}/oqtans_tools/GFFtools/0.1
PERL_BIN_DIR=${PACKAGE_DIR}/perl-5.30.2/bin
PYTHON_DIR=${PACKAGE_DIR}/python-3.8.2
PERL_BIN_DIR=${PACKAGE_DIR}/perl/bin

## Output directories
PROJECT_DIR=${OUTPUT_DIRECTORY}

RNAMMER_OUTPUT_DIR=${PROJECT_DIR}/rnammer
TRNASCAN_OUTPUT_DIR=${PROJECT_DIR}/trnascan
TMHMM_OUTPUT_DIR=${PROJECT_DIR}/tmhmm
HMMSCAN_OUTPUT_DIR=${PROJECT_DIR}/hmm_scan
RAPSEARCH_OUTPUT_DIR=${PROJECT_DIR}/rapsearch
ATTRIBUTOR_OUTPUT_DIR=${PROJECT_DIR}/attributor

## Input Files
AUGUSTUS_HINTS_AA=${AA_FASTA} 
AUGUSTUS_HINTS_GFF3=${GTF}
ATTRIBUTOR_CONFIG_FILE=${ATTRIBUTOR_OUTPUT_DIR}/config.txt

## DATABASES AND INDEXES
CODING_HMM_LIB=${DBS_DIR}/coding_hmm
UNIREF100=${DBS_DIR}/uniprot_sprot_all_and_trembl_characterized.20180904

## HMMSCAN Output Files
HMMSCAN_FILE=${HMMSCAN_OUTPUT_DIR}/hmmscan.file
HMMSCAN_HTAB=${HMMSCAN_OUTPUT_DIR}/hmmscan.htab

## RAPSEARCH Output Files
RAPSEARCH_FILE=${RAPSEARCH_OUTPUT_DIR}/rapsearch.file

## TMHMM Output Files
TMHMM_FILE=${TMHMM_OUTPUT_DIR}/tmhmm.txt

## RNAMMER Output Files
RNAMMER_GFF=${RNAMMER_OUTPUT_DIR}/rnammer.gff
RNAMMER_FA=${RNAMMER_OUTPUT_DIR}/rnammer.fa

## TRNASCAN Output Files
TRNASCAN_GFF3=${TRNASCAN_OUTPUT_DIR}/trnascan.gff3

## ATTRIBUTOR Output Files
ATTRIBUTOR_FASTA_NAME=attributor_fasta
ATTRIBUTOR_GFF3_NAME=attributor_gff
ATTRIBUTOR_GFF3=${ATTRIBUTOR_OUTPUT_DIR}/attributor_gff.gff3

## Set Environment Variables
export LD_LIBRARY_PATH=${PYTHON_DIR}/lib/:${PACKAGE_BIN_DIR}/gcc/lib64
export PYTHONPATH=${USER_PACKAGE_DIR}/biocode/lib:$PYTHONPATH

## Create all output directories
mkdir -p ${PROJECT_DIR} ${RNAMMER_OUTPUT_DIR} ${TRNASCAN_OUTPUT_DIR} ${TMHMM_OUTPUT_DIR} ${HMMSCAN_OUTPUT_DIR} ${RAPSEARCH_OUTPUT_DIR} ${ATTRIBUTOR_OUTPUT_DIR}

## Set threads
THREADS=10
#--------------------------------------------------------
# 1. RNAMMER - IDENTIFY rRNA's
#--------------------------------------------------------
echo -e " 
${PERL_BIN_DIR}/perl ${RNAMMER_DIR}/rnammer ${REFERENCE_FASTA} -S euk -gff ${RNAMMER_GFF} -f ${RNAMMER_FA} 
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N rnammer -l mem_free=30G -wd ${RNAMMER_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q


#--------------------------------------------------------
# 2. tRNASCAN - IDENTIFY tRNA's
#--------------------------------------------------------

echo -e " 
${TRNASCAN_BIN_DIR}/tRNAscan-SE --thread ${THREADS} -E -o ${TRNASCAN_OUTPUT_DIR}/output.txt -b trnascan.bed -a trnascan.fa -l trnascan.log ${REFERENCE_FASTA} 
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N trnascan -l mem_free=30G -wd ${TRNASCAN_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q

## Convert trnascan bed file to a gff3 file

## Need to use python 2.7 for this specific tool, reset to python 3.8 afterwards
module switch python/2.7

echo -e " 
${GFFTOOLS_DIR}/bed_to_gff_conv.py ${TRNASCAN_OUTPUT_DIR}/trnascan.bed > ${TRNASCAN_GFF3}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N conv_bed_to_gff3 -l mem_free=30G -wd ${TRNASCAN_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q -hold_jid trnascan

module switch python/3.8

#--------------------------------------------------------
# 3. TMHMM - Predict transmembrane helices in protein
#--------------------------------------------------------

echo -e " 
${TMHMM_BIN_DIR}/tmhmm ${AUGUSTUS_HINTS_AA}  > ${TMHMM_FILE}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N tmhmm -l mem_free=30G -wd ${TMHMM_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q
 
 
#--------------------------------------------------------
# 4. HMM_UNIPROT - Search Uniprot db using HMM's to identify orthologs
#--------------------------------------------------------

echo -e " 
${HMM_BIN_DIR}/hmmscan --cpu ${THREADS} --acc --cut_ga ${CODING_HMM_LIB}.lib.bin ${AUGUSTUS_HINTS_AA} > ${HMMSCAN_FILE}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N hmm_scan -l mem_free=50G -wd ${HMMSCAN_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q


ERGARTIS_DIR=/usr/local/projects
echo -e " 
${ERGARTIS_DIR}/ergatis/package-nightly/bin/hmmpfam32htab --donotrun=0 --input_file=${HMMSCAN_FILE} \
    --output_htab=${HMMSCAN_HTAB} --mldbm_file='${CODING_HMM_LIB}.lib.db' 
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N htab -l mem_free=30G -wd ${HMMSCAN_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q -hold_jid hmm_scan



#--------------------------------------------------------
# 5. RAPSEARCH2 - Search Uniprot db using rapsearch 
#--------------------------------------------------------

echo -e " 
${RAPSEARCH2_BIN_DIR}/rapsearch -z ${THREADS} -q ${AUGUSTUS_HINTS_AA} -d ${UNIREF100}.rapsearch2.db \
    -o ${RAPSEARCH_FILE} -s F -e 1e-5 -l 0 -b 10 -v 10 -s F 
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N rapsearch -l mem_free=30G -wd ${RAPSEARCH_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q

#--------------------------------------------------------
# 6. ATTRIBUTOR CONFIG FILE - Generate Config File for Attributor
#--------------------------------------------------------

echo -e "
sh ${SCRIPT_DIR}/create_attributor_config_file.sh -p ${AUGUSTUS_HINTS_AA} -g ${AUGUSTUS_HINTS_GFF3} \
    -s ${HMMSCAN_HTAB} -r ${RAPSEARCH_FILE}.m8 -t ${TMHMM_FILE} -o ${ATTRIBUTOR_CONFIG_FILE} -l ${CODING_HMM_LIB} -u ${UNIREF100}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N create_config_file -l mem_free=5G -wd ${ATTRIBUTOR_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q \
    -hold_jid rnammer,trnascan,cov_bed_to_gff3,tmhmm,hmm_scan,htab,rapsearch
#--------------------------------------------------------
# 7. ATTRIBUTOR - Assign functional Annotation based on Database Searches 
#--------------------------------------------------------
## Specifying python 3.8 for attributor
module switch python/3.8


echo -e " 
${ATTRIBUTOR_DIR}/assign_functional_annotation.py  -c ${ATTRIBUTOR_CONFIG_FILE} -o ${ATTRIBUTOR_OUTPUT_DIR}/${ATTRIBUTOR_FASTA_NAME} -f fasta >> ${ATTRIBUTOR_OUTPUT_DIR}/${ATTRIBUTOR_FASTA_NAME}.log 
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N ${ATTRIBUTOR_FASTA_NAME} -l mem_free=70G -wd ${ATTRIBUTOR_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q \
    -hold_jid rnammer,trnascan,cov_bed_to_gff3,tmhmm,hmm_scan,htab,rapsearch,create_config_file


echo -e " 
${ATTRIBUTOR_DIR}/assign_functional_annotation.py -c ${ATTRIBUTOR_CONFIG_FILE} -o ${ATTRIBUTOR_OUTPUT_DIR}/${ATTRIBUTOR_GFF3_NAME} -f gff3 >> ${ATTRIBUTOR_OUTPUT_DIR}/${ATTRIBUTOR_GFF3_NAME}.log 
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N ${ATTRIBUTOR_GFF3_NAME} -l mem_free=50G -wd ${ATTRIBUTOR_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q \
    -hold_jid rnammer,trnascan,cov_bed_to_gff3,tmhmm,hmm_scan,htab,rapsearch,create_config_file

#--------------------------------------------------------
# 8. Combine GFF3 Files - Combines the tRNAScan, RNAmmer, and attributor gff3 file into one
#--------------------------------------------------------
echo -e " 
cat ${ATTRIBUTOR_GFF3} ${RNAMMER_GFF} ${TRNASCAN_GFF3} > ${FINAL_GFF3_FILE}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N combine_gff_files -l mem_free=50G -wd ${PROJECT_DIR} -pe thread ${THREADS} -q threaded.q -hold_jid ${ATTRIBUTOR_GFF3_NAME},${ATTRIBUTOR_FASTA_NAME}


#--------------------------------------------------------
# END OF SCRIPT
#--------------------------------------------------------
