#!/bin/bash

#--------------------------------------------------------
# Set Usage Parametres
#--------------------------------------------------------
## This script will create the config file needed
## for Biocode's attributor functional annotation
## scripts to work



programmename=$0
function usage {
    
    echo "
    $(basename $programmename) [-p <polypeptide_fasta>] [-g <braker_gff3>] [-s <hmmscan_htab>] [-r <rapsearch_m8>] [-t <tmhmm_txt>] [-o <config_file>] [-l <coding_hmm_lib>] [-u <uniref100>] [-h <help>]
    
    This script will generate the config file needed for Attributor to run with the correct file paths


    ** FULL PATHS ARE REQUIRED FOR EACH ARGUMENT **

    "


}




if [ "$#" -eq 0 ]; then
    usage
fi
if [ "$#" -eq 1 ] && [ "$1" == "-h" ]; then
    usage
fi


while getopts p:g:s:r:t:o:l:u:h flag
do
    case "${flag}" in
        p) POLYPEPTIDE_FASTA=${OPTARG};;
        g) BRAKER_GFF3=${OPTARG};;
        s) HMMSCAN_HTAB=${OPTARG};;
        r) RAPSEARCH_M8=${OPTARG};;
        t) TMHMM_TXT=${OPTARG};;
        o) CONFIG_FILE=${OPTARG};;
        l) CODING_HMM_LIB=${OPTARG};;
        u) UNIREF100=${OPTARG};;
    esac
done

echo -e "---
#Config file for IGS Prokaryotic Functional Annotation
#According to guidlines here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3111993/pdf/sigs.1223234.pdf

general:
   # All proteins will start out with this product name, before evidence
   #  is evaluated
   default_product_name: hypothetical protein

   # If Yes, each annotation attribute is allowed from a different source.  For
   #  example, if there is a high-quality HMM hit that provides a gene product
   #  name and GO terms but lacks an EC number, another source (such as a BLAST
   #  hit) may provide the EC number.
   allow_attributes_from_multiple_sources: No

   # This is purely for development purposes and should usually best set to 0.  If
   #  you set this to any other integer, it will limit the number of polypeptides
   #  for which results are considered.  This allows for more rapid testing of
   #  the parser on larger datasets.
   debugging_polypeptide_limit: 0

indexes:
   coding_hmm_lib: ${CODING_HMM_LIB}_lib.sqlite3
   uniref100: ${UNIREF100}.sqlite3

input:
   # These are the files on which annotation will be applied.  At a minimum, the 'fasta'
   #  option must be defined.  GFF3 output cannot be specified unless 'gff3' input is
   #  also provided.
   # These test data are for SRS15430
   polypeptide_fasta: ${POLYPEPTIDE_FASTA}
   gff3: ${BRAKER_GFF3}
order:
   - coding_hmm_lib__equivalog
   - rapsearch2__uniref100__trusted_full_full
   - coding_hmm_lib__equivalog_domain
   - rapsearch2__uniref100__trusted_partial_full
   - coding_hmm_lib__subfamily
   - coding_hmm_lib__superfamily
   - coding_hmm_lib__subfamily_domain
   - coding_hmm_lib__domain
   - coding_hmm_lib__pfam
   - rapsearch2__uniref100__trusted_full_partial
   - tmhmm
   - coding_hmm_lib__hypothetical_equivalog

evidence:
   - label: coding_hmm_lib__equivalog
     type: HMMer3_htab
     path: ${HMMSCAN_HTAB}
     class: equivalog
     index: coding_hmm_lib

   - label: coding_hmm_lib__equivalog_domain
     type: HMMer3_htab
     path: ${HMMSCAN_HTAB}
     class: equivalog_domain
     index: coding_hmm_lib

   - label: coding_hmm_lib__subfamily
     type: HMMer3_htab
     path: ${HMMSCAN_HTAB}
     class: subfamily
     index: coding_hmm_lib
     append_text: family protein

   - label: coding_hmm_lib__superfamily
     type: HMMer3_htab
     path: ${HMMSCAN_HTAB}
     class: superfamily
     index: coding_hmm_lib
     append_text: family protein

   - label: coding_hmm_lib__subfamily_domain
     type: HMMer3_htab
     path: ${HMMSCAN_HTAB}
     class: subfamily_domain
     index: coding_hmm_lib
     append_text: domain protein

   - label: coding_hmm_lib__domain
     type: HMMer3_htab
     path: ${HMMSCAN_HTAB}
     class: domain
     index: coding_hmm_lib
     append_text: domain protein

   - label: coding_hmm_lib__pfam
     type: HMMer3_htab
     path: ${HMMSCAN_HTAB}
     class: pfam
     index: coding_hmm_lib
     append_text: family protein

   - label: coding_hmm_lib__hypothetical_equivalog
     type: HMMer3_htab
     path: ${HMMSCAN_HTAB}
     class: hypoth_equivalog
     index: coding_hmm_lib

   - label: rapsearch2__uniref100__trusted_full_full
     type: RAPSearch2_m8
     path: ${RAPSEARCH_M8}
     class: trusted
     index: uniref100
     query_cov: 80%
     match_cov: 80%
     percent_identity_cutoff: 50%

   - label: rapsearch2__uniref100__trusted_partial_full
     type: RAPSearch2_m8
     path: ${RAPSEARCH_M8}
     class: trusted
     index: uniref100
     match_cov: 80%
     percent_identity_cutoff: 50%
     append_text: domain protein

   - label: rapsearch2__uniref100__trusted_full_partial
     type: RAPSearch2_m8
     path: ${RAPSEARCH_M8}
     class: trusted
     index: uniref100
     query_cov: 80%
     percent_identity_cutoff: 50%
     append_text: domain protein

   - label: rapsearch2__uniref100__all_full_full
     type: RAPSearch2_m8
     path: ${RAPSEARCH_M8}
     index: uniref100
     query_cov: 80%
     match_cov: 80%
     percent_identity_cutoff: 50%
     prepend_text: putative

   - label: rapsearch2__uniref100__all_partial_full
     type: RAPSearch2_m8
     path: ${RAPSEARCH_M8}
     index: uniref100
     match_cov: 80%
     percent_identity_cutoff: 50%
     prepend_text: putative
     append_text: domain protein

   - label: rapsearch2__uniref100__all_full_partial
     type: RAPSearch2_m8
     path: ${RAPSEARCH_M8}
     index: uniref100
     query_cov: 80%
     percent_identity_cutoff: 50%
     prepend_text: putative
     append_text: domain protein

   - label: tmhmm
     type: TMHMM
     # this is the product name that will be assigned given a positive TMHMM match
     product_name: putative integral membrane protein
     # minimum required predicted helical spans across the membrane required to apply evidence
     min_helical_spans: 5
     path: ${TMHMM_TXT} " > ${CONFIG_FILE}
