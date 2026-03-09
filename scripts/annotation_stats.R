#-------------------------------------------------------
## This script reads in a GFF file and will print out:
##      The number of genes with multiple transcripts
##      Total number of genes
##      Average number of CDS per gene
##      Number of inferred exons per gene
##      Average CDS length
##      Average exon length
#-------------------------------------------------------

## Loading in libraries
library(readr)
library(tidyverse)

format_gff_file <- function(annotations){
  annotations <- annotations %>% separate('X9', into = c("Name", "X9"), sep = ';')
  annotations$X9 <- gsub(".*:", "", annotations$X9)
  annotations$X9 <- gsub(".*=", "", annotations$X9)
  annotations$Name <- gsub(".*:", "", annotations$Name)
  annotations$Name <- gsub(".*=", "", annotations$Name)
  
  return(annotations)
}

num_rows <- function(annotations, feature){
  return(annotations %>% filter(X3 == paste0(feature)) %>% nrow(.))
}

num_feature_per_gene <- function(annotations, feature){
  (num_rows(annotations, paste0(feature)))/num_rows(annotations, 'gene')
}

average_length <- function(annotations,feature){
  return((annotations %>% filter(X3 == paste0(feature)) %>% mutate(length = X5 - X4) %>% 
            summarise(x = sum(length)))/(annotations %>% filter(X3 == paste0(feature)) %>% nrow(.)))
}

summary_statistics <- function(annotations){
  message('Number of Genes with Multiple Transcripts: ', as.data.frame(table(annotations[annotations$X3 == "mRNA",]$X10)) %>% filter(Freq >1) %>% nrow(.))
  message('total number of genes: ', num_rows(annotations, 'gene'))
  message('Avg Number CDS/gene: ', num_feature_per_gene(annotations, 'CDS'))
  message('Number of Inferred Introns per Gene: ', num_feature_per_gene(annotations, "exon") - num_feature_per_gene(annotations, 'mRNA'))
  message('Avg CDS Length: ', average_length(annotations, 'CDS'))
  ## Avg Exon size
  message('Avg Exon Length: ', average_length(annotations, 'exon'))
  
}


## Bpahangi from NCBI
ANNOTATION_FILE <- "GCA_012070555.1_ASM1207055v1_genomic.no_mitochondria.gff"
## Bmalayi from NCBI
ANNOTATION_FILE <- "GCF_000002995.4_B_malayi-4.0_genomic.gff"



## Read in Annotation File
annotations <- readr::read_delim(ANNOTATION_FILE, "\t", escape_double = FALSE, col_names = FALSE,
                                 trim_ws = TRUE, show_col_types = F, comment = "#")
annotations <- format_gff_file(annotations)


mRNA_features <- annotations %>% filter(X3 == "mRNA")
mRNA_features <- unique(mRNA_features$X9)

annotations <- annotations %>% mutate(X10 = ifelse(.$X3 == "gene", .$Name, NA))
annotations <- tidyr::fill(annotations, X10)
annotations <- annotations[annotations$X10 %in% mRNA_features,]

## Bpahangi
annotations_on_chromosomes <- annotations %>% filter(X1 %in% c("JAAVKF010000001.1", "JAAVKF010000002.1", "JAAVKF010000003.1", "JAAVKF010000004.1",
                                                "JAAVKF010000005.1", "JAAVKF010000006.1",
                                                "JAAVKF010000007.1", "JAAVKF010000008.1",
                                                "CM022444.1",
                                                "JAAVKF010000010.1", "JAAVKF010000011.1", "JAAVKF010000012.1"))


## Bmalayi
## Chr 1: NW_025062475.1
## Chr 3: NW_025062477.1
## Chr 2: NW_025062476.1
## Chr 4: NW_025062478.1
## Chr X: NW_025062479.1
annotations_on_chromosomes <- annotations %>% filter(X1 %in% c("NW_025062475.1",
                                                               "NW_025062477.1",
                                                               "NW_025062476.1",
                                                               "NW_025062478.1",
                                                               "NW_025062479.1"))


summary_statistics(annotations)
summary_statistics(annotations_on_chromosomes)

