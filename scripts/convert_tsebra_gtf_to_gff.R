#!/usr/local/packages/r-4.1.2/bin/Rscript
#-------------------------------------------------------
## This script will take the GTF output from 
## BRAKER's TSEBRA function and convert it
## to a GFF file
##
## How to run: ./convert_tsebra_gtf_to_gff.R --gtf /path/to/gtf/from/tsebra --gff /path/to/new/gff
##
## usage: ./convert_tsebra_gtf_to_gff.R [-h] --gtf GTF --gff GFF
##
## optional arguments:
##   -h, --help  show this help message and exit
##
## required arguments:
##   --gtf GTF   File path to the GTF file from TSEBRA
##   --gff GFF   File path to the output GFF3 formatted file
#-------------------------------------------------------


library(tibble)
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(data.table)
library(argparse)

parser = argparse::ArgumentParser()
parser$add_argument("--gtf", help = "File path to the GTF file from TSEBRA, Required")
parser$add_argument("--gff", help = "File path to the output GFF3 formatted file, Required")

args <- parser$parse_args()

## Read in TSEBRA file
gtf.path <- args$gtf 
gtf <- readr::read_delim(gtf.path, delim = "\t", col_names = F)
gtf$X9 <- paste0("ID=", gtf$X9)
gtf$X9 <- gsub("\"", "", gtf$X9)
gtf$X9 <- gsub("transcript_id ", "", gtf$X9)
gtf$X9 <- gsub("; gene_id ", ";gene_id=", gtf$X9)
gtf <- gtf %>% separate(X9, into = c("transcript", 'gene'), sep = ';')
gtf <- gtf %>% mutate(gene = ifelse(.$X3 == 'gene', paste0(.$transcript), NA))
gtf <- gtf %>% fill(gene)

gtf <- gtf %>% mutate(X1 = ifelse(grepl("_pilon_LN",.$X1), str_remove(.$X1, "_LN.*"), .$X1))


format_feature_name <- function(feature){
  feature$transcript <- feature$transcript_new
  feature <- feature %>% dplyr::mutate(gene_new = ifelse(.$X3 %in% c("gene", "transcript"), .$gene, .$transcript))
  feature$gene <- feature$gene_new
  feature$transcript_new <- NULL
  feature$gene_new <- NULL
  
  return(feature)
  
}



genes <- unique(gtf$gene)
gff <- data.frame()
for(i in 1:length(genes)){
  gene_gtf <- gtf %>% filter(gene == genes[i])
  if(nrow(gene_gtf[gene_gtf$X3 == "transcript",]) == 1){
    gene_gtf <- gene_gtf %>% dplyr::mutate(transcript_new = ifelse(.$X3 != "gene", paste0(gene_gtf$gene, '.', purrr::map(strsplit(gene_gtf$transcript, split = "\\."), 3), ';'), 
                                                            ifelse(.$X3 == "gene", paste0(gene_gtf$gene, ';'), NA)))
    
    gene_gtf <- format_feature_name(gene_gtf)
    
    gff <- rbind(gff, gene_gtf)
  }else{
    transcripts <- unique(gene_gtf[gene_gtf$X3 == "transcript",]$transcript)
    final_transcript_gtf <- data.frame()
    final_transcript_gtf <- rbind(final_transcript_gtf, gene_gtf[gene_gtf$X3 == "gene",])
    for(j in 1:length(transcripts)){
      transcript_gtf <- gene_gtf %>% dplyr::filter(transcript == transcripts[j])
      transcript_gtf <- transcript_gtf %>% dplyr::mutate(transcript_new = paste0(.$gene, ".t", j))
      
      transcript_gtf <- format_feature_name(transcript_gtf)

      final_transcript_gtf <- rbind(final_transcript_gtf, transcript_gtf)
    }

    gff <- rbind(gff, final_transcript_gtf)
    
  }
}

gff$transcript <- stringr::str_remove(gff$transcript, ';')
gff$gene <- stringr::str_remove(gff$gene, ';')

gff <- gff %>% dplyr::mutate(transcript = ifelse(.$X3 %in% c("CDS", "exon", "intron"), paste0(.$transcript, '.', .$X3),
                                   ifelse(.$X3 %in% c("start_codon"), paste0(.$transcript, '.start1'), 
                                          ifelse(.$X3 %in% c("stop_codon"), paste0(.$transcript, '.stop1'), .$transcript))))

gff <- tibble::rownames_to_column(gff, "row_num")
gff_to_format <- gff %>% dplyr::filter(X3 %in% c('CDS', 'intron', 'exon') )
gff_remain <- gff[!(gff$row_num %in% gff_to_format$row_num),]

gff_to_format <- data.table::data.table(gff_to_format)
gff_to_format <- gff_to_format[, id := seq_len(.N), by = transcript]
gff_to_format <- gff_to_format[, id := rowid(transcript)]

gff_to_format$transcript <- paste0(gff_to_format$transcript, gff_to_format$id)
gff_to_format$id <- NULL

gff <- rbind(gff_remain, gff_to_format)
gff$row_num <- as.double(gff$row_num)
gff <- gff[order(gff$row_num),]
gff$row_num <- NULL


gff$transcript <- paste0(gff$transcript, ';')
gff$gene <- paste0(gff$gene, ';')


gff$gene <- gsub("ID=", "Parent=", gff$gene)
gff <- gff %>% dplyr::mutate(X9 = ifelse(.$X3 != 'gene', paste0(.$transcript, .$gene), .$transcript))
gff$transcript <- NULL
gff$gene <- NULL

gff$X3 <- stringr::str_replace(gff$X3, "transcript", "mRNA")


readr::write_delim(x = gff, 
            file = args$gff,
            delim = '\t', col_names = F, quote = 'none', escape = "backslash")

