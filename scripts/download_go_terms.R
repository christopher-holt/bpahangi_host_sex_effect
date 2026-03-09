#!/usr/local/packages/r-4.1.2/bin/Rscript

#-------------------------------------------------------
## This script will download a file that contains
## GO terms and the associated definition. It will
## be required when generating a geneinfo file
##
## How to run: ./download_go_terms.R --output /path/to/gomap/file
##
## usage: ./download_go_terms.R [-h] --output OUTPUT
##
## optional arguments:
##  -h, --help       show this help message and exit
## 
## required arguments:
##  --output OUTPUT  Output file name for go id map file
#-------------------------------------------------------


## Import libraries
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(GO.db))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(argparse))

## Define arguments
parser <- ArgumentParser()
parser$add_argument("--output", help = "Output file name for go id map file")
args <- parser$parse_args()

## Create a map file
### list of go terms and descriptions
goterms.list <- as.list(GOTERM)
gomap <- tibble::tibble('goterm' = as.character(), 'gocatergory' = as.character(), 'godescription' = as.character())
### Assign go terms their catergory (biological process, molecular function, cellular component, or universal)  and functional description
for(goterm in 1:length(goterms.list)){
  gomap <- gomap %>% tibble::add_row(goterm = as.character(AnnotationDbi::GOID(goterms.list[[goterm]])),
                             gocatergory = as.character(AnnotationDbi::Ontology(goterms.list[[goterm]])),
                             godescription = as.character(AnnotationDbi::Term(goterms.list[[goterm]])))
  
}

## Format catergory
gomap <- gomap %>% dplyr::mutate(gocatergory = ifelse(grepl("BP", gomap$gocatergory), "go_biologicalprocess",
                                               ifelse(grepl("MF", gomap$gocatergory), "go_molecularfunction",
                                                      ifelse(grepl("CC", gomap$gocatergory), "go_cellularcomponent", "go_universal"))))


## Write out to a text file
readr::write_delim(gomap, args$output,delim = "\t", quote = "none")
