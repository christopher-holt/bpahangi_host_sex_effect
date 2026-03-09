#!/usr/local/packages/r-4.1.2/bin/Rscript
#-------------------------------------------------------
## This script reads:
##      formatted interproscan tsv
##      htseq counts file (or file of gene names)
##      go map file
##      ipr map file
## And will write out:
##      geneinfo file
##
## How to run: ./interproscan_to_geneinfo.R --interproscan /path/to/formatted/interproscan/tsv  \
##                                          --counts /path/to/counts/file \
##                                          --gomap /path/to/gomap/file \
##                                          --iprmap /path/to/iprmap/file \
##                                          --out /path/to/output/geneinfo/file
##
## usage: ./interproscan_to_geneinfo.R [-h] --interproscan INTERPROSCAN
##                                     --counts COUNTS --gomap GOMAP
##                                     --iprmap IPRMAP --out OUT
##
## optional arguments:
##   -h, --help            show this help message and exit
## 
## required arguments:
##   --interproscan INTERPROSCAN
##                         File path to the tsv formatted interproscan file
##   --counts COUNTS       File path to the output of HTSeq Counts file or list of gene names
##   --gomap GOMAP         File path to the goid map file
##   --iprmap IPRMAP       File path to the ipr map file
##   --out OUT             File path to write out the geneinfo file
##
#-------------------------------------------------------

library(readr)
library(argparse)
library(stringr)
doesFileExist <- function(file){
  ## Check to see if a file exists, if not then quit with status = -1
  if (file.access(file) == -1){
    message(paste0(file, " does not exist"))
    q(status=1)
  }
}

## Prevent URL scraping from stalling
options(timeout = 1000)

parser = argparse::ArgumentParser()
parser$add_argument("--interproscan", help = "File path to the tsv formatted interproscan file, Required")
parser$add_argument("--counts", help = "File path to the output of HTSeq Counts file, Required")
parser$add_argument("--gomap", help = "File path to the goid map file, Required")
parser$add_argument("--iprmap", help = "File path to the ipr map file, Required")
parser$add_argument("--out", help = "File path to write out the geneinfo file")

args <- parser$parse_args()

## Check to see if the four required files exist
doesFileExist(args$interproscan)
doesFileExist(args$counts)
doesFileExist(args$gomap)
doesFileExist(args$iprmap)

interproscan.path = args$interproscan
counts.path = args$counts
gomap.path = args$gomap
iprmap.path = args$iprmap

interproscan <- readr::read_delim(interproscan.path,
                           col_names = seq(1,15),
                           delim = "\t",
                           show_col_types = FALSE)

counts <- readr::read_delim(counts.path,
                    delim = '\t',
                    col_names = F,
                    show_col_types = FALSE)

## Remove the __ lines from HTSeq file
counts <- counts[!grepl("__", counts$X1),]

gomap <- readr::read_delim(gomap.path,
                    delim = "\t",
                    show_col_types = FALSE)

gomap$goterm <- as.character(gomap$goterm)
gomap$gocatergory <- as.character(gomap$gocatergory)
gomap$godescription <- as.character(gomap$godescription)

iprmap <- readr::read_delim(iprmap.path,
                     col_names = F,
                     delim = "\t",
                     show_col_types = FALSE)
iprmap <- iprmap[1:2]

iprmap$X1 <- as.character(iprmap$X1)
iprmap$X2 <- as.character(iprmap$X2)

geneinfo <- as.data.frame(matrix(nrow=length(unique(counts$X1)),
                                 ncol=6))
colnames(geneinfo) <- c("gene","interpro_description","go_biologicalprocess","go_cellularcomponent","go_molecularfunction", "go_universal")



goterms <- sort(unique(unlist(strsplit(paste(interproscan$X14,collapse="|"),split="[|]"))))
iprterms <- sort(unique(unlist(strsplit(paste(interproscan$X12,collapse="|"),split="[|]"))))

goterms <- goterms[!(goterms %in% c("", "-", "NA"))]
iprterms <- iprterms[!(iprterms %in% c("", "-"))]


missinggoterms <- goterms[!(goterms %in% gomap$goterm)]
missinggoterms <- missinggoterms[!(is.na(missinggoterms))]
missingiprterms <- iprterms[!(iprterms %in% iprmap$X1)]

if(length(missinggoterms) > 0){
  for(missinggoterm in missinggoterms){
    tryCatch(
      expr = {
        go <- readLines(paste0("https://gowiki.tamu.edu/wiki/index.php/Category:",missinggoterm),warn=F)
        gocategory <- grep("namespace:",go,value=T)
        gocategory <- gsub(".* ! ","",gocategory)
        gocategory <- gsub("\".*","",gocategory)
        
        godescription <- grep("wgTitle",go,value=T)
        godescription <- gsub(".* ! ","",godescription)
        godescription <- gsub("\",","",godescription)
        
        gomap <- as.data.frame(rbind(gomap,
                                     c(missinggoterm,gocategory,godescription)))
        
      },
      error = function(e){
        message(paste0(missinggoterm, " was not found"))
      })
  }
}else{
  message("No missing go terms")
}


if(length(missingiprterms) > 0){
  for(missingiprterm in missingiprterms){
    tryCatch(
      expr = {
        interprodescription <- readLines(paste0("https://www.ebi.ac.uk/interpro/entry/",missingiprterm))
        interprodescription <- grep("h2 class", interprodescription,value = T)
        interprodescription <- gsub(".*<h2 class=\"strapline\">","",interprodescription)
        interprodescription <- gsub(" <span>.*","",interprodescription)
        iprmap <- as.data.frame(rbind(iprmap,
                                      c(missingiprterm,interprodescription)))
      },
      error = function(e){
        message(paste0(missingiprterm, " was not found"))
      }
    )
  }
}else{
  message("No missing ipr terms")
}

## Terms that are missing from the geneinfo file
missing_from_map <- c()

geneinfo$gene <- unique(counts$X1)
for(i in 1:length(unique(counts$X1))){
  interproscan.subset <- interproscan[interproscan$X1 == unique(counts$X1)[i],]
  goterms <- sort(unique(unlist(strsplit(paste(interproscan.subset$X14,collapse="|"),split="[|]"))))
  iprterms <- sort(unique(unlist(strsplit(paste(interproscan.subset$X12,collapse="|"),split="[|]"))))
  
  goterms <- goterms[!(goterms %in% c("", "-", "NA"))]
  iprterms <- iprterms[!(iprterms %in% c("", "-"))]
  while(length(iprterms) > 0){
    geneinfo$interpro_description[i] <- paste0(geneinfo$interpro_description[i],"|",iprterms[1],":",iprmap[iprmap$X1 == iprterms[1],2])
    if(length(iprmap[iprmap$X1 == iprterms[1],2]) == 0){
      missing_from_map <- c(missing_from_map,iprterms[1])
    }
    iprterms <- iprterms[-1]
  }
  while(length(goterms) > 0){
    if(gomap[gomap$goterm == goterms[1],2] == "go_biologicalprocess"){
      geneinfo$go_biologicalprocess[i] <- paste0(geneinfo$go_biologicalprocess[i],"|",goterms[1],":",gomap[gomap[,1] == goterms[1],3])
    }else if(gomap[gomap$goterm == goterms[1],2] == "go_cellularcomponent"){
      geneinfo$go_cellularcomponent[i] <- paste0(geneinfo$go_cellularcomponent[i],"|",goterms[1],":",gomap[gomap[,1] == goterms[1],3])
    }else if(gomap[gomap$goterm == goterms[1],2] == "go_molecularfunction"){
      geneinfo$go_molecularfunction[i] <- paste0(geneinfo$go_molecularfunction[i],"|",goterms[1],":",gomap[gomap[,1] == goterms[1],3])
    }else if(gomap[gomap$goterm == goterms[1],2] == "go_universal"){
      geneinfo$go_universal[i] <- paste0(geneinfo$go_universal[i],"|",goterms[1],":",gomap[gomap$goterm == goterms[1],3])
    }else{
      missing_from_map <- c(missing_from_map,goterms[1])
    }
    goterms <- goterms[-1]
  }
}

missing_from_map <- unique(missing_from_map)
print(missing_from_map)

for(i in 1:ncol(geneinfo)){
  geneinfo[,i] <- gsub("NA[|]","",geneinfo[,i])
}
## Fix NA Values
geneinfo$interpro_description[which(is.na(geneinfo$interpro_description))] <- "No InterPro entry"
geneinfo$go_biologicalprocess[which(is.na(geneinfo$go_biologicalprocess))] <- "No GO terms for biological process"
geneinfo$go_cellularcomponent[which(is.na(geneinfo$go_cellularcomponent))] <- "No GO terms for cellular component"
geneinfo$go_molecularfunction[which(is.na(geneinfo$go_molecularfunction))] <- "No GO terms for molecular function"
geneinfo$go_universal[which(is.na(geneinfo$go_universal))] <- "No GO terms for universal function"

if(!is.null(args$out)){
  write.table(geneinfo,
              file = paste0(args$out),
              quote = F,
              col.names = T,
              row.names = F,
              sep = "\t")
}else{
  print(geneinfo)
}



