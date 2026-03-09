library(tidyverse)
## Convert BmID to WBGene

ANNOTATION_URL <- "b_malayi.PRJNA10729.WS276.annotations.wormbase.gff3"
annotations <- readr::read_delim(ANNOTATION_URL, "\t", escape_double = FALSE, col_names = FALSE,
                                 trim_ws = TRUE, show_col_types = F)

annotations <- annotations %>% separate('X9', into = c("Name", "X9"), sep = ';')
annotations$X9 <- gsub(".*:", "", annotations$X9)
annotations$X9 <- gsub(".*=", "", annotations$X9)
annotations$Name <- gsub(".*:", "", annotations$Name)
annotations$Name <- gsub(".*=", "", annotations$Name)


Bm_to_WB <- annotations %>% filter(X3 == "mRNA") %>% select("Name", "X9")

write_delim(Bm_to_WB, "Bm_to_WB.tsv", delim = "\t",
            col_names = F)


interproscan <- readr::read_delim("b_malayi.PRJNA10729.WS276.genomic.aa.formatted.tsv",
                                  col_names = seq(1,15),
                                  delim = "\t",
                                  show_col_types = FALSE)


interproscan <- interproscan %>% dplyr::inner_join(Bm_to_WB, by = c("X1" = "Name")) %>% select("X9.y", everything()) %>% select(-"X1")

write_delim(interproscan, "b_malayi.PRJNA10729.WS276.genomic.aa.formatted.names_updated.tsv", delim = "\t",
            col_names = F)
