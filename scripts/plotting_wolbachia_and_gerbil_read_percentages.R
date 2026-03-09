#-------------------------------------------------------
## This script will read in the idxstats results
## for brugia pahangi and gerbil genomes. A plot showing
## the percentage of reads mapping to each genome
## for each sample will be generated
#-------------------------------------------------------

## Load libraries
library(tidyverse)

## Read in idxstats results for brugia pahangi genomes
idxstat_url <- "https://raw.githubusercontent.com/christopher-holt/bpahangi_host_sex_effect/refs/heads/main/input_files/combined.brugia_pahangi.idxstat.txt"
idxstats_brugia_pahangi <- readr::read_delim(idxstat_url,
                  "\t", escape_double = FALSE, col_names = FALSE, 
                  trim_ws = TRUE, show_col_types = F)


## Name columns and assign each contig to either Wolbachia (wBp), mitochondria, or nuclear genome
colnames(idxstats_brugia_pahangi) <- c("contig", "contig_length", "mapped_reads", "unmapped_reads", "sample")
idxstats_brugia_pahangi <- idxstats_brugia_pahangi %>% mutate(species = ifelse(contig == "CP050521.1", "wBp",
                                                 ifelse(contig == "CM022469.1", "mito",
                                                        ifelse(contig == "*", "unknown", "nuclear"))))


## Read on idxstats for gerbil genome
idxstat_url <- "https://raw.githubusercontent.com/christopher-holt/bpahangi_host_sex_effect/refs/heads/main/input_files/combined.gerbil.idxstat.txt"
idxstats_gerbil <- readr::read_delim(idxstat_url,
                              "\t", escape_double = FALSE, col_names = FALSE, 
                              trim_ws = TRUE, show_col_types = F)
## Name each column
colnames(idxstats_gerbil) <- c("contig", "contig_length", "mapped_reads", "unmapped_reads", "sample")
idxstats_gerbil$species <- "gerbil"


## Combining idxstats and idxstats_gerbil
idxstats_combined <- rbind(idxstats_brugia_pahangi, idxstats_gerbil)

idxstats_combined <- idxstats_combined %>%
  group_by(sample, species) %>% ## Group by sample and species (contig)
  summarise(sum_mapped_reads = sum(mapped_reads)) %>% ## Create sum of total number reads for each contig type (nuclear, mito, wBp or gerbil) for each sample
  filter(species != "unknown") %>% ## remove contigs classified as unknown
  ungroup() %>%
  group_by(sample) %>% ## regroup by sample name only
  pivot_wider(names_from = species,
              values_from = sum_mapped_reads) %>% ## create a wide table with sample name in col 1, col 2 is # gerbil reads, col 3 is # mito reads ...
  group_by(sample) %>%
  mutate(total_reads = mito + nuclear + wBp) %>% ## new column with total numbers of all reads
  mutate(total_reads_including_gerbil = total_reads + gerbil) %>% ## total number of reads including the gerbil
  ## Next 7 lines are just calculating percentages for each read type based on total # of reads
  mutate(gerbil_percentage = (gerbil/total_reads_including_gerbil)*100) %>%
  mutate(mito_gerbil_percentage = (mito/total_reads_including_gerbil)*100) %>%
  mutate(nuclear_gerbil_percentage = (nuclear/total_reads_including_gerbil)*100) %>%
  mutate(wBp_gerbil_percentage = (wBp/total_reads_including_gerbil)*100) %>%
  mutate(mito_percentage = (mito/total_reads)*100) %>%
  mutate(nuclear_percentage = (nuclear/total_reads)*100) %>%
  mutate(wBp_percentage = (wBp/total_reads)*100) %>%
  ## Select specific columns
  select("sample", "gerbil_percentage","mito_gerbil_percentage","nuclear_gerbil_percentage","wBp_gerbil_percentage",
         "mito_percentage","nuclear_percentage","wBp_percentage") %>% 
  ## Make into a longer table for plotting purposes
  pivot_longer(cols = c("gerbil_percentage","mito_gerbil_percentage","nuclear_gerbil_percentage","wBp_gerbil_percentage",
                        "mito_percentage","nuclear_percentage","wBp_percentage"),
               names_to = "variable",
               values_to = "values") %>%
  ungroup() %>%
  mutate(contains_gerbil = ifelse(grepl("gerbil", .$variable), "contains_gerbil", "no_gerbil"))

## Separate into percentages with and without gerbil, 
idxstats_combined_gerbil <- idxstats_combined %>% filter(contains_gerbil == "contains_gerbil")
idxstats_combined_no_gerbil <- idxstats_combined %>% filter(contains_gerbil == "no_gerbil")

## Reorder so samples are in decreasing order of gerbil percentage
sample_order <- idxstats_combined_gerbil %>%
  filter(variable == "gerbil_percentage") %>%
  arrange(desc(values)) %>%
  pull(sample)

## Refactor
idxstats_combined_gerbil$sample <- factor(idxstats_combined_gerbil$sample, levels = sample_order)

## Plot of read percentages including the gerbil
ggplot(idxstats_combined_gerbil, aes(x = sample, y = values, fill = variable)) +
  geom_bar(stat = 'identity', position = "stack") + 
  theme_bw() +   theme(axis.text.x = element_text(angle = 90)) + labs(x = "Sample")


## Redo but for only Bphangi genomes, no gerbil
sample_order <- idxstats_combined_no_gerbil %>%
  filter(variable == "nuclear_percentage") %>%
  arrange(values) %>%
  pull(sample)

idxstats_combined_no_gerbil$sample <- factor(idxstats_combined_no_gerbil$sample, levels = sample_order)

ggplot(idxstats_combined_no_gerbil, aes(x = sample, y = values, fill = variable)) +
  geom_bar(stat = 'identity', position = "stack") + 
  theme_bw() +   theme(axis.text.x = element_text(angle = 90)) + labs(x = "Sample")





















