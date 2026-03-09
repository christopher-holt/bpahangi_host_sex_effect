#!/usr/local/bin/Rscript
#-------------------------------------------------------
## This script will download the worm recovery table
## and generate a presence/absence table for different
## life stages of B. pahangi isolated per SQ infected
## gerbil (S1 Fig)
#-------------------------------------------------------


library(tidyverse)
library(tidyr)
library(ggpubr)
library(pwr)
devtools::source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

## Set paths
cwd <- getwd()
OUTPUT_DIR <- paste0(cwd, "/output_files/")
FIGURE_DIR <- paste0(OUTPUT_DIR, '/figures/')

## Tidy Data
WORM_RECOVERY_URL <- "https://raw.githubusercontent.com/christopher-holt/bpahangi_host_sex_effect/refs/heads/main/input_files/Worm_Recovery_Table.csv"

df <- read_csv(WORM_RECOVERY_URL)


df$`# Male` <- as.double(df$`# Male`)
df$`# Female` <- as.double(df$`# Female`)
df$`# Damaged` <- as.double(df$`# Damaged`)
df$Total <- as.double(df$Total)
df$`Total # Males + Females` <- as.double(df$`Total # Males + Females`)

## List of five SQ samples that should be removed
df[!(grepl("dpi",df$`Days Post Infection at Necropsy`)),]
## Remove the five SQ Samples, leaving total of 63 Samples
df <- df[!(df$`Gerbil ID` %in% c("M058_SQ_Bp", "M059_SQ_Bp", "F014_SQ_Bp", "M023_SQ_Bp", "F027_SQ_Bp")),]
## Remove Poorly Infected IP Sample, leaving total of 62
#df <- df[!(df$`Gerbil ID` %in% c("M201_IP_Bp")),]
df[df$`Gerbil ID` == "M201_IP_Bp",]$Total <- 3
df[df$`Gerbil ID` == "M201_IP_Bp",]$`Total # Males + Females` <- 3


## Number of SQ with MF but not adult worms
df %>% filter(`Method of Exposure` == "Subcutaneous injection" & Total == 0 & `mf/20 μL` > 0)

## SQ with all life stages
as.data.frame(table(df %>% filter(`Method of Exposure` == "Subcutaneous injection" & `mf/20 μL` > 0 & `# Male` > 0 & `# Female` > 0) %>% select('Sex')))


SQ <- df[df$`Method of Exposure` == "Subcutaneous injection",]
SQ <- SQ[!SQ$`Gerbil ID` %in% c("FBpSQ014", "FBpSQ023", "FBpSQ027"),]

SQ <- SQ %>% select('Gerbil ID', "Sex","# Male", "# Female", "mf/20 μL", "Litter")

SQ.final <- data.frame()
## Make SQ tidy
SQ.F.female <- SQ %>% filter(Sex == "F") %>% select(-'# Male', -"mf/20 μL")
SQ.F.female$wormstage <- 'female_worm'
colnames(SQ.F.female)[3] <- 'counts'
SQ.final <- rbind(SQ.final, SQ.F.female)

SQ.F.male <- SQ %>% filter(Sex == "F") %>% select(-'# Female', -"mf/20 μL")
SQ.F.male$wormstage <- 'male_worm'
colnames(SQ.F.male)[3] <- 'counts'
SQ.final <- rbind(SQ.final, SQ.F.male)


SQ.F.mf <- SQ %>% filter(Sex == "F") %>% select(-'# Female', -"# Male")
SQ.F.mf$wormstage <- 'mf'
colnames(SQ.F.mf)[3] <- 'counts'
SQ.final <- rbind(SQ.final, SQ.F.mf)

SQ.M.female <- SQ %>% filter(Sex == "M") %>% select(-'# Male', -"mf/20 μL")
SQ.M.female$wormstage <- 'female_worm'
colnames(SQ.M.female)[3] <- 'counts'
SQ.final <- rbind(SQ.final, SQ.M.female)

SQ.M.male <- SQ %>% filter(Sex == "M") %>% select(-'# Female', -"mf/20 μL")
SQ.M.male$wormstage <- 'male_worm'
colnames(SQ.M.male)[3] <- 'counts'
SQ.final <- rbind(SQ.final, SQ.M.male)

SQ.M.mf <- SQ %>% filter(Sex == "M") %>% select(-'# Female', -"# Male")
SQ.M.mf$wormstage <- 'mf'
colnames(SQ.M.mf)[3] <- 'counts'
SQ.final <- rbind(SQ.final, SQ.M.mf)

SQ.final <- SQ.final %>% mutate(prescence_of_worms = ifelse(counts > 0, "present", "not_present"))

SQ.final <- SQ.final[order(SQ.final$counts, SQ.final$wormstage),]
SQ.final$`Gerbil ID` <- as.factor(SQ.final$`Gerbil ID`)
SQ.final$prescence_of_worms <- as.factor(SQ.final$prescence_of_worms)
SQ.final$wormstage <- as.factor(SQ.final$wormstage)

SQ.final <- SQ.final %>% mutate(binary = ifelse(prescence_of_worms == "not_present", 0, 1))

SQ.final <- SQ.final %>% arrange(-binary, wormstage)

pdf(paste0(FIGURE_DIR, "/SQ_worm_recovery_heatmap.pdf"),
    width = 10,
    height = 10)
  print(
  ggplot(SQ.final, aes(x = `Gerbil ID`, y = wormstage, fill = prescence_of_worms) ) + 
    geom_tile(colour = "black") + theme_bw()  + theme(axis.text.x = element_text(angle = 90), line = element_blank()) +
    scale_fill_manual(breaks = c("not_present", "present"),
                      values = c("white", "salmon")) +
    labs(title = "Presence of worms in SQ-injected Gerbils") + coord_flip()
  )
dev.off()










