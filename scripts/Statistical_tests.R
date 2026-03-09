#-------------------------------------------------------
## This script will read in the worm recovery table
## and perform statistical tests on a variety
## of variables as well as a power analysis
#-------------------------------------------------------
## Load in libraries
library(tidyverse)
library(tidyr)
library(ggpubr)
library(pwr)

## Assess number of worms injected, Brugia pahangi
## Tidy Data
WORM_RECOVERY_URL <- "https://raw.githubusercontent.com/christopher-holt/bpahangi_host_sex_effect/refs/heads/main/input_files/Worm_Recovery_Table.csv"

df <- read_csv(WORM_RECOVERY_URL)

## Format data types
df$`# Male` <- as.double(df$`# Male`)
df$`# Female` <- as.double(df$`# Female`)
df$`# Damaged` <- as.double(df$`# Damaged`)
df$`mf/20 μL` <- as.double(df$`mf/20 μL`)
df$Total <- as.double(df$Total)


## List of five SQ samples that should be removed
df[!(grepl("dpi",df$`Days Post Infection at Necropsy`)),]


## Remove the five SQ Samples, leaving total of 63 Samples
df <- df[!(df$`Gerbil ID` %in% c("M058_SQ_Bp", "M059_SQ_Bp", "F014_SQ_Bp", "M023_SQ_Bp", "F027_SQ_Bp")),]
## Remove Poorly Infected IP Sample, leaving total of 62
df <- df[!(df$`Gerbil ID` %in% c("M201_IP_Bp")),]

## Add in location
df <- df %>% mutate(location = ifelse(Litter == "CRL", "CRL", "in-house"))

## Subset Methods
SQ <- df %>% filter(`Method of Exposure` == "Subcutaneous injection")

## This detectes the difference between male gerbils and female gerbils re worm recovery (Infection Intensity)
## Power Calculation
power = 0.8
alpha = 0.05

M_SQ_mean <- mean(na.omit(SQ[grepl('M', SQ$Sex),]$Total))
F_SQ_mean <- mean(na.omit(SQ[grepl('F', SQ$Sex),]$Total))
SQ_SD <- (sd(na.omit(SQ$Total)))
pwr.t.test(d = (M_SQ_mean-F_SQ_mean)/SQ_SD, power = power, sig.level = alpha, type = "one.sample") ## SQ



## Chi Square of microfilariae prescence
SQ <- SQ %>% mutate(successful_infection = ifelse(`mf/20 μL` > 0, "yes", "no"))
SQ <- SQ %>% select("Sex", "successful_infection")
chisq.test(table(SQ$Sex, SQ$successful_infection))





### Using all passing SQ gerbils
## Male Gerbils had more female worms recovered than female gerbils, but fewer male worms as a raw sum
median(SQ[SQ$Sex == "F",]$`# Male`)
median(SQ[SQ$Sex == "M",]$`# Male`)

median(SQ[SQ$Sex == "F",]$`# Female`)
median(SQ[SQ$Sex == "M",]$`# Female`)

## What about those with successful infection, male gerbils had more male and female worms (sum), not for median, more female worms as median (very small difference)
mean(SQ[SQ$Sex == "F" & SQ$successful_infection == "yes",]$`# Male`)
mean(SQ[SQ$Sex == "M" & SQ$successful_infection == "yes",]$`# Male`)

mean(SQ[SQ$Sex == "F" & SQ$successful_infection == "yes",]$`# Female`)
mean(SQ[SQ$Sex == "M" & SQ$successful_infection == "yes",]$`# Female`)


## Plot number days infected as histogram
### Days infected
SQ$`Days Post Infection at Necropsy` <- as.double(stringr::str_replace_all(SQ$`Days Post Infection at Necropsy`, " dpi", ""))
SQ <- SQ %>% mutate(`Days Post Infection at Necropsy Binned` = ifelse(.$`Days Post Infection at Necropsy` < 100, "0-99",
                                                          ifelse(.$`Days Post Infection at Necropsy` >= 100 & .$`Days Post Infection at Necropsy` < 200, "100-199",
                                                                 ifelse(.$`Days Post Infection at Necropsy` >= 200, "200_or_more", NA))))



ggplot(SQ, aes(x = `Days Post Infection at Necropsy Binned`, y = `# Female`)) + geom_boxplot() + facet_wrap(successful_infection~Sex) + ylim(0, 8)


## Statistical test for number of days infected
table(SQ$`Days Post Infection at Necropsy Binned`, SQ$Sex)
chisq.test(table(SQ$Sex, SQ$`Days Post Infection at Necropsy Binned`))



## mf levels vs worm recovery
SQ$Age_of_infection_when_mf_tested <- difftime(as.Date(SQ$`Date mf tested`, format = "%m/%d/%Y"), as.Date(SQ$`Date exposed`,format = "%m/%d/%Y"), units = "days")
SQ$Age_of_infection_when_mf_tested <- as.double(SQ$Age_of_infection_when_mf_tested)

SQ <- SQ %>% mutate(Age_of_infection_when_mf_tested_binned = ifelse(SQ$Age_of_infection_when_mf_tested < 150, "123-150",
                                                                    ifelse(SQ$Age_of_infection_when_mf_tested >= 150, "150_or_more", "NA")))





table(SQ$Age_of_infection_when_mf_tested_binned, SQ$Sex)
chisq.test(table(SQ$Sex, SQ$Age_of_infection_when_mf_tested_binned))


## Also age when mf tested and mf levels
ggplot(SQ, aes(x = Age_of_infection_when_mf_tested_binned, y = `mf/20 μL`)) + geom_boxplot() + facet_wrap(successful_infection~Sex)
