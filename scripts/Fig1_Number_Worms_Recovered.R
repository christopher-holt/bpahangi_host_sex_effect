#-------------------------------------------------------
## This script reads in Worm Recovery Table and 
## will generate Figures 1C through 1F and the 
## numbers used in Fig 1A
#-------------------------------------------------------


library(tidyverse)
library(tidyr)
library(ggpubr)

cwd <- getwd()
OUTPUT_DIR <- paste0(cwd, "/output_files/")
FIGURE_DIR <- paste0(OUTPUT_DIR, '/figures/')

dir.create(OUTPUT_DIR, showWarnings = F, recursive = T)
dir.create(FIGURE_DIR, showWarnings = F, recursive = T)


## Assess number of worms injected, Brugia pahangi
## Tidy Data
WORM_RECOVERY_URL <- "https://raw.githubusercontent.com/christopher-holt/bpahangi_host_sex_effect/refs/heads/main/input_files/Worm_Recovery_Table.csv"

df <- read_csv(WORM_RECOVERY_URL)

gerbil_id <- df$`Gerbil ID`



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
df[df$`Gerbil ID` == "M201_IP_Bp",]$Total <- 3
df[df$`Gerbil ID` == "M201_IP_Bp",]$`Total # Males + Females` <- 3



df$Sex <- factor(df$Sex)
df$`Method of Exposure` <- factor(df$`Method of Exposure`)
df$Age_when_infected <- difftime(as.Date(df$`Date exposed`, format = "%m/%d/%Y"), as.Date(df$DOB,format = "%m/%d/%Y"), units = "days")
df$Age_of_infection_when_mf_tested <- difftime(as.Date(df$`Date mf tested`, format = "%m/%d/%Y"), as.Date(df$`Date exposed`,format = "%m/%d/%Y"), units = "days")

## Age Range when infected
min(difftime(as.Date(df$`Date exposed`, format = "%m/%d/%Y"), as.Date(df$DOB,format = "%m/%d/%Y"), units = "days"))
max(difftime(as.Date(df$`Date exposed`, format = "%m/%d/%Y"), as.Date(df$DOB,format = "%m/%d/%Y"), units = "days"))



df.plot <- df %>% select('Sex', "Method of Exposure", "Total # Males + Females", "Days Post Infection at Necropsy", "mf/20 μL")
df.plot$`mf/20 μL` <- as.double(df.plot$`mf/20 μL`)
df.plot$`Method of Exposure` <- str_replace_all(df.plot$`Method of Exposure`, "Intraperitoneal injection", "IP")
df.plot$`Method of Exposure` <- str_replace_all(df.plot$`Method of Exposure`, "Subcutaneous injection", "SQ")
df.plot$Exposure <- df.plot$`Method of Exposure`
#df.plot <- df.plot %>% filter(`Age of Infection at Harvest` != "Omit from study")


## Panel 1F
pdf(paste0(FIGURE_DIR, "/Fig1_F.pdf"),
    width = 10,
    height = 10)
    print(
      ggplot(df.plot %>% filter(!(is.na(`Total # Males + Females`)) & `Method of Exposure` == "SQ"), aes(x = Sex, y = `mf/20 μL`)) + 
        geom_boxplot() + geom_jitter(height=0) +
        theme_bw() + stat_compare_means(method = 't.test', label.x = 1.5) + labs(y = "mf/20 microL")
    )
dev.off()

mean(df.plot[df.plot$Sex == "M" & df.plot$`Method of Exposure` == "SQ",]$`mf/20 μL`)
median(df.plot[df.plot$Sex == "M" & df.plot$`Method of Exposure` == "SQ",]$`mf/20 μL`)

mean(df.plot[df.plot$Sex == "F" & df.plot$`Method of Exposure` == "SQ",]$`mf/20 μL`)
median(df.plot[df.plot$Sex == "F" & df.plot$`Method of Exposure` == "SQ",]$`mf/20 μL`)



## Panel 1C
pdf(paste0(FIGURE_DIR, "/Fig1_C.pdf"),
    width = 10,
    height = 10)
plot(
  ggplot(df.plot, aes(x = Exposure, y = `Total # Males + Females`)) + 
    geom_boxplot() + geom_jitter(height=0) + theme(axis.text.x = element_text(angle = 45)) + stat_compare_means(method = 't.test', label.x = 1.5) + 
    ylim(-0.5, 250) + theme_bw()
)
dev.off()



## Panel 1D
pdf(paste0(FIGURE_DIR, "/Fig1_D.pdf"),
    width = 10,
    height = 10)
plot(
  ggplot(df.plot[grepl("IP", df.plot$Exposure),], aes(x = Sex, y = `Total # Males + Females`)) + 
    geom_boxplot() + geom_jitter(height=0) + theme(axis.text.x = element_text(angle = 45)) + stat_compare_means(method = 't.test', label.x = 1.5) + 
    ylim(0, 250) + theme_bw()
)
dev.off()


## Panel 1E
pdf(paste0(FIGURE_DIR, "/Fig1_E.pdf"),
    width = 10,
    height = 10)
plot(
  ggplot(df.plot[grepl("SQ", df.plot$Exposure),], aes(x = Sex, y = `Total # Males + Females`)) + 
    geom_boxplot() + geom_jitter(height = 0) + theme(axis.text.x = element_text(angle = 45)) + 
    stat_compare_means(method = 't.test', label.x = 1.3) + theme_bw() + ylim(-0.5,10) + scale_y_continuous(breaks = c(0,2,4,6,8,10))
)
dev.off()





## Avg Number of Worms Recovered by Exposure Method
mean(df.plot[df.plot$`Method of Exposure` == "IP",]$`Total # Males + Females`)
mean(df.plot[df.plot$`Method of Exposure` == "SQ",]$`Total # Males + Females`)

## Avg Number of Worms Recovered by Exposure Method and Gerbil Sex
## IP
mean(df.plot[df.plot$Sex == "M" & df.plot$`Method of Exposure` == "IP",]$`Total # Males + Females`)
mean(df.plot[df.plot$Sex == "F" & df.plot$`Method of Exposure` == "IP",]$`Total # Males + Females`)
## SQ
mean(df.plot[df.plot$Sex == "M" & df.plot$`Method of Exposure` == "SQ",]$`Total # Males + Females`)
mean(df.plot[df.plot$Sex == "F" & df.plot$`Method of Exposure` == "SQ",]$`Total # Males + Females`)



## Median Number of Worms Recovered by Exposure Method
median(df.plot[df.plot$`Method of Exposure` == "IP",]$`Total # Males + Females`)
median(df.plot[df.plot$`Method of Exposure` == "SQ",]$`Total # Males + Females`)

## Avg Number of Worms Recovered by Exposure Method and Gerbil Sex
## IP
median(df.plot[df.plot$Sex == "M" & df.plot$`Method of Exposure` == "IP",]$`Total # Males + Females`)
median(df.plot[df.plot$Sex == "F" & df.plot$`Method of Exposure` == "IP",]$`Total # Males + Females`)
## SQ
median(df.plot[df.plot$Sex == "M" & df.plot$`Method of Exposure` == "SQ",]$`Total # Males + Females`)
median(df.plot[df.plot$Sex == "F" & df.plot$`Method of Exposure` == "SQ",]$`Total # Males + Females`)
















