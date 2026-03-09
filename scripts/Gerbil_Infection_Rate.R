#-------------------------------------------------------
## This script will download the worm recovery table
## and generate table/plot of the number of gerbils
## where an infection took hold (Fig 1B)
#-------------------------------------------------------


library(tidyverse)

WORM_RECOVERY_URL <- "https://raw.githubusercontent.com/christopher-holt/bpahangi_host_sex_effect/refs/heads/main/input_files/Worm_Recovery_Table.csv"

df <- read_csv(WORM_RECOVERY_URL)

cwd <- getwd()
OUTPUT_DIR <- paste0(cwd, "/output_files/")
FIGURE_DIR <- paste0(OUTPUT_DIR, '/figures/')

dir.create(OUTPUT_DIR, showWarnings = F, recursive = T)
dir.create(FIGURE_DIR, showWarnings = F, recursive = T)


## List of five SQ samples that should be removed
df[!(grepl("dpi",df$`Days Post Infection at Necropsy`)),]

## Remove the five SQ Samples, leaving total of 63 samples
df <- df[!(df$`Gerbil ID` %in% c("M058_SQ_Bp", "M059_SQ_Bp", "F014_SQ_Bp", "M023_SQ_Bp", "F027_SQ_Bp")),]

## Infection Table
df <- df %>% mutate(Exposure_Status = ifelse(`Method of Exposure` == "Intraperitoneal injection" & !is.na(Total), 'Successful_Infection', 
                                             ifelse(`Method of Exposure` == "Subcutaneous injection" & `mf/20 μL` > 0, 'Successful_Infection', 'Injected')))



gerbil_infection_summary <- tibble('Sex' = as.character() ,
                                   'Exposure_Method' = as.character(),
                                   'Exposure_Status' = as.character(),
                                   'Number_of_Gerbils' = as.numeric())

gerbil_infection_summary <- gerbil_infection_summary %>% add_row(Sex = 'F', Exposure_Method = 'Subcutaneous', 
                                                                 Exposure_Status = "Injected",
                                                                 Number_of_Gerbils = nrow(df %>% filter(Sex == "F" & 
                                                                                                          `Method of Exposure` == "Subcutaneous injection")))

gerbil_infection_summary <- gerbil_infection_summary %>% add_row(Sex = 'F', Exposure_Method = 'Subcutaneous', 
                                                                 Exposure_Status = "Successful_Infection",
                                                                 Number_of_Gerbils = nrow(df %>% filter(Sex == "F" & 
                                                                                                          `Method of Exposure` == "Subcutaneous injection" & 
                                                                                                          Exposure_Status == "Successful_Infection")))

gerbil_infection_summary <- gerbil_infection_summary %>% add_row(Sex = 'M', Exposure_Method = 'Subcutaneous', 
                                                                 Exposure_Status = "Injected",
                                                                 Number_of_Gerbils = nrow(df%>% filter(Sex == "M" & 
                                                                                                         `Method of Exposure` == "Subcutaneous injection")))

gerbil_infection_summary <- gerbil_infection_summary %>% add_row(Sex = 'M', Exposure_Method = 'Subcutaneous', 
                                                                 Exposure_Status = "Successful_Infection",
                                                                 Number_of_Gerbils = nrow(df %>% filter(Sex == "M" & 
                                                                                                          `Method of Exposure` == "Subcutaneous injection" & 
                                                                                                          Exposure_Status == "Successful_Infection")))

gerbil_infection_summary <- gerbil_infection_summary %>% add_row(Sex = 'F', Exposure_Method = 'Intraperitoneal', 
                                                                 Exposure_Status = "Injected",
                                                                 Number_of_Gerbils = nrow(df %>% filter(Sex == "F" & 
                                                                                                          `Method of Exposure` == "Intraperitoneal injection")))

gerbil_infection_summary <- gerbil_infection_summary %>% add_row(Sex = 'F', Exposure_Method = 'Intraperitoneal', 
                                                                 Exposure_Status = "Successful_Infection",
                                                                 Number_of_Gerbils = nrow(df %>% filter(Sex == "F" & 
                                                                                                          `Method of Exposure` == "Intraperitoneal injection" & 
                                                                                                          Exposure_Status == "Successful_Infection")))

gerbil_infection_summary <- gerbil_infection_summary %>% add_row(Sex = 'M', Exposure_Method = 'Intraperitoneal', 
                                                                 Exposure_Status = "Injected",
                                                                 Number_of_Gerbils = nrow(df %>% filter(Sex == "M" & 
                                                                                                          `Method of Exposure` == "Intraperitoneal injection")))

## 7th gerbil was successfully infected, but was removed due to low infection numbers
gerbil_infection_summary <- gerbil_infection_summary %>% add_row(Sex = 'M', Exposure_Method = 'Intraperitoneal', 
                                                                 Exposure_Status = "Successful_Infection",
                                                                 Number_of_Gerbils = nrow(df %>% filter(Sex == "M" & 
                                                                                                          `Method of Exposure` == "Intraperitoneal injection")))


gerbil_infection_summary <- gerbil_infection_summary %>% mutate(total_number_of_gerbils = ifelse(Exposure_Status == "Injected", .$Number_of_Gerbils, NA))
gerbil_infection_summary <- gerbil_infection_summary %>% fill(total_number_of_gerbils)
gerbil_infection_summary <- gerbil_infection_summary %>% mutate(percentage = ifelse(Exposure_Status == "Successful_Infection", 
                                                                                    Number_of_Gerbils/total_number_of_gerbils, NA)) %>% 
  fill(percentage, .direction = "up") %>% mutate(percentage = ifelse(Exposure_Status == "Injected", 1-percentage, .$percentage))


pdf(paste0(FIGURE_DIR,"Fig1_B.pdf"),
    width = 10,
    height = 10)
    print(
      ggplot(gerbil_infection_summary, aes(x = Sex, y = percentage, fill = Exposure_Status)) + 
        geom_bar(stat = 'identity', position = "stack", colour = "black") + 
        facet_grid(~Exposure_Method) + theme_bw() + labs(x = "Gebril Sex") + guides(fill = guide_legend(title = 'Exposure Status'))  +
        scale_fill_manual(breaks = c('Injected', 'Successful_Infection'),
                          values = c('white', 'black'))
    )
dev.off()

