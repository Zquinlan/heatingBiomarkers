## Plan-c analysis pipeline from Mo'orea 2019 thermal bleaching experiment
## Script created Oct. 24, 2019 - Zach Quinlan

# READING -- Libraries ----------------------------------------------------
#Data mungering
library(tidyverse)
library(data.table)
library(DescTools)
library(broom)
library(readxl)
library(multcomp)
library(CHNOSZ)
library(furrr)
library(future)
library(readxl)

#PCoA, PERMANOVA
library(vegan)
library(ape)
library(wesanderson)
library(RColorBrewer)

select <- dplyr::select

nameSamples <- function(x) { 
  separate(x, Code, c("Experiment", "Health", "Day"), sep = "_", remove = FALSE)%>%
    separate(Health, c("Health", "Replicate"), sep = 2)%>%
    separate(Day, c("Day", "Temperature"), sep = 2)%>%
    mutate(Health = case_when(Health == "HE" ~ "Healthy",
                              Health == "PB" ~ "Partially bleached",
                              Health == "BL" ~ "Bleached",
                              Health == "WA" ~ "Water control",
                              TRUE ~ as.character(Health)),
           Day = gsub("T", "", Day),
           Temperature = case_when(Temperature == "C" ~ "Ambient",
                                   Temperature == "H" ~ "Heated",
                                   TRUE ~ as.character(Temperature)))
}

zscore <- function(x) {
  (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}

genTheme <- function(x) {
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 25),
    legend.key = element_rect(fill = "transparent"),
    legend.key.size = unit(3, 'line'),
    legend.box= "vertical",
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title = element_text(size = 30),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank()) # get rid of minor grid
}

# READING -- dataframes ---------------------------------------------------
# FCM data
setwd("~/Documents/GitHub/Plan_C/data/raw/fcm/")
fcmFilesRaw <- dir(path = "~/Documents/GitHub/Plan_C/data/raw/fcm/", pattern = "*.csv")%>%
  map(read_csv)%>%
  map(~ dplyr::select(., c("Plate Name", "Sample Name", "Gate Name", "Event Count")))%>%
  reduce(bind_rows)%>%
  rename("plate" = "Plate Name",
         "Well" = "Sample Name")

fcmNamesRaw <- dir(path = "~/Documents/GitHub/Plan_C/data/raw/fcm/", pattern = "*.xlsx")

fcmNames <- fcmNamesRaw%>%
  map(read_xlsx)%>%
  set_names(fcmNamesRaw)

## Surface area
setwd("~/Documents/GitHub/Plan_C/data")
surfaceAreaRaw <- read_csv("~/Documents/GitHub/Plan_C/data/raw/physiology/Nubbin_sizing_notes.csv")

docExclude <- read_xlsx("~/Documents/GitHub/Plan_C/data/raw/doc/DOC_exlude_variance.xlsx")$`Sample Code`%>%
  as.vector()

# DOC data
pcDocRaw <- read_xlsx("~/Documents/GitHub/Plan_C/data/raw/doc/PlanC_DOC_upload.xlsx")%>%
  filter(!`Sample Code` %in% docExclude)%>%
  rename(Code = `Sample Code`)%>%
  nameSamples()
  
# Symbiodinium fcm data
symRaw <- read_csv('~/Documents/GitHub/Plan_C/data/raw/physiology/PLANC_Blastate_FCM.csv')
symNames <- read_csv('~/Documents/GitHub/Plan_C/data/raw/physiology/PLANC_Blastate_FCM_Sample_Map.csv')%>%
  mutate(Plate = gsub('0', '', Plate))

# fDOM data
fdomRaw <- read_xlsx('~/Documents/GitHub/Plan_C/data/raw/fdom/PCfDOM.xlsx')

# 16s microbial data
microbialCompositionRaw <- read_tsv('~/Documents/GitHub/Plan_C/data/raw/16s/abundance_table_100.database.tsv')%>%
  select(-c(repSeqName, repSeq))%>%
  select(OTUNumber, OTUConTaxonomy, everything())%>%
  gather(sample, reads, 3:ncol(.))

microbialSampleNames <- read_csv('~/Documents/GitHub/Plan_C/data/raw/16s/sequenceLabels.csv')%>%
  select(`Sampe ID`, `Sample name`)%>%
  rename(sample = 1,
         Code = 2)%>%
  filter(Code %like% 'PC%')
  

# CLEANING -- FCM ---------------------------------------------------------
fcmRaw <- fcmNames%>%
  map_df( ~as.data.frame(.), .id="plate")%>%
  mutate(plate = gsub(".Map.xlsx", "", plate),
         plate = gsub("OCT.2019", "OCT.18.2019", plate))%>%
  full_join(., fcmFilesRaw, by = c("plate", "Well"))%>%
  separate(plate, c("experiment_plate", "dye", "month", "day", "year", "plate", "number"))%>%
  unite(date, c("month", "day", "year"), sep = "-")%>%
  unite(plate_numnber, c("plate", "number"), sep = "_")%>%
  # filter(!Sample %in% docExclude)%>%
  # filter(Sample != 'PC_HE2_T1H',
  #        Sample != 'PC_PB2_T1H',
  #        Sample != 'PC_BL3_T2H',
  #        Sample != 'PC_HE3_T2H',
  #        Sample != 'PC_BL3_T3H',
  #        Sample != 'PC_PB2_T3H',
  #        Sample != 'PC_HE2_T3H',
  #        Sample != 'PC_WA2_T4H',
  #        Sample != 'PC_HE1_T5H',
  #        Sample != 'PC_BL1_T6H',
  #        Sample != 'PC_PB3_T6H',
  #        Sample != 'PC_BL1_T7H')%>%
  separate(Sample, c("Experiment", "Health", "Day"), sep = "_", remove = FALSE)%>%
  separate(Health, c("Health", "Replicate"), sep = 2)%>%
  separate(Day, c("Day", "Temperature"), sep = 2)%>%
  mutate(Health = case_when(Health == "HE" ~ "Healthy",
                            Health == "PB" ~ "Partially bleached",
                            Health == "BL" ~ "Bleached",
                            Health == "WA" ~ "Water control",
                            TRUE ~ as.character(Health)),
         Day = gsub("T", "", Day),
         Temperature = case_when(Temperature == "C" ~ "Ambient",
                                 Temperature == "H" ~ "Heated",
                          TRUE ~ as.character(Temperature)))

fcmSyber <- fcmRaw%>%
  filter(Experiment == "PC",
         dye == "SYBR",
         `Gate Name` == "SYBR Polygon")%>%
  # filter(bath == "hot")%>%
  group_by(Day, Health, Temperature)%>%
  mutate(cells_mL = `Event Count`/0.075,
         cells_ul = cells_mL/1000,
         log10_ul = log10(cells_ul),
         std = sd(log10_ul))

fcmMeans <- fcmSyber%>%
  select(-c(1:5, Replicate))%>%
  group_by(Day, Health, Temperature)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()




# CLEANING -- Symbiodinium counts -----------------------------------------
symBad <- c('PC_BL2_T0C_A',
             'PC_BL3_T0C_D',
             'PC_HE3_T0C_A',
             'PC_HE3_T0C_E',
             'PC_PB3_T7C_D',
             'PC_PB3_T7H_E',
             'PC_BL1_T7C_A',
             'PC_BL1_T7H_E',
             'PC_BL2_T7C_A',
             'PC_BL3_T7C_D',
             'PC_BL3_T7C_E',
             'PC_BL3_T7H_A',
             'PC_HE1_T7C_B',
             'PC_HE1_T7C_E',
             'PC_HE1_T7H_D',
             'PC_HE2_T7C_F',
             'PC_BL1_T7C_E_v1')

symStats <- symNames%>%
  left_join(symRaw%>%
              mutate(Plate = as.character(Plate)), by = c('Plate', 'Well'))%>%
  mutate(count = case_when(is.na(Notes) ~ `Sym-FSC Events`,
                              !is.na(Notes) ~ `Sym-SSC Events`))%>%
  select(Sample, count)%>%
  rename(Code = Sample)%>%
  filter(!Code %in% symBad)%>%
  left_join(surfaceAreaRaw%>%
              select(Code, SA, Notes),
            by = 'Code')%>%
  mutate(sym_sa = count/SA)%>%
  nameSamples()%>%
  filter(!Code %like% 'mil%')%>%
  separate(Code, c('Code', 'nubbin'), sep = -1)%>%
  mutate(log = log10(sym_sa))

  



# CLEANING -- fDOM --------------------------------------------------------
fdomClean <- fdomRaw%>%
  rename(Code = 1)%>%
  filter(!Code %like% 'BLANK%')%>%
  select(-c(2:6, `Optical Brighteners`:`Lignin-like`))%>%
  nameSamples()%>%
  gather(fluor, rawVal, 7:ncol(.))%>%
  mutate(log10 = log10(rawVal + 1))


# CLEANING -- 16s microbial composition ---------------------------------------
microbialComposition <- microbialCompositionRaw%>%
  inner_join(microbialSampleNames, by = 'sample')%>%
  nameSamples()%>%
  group_by(Code)%>%
  mutate(numOtus = nrow(.),
         reads = case_when(reads == 0 ~ 1/numOtus,
                           TRUE ~ as.numeric(reads)),
         ra = reads/sum(reads),
         asin = asin(sqrt(ra)),
         log = log10(ra),
         logit = Logit(ra))%>%
  ungroup()%>%
  group_by(OTUNumber)%>%
  mutate(abundant = case_when(max(ra) > 0.05 | sum(ra > 0.001) >= 3 ~ "abundant",
                              TRUE ~ "rare"))%>%
  filter(abundant == 'abundant')%>%
  ungroup()



greaterThan70 <- microbialComposition%>%
  group_by(OTUNumber, OTUConTaxonomy, Health, Temperature, Day)%>%
  summarize_if(is.numeric, mean)%>%
  arrange(-ra)%>%
  ungroup()%>%
  group_by(Health, Temperature, Day)%>%
  nest()%>%
  mutate(data = map(data, ~ mutate(.x, cum = cumsum(ra))))%>%
  unnest(data)%>%
  ungroup()%>%
  # filter(ra >= 0.05)
  filter(cum <= 0.7)%>%
  select(OTUNumber)%>%
  unique()

# PRE-STATS -- qqplots and filtering ----------------------------------------------------
pdf('plots/qqplots.pdf', width = 15, height = 10)
car::qqPlot(microbialComposition$ra, 
            ylab = "Relative Abundance", xlab = "Normal quantiles",
            main = 'QQ-plot: Relative Abundance') 

car::qqPlot(microbialComposition$asin, 
            ylab = "Angularly transformed relative abundance", xlab = "Normal quantiles",
            main = 'QQ-plot: Angularly transformation') 

car::qqPlot(microbialComposition$log, 
            ylab = "Log10 relative abundance", xlab = "Normal quantiles",
            main = 'QQ-plot: Log10 transformation') 

car::qqPlot(microbialComposition$logit, 
            ylab = "Logit relative abundance", xlab = "Normal quantiles",
            main = 'QQ-plot: Logit transformation') 
dev.off()


# STATS -- 16s Two-Way ANOVA ----------------------------------------------
threeWayAnova <- microbialComposition%>%
  select(OTUNumber, Health, Replicate, Day, Temperature, asin)%>%
  group_by(OTUNumber)%>%
  nest()%>%
  mutate(anova = map(data, ~ aov(asin ~ Health*Temperature*Day, .x)%>%
                       tidy()%>%
                       filter(!term == "Residuals")%>%
                       select(term, p.value)))%>%
  select(-data)%>%
  unnest(anova)%>%
  mutate(fdr = p.adjust(p.value, method = 'BH'))%>%
  filter(fdr <= 0.05)

differentMicrobes <- threeWayAnova%>%
  filter(term == 'Temperature' | term == 'Temperature:Day' | term == 'Health:Temperature:Day')%>%
  select(OTUNumber)%>%
  unique()

allDifferentMicrobes <- threeWayAnova%>%
  select(OTUNumber)%>%
  unique()

# GRAPHING -- Color Pallette and theme----------------------------------------------
healthColors <- c("#FD6467", "#E1AF00", "#0B775E", "#78B7C5")

# GRAPHING -- FCM ---------------------------------------------------------
pdf("~/Documents/GitHub/Plan_C/data/plots/Plan_c_FCM.pdf", width = 15, height = 10)
# Ambient for presentations
fcmSyber%>%
  unite(cb, c(Health, Temperature), sep = "_", remove = FALSE)%>%
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  rename(`Initial Health` = Health,)%>%
  filter(Temperature == 'Ambient')%>%
  ggplot(aes(x = Day, y = cells_ul, color = `Initial Health`)) +
  geom_point(aes(group = cb), stat = "identity", size = 2) + 
  # geom_errorbar(aes(ymin = mean(cells_mL) - std, ymax = mean(cells_mL) + std)) +
  geom_line(aes(group = cb, linetype = Temperature), stat = "summary", fun.y = "mean", lwd = 1.4) +
  ylab("Cells µL-1") +
  scale_linetype_manual(values=c("dashed", "solid")) +
  scale_color_manual(values = healthColors) +
  scale_y_log10(limits = c(100,2600)) +
  ylab(bquote(Cells ~µL^-1)) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 25),
    legend.key = element_rect(fill = "white",  color = "white"),
    legend.key.size = unit(3, 'line'),
    legend.box= "vertical",
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title = element_text(size = 30),
    panel.background = element_rect(fill = "white",  color = "white"), # bg of the panel
    plot.background = element_rect(fill = "white",  color = "white"), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    # legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    # legend.box.background = element_rect(fill = "transparent"))
  ) 

fcmSyber%>%
  unite(cb, c(Health, Temperature), sep = "_", remove = FALSE)%>%
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  rename(`Initial Health` = Health)%>%
  filter(Temperature != 'Ambient')%>%
  ggplot(aes(x = Day, y = cells_ul, color = `Initial Health`)) +
  geom_point(aes(group = cb), stat = "identity", size = 2) + 
  # geom_errorbar(aes(ymin = mean(cells_mL) - std, ymax = mean(cells_mL) + std)) +
  geom_line(aes(group = cb, linetype = Temperature), stat = "summary", fun.y = "mean", lwd = 1.4) +
  ylab("Cells µL-1") +
  scale_linetype_manual(values=c("solid")) +
  scale_color_manual(values = healthColors) +
  scale_y_log10(limits = c(100,2600)) +
  ylab(bquote(Cells ~µL^-1)) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 25),
    legend.key = element_rect(fill = "white",  color = "white"),
    legend.key.size = unit(3, 'line'),
    legend.box= "vertical",
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title = element_text(size = 30),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent"), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    # legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    # legend.box.background = element_rect(fill = "transparent"))
  ) 

fcmSyber%>%
  unite(cb, c(Health, Temperature), sep = "_", remove = FALSE)%>%
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  rename(`Initial Health` = Health,)%>%
  # filter(Temperature == 'Ambient')%>%
  ggplot(aes(x = Day, y = cells_ul, color = `Initial Health`)) +
  geom_point(aes(group = cb), stat = "identity", size = 2) + 
  # geom_errorbar(aes(ymin = mean(cells_mL) - std, ymax = mean(cells_mL) + std)) +
  geom_line(aes(group = cb, linetype = Temperature), stat = "summary", fun.y = "mean", lwd = 1.4) +
  ylab("Cells µL-1") +
  scale_linetype_manual(values=c("dashed", "solid")) +
  scale_color_manual(values = healthColors) +
  scale_y_log10(limits = c(100,2600)) +
  ylab(bquote(Cells ~µL^-1)) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 25),
    legend.key = element_rect(fill = "white",  color = "white"),
    legend.key.size = unit(3, 'line'),
    legend.box= "vertical",
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title = element_text(size = 30),
    panel.background = element_rect(fill = "white",  color = "white"), # bg of the panel
    plot.background = element_rect(fill = "white",  color = "white"), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    # legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    # legend.box.background = element_rect(fill = "transparent"))
  ) 
dev.off()

pdf("~/Documents/GitHub/Plan_C/data/plots/Plan_c_FCM.pdf", width = 15, height = 10)
fcmMeans%>%
  # filter(Health != 'Bleached')%>%
  rename(`Starting Health` = Health)%>%
  ggplot(aes(x = Day, y = cells_ul, color = `Starting Health`))+
  geom_point(stat = "summary", fun.y = "mean") + 
  # geom_errorbar(aes(ymin = log10_ul - std, ymax =  log10_ul + std)) +
  geom_line(aes(group = `Starting Health`), stat = "summary", fun.y = "mean") +
  ylab(bquote(Cells ~µL^-1)) +
  theme(
    legend.position = "top",
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"))
# facet_wrap(~bath)
dev.off()



# GRAPHING -- Symbiodinium counts -----------------------------------------

pdf('~/Documents/GitHub/Plan_C/data/plots/Plan_c_SymbiodiniumCounts.pdf', width = 15, height = 10)
symStats%>%
  filter(!Code %like% '%T7C%')%>%
  unite(dayHealth, c('Day', 'Health'), sep = "_", remove = FALSE)%>%
  mutate(Day = as.factor(Day))%>% 
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  ggplot(aes(x = Day, y = log, color = Health)) + 
  geom_boxplot(aes(group = dayHealth)) +
  scale_color_manual(values = healthColors) +
  facet_wrap(~ Health) +
  labs(y = bquote(Log[10] ~(Symbiodinium ~cm^-2)), title = 'Symbiodinium counts heated treatment') +
  genTheme()

symStats%>%
  filter(!Code %like% '%T7H%')%>%
  unite(dayHealth, c('Day', 'Health'), sep = "_", remove = FALSE)%>%
  mutate(Day = as.factor(Day))%>% 
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  ggplot(aes(x = Day, y = log, color = Health)) + 
  geom_boxplot(aes(group = dayHealth)) +
  scale_color_manual(values = healthColors) +
  facet_wrap(~ Health) +
  labs(y = bquote(Log[10] ~(Symbiodinium ~cm^-2)), title = 'Symbiodinium counts Ambient treatment') +
  genTheme()
dev.off()


# GRAPHING -- DOC ---------------------------------------------------------
pdf("~/Documents/GitHub/Plan_C/data/plots/Plan_c_DOC.pdf", width = 15, height = 10)
pcDocRaw%>%
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  unite(cb, c(Health, Temperature), sep = "_", remove = FALSE)%>%
  rename(`Initial Health` = Health)%>%
  filter(Temperature == 'Ambient')%>%
  ggplot(aes(Day, `DOC (µMC)`, color = `Initial Health`)) +
  geom_point(stat = 'identity') +
  geom_line(aes(group = cb, linetype = Temperature), stat = "summary", fun.y = "mean", lwd = 1.4) +
  scale_linetype_manual(values=c("dashed", "solid")) +
  scale_color_manual(values = healthColors) +
  scale_x_discrete(limits = c(1,2,3,4,5,6,7)) +
  labs(x = 'Day') +
  ylim(65,176) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 25),
    legend.key = element_rect(fill = "transparent"),
    legend.key.size = unit(3, 'line'),
    legend.box= "vertical",
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title = element_text(size = 30),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank()) # get rid of minor grid 


pcDocRaw%>%
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  unite(cb, c(Health, Temperature), sep = "_", remove = FALSE)%>%
  rename(`Initial Health` = Health)%>%
  filter(Temperature != 'Ambient')%>%
  ggplot(aes(Day, `DOC (µMC)`, color = `Initial Health`)) +
  geom_point(stat = 'identity') +
  geom_line(aes(group = cb, linetype = Temperature), stat = "summary", fun.y = "mean", lwd = 1.4) +
  scale_linetype_manual(values=c("solid")) +
  scale_color_manual(values = healthColors) +
  scale_x_discrete(limits = c(1,2,3,4,5,6,7)) +
  labs(x = 'Day', y = "Dissolved Orgnic Carbon (µMC)") +
  ylim(65,176) +
  genTheme()

pcDocRaw%>%
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  unite(cb, c(Health, Temperature), sep = "_", remove = FALSE)%>%
  rename(`Initial Health` = Health)%>%
  # filter(Temperature != 'Ambient')%>%
  ggplot(aes(Day, `DOC (µMC)`, color = `Initial Health`)) +
  geom_point(stat = 'identity') +
  geom_line(aes(group = cb, linetype = Temperature), stat = "summary", fun.y = "mean", lwd = 1.4) +
  scale_linetype_manual(values=c("dashed", "solid")) +
  scale_color_manual(values = healthColors) +
  scale_x_discrete(limits = c(1,2,3,4,5,6,7)) +
  labs(x = 'Day', y = "Dissolved Orgnic Carbon (µMC)") +
  ylim(65,176) +
  genTheme()
dev.off()

# GRAPHING -- combine DOC FCM bad -----------------------------------------------------
fcmDoc <- fcmSyber%>%
  full_join(pcDocRaw%>%
              rename(Sample = Code), by = c('Sample', 'Health', 'Day', 'Temperature'))

fcmDoc%>%
  ungroup()%>%
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  unite(cb, c(Health, Temperature), sep = "_", remove = FALSE)%>%
  rename(`Initial Health` = Health)%>%
  filter(Temperature != 'Ambient')%>%
  ggplot(aes(x = Day, color = `Initial Health`)) +
  # geom_point(aes(y = `DOC (µMC)`), stat = 'identity') +
  geom_line(aes(y = `DOC (µMC)`, group = cb, linetype = "dashed"), stat = "summary", fun.y = "mean", lwd = 1.4) +
  # geom_point(aes(y = `cells_ul`)) +
  geom_line(aes(y = `cells_ul`, group = cb, linetype = "solid"), stat = "summary", fun.y = "mean", lwd = 1.4) +
  scale_y_continuous(name = "DOC (µMC)",
  sec.axis = sec_axis(~.*100, name= bquote(Cells ~µL^-1))) +
  scale_color_manual(values = healthColors) +
  genTheme()



# GRAPHING -- fDOM --------------------------------------------------------
pdf('plots/plan_c_fDOM.pdf', width = 15, height = 10)
fdomClean%>%
  filter(!fluor %like any% c('M:C', 'BIX', 'HIX', 'FI', 'PARAFAC%'))%>%
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  unite(cb, c(Health, Temperature), sep = "_", remove = FALSE)%>%
  rename(`Initial Health` = Health)%>%
  filter(!Temperature == 'Ambient')%>%
  ggplot(aes(Day, log10, color = `Initial Health`)) +
  geom_point(stat = 'identity') +
  geom_line(aes(group = cb, linetype = Temperature), stat = "summary", fun.y = "mean", lwd = 1.4) +
  scale_linetype_manual(values=c("solid")) +
  scale_color_manual(values = healthColors) +
  scale_x_discrete(limits = c(1,2,3,4,5,6,7)) +
  labs(x = 'Day') +
  facet_wrap(~fluor, scales = 'free_y') +
  genTheme()

fdomClean%>%
  filter(!fluor %like any% c('M:C', 'BIX', 'HIX', 'FI', 'PARAFAC%'))%>%
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  unite(cb, c(Health, Temperature), sep = "_", remove = FALSE)%>%
  rename(`Initial Health` = Health)%>%
  filter(Temperature == 'Ambient')%>%
  ggplot(aes(Day, log10, color = `Initial Health`)) +
  geom_point(stat = 'identity') +
  geom_line(aes(group = cb, linetype = Temperature), stat = "summary", fun.y = "mean", lwd = 1.4) +
  scale_linetype_manual(values=c("dashed")) +
  scale_color_manual(values = healthColors) +
  scale_x_discrete(limits = c(1,2,3,4,5,6,7)) +
  labs(x = 'Day') +
  facet_wrap(~fluor, scales = 'free_y') +
  genTheme()

fdomClean%>%
  filter(!fluor %like any% c('M:C', 'BIX', 'HIX', 'FI', 'PARAFAC%'))%>%
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  unite(cb, c(Health, Temperature), sep = "_", remove = FALSE)%>%
  rename(`Initial Health` = Health)%>%
  # filter(Temperature == 'Ambient')%>%
  ggplot(aes(Day, log10, color = `Initial Health`)) +
  geom_point(stat = 'identity') +
  geom_line(aes(group = cb, linetype = Temperature), stat = "summary", fun.y = "mean", lwd = 1.4) +
  scale_linetype_manual(values=c("dashed", "solid")) +
  scale_color_manual(values = healthColors) +
  scale_x_discrete(limits = c(1,2,3,4,5,6,7)) +
  labs(x = 'Day') +
  facet_wrap(~fluor, scales = 'free_y') +
  genTheme()

dev.off()


# GRAPHING -- 16s change --------------------------------------------------

pdf('plots/abundantDiversity.pdf', width = 15, height = 10)
microbialComposition%>%
  inner_join(greaterThan70, by = 'OTUNumber')%>%
  select(OTUNumber, OTUConTaxonomy, Code, Health, Replicate, Day, Temperature, ra)%>%
  separate(OTUConTaxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";")%>%
  group_by(Class, Health, Day, Temperature, Replicate)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(Class, Health, Day, Temperature)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  # filter(ra > 0.07)%>%
  unite(treatment, c(Temperature, Health), sep = '  ')%>%
  ggplot(aes(Day, ra, fill = Class)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~treatment, nrow = 2) +
  gen_theme()
dev.off()  

# pdf('plots/abundantDiversity.pdf', width = 15, height = 10)
microbialComposition%>%
  inner_join(greaterThan70, by = 'OTUNumber')%>%
  select(OTUNumber, OTUConTaxonomy, Code, Health, Replicate, Day, Temperature, ra)%>%
  separate(OTUConTaxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";")%>%
  filter(Family == 'Rhodobacteraceae',
         Temperature == 'Ambient',
         Day == '1' | Day == '3',
         Health == 'Bleached' | Health == 'Healthy')%>%
  group_by(Genus, Health, Day, Temperature, Replicate)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(Genus, Health, Day, Temperature)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  unite(treatment, c(Temperature, Health), sep = '  ')%>%
  ggplot(aes(Day, ra, fill = Genus)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~treatment, nrow = 2) +
  gen_theme()
# dev.off()  

# GRAPHING -- 16S bicluster ------------------------------------
biclusterDf <- microbialComposition%>%
  inner_join(greaterThan70, by = 'OTUNumber')%>%
  inner_join(allDifferentMicrobes, by = 'OTUNumber')%>%
  filter(Day != '7')%>%
  select(OTUNumber, OTUConTaxonomy, Health, Replicate, Day, Temperature, asin)%>%
  unite(sample, c(Temperature, Health, Day, Replicate), sep = '_')%>%
  group_by(OTUConTaxonomy)%>%
  mutate(asin = zscore(asin))%>%
  ungroup()%>%
  separate(OTUConTaxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";")%>%
  unite(OFGO, c("Order", "Family", "Genus", "OTUNumber"), sep = ";")%>%
  select(-c("Kingdom", "Phylum", "Class", "OTU_id"))%>%
  spread(OFGO, asin)

write_csv(biclusterDf, 'analysis/biclusterMicrobial.csv')


# GRAPHING -- 16s PCoA ----------------------------------------------------
pcoaDf <- microbialComposition%>%
  filter(Day != 7,
         Temperature == 'Heated')%>%
  inner_join(allDifferentMicrobes, by = 'OTUNumber')%>%
  select(OTUNumber, OTUConTaxonomy, Health, Replicate, Day, Temperature, asin)%>%
  unite(sample, c(Temperature, Health, Day, Replicate), sep = '_')%>%
  separate(OTUConTaxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";")%>%
  unite(OFGO, c("Order", "Family", "Genus", "OTUNumber"), sep = ";")%>%
  select(-c("Kingdom", "Phylum", "Class", "OTU_id"))%>%
  spread(OFGO, asin)%>%
  column_to_rownames(var = "sample")%>%
  vegdist(na.rm = TRUE)%>%
  pcoa()


pdf('plots/pcoaMicrobes.pdf', width = 15, height = 10)
pcoaDf$vectors%>%
  as.data.frame()%>%
  rownames_to_column(var = "sample")%>%
  separate(sample, c('Temperature', 'Health', 'Day', 'Replicate'), sep = '_')%>%
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  ggplot(., aes(x = Axis.1, y = Axis.2, color = Health, shape = Day)) +
  geom_text(aes(label = Day), cex = 13) +
  scale_color_manual(values = healthColors) +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(face = "italic"),
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 25)) +
  xlab(str_c("Axis 1", " (", round(pcoaDf$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
  ylab(str_c("Axis 2", " (", round(pcoaDf$values$Relative_eig[2], digits = 4)*100, "%)", sep = ""))
dev.off()




