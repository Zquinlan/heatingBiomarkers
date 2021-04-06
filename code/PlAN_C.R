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


# READING -- dataframes ---------------------------------------------------
# FCM data
setwd("~/Documents/GitHub/Plan_C/data/raw/fcm/")
fcm_files_raw <- dir(path = "~/Documents/GitHub/Plan_C/data/raw/fcm/", pattern = "*.csv")%>%
  map(read_csv)%>%
  map(~ dplyr::select(., c("Plate Name", "Sample Name", "Gate Name", "Event Count")))%>%
  reduce(bind_rows)%>%
  rename("plate" = "Plate Name",
         "Well" = "Sample Name")

fcm_names_raw <- dir(path = "~/Documents/GitHub/Plan_C/data/raw/fcm/", pattern = "*.xlsx")

fcm_names <- fcm_names_raw%>%
  map(read_xlsx)%>%
  set_names(fcm_names_raw)

## Surface area
setwd("~/Documents/GitHub/Plan_C/")
surface_area_raw <- read_csv("~/Documents/GitHub/Plan_C/data/raw/physiology/Nubbin_sizing_notes.csv")

doc_exclude <- read_xlsx("~/Documents/GitHub/Plan_C/data/raw/doc/DOC_exlude_variance.xlsx")$`Sample Code`%>%
  as.vector()

name_samples <- function(x) { 
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

pc_doc_raw <- read_xlsx("~/Documents/GitHub/Plan_C/data/raw/doc/PlanC_DOC_upload.xlsx")%>%
  filter(!`Sample Code` %in% doc_exclude)%>%
  rename(Code = `Sample Code`)%>%
  name_samples()
  
  
sym_raw <- read_csv('~/Documents/GitHub/Plan_C/data/raw/physiology/PLANC_Blastate_FCM.csv')
sym_names <- read_csv('~/Documents/GitHub/Plan_C/data/raw/physiology/PLANC_Blastate_FCM_Sample_Map.csv')%>%
  mutate(Plate = gsub('0', '', Plate))

# CLEANING -- FCM ---------------------------------------------------------
fcm_raw <- fcm_names%>%
  map_df( ~as.data.frame(.), .id="plate")%>%
  mutate(plate = gsub(".Map.xlsx", "", plate),
         plate = gsub("OCT.2019", "OCT.18.2019", plate))%>%
  full_join(., fcm_files_raw, by = c("plate", "Well"))%>%
  separate(plate, c("experiment_plate", "dye", "month", "day", "year", "plate", "number"))%>%
  unite(date, c("month", "day", "year"), sep = "-")%>%
  unite(plate_numnber, c("plate", "number"), sep = "_")%>%
  # filter(!Sample %in% doc_exclude)%>%
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

fcm_sybr <- fcm_raw%>%
  filter(Experiment == "PC",
         dye == "SYBR",
         `Gate Name` == "SYBR Polygon")%>%
  # filter(bath == "hot")%>%
  group_by(Day, Health, Temperature)%>%
  mutate(cells_mL = `Event Count`/0.075,
         cells_ul = cells_mL/1000,
         log10_ul = log10(cells_ul),
         std = sd(log10_ul))

# fcm_means <- fcm_sybr%>%
#   select(-c(1:5, replicate))%>%
#   group_by(Day, health, bath)%>%
#   summarize_if(is.numeric, mean)%>%
#   ungroup()




# CLEANING -- Symbiodinium counts -----------------------------------------
sym_bad <- c('PC_BL2_T0C_A',
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

sym_stats <- sym_names%>%
  left_join(sym_raw%>%
              mutate(Plate = as.character(Plate)), by = c('Plate', 'Well'))%>%
  mutate(count = case_when(is.na(Notes) ~ `Sym-FSC Events`,
                              !is.na(Notes) ~ `Sym-SSC Events`))%>%
  select(Sample, count)%>%
  rename(Code = Sample)%>%
  filter(!Code %in% sym_bad)%>%
  left_join(surface_area_raw%>%
              select(Code, SA, Notes),
            by = 'Code')%>%
  mutate(sym_sa = count/SA)%>%
  name_samples()%>%
  filter(!Code %like% 'mil%')%>%
  separate(Code, c('Code', 'nubbin'), sep = -1)%>%
  mutate(log = log10(sym_sa))

  


# GRAPHING -- Color Pallette and theme----------------------------------------------
health_colors <- c("#FD6467", "#E1AF00", "#0B775E","#78B7C5")

gen_theme <- function(x) {
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

# GRAPHING -- FCM ---------------------------------------------------------
pdf("~/Documents/GitHub/Plan_C/data/plots/Plan_c_FCM.pdf", width = 15, height = 10)
# Ambient for presentations
fcm_sybr%>%
  unite(cb, c(health, bath), sep = "_", remove = FALSE)%>%
  mutate(health = as.factor(health),
         health = fct_relevel(health, "Bleached", "Partially bleached"))%>%
  rename(`Initial Health` = health,
         Temperature = bath)%>%
  filter(Temperature == 'Ambient')%>%
  ggplot(aes(x = Day, y = cells_ul, color = `Initial Health`)) +
  geom_point(aes(group = cb), stat = "identity", size = 2) + 
  # geom_errorbar(aes(ymin = mean(cells_mL) - std, ymax = mean(cells_mL) + std)) +
  geom_line(aes(group = cb, linetype = Temperature), stat = "summary", fun.y = "mean", lwd = 1.4) +
  ylab("Cells µL-1") +
  scale_linetype_manual(values=c("dashed", "solid")) +
  scale_color_manual(values = health_colors) +
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

fcm_sybr%>%
  unite(cb, c(health, bath), sep = "_", remove = FALSE)%>%
  mutate(health = as.factor(health),
         health = fct_relevel(health, "Bleached", "Partially bleached"))%>%
  rename(`Initial Health` = health,
         Temperature = bath)%>%
  filter(Temperature != 'Ambient')%>%
  ggplot(aes(x = Day, y = cells_ul, color = `Initial Health`)) +
  geom_point(aes(group = cb), stat = "identity", size = 2) + 
  # geom_errorbar(aes(ymin = mean(cells_mL) - std, ymax = mean(cells_mL) + std)) +
  geom_line(aes(group = cb, linetype = Temperature), stat = "summary", fun.y = "mean", lwd = 1.4) +
  ylab("Cells µL-1") +
  scale_linetype_manual(values=c("solid")) +
  scale_color_manual(values = health_colors) +
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
dev.off()

pdf("~/Documents/GitHub/Plan_C/data/plots/Plan_c_FCM.pdf", width = 15, height = 10)
fcm_means%>%
  filter(health != 'bleached')%>%
  rename(`Starting Health` = health)%>%
  ggplot(aes(x = Day, y = , color = `Starting Health`))+
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
sym_stats%>%
  filter(!Code %like% '%T7C%')%>%
  unite(dayHealth, c('Day', 'Health'), sep = "_", remove = FALSE)%>%
  mutate(Day = as.factor(Day))%>% 
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  ggplot(aes(x = Day, y = log, color = Health)) + 
  geom_boxplot(aes(group = dayHealth)) +
  scale_color_manual(values = health_colors) +
  facet_wrap(~ Health) +
  labs(y = bquote(Log[10] ~(Symbiodinium ~cm^-2)), title = 'Symbiodinium counts heated treatment') +
  gen_theme()

sym_stats%>%
  filter(!Code %like% '%T7H%')%>%
  unite(dayHealth, c('Day', 'Health'), sep = "_", remove = FALSE)%>%
  mutate(Day = as.factor(Day))%>% 
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  ggplot(aes(x = Day, y = log, color = Health)) + 
  geom_boxplot(aes(group = dayHealth)) +
  scale_color_manual(values = health_colors) +
  facet_wrap(~ Health) +
  labs(y = bquote(Log[10] ~(Symbiodinium ~cm^-2)), title = 'Symbiodinium counts Ambient treatment') +
  gen_theme()
dev.off()


# GRAPHING -- DOC ---------------------------------------------------------
pdf("~/Documents/GitHub/Plan_C/data/plots/Plan_c_DOC.pdf", width = 15, height = 10)
pc_doc_raw%>%
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  unite(cb, c(Health, Temperature), sep = "_", remove = FALSE)%>%
  rename(`Initial Health` = Health)%>%
  filter(Temperature == 'Ambient')%>%
  ggplot(aes(Day, `DOC (µMC)`, color = `Initial Health`)) +
  geom_point(stat = 'identity') +
  geom_line(aes(group = cb, linetype = Temperature), stat = "summary", fun.y = "mean", lwd = 1.4) +
  scale_linetype_manual(values=c("dashed", "solid")) +
  scale_color_manual(values = health_colors) +
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


pc_doc_raw%>%
  mutate(Health = as.factor(Health),
         Health = fct_relevel(Health, "Bleached", "Partially bleached"))%>%
  unite(cb, c(Health, Temperature), sep = "_", remove = FALSE)%>%
  rename(`Initial Health` = Health)%>%
  filter(Temperature != 'Ambient')%>%
  ggplot(aes(Day, `DOC (µMC)`, color = `Initial Health`)) +
  geom_point(stat = 'identity') +
  geom_line(aes(group = cb, linetype = Temperature), stat = "summary", fun.y = "mean", lwd = 1.4) +
  scale_linetype_manual(values=c("solid")) +
  scale_color_manual(values = health_colors) +
  scale_x_discrete(limits = c(1,2,3,4,5,6,7)) +
  labs(x = 'Day', y = "Dissolved Orgnic Carbon (µMC)") +
  ylim(65,176) +
  gen_them()
dev.off()

# GRAPHONG -- combine DOC FCM bad -----------------------------------------------------
fcm_doc <- fcm_sybr%>%
  full_join(pc_doc_raw%>%
              rename(Sample = `Sample Code`), by = c('Sample', 'Health', 'Day', 'Temperature'))

fcm_doc%>%
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
  sec.axis = sec_axis(~.*100, name= bquote(Cells ~µL^-1)))

