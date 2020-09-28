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

# # SETTING CORES -- for parrelelizing code ---------------------------------
# # sets the number of cores to one fewer than the max available
# num_cores <- availableCores() -1
# # actually sets the planning algoreithm for how the code will be mapped across the cores.
# plan(multiprocess, workers = num_cores) 


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

# CLEANING -- FCM ---------------------------------------------------------
fcm_raw <- fcm_names%>%
  map_df( ~as.data.frame(.), .id="plate")%>%
  mutate(plate = gsub(".Map.xlsx", "", plate),
         plate = gsub("OCT.2019", "OCT.18.2019", plate))%>%
  full_join(., fcm_files_raw, by = c("plate", "Well"))%>%
  separate(plate, c("experiment_plate", "dye", "month", "day", "year", "plate", "number"))%>%
  unite(date, c("month", "day", "year"), sep = "-")%>%
  unite(plate_numnber, c("plate", "number"), sep = "_")%>%
  separate(Sample, c("Experiment", "health", "Timepoint"), sep = "_")%>%
  separate(health, c("health", "replicate"), sep = 2)%>%
  separate(Timepoint, c("timepoint", "bath"), sep = 2)%>%
  mutate(health = case_when(health == "HE" ~ "healthy",
                            health == "PB" ~ "partially bleached",
                            health == "BL" ~ "bleached",
                            health == "WA" ~ "water control",
                            TRUE ~ as.character(health)),
         timepoint = gsub("T", "", timepoint),
         bath = case_when(bath == "C" ~ "control",
                          bath == "H" ~ "hot",
                          TRUE ~ as.character(bath)))

fcm_sybr <- fcm_raw%>%
  filter(Experiment == "PC",
         dye == "SYBR",
         `Gate Name` == "SYBR Polygon")%>%
  mutate(cells_mL = `Event Count`/0.075,
         cells_ul = cells_mL/1000)%>%
  group_by(timepoint, health, bath)%>%
  mutate(std = sd(cells_mL))%>%
  rename('Day' = timepoint)%>%
  filter(bath == "hot")%>%
  ungroup()



# GRAPHING -- FCM ---------------------------------------------------------
pdf("~/Documents/GitHub/Plan_C/data/plots/Plan_c_FCM.pdf", width = 15, height = 10)
fcm_sybr%>%
  ggplot(aes(x = Day, y = cells_ul, color = health))+
  geom_point(stat = "summary", fun.y = "mean") + 
  # geom_errorbar(aes(ymin = mean(cells_mL) - std, ymax = mean(cells_mL) + std)) +
  geom_line(aes(group = health), stat = "summary", fun.y = "mean") +
  theme(
    legend.position = "top",
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"))
  # facet_wrap(~bath)
dev.off()





