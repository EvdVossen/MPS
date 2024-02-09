rm(list=ls())

library(dplyr)

path_data <- "~/Data_files/MPS/Eduard/Data_transfer"
setwd(path_data)

#SGB data from metaphlan
sp <- rio::import("Pipeline_output/Metaphlan_merged_abundance_table_SGB.txt") %>% 
  tibble::column_to_rownames("clade_name")

#getting the triads
pd <- rio::import("patient_donor_mapping_w_febaligo.xlsx") %>%
  dplyr::mutate(Pre_FMT = gsub("MP-M", "MPS_M", .$Pt_ID),
                Post_FMT = gsub("MP-M","MPP_M", .$Pt_ID),
                Donor_Sample_ID = gsub("MP-H", "MPS_H", .$Donor_Sample_ID)) %>% 
  dplyr::select(Study_origin, Pt_ID, Donor_Sample_ID, Pre_FMT, Post_FMT, Intervention) %>% 
  dplyr::filter(!grepl("NA",.$Donor_Sample_ID),
                !grepl("FATLOSE2",.$Study_origin),
                !.$Pt_ID=="MP-M11",
                .$Intervention == 1) %>% 
  dplyr::mutate(Triad = paste0("Triad_",rownames(.)))

species <- sp[,names(sp) %in% c(unique(pd$Donor_Sample_ID), pd$Pre_FMT, pd$Post_FMT)] %>%
  t(.) %>%
  as.data.frame(.)

#extracting the data
pd.long <- pd %>% 
  dplyr::select(-Triad, -Pt_ID) %>% 
  tidyr::gather("Triad", "Subject_ID",2:4) %>% 
  dplyr::mutate(Subject = gsub("MPP_","",gsub("MPS_","", Subject_ID)),
                Study_origin = ifelse(Triad == "Donor_Sample_ID", "", Study_origin)) %>% 
  distinct(Study_origin,Triad,Subject_ID,Subject) %>% 
  mutate(Index = paste0("MP-",.$Subject))

dat <- as.data.frame(rownames(species)) %>%
  setNames("Subject_ID") %>%
  mutate(Subject = gsub("MPP_","", gsub("MPS_","", .$Subject_ID))) %>%
  mutate(Triad = ifelse(grepl("MPP",.$Subject_ID), "Post-FMT", ifelse(grepl("MPS_H", .$Subject_ID), "Donor", "Pre-FMT"))) %>%
  mutate(Species = sp$SGB385[match(.$Subject_ID, rownames(sp))]) %>%
  mutate(Index = paste0("MP-",.$Subject))

dir.create(paste0(path_data, "/Intermediate_files/"))
rio::export(pd.long, paste0(path_data, "/Intermediate_files/", "Triad_data_long.xlsx"))
rio::export(pd, paste0(path_data,"/Intermediate_files/", "Triad_numbers.xlsx"))
rio::export(species, paste0(path_data,"/Intermediate_files/", "Triads_MP4.xlsx"), rowNames=T)
