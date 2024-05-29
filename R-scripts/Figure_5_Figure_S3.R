rm(list=ls())

#Load packages and functions
source("~/Data_files/MPS/Eduard/Data_transfer/R-scripts/functions.R")

#Set working directory
setwd(path_data)

#import files (imports are done in the functions script)
d <- get_data(st = T, triad = T, triad_long = T, met=T, mp4 = T, metab = T)
base::list2env(d,envir=.GlobalEnv)

#Imports + standard dataprocessing (somewhat sub-optimal programming)
top_list_engraftment <- rio::import(paste0(path_data, "/Intermediate_files/top_list_strain_engraftment.xlsx")) %>% 
  tibble::column_to_rownames("...1")

#strain binary
st1 <- st %>% dplyr::filter(paste0(X1,X2) %in% paste0(triad$Post_FMT, triad$Donor_Sample_ID)) %>%
  dplyr::select(X1, X2, all_of(rownames(top_list_engraftment))) %>%
  dplyr::mutate_at(vars(-X1, -X2), ~ifelse(is.na(.) | . == "" | !., 0, 1)) %>%
  dplyr::mutate(X1 = gsub("MPP_","", X1)) %>%
  dplyr::select(-X2) %>%
  tibble::column_to_rownames("X1") %>%
  .[order(as.numeric(gsub("M","", rownames(.)))),] %>%
  setNames(top_list_engraftment$plot_name[match(names(.), rownames(top_list_engraftment))])

#Fix metadata
met_post_comp <- met_post %>% 
  dplyr::select(-Study_origin, -Weight)
met_pre_comp <- met_pre[rownames(met_pre) %in% rownames(met_post_comp), names(met_pre) %in% names(met_post_comp),]

met_delta <- met_post_comp - met_pre_comp
met_delta_comp_diast <- met_delta %>% 
  filter(!is.na(Diast)) %>% 
  select_if(~ !any(is.na(.))) %>% 
  dplyr::select(c(Rd, Diast,EGPsupp, Insulin, HOMA))

# Add publication names & threshold levels according to literature (or clinicians)
df_names_met <- as.data.frame(names(met_delta)) %>% 
  stats::setNames("Abbreviation") %>%
  dplyr::mutate(Pub_names = c("BMI (kg/m^2)", "Systolic blood pressure (mm Hg)", "Diastolic blood pressure (mm Hg)", 
                              "Glucose levels (mmol/L)", "Insulin levels (pmol/L)", "HbA1c levels (mmol/mol)", "HOMA-IR",
                              "Total cholesterol (mmol/L)", "LDL cholesterol (mmol/L)", "HDL cholesterol (mmol/L)", 
                              "Triglyceride levels (mmol/L)", "ALAT (units/L)", "C-reactive protein (CRP) levels (mmol/L)",
                              "Rate of glucose disappearance (μmol/kg/min)", "Suppresion of endogenous glucose production (μmol/kg/min)",
                              "Resting energy expenditure (kcal/day"),
                Threshold_levels = c(25, 120, 80,
                                     5.6, 174, 42, 1.9,
                                     5.17, 2.6, 1.55,
                                     1.47, 56, 10,
                                     37.3, 46.5,
                                     NA))

#colors
colors <- c("#440154FF","#21908CFF", "#35B779FF")
date <- Sys.Date()

#metadata_delta x strain_bin
st_bin_diast <- st1 %>% dplyr::filter(rownames(.) %in% rownames(met_delta_comp_diast))
cca_res <- CCA_fun(X = met_delta_comp_diast, Y = st_bin_diast, 
                   x_data_name = "Metadata", y_data_name = "Strain engraftment",
                   color_x_data = colors[1], color_y_data = colors[2],
                   grid1 = seq(0.1, 1, length = 5),
                   grid2 = seq(0.1, 1, length = 5),
                   n_comp = 10,validation = "loo",
                   name_map = paste0("diast_", date),
                   name_submap = "delta_met_strainsbin_no_na_variables_fin",
                   cor.method = "spearman",
                   output_path = paste0(path_data, "/plots/CCA"))

### Downstream for LMER
all(rownames(met_pre) == rownames(met_post))

met_post_comp <- met_post %>% 
  dplyr::select(-c(Study_origin, Weight)) %>%
  dplyr::mutate(Time = "Post", Subject_ID = rownames(.))
met_pre_comp <- met_pre[rownames(met_pre) %in% rownames(met_post_comp), names(met_pre) %in% names(met_post_comp),] %>% 
  dplyr::mutate(Time = "Pre", Subject_ID = rownames(.))

met_full_comp <- rbind(met_post_comp, met_pre_comp)
rownames(met_full_comp)<-NULL

pairs_oi <- cca_res$pairs_oi %>% head(10)

#linear mixed effect models
lmer_df <- data.frame(
  var_dataset1 = numeric(0),
  var_dataset2 = numeric(0),
  estimate = numeric(0),
  `p-value` = numeric(0))

comb_for_lmer <- list()
for (idx in 1:nrow(pairs_oi)){
  name_clinparam <- pairs_oi$var_dataset1[idx]
  name_strain <- pairs_oi$var_dataset2[idx]
  comb_for_lmer[[idx]] <- merge(st1 %>% dplyr::select(all_of(name_strain)), met_full_comp %>% dplyr::select(all_of(name_clinparam), Time, Subject_ID), 
                         by.x = "row.names", by.y = "Subject_ID") %>% 
    dplyr::rename(Subject_ID = Row.names) %>% 
    dplyr::mutate(Time = ifelse(Time == "Pre","Pre FMT", "Post FMT")) %>% 
    dplyr::mutate(Time = factor(Time, levels = c("Pre FMT", "Post FMT"))) 
  lmer_test <- summary(lmerTest::lmer(comb_for_lmer[[idx]][,name_clinparam] ~ comb_for_lmer[[idx]][,name_strain] * Time + (1|Subject_ID), data = comb_for_lmer[[idx]]))
  
  lmer_df[idx,] <- c(name_clinparam, 
                     name_strain, 
                     lmer_test$coefficients[4], 
                     lmer_test$coefficients[20]) 
}

#Results lmer strain engraftment and clinical parameters
lmer_df$p.adj <- p.adjust(lmer_df$p.value, method = "fdr")
print(lmer_df)

#Supplementary Figure S3
print(cca_res$p_s3)

#Figure 5
print(cca_res$pf5)
