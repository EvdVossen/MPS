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

#strain engraftment
st1 <- st %>% dplyr::filter(paste0(X1,X2) %in% paste0(triad$Post_FMT, triad$Donor_Sample_ID)) %>%
  dplyr::select(X1, X2, all_of(rownames(top_list_engraftment))) %>%
  dplyr::mutate_at(vars(-X1, -X2), ~ifelse(is.na(.) | . == "" | !., 0, 1)) %>%
  dplyr::mutate(X1 = gsub("MPP_","", X1)) %>%
  dplyr::select(-X2) %>%
  tibble::column_to_rownames("X1") %>%
  .[order(as.numeric(gsub("M","", rownames(.)))),] %>%
  setNames(top_list_engraftment$plot_name[match(names(.), rownames(top_list_engraftment))]) %>% 
  dplyr::select("C. aerofaciens (SGB14546 group)", "F. saccharivorans (SGB4874)")

#metab_delta
m_delta <- m_post - m_pre

#2-oxoarginine comes out -> within the ADH pathway(?); 4-guanidinobutanoate comes after (https://www.genome.jp/pathway/map00330+C03771)
metab_for_lmer <- metab %>% 
  dplyr::mutate(Time = ifelse(grepl("MPP", rownames(.)), "Post FMT", "Pre FMT"),
                Subject_ID = gsub("MPS_","",gsub("MPP_","", rownames(.)))) %>% 
  dplyr::mutate(subject_ID_time = paste0(Subject_ID, "_", Time)) %>% 
  dplyr::select(Time, Subject_ID, subject_ID_time, `4-guanidinobutanoate`, `2-oxoarginine*`)

## For LMER
met_post_comp <- met_post %>% 
  dplyr::select(-c(Study_origin, Weight)) %>%
  dplyr::mutate(Time = "Post FMT", Subject_ID = rownames(.))
met_pre_comp <- met_pre[rownames(met_pre) %in% rownames(met_post_comp), names(met_pre) %in% names(met_post_comp),] %>% 
  dplyr::mutate(Time = "Pre FMT", Subject_ID = rownames(.))

met_full_comp <- rbind(met_post_comp, met_pre_comp) %>% 
  `rownames<-`(NULL)
lmer_dat_diast <- merge(st1, met_full_comp %>% dplyr::select(Diast, Time, Subject_ID), 
                              by.x = "row.names", by.y = "Subject_ID") %>% 
  dplyr::rename(Subject_ID = Row.names) %>% 
  dplyr::mutate(subject_ID_time = paste0(Subject_ID, "_", Time)) %>% 
  dplyr::select(Subject_ID, subject_ID_time, Time, Diast, everything())

p_ox <- list()
p_gua <- list()
for (k in 5:ncol(lmer_dat_diast)){
  n_p <- k-4
  var_name_st <- names(lmer_dat_diast)[k]
  
  for(l in 4:ncol(metab_for_lmer)){
    var_name_metab <- names(metab_for_lmer)[l]
    comb_for_lmer <- merge(metab_for_lmer %>% dplyr::select(all_of(var_name_metab), subject_ID_time, Time, Subject_ID), lmer_dat_diast %>% dplyr::select(all_of(var_name_st), subject_ID_time), by = "subject_ID_time") %>% 
      mutate(Time = factor(Time, levels = c("Pre FMT", "Post FMT")))

    lmer_test <- summary(lmerTest::lmer(comb_for_lmer[,var_name_metab] ~ comb_for_lmer[,var_name_st] * Time + (1|Subject_ID), data = comb_for_lmer))
    
    # for the boxplot
    comb_for_plot <- comb_for_lmer %>% 
      dplyr::select(Subject_ID, all_of(var_name_metab), all_of(var_name_st), Time) %>% 
      setNames(c("Subject_ID","var_name_metab", "var_name_st", "Time")) %>% 
      dplyr::mutate(Strain_engraftment = ifelse(var_name_st == 0, "No engraftment", "Engraftment")) %>%
      dplyr::mutate(Strain_engraftment = factor(Strain_engraftment, levels = c("No engraftment","Engraftment"))) %>% 
      dplyr::mutate(Group = paste0(Time, " ", Strain_engraftment)) %>% 
      dplyr::mutate(Group = factor(Group, levels = c("Pre FMT No engraftment", "Post FMT No engraftment", "Pre FMT Engraftment", "Post FMT Engraftment")))
    
    # for the barplot
    df2 <- merge(m_delta %>% dplyr::select(all_of(var_name_metab)), st1 %>% dplyr::select(all_of(var_name_st)), by = "row.names") %>% 
      tibble::column_to_rownames("Row.names") %>% 
      setNames(c("var_name_metab", "var_name_st")) %>% 
      data_summary(., varname="var_name_metab", 
                   groupnames=c("var_name_st")) %>% 
      dplyr::mutate(time = c("Delta")) %>%  
      dplyr::rename("delta_var" = median) %>% 
      dplyr::mutate(Strain_engraftment = ifelse(var_name_st==0, "No engraftment", "Engraftment")) %>% 
      dplyr::mutate(Strain_engraftment = factor(Strain_engraftment, levels = c("No engraftment", "Engraftment"))) %>% 
      dplyr::mutate(loc_plot = c(1.5, 3.5))
    
    # combination plot of boxplot with barplot with p-value located at the lines
  p <- ggplot() +
    geom_boxplot(data = comb_for_plot,
                 aes(x = Group, 
                     y = var_name_metab, 
                     fill = Strain_engraftment),
                 alpha = 0.5,
                 width = 0.5,
                 show.legend = F) +
    geom_point(data = comb_for_plot,
               aes(x = Group, 
                   y = var_name_metab, 
                   col = Strain_engraftment),
               size = 1.5, color = "black") +
    geom_line(data = comb_for_plot,
              aes(x = Group, 
                  y = var_name_metab, 
                  col = Strain_engraftment, 
                  group = Subject_ID),
              color = "gray30") +
    geom_hline(yintercept = 0) +
    geom_bar(data = df2,
             aes(x = loc_plot,
                 y = delta_var,
                 fill = Strain_engraftment),
             stat = "identity",
             alpha = 0.5,
             color = "black", 
             width = 0.5) +
    geom_errorbar(data = df2,
                  aes(x = loc_plot,
                      y = delta_var,
                      xmin = loc_plot, 
                      xmax = loc_plot,
                      ymax = delta_var + (delta_var > 0)*sd,
                      ymin = delta_var - (delta_var < 0)*sd),
                  position = position_dodge(), 
                  width = 0.2)  +
    stat_pvalue_manual(data = data.frame(group1 = 1.5, group2 = 3.5, p.adj = round(as.numeric(lmer_test$coefficients[20]),3)),
                       label = "p.adj",
                       y.position = max(comb_for_plot$var_name_metab, na.rm = T)*1.05,
                       tip.length = .01) +
    ylab(label = var_name_metab) +
    xlab(label = "Time") +
    scale_x_discrete(labels = c(levels(comb_for_plot$Time), levels(comb_for_plot$Time))) +
    theme_Publication() +
    scale_fill_manual(name = var_name_st, labels = c("No engraftment", "Engraftment"),
                      values=c("firebrick1", "dodgerblue1","firebrick1", "dodgerblue1")) +
    theme(legend.text.align = 0,
          strip.background =  element_rect(fill = NA, 
                                           colour = NA),
          plot.margin = unit(c(1, 1, 1, 1), "cm"))
  if(l == 4){p_gua[[n_p]] <- p}
  if(l == 5){p_ox[[n_p]] <- p + theme(legend.position ="none")}
  }
}

library(patchwork)
p_s6 <- p_ox[[1]] + p_gua[[1]] +
  p_ox[[2]] + p_gua[[2]]

ggsave(filename = "Supplementary_Figure_S6.pdf", plot = p_s6, device = "pdf", path = "Manuscript/Supplementary_information/", height = 9, width = 12)
