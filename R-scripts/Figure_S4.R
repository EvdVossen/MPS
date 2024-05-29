rm(list=ls())

library(tidyverse)
library(ggplot2)
library(phyloseq)
library(vegan)
library(mixOmics)
library(ggplot2)
library(ggpubr)
source("~/Data_files/MPS/Eduard/Data_transfer/R-scripts/functions.R")

setwd(path_data)

#import files (imports are done in the functions script)
d <- get_data(triad = T, triad_long = T, met=T, metab = T)
base::list2env(d,envir=.GlobalEnv)

# Intermediate file on global strainsharing information
df_fractionfig <- rio::import("Intermediate_files/Strainsharing_global_information.xlsx")
  

############ Mark Multilevel PCA - metabolites
metab_pca <- metab %>% 
  mutate(Subject_ID = triad_long$Index[match(rownames(.), triad_long$Subject_ID)]) %>% 
  mutate(Time = factor(triad_long$Triad[match(rownames(.), triad_long$Subject_ID)], levels=c("Pre_FMT", "Post_FMT")))

#clr-scaling
clrdat <- as.data.frame(compositions::clr(as.data.frame((metab))))

#PCA
pca.result <- mixOmics::pca(clrdat, scale = TRUE, multilevel = as.factor(metab_pca$Subject_ID))

#Creating df with visit
df <- data.frame(pca.result$variates) %>% 
  dplyr::mutate(Visit = metab_pca$Time[match(rownames(.), rownames(metab_pca))],
                sample = rownames(.)) %>%
  dplyr::mutate(Subject_ID = gsub("MPS_", "MP-", gsub("MPP_", "MP-", sample))) %>% 
  merge(., aggregate(.$X.PC1 ~ Visit, data = ., mean), by="Visit") %>%   ### include centroids
  merge(., aggregate(.$X.PC2 ~ Visit, data = ., mean), by="Visit")  %>% ### include centroids
  setNames(c("Visit", "X","Y", "Sample", "Subject_ID", "X.centroid", "Y.centroid")) %>% 
  dplyr::mutate(Donor_sharing_rate = ifelse(.$Visit == "Post_FMT", df_fractionfig$Donor_sharing_rate[match(.$Sample, df_fractionfig$Post_FMT)], 0))


#adding envfit for diastolic bloodpressure
met_prepost <- rbind(met_pre %>% dplyr::mutate(Subject_ID = paste0("MPS_", rownames(.))), 
                     met_post %>% dplyr::mutate(Subject_ID = paste0("MPP_", rownames(.)))) %>%
  arrange(match(.$Subject_ID, rownames(metab_pca))) %>% 
  tibble::rownames_to_column(var = "RowID") %>%  # Convert rownames to a column named "RowID"
  dplyr::select(-RowID) %>% 
  tibble::column_to_rownames("Subject_ID") %>% 
  dplyr::mutate(donor_sharing_rate = df_fractionfig$Donor_sharing_rate[match(rownames(.), df_fractionfig$Post_FMT)])
met_prepost_diast <- met_prepost %>% 
  dplyr::filter(!is.na(Diast))

#Envfit
set.seed(0)
en_diast <- envfit(pca.result, met_prepost %>% dplyr::select(Diast), na.rm=T) #contains NA
en_rd <- envfit(pca.result, met_prepost %>% dplyr::select(Rd), na.rm=F) #contains no NA

#figure attributes
p <- ggplot(data=df, 
       aes_string(x = df$X,
                  y = df$Y)) +
  geom_point(aes(shape = Visit, 
                 fill = Donor_sharing_rate), 
             size = 3, 
             stroke = .1) + 
  geom_segment(aes(x = 0, y = 
                     0, 
                   xend = en_diast$vectors$arrows[1]*10, 
                   yend = en_diast$vectors$arrows[2]*10),
               arrow = arrow(length = unit(0.1, "inches")), 
               color = "grey60") +
  geom_segment(aes(x = 0, y = 
                     0, 
                   xend = en_rd$vectors$arrows[1]*10, 
                   yend = en_rd$vectors$arrows[2]*10),
               arrow = arrow(length = unit(0.1, "inches")), 
               color = "grey60") +
  # ggrepel::geom_label_repel(aes(x=df$X, y=df$Y, label = gsub("MP-","",df$Subject_ID))) +
  scale_shape_manual(values = c(21, 24)) +
  geom_segment(aes_string(x=df$X.centroid, 
                          y=df$Y.centroid, 
                          xend=df$X, 
                          yend=df$Y), 
               alpha=0.3, 
               color="black") +
  geom_text(aes(x = en_diast$vectors$arrows[1]*10, y = en_diast$vectors$arrows[2]*10, label = "Diast*"), vjust = -0.5, color = "black") +
  geom_text(aes(x = en_rd$vectors$arrows[1]*10, y = en_rd$vectors$arrows[2]*10, label = "Rd"), vjust = -0.5, color = "black") +
  theme_bw() +
  labs(x=paste0("PC1: ",round(pca.result$prop_expl_var$X[1], 3)*100, "% variance explained"),
       y=paste0("PC2: ",round(pca.result$prop_expl_var$X[2], 3)*100, "% variance explained"),
       fill = "Fraction of strains  \nshared with the donor") +
  # facet_wrap(~Visit) +
  theme(panel.spacing.x = unit(6, "mm")) +
  scale_fill_gradient(low = "black", high = "red")

ggsave(filename = "Supplementary_Figure_S4.pdf", 
       plot = p, 
       device = "pdf", 
       path = "Manuscript/Supplementary_information/", 
       width = 14, height = 8)
