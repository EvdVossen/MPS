rm(list=ls())

library(tidyverse)
library(ggplot2)
source("~/Data_files/MPS/Eduard/Data_transfer/R-scripts/functions.R")

setwd(path_data)

#import files (imports are done in the functions script)
d <- get_data(st = T, triad = T, triad_long = T, met=T, mp4 = T, path_dat = T,  metab = T)
base::list2env(d,envir=.GlobalEnv);

#complete cases; otherwise skewed results on the pca
met_pre <- met_pre %>% select_if(~ !any(is.na(.))) 
met_post <- met_post %>% select_if(~ !any(is.na(.)))

## check if all rownames are aligned
all(c(
  rownames(mp4_pre) == rownames(mp4_post),
  rownames(p_pre) == rownames(p_post),
  rownames(m_pre) == rownames(m_post),
  rownames(met_pre) == rownames(met_post),
  rownames(mp4_pre) == rownames(m_post),
  rownames(mp4_pre) == rownames(p_post),
  rownames(mp4_pre) == rownames(met_post)
))

## Lists for the loops.
dat_list <- list(mp4_pre, p_pre, m_pre, met_pre,
                 mp4_post, p_post, m_post, met_post)
dist_list <- vector(mode = "list", length = length(dat_list))
names(dist_list)  <- c("mp4_pre", "p_pre", "m_pre", "met_pre",
                       "mp4_post", "p_post", "m_post", "met_post")
pcoa_list <- dist_list
vector_list <- pcoa_list

for(i in 1:length(dist_list)){
  dist_list[[i]] <- vegan::vegdist(as.data.frame(compositions::clr(dat_list[[i]])), method = "euclidean")
  pcoa_list[[i]] <- ape::pcoa(as.matrix(dist_list[[i]]))
  vector_list[[i]] <- pcoa_list[[i]]$vectors
}

## All combinations of a list of 8 dataframes
all_combinations <- utils::combn(length(vector_list), 2)
proc_list <- vector(mode= "list", length = ncol(all_combinations))
protest_list <- proc_list
pl <- proc_list

names_plot <- data.frame(names(dist_list), 
                         c("Metaphlan pre-FMT","Pathways pre-FMT","Metabolites pre-FMT","Metadata Pre-FMT",
                           "Metaphlan post-FMT","Pathways post-FMT","Metabolites post-FMT","Metadata Post-FMT")) %>% 
  setNames(c("short_name", "full_name"))
# dir.create("plots/Procrustes/")

df_sup <- data.frame(Dataset_1 = numeric(0), 
                     Dataset_2 = numeric(0), 
                     Sum_of_squares = numeric(0), 
                     rho = numeric(0), 
                     p.value = numeric(0))

for (i in 1:ncol(all_combinations)) {
  set.seed(1)
  comb <- all_combinations[,i]
  name_comb1 <- names(vector_list)[comb[1]]
  name_comb2 <- names(vector_list)[comb[2]]
  name1_plot <- names_plot[names_plot$short_name %in% name_comb1,]$full_name
  name2_plot <- names_plot[names_plot$short_name %in% name_comb2,]$full_name
  proc_list[[i]] <- vegan::procrustes(X = vector_list[[name_comb1]], Y = vector_list[[name_comb2]], symmetric=TRUE)
  
  protest_list[[i]] <- vegan::protest(X = vector_list[[name_comb1]], Y = vector_list[[name_comb2]], permutations = 999)
  
  proc.df <- rbind(
    data.frame(PC1=proc_list[[i]]$X[,1],PC2=proc_list[[i]]$X[,2], 
               type=name1_plot,
               Sample.Group=triad$Study_origin[match(rownames(vector_list[[name_comb1]]),gsub("MP-","",triad$Pt_ID))],
               Sample.ID=rownames(vector_list[[name_comb1]]), stringsAsFactors=F),
    data.frame(PC1=proc_list[[i]]$Yrot[,1],PC2=proc_list[[i]]$Yrot[,2], 
               type=name2_plot, 
               Sample.Group=triad$Study_origin[match(rownames(vector_list[[name_comb2]]),gsub("MP-","",triad$Pt_ID))], 
               Sample.ID=rownames(vector_list[[name_comb2]]), stringsAsFactors=F))
  
  pl[[i]] <-   ggplot(proc.df) +
    geom_point(alpha=.7, aes(x = PC1, y = PC2, color = Sample.Group, shape=type), size = 3) +
    geom_line(aes(x = PC1, y = PC2, group=Sample.ID, color=Sample.Group), alpha=1) +
    ggrepel::geom_label_repel(aes(x = PC1, y = PC2, label = gsub("MPP_","",Sample.ID)))+
    xlab("PC1") + ylab("PC2") + 
    theme(legend.title = element_blank()) +
    labs(x = "PC1", y = "PC2",
         title=paste0("Combination: ", name1_plot, " with ", name2_plot),
         subtitle=paste0("sum of squares: ", round(protest_list[[i]]$ss,2), "; r = ", round(protest_list[[i]]$t0,3), "; p = ", protest_list[[i]]$signif))
  
  ggsave(filename = paste0("cc_procrustes_", name1_plot, "_", name2_plot,".png"),plot = print(pl[[i]]),device = "png", path = "plots/Procrustes/",height = 10, width = 15)

  df_sup[i, ] <- c(name1_plot,
                   name2_plot,
                   round(protest_list[[i]]$ss,2),
                   round(protest_list[[i]]$t0,3), 
                   protest_list[[i]]$signif)  
}