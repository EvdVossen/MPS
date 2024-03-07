rm(list=ls())

library(tidyverse)
library(igraph)
library(ggplot2)
source("~/Data_files/MPS/Eduard/Data_transfer/R-scripts/functions.R")

setwd(path_data)

#import files (imports are done in the functions script)
d <- get_data(st = T, triad = T, triad_long = T, met=T, mp4 = T, metab = F)
base::list2env(d,envir=.GlobalEnv)

#### per study
fatmed <- triad[triad$Study_origin=="FATMED",]
appetite <- triad[triad$Study_origin=="APPETITE" & triad$Intervention==1,]
febaligo <- triad[triad$Study_origin=="FEBALIGO",]

fatmed_met <- fatmed %>% 
  mutate(donor_name_plot=LETTERS[match(.$Donor_Sample_ID, unique(.$Donor_Sample_ID))]) %>% 
  group_by(donor_name_plot) %>%
  mutate(subject_name_plot = paste0(donor_name_plot, row_number())) %>% 
  gather("Sample_type", "Subject_ID", 3:5) %>% 
  mutate(plot_names = ifelse(Sample_type=="Donor_Sample_ID", donor_name_plot, subject_name_plot)) %>% 
  ungroup() %>% 
  dplyr::select(-c(Pt_ID, Triad, Study_origin, subject_name_plot,donor_name_plot)) %>% 
  distinct()

appetite_met <- appetite %>% 
  mutate(donor_name_plot=LETTERS[match(.$Donor_Sample_ID, unique(.$Donor_Sample_ID))]) %>% 
  group_by(donor_name_plot) %>%
  mutate(subject_name_plot = paste0(donor_name_plot, row_number())) %>% 
  gather("Sample_type", "Subject_ID", 3:5) %>% 
  mutate(plot_names = ifelse(Sample_type=="Donor_Sample_ID", donor_name_plot, subject_name_plot)) %>% 
  ungroup() %>% 
  dplyr::select(-c(Pt_ID, Triad, Study_origin, subject_name_plot,donor_name_plot)) %>% 
  distinct()

febaligo_met <- febaligo %>% 
  mutate(donor_name_plot=LETTERS[match(.$Donor_Sample_ID, unique(.$Donor_Sample_ID))]) %>% 
  group_by(donor_name_plot) %>%
  mutate(subject_name_plot = paste0(donor_name_plot, row_number())) %>% 
  gather("Sample_type", "Subject_ID", 3:5) %>%
  mutate(plot_names = ifelse(Sample_type=="Donor_Sample_ID", donor_name_plot, subject_name_plot)) %>% 
  ungroup() %>% 
  dplyr::select(-c(Pt_ID, Triad, Study_origin, subject_name_plot,donor_name_plot)) %>% 
  distinct()

#profiled strains post_fmt
strains_post <- st %>% 
  dplyr::filter(grepl("MPP", .$X1)) %>% 
  dplyr::select(-X2) %>% 
  dplyr::mutate(across(everything(), ~replace(., . %in% c("N.A.", "NA", "N/A", ""), NA)))

strains_pre <- st %>% 
  dplyr::filter(grepl("MPS_M", .$X1)) %>% 
  dplyr::select(-X2) %>% 
  dplyr::mutate(across(everything(), ~replace(., . %in% c("N.A.", "NA", "N/A", ""), NA)))

strains_donor <- st %>% 
  dplyr::filter(grepl("MPS_H", .$X1)) %>% 
  dplyr::select(-X2) %>% 
  dplyr::mutate(across(everything(), ~replace(., . %in% c("N.A.", "NA", "N/A", ""), NA)))

#for-loop to identify the total SGBs that were used in strainphlan for each respective post-FMT sample
triad$total_strains_postfmt <- NA
triad$total_strains_prefmt <- NA
triad$total_strains_donor <- NA
for (i in 1:nrow(triad)){
  print(i)
  Post_FMT_ID = triad$Post_FMT[i]
  Pre_FMT_ID = triad$Pre_FMT[i]
  Donor_FMT_ID = triad$Donor_Sample_ID[i]
  n_sgb_post <- strains_post %>% 
    filter(X1 %in% Post_FMT_ID) %>%
    dplyr::select_if(~sum(!is.na(.)) > 0) %>% 
    dplyr::select(-1) %>% 
    ncol(.)
  n_sgb_pre <- strains_pre %>% 
    filter(X1 %in% Pre_FMT_ID) %>%
    dplyr::select_if(~sum(!is.na(.)) > 0) %>% 
    dplyr::select(-1) %>% 
    ncol(.)
  n_sgb_donor <- strains_donor %>% 
    filter(X1 %in% Donor_FMT_ID) %>%
    dplyr::select_if(~sum(!is.na(.)) > 0) %>% 
    dplyr::select(-1) %>% 
    ncol(.)
  triad[triad$Post_FMT==Post_FMT_ID,]$total_strains_postfmt <- n_sgb_post
  triad[triad$Pre_FMT==Pre_FMT_ID,]$total_strains_prefmt <- n_sgb_pre
  triad[i,]$total_strains_donor <- n_sgb_donor
}

# definition strain-sharing rate (Here called shared_strain_percentage) (see natmed article; DOI: 10.1038/s41591-022-01964-3):
# -> total number of shared strains between two samples divided by the number of species profiled by strainphlan in common between the two samples.
st_feb <- st %>%
  dplyr::filter(X1 %in% febaligo_met$Subject_ID & X2 %in% febaligo_met$Subject_ID) %>%
  dplyr::select(1,2, where(is.logical)) %>% 
  rowwise() %>% 
  mutate(true_counts = sum(c_across(3:ncol(.)),na.rm = T),
         false_counts = sum(!c_across(3:ncol(.)),na.rm = T)) %>% 
  ungroup() %>% 
  mutate(shared_strain_percentage = round((true_counts / (true_counts+false_counts))*100),0) %>% 
  dplyr::select(X1,X2, shared_strain_percentage) %>% 
  tidyr::pivot_wider(., names_from = X1, values_from = shared_strain_percentage) %>% 
  column_to_rownames("X2") %>% 
  dplyr::select(unique(febaligo_met$Subject_ID)) %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  dplyr::select(unique(febaligo_met$Subject_ID)) 

st_fat <- st %>%
  dplyr::filter(X1 %in% fatmed_met$Subject_ID & X2 %in% fatmed_met$Subject_ID) %>%
  dplyr::select(1,2, where(is.logical)) %>% 
  rowwise() %>% 
  mutate(true_counts = sum(c_across(3:ncol(.)),na.rm = T),
         false_counts = sum(!c_across(3:ncol(.)),na.rm = T)) %>% 
  ungroup() %>% 
  mutate(shared_strain_percentage = round((true_counts / (true_counts+false_counts))*100),0) %>% 
  dplyr::select(X1,X2, shared_strain_percentage) %>% 
  tidyr::pivot_wider(., names_from = X1, values_from = shared_strain_percentage) %>% 
  column_to_rownames("X2") %>% 
  dplyr::select(unique(fatmed_met$Subject_ID)) %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  dplyr::select(unique(fatmed_met$Subject_ID))

st_app <-   st %>%
  dplyr::filter(X1 %in% appetite_met$Subject_ID & X2 %in% appetite_met$Subject_ID) %>%
  dplyr::select(1,2, where(is.logical)) %>% 
  rowwise() %>% 
  mutate(true_counts = sum(c_across(3:ncol(.)),na.rm = T),
         false_counts = sum(!c_across(3:ncol(.)),na.rm = T)) %>% 
  ungroup() %>% 
  mutate(shared_strain_percentage = round((true_counts / (true_counts+false_counts))*100),0) %>% 
  dplyr::select(X1,X2, shared_strain_percentage) %>% 
  tidyr::pivot_wider(., names_from = X1, values_from = shared_strain_percentage) %>% 
  column_to_rownames("X2") %>% 
  dplyr::select(unique(appetite_met$Subject_ID)) %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  dplyr::select(unique(appetite_met$Subject_ID))


###Febaligo study figure
g_feb <- graph.adjacency(as.matrix(st_feb),
                     mode="undirected",
                     weighted=T,
                     diag = F)
g_feb <- simplify(g_feb, remove.multiple=TRUE, remove.loops=TRUE)
E(g_feb)$color <- paste0("gray",100-E(g_feb)$weight)
E(g_feb)$weight <- abs(E(g_feb)$weight)

# Remove any vertices remaining that have no edges
g_feb <- delete.vertices(g_feb, degree(g_feb)==0)

# Change shape of graph vertices
V(g_feb)$shape <- "sphere"

# Change colour of graph vertices
V(g_feb)$color <- ifelse(grepl("MPS_M", V(g_feb)$name), "darkorange", ifelse(grepl("MPP", V(g_feb)$name), "#00b200", "#B1D4E0"))

# Change colour of vertex frames
V(g_feb)$vertex.frame.color <- "white"

# Assign names to the graph vertices (optional)
V(g_feb)$name <- febaligo_met$plot_names[match(V(g_feb)$name, febaligo_met$Subject_ID)]

# Put the sharing lines on top (set lines to NA that have a sharing rate below 15)
g1_feb <- g_feb
E(g1_feb)[E(g1_feb)$weight <=15]$color <- NA

#### FATMED study figure
g_fat <- graph.adjacency(as.matrix(st_fat),
                     mode="undirected",
                     weighted=T,
                     diag = F)
g_fat <- simplify(g_fat, remove.multiple=TRUE, remove.loops=TRUE)
E(g_fat)$color <- paste0("gray",100-E(g_fat)$weight)
E(g_fat)$weight <- abs(E(g_fat)$weight)

# Remove any vertices remaining that have no edges
g_fat <- delete.vertices(g_fat, degree(g_fat)==0)

# Change shape of graph vertices
V(g_fat)$shape <- "sphere"

# Change colour of graph vertices
V(g_fat)$color <- ifelse(grepl("MPS_M", V(g_fat)$name), "darkorange", ifelse(grepl("MPP", V(g_fat)$name), "#00b200", "#B1D4E0"))

# Change colour of vertex frames
V(g_fat)$vertex.frame.color <- "white"

# Assign names to the graph vertices (optional)
V(g_fat)$name <- fatmed_met$plot_names[match(V(g_fat)$name, fatmed_met$Subject_ID)]

# Put the sharing lines on top (set lines to NA that have a sharing rate below 15)
g1_fat <- g_fat
E(g1_fat)[E(g1_fat)$weight <=15]$color <- NA

#### APPETITE study figure
g_app <- graph.adjacency(as.matrix(st_app),
                     mode="undirected",
                     weighted=T,
                     diag = F)
g_app <- simplify(g_app, remove.multiple=TRUE, remove.loops=TRUE)
E(g_app)$color <- paste0("gray",100-E(g_app)$weight)
E(g_app)$weight <- abs(E(g_app)$weight)

# Remove any vertices remaining that have no edges
g_app <- delete.vertices(g_app, degree(g_app)==0)

# Change shape of graph vertices
V(g_app)$shape <- "sphere"

# Change colour of graph vertices
V(g_app)$color <- ifelse(grepl("MPS_M", V(g_app)$name), "darkorange", ifelse(grepl("MPP", V(g_app)$name), "#00b200", "#B1D4E0"))

# Change colour of vertex frames
V(g_app)$vertex.frame.color <- "white"

# Assign names to the graph vertices (optional)
V(g_app)$name <- appetite_met$plot_names[match(V(g_app)$name, appetite_met$Subject_ID)]

#Group for legend
legend_stuff <- data.frame(
  group = factor(c("Pre FMT", "Post FMT", "Donor")),
  color = factor(unique(V(g_app)$color), levels = c("#00b200", "darkorange", "#B1D4E0")))

# Put the sharing lines on top (set lines to NA that have a sharing rate below 15)
g1_app <- g_app
E(g1_app)[E(g1_app)$weight <=15]$color <- NA

#Saving the figure
cairo_pdf(filename ="Manuscript/Main_Figures/Figure_1.pdf", width = 34, height = 8)

# Set up the layout with 1 row and 3 columns
par(mfrow = c(1, 4), oma = c(0, 0, 0, 0), mar = c(4, 4, 2, 1), cex.main = 2.0)

# Plot 1: Appetite study
set.seed(1)
plot(g_app,
     main = "",
     vertex.label.color = "black",
     vertex.label.dist = 1.2,
     vertex.label.cex = 2,
     vertex.size = 7,
     layout = layout.fruchterman.reingold,
     edge.width = 3)
set.seed(1)
plot(g1_app,
     vertex.label.color = "black",
     vertex.label.dist = 1.2,
     vertex.size = 7,
     vertex.label.cex = 2,
     layout = layout.fruchterman.reingold,
     add = TRUE,
     edge.width = 3)

# Add a custom main title with reduced line height
title(main = "Appetite study")

# Plot 2: Fatmed study
set.seed(1)
plot(g_fat,
     main = "",
     vertex.label.color = "black",
     vertex.label.dist = 1.2,
     vertex.size = 7,
     vertex.label.cex = 2,
     layout = layout.fruchterman.reingold,
     edge.width = 3)
set.seed(1)
plot(g1_fat,
     vertex.label.color = "black",
     vertex.label.dist = 1.2,
     vertex.size = 7,
     vertex.label.cex = 2,
     layout = layout.fruchterman.reingold,
     add = TRUE,
     edge.width = 3)
title(main = "Fatmed study")

# Plot 3: Febaligo study
set.seed(1)
plot(g_feb,
     main = "",
     vertex.label.color = "black",
     vertex.label.dist = 1.2,
     vertex.size = 7,
     vertex.label.cex = 2,
     layout = layout.fruchterman.reingold,
     edge.width = 3)
set.seed(1)
plot(g1_feb,
     vertex.label.color = "black",
     vertex.label.dist = 1.2,
     vertex.size = 7,
     vertex.label.cex = 2,
     layout = layout.fruchterman.reingold,
     add = TRUE,
     edge.width = 3)
title(main = "Febaligo study")

# Plot 4: Legend
plot(NULL)
legend(x=.1, y=.9,
       legend = c("Pre FMT", "Post FMT", "Donor"),
       fill = c("darkorange", "#00b200", "#B1D4E0"),
       border = TRUE,
       title = "Sample type",
       cex = 2.5) 
legend(NULL, x=0.1,y=.6 ,lty=1, col=c("gray95","gray50","gray1"),lwd=2,
       legend=c("0 %", "50 %", "100 %"),
       border=NA,
       title="Shared strains",
       cex = 2.5)
dev.off()

##dataframe for x/y axis fraction sharing with donor and fraction sharing with pre-fmt sample
#x-axis -> fraction of strains post-fmt strains shared with pre-FMT | y-axis -> fraction of post-fmt strains shared with donor
#Formal definition (see natmed article: DOI: 10.1038/s41591-022-01964-3):
# -> The fraction of retained strains is defined as: the fraction of post-FMT strains shared with pre-FMT (shared strains between
# post-FMT and pre-FMT divided by the number of strains profiled at post-FMT)
# -> The fraction of donor strains as the fraction of post-FMT strains shared with
# the donor (shared strains between post-FMT and donor divided by the number of strains profiled at post-FMT).

triad_ff <- triad %>% 
  dplyr::select(-c(total_strains_prefmt, total_strains_donor)) %>% 
  dplyr::select(3,4, everything()) %>% 
  gather("don_or_pre_type", value="don_or_pre_sample", 1:2)

st_fractionfig <- st %>%
  dplyr::mutate(X_comb = paste0 (X1, X2)) %>% 
  dplyr::filter(X_comb %in% paste0(triad_ff$don_or_pre_sample, triad_ff$Post_FMT)) %>% 
  dplyr::select(1,2, where(is.logical)) %>%
  dplyr::select_if(~sum(!is.na(.)) > 0) %>% 
  rowwise() %>% 
  dplyr::mutate(true_counts = sum(c_across(3:ncol(.)),na.rm = T),
                false_counts = sum(!c_across(3:ncol(.)),na.rm = T)) %>% 
  ungroup() %>% 
  dplyr::select(X1,X2, true_counts) %>%
  dplyr::mutate(don_or_pre_type = ifelse(grepl("MPS_M",X1), "Pre_FMT_sharing_rate", "Donor_sharing_rate"))

df_fractionfig <- triad %>% 
  dplyr::select(Study_origin, Donor_Sample_ID, Pre_FMT, Post_FMT, total_strains_postfmt) %>% 
  dplyr::mutate(donor_comb = paste0(Donor_Sample_ID, Post_FMT),
                pre_comb = paste0(Pre_FMT, Post_FMT)) %>% 
  dplyr::mutate(Pre_FMT_sharing_counts = st_fractionfig$true_counts[match(.$pre_comb, paste0(st_fractionfig$X1, st_fractionfig$X2))],
                Donor_sharing_counts = st_fractionfig$true_counts[match(.$donor_comb, paste0(st_fractionfig$X1, st_fractionfig$X2))]) %>% 
  dplyr::mutate(Pre_FMT_sharing_rate = Pre_FMT_sharing_counts / total_strains_postfmt,
                Donor_sharing_rate = Donor_sharing_counts / total_strains_postfmt) %>% 
  dplyr::select(-c(donor_comb, pre_comb))

# Assuming df_centroids is your centroid data frame
df_centroids <- df_fractionfig %>%
  group_by(Study_origin) %>%
  summarize(
    centroid_x = mean(Pre_FMT_sharing_rate),
    centroid_y = mean(Donor_sharing_rate))

p2 <- ggplot(data = df_fractionfig, aes(x = Pre_FMT_sharing_rate, y = Donor_sharing_rate)) +
  geom_point(aes(color = Study_origin), size = 2) +
  geom_point(data = df_centroids, aes(x = centroid_x, y = centroid_y, color = Study_origin), size = 4, shape = 17, show.legend = F) +
  theme_Publication() +
  scale_x_continuous(limits = c(0, 0.75), breaks = c(0.0, 0.25, 0.5, 0.75), labels = c(0, 0.25, 0.5, 0.75)) +
  scale_y_continuous(limits = c(0, 0.75), breaks = c(0.0, 0.25, 0.5, 0.75), labels = c(0, 0.25, 0.5, 0.75)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(color = 'Study') +
  xlab("Fraction of post-FMT strains shared with pre-FMT") +
  ylab("Fraction of Post-FMT strains shared with donor") +
  scale_color_manual(
    name = "Study",
    labels = c("Appetite", "Fatmed", "Febaligo"),
    values = c("#F8766D", "#00BA38", "#619CFF")
  ) +
  guides(size = "none")
ggsave(filename = "Figure_2.pdf",plot = p2,
path = "Manuscript/Main_Figures/", height = 7, width = 8)

#### Permanova on strain-sharing-based dissimilarities
# strain-sharing-based dissimilarity -> within each dataset
# Ianiro (Segata) 2022 -> dissimilarity = 1-(n/M);
# n = number of shared strains
# M = maximum number of shared strains -> assuming this means the SGBs assessed by strainphlan in which both subjects occur
st_dissimilarity <- st %>%
  dplyr::mutate(X_comb = paste0(X1, X2)) %>% 
  dplyr::select(1,2, where(is.logical)) %>%
  rowwise() %>% 
  dplyr::mutate(true_counts = sum(c_across(3:ncol(.)),na.rm = T),
                false_counts = sum(!c_across(3:ncol(.)),na.rm = T)) %>% 
  ungroup() %>% 
  dplyr::select(X1,X2, true_counts, false_counts) %>%
  dplyr::mutate(max_shared_strains = true_counts + false_counts) %>% 
  dplyr::mutate(Dissimilarity = 1 - (true_counts / max_shared_strains)) %>% 
  dplyr::select(X1, X2, Dissimilarity)

st_dis_feb <- st_dissimilarity %>%
  dplyr::filter(X1 %in% febaligo_met$Subject_ID & X2 %in% febaligo_met$Subject_ID) %>% 
  tidyr::pivot_wider(., names_from = X1, values_from = Dissimilarity) %>% 
  column_to_rownames("X2") %>% 
  mutate_all(~ replace_na(., 0)) %>%
  dplyr::select(unique(febaligo_met$Subject_ID)) %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  dplyr::select(unique(febaligo_met$Subject_ID))

st_dis_fat <- st_dissimilarity %>%
  dplyr::filter(X1 %in% fatmed_met$Subject_ID & X2 %in% fatmed_met$Subject_ID) %>% 
  tidyr::pivot_wider(., names_from = X1, values_from = Dissimilarity) %>% 
  column_to_rownames("X2") %>% 
  mutate_all(~ replace_na(., 0)) %>%
  dplyr::select(unique(fatmed_met$Subject_ID)) %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  dplyr::select(unique(fatmed_met$Subject_ID))

st_dis_appetite <- st_dissimilarity %>%
  dplyr::filter(X1 %in% appetite_met$Subject_ID & X2 %in% appetite_met$Subject_ID) %>% 
  tidyr::pivot_wider(., names_from = X1, values_from = Dissimilarity) %>% 
  column_to_rownames("X2") %>% 
  mutate_all(~ replace_na(., 0)) %>%
  dplyr::select(unique(appetite_met$Subject_ID)) %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  dplyr::select(unique(appetite_met$Subject_ID))

## Febaligo estimates & p-values
febaligo_met <- febaligo_met %>% 
  dplyr::mutate(donor_relation = gsub("[0-9]+","", febaligo_met$plot_names),
                triad = ifelse(grepl("MPS_H", febaligo_met$Subject_ID), "", 
                               ifelse(grepl("MPS_M", febaligo_met$Subject_ID), 
                                      triad$Triad[match(febaligo_met$Subject_ID, triad$Pre_FMT)],
                                      triad$Triad[match(febaligo_met$Subject_ID, triad$Post_FMT)])))

## LM
febaligo_pp <- febaligo %>% gather("Time", "Subject", 4:5)
donors <- colnames(as.matrix(st_dis_feb)) %in% febaligo_met$Subject_ID[1:2]
subjects <- !colnames(as.matrix(st_dis_feb)) %in% febaligo_met$Subject_ID[1:2]
df.long <- st_dis_feb %>% 
  mutate(rownames = rownames(.))
df.long <- reshape2::melt(df.long[subjects, donors]) %>%
  setNames(c("Subject","Donor", "Dissimilarity")) %>% 
  mutate(TF_donor = ifelse(paste0(.$Subject, .$Donor) %in% paste0(febaligo_pp$Subject, febaligo_pp$Donor_Sample_ID), T, F), 
         Sample_type = ifelse(grepl("MPS_M", .$Subject), "Pre-FMT", "Post-FMT")) %>% 
  mutate(subject_rep = gsub("MPS_","",gsub("MPP_","", Subject)))
summary(lm(Dissimilarity ~ TF_donor+Sample_type, data=df.long))

## APPETITE estimates & p-values
appetite_met <- appetite_met %>% 
  dplyr::mutate(donor_relation = gsub("[0-9]+","", appetite_met$plot_names),
                triad = ifelse(grepl("MPS_H", appetite_met$Subject_ID), "", 
                               ifelse(grepl("MPS_M", appetite_met$Subject_ID), 
                                      triad$Triad[match(appetite_met$Subject_ID, triad$Pre_FMT)],
                                      triad$Triad[match(appetite_met$Subject_ID, triad$Post_FMT)])))

## LM
appetite_pp <- appetite %>% gather("Time", "Subject", 4:5)
donors <- colnames(as.matrix(st_dis_appetite)) %in% appetite_met$Subject_ID[1:2]
subjects <- !colnames(as.matrix(st_dis_appetite)) %in% appetite_met$Subject_ID[1:2]
df.long <- st_dis_appetite %>% 
  mutate(rownames = rownames(.))
df.long <- reshape2::melt(df.long[subjects, donors]) %>%
  setNames(c("Subject","Donor", "Dissimilarity")) %>% 
  mutate(TF_donor = ifelse(paste0(.$Subject, .$Donor) %in% paste0(appetite_pp$Subject, appetite_pp$Donor_Sample_ID), T, F),
         Sample_type = ifelse(grepl("MPS_M", .$Subject), "Pre-FMT", "Post-FMT")) %>% 
  mutate(subject_rep = gsub("MPS_","",gsub("MPP_","", Subject)))
summary(lm(Dissimilarity ~ TF_donor+Sample_type, data=df.long)) 

## FATMED estimates & p-values
fatmed_met <- fatmed_met %>% 
  dplyr::mutate(donor_relation = gsub("[0-9]+","", fatmed_met$plot_names),
                triad = ifelse(grepl("MPS_H", fatmed_met$Subject_ID), "", 
                               ifelse(grepl("MPS_M", fatmed_met$Subject_ID), 
                                      triad$Triad[match(fatmed_met$Subject_ID, triad$Pre_FMT)],
                                      triad$Triad[match(fatmed_met$Subject_ID, triad$Post_FMT)])))

## LM
fatmed_pp <- fatmed %>% gather("Time", "Subject", 4:5)
donors <- colnames(as.matrix(st_dis_fat)) %in% fatmed_met$Subject_ID[1:2]
subjects <- !colnames(as.matrix(st_dis_fat)) %in% fatmed_met$Subject_ID[1:2]
df.long <- st_dis_fat %>% 
  mutate(rownames = rownames(.))
df.long <- reshape2::melt(df.long[subjects, donors]) %>%
  setNames(c("Subject","Donor", "Dissimilarity")) %>% 
  mutate(TF_donor = ifelse(paste0(.$Subject, .$Donor) %in% paste0(fatmed_pp$Subject, fatmed_pp$Donor_Sample_ID), T, F),
         Sample_type = ifelse(grepl("MPS_M", .$Subject), "Pre-FMT", "Post-FMT")) %>% 
  mutate(subject_rep = gsub("MPS_","",gsub("MPP_","", Subject)))
summary(lm(Dissimilarity ~ TF_donor+Sample_type, data=df.long))



#### Supplementary table 2
met_delta <- (met_post %>% dplyr::select(-Study_origin) - met_pre %>% dplyr::select(-Study_origin)) %>% 
  .[order(as.numeric(gsub("M","", rownames(.)))),] %>% 
  dplyr::mutate(Donor_sharing_rate = df_fractionfig$Donor_sharing_rate[match(rownames(.), gsub("MPP_","", df_fractionfig$Post_FMT))])


#create empty df
df <- data.frame(
  `Clinical parameter` = character(),
  Estimate = numeric(),
  p_value = numeric())

#loop
for (i in 1:(ncol(met_delta)-1)){
  name_clinparam <- names(met_delta)[i]
  df1 <- met_delta %>% dplyr::select(all_of(name_clinparam), Donor_sharing_rate) %>% 
    setNames(c("clinparam", "Donor_sharing_rate"))
  cor_res <- cor.test(df1$clinparam, df1$Donor_sharing_rate, method = "spearman")
  df[i,] <- c(name_clinparam, round(cor_res$estimate,2), round(cor_res$p.value,3))
}

#of interest -> Diast & RD
df <- df[order(as.numeric(df$p_value)),]
rio::export(df,"Manuscript/Supplementary_information/Supplementary_table_S2.xlsx")

