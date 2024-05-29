rm(list=ls())

#Load packages and functions
source("~/Data_files/MPS/Eduard/Data_transfer/R-scripts/functions.R")

#Set working directory
setwd(path_data)

#import files (imports are done in the functions script)
d <- get_data(st = F, triad = T, triad_long = T, met=T, mp4 = T, metab = F)
base::list2env(d,envir=.GlobalEnv)

#Bray-curtis PCOA -> based on SGB level
bray <- vegan::vegdist(mp4, method = 'bray')

mds.stuff <- cmdscale(bray, eig = TRUE, x.ret = TRUE) #eig = true -> returns eigenvalues,
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100,1)
mds.values <- mds.stuff$points

#Dataframe for the plot
mds.data <- data.frame(Sample=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2]) %>% 
  dplyr::mutate(Group = gsub("_","-", gsub("_Sample_ID","", (triad_long$Triad))[match(Sample, triad_long$Subject_ID)]),
                Triad_donor = triad$Donor_Sample_ID[match(gsub("MPS_","", gsub("MPP_","",Sample)), gsub("MP-","", triad$Pt_ID))]) %>% 
  dplyr::mutate(Group = factor(Group, levels = c("Donor","Pre-FMT", "Post-FMT")),
                Triad_donor = ifelse(is.na(Triad_donor), gsub("MPS_","", Sample), gsub("MPS_","",Triad_donor)))

#The plot
p <-
  ggplot(mds.data, aes(x=X, y=Y))+
  geom_point(aes(shape = Group, color = Group), size=2)+
  stat_ellipse(aes(color=Group), level = 0.8, size=.75, show.legend = F)+
  scale_color_manual(values = c("#B1D4E0", "darkorange", "#00b200"))+
  xlab(paste0('PCo1 [',mds.var.per[1],'%]'))+
  ylab(paste0('PCo2 [',mds.var.per[2],'%]'))+
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlim(-0.4, 0.4) +
  ylim(-0.3, 0.3) +
  labs(title = "Overall PCoA") +
  theme_Publication()

#check per donor how the pcoa looks (forloop)
donors <- unique(triad$Donor_Sample_ID)[order(as.numeric(gsub("MPS_H","", unique(triad$Donor_Sample_ID))))]
ppd <- list()
ppd[[1]] <- p
for(idx in 1:length(donors)){
  donor = gsub("MPS_","", donors[[idx]])
  mds.data_don <- mds.data[mds.data$Triad_donor == donor,] %>% 
    dplyr::mutate(Subject = ifelse(.$Group == "Donor", "Donor", gsub("MPS_","", gsub("MPP_","", Sample))))
  
  ppd[[idx + 1]] <-
    ggplot(mds.data_don, aes(x=X, y=Y))+
    geom_point(aes(shape = Group, color = Subject), size=3)+
    xlab(paste0('PCo1 [',mds.var.per[1],'%]'))+
    ylab(paste0('PCo2 [',mds.var.per[2],'%]'))+
    xlim(-0.4, 0.4) +
    ylim(-0.3, 0.3) +
    labs(title = paste0("Donor sample: ",donor))+
    guides(fill = guide_legend(override.aes=list(shape=21))) +
    theme_Publication()
}

#Combine the plots using patchwork::wrap_plots()
comb_plot <- wrap_plots(ppd)

ggsave(filename = "Supplementary_Figure_S1.pdf", plot = comb_plot, device = "pdf", path = "Manuscript/Supplementary_information/", width = 20, height = 12)

print(comb_plot)