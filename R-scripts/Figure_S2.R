rm(list=ls())

#Load packages and functions
source("~/Data_files/MPS/Eduard/Data_transfer/R-scripts/functions.R")

#Set working directory
setwd(path_data)

#import files (imports are done in the functions script)
d <- get_data(st = T, triad = T, triad_long = F, met=F, mp4 = T, metab = F, path_dat = F)
base::list2env(d,envir=.GlobalEnv);

#Imports + standard dataprocessing
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

# mp4 pre + post together
mp4_pre <- mp4_pre %>% 
  mutate(Subject_ID = rownames(.)) %>% 
  `rownames<-`(paste0("MPS_", rownames(.)))
mp4_post <- mp4_post %>% 
  mutate(Subject_ID = rownames(.)) %>% 
  `rownames<-`(paste0("MPP_", rownames(.)))

mp4 <- bind_rows(mp4_pre, mp4_post) %>% 
  mutate(Time = ifelse(grepl("MPP_", rownames(.)), "Post FMT", "Pre FMT")) %>% 
  mutate(Time = factor(Time, levels = c("Pre FMT", "Post FMT"))) 

mp4_delta <- (mp4_post %>% dplyr::select(all_of(top_list_engraftment$SGB_ID)) - mp4_pre %>% dplyr::select(all_of(top_list_engraftment$SGB_ID))) %>% 
  `rownames<-`(gsub("MPP_","", rownames(.))) %>% 
  merge(., st1, by = "row.names") %>% 
  tibble::column_to_rownames("Row.names")


p <- list()
for (idx in 1:nrow(top_list_engraftment)){
  name_sgb <- top_list_engraftment$SGB_ID[idx]
  name_engraftment_sgb <- top_list_engraftment$plot_name[idx]
  parts_name <- str_split(name_engraftment_sgb, "\\(")[[1]]
  part_name <- ifelse(grepl("clade", parts_name[1]), "P. copri", parts_name[1])
  part_sgb <- ifelse(grepl("clade", parts_name[1]),  paste0(gsub("P. copri ","", parts_name[[1]]), "(",parts_name[2]), paste0("(",parts_name[2]))
  
  
  comb_for_lmer <- mp4 %>% 
    dplyr::select(all_of(name_sgb), Time, Subject_ID) %>%
    stats::setNames(c("name_sgb", "Time", "Subject_ID")) %>% 
    mutate(engraftment = st1[,name_engraftment_sgb][match(.$Subject_ID, rownames(st1))]) %>% 
    dplyr::mutate(Strain_engraftment = ifelse(engraftment == 0, "No engraftment", "Engraftment")) %>% 
    dplyr::mutate(Strain_engraftment = factor(Strain_engraftment, levels = c("No engraftment","Engraftment"))) %>% 
    dplyr::mutate(Group = paste0(Time, " ", Strain_engraftment)) %>% 
    dplyr::mutate(Group = factor(Group, levels = c("Pre FMT No engraftment", "Post FMT No engraftment", "Pre FMT Engraftment", "Post FMT Engraftment")))

lmer_test <- summary(lmerTest::lmer(name_sgb ~ engraftment * Time + (1|Subject_ID), data = comb_for_lmer))

# for the barplot
mp4_delta_fp <- mp4_delta %>% 
  dplyr::select(all_of(name_sgb), all_of(name_engraftment_sgb)) %>% 
  setNames(c("name_sgb", "name_engraftment_sgb"))

df2 <- 
  mp4_delta_fp %>% 
  data_summary(., varname="name_sgb", 
               groupnames=c("name_engraftment_sgb")) %>% 
  mutate(time = c("Delta")) %>%  
  dplyr::rename("delta_var" = median) %>% 
  dplyr::mutate(Strain_engraftment = ifelse(name_engraftment_sgb==0, "No engraftment", "Engraftment")) %>% 
  dplyr::mutate(Strain_engraftment = factor(Strain_engraftment, levels = c("No engraftment", "Engraftment"))) %>% 
  dplyr::mutate(loc_plot = c(1.5, 3.5))


p[[idx]] <-
  ggplot() +
  geom_boxplot(data = comb_for_lmer,
               aes(x = Group, 
                   y = name_sgb, 
                   fill = Strain_engraftment),
               alpha = 0.5,
               width = 0.4,
               show.legend = F) +
  geom_point(data = comb_for_lmer,
             aes(x = Group, 
                 y = name_sgb, 
                 col = Strain_engraftment),
             size = 1.5, color = "black") +
  geom_line(data = comb_for_lmer,
            aes(x = Group, 
                y = name_sgb, 
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
           width = 0.4) +
  geom_errorbar(data = df2,
                aes(x = loc_plot,
                    y = delta_var,
                    xmin = loc_plot, 
                    xmax = loc_plot,
                    ymax = delta_var + (delta_var > 0)*sd,
                    ymin = delta_var - (delta_var < 0)*sd),
                position = position_dodge(), 
                width = 0.2)  +
  stat_pvalue_manual(data = data.frame(group1 = 1.5, group2 = 3.5, p.adj = round(as.numeric(lmer_test$coefficients[20]),4)),
                     label = "p.adj",
                     y.position = max(comb_for_lmer$name_sgb, na.rm = T)*1.03,
                     tip.length = .01) +
  ylab(label = bquote("Relative abundance -"~italic(.(part_name))~.(part_sgb))) +
  xlab(label = "Time") +
  theme_Publication() +
  scale_fill_manual(name = "FMT group", labels = c("No engraftment", "Engraftment"),
                    values=c("firebrick1", "dodgerblue1","firebrick1", "dodgerblue1")) +
  theme(legend.text.align = 0,
        strip.background =  element_rect(fill = NA, colour = NA))
} 

library(patchwork)
p1 <- p[[2]] + theme(axis.text.x = element_blank())+ xlab("")
p3 <- p[[3]] + theme(axis.text.x = element_blank())+ xlab("")
p_s2 <- p1 + p3 + p[[10]] + plot_layout(nrow = 3, guides = 'collect')

ggsave(filename = paste0("Supplementary_Figure_S2.pdf"), plot = p_s2, device = "pdf",
       path =  "Manuscript/Supplementary_information/",
       height = 18, width = 12)
