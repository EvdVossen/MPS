# Script to find the most engrafting SGBs with their names.
rm(list=ls())

#Load packages and functions
source("~/Data_files/MPS/Eduard/Data_transfer/R-scripts/functions.R")

#Set working directory
setwd(path_data)

#import files (imports are done in the functions script)
d <- get_data(triad = T, triad_long = T, mp4 = T, st=T, met=T, metab = T)
base::list2env(d,envir=.GlobalEnv)

mp4_post <- mp4 %>% filter(rownames(.) %in% triad$Post_FMT)

mp4_stats <- apply(mp4,2,summary) %>% 
  t() %>%  
  as.data.frame()
mp4_stats_post_FMT <- mp4_post %>% 
  apply(.,2,summary) %>% 
  t() %>% 
  as.data.frame()

### total strain engrafting rates (first without looking at substudies)
stt <- st[grepl("MPS_H", st$X1) & grepl("MPP_M", st$X2),] %>% 
  dplyr::filter(paste0(.$X1, .$X2) %in% paste0(triad$Donor_Sample_ID, triad$Post_FMT)) %>% 
  dplyr::select(where(is.logical)) %>% 
  colSums(., na.rm = T) %>% 
  as.data.frame(.) %>% 
  setNames("n_true") %>% 
  filter(n_true>=5) %>% 
  mutate(SGB_ID = rownames(.),
         median_abundance_all = mp4_stats$Median[match(rownames(.), rownames(mp4_stats))],
         mean_abundance_all =  mp4_stats$Mean[match(rownames(.), rownames(mp4_stats))],
         median_abundance_post_FMT = mp4_stats_post_FMT$Median[match(rownames(.), rownames(mp4_stats_post_FMT))],
         mean_abundance_post_FMT =  mp4_stats_post_FMT$Mean[match(rownames(.), rownames(mp4_stats_post_FMT))]) %>%
  mutate(plot_name = mp4_db$plot_name[match(rownames(.), mp4_db$SGB_ID)]) %>% 
  .[order(as.numeric(-.$n_true)),]

top_n_feat = 30
  
mp4_post_top <- mp4_post[,stt$SGB_ID[1:top_n_feat]] %>%
  gather(., "SGB_ID", "relative_abundance", 1:top_n_feat) %>%
  mutate(Phylum=mp4_db$phylum[match(.$SGB_ID, mp4_db$SGB_ID)])

top_feat_30 <- stt[stt$SGB_ID %in% mp4_post_top$SGB_ID,] %>%
  arrange(factor(SGB_ID, levels = unique(mp4_post_top))) %>%
  `rownames<-`(.$SGB_ID)
  
my_palette <- colorRampPalette(c("black", "darkgrey"))(n = top_n_feat)
  
p2 <-
ggplot(top_feat_30, aes(x = reorder(plot_name, as.double((n_true/nrow(mp4_post))*100)), y = as.double((n_true/nrow(mp4_post))*100))) +
  geom_segment(
    aes(x = reorder(plot_name, as.double((n_true/nrow(mp4_post))*100)), xend = plot_name, y = 0.0, yend = as.double((n_true/nrow(mp4_post))*100)),
    color = my_palette, size = 1) +
  geom_point(shape = 21, colour = my_palette, fill = "white", size = 2, stroke = 2) +
  theme_Publication()+
  xlab(" ") + ylab("Strain transfer (%)") +
  scale_y_continuous(limits = c(0,50)) +
  coord_flip()

#### addition heatmap
st_hm <- st[grepl("MPS_H", st$X1) & grepl("MPP_M", st$X2),] %>% 
  dplyr::filter(paste0(.$X1, .$X2) %in% paste0(triad$Donor_Sample_ID, triad$Post_FMT)) %>% 
  dplyr::mutate(Subject_ID = gsub("MPP_","", X2)) %>% 
  tibble::column_to_rownames("Subject_ID") %>% 
  dplyr::select(all_of(rownames(top_feat_30))) %>%
  dplyr::mutate(across(everything(), ~ replace(., is.na(.), FALSE))) %>% #comment this partif also white should be in (similar to figure 1 "analyses")
  setNames(top_feat_30$plot_name[match(names(.), top_feat_30$SGB_ID)]) %>% 
  dplyr::mutate(`Study origin` = factor(str_to_title(triad$Study_origin)[match(rownames(.), gsub("MP-","", triad$Pt_ID))], levels = c("Appetite","Fatmed","Febaligo")),
                Subject_ID = rownames(.)) %>% 
  dplyr::arrange(`Study origin`, as.numeric(gsub("M","", rownames(.)))) %>%
  dplyr::mutate(Subject_ID = factor(Subject_ID, levels = unique(.$Subject_ID))) %>%
  ungroup()

st_hm_long <- st_hm %>% 
  tidyr::pivot_longer(data = ., cols = 1:30, values_to = "FMT Engraftment") %>% 
  dplyr::mutate(name = factor(name, levels = c(rev(top_feat_30$plot_name))))

#"firebrick1", "dodgerblue1"
p_hm <-
ggplot(st_hm_long, aes(x = Subject_ID, y = name, fill = `FMT Engraftment`)) +
  geom_tile(alpha = 0.6) +
  labs(x = "", y = "", fill = "FMT Engraftment") +
  theme_minimal() +
        # panel.grid = element_blank()) +
        # axis.text.x.top = element_text(angle = 45, hjust = 0)) +  # Adjust top x-axis labels
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = c("firebrick1", "dodgerblue1"), labels = c("No engraftment", "Engraftment")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
p_cb <-
  ggplot(st_hm)+
  geom_bar(mapping = aes(x = Subject_ID, y = 1, fill = `Study origin`), 
           stat = "identity", 
           width = 1)+
  theme_void()+
  labs(x = "") +
  scale_fill_manual(values = c("Appetite" = "#F8766D", "Fatmed" = "#00BA38", "Febaligo" = "#619CFF"),
                    labels = c("Appetite" = "Appetite (n = 12)", "Fatmed" = "Fatmed (n = 12)", "Febaligo" = "Febaligo (n = 5)"))+
  theme(panel.spacing.x = unit(1, "mm"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
p_blank <- 
ggplot(st_hm) +
  geom_blank(mapping = aes(x = Subject_ID, y = 1)) +
  theme_void() +
  theme(panel.spacing.x = unit(1, "mm")) +
  labs(x = "")
####



mp4_post_top <- mp4_post_top %>% 
  mutate(plot_name = top_feat_30$plot_name[match(.$SGB_ID, top_feat_30$SGB_ID)]) %>% 
  mutate(plot_name = factor(plot_name, levels = c(levels(reorder(top_feat_30$plot_name, as.double((top_feat_30$n_true/nrow(mp4_post)) * 100))))))

### post abundance only engrafters
st1 <- st %>% dplyr::filter(paste0(X1,X2) %in% paste0(triad$Post_FMT, triad$Donor_Sample_ID)) %>%
  dplyr::select(X1, X2, c(names(.)[names(.) %in% stt$SGB_ID[1:top_n_feat]])) %>% 
  tibble::column_to_rownames("X1") %>% dplyr::select(-X2) %>% 
  .[order(as.numeric(gsub("MPP_M","",rownames(.)))),]

mp4_post_top30 <- mp4_post %>% dplyr::select(names(st1)) %>% 
  .[order(as.numeric(gsub("MPP_M","",rownames(.)))),]
mp4_post_top30[(is.na(st1) | st1 == FALSE)] <- NA

mp4_post_top30_engrafters<- mp4_post_top30[,stt$SGB_ID[1:top_n_feat]] %>%
  gather(., "SGB_ID", "relative_abundance", 1:top_n_feat) %>% 
  mutate(Phylum=mp4_db$phylum[match(.$SGB_ID, mp4_db$SGB_ID)])%>% 
  mutate(plot_name = top_feat_30$plot_name[match(.$SGB_ID, top_feat_30$SGB_ID)]) %>% 
  mutate(plot_name = factor(plot_name, levels = c(levels(reorder(top_feat_30$plot_name, as.double((top_feat_30$n_true/nrow(mp4_post)) * 100)))))) %>% 
  dplyr::filter(!is.na(relative_abundance))
 
p3 <-
ggplot(mp4_post_top30_engrafters, aes(y=relative_abundance, x=plot_name)) +
  geom_boxplot(aes(fill=Phylum),outlier.shape = NA) +
  geom_point(position = "jitter", size = 0.1, na.rm = T)+
  theme_Publication() +
  scale_y_continuous(limits = c(0, 20),breaks = seq(0, 20, by = 2)) +
  xlab(" ") + ylab("Relative abundance (%)") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines")) +
  labs(fill = "Phylum") +
  coord_flip(ylim = c(0, 10))

## Final part with plot heatmap
legend <- plot_grid(get_legend(p_cb), get_legend(p_hm), get_legend(p3), 
                    get_legend(p_blank),get_legend(p_blank),get_legend(p_blank), #Quick & not so nice method to condense the legends 
                    ncol = 1)
p_hm <- p_hm + theme(legend.position = "none")
p_cb <- p_cb + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
plot <- plot_grid(p_cb, p_hm, p_blank, align = "v", ncol = 1, axis = "tb", rel_heights = c(0.5, 15, 1.))
plot2 <- plot_grid(p_blank, p3, align = "v", ncol = 1, axis = "tb", rel_heights = c(0.3, 15))

pf4 <- plot + plot2 + legend + plot_layout(widths = c(4, 3,1))
ggsave(filename = "Figure_4.pdf",plot = pf4, device = "pdf",
       path = "Manuscript/Main_Figures", height = 10, width = 18)
##

top_list_cca <- top_feat_30 %>% filter(median_abundance_post_FMT > 1 | mean_abundance_post_FMT > 1 | n_true >= 10) %>%
  `rownames<-`(.$SGB_ID) %>% 
  rio::export(., paste0(path_data, "/Intermediate_files/top_list_strain_engraftment.xlsx"), rowNames = T)
