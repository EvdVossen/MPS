rm(list=ls())

source("~/Data_files/MPS/Eduard/Data_transfer/R-scripts/functions.R")

setwd(path_data)

#import files (imports are done in the functions script)
d <- get_data(triad = T, triad_long = T, met=T, metab = T)
base::list2env(d,envir=.GlobalEnv)

############ Mark Multilevel PCA - metabolites
m_delta <- m_post - m_pre

#Fix metadata
met_post_comp <- met_post %>% 
  dplyr::select(-Study_origin, -Weight)
met_pre_comp <- met_pre[rownames(met_pre) %in% rownames(met_post_comp), names(met_pre) %in% names(met_post_comp),]
met_delta <- met_post_comp - met_pre_comp
met_delta_comp_diast <- met_delta %>% 
  filter(!is.na(Diast)) %>% 
  dplyr::select_if(~ !any(is.na(.)))


m_delta_comp_diast <- m_delta %>% filter(rownames(.) %in% rownames(met_delta_comp_diast))

all(rownames(met_delta_comp_diast) == rownames(m_delta_comp_diast))
cor.test(met_delta_comp_diast$Diast, m_delta_comp_diast$`1-(1-enyl-oleoyl)-GPE (P-18:1)*`, method="spearman")

df_p <- data.frame(metabolites = numeric(0),
                   estimate = numeric(0),
                   `p-value` = numeric(0))

for (idx in 1:ncol(m_delta_comp_diast)){
  # idx=1
  ct <- suppressWarnings(cor.test(met_delta_comp_diast$Diast, m_delta_comp_diast[,idx], method="spearman"))
  df_p[idx,] <- c(names(m_delta_comp_diast)[idx], ct$estimate, ct$p.value)
}

df_sig <- df_p[df_p$p.value<=0.05,] %>% 
  .[order(as.numeric(.$p.value)),]
# rio::export(df_sig,"CCA/metab_sig_diast.xlsx", rowNames=F)

p <- list()
for (j in 1:nrow(df_sig)){
  # j=1
  var_name <- df_sig$metabolites[j]
  var_oi <- m_delta_comp_diast %>% 
    dplyr::select(all_of(var_name)) %>% 
    stats::setNames("metabolite")
  
  p[[j]] <-
    ggscatter(merge(var_oi, met_delta_comp_diast, by = "row.names"),
              x = "Diast", 
              y = "metabolite",
              # xlab = name_dataset1,
              ylab = var_name,
              add = "reg.line",
              conf.int = TRUE,
              cor.coef = TRUE,
              cor.method = "spearman",
              p.title = "P-Value: {stat.test.p.value}")
}
p_cor <- wrap_plots(p)
ggsave(filename = "Figure_S5.pdf", plot = p_cor, device = "pdf",
       path = "Manuscript/Supplementary_information/", width = 16, height = 24)