library(tidyverse)

path_data <- "~/Data_files/MPS/Eduard/Data_transfer/"

get_data <- function(st = F, triad = F, triad_long = F, met = F, mp4 = F, path_dat=F, metab = F){
  new_data <- "Pipeline_output/"
  t_long <- rio::import(paste0(path_data, "Intermediate_files/Triad_data_long.xlsx"))
  
  if(st == T){
    st <- rio::import(paste0(path_data,"Intermediate_files/table_strainsharing_checked.txt"))
  }
  if(triad == T){
    triad  <- rio::import(paste0(path_data, "Intermediate_files/Triad_numbers.xlsx"))
  }
  if(triad_long == T){
    triad_long = t_long
  }
  if(met == T){
    # metadata full
    met <- xlsx::read.xlsx(paste0(path_data, "MPS_DATABASE_14-03-19.xlsx"), sheetName = "PRE & POST") %>%
      dplyr::filter(.$Index %in% t_long$Index)
    
    #metadata pre
    met_pre <- met %>% 
      dplyr::select(Index, !grep("post", names(.))) %>% 
      dplyr::mutate(Index = gsub("MP-","", Index)) %>% 
      tibble::column_to_rownames("Index") %>% 
      .[order(as.numeric(gsub("M","",rownames(.)))),] %>% 
      dplyr::mutate(CRP = as.numeric(CRP)) %>% 
      dplyr::select(-c(ASAT, GGT, Leuko, Study_origin_ID, Study_group, Intervention, Gender, Age, Length))
    
    # metadata post
    met_post <- met %>% 
      dplyr::select(Index, Study_origin, grep("post", names(.))) %>% 
      dplyr::mutate(Index = gsub("MP-","", Index)) %>% 
      tibble::column_to_rownames("Index") %>% 
      .[order(as.numeric(gsub("M","",rownames(.)))),] %>% 
      dplyr::mutate(CRP_post = as.numeric(CRP_post)) %>% 
      dplyr::select(-c(ASAT_post, GGT_post, Leuko_post)) %>% 
      setNames(gsub("_post", "", names(.))) %>% 
      dplyr::select(names(met_pre))
    
    met_bl <- xlsx::read.xlsx(paste0(path_data, "MPS_DATABASE_14-03-19.xlsx"), sheetName = "Export_bl") %>% 
      dplyr::filter(Pt_ID %in% t_long$Index & Study_group == "Healthy")
  } else{
    met_pre <- F
    met_post <- F
    met_bl <- F
  }
  if(mp4 == T){
    mp4 <- rio::import(paste0(path_data, new_data, "Metaphlan_merged_abundance_table_SGB.txt")) %>% 
      tibble::column_to_rownames("clade_name") %>% 
      t(.) %>% as.data.frame(.) %>%
      `rownames<-`(ifelse(rownames(.)=="MPS_M05","MPS_M09",
                                    ifelse(rownames(.) =="MPS_M09","MPS_M05", rownames(.)))) %>% 
      dplyr::filter(rownames(.) %in% t_long$Subject_ID) %>% 
      .[,order(colSums(-.,na.rm=TRUE))] %>% 
      dplyr::select_if(colSums(.) != 0) %>% 
      dplyr::select(-UNCLASSIFIED)
    
    ## Pre-FMT
    mp4_pre <- mp4[grep("MPS_M", rownames(mp4)),] %>% 
      `rownames<-`(gsub("MPS_","",rownames(.))) %>% 
      .[order(as.numeric(gsub("M","",rownames(.)))),]
    
    ## Post-FMT
    mp4_post <- mp4[grep("MPP_M", rownames(mp4)),] %>% 
      `rownames<-`(gsub("MPP_","",rownames(.))) %>%
      .[order(as.numeric(gsub("M","",rownames(.)))),order(-colSums(.))]
    
    ## Database
    mp4_db <- rio::import(paste0(path_data, new_data, "Bugs_list_joined.tsv")) %>% 
      dplyr::select(1) %>% 
      dplyr::filter(grepl("t__", .$clade_name)) %>% 
      dplyr::mutate(SGB_ID = gsub(".*t__","", .$clade_name),
             species = gsub(".{1}$","", gsub(".*s__", "", gsub("t__.*","", .$clade_name))),
             phylum = gsub(".{1}$","", gsub(".*p__", "", gsub("c__.*","", .$clade_name))),
             family = gsub(".{1}$","", gsub(".*f__", "", gsub("g__.*","", .$clade_name)))) %>% 
      dplyr::mutate(names_abb = ifelse(grepl("sp", species), species, 
                                ifelse(grepl("copri", species), paste0("P. copri", gsub("Prevotella_copri_", " ", .$species)),
                                       ifelse(grepl("GGB", species), family,
                                              ifelse(grepl("GGB", species) & grepl("FGB", family), phylum, 
                                                     ifelse(grepl("KLE", species), family,
                                                            ifelse(grepl("SGB", species), gsub("_SGB.*","",species), 
                                                                   gsub("^(\\w).*_", "\\1. ",species)))))))) %>% 
      dplyr::mutate(plot_name = paste0(names_abb, " (", SGB_ID, ")")) %>% 
      dplyr::mutate(plot_name = gsub("_"," ", gsub("unclassified ","", plot_name))) %>% 
      dplyr::mutate(plot_name = ifelse(grepl("Lacrimispora", plot_name), gsub("Lacrimispora", "L.", plot_name), plot_name))
  } else{
    mp4_pre <- F
    mp4_post <- F
    mp4_db <- F
  }
  if(path_dat == T){
    path_dat <- rio::import(paste0(path_data, new_data, "Humann_merged_path_abundance_cpm_unstratified.txt")) %>% 
      column_to_rownames("# Pathway") %>% 
      t(.) %>% as.data.frame(.) %>% 
      `rownames<-`(ifelse(rownames(.)=="MPS_M05","MPS_M09",
                          ifelse(rownames(.) =="MPS_M09","MPS_M05", rownames(.)))) %>% 
      dplyr::filter(rownames(.) %in% t_long$Subject_ID) %>%
      .[,order(colSums(-., na.rm=T))] %>% 
      dplyr::select_if(colSums(.) != 0) %>% 
      dplyr::select(c(-UNINTEGRATED, -UNMAPPED))

    # Pre FMT
    p_pre <- path_dat[grep("MPS_M", rownames(path_dat)),] %>% 
      select_if(colSums(.) != 0) %>% 
      `rownames<-`(gsub("MPS_","",rownames(.))) %>% 
      .[order(as.numeric(gsub("M","",rownames(.)))),]
    
    # Post FMT  
    p_post <- path_dat[grep("MPP_M", rownames(path_dat)),] %>% 
      select_if(colSums(.) != 0) %>% 
      `rownames<-`(gsub("MPP_","",rownames(.))) %>%
      .[order(as.numeric(gsub("M","",rownames(.)))),order(-colSums(.))]
  } else{
    p_post <- F
    p_pre <- F
  }
  
  if(metab == T){
    metab <- rio::import(paste0(path_data,"Metabolon/normalized_UNGO_and_UVAM_data_all.csv")) %>%
      tibble::column_to_rownames("ID") %>% 
      dplyr::filter(rownames(.) %in% t_long$Subject_ID)
    metab_info1 <- rio::import("Metabolon/prepped/Batch_1_metabolite_info.csv") 
    metab_info2 <- rio::import("Metabolon/prepped/Batch_2_metabolite_info.csv")
    metab_info <- base::rbind(metab_info1, metab_info2) %>% 
      dplyr::select(-c(V1, MASS, `PATHWAY SORTORDER`, RI)) %>% 
      dplyr::distinct()
    
    #pre data
    m_pre <- metab %>% 
      dplyr::filter(grepl("MPS", rownames(.))) %>%
      `rownames<-`(gsub("MPS_","",rownames(.))) %>% 
      .[order(as.numeric(gsub("M","",rownames(.)))),]
    
    #post data
    m_post <- metab %>%
      dplyr::filter(grepl("MPP", rownames(.))) %>%
      `rownames<-`(gsub("MPP_","",rownames(.))) %>% 
      .[order(as.numeric(gsub("M","",rownames(.)))),]
  } else{
    metab_info <- F
    m_pre <- F
    m_post <- F
  }
  
  d <- list(st, triad, triad_long, 
            met, met_pre, met_post, met_bl,
            mp4, mp4_db, mp4_pre, mp4_post,
            path_dat, p_pre, p_post,
            metab, metab_info, m_pre, m_post)
  names(d) <- c("st", "triad", "triad_long", 
                "met", "met_pre", "met_post", "met_bl",
                "mp4", "mp4_db", "mp4_pre", "mp4_post",
                "path_dat", "p_pre", "p_post",
                "metab", "metab_info", "m_pre", "m_post")
  
  d <-  base::Filter(function(x) !is.logical(x) || x, d)
  
  suppressWarnings(dir.create(paste0(path_data, "Manuscript/Main_Figures"), recursive = T))
  suppressWarnings(dir.create(paste0(path_data, "Manuscript/Supplementary_information"), recursive = T))
  
  return(d)
}

# publication theme
theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.spacing  = unit(0, "cm"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

# Python's standardscaler
standard_scaler <- function(data, ddof=0) {
  
  #standard deviation with the delta degrees of freedom (ddof) -> similar to Sklearn's StandardScaler
  sd_fun_ddof <- function(data, ddof=ddof){
    stds <- apply(data, 2, function(x){sqrt(abs(sum((x - mean(x))^2)) / (length(x) - ddof))})
    return(stds)
  }
  
  # Calculate the mean and standard deviation for each column
  means <- colMeans(data)
  
  # Calculate the standard deviation (included the ddof)
  stds <- sd_fun_ddof(data = data, ddof = ddof)
  
  # Scale the data using the mean and standard deviation
  data_list_scaled <- vector(mode = "list", length = ncol(data))
  for (idxi in 1:ncol(data)){
    data_list_scaled[[idxi]] = (data[,idxi]-means[idxi])/stds[idxi]
  }
  scaled_data <- as.data.frame(do.call("cbind", data_list_scaled)) %>% setNames(names(data))
  
  return(scaled_data)
}

## Function CCA -> easy output per map to screen for potential interesting hits to do statistics on.
CCA_fun <- function(X, Y, 
                    x_data_name = "dataset_1", y_data_name = "dataset_2", 
                    color_x_data = "red", color_y_data = "blue",
                    n_comp = 10, validation = "loo", 
                    grid1 = seq(0.001, 1, length = 5), 
                    grid2 = seq(0.001, 1, length = 5),
                    name_map = "example_map", 
                    name_submap = "submap",
                    threshold_cor = 0.8,
                    rad.in = 0.5,
                    cor.method = "pearson",
                    output_path = getwd()){
  
  #necessary libraries
  library(mixOmics)
  library(patchwork)
  
  #create the map for the output
  file_path <- file.path(output_path, name_map, name_submap)
  suppressWarnings(dir.create(file_path, recursive = T))
  
  cat(str_glue("Saving inputs: {x_data_name} & {y_data_name} \n"))
  rio::export(X, paste0(file.path(file_path, x_data_name), ".xlsx"), format = "xlsx", rowNames=T)
  rio::export(Y, paste0(file.path(file_path, y_data_name), ".xlsx"), format = "xlsx", rowNames=T)
  
  #scale using Python's standardscaler (rewritten in R; delta degrees of freedom included)
  x1_sc <- standard_scaler(X, ddof = 0)
  x2_sc <- standard_scaler(Y, ddof = 0)
  
  # Tune lambda parameters for regularized CCA
  cat("\nTuning lambda for the regularised CCA using the prespecified grids.
If no grid is specified, than the default is taken: seq(0.001, 1, length = 5) \n")
  tune_res <- tune.rcc(X = x1_sc,
                       Y = x2_sc,
                       validation = validation,
                       grid1 = grid1,
                       grid2 = grid2,
                       plot = F)
  
  lambda.opt1 = round(tune_res$opt.lambda1,3)
  lambda.opt2 = round(tune_res$opt.lambda2,3)   
  
  #check for n_comp
  if(n_comp >= min(ncol(x1_sc), ncol(x2_sc))) {
    n_comp <- min(ncol(x1_sc), ncol(x2_sc))
    cat(str_glue("\n\n minimal columns in dataset is {min(ncol(x1_sc), ncol(x2_sc))}. Therefore, the number of components is set to \n this number.\n\n"))
  }
  
  
  # Regularized CCA
  cca_res <- mixOmics::rcc(x1_sc, 
                           x2_sc, 
                           lambda1 = tune_res$opt.lambda1, 
                           lambda2 = tune_res$opt.lambda2,
                           ncomp = n_comp)
  
  # lambda plot
  df_long <- as.data.frame(tune_res$mat) %>% `rownames<-`(tune_res$grid1) %>%
    setNames(tune_res$grid2) %>%
    as.matrix() %>% 
    reshape2::melt(varnames = c("grid1", "grid2"), value.name = "value")
  
  p1 <- ggplot(df_long, aes(grid1, grid2, fill = value)) + 
    geom_tile()+
    scale_fill_gradient2(low = "red", mid = "yellow", high = "white", midpoint = max(df_long$value)/1.065) +
    scale_x_continuous(breaks = unique(df_long$grid1), name = expression(lambda[1])) +
    scale_y_continuous(breaks = unique(df_long$grid2), name = expression(lambda[2])) +
    theme_minimal() + 
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10, vjust = 0.5, hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = -20)),
          axis.text.y = element_text(margin = margin(r = -20)),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, margin = margin(b = -25)),
          legend.title = element_blank()) +
    ggtitle(substitute(paste("CV(", lambda[1], ", ", lambda[2], ");   ",
                             lambda[1*",opt"], " = ", j, ";   ", lambda[2*",opt"], " = ", k), 
                       list(j = lambda.opt1, k = lambda.opt2)))
  
  #bar plot
  p2 <- ggplot(data.frame(index = names(cca_res$cor), cor = cca_res$cor), aes(x = index, y = cor)) +
    geom_col(fill = "grey", color = "black", linewidth = .9) +
    labs(x = "Canonical Variates",
         y = "Correlation") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ceiling(max(cca_res$cor) * 10) / 10))
  
  #combination
  p_hc <- p1 + p2 + plot_annotation(tag_levels = 'A')
  ggsave(filename = "Lambda_est_CCA_correlations.pdf", plot = p_hc, device = "pdf", path = file_path, height = 10, width = 16)
  ggsave(filename = "Supplementary_figure_S2.pdf", plot = p_hc, device = "pdf", 
         path = paste0(getwd(), "/Manuscript/Supplementary_information/"), height = 10, width = 16)
  cca_res[["p_s2"]] <- p_hc #Small line to add supplemental figure in the cca_res object
  
  
  # Add original dataframes to the CCA_RES object
  cca_res$X_orig_scale <- X
  cca_res$Y_orig_scale <- Y
  
  cat("Plotting the correlation circle plot and the most interesting correlations derived \nfrom the rCCA. \n")
  cca_res <- plotrCCA(object = cca_res, 
                      x_data_name = x_data_name, 
                      y_data_name = y_data_name, 
                      threshold_cor = threshold_cor, 
                      rad.in = rad.in, 
                      file_path = file_path, 
                      color_x_data = color_x_data,
                      color_y_data = color_y_data,
                      cor.method = cor.method)
  
  cca_res[["orig_data"]] <- merge(X, Y, by = "row.names") %>% 
    rename_all(~gsub("__","_", gsub("([[:space:][:punct:]])+$", "", gsub("[[:space:][:punct:]]", "_", .))))
  saveRDS(cca_res, paste0(file_path, "/cca_result.RDS"))
  return(cca_res)
  cat("\nDONE\n\n")
}


plotrCCA <- function(object, x_data_name, y_data_name, threshold_cor = 0.8, 
                     rad.in = rad.in, file_path = file_path, 
                     color_x_data = "red", color_y_data = "blue", 
                     cor.method = "pearson"){

  library(gtools)
  library(ggplot2)
  library(ggrepel)
  library(tidyverse)
  library(patchwork)
  library(ellipse)
  
  #separate datasets (necessary for single correlations)
  dataset_1 <- as.data.frame(object$X_orig_scale)
  dataset_2 <- as.data.frame(object$Y_orig_scale)
  
  #colors for the plot
  colors_plot <- c(color_x_data, color_y_data)
  
  #code from mixOmics::plotVar function
  circle = list()
  circle_center = 0
  circle[[1]] = ellipse::ellipse(circle_center, levels = 1, t = 1)
  circle[[2]] = ellipse::ellipse(circle_center, levels = 1, t = rad.in)
  circle = data.frame(do.call("rbind", circle), Circle = c(rep("Main circle", 
                                                               100), rep("Inner circle", 100)))
  
  if(length(object$cor[object$cor>threshold_cor]) <2){
    cat(str_glue("Less than 2 Canonical variates have a correlation > {threshold_cor}; therefore only plotting the\n first two components. \n"))
    comps = object$cor[1:2]
    df_comb = mixOmics::plotVar(object, legend = c(x_data_name, y_data_name), style = "ggplot2", comp = c(1,2), plot = F) %>%
      dplyr::mutate(distance_to_inner_circle = sqrt((x - circle_center)^2 + (y - circle_center)^2)) %>% 
      dplyr::mutate(outer_circle = ifelse(distance_to_inner_circle > rad.in, rownames(.), "")) 
    
    p = ggplot(df_comb, aes(x = x, y = y, color = Block))
    
    for (j in c("Main circle", "Inner circle")) {
      p = p + geom_path(data = subset(circle, Circle == 
                                        j), aes_string(x = "x", y = "y"), color = "grey30")
    }
    
    p = p + 
      geom_point() +
      scale_color_manual(values=colors_plot) +
      geom_text_repel(aes(label = outer_circle)) +
      labs(title = "Correlation Circle Plot", 
           x = substitute(paste("Component 1: ",rho, " = ", j), 
                          list(j = round(object$cor[1],2))), 
           y = substitute(paste("Component 2: ",rho, " = ", j), 
                          list(j = round(object$cor[2],2)))) + 
      theme_Publication()
    
    #Credits: https://biodatamining.biomedcentral.com/articles/10.1186/1756-0381-5-19
    grDevices::cairo_pdf(filename = paste0(file_path, "/Circleplot_dim_1_2.pdf"), width = 16, height = 12)
    base::plot(p)
    grDevices::dev.off()
    
    grDevices::cairo_pdf(filename = paste0(getwd(), "/Manuscript/Main_Figures/Figure_4.pdf"), width = 16, height = 12)
    base::plot(p)
    grDevices::dev.off()  
    object[["pf4"]] <- p #Small line to add Figure 4 in the cca_res object
    # Filter points from dataset_1 (e.g., parameter A)
    dataset_1_points <- df_comb %>%
      dplyr::filter(Block == x_data_name & distance_to_inner_circle > rad.in) %>% 
      dplyr::select(x, y, names) %>% 
      dplyr::rename(x1 = x, y1 = y, var_dataset1 = names) 
    
    # Filter points from dataset_2 (e.g., parameter C)
    dataset_2_points <- df_comb %>%
      dplyr::filter(Block == y_data_name & distance_to_inner_circle > rad.in) %>% 
      dplyr::select(x, y, names) %>% 
      dplyr::rename(x2 = x, y2 = y, var_dataset2 = names)
    
    # Create a pairwise combination of points from different datasets that are closely positioned
    pairs_oi <- dataset_1_points %>%
      tidyr::crossing(dataset_2_points) %>%
      dplyr::mutate(cosine_angle = (x1 * x2 + y1 * y2) / (sqrt(x1^2 + y1^2)*sqrt(x2^2 + y2^2))) %>%
      dplyr::mutate(inner_product = sqrt(x1^2 + y1^2)*sqrt(x2^2 + y2^2)*cosine_angle) %>% 
      dplyr::mutate(abs_inner_product = abs(inner_product)) %>%
      .[order(as.numeric(-.$abs_inner_product)),] %>% 
      dplyr::select(var_dataset1, x1, y1, var_dataset2, x2, y2, cosine_angle, inner_product, abs_inner_product)
    
    #save the candidate combinations between the two panels
    rio::export(pairs_oi, paste0(file_path, "/Candidate_combinations_rCCA.xlsx"), rowNames=T)
    rio::export(pairs_oi, paste0(getwd(), "/Intermediate_files/Candidate_combinations_rCCA.xlsx"), rowNames=T)

  } else{
    comps = object$cor[object$cor>threshold_cor]
    total_number <- length(comps)
    combination_size <- 2
    cat(str_glue("The number of canonical variates that have a correlation > {threshold_cor}, is {length(comps)}; \nPlotting these components. \n"))
    
    # Generate and loop through combinations
    comb_canvar <- gtools::combinations(total_number, combination_size)
    p_comb = pairs_l = list()
    for (i in 1:nrow(comb_canvar)) {
      # i=1
      current_combination <- comb_canvar[i, ]
      comp_cor <- object$cor[current_combination]
      cat(str_glue("\nCombination", i, ": ", paste(current_combination, collapse = " "), ". \n Correlations are: {round(comp_cor[1],2)} and {round(comp_cor[2],2)} ", "\n\n"))
      
      df_comb = mixOmics::plotVar(object, legend = c(x_data_name, y_data_name), style = "ggplot2", comp = current_combination, plot = F) %>%
        dplyr::mutate(distance_to_inner_circle = sqrt((x - circle_center)^2 + (y - circle_center)^2)) %>% 
        dplyr::mutate(outer_circle = ifelse(distance_to_inner_circle > rad.in, rownames(.), "")) 
      
      p = ggplot(df_comb, aes(x = x, y = y, color = Block))
      
      
      for (j in c("Main circle", "Inner circle")) {
        p = p + geom_path(data = subset(circle, Circle == 
                                          j), aes_string(x = "x", y = "y"), color = "grey30")
      }
      
      p = p + 
        geom_point() +
        scale_color_manual(values=colors_plot) +
        geom_text_repel(aes(label = outer_circle)) +
        labs(title = "Correlation Circle Plot", 
             x = bquote("Component"~.(current_combination[1])*":"~italic(rho) == .(round(comp_cor[1],2))), 
             y = bquote("Component"~.(current_combination[2])*":"~italic(rho) == .(round(comp_cor[2],2)))) + 
        theme_Publication()
      
      #saving the plot with the dimensions
      grDevices::cairo_pdf(filename = paste0(file_path, "/Circleplot_dim_",paste(current_combination, collapse = "_"), ".pdf"), width = 16, height = 12)
      base::plot(p)
      grDevices::dev.off()
      
      object[["pf4"]] <- p #Small line to add Figure 4 in the cca_res object
      
      # Filter points from dataset_1 (e.g., parameter A)
      dataset_1_points <- df_comb %>%
        dplyr::filter(Block == x_data_name & distance_to_inner_circle > rad.in) %>% 
        dplyr::select(x, y, names) %>% 
        dplyr::rename(x1 = x, y1 = y, var_dataset1 = names)# %>% 
      
      # Filter points from dataset_2 (e.g., parameter C)
      dataset_2_points <- df_comb %>%
        dplyr::filter(Block == y_data_name & distance_to_inner_circle > rad.in) %>% 
        dplyr::select(x, y, names) %>% 
        dplyr::rename(x2 = x, y2 = y, var_dataset2 = names)# %>% 
      
      # Create a pairwise combination of points from different datasets that are closely positioned
      pairs_l[[i]] <- dataset_1_points %>%
        tidyr::crossing(dataset_2_points) %>%
        dplyr::mutate(cosine_angle = (x1 * x2 + y1 * y2) / (sqrt(x1^2 + y1^2)*sqrt(x2^2 + y2^2))) %>%
        dplyr::mutate(inner_product = sqrt(x1^2 + y1^2)*sqrt(x2^2 + y2^2)*cosine_angle) %>% 
        dplyr::mutate(abs_inner_product = abs(inner_product)) %>%
        .[order(as.numeric(-.$abs_inner_product)),] %>% 
        dplyr::mutate(dim_pair = paste(current_combination, collapse = "_"),
                      comb_cor = comp_cor[1] + comp_cor[2])
    }
    
    
    #Combinations on which a scatter plot is done with p-value calculation
    pairs_oi <- do.call(rbind, pairs_l) %>% 
      dplyr::arrange(var_dataset1, var_dataset2, desc(comb_cor), desc(abs_inner_product)) %>% 
      dplyr::distinct(var_dataset1, var_dataset2, .keep_all = T) %>% 
      .[order(as.numeric(-.$abs_inner_product)),] %>% 
      dplyr::select(var_dataset1, x1, y1, var_dataset2, x2, y2, cosine_angle, inner_product, abs_inner_product)

    #save the candidate combinations between the two panels
    rio::export(pairs_oi, paste0(file_path, "/Candidate_combinations.xlsx"), rowNames=F)
  }
  object$pairs_oi <- pairs_oi
  return(object)
}

# summarize data
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  return(data_sum)
}
