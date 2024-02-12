rm(list=ls())

library(tidyverse)

path_data <- "~/Data_files/MPS/Eduard/Data_transfer/"
setwd(path_data)

st <- rio::import(paste0(path_data, "Pipeline_output/Table_strainsharing.txt"))
st2 <- st %>% dplyr::rename("X2"= "X1","X1" = "X2") %>% dplyr::select(X1, X2, everything())#easiest solution for later issues with one missing value in the matrix
st <- rbind(st, st2); rm(st2)

st <- st %>% dplyr::filter(!grepl(paste(c("V1", "W2", "W3", "W12"), collapse = "|"), X1) & !grepl(paste(c("V1", "W2", "W3", "W12"), collapse = "|"), X2))

### replace the sample swap (MPS_M05 should be MPS_M09 and vise versa)
st$`X1` <- ifelse(st$`X1`=="MPS_M05","MPS_M09",
                  ifelse(st$`X1`=="MPS_M09","MPS_M05", st$`X1`))
st$`X2` <- ifelse(st$`X2`=="MPS_M05","MPS_M09",
                  ifelse(st$`X2`=="MPS_M09","MPS_M05", st$`X2`))

rownames(st) <- paste0 (st$X1, " ", st$X2)

ngd <- rio::import(paste0(path_data, "Pipeline_output/Table_nGD.txt"))
ngd2 <- ngd %>% dplyr::rename("X2"= "X1","X1" = "X2") 
ngd <- rbind(ngd, ngd2); rm(ngd2)

ngd <- ngd %>% dplyr::filter(!grepl(paste(c("V1", "W2", "W3", "W12"), collapse = "|"), X1) & !grepl(paste(c("V1", "W2", "W3", "W12"), collapse = "|"), X2))

### replace the sample swap (MPS_M05 should be MPS_M09 and vise versa)
ngd$`X1` <- ifelse(ngd$`X1`=="MPS_M05","MPS_M09",
                  ifelse(ngd$`X1`=="MPS_M09","MPS_M05", ngd$`X1`))
ngd$`X2` <- ifelse(ngd$`X2`=="MPS_M05","MPS_M09",
                  ifelse(ngd$`X2`=="MPS_M09","MPS_M05", ngd$`X2`))

triad  <- rio::import(paste0(path_data, "Intermediate_files/Triad_numbers.xlsx"))

st1 <- st
#### For loop to do it sample-wise
for (idx in 1: nrow(triad)){
  cat(paste0("Big loop number: ", (idx), "\n"))
  # idx = 13
  triad_oi <- triad[idx, ]
  pre_fmt <- triad_oi$Pre_FMT
  post_fmt <- triad_oi$Post_FMT
  don_fmt <- triad_oi$Donor_Sample_ID
  ngd_d <- ngd %>% dplyr::filter(X1 == post_fmt & X2 == don_fmt) %>% select_if(~all(!is.na(.)))
  ngd_p <- ngd %>% dplyr::filter(X1 == post_fmt & X2 == pre_fmt) %>% select_if(~all(!is.na(.)))
  st_d <- st %>% dplyr::filter(X1 == post_fmt & X2 == don_fmt) %>% select_if(~all(!is.na(.) & . != ""))
  st_p <- st %>% dplyr::filter(X1 == post_fmt & X2 == pre_fmt) %>% select_if(~all(!is.na(.) & . != ""))
  
  ngd_dp <- bind_rows(ngd_d, ngd_p)
  st_dp <- bind_rows(st_d,st_p)
  
  don_pre_T <- names(st_dp %>% select_if(~all(. == TRUE)))
  if (length(don_pre_T) == 0) {
    next
  }
  
  for (idxj in 1:length(don_pre_T)){
    cat(paste0("Small loop number: ", (idxj), "\n"))
    name_sgb <- don_pre_T[idxj]
    sgb_oi <- ngd_dp %>% select(X1, X2,all_of(name_sgb)) %>% 
      dplyr::filter(get({{name_sgb}}) == max(get({{name_sgb}})))
    
    rows_to_update <- row.names(st) %in% c(paste0(sgb_oi$X1, " ", sgb_oi$X2), paste0(sgb_oi$X2, " ", sgb_oi$X1))
    st[rows_to_update, name_sgb] <- FALSE
  }
}
# Also note that the sample swap was fixed in this dataset.
rio::export(st, "Intermediate_files/table_strainsharing_checked.txt")
