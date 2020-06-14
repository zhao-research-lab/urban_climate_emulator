#================================================#
# Prediction: using emulator to predict across 
# CMIP5 models
#================================================#
# This script is for the use in Cheyenne
library(abind)
library(readr)
library(purrr)
library(tidyverse)
library(plyr)
library(stringr)

#=======================#
wd<-getwd()  #directory where the script is in on Cheyenne
cat(paste("********* working directory:",wd,"******************"),sep='\n')
#************************************
RCP <- c('RCP85', 'RCP45')
dataDir <- '/[YOUR DATA DIR]]/'
ESMlist <- read.table(paste0(dataDir,"CMIP5ESMlist.txt"))  #your CIMP5 model name list
for (rcp in RCP) {
  cat(paste("*********work on scenario",rcp,"****************"),sep='\n')
  cat(paste("*********loading emulator......"),sep='\n')
  load(paste0(dataDir, rcp, '/CESM/emulator_tsa_', rcp, '.RData'))    #loading emulator
  for (ESM in ESMlist[,1]) {
    ncdir <-paste0(dataDir, rcp, '/', ESM,'/')
    cat(paste("*********processing",ESM,"****************"),sep='\n')
    cat(paste("*********loading",ESM," new data for prediction****************"),sep='\n')
    load(paste0(ncdir,'urban_newdata_',ESM,'.RData'))    #loading newdata for ESM to predict
    #******************************************************************************
    #        prediction section
    #******************************************************************************
    cat(paste("*********Predicting on",ESM,"****************"),sep='\n')
    # predict TSA_U
    TSA_U_pred_list <- purrr::map2(.x = m_full_list, .y = urban_newdata, ~predict(.x, newdata = .y))  #urban_newdata is your regridded CMIP5 data
    TSA_U_pred_mat <- unlist(TSA_U_pred_list) %>%
      matrix(ncol = 1140, byrow = TRUE)
    # output
    save(TSA_U_pred_mat,file=paste0(dataDir,'TSA_U_pred_',ESM,'_',rcp,'.RData'))
    cat("************** TSA_U prediction done and results saved ******************",sep='\n')
    # predict RH2M_U
    RH2MU_pred_list <- purrr::map2(.x = m_RH2MU_full_list, .y = urban_newdata, ~predict(.x, newdata = .y))
    RH2MU_pred_mat <- unlist(RH2MU_pred_list) %>%
      matrix(ncol = 1140, byrow = TRUE)
    # output
    save(RH2MU_pred_mat,file=paste0(dataDir,'RH2MU_pred_',ESM,'_',rcp,'.RData'))
    cat("************** RH2M_U prediction done and results saved ******************",sep='\n')
    # predict Tmax_U
    TmaxU_pred_list <- purrr::map2(.x = m_TmaxU_full_list, .y = urban_newdata, ~predict(.x, newdata = .y))
    TmaxU_pred_mat <- unlist(TmaxU_pred_list) %>%
      matrix(ncol = 1140, byrow = TRUE)
    # output
    save(TmaxU_pred_mat,file=paste0(dataDir,'TmaxU_pred_',ESM,'_',rcp,'.RData'))
    cat("************** Tmax_U prediction done and results saved ******************",sep='\n')
    # predict Tmin_U
    TminU_pred_list <- purrr::map2(.x = m_TminU_full_list, .y = urban_newdata, ~predict(.x, newdata = .y))
    TminU_pred_mat <- unlist(TminU_pred_list) %>%
      matrix(ncol = 1140, byrow = TRUE)
    # output
    save(TminU_pred_mat,file=paste0(dataDir,'TminU_pred_',ESM,'_',rcp,'.RData'))
    cat("************** Tmin_U prediction done and results saved ******************",sep='\n')
    # memory collection
    rm(urban_newdata, TSA_U_pred_list, TSA_U_pred_mat, TmaxU_pred_list, TmaxU_pred_mat, TminU_pred_list, TminU_pred_mat, RH2MU_pred_list, RH2MU_pred_mat)
    gc()
  }
  # memory collection
  rm(m_full_list, m_TmaxU_full_list, m_TminU_full_list, m_RH2MU_full_list)
  gc()
}
cat("************** ALL DONE !!! ******************",sep='\n')
