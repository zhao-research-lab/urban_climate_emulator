#================================================#
# Model emulator training: building emulator
# and Cross-validation
#================================================#
# This script is for the use in Cheyenne
library(abind)
library(readr)
library(purrr)
library(tidyverse)
library(plyr)
library(cvTools)
library(stringr)
library(ggplot2)
library(maps)

#=======================#
setwd("[YOUR WORKING DIR]")
wd<-getwd()  #directory where the script is in on cluster
cat(paste("********* working directory:",wd,"******************"),sep='\n')
#************** Model Fitting ****************
RCP <- c('RCP85', 'RCP45')
get_TSAU_rmse <- function(x, formula, K, R) {
  call <- call("lm", formula = formula)
  folds <- cvFolds(nrow(x), K = K, R = R)
  return(mean(cvTool(call, data = x, y = x$TSA_U, folds = folds)))
}
get_TmaxU_rmse <- function(x, formula, K, R) {
  call <- call("lm", formula = formula)
  folds <- cvFolds(nrow(x), K = K, R = R)
  return(mean(cvTool(call, data = x, y = x$TREFMXAV_U, folds = folds)))
}
get_TminU_rmse <- function(x, formula, K, R) {
  call <- call("lm", formula = formula)
  folds <- cvFolds(nrow(x), K = K, R = R)
  return(mean(cvTool(call, data = x, y = x$TREFMNAV_U, folds = folds)))
}
get_RH2MU_rmse <- function(x, formula, K, R) {
  call <- call("lm", formula = formula)
  folds <- cvFolds(nrow(x), K = K, R = R)
  return(mean(cvTool(call, data = x, y = x$RH2M_U, folds = folds)))
}
for (rcp in RCP) {
  cat(paste("********************work on scenario", rcp, "**********************"), sep='\n')
  # Load training data
  cat("************** loading CESM training data ******************", sep='\n')
  load(paste0(wd, '/urban_data_', rcp, '.RData'))  #urban_list, urban_lat, urban_lon
  cat("************** Done: training data loaded ******************",sep='\n')
  # Fit model: full variables plus all two-way interactions
  cat("************** model fitting ******************",sep='\n')
  m_full_list <- purrr::map(urban_list, ~lm(TSA_U ~ TBOT + FSDS  + FLDS + RAIN + PBOT + QBOT + U + V  + ZBOT + as.factor(MONTH) +
                                              as.factor(MONTH):TBOT + as.factor(MONTH):FSDS + as.factor(MONTH):FLDS + 
                                              as.factor(MONTH):RAIN + as.factor(MONTH):PBOT +
                                              as.factor(MONTH):QBOT + as.factor(MONTH):U + as.factor(MONTH):V +
                                              as.factor(MONTH):ZBOT, data = as.data.frame(.x)))
  cat("************** TSA_U model fitting done!******************",sep='\n')
  m_TmaxU_full_list <- purrr::map(urban_list, ~lm(TREFMXAV_U ~ TBOT + FSDS  + FLDS + RAIN + PBOT + QBOT + U + V  + ZBOT + as.factor(MONTH) +
                                              as.factor(MONTH):TBOT + as.factor(MONTH):FSDS + as.factor(MONTH):FLDS + 
                                              as.factor(MONTH):RAIN + as.factor(MONTH):PBOT +
                                              as.factor(MONTH):QBOT + as.factor(MONTH):U + as.factor(MONTH):V +
                                              as.factor(MONTH):ZBOT, data = as.data.frame(.x)))
  cat("************** TREFMXAV_U model fitting done!******************",sep='\n')
  m_TminU_full_list <- purrr::map(urban_list, ~lm(TREFMNAV_U ~ TBOT + FSDS  + FLDS + RAIN + PBOT + QBOT + U + V  + ZBOT + as.factor(MONTH) +
                                                    as.factor(MONTH):TBOT + as.factor(MONTH):FSDS + as.factor(MONTH):FLDS + 
                                                    as.factor(MONTH):RAIN + as.factor(MONTH):PBOT +
                                                    as.factor(MONTH):QBOT + as.factor(MONTH):U + as.factor(MONTH):V +
                                                    as.factor(MONTH):ZBOT, data = as.data.frame(.x)))
  cat("************** TREFMNAV_U model fitting done!******************",sep='\n')
  m_RH2MU_full_list <- purrr::map(urban_list, ~lm(RH2M_U ~ TBOT + FSDS  + FLDS + RAIN + PBOT + QBOT + U + V  + ZBOT + as.factor(MONTH) +
                                                    as.factor(MONTH):TBOT + as.factor(MONTH):FSDS + as.factor(MONTH):FLDS + 
                                                    as.factor(MONTH):RAIN + as.factor(MONTH):PBOT +
                                                    as.factor(MONTH):QBOT + as.factor(MONTH):U + as.factor(MONTH):V +
                                                    as.factor(MONTH):ZBOT, data = as.data.frame(.x)))
  cat("************** RH2M_U model fitting done!******************",sep='\n')
  cat("************** Saving the fitted models ........",sep='\n')
  save(m_full_list, m_TmaxU_full_list, m_TminU_full_list, m_RH2MU_full_list, file=paste0('emulator_tsa_', rcp, '.RData'))
  cat("************** Model saving done ******************",sep='\n')
  # Cross-validation
  cat("************** Cross validation ..........", sep='\n')
  #********************For TSA_U***********************************
  m_full_cv <- map_dbl(urban_list, ~get_TSAU_rmse(as.data.frame(.x),
                                                 TSA_U ~ TBOT + FSDS  + FLDS + RAIN + PBOT + QBOT + U + V  + ZBOT + as.factor(MONTH) +
                                                   as.factor(MONTH):TBOT + as.factor(MONTH):FSDS + as.factor(MONTH):FLDS + 
                                                   as.factor(MONTH):RAIN + as.factor(MONTH):PBOT +
                                                   as.factor(MONTH):QBOT + as.factor(MONTH):U + as.factor(MONTH):V +
                                                   as.factor(MONTH):ZBOT,
                                                   K = 10, R = 1))
  df_m_full_cv <- data.frame(LAT = urban_lat,
                             LON = urban_lon,
                             RMSE = c(m_full_cv ))
  #********************For TREFMXAV_U***********************************
  m_TmaxU_full_cv <- map_dbl(urban_list, ~get_TmaxU_rmse(as.data.frame(.x),
                                                  TREFMXAV_U ~ TBOT + FSDS  + FLDS + RAIN + PBOT + QBOT + U + V  + ZBOT + as.factor(MONTH) +
                                                    as.factor(MONTH):TBOT + as.factor(MONTH):FSDS + as.factor(MONTH):FLDS + 
                                                    as.factor(MONTH):RAIN + as.factor(MONTH):PBOT +
                                                    as.factor(MONTH):QBOT + as.factor(MONTH):U + as.factor(MONTH):V +
                                                    as.factor(MONTH):ZBOT,
                                                  K = 10, R = 1))
  df_m_TmaxU_full_cv <- data.frame(LAT = urban_lat,
                             LON = urban_lon,
                             RMSE = c(m_TmaxU_full_cv ))
  #********************For TREFMNAV_U***********************************
  m_TminU_full_cv <- map_dbl(urban_list, ~get_TminU_rmse(as.data.frame(.x),
                                                   TREFMNAV_U ~ TBOT + FSDS  + FLDS + RAIN + PBOT + QBOT + U + V  + ZBOT + as.factor(MONTH) +
                                                     as.factor(MONTH):TBOT + as.factor(MONTH):FSDS + as.factor(MONTH):FLDS + 
                                                     as.factor(MONTH):RAIN + as.factor(MONTH):PBOT +
                                                     as.factor(MONTH):QBOT + as.factor(MONTH):U + as.factor(MONTH):V +
                                                     as.factor(MONTH):ZBOT,
                                                   K = 10, R = 1))
  df_m_TminU_full_cv <- data.frame(LAT = urban_lat,
                                   LON = urban_lon,
                                   RMSE = c(m_TminU_full_cv ))
  #********************For RH2M_U***********************************
  m_RH2MU_full_cv <- map_dbl(urban_list, ~get_RH2MU_rmse(as.data.frame(.x),
                                                         RH2M_U ~ TBOT + FSDS  + FLDS + RAIN + PBOT + QBOT + U + V  + ZBOT + as.factor(MONTH) +
                                                           as.factor(MONTH):TBOT + as.factor(MONTH):FSDS + as.factor(MONTH):FLDS + 
                                                           as.factor(MONTH):RAIN + as.factor(MONTH):PBOT +
                                                           as.factor(MONTH):QBOT + as.factor(MONTH):U + as.factor(MONTH):V +
                                                           as.factor(MONTH):ZBOT,
                                                         K = 10, R = 1))
  df_m_RH2MU_full_cv <- data.frame(LAT = urban_lat,
                                   LON = urban_lon,
                                   RMSE = c(m_RH2MU_full_cv ))
  rm(list=c("m_full_list", "urban_list", "m_full_cv", "urban_lat", "urban_lon"))
  gc()
}




