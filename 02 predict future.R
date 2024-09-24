
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, spocc, gtools, geodata, pROC, Hmisc, randomForest, cclust, corrplot, raptr, sf, outliers, glue, CoordinateCleaner, tidyverse, rgbif, readx, openxlsx, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Ghana 
ghan <- gadm(country = 'GHA', level = 0, path = './tmpr')

## Bioclimatic Future
fles <- as.character(dir_ls('../01 PREPARE CLIMATE/data/tif/ssp370_gcms'))
fles <- grep('bios', fles, value = T)
gcms <- str_split(basename(fles), '_') %>% map_chr(2) %>% unique()
prds <- c('2021-2040', '2041-2060')

## Species 
dirs <- as.character(dir_ls('./rData'))
lbls <- tibble(specie = c('Theobroma cacao', 'Anacardium occidentale', 'Cocos nucifera', 'Elaeis guineensis', 'Hevea brasiliensis', 'Vitellaria paradoxa', 'Mangifera indica'), common = c('Cocoa', 'Cashew', 'Coconut', 'Oil palm', 'Rubber', 'Shea', 'Mango'))

## Points vif 
# pnts.vif <- read_csv('./rData/')

# Function to predic the future -------------------------------------------
prdc.ftre <- function(dir, gcm, prd){
  
  # dir <- dirs[1]
  # gcm <- gcms[1]
  # prd <- prds[1]
  
  ## Model 
  cat('To process: ', dir, gcm, prd, '\n')
  spc <- basename(dir)
  fls <- dir_ls(dir) %>% dir_ls() %>% grep('run_', ., value = T) %>% as.character()
  # fls <- dir_ls(dir) %>% as.character()
  rff <- grep('rff_dist', fls, value =T )
  load(rff)
  load(fls[grep('/allclasses_swd', fls, value = F)])
  
  ## Cluster numbers
  NumberOfClusters <- max(as.numeric(allclasses_swd$pb)) - 2
  
  ## Get the variables
  vrs <- colnames(dplyr::select(allclasses_swd, starts_with('bio')))
  
  ## Climate layers
  bio <- grep(gcm, fles, value = T) 
  bio <- grep(prd, bio, value = T)
  bio <- rast(bio)
  bio <- bio[[vrs]]
  
  ## Make the climate matrix 
  vls <- as.data.frame(bio, xy = TRUE, na.rm = FALSE)
  vls <- drop_na(vls)
  bio <- rast(vls, type = 'xyz')
  vls <- as.data.frame(bio, xy = FALSE, na.rm = FALSE)
  
  ## Predict
  rasterProbs <- predict(rff, vls, type = 'prob')
  
  ## Fix the rasters
  rasterRF <- rowSums(rasterProbs[,c(3:(NumberOfClusters+2))])
  uncertainty <- apply(rasterProbs, 1, max)  
  
  rasterRFprob <- bio[[1]]
  values(rasterRFprob) <- rasterRF 
  
  rasterRFuncertainty <- bio[[1]]
  values(rasterRFuncertainty) <- uncertainty 
  
  rasterRF <- max.col(rasterProbs, 'first')
  rasterRFclass <- bio[[1]]
  values(rasterRFclass) <- rasterRF
  
  ## Plot the rasters
  plot(rasterRFclass)
  plot(rasterRFprob)
  plot(rasterRFuncertainty)
  
  cat('To write the final rasters!\n')
  terra::writeRaster(x = rasterRFclass, filename = glue('./rf/output/{spc}/run_1/results/RF_5Clust_ssp370_{gcm}_{prd}.tif'), overwrite = TRUE)
  terra::writeRaster(x = rasterRFprob, filename = glue('./rf/output/{spc}/run_1/results/RF_5Prob_ssp370_{gcm}_{prd}.tif'), overwrite = TRUE)
  terra::writeRaster(x = rasterRFuncertainty, filename = glue('./rf/output/{spc}/run_1/results/RF_5Unc_ssp370_{gcm}_{prd}.tif'), overwrite = TRUE)
  cat('Done!\n')
  
}
make.uncr <- function(dir){
  
  cat('To process: ', dir, '\n')
  fls <- dir_ls(dir) %>% as.character()
  spc <- unique(basename(dirname(dirname(dirname(fls)))))
  thr <- thrs %>% filter(crop == spc)
  
  thr.prb <- thr %>% filter(percentile == 5 & type == 'world') %>% pull(prob)
  thr.unc <- thr %>% filter(percentile == 5 & type == 'world') %>% pull(uncr)
  
  # thr.prb <- pull(thr, prob)
  # thr.unc <- pull(thr, uncr) 
  
  rst <- map(.x = 1:length(gcms), .f = function(g){
    
    rst <- map(.x = 1:length(prds), .f = function(p){
      
      ## To read as a raster file
      cat('Time: ', prds[p], '\n')
      fll <- grep(prds[p], fls, value = T) %>% grep(gcms[g], ., value = T)
      cls <- rast(grep('Clust', fll, value = T))
      prb <- rast(grep('Prob', fll, value = T))
      unc <- rast(grep('Unc', fll, value = T))
      
      ## To change the names
      names(cls) <- 'clst'
      names(prb) <- 'prob'
      names(unc) <- 'uncr'
      
      ## To classify - limitations
      no.clusters <- 5
      
      ## Matrix  
      mtx.prb <- matrix(c(0, thr.prb, 0, thr.prb, 1, 2), ncol = 3, byrow = T)
      mtx.cls <- matrix(c(0.5, 2 + 0.5, 0, 2 + 0.5, 2 + no.clusters + 0.5, 1), nrow = 2, byrow = T)
      
      ## To classify the cluster raster
      prb.rcl <- terra::classify(prb, mtx.prb)
      cls.rcl <- terra::classify(cls, mtx.cls)
      
      ## To make the stack 
      names(prb.rcl) <- 'prob_rcl'
      names(cls.rcl) <- 'clst_rcl'
      stk <- c(prb, cls, unc, prb.rcl, cls.rcl)
      
      # To calculate the difference [probability]
      tbl <- as_tibble(terra::as.data.frame(stk, xy = T))
      tbl <- mutate(tbl, dfrn =  prob_rcl - clst_rcl)
      tbl <- mutate(tbl, clst_fnl = ifelse(dfrn %in% c(-1, 2), 8, clst))
      
      tbl %>% dplyr::select(x, y, clst, dfrn, clst_fnl) %>% filter(!clst %in% c(1, 2))
      
      # To calculate the difference [uncertainty]
      tbl <- mutate(tbl, uncr_rcl = ifelse(uncr < thr.unc, 0, 2))
      tbl <- mutate(tbl, clst_unc = ifelse(uncr_rcl == 2 & prob_rcl == 2, 9, clst))
      
      ## To create the final raster
      rsl <- terra::rast(tbl[,c(1, 2, 11)])
      
      ## To finish 
      names(rsl) <- glue('{spc}_{gcms[g]}_{prds[p]}')
      return(rsl)
      
    }) %>% 
      reduce(., c)
    
  }) %>% 
    reduce(., c)
  
  return(rst)
  
}

# To predict  -------------------------------------------------------------
map(.x = 1, .f = function(d){
  map(.x = 1:length(gcms), .f = function(g){
    map(.x = 1:length(prds), .f = function(p){
      prdc.ftre(dir = dirs[d], gcm = gcms[g], prd = prds[p])
    })
  })
})

# Thresholds --------------------------------------------------------------
thrs.fles <- dir_ls('./rf/output', type = 'directory') %>% map(dir_ls) %>% unlist() %>%  grep('run_1', ., value = T) %>% as.character() %>% map(dir_ls) %>% unlist() %>% grep('thresholds_gha-world', ., value = T) %>% as.character()
thrs <- map(thrs.fles, read_csv, show_col_types = F)
thrs <- map_dfr(1:length(thrs), .f = function(i) thrs[[i]] %>% mutate(crop = basename(dirname(dirname(thrs.fles)))[i]))
# thrs <- filter(thrs, type == 'world' & percentile == 5)

# Function to make limitations / mixed ------------------------------------
dirs <- as.character(dir_ls('./rf/output', type = 'directory'))
dirs <- glue('{dirs}/run_1') %>% as.character()
dirs <- glue('{dirs}/results') %>% as.character()
gcms <- c('ACCESS-ESM1', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')
prds <- c('2021-2040', '2041-2060')

# 07 to make limitations / mixed ------------------------------------------
ftre <- map(.x = 1:length(dirs), .f = function(d){make.uncr(dir = dirs[d])})
ftre <- reduce(ftre, c)
terra::writeRaster(x = ftre, filename = './rf/output/predict_future_all-run1_percentile_5.tif', overwrite = TRUE)

# To list the raw files  --------------------------------------------------

fles <- './rf/output/Cocoa/run_1/results' %>% dir_ls(., regexp = '.tif$') %>% as.character()
rs30 <- fles %>% grep('2021', ., value = T) %>% grep('Clust', ., value = T) %>% .[2:6] %>% rast()
rs50 <- fles %>% grep('2040', ., value = T) %>% grep('Clust', ., value = T) %>% .[2:6] %>% rast()
rs30.mdl <- terra::modal(rs30)
rs50.mdl <- terra::modal(rs50)

