
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, usdm, raster, gtools, spocc, sf, outliers, glue, CoordinateCleaner, tidyverse, rgbif, readxl, xlsx, openxlsx, rnaturalearthdata, rnaturalearth, geodata)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Ghana 
ghan <- gadm(country = 'GHA', level = 0, path = './tmpr')

## Bioclimatic Future
fles <- as.character(dir_ls('../v2_allCrops/tif/c6_5m/ssp370/2041-2060'))
fles <- grep('gha', fles, value = T)
fles <- grep('bioc-alls', fles, value = T)
gcms <- str_split(basename(fles), '_') %>% map_chr(2) %>% unique()
prds <- '2041-2060'

## Species 
dirs <- as.character(dir_ls('./rData'))
spcs <- basename(dirs)

# Function to predic the future -------------------------------------------

##
prdc.ftre <- function(dir, gcm, prd){
  
  # dir <- dirs[2]
  # gcm <- gcms[1]
  # prd <- prds[1]
  
  ## Model 
  cat('To process: ', dir, gcm, prd, '\n')
  spc <- basename(dir)
  fls <- dir_ls(dir) %>% dir_ls() %>% grep('run_1', ., value = T) %>% as.character()
  # rff <- grep('rff_dist', fls, value = T)
  rff <- grep('/rflist', fls, value = T)
  load(rff)
  rff <- do.call(randomForest::combine, rflist)
  
  load(fls[grep('/allclasses_swd', fls, value = F)])
  table(allclasses_swd$pb)
  
  ## Cluster numbers
  NumberOfClusters <- max(as.numeric(allclasses_swd$pb)) - 2
  vrs <- allclasses_swd %>% dplyr::select(starts_with('bio')) %>% colnames()
  
  ## Climate layers
  bio <- grep(paste0(gcm, '_'), fles, value = T) 
  bio <- grep(prd, bio, value = T)
  bio <- grep('gha', bio, value = T)
  bio <- grep('bioc-alls', bio, value = T)
  bio <- rast(bio)
  bio <- bio[[vrs]]
  
  ## To extract by mask 
  bio <- terra::crop(bio, ghan)
  bio <- terra::mask(bio, ghan)
  
  ## Make the climate matrix 
  vls <- as.data.frame(bio, xy = TRUE, na.rm = FALSE)
  vls <- drop_na(vls)
  bio <- rast(vls, type = 'xyz')
  vls <- as.data.frame(bio, xy = FALSE, na.rm = FALSE)
  
  ## Predict
  rasterProbs <- terra::predict(rff, vls, type = 'prob')
  
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
  # plot(rasterRFclass)
  # plot(rasterRFprob)
  # plot(rasterRFuncertainty)
  
  cat('To write the final rasters!\n')
  dir_create(glue('./rf/output/{spc}/run_1/results/v2'))
  terra::writeRaster(x = rasterRFclass, filename = glue('./rf/output/{spc}/run_1/results/v2/gha_RF_5Clust_ssp370_{gcm}_{prd}.tif'), overwrite = TRUE)
  terra::writeRaster(x = rasterRFprob, filename = glue('./rf/output/{spc}/run_1/results/v2/gha_RF_5Prob_ssp370_{gcm}_{prd}.tif'), overwrite = TRUE)
  terra::writeRaster(x = rasterRFuncertainty, filename = glue('./rf/output/{spc}/run_1/results/v2/gha_RF_5Unc_ssp370_{gcm}_{prd}.tif'), overwrite = TRUE)
  cat('Done!\n')
  
}

##
map(.x = 1:length(dirs), .f = function(d){
  map(.x = 1:length(gcms), .f = function(g){
    map(.x = 1:length(prds), .f = function(p){
      prdc.ftre(dir = dirs[d], gcm = gcms[g], prd = prds[p])
    })
  })
})

## Mango
map(.x = 4, .f = function(d){
  map(.x = 1:length(gcms), .f = function(g){
    map(.x = 1:length(prds), .f = function(p){
      prdc.ftre(dir = dirs[d], gcm = gcms[g], prd = prds[p])
    })
  })
})
dirs


# Thresholds --------------------------------------------------------------

# Limitations and mixed
dirs <- as.character(dir_ls('./rf/output', type = 'directory'))
dirs <- glue('{dirs}/run_1') %>% as.character()
dirs <- glue('{dirs}/results/v2') %>% as.character()
fles <- dir_ls(dirs) %>% grep('.tif$', ., value = T) %>% as.character() %>% grep('ssp370', ., value = T)
gcms <- fles %>% str_split('_') %>% map_chr(6) %>% unique()
prds <- c('2041-2060')

# Make limitations and mixed ----------------------------------------------

##
make.uncr <- function(spc, gcm){
  
  # spc <- spcs[5]
  # gcm <- gcms[1]
  
  cat('To process: ', spc, ' ', gcm, '\n')
  fls <- grep(paste0('_', gcm, '_'), fles, value = T)
  fls <- grep('gha', fls, value = T)
  fls <- grep(spc, fls, value = T)
  
  ## To read as a raster file
  cls <- rast(grep('Clust', fls, value = T))
  prb <- rast(grep('Prob', fls, value = T))
  unc <- rast(grep('Unc', fls, value = T))
  
  ## To change the names
  names(cls) <- 'clst'
  names(prb) <- 'prob'
  names(unc) <- 'uncr'
  
  ## To classify - limitations
  # no.clusters <- 5
  
  ## To read the thresholds 
  thrs <- grep('thr', as.character(dir_ls(glue('./rData/{spc}/run_1'))), value = T)
  thr.prb <- grep('prb_1.rds', thrs, value = T) %>% readRDS()
  thr.unc <- grep('unc_', thrs, value = T) %>% readRDS()
  
  ## Matrix  
  mtx.prb <- matrix(c(0, thr.prb, 0, thr.prb, 1, 2), ncol = 3, byrow = T)
  mtx.cls <- matrix(c(0.5, 2 + 0.5, 0, 2 + 0.5, 2 + 5 + 0.5, 1), nrow = 2, byrow = T)
  
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
  names(rsl) <- glue('{spc}_{gcm}_2041-2060')
  return(rsl)
  
}



# Cashew ------------------------------------------------------------------
cshw <- map(.x = 1:length(gcms), .f = function(g){
  make.uncr(gcm = gcms[g], spc = 'Cashew')
}) 
cshw <- map(.x = 1:length(cshw), .f = function(g){
  terra::resample(cshw[[g]], cshw[[1]], method = 'near')
}) %>% 
  reduce(., c)

# Cocoa -------------------------------------------------------------------
coco <- map(.x = 1:length(gcms), .f = function(g){
  make.uncr(gcm = gcms[g], spc = 'Cocoa')
}) 
coco <- map(.x = 1:length(coco), .f = function(g){
  terra::resample(coco[[g]], coco[[1]], method = 'near')
}) %>% 
  reduce(., c)

# Coconut -----------------------------------------------------------------
cocn <- map(.x = 1:length(gcms), .f = function(g){
  make.uncr(gcm = gcms[g], spc = 'Coconut')
}) 
cocn <- map(.x = 1:length(cocn), .f = function(g){
  terra::resample(cocn[[g]], cocn[[1]], method = 'near')
}) %>% 
  reduce(., c)

# Oil palm ----------------------------------------------------------------
oilp <- map(.x = 1:length(gcms), .f = function(g){
  make.uncr(gcm = gcms[g], spc = 'Oil palm')
}) 
oilp <- map(.x = 1:length(oilp), .f = function(g){
  terra::resample(oilp[[g]], oilp[[1]], method = 'near')
}) %>% 
  reduce(., c)

# Rubber ------------------------------------------------------------------
rbbr <- map(.x = 1:length(gcms), .f = function(g){
  make.uncr(gcm = gcms[g], spc = 'Rubber')
}) 
rbbr <- map(.x = 1:length(rbbr), .f = function(g){
  terra::resample(rbbr[[g]], rbbr[[1]], method = 'near')
}) %>% 
  reduce(., c)

# Shea --------------------------------------------------------------------
shea <- map(.x = 1:length(gcms), .f = function(g){
  try(expr = {make.uncr(gcm = gcms[g], spc = 'Shea')})
}) 
shea <- map(.x = 1:length(shea), .f = function(g){
  terra::resample(shea[[g]], shea[[1]], method = 'near')
}) %>% 
  reduce(., c)

# Mango -------------------------------------------------------------------
mngo <- map(.x = 1:length(gcms), .f = function(g){
  try(expr = {make.uncr(gcm = gcms[g], spc = 'Mango')})
}) 
mngo <- map(.x = 1:length(mngo), .f = function(g){
  terra::resample(mngo[[g]], mngo[[1]], method = 'near')
}) %>% 
  reduce(., c)

# Compile all the crops ---------------------------------------------------
crps <- list(cshw, coco, cocn, coco, mngo, oilp, rbbr, shea)
crps <- reduce(crps, c)
terra::writeRaster(x = crps, filename = './rf/output/future_crops-gcms.tif', overwrite = TRUE)

# To make the modal  ------------------------------------------------------
nmes <- names(crps)
spcs
calc.mdal <- function(spce){
  cat('To process: ', spce, '\n')
  ftr <- crps[[grep(spce, names(crps))]]
  ftr <- terra::modal(ftr)
  names(ftr) <- glue('{spce}_ftr')
  return(ftr)
}
ftre.mdal <- map(spcs, calc.mdal)
ftre.mdal <- reduce(ftre.mdal, c)
terra::writeRaster(x = ftre.mdal, filename = './rf/output/future_modal.tif', overwrite = TRUE)

