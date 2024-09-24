

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, spocc, gtools, geodata, hrbrthemes, pROC, Hmisc, randomForest, cclust, corrplot, raptr, sf, outliers, glue, CoordinateCleaner, tidyverse, rgbif, readx, openxlsx, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Ghana 
ghan <- gadm(country = 'GHA', level = 0, path = './tmpr')

## List the results 
dirs <- dir_ls('./rf/output') %>% map(., dir_ls) %>% unlist() %>% as.character()

## Species 
spcs <- basename(as.character(dir_ls('./rData')))
lbls <- tibble(specie = c('Theobroma cacao', 'Anacardium occidentale', 'Cocos nucifera', 'Elaeis guineensis', 'Hevea brasiliensis', 'Vitellaria paradoxa', 'Mangifera indica'), common = c('Cocoa', 'Cashew', 'Coconut', 'Oil palm', 'Rubber', 'Shea', 'Mango'))

# Function to use ---------------------------------------------------------

make.uncr <- function(dir){
  
  # dir <- dirs[1]
  
  ## To list all the files and get the specie name
  cat('To process: ', dir, '\n')
  fls <- dir_ls(dir) %>% dir_ls() %>% unlist() %>% as.character() %>% grep('', ., value = T) 
  spc <- dirname(dir) %>% basename()
  
  ## Current rasters 
  fls.crn <- grep('current', fls, value = T)
  crn.prb <- grep('Prob', fls.crn, value = T) %>% rast() %>% setNames('Current_prob')
  crn.cls <- grep('Clust', fls.crn, value = T) %>% rast() %>% setNames('Current_clst')
  crn.unc <- grep('Unc', fls.crn, value = T) %>% rast() %>% setNames('Current_uncr')
  
  ## To read the points
  pnts <- dir_ls('./rData') %>% grep(spc, ., value = T) %>% dir_ls() %>% dir_ls()
  load(grep('/clustereddata', pnts, value = T))
  
  ## To get the general points
  pnt <- './rf/output' %>% dir_ls(., regexp = spc) %>% dir_ls(., regexp = 'run_1') %>% dir_ls(regexp = 'results') %>% dir_ls(., regexp = '.csv$') %>% read_csv(., show_col_types = FALSE)
  pnt <- as_tibble(cbind(occ.clp[,1:4],pnt))
  pnt <- filter(pnt, pb == 1)
  
  ## Cluster number
  # no.clusters <- max(crn.cls[], na.rm = T) - 2
  
  ## To extract the probability
  qnt <- as_tibble(rownames_to_column(as.data.frame(quantile(pull(pnt[,5]), seq(0, 1, 0.01))))) %>% mutate(rowname = parse_number(rowname)) %>% setNames(c('rowname', 'value'))
  qnt.prb <- qnt
  qnt.unc <- as_tibble(rownames_to_column(as.data.frame(quantile(pull(pnt[,6]), seq(0, 1, 0.01))))) %>% mutate(rowname = parse_number(rowname))
  colnames(qnt.unc) <- c('rowname', 'value')
  
  ## To select the thresholds
  thr.prb <- filter(qnt, rowname == 10) %>% pull(., 2)
  thr.unc <- filter(qnt, rowname == 10) %>% pull(., 2)
  
  ## To make a stack 
  stk <- c(crn.prb, crn.cls, crn.unc)
  
  ## Raster to table
  tbl <- as_tibble(terra::as.data.frame(stk, xy = T))
  tbl <- mutate(tbl, Current_prob_rcl = ifelse(Current_prob > thr.prb, 2, 0))
  tbl <- mutate(tbl, Current_clst_rcl = ifelse(Current_clst > 2, 1, 0))
  tbl <- mutate(tbl, dfrn = Current_prob_rcl - Current_clst_rcl)  

  ## To reclassify
  tbl <- mutate(tbl, dfrn = ifelse(dfrn %in% c(-1, 2), 8, Current_clst))
  
  ## Mixed
  tbl <- mutate(tbl, Current_unc_rcl = ifelse(Current_uncr > thr.unc, 2, 0))
  tbl <- tbl %>% mutate(final = ifelse(Current_uncr < thr.unc & Current_prob > thr.prb, 9, dfrn))
  
  ## Table to raster 
  fnl <- terra::rast(tbl[,c(1, 2, 9)])
  
  # fnl[crn.unc < thr.unc & crn.prb > thr.prb] <- mxs + 1
  
  ## To write the raster
  terra::writeRaster(x = fnl, filename = glue('{dir}/results/RF_5Mixed_current.tif'), overwrite = TRUE)
  cat('Done!\n')
  
}



