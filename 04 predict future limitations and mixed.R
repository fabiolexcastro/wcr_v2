


# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, spocc, gtools, RColorBrewer, geodata, pROC, Hmisc, randomForest, cclust, corrplot, raptr, sf, outliers, glue, CoordinateCleaner, tidyverse, rgbif, readxl, openxlsx, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
fles <- as.character(dir_ls(as.character(dir_ls('./rf/output/Cocoa/run_1', regexp = 'results')), regexp = '.tif$'))
dirs <- dir_ls('./rf/output') %>% as.character()

## Vector data
gha0 <- gadm(country = 'GHA', level = 0, path = './tmpr')
gha1 <- gadm(country = 'GHA', level = 1, path = './tmpr')

# Function to use ---------------------------------------------------------
make.uncr <- function(dir){
  
  ## Proof
  # dir <- dirs[1]
  
  ## To list the files
  cat('To process: ', dir, '\n')
  fls <- dir_ls(dir) %>% as.character() %>% dir_ls()
  spc <- basename(dir)
  
  ## Threshold 
  thr <- fls %>% grep('results', ., value = T) %>% dir_ls(., regexp = '.csv$') %>% read_csv()
  
  ## To calcualte the quantile
  qnt.prb <- as_tibble(rownames_to_column(as.data.frame(quantile(pull(thr[,1]), seq(0, 1, 0.01))))) %>% mutate(rowname = parse_number(rowname)) %>% setNames(c('rowname', 'value'))
  qnt.unc <- as_tibble(rownames_to_column(as.data.frame(quantile(pull(thr[,2]), seq(0, 1, 0.01))))) %>% mutate(rowname = parse_number(rowname)) %>% setNames(c('rowname', 'value'))
  
  thr.prb <- qnt.prb %>% filter(rowname == 5) %>% pull(2)
  thr.unc <- qnt.unc %>% filter(rowname == 5) %>% pull(2)
  
  ## GCMs / Periods
  fles <- dir_ls(grep('results', fls, value = T)) %>% grep('ssp370', ., value = T) %>% as.character()
  gcms <- c("ACCESS-ESM1-5", "EC-Earth3", "INM-CM5-0", "MPI-ESM1-2-HR", "MRI-ESM2-0")
  prds <- c('2021-2040', '2041-2060')
  
  rst <- map(.x = 1:length(gcms), .f = function(g){
    
    rst <- map(.x = 1:length(prds), .f = function(p){
      
      ## To read as a raster file
      cat('Time: ', gcms[g], prds[p], '\n')
      fll <- grep(prds[p], fles, value = T) %>% grep(gcms[g], ., value = T)
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

# To apply the function ---------------------------------------------------
ftre.ccoa <- make.uncr(dir = dirs[1])
ftre.ccoa
terra::writeRaster(x = ftre.ccoa, filename = './rf/output/Cocoa/run_1/cocoa_future-unc.tif', overwrite = TRUE)

ftre.ccoa

# To make the maps --------------------------------------------------------

## Function
make.map <- function(rst, prd, crp){
  
  # rst <- ftre.ccoa
  # prd <- prds[1]
  # crp <- 'Cocoa'
  
  ## To filter the raster
  cat('To start the analysis: ', prd, ' ', crp, '\n')
  ftr <- rst[[grep(prd, names(rst))]]
  
  ## Raster to table 
  tbl <- terra::as.data.frame(ftr, xy = T) %>% as_tibble() %>% gather(var, value, -c(x, y))
  tbl <- tbl %>% separate(data = ., col = 'var', into = c('crop', 'gcm', 'prd'), sep = '_')
  tbl <- mutate(tbl, gcm = factor(gcm, levels = gcms))
  lbl <- tibble(value = c(1:9), class = c('Unsuitable', 'Unsuitable', paste0('Type ', 1:5), 'Limitations', 'Mixed'))
  tbl <- inner_join(tbl, lbl, by = c('value'))
  tbl <- mutate(tbl, class = factor(class, levels = unique(lbl$class)))
  
  ## To make the maps
  clrs <- c('white', brewer.pal(n = 5, name = 'Set2'), 'grey50', '#ffffc8')
  names(clrs) <- c('Unsuitable', paste0('Type ', 1:5), 'Limitations', 'Mixed')
  
  gmp <- ggplot() + 
    geom_tile(data = tbl, aes(x = x, y = y, fill = class)) + 
    facet_wrap(.~gcm) +
    scale_fill_manual(values = clrs) +
    geom_sf(data = st_as_sf(gha1), fill = NA, col = 'grey50') +
    coord_sf() + 
    ggtitle(label = paste0(str_to_title(crp), ' ', prd)) +
    labs(x = 'Lon', y = 'Lat', fill = 'AEZ') +
    theme_minimal() + 
    theme(
      strip.text = element_text(face = 'bold', hjust = 0.5), 
      axis.text.y = element_text(angle = 90, hjust = 0.5, size = 5), 
      axis.text.x = element_text(size = 5), 
      plot.title = element_text(face = 'bold', hjust = 0.5)
    )
  
  ## To save the maps
  ggsave(plot = gmp, filename = glue('./png/maps/run_1/{crp}_future_{prd}.jpg'), units = 'in', width = 9, height = 7, dpi = 300)
  cat('Done!\n')
  
}

## To apply the function
make.map(rst = ftre.ccoa, prd = prds[1], crp = 'Cocoa')
make.map(rst = ftre.ccoa, prd = prds[2], crp = 'Cocoa')


# Current map -------------------------------------------------------------

rstr <- rast(as.character(dir_ls('./rf/output/Cocoa/run_1/results', regexp = 'Mixed_current')))
tble <- terra::as.data.frame(rstr, xy = T)
tble <- as_tibble(tble)
lbls <- tibble(value = c(1:9), class = c('Unsuitable', 'Unsuitable', paste0('Type ', 1:5), 'Limitations', 'Mixed'))
tble <- inner_join(tble, lbls, by = c('final' = 'value'))
tble <- mutate(tble, class = factor(class, levels = unique(lbl$class)))

clrs <- c('white', brewer.pal(n = 5, name = 'Set2'), 'grey50', '#ffffc8')
names(clrs) <- c('Unsuitable', paste0('Type ', 1:5), 'Limitations', 'Mixed')

gmp <- ggplot() + 
  geom_tile(data = tble, aes(x = x, y = y, fill = class)) + 
  scale_fill_manual(values = clrs) +
  geom_sf(data = st_as_sf(gha1), fill = NA, col = 'grey50') +
  coord_sf() + 
  ggtitle(label = paste0(str_to_title(crp), ' ', prd)) +
  labs(x = 'Lon', y = 'Lat', fill = 'AEZ') +
  theme_minimal() + 
  theme(
    strip.text = element_text(face = 'bold', hjust = 0.5), 
    axis.text.y = element_text(angle = 90, hjust = 0.5, size = 5), 
    axis.text.x = element_text(size = 5), 
    plot.title = element_text(face = 'bold', hjust = 0.5)
  )

ggsave(plot = gmp, filename = './png/maps/run_1/Cocoa_baseline_baseline.jpg', units = 'in', width = 9, height = 7, dpi = 300)

load('./rData/Cocoa/run_1/clustereddata.rData')
