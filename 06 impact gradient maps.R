
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, gtools, RColorBrewer)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Load raster and vector data

### Current
crnt <- dir_ls('./rf/output', type = 'directory') %>% as.character() %>% dir_ls() %>% dir_ls() %>% grep('results', ., value = T) %>%
  paste0(., '/v1') %>% dir_ls() %>% grep('.tif', ., value = T) %>% as.character() %>% grep('Mixed', ., value = T)
crnt <- crnt[-2]
nmes <- basename(dirname(dirname(dirname(dirname(crnt)))))
crnt <- rast(crnt)
names(crnt) <- glue('{nmes}_crn')

### Future
ftre <- terra::rast('./rf/output/future_modal.tif')

### Vector data
ghan <- gadm(country = 'GHA', level = 0, path = './tmpr')
gha1 <- gadm(country = 'GHA', level = 1, path = './tmpr')

## Options
all_options <- read_csv('../v1/tbl/classesImpGraLimMix.csv', show_col_types = FALSE); unique(all_options$category) 
labelss <- data.frame(value = c(0, 1, 2, 3, 4, 5), category = c('Unsuit', 'cope', 'adjust', 'transform', 'opportunity', 'resilience'))

## Species
spcs <- nmes
names(crnt) <- glue('{spcs}_crn')

# Impact gradient ---------------------------------------------------------

##
make.impg <- function(spc){
  
  ## Specie 
  # spc <- 'Cocoa'
  
  ## To filter the rasters
  crn <- crnt[[grep(spc, names(crnt))]]
  ftr <- ftre[[grep(spc, names(ftre))]]
  
  ## To raster 
  crn <- raster(crn)
  ftr <- raster(ftr)
  
  ## To start
  msk <- crn * 0
  crd_df <- coordinates(crn)
  
  ## To extract the values
  x <- raster::extract(crn, crd_df, cellnumbers = TRUE) %>% as_data_frame()
  ncell <- dplyr::select(x, cells)
  x <- select_(x, names(crn))
  colnames(x) <- 'current'
  
  y <- raster::extract(ftr, crd_df[,c('x', 'y')], cellnumbers = TRUE) %>% as_data_frame()
  y <- select_(y, names(ftr))
  colnames(y) <- 'future'
  
  z <- data.frame(x, y, ncell) %>% as_tibble()
  
  rslts <- left_join(z, all_options, by = c('current', 'future'))
  labls <- as_tibble(labelss) %>% mutate(category = as.character(category))
  final <- left_join(rslts, labls, by = 'category') %>% dplyr::select(value) %>% pull(1)
  
  length(final)
  length(msk)
  hist(final)
  
  rst <- raster::setValues(msk, final)
  rst <- rast(rst)
  names(rst) <- glue('impact-gradient_{spc}')
  plot(rst)
  terra::writeRaster(x = rst, filename = glue('./rf/output/{spc}/run_1/results/v2/impactGradient_mdlSSP370.tif'), overwrite = TRUE)  
  cat('Done!\n')
  return(rst)
  
}

##
impr <- map(spcs, make.impg)
impr <- reduce(impr, c)

# Colors ------------------------------------------------------------------
# Newly suitable: 92 137 88 Expansion (4) #5c8a4d
# AEZ unchanged: 128 188 122 Incremental adaptation (1) #80bc8d   
# AEZ change 229 213 28 Systemic adaptation (2) #e5d51c
# AEZ uncertain 226 151 26 Systemic resilience (5) #e2971a
# Newly unsuitable 201 88 85 Transform (3)  #c95855

# To make the maps --------------------------------------------------------
tble <- terra::as.data.frame(impr, xy = T) %>% 
  as_tibble() %>% 
  gather(var, value, -c(x, y)) %>% 
  separate(data = ., col = 'var', into = c('imp', 'crop'), sep = '_') %>% 
  inner_join(., tibble(value = c(0, 1, 2, 3, 4, 5), class = c('Unsuitable', 'Incremental adaptation', 'Systemic adaptation', 'Transform', 'Opportunities', 'Systemic resilience')), by = 'value') %>% 
  mutate(class = factor(class, levels = c('Unsuitable', 'Incremental adaptation', 'Systemic adaptation', 'Transform', 'Opportunities', 'Systemic resilience')))

## Colors tibble
clrs <- c('white', '#80bc8d', '#e5d51c', '#c95855', '#5c8a4d', '#e2971a')
names(clrs) <- c('Unsuitable', 'Incremental adaptation', 'Systemic adaptation', 'Transform', 'Opportunities', 'Systemic resilience')

## To make the map 

##
make.map <- function(spc){
  
  ## Filter
  # spc <- 'Cashew'
  cat('To process: ', spc, '\n')
  tbl <- filter(tble, crop == spc)
  
  ## To make the map
  g.imp <- ggplot() +
    geom_tile(data = tbl, aes(x = x, y = y, fill = class)) +
    scale_fill_manual(values = clrs) +
    geom_sf(data = st_as_sf(ghan), fill = NA, col = 'grey30') + 
    geom_sf(data = st_as_sf(gha1), fill = NA, col = 'grey30') + 
    labs(x = 'Lon', y = 'Lat', fill = 'AEZ') +
    ggtitle(label = glue('{spc} impact gradient')) +
    coord_sf() + 
    theme_minimal() +
    theme(
      strip.text = element_text(face = 'bold', hjust = 0.5),
      axis.text.x = element_text(size = 6), 
      axis.text.y = element_text(size = 6, angle = 90),
      plot.title = element_text(size = 12, face = 'bold', hjust = 0.5),
      legend.title = element_text(size = 12, face = 'bold'),
      legend.text = element_text(size = 11),
      axis.title = element_text(size = 7, face = 'bold')
    ) 
  
  ## To save the map 
  ggsave(plot = g.imp, filename = glue('./png/maps/run_1/impact/impact_{spc}.jpg'), units = 'in', width = 5, height = 5, dpi = 300)
  
}

## 
map(spcs, make.map)

