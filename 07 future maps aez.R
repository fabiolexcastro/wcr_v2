
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, usdm, RColorBrewer, raster, spocc, sf, outliers, glue, CoordinateCleaner, tidyverse, rgbif, readxl, xlsx, openxlsx, rnaturalearthdata, rnaturalearth, geodata)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Ghana 
gha0 <- gadm(country = 'GHA', level = 0, path = './tmpr') %>% st_as_sf()
gha1 <- gadm(country = 'GHA', level = 1, path = './tmpr') %>% st_as_sf()

## Rasters
ftre <- terra::rast('./rf/output/future_crops-gcms.tif')

## Species 
spcs <- basename(as.character(dir_ls('./rData')))
lbls <- tibble(specie = c('Theobroma cacao', 'Anacardium occidentale', 'Cocos nucifera', 'Elaeis guineensis', 'Hevea brasiliensis', 'Vitellaria paradoxa', 'Mangifera indica'), common = c('Cocoa', 'Cashew', 'Coconut', 'Oil palm', 'Rubber', 'Shea', 'Mango'))

## Colors 
clrs <- read_csv('./tbl/Types clusters TCG Three Crops Ghana.csv', show_col_types = FALSE)

# Function to make the map  -----------------------------------------------

##
make.map <- function(spce){
  
  # spce <- 'Coconut'

  ## Read the raster
  cat('To process: ', spce, '\n')
  gcms <- c('GISS-E2-1-H', 'INM-CM4-8')
  rstr <- ftre[[grep(paste0(gcms, collapse = '|'), names(ftre))]]
  rstr <- rstr[[grep(spce, names(rstr))]]
  rstr <- rstr[[1:2]]
  
  ## Raster to table 
  tble <- terra::as.data.frame(rstr, xy = T) %>% as_tibble() 
  clor <- filter(clrs, Crop == spce)
  clor <- drop_na(clor)
  
  tble <- tble %>% gather(var, value, -c(x, y)) %>% separate(data = ., col = 'var', into = c('crop', 'gcm', 'year'), sep = '_')
  tble <- tble %>% mutate(gcm = factor(gcm, levels = gcms))
  
  ## Join 
  tble <- full_join(tble, clor, by = c('value' = 'final'))
  tble <- mutate(tble, Type = factor(Type, levels = unique(clor$Type)))
  tble <- drop_na(tble)
  
  ## Colors 
  clrs <- unique(clor$Color)
  names(clrs) <- unique(clor$Type)
  
  ## To make the map
  g.map <- ggplot() +
    geom_tile(data = tble, aes(x = x, y = y, fill = Type)) +
    scale_fill_manual(values = clrs) +
    facet_wrap(~gcm) +
    geom_sf(data = gha0, fill = NA, col = 'grey30') + 
    geom_sf(data = gha1, fill = NA, col = 'grey30') + 
    labs(x = 'Lon', y = 'Lat', fill = 'AEZ') +
    ggtitle(label = glue('{spce}')) +
    coord_sf() + 
    theme_minimal() +
    theme(
      strip.text = element_text(face = 'bold', hjust = 0.5),
      axis.text.x = element_text(size = 6), 
      axis.text.y = element_text(size = 6, angle = 90),
      plot.title = element_text(size = 12, face = 'bold', hjust = 0.5),
      legend.title = element_text(size = 9, face = 'bold'),
      axis.title = element_text(size = 7, face = 'bold'), 
      legend.position = 'bottom'
    ) 
  
  g.map
  out <- glue('./png/maps/run_1/aez/gcms')
  ggsave(plot = g.map, filename = glue('{out}/gcms_map_{spce}_1_5.jpg'), units = 'in', width = 6, height = 5.0, dpi = 300)
  cat('Done!\n')
  
}

## 
map(spcs, make.map)


## GISS-E2-1-H: very hot - dry
## INM-CM4-8: hot - wet




