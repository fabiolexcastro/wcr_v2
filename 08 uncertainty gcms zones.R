



# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, gtools, RColorBrewer)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

### Vector data
gha0 <- gadm(country = 'GHA', level = 0, path = './tmpr') %>% st_as_sf()
gha1 <- gadm(country = 'GHA', level = 1, path = './tmpr') %>% st_as_sf()

### Future
ftre <- terra::rast('./rf/output/future_modal.tif')

## Species 
spcs <- basename(as.character(dir_ls('./rData')))
lbls <- tibble(specie = c('Theobroma cacao', 'Anacardium occidentale', 'Cocos nucifera', 'Elaeis guineensis', 'Hevea brasiliensis', 'Vitellaria paradoxa', 'Mangifera indica'), common = c('Cocoa', 'Cashew', 'Coconut', 'Oil palm', 'Rubber', 'Shea', 'Mango'))

## Colors 
clrs <- read_csv('./tbl/Types clusters TCG Three Crops Ghana.csv', show_col_types = FALSE)


