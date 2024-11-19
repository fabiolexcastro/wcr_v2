
## Fabio Alexander Castro Llanos 
## Alliance Bioversity - CIAT 
## November 12th 2024

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, data.table, usdm, classInt, raster, gtools, spocc, ggpubr, sf, outliers, glue, CoordinateCleaner, tidyverse, rgbif, readxl, xlsx, openxlsx, rnaturalearthdata, rnaturalearth, geodata, MetBrewer)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data --------------------------------------------------------------

## Current and future data
fles.crnt <- as.character(dir_ls('../v2_allCrops/tif/wc_5m', regexp = '.tif$'))
fles.crnt <- grep('bioc-all', fles.crnt, value = T)

fles.ftre <- as.character(dir_ls('../v2_allCrops/tif/c6_5m/ssp370/2041-2060', regexp = '.tif$'))
fles.ftre <- grep('gha', fles.ftre, value = T)
fles.ftre <- grep('bioc-alls', fles.ftre, value = T)

##
gcms <- basename(fles.ftre) %>% str_split(., '_') %>% map_chr(2) %>% unique()
# fwrite(as.list(gcms), './gcm.csv')

##
rst.crn <- rast(fles.crnt)

## Vector data 
wrld <- ne_countries(returnclass = 'sf', scale = 50)
ghan <- gadm(country = 'GHA', level = 0, path = './tmpr')
gha1 <- gadm(country = 'GHA', level = 1, path = './tmpr')

# To calculate the ensemble  ----------------------------------------------

##
rst.crn <- rast(fles.crnt)
rst.crn <- crop(rst.crn, ghan)
rst.crn <- mask(rst.crn, ghan)

##
rst.ftr <- map(.x = 1:19, .f = function(i){
  cat('Bioclimatic ', i, '\n')
  rst <- reduce(map(fles.ftre, rast, lyr = i), c)
  rst <- mean(rst)
  names(rst) <- glue('bioc_{i}')
  return(rst)
})
rst.ftr <- reduce(rst.ftr, c)

# To make the difference --------------------------------------------------
rst.crn <- rst.crn[[c(1, 12)]]
rst.ftr <- rst.ftr[[c(1, 12)]]

### Temperature 
dfr.tmp <- rst.ftr[[1]] - rst.crn[[1]]

### Precipitation
dfr.ppt <- rst.ftr[[2]] - rst.crn[[2]]
dfr.ppt <- (dfr.ppt / rst.crn[[2]]) * 100

plot(rst.crn[[2]])
plot(rst.ftr[[2]])

### Final Stack
tbl <- as_tibble(as.data.frame(c(dfr.tmp, dfr.ppt), xy = T))
tbl <- gather(tbl, var, value, -c(x, y))
tbl <- mutate(tbl, var = ifelse(var == 'bioc_1', 'Temperature (°C)', 'Precipitation (mm)'))
tbl <- mutate(tbl, var = factor(var, levels = c('Temperature (°C)', 'Precipitation (mm)')))

# To make the maps for the change -----------------------------------------

## Temperature
tbl.tmp <- filter(tbl, var == 'Temperature (°C)')
brk.tmp <- classInt::classIntervals(var = pull(tbl.tmp, 4), n = 5, style = 'pretty')$brks
tbl.tmp <- tbl.tmp %>% mutate(value_class = findInterval(value, brk.tmp))
tbl.tmp <- inner_join(tbl.tmp, tibble(value_class = 1:5, intr = brk.tmp[1:5]), by = 'value_class')
tbl.tmp <- mutate(tbl.tmp, intr = factor(intr, levels = brk.tmp))

g.tmp <- ggplot() + 
  geom_tile(data = filter(tbl.tmp, var == 'Temperature (°C)'), aes(x = x, y = y, fill = intr)) + 
  scale_fill_manual(values = brewer.pal(n = 5, name = 'YlOrRd')) +
  geom_sf(data = st_as_sf(ghan), fill = NA, col = 'grey30') + 
  geom_sf(data = st_as_sf(gha1), fill = NA, col = 'grey30') + 
  labs(fill = 'Difference\nTemperature (°C)') +
  coord_sf() + 
  theme_void() + 
  theme(
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 14),
    legend.position = 'bottom'
  )

## Precipitation
tbl.ppt <- filter(tbl, var != 'Temperature (°C)')
brk.ppt <- classInt::classIntervals(var = pull(tbl.ppt, 4), n = 5, style = 'pretty')$brks
brk.ppt <- c(-2, -1, 1, 5, 10, 15, 20)
tbl.ppt <- tbl.ppt %>% mutate(value_class = findInterval(value, brk.ppt))
tbl.ppt <- inner_join(tbl.ppt, tibble(value_class = 1:6, intr = brk.ppt[1:6]), by = 'value_class')
tbl.ppt <- mutate(tbl.ppt, intr = factor(intr, levels = brk.ppt))

clrs.ppt <- c("#8C510A", "#D8B365",  'white', "#C7EAE5", "#5AB4AC", "#01665E")

g.ppt <- ggplot() + 
  geom_tile(data = tbl.ppt, aes(x = x, y = y, fill = intr)) + 
  scale_fill_manual(values = clrs.ppt, labels = c('-2', '-1', '~0', '5', '10', '15')) +
  # scale_fill_gradientn(colors = met.brewer(palette_name = 'Isfahan1', n = 8)) +
  geom_sf(data = st_as_sf(ghan), fill = NA, col = 'grey30') + 
  geom_sf(data = st_as_sf(gha1), fill = NA, col = 'grey30') + 
  labs(fill = 'Difference\nPrecipitation (%)') +
  coord_sf() + 
  theme_void() + 
  theme(
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 14),
    legend.position = 'bottom'
  )

g.ppt

## Join both maps into only one 
g.all <- ggpubr::ggarrange(g.tmp, g.ppt, ncol = 2, nrow = 1, labels = c('     A', 'B'), hjust = 0.5, font.label = list(size = 13, color = 'grey30', face = 'bold'))
ggsave(plot = g.all, filename = './png/maps/climate/difference_prec-tasm_v2.jpg', units = 'in', width = 9, height = 5, dpi = 300, create.dir = TRUE)

# To check the dry and the wet GCMs ---------------------------------------
fles.ftre

dfrn <- map(.x = 1:length(fles.ftre), .f = function(i){
  
  ## To read the future raster 
  cat('To process: ', i, '\n')
  fls.ftr <- fles.ftre[i]
  gcm <- basename(fls.ftr) %>% str_split('_') %>% map_chr(2)
  rst.ftr <- rast(fls.ftr, lyr = c(1, 12))
  
  ## To calculate the difference
  tmpr.dfrn <- rst.ftr[[1]] - rst.crn[[1]]
  prec.dfrn <- ((rst.ftr[[2]] - rst.crn[[2]]) / rst.crn[[2]]) * 100
  dfrn <- c(tmpr.dfrn, prec.dfrn)
  names(dfrn) <- c(glue('bioc_1_{gcm}'), glue('bioc_12_{gcm}'))
  return(dfrn)
  
}) %>% 
  reduce(., c)

dfrn.prec <- dfrn[[grep('bioc_12_', names(dfrn), value = FALSE)]]

# Get the percentiles -----------------------------------------------------
prec.05 <- terra::app(x = dfrn.prec, fun = stats::quantile, probs = 0.05, na.rm = T)
prec.95 <- terra::app(x = dfrn.prec, fun = stats::quantile, probs = 0.95, na.rm = T)

prec.pr <- c(prec.05, prec.95)
names(prec.pr) <- c('prec_05', 'prec_95')

## To make the map 

### Tidy the table
tble.pr <- terra::as.data.frame(prec.pr, xy = T, na.rm = T)
tble.pr <- tble.pr %>% as_tibble() %>% gather(var, value, -c(x, y))
tble.pr <- tble.pr %>% mutate(percentile = ifelse(var == 'prec_05', 'Percentile 5', 'Percentile 95'))
tble.pr <- tble.pr %>% mutate(percentile = factor(percentile, levels = c('Percentile 5', 'Percentile 95')))

### To classify the precipitation values
tble.pr

brk.ppt <- classInt::classIntervals(var = pull(tble.pr, 4), n = 5, style = 'pretty')$brks
brk.ppt <- classInt::classIntervals(var = pull(tble.pr, 4), n = 5, style = 'equal')$brks

brk.ppt <- c(-12, -7, -1, 1, 10, 20, 40)

tble.pr <- tble.pr %>% mutate(value_class = findInterval(value, brk.ppt))
tble.pr <- inner_join(tble.pr, tibble(value_class = 1:7, intr = brk.ppt[1:7]), by = 'value_class')
tble.pr <- mutate(tble.pr, intr = factor(intr, levels = brk.ppt))

clrs.ppt <- c("#8C510A", "#D8B365",  'white', "#80CDC1", "#35978F", "#01665E")
# "#8C510A" "#BF812D" "#DFC27D" "#F6E8C3" "#F5F5F5" "#C7EAE5" "#80CDC1" "#35978F" "#01665E"

### Mapping 
gprc.prec <- ggplot() + 
  geom_tile(data = tble.pr, aes(x = x, y = y, fill = intr)) + 
  facet_wrap(~percentile) + 
  scale_fill_manual(values = clrs.ppt, labels = c('< -12', '-7', '~0', '10', '>20')) +
  # scale_fill_gradientn(colors = brewer.pal(n = 5, name = 'BrBG')) +
  geom_sf(data = st_as_sf(ghan), fill = NA, col = 'grey30') + 
  geom_sf(data = st_as_sf(gha1), fill = NA, col = 'grey30') + 
  labs(fill = 'Difference\nPrecipitation (%)') +
  coord_sf() + 
  theme_void() + 
  theme(
    strip.text = element_text(face = 'bold', hjust = 0.5, size = 14),
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 14),
    legend.position = 'bottom'
  )

gprc.prec

## To save the map
ggsave(plot = gprc.prec, filename = './png/maps/climate/difference_prec_percentiles.jpg', units = 'in', width = 9, height = 7, dpi = 300)












