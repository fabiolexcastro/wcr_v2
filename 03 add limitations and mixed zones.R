
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, usdm, RColorBrewer, raster, spocc, sf, outliers, glue, CoordinateCleaner, tidyverse, rgbif, readxl, xlsx, openxlsx, rnaturalearthdata, rnaturalearth, geodata)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Ghana 
ghan <- gadm(country = 'GHA', level = 0, path = './tmpr') %>% st_as_sf()
gha1 <- gadm(country = 'GHA', level = 1, path = './tmpr') %>% st_as_sf()

## List the results 
dirs <- dir_ls('./rf/output', type = 'directory') %>% dir_ls() %>% dir_ls() %>% grep('/results', ., value = T) %>% as.character()

## Species 
spcs <- basename(as.character(dir_ls('./rData')))
lbls <- tibble(specie = c('Theobroma cacao', 'Anacardium occidentale', 'Cocos nucifera', 'Elaeis guineensis', 'Hevea brasiliensis', 'Vitellaria paradoxa', 'Mangifera indica'), common = c('Cocoa', 'Cashew', 'Coconut', 'Oil palm', 'Rubber', 'Shea', 'Mango'))

# Function to use ---------------------------------------------------------

##
make.uncr <- function(dir, perc.pr, perc.un){
  
  dir <- dirs[4]
  perc.pr <- 1
  perc.un <- 5
  
  ## To list all the files and get the specie name
  cat('To process: ', dir, '\n')
  fls <- dir_ls(dir) %>% as.character()
  spc <- dirname(dir) %>% dirname() %>% basename()
  
  ## Current rasters 
  fls.crn <- grep('current', fls, value = T)
  crn.prb <- grep('Prob', fls.crn, value = T) %>% rast() %>% setNames('Current_prob')
  crn.cls <- grep('Clust', fls.crn, value = T) %>% rast() %>% setNames('Current_clst')
  crn.unc <- grep('Unc', fls.crn, value = T) %>% rast() %>% setNames('Current_uncr')
  
  plot(crn.prb)
  
  ## To read the points
  pnts <- grep('.csv', fls, value = T) %>% read_csv(., show_col_types = FALSE)
  
  ## To extract the probability and uncertainty values
  occ.clp <- pnts
  
  ## Cluster number
  no.clusters <- max(crn.cls[], na.rm = T) - 2
  occ.clp <- drop_na(occ.clp)
  
  ## To extract the probability
  qnt <- as_tibble(rownames_to_column(as.data.frame(quantile(pull(occ.clp[,'prob.pnts']), seq(0, 1, 0.01))))) %>% mutate(rowname = parse_number(rowname)) %>% setNames(c('rowname', 'value'))
  qnt.prb <- qnt
  qnt.unc <- as_tibble(rownames_to_column(as.data.frame(quantile(pull(occ.clp[,'uncr.pnts']), seq(0, 1, 0.01))))) %>% mutate(rowname = parse_number(rowname))
  colnames(qnt.unc) <- c('rowname', 'value')
  
  ## Quantile table
  qnt <- inner_join(rename(qnt.prb, prb = value), rename(qnt.unc, unc = value))
  qnt <- qnt %>% gather(var, value, -rowname)
  
  ## Make a simple graph 
  g.qnt <- ggplot(data = qnt, aes(x = rowname, y = value, col = var, group = var)) + 
    geom_line(size = 1.2) + 
    labs(x = 'Percentile', y = 'Value', col = '') +
    ggtitle(label = spc) +
    theme_light() + 
    theme(
      plot.title = element_text(size = 16, face = 'bold', hjust = 0.5),
      legend.position = c(0.9, 0.2),
      legend.text = element_text(size = 12), 
      legend.title = element_text(size = 14), 
      axis.text = element_text(size = 11), 
      axis.text.y = element_text(angle = 90, hjust = 0.5),
      axis.title = element_text(size = 12),
      legend.spacing.y = unit(0.5, "cm"),  # Increase vertical spacing between legend items
      legend.key.height = unit(1, "cm") 
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 2)  # Increase line size in legend
    ))
  
  g.qnt
  
  ggsave(plot = g.qnt, filename = glue('./png/graphs/prob_uncr/{spc}.jpg'), units = 'in', width = 5, height = 4, dpi = 300, create.dir = T)
  
  ## To select the thresholds
  perc.pr <- 1
  perc.un <- 5
  
  thr.prb <- filter(qnt.prb, rowname == perc.pr) %>% pull(., 2)
  thr.unc <- filter(qnt.unc, rowname == perc.un) %>% pull(., 2)
  
  saveRDS(object = thr.prb, file = glue('./rData/{spc}/run_1/thr_prb_{perc.pr}.rds'))
  saveRDS(object = thr.unc, file = glue('./rData/{spc}/run_1/thr_unc_{perc.un}.rds'))
  
  ## To make a stack 
  stk <- c(crn.prb, crn.cls, crn.unc)
  stk <- terra::crop(stk, ghan) %>% terra::mask(., ghan)
  
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
  
  ## Label classes
  
  ## Table to raster 
  fnl <- terra::rast(tbl[,c(1, 2, 10)])
  plot(fnl)
  
  ## To write the raster
  dir_create(glue('{dir}/v1'))
  terra::writeRaster(x = fnl, filename = glue('{dir}/v1/RF_5Mixed_current_{perc.pr}_{perc.un}.tif'), overwrite = TRUE)
  cat('Done!\n')
  
}

## 
make.uncr(dir = dirs[1])
make.uncr(dir = dirs[2], perc.pr = 1, perc.un = 5)




