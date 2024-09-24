


# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, usdm, spocc, sf, outliers, glue, CoordinateCleaner, tidyverse, rgbif, readxl, xlsx, openxlsx, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

source('../02 MAKE MODEL/FunctionsRFclustering.R')
lbls <- tibble(specie = c('Theobroma cacao', 'Anacardium occidentale', 'Cocos nucifera', 'Elaeis guineensis', 'Hevea brasiliensis', 'Vitellaria paradoxa', 'Mangifera indica'), common = c('Cocoa', 'Cashew', 'Coconut', 'Oil palm', 'Rubber', 'Shea', 'Mango'))

# Functions to use --------------------------------------------------------
make.clea <- function(spce){
  
  # spce <- 'Theobroma cacao'
  
  ## To filter the specie
  cat('To process: ', spce, '\n')
  pnt <- filter(pnts, nombre == spce & pb == 1)
  bck <- filter(pnts, nombre == spce & pb == 0)
  
  ## To clean the coordinates
  cln <- clean_coordinates(x = as.data.frame(pnt), lon = 'Longitude', lat = 'Latitude', species = 'nombre', tests = c('capitals', 'centroids', 'equal', 'zeros', 'institutions', 'seas'))
  cla <- cln[cln$.summary,]
  cla <- dplyr::select(cla, pb:bioc_29)
  cla <- as_tibble(cla)
  
  ## Join with the background 
  cla <- rbind(cla, bck[,-33])
  
  ## Return
  cat('Done!\n')
  return(cla)
  
  
}
make.vifs <- function(tble, spce){
  
  # tble <- pnts.clea
  # spce <- 'Theobroma cacao'
  
  ## To filter the specie
  cat('To process: ', spce, '\n')
  pnt <- filter(tble, nombre == spce)
  occ <- filter(pnt, pb == 1)
  bck <- filter(pnt, pb == 0)
  
  ## To make the VIF 
  vif <- usdm::vifstep(x = as.data.frame(occ[,5:32]), th = 10)
  vrs <- vif@results$Variables
  
  ## To select the variables 
  rsl <- dplyr::select(pnt, pb:Latitude, all_of(vrs))
  
  ## Finish
  cat('Done!\n')
  return(rsl)
  
}
make.clst.occr <- function(tble, spce){
  
  ## To start the process
  cat('To star the process: ', spce, '\n')
  nme <- filter(lbls, specie == spce) %>% pull(2)
  pnt <- filter(tble, nombre == spce & pb == 1)
  
  ## No Forest / No Trees
  no.forest <- 25
  no.trees <- 100
  nVars <- 8
  
  ## Clustering presences 
  occ <- pnt
  occ.mtx <- occ[,5:ncol(occ)]
  occ.dst <- RFdist(occ.mtx, mtry1 = nVars, no.trees, no.forest, addcl1 = T, addcl2 = F, imp = T, oob.prox1 = T)
  occ.ncl <- 5
  occ.lrf <- pamNew(occ.dst$cl1, occ.ncl)
  occ.cls <- hclust(as.dist(occ.dst$cl1), method = 'ward.D2')
  occ.cld <- cbind(pb = as.factor(occ.lrf), occ[2:ncol(occ)])
  occ.clp <- cbind(occ, cluster = occ.lrf) %>% na.omit() %>% as_tibble()
  
  ## To save the results
  dir <- glue('./rData/{nme}/run_1'); dir_create(dir)
  save(occ.mtx, file = glue('{dir}/datRF.rData'))
  save(occ.cls, file = glue('{dir}/clusterdata.rData'))
  save(occ, occ.clp, occ.ncl, occ.lrf, file = glue('{dir}/clustereddata.rData'))
  save(occ.cld, file = glue('{dir}/occ_cld.rData'))
  cat('Done!\n')
  
}
make.rfrs.mdel <- function(tble, crop){
  
  tble <- pnts.vifs
  crop <- 'Theobroma cacao'
  
  ## To start the process
  cat('To start the process: ', crop, '\n')
  nme <- filter(lbls, specie == spce) %>% pull(2)
  
  ## No Forest / No Trees
  no.forest <- 25
  no.trees <- 100
  nVars <- 8
  
  ## Clustering pseudo-absences 
  bck <- filter(tble, pb == 0)
  bck.mtx <- bck[,5:ncol(bck)]
  bck.dst <- RFdist(bck.mtx, mtry1 = nVars, no.trees, no.forest, addcl1 = T, addcl2 = F, imp = T, oob.prox1 = T)
  bck.ncl <- 2
  bck.lrf <- pamNew(bck.dst$cl1, bck.ncl)
  bck.cls <- hclust(as.dist(bck.dst$cl1), method = 'ward.D2')
  bck.cld <- cbind(pb = as.factor(bck.lrf), bck[2:ncol(bck)])
  bck.clp <- cbind(bck, cluster = bck.lrf) %>% na.omit() %>% as_tibble()
  
  ## Read the results of presences 
  fls.occ <- as.character(dir_ls(glue('./rData/{nme}/run_1'), regexp = '.rData'))
  load(grep('/clusterdata', fls.occ, value = T))
  load(grep('/clustereddata', fls.occ, value = T))
  load(grep('/datRF', fls.occ, value = T))
  load(grep('occ_cld', fls.occ, value = T))
  
  clusteredpresdata <- occ.clp
  no.absenceclasses <- 2
  
  presvalue_swd  <- clusteredpresdata[,3:ncol(clusteredpresdata)] %>%
    cbind(pb = (clusteredpresdata$cluster + no.absenceclasses), .) %>%
    na.omit() %>%
    as.data.frame() %>%
    mutate(cluster = cluster + no.absenceclasses)
  presvalue_swd <- dplyr::select(presvalue_swd, pb, bioc_1:bioc_19, bioc_21:bioc_29)
  presvalue_swd <- mutate(presvalue_swd, pb = as.factor(pb))
  classdata <- occ.cld
  
  classdata_2 <- cbind(pb = as.data.frame(classdata)$pb, classdata[,2:ncol(classdata)]) # Background
  
  dim(classdata_2); dim(presvalue_swd)
  # presvalue_swd <- presvalue_swd %>% dplyr::select(-cluster)
  classdata_2 <- dplyr::select(classdata_2, -Longitude, -Latitude)
  
  allclasses_swd <- rbind(classdata_2, presvalue_swd[,1:ncol(classdata_2)])
  unique(allclasses_swd$pb)
  
  ## To save the results
  dir <- glue('./rData/{nme}/run_1'); dir_create(dir)
  save(bck.mtx, file = glue('{dir}/back_datRF.rData'))
  save(bck.cls, file = glue('{dir}/back_clusterdata.rData'))
  save(bck, bck.clp, bck.ncl, bck.lrf, file = glue('{dir}/back_clustereddata.rData'))
  save(allclasses_swd, file = glue('{dir}/allclasses_swd.rData'))
  
  ## To make the random forest
  
  # To make the random forest analysis --------------------------------------
  vrs <- c(paste0('bioc_', 1:19), paste0('bioc_', 21:29))
  model1 <- as.formula(paste('factor(pb) ~', paste(paste(vrs), collapse = '+', sep =' ')))
  rflist <- vector('list', 50) 
  auc <- vector('list', 50)
  NumberOfClusters <- 5
  
  clusteredpresdata
  allclasses_swd
  allclasses_swd <- sample_n(tbl = allclasses_swd, size = nrow(allclasses_swd) * 0.5, replace = FALSE)
  
  samplesize <- round(min(summary(as.factor(allclasses_swd$pb))) / 2, 0) 
  
  for(repe in 1:50){ # 50 bosques
    
    print(repe)
    pressample <- list()
    
    for (i in 1:(NumberOfClusters+no.absenceclasses)){
      
      if(any(i==c(1:no.absenceclasses))) { 
        
        rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), 
                       size = samplesize*NumberOfClusters/2/no.absenceclasses)
      } else {
        rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), size=samplesize)
      }
      pressample[[i]] <- allclasses_swd[rows,] 
    }
    
    species <- na.omit(do.call(rbind, pressample)) 
    head(species)
    Samplesplit <- sample(rownames(species)) 
    
    envtrain <- species[Samplesplit[1:(0.8*nrow(species))],] 
    envtest <- species[Samplesplit[(0.8*nrow(species)):nrow(species)],] 
    
    rfmodel <- randomForest(model1, data = envtrain, ntree = 500, na.action = na.omit, nodesize = 2) 
    
    dir_create(glue('./rf/output/{nme}/run_1/models'))
    save(rfmodel, file = glue('./rf/output/{nme}/run_1/models/RF_', NumberOfClusters, 'Prob_' , 'rep_' ,repe, '.rdata' ,sep=''))
    rflist[[repe]] <- rfmodel
    
    # AUC 
    predicted <- as.numeric(predict(rfmodel, envtest))
    observed <- as.vector(envtest[,'pb'])
    auc[[repe]] <- auc(observed, predicted) 
    rm(rfmodel)
    
    cat(auc[[repe]] ,'\n')
    
  }
  
  auc <- unlist(auc)
  rff <- do.call(randomForest::combine, rflist)
  importance <- as.data.frame(rff$importance)
  
  dir.out <- glue('./rData/{nme}/run_1')
  dir_create(dir.out)
  
  save(rflist, file = paste(dir.out, '/rflist_', NumberOfClusters, '.rdata', sep = ''))
  save(importance, file = paste0(dir.out, '/importanceRF.rData'))
  save(auc, file = paste0(dir.out, '/aucRF_dist.rData'))
  save(rff, file = paste0(dir.out, '/rff_dist.rData'))
  
  # To extract by mask 
  lyr <- rast('../01 PREPARE CLIMATE/data/tif/bioc_all.tif')
  lyr <- terra::crop(lyr, limt)
  lyr <- terra::mask(lyr, limt)
  
  # Predict modell
  climatevalues <- as.data.frame(lyr, xy = T, na.rm = F)
  NumberOfClusters <- 5
  
  pnts.df <- allclasses_swd[,2:ncol(allclasses_swd)]
  
  rasterProbs <- predict(rff, climatevalues[,3:ncol(climatevalues)], type = 'prob') # proximity = T
  ## Hacer el predict para todos los puntos, usando la matriz all_classesswd; y así generar la curva de los percentiles; 
  ## Primer percentile..., y te comparto nuevamente las figuras para que le des una revisada 
  ## Hacer dos curvas, una para todo el Mundo y una para Ghana (dos curvas en el mismo gráfico)
  
  rasterProbs_na <- na.omit(rasterProbs)
  sum_rasterProbs_na <- apply(rasterProbs_na, 1, sum)
  
  rasterRF <- rowSums(rasterProbs[,c(3:(NumberOfClusters+2))])
  uncertainty <- apply(rasterProbs, 1, max)  
  
  rasterRFprob <- lyr[[1]]
  values(rasterRFprob) <- rasterRF 
  
  rasterRFuncertainty <- lyr[[1]]
  values(rasterRFuncertainty) <- uncertainty 
  
  rasterRF <- max.col(rasterProbs, 'first')
  rasterRFclass <- lyr[[1]]
  values(rasterRFclass) <- rasterRF
  
  ## For the points
  rasterProbs_points <- predict(rff, pnts.df, type = 'prob')
  prob.pnts <- rowSums(rasterProbs_points[,3:(NumberOfClusters+2)])
  uncr.pnts <- apply(rasterProbs_points, 1, max)
  prob.pnts <- cbind(prob.pnts, uncr.pnts)
  
  ## To write the raster
  plot(rasterRFclass)
  dir.out <- glue('./rf/output/{nme}/run_2/results')
  dir_create(dir.out)
  terra::writeRaster(rasterRFclass, paste0(dir.out, '/RF_5Clust_current.tif'), overwrite = T)
  terra::writeRaster(rasterRFprob, paste0(dir.out, '/RF_5Prob_current.tif'), overwrite = T)
  terra::writeRaster(rasterRFuncertainty, paste0(dir.out, '/RF_5Unc_current.tif'), overwrite = T)
  write.csv(prob.pnts, paste0(dir.out, '/pnts_prob-clust.csv'), row.names = FALSE)
  cat('Done!\n')
  
}


# Load data ---------------------------------------------------------------
pnts <- read_csv('../02 MAKE MODEL/tbl/points/processed/03 points background.csv')

# To cleaning the points  -------------------------------------------------
pnts.clea <- make.clea(spce = 'Theobroma cacao')

# To make the VIF ---------------------------------------------------------
pnts.vifs <- make.vifs(tble = pnts.clea, spce = 'Theobroma cacao')

# To make the Random Forest -----------------------------------------------






