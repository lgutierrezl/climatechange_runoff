# ====================================================================================
# =============== Extract runoff data and plot by climate change model ===============
# by: Gutierrez Lope Leonardo | lgutierrezlf@gmail.com ===============================
# ====================================================================================

library(ggplot2)
library(raster)
library(cowplot)
library(fields)

source('src/01_plot.R')

# 01. Import data ----------
list_nc <- list.files('data/raw/nc/', pattern = 'median', full.names = T)[-5]
mod_cc      <- raster::stack(list_nc)
names(mod_cc) <- c('PRS', 'ACC','HAD','MPI')
sbs_shp <- shapefile('data/raw/shp/C_aporte.shp')
nom_eps<- sort(unique((sbs_shp$nomeps)))

# 02. Initial extract data and boxplot ----------
tab_stk <- list()
for (i in seq_along(nom_eps)) {
  DEP_sel  <- sbs_shp[sbs_shp$nomeps == nom_eps[i],]
  mod_crop <- crop(mod_cc,rgeos::gBuffer(DEP_sel, width=0.1, byid=TRUE))
  mod_mask <- crop(mask(mod_cc,DEP_sel),DEP_sel)
  REG_DEPf <- na.omit(as.data.frame(mod_mask, xy = TRUE))
  row.names(REG_DEPf) <- NULL
  REG_DEPc    <- data.frame(tidyr::pivot_longer(REG_DEPf, !c(x,y), names_to = "MOD", values_to = "RNF"))
  
  REG_DEPc[,3] <- factor(REG_DEPc[,3], levels = c('PRS', 'ACC', 'HAD', 'MPI'), 
                         labels = c('Presente (1981-2019)','ACCESS 1.0 (2035-2060)', 'HadGEM2-ES (2035-2060)', 'MPI-ESM-LR (2035-2060)'))
  
  list_rsmp <- list()
  for (k in 1:nlayers(mod_crop)) {
    mod_ccsel <- mod_crop[[k]]
    
    x1 <- extent(DEP_sel)[1]
    x2 <- extent(DEP_sel)[2]
    
    zone1 <- ifelse(x2<=-78&x1>-84 | abs(-81-mean(c(x1, x2))) < 3, 17, NA)
    zone2 <- ifelse(x2<=-72&x1>-78 | abs(-75-mean(c(x1, x2))) < 3, 18, NA)
    zone3 <- ifelse(x2<=-66&x1>-72 | abs(-69-mean(c(x1, x2))) < 3, 19, NA)
    utm_z = CRS(paste0('+proj=utm +zone=', max(na.omit(c(zone1, zone2, zone3))[1]),' +south +datum=WGS84 +units=m +no_defs'))
    
    mod_cs_utm <- projectRaster(mod_ccsel, crs=utm_z)
    mod_cs_agg <- aggregate(mod_cs_utm, fact=12, fun=median)
    RPCOP_DF <- data.frame(rasterToPoints(mod_cs_agg))
    EST_SP <- SpatialPointsDataFrame(coords = RPCOP_DF[,1:2], data = RPCOP_DF) 
    proj4string(EST_SP) = utm_z
    
    x1 <- floor(extent(EST_SP)[1]/1000)*1000
    x2 <- ceiling(extent(EST_SP)[2]/1000)*1000
    y1 <- floor(extent(EST_SP)[3]/1000)*1000
    y2 <- ceiling(extent(EST_SP)[4]/1000)*1000
    
    grd_r <-  raster(xmn= x1, xmx= x2, ymn= y1, ymx= y2, resolution=1000,crs=as.character(utm_z))
    
    ESTAC_UTM <- na.omit(EST_SP)
    TPS_LIST  <- Tps(coordinates(ESTAC_UTM), data.frame(ESTAC_UTM)[,3])
    TPS_GRID <- raster::interpolate(grd_r, TPS_LIST)
    resamp_wgs <- projectRaster(TPS_GRID, crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    resamp_wgs <- crop(mask(resamp_wgs,DEP_sel),DEP_sel)
    list_rsmp[[k]] <- resamp_wgs
    rm(grd_r, TPS_GRID, resamp_wgs, TPS_LIST, ESTAC_UTM, EST_SP)
  }
  list_mnc <- raster::stack(list_rsmp)
  names(list_mnc) <- names(mod_crop)
  
  writeRaster(list_mnc[[1]], paste0('data/processed/tif/',sprintf('%02d',i),'_', nom_eps[i],'_PRS.tif'), 'GTiff', overwrite=TRUE)
  writeRaster(list_mnc[[2]], paste0('data/processed/tif/',sprintf('%02d',i),'_', nom_eps[i],'_ACC.tif'), 'GTiff', overwrite=TRUE)
  writeRaster(list_mnc[[3]], paste0('data/processed/tif/',sprintf('%02d',i),'_', nom_eps[i],'_HAD.tif'), 'GTiff', overwrite=TRUE)
  writeRaster(list_mnc[[4]], paste0('data/processed/tif/',sprintf('%02d',i),'_', nom_eps[i],'_MPI.tif'), 'GTiff', overwrite=TRUE)
  
  shapefile(DEP_sel, paste0('data/processed/shp/',sprintf('%02d',i),'_', nom_eps[i],'.shp'), overwrite=TRUE)
  REG_DEPc$MOD <- factor(REG_DEPc$MOD, levels = c('Presente (1981-2019)','ACCESS 1.0 (2035-2060)', 'HadGEM2-ES (2035-2060)', 'MPI-ESM-LR (2035-2060)'),
                         labels = c('Presente','ACCESS 1.0', 'HadGEM2-ES', 'MPI-ESM-LR'))
  
  barras1 <- fig_abs(datab = REG_DEPc)

  tab_val <- data.frame(mod = aggregate(REG_DEPc[,'RNF'], by = list(MOD = REG_DEPc$MOD), FUN = mean)[,1],
                        mean = aggregate(REG_DEPc[,'RNF'], by = list(MOD = REG_DEPc$MOD), FUN = mean)[,2],
                        median = aggregate(REG_DEPc[,'RNF'], by = list(MOD = REG_DEPc$MOD), FUN = median)[,2],
                        var = 'abs') 

  for (j in 4:6) {REG_DEPf[,j+3] <- (REG_DEPf[,j]-(REG_DEPf[,3]+1))/(REG_DEPf[,3]+1)*100} 
  
  rnf_cc <- REG_DEPf[,c(1:2,7:9)]
  names(rnf_cc)[3:5] <- names(REG_DEPf)[4:6]
  rnf_cc    <- data.frame(tidyr::pivot_longer(rnf_cc, !c(x,y), names_to = "MOD", values_to = "RNF"))
  
  rnf_cc[,3] <- factor(rnf_cc[,3], levels = c('ACC', 'HAD', 'MPI'), 
                       labels = c('ACCESS 1.0', 'HadGEM2-ES', 'MPI-ESM-LR'))

  barras2 <- fig_chn(datab = rnf_cc)

  map_hist <- ggdraw() +
    draw_plot(barras2, x = 0, y = 0, width = 1, height = 0.4) +
    draw_plot(barras1, x = 0, y = 0.4, width = 1, height = .6)
  
  tab_cc <- data.frame(mod = aggregate(rnf_cc[,'RNF'], by = list(MOD = rnf_cc$MOD), FUN = mean)[,1],
                       mean = aggregate(rnf_cc[,'RNF'], by = list(MOD = rnf_cc$MOD), FUN = mean)[,2],
                       median = aggregate(rnf_cc[,'RNF'], by = list(MOD = rnf_cc$MOD), FUN = median)[,2], 
                       var = 'cambio')
  tab_exp <- rbind(tab_val, tab_cc)
  tab_exp$eps <- nom_eps[i]
  tab_stk[[i]] <- tab_exp

  ggsave(plot = map_hist, paste0("figs/",sprintf('%02d',i),'_', nom_eps[i],".png"),
         units = "mm", width = 97, height = 140, dpi = 900)
}

tab_df <- do.call(rbind, tab_stk)
write.csv(tab_df[,c(5,1,4,2:3)], 'data/processed/tab2boxplot.csv', row.names = F)
