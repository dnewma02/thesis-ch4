library(raster)
library(caret)
library(tidyverse)
library(ithir)
library(sp)
library(rgdal)
library(doParallel)
library(onsoilsurvey)
library(kernlab)
library(snowfall)
library(dplyr)
library(quantregForest)
library(compositions)

## configure work space
wd <- "path/to/data"
source("oss_seq_vif.R")
out_dir <- paste(wd,"Outputs/", sep="")
data_dir <- paste0(wd,'data/', sep='')
input_dir <- paste0(wd, 'covariates/')
setwd(data_dir)

spatial_prediction = F

# Containers for performance metrics
if(file.exists(paste0(out_dir, "performance.csv"))){
  performance_data = read.csv(file=paste0(out_dir, "performance.csv"), header=TRUE, sep=",")
} else {
  performance_data = data.frame()
}

# Set the depth increments -> GlobalSoilMap.net standard depths
depths = c(0, 5, 15, 30, 60, 100)
vars = c("Sand", "Silt", "Clay")
LSPs = c('Anisotropy', 'DME','Eastness','MeanCurvature','Northness',
         'ProfileCurvature','Roughness','Slope','TanCurvature')

## Start experiment
inp_dir_methods = c(
  paste0(input_dir, 'mono/'),            # 1
  paste0(input_dir, 'heterogeneous/'),   # 2
  paste0(input_dir, 'heterogeneous/'),   # 3
  paste0(input_dir, 'homogeneous/'),     # 4
  paste0(input_dir, 'homogeneous/'),     # 5
  paste0(input_dir, 'homogeneous/'),     # 6
  paste0(input_dir, 'conceptual/')       # 7
)

# loop through mono scale and heterogeneous scale input data
for(method in 1:length(inp_dir_methods)){
  if(method==5 | method==6 | method==7 | method==2 | method==3){next} # skip v, v+g, PCA, RFE, and VIF stuff
  method_dir <- inp_dir_methods[method]
  raster_paths = list.files(path = method_dir, pattern = ".tif", full.names = T)
  if(method==1){method_flag <- "mono"}
  if(method==2){
    method_flag <- "het-z"
    raster_paths = raster_paths[endsWith(raster_paths, 'z.tif')]
  }
  if(method==3){
    method_flag <- "het-zg"
    raster_paths = raster_paths[endsWith(raster_paths, 'z.tif') | endsWith(raster_paths, 'g.tif')]
  }
  if(method==4 | method==5 | method==6){
    setwd(method_dir)
    raster_paths = list.files(path = method_dir, pattern = ".tif", full.names = F)
    hom_lsp_paths = data.frame(
      Anisotropy = raster_paths[startsWith(raster_paths, LSPs[1])],
      DME = raster_paths[startsWith(raster_paths, LSPs[2])],
      Eastness = raster_paths[startsWith(raster_paths, LSPs[3])],
      MeanCurvature = raster_paths[startsWith(raster_paths, LSPs[4])],
      Northness = raster_paths[startsWith(raster_paths, LSPs[5])],
      ProfileCurvature = raster_paths[startsWith(raster_paths, LSPs[6])],
      Roughness = raster_paths[startsWith(raster_paths, LSPs[7])],
      Slope = raster_paths[startsWith(raster_paths, LSPs[8])],
      TanCurvature = raster_paths[startsWith(raster_paths, LSPs[9])]
    )
    if(method==4){
      method_flag <- "hom-VIF"
      print("Running VIF elimination...")
      raster_paths = c()
      for(lsp in 1:length(LSPs)){
        lsp_rasters = hom_lsp_paths[,lsp]
        set.seed(17)
        cov_data <- raster::sampleRandom(x=stack(lsp_rasters), size=1000, na.rm=TRUE, sp=TRUE)
        vif_results <- oss.seq_vif(cov_data@data, thresh=1, trace=FALSE, show.R2.vals=TRUE)
        raster_paths = c(raster_paths, paste0(method_dir, vif_results$CovariatesRetained,'.tif'))
        rm(cov_data, vif_results)#, VIF_threshold, cov_df)
      }
      closeAllConnections()
    }
    if(method==5){method_flag <- "hom-T1"}
    if(method==6){method_flag <- "hom-T2"}
  }
  if(method==7){
    method_flag <- "het-c"
    raster_paths = raster_paths[endsWith(raster_paths, 'z.tif')]
  }

  CovariateStack <- stack(raster_paths)

  # Extract covariate values to points
  # Read in soil data csv for SiteIDs, x, and y
  Sites <- unique.data.frame(read.csv(file=paste0(data_dir, "data.csv"), header=TRUE, sep=",")[,c("SiteID","x","y")])

  # we extract 2 times here, once for CU and RF data, once for SVM, GBM and kNN data
  print("Extracting covariate values...")
  coordinates(Sites) <- ~ x + y
  cores <- (detectCores())
  sfInit(parallel=TRUE,cpus=cores)  # This specifies the number of CPUs to use in the virtual cluster.
  sfLibrary(raster)                 # This will force the snowfall package to use the 'raster' package
  sfLibrary(rgdal)

  inSlist <- unstack(CovariateStack)
  names(inSlist) <- names(CovariateStack)
  values <- sfSapply(inSlist,raster::extract,y=Sites)

  sfStop()  # This terminates the virtual cluster

  ## Prepare the regression matrices by joining harmonized data to covariate values
  SiteData<- cbind(Sites@data, values)
  rm(values, inSlist)

  #Removing rows where covariate values are NA (where sampling may be outside the DEM, oops)
  SiteData <- SiteData[complete.cases(SiteData[,c(which(colnames(SiteData) %in% names(CovariateStack)))]),]

  # apply ILR to textures
  ilrvars = c('ILRtex1','ILRtex2')
  csv_Data <- read.csv(file=paste0(data_dir, "data.csv"), header=TRUE, sep=",")
  csv_Data = csv_Data[complete.cases(csv_Data[, c(which(colnames(csv_Data) %in% vars))]),]
  ilrdata = ilr(csv_Data[colnames(csv_Data) %in% vars])
  colnames(ilrdata) = ilrvars
  csv_Data = cbind(csv_Data, ilrdata)

  Soil_Data = SiteData
  for(v in 1:length(ilrvars)){
    var = ilrvars[v]
    # adjust the min and max if required
    vlow <- min(csv_Data[,var],na.rm=TRUE)*1.0
    vhigh <- max(csv_Data[,var],na.rm=TRUE)*1.0
    # remove rows with NAs for target variable
    VarData <- csv_Data[!is.na(csv_Data[,var]), ]

    ## Run spline on data with lambda 0.1 (default and recommended) to generate data for model building
    print("Running EA Spline...")
    DataSpline <- ea_spline(VarData, var.name = var, d = t(depths),lam = 0.1, vlow = vlow, vhigh = vhigh, show.progress=F)
    id <- DataSpline$harmonised[,1]

    # extract the harmonised soil profile data
    df <- DataSpline$harmonised[,c(2:(ncol(DataSpline$harmonised)-1))]; colnames(df)<- paste0(var,colnames(df))
    # ea_spline returns -9999 in some instances, and we need to convert these to NA
    df[df==-9999]<-NA
    # rename column headers
    colnames(df) <- gsub(x=colnames(df), "-", "_")
    colnames(df) <- gsub(x=colnames(df), " cm", "")
    colnames(df) <- gsub(x=colnames(df), var, paste0(var,'_'))

    # Bind with site id so we can attached coordinates
    VarData <- cbind(id,df)
    # Prepare the regression matrices by joining harmonized data to covariate values
    colnames(VarData)[1] = 'SiteID'
    Soil_Data = merge(x=Soil_Data, y=VarData, by='SiteID')
  }
  Soil_Data <- Soil_Data[complete.cases(Soil_Data[,c(which(colnames(Soil_Data) %in% names(CovariateStack)))]),]
  rm(csv_Data, DataSpline, VarData, id)

  # Set Output Directory for method and variable
  if(!dir.exists(paste0(out_dir, method_flag, '/TEX'))){
    dir.create(paste0(out_dir, method_flag, '/TEX'), recursive = T)
  }
  output_dir <- paste0(out_dir, method_flag, '/TEX/')

  ## Set up modelling parameters for caret
  set.seed(17)
  fitControl <- trainControl(
    method = "repeatedcv",
    number=10,
    repeats=5,
    allowParallel = TRUE,
    returnResamp = "all",
    savePredictions = TRUE
  )

  # Set the cores caret is allowed to access for parallel processing
  registerDoParallel(as.integer(cores*0.75)) # memory consumption is concern

  for(i in 1:(length(depths)-1)){
    ## Start with ILR texture 2

    # Column number of the soil variable
    col_number <- which(colnames(Soil_Data) %in% colnames(df)[i])
    depth_str = gsub(x=substring(gsub(x=colnames(df)[i], var, "",), 2), "_", "-")

    # Remove NAs in response variable and create training regression matrices
    Train_Data = Soil_Data[!is.na(Soil_Data[,col_number]),]

    CovariateStack2 = CovariateStack
    if(method==5 || method==6){
      print("Running T1 or T2...")
      top_lsps = c()
      for(lsp in 1:length(LSPs)){
        lsp_rasters = c(hom_lsp_paths[,lsp])
        lsp_rasters = gsub(x=lsp_rasters, '.tif','')
        RF_tune <- data.frame(mtry=c(seq(2,length(lsp_rasters), by = 2)))
        ML_Equation <- as.formula(paste(names(df)[i], "~", paste(lsp_rasters, collapse="+"),sep=""))
        set.seed(17)
        rfo <- train(
          data=Train_Data, ML_Equation,
          method = "qrf",
          tuneGrid = RF_tune,
          importance = TRUE,
          metric='Rsquared',
          maximize=TRUE,
          trControl = fitControl,
          allowParallel=T
        )
        vimp = varImp(rfo, useModel = F, nonpara = T)$importance
        vimp$Covariate <- rownames(vimp)
        vimp <- vimp[order(-vimp$Overall),]
        if(method==5){
          top_lsps = c(top_lsps, vimp$Covariate[1])
        } else {
          top_lsps = c(top_lsps, vimp$Covariate[1:2])
        }
      }
       CovariateStack2 = stack(paste0(method_dir, top_lsps, '.tif'))
    } # New  raster stack

    # Set up the tuning grid for Random Forest/QRF
    RF_tune <- data.frame(mtry=c(seq(2,nlayers(CovariateStack2), by = 2)))

    # Create an equation that will be fed into the train function
    ML_Equation2 <- as.formula(paste(names(df)[i], "~", paste(names(CovariateStack2), collapse="+"),sep=""))

    # Quantile Regression Forest training
    print('Optimizing quantile regression forest...')
    set.seed(17)
    rfo2 <- train(
      data=Train_Data, ML_Equation2,
      method = "qrf",
      tuneGrid = RF_tune,
      importance = TRUE,
      metric='Rsquared',
      maximize=TRUE,
      trControl = fitControl,
      allowParallel=T
    )
    ccc_rf2 <- oss.getCCC(rfo2)

    ## Train the final qrf model using quantreg package to be used for predictions
    print("Training final QRF...")
    set.seed(17)
    rf_final2 <- quantregForest(
      x=Train_Data[,names(CovariateStack2)],
      y=Train_Data[,col_number],
      mtry=ccc_rf2$ccc_best_tune$mtry
    )

    ## Now ILR texture 1

    #Column number of the soil variable
    col_number <- which(colnames(Soil_Data) %in% gsub(x=colnames(df)[i], ilrvars[2], ilrvars[1]))

    #Remove NAs in response variable and create training regression matrices
    Train_Data = Soil_Data[!is.na(Soil_Data[,col_number]),]

    CovariateStack1 = CovariateStack
    if(method==5 || method==6){
      print("Running second T1 or T2")
      top_lsps = c()
      for(lsp in 1:length(LSPs)){
        lsp_rasters = c(hom_lsp_paths[,lsp])
        lsp_rasters = gsub(x=lsp_rasters, '.tif','')
        RF_tune <- data.frame(mtry=c(seq(2,length(lsp_rasters), by = 2)))
        ML_Equation <- as.formula(paste(gsub(x=names(df)[i], ilrvars[2], ilrvars[1]), "~", paste(lsp_rasters, collapse="+"),sep=""))
        set.seed(17)
        rfo <- train(
          data=Train_Data, ML_Equation,
          method = "qrf",
          tuneGrid = RF_tune,
          importance = TRUE,
          metric='Rsquared',
          maximize=TRUE,
          trControl = fitControl,
          allowParallel=T
        )
        vimp = varImp(rfo, useModel = F, nonpara = T)$importance
        vimp$Covariate <- rownames(vimp)
        vimp <- vimp[order(-vimp$Overall),]
        if(method==5){
          top_lsps = c(top_lsps, vimp$Covariate[1])
        } else {
          top_lsps = c(top_lsps, vimp$Covariate[1:2])
        }
      }
      CovariateStack1 = stack(paste0(method_dir, top_lsps, '.tif'))
    } # New  raster stack

    # Set up the tuning grid for Random Forest/QRF
    RF_tune <- data.frame(mtry=c(seq(2,nlayers(CovariateStack1), by = 2)))

    # Create an equation that will be fed into the train function
    ML_Equation1 <- as.formula(paste(gsub(x=names(df)[i], ilrvars[2], ilrvars[1]), "~", paste(names(CovariateStack1), collapse="+"),sep=""))

    # Quantile Regression Forest training
    print('Optimizing quantile regression forest...')
    set.seed(17)
    rfo1 <- train(
      data=Train_Data, ML_Equation1,
      method = "qrf",
      tuneGrid = RF_tune,
      importance = TRUE,
      metric='Rsquared',
      maximize=TRUE,
      trControl = fitControl,
      allowParallel=T
    )
    ccc_rf1 <- oss.getCCC(rfo1)

    ## Train the final qrf model using quantreg package to be used for predictions
    print("Training final QRF...")
    set.seed(17)
    rf_final1 <- quantregForest(
      x=Train_Data[,names(CovariateStack1)],
      y=Train_Data[,col_number],
      mtry=ccc_rf1$ccc_best_tune$mtry
    )

    ## Predict ILR textures 1 and 2, and invert back to original value
    # Get ILR predictions from model and filter by the optimal mtry
    rfo_pred2 = rfo2$pred[rfo2$pred$mtry == ccc_rf2$ccc_best_tune$mtry,]
    rfo_pred2 = rfo_pred2[order(rfo_pred2$rowIndex, rfo_pred2$Resample),]
    colnames(rfo_pred2) = paste0(colnames(rfo_pred2), '2')
    rfo_pred1 = rfo1$pred[rfo1$pred$mtry == ccc_rf1$ccc_best_tune$mtry,]
    rfo_pred1 = rfo_pred1[order(rfo_pred1$rowIndex, rfo_pred1$Resample),]
    colnames(rfo_pred1) = paste0(colnames(rfo_pred1), '1')
    rfo_pred = cbind(rfo_pred1, rfo_pred2)
    # back transform the observations and predictions
    ilr_obs = as.data.frame(ilrInv(data.matrix(rfo_pred[,c(2,7)])))
    colnames(ilr_obs) = paste0(vars, 'Obs')
    ilr_pred = as.data.frame(ilrInv(data.matrix(rfo_pred[,c(1,6)])))
    colnames(ilr_pred) = paste0(vars, 'Pred')
    rfo_pred = cbind(rfo_pred, ilr_obs, ilr_pred)
    # compute quantiles from observed vs predicted data
    fl_sand<- quantreg::rq(rfo_pred$SandObs ~ rfo_pred$SandPred, tau=c(0.05, 0.95))
    fl_silt<- quantreg::rq(rfo_pred$SiltObs ~ rfo_pred$SiltPred, tau=c(0.05, 0.95))
    fl_clay<- quantreg::rq(rfo_pred$ClayObs ~ rfo_pred$ClayPred, tau=c(0.05, 0.95))

    # now get performance metrics
    perf = ccc_rf2$ccc_optimal_model[,-1]
    perf[1,] = NA
    sand = cbind(goof(rfo_pred$SandObs, rfo_pred$SandPred, type = 'DSM'), NA)
    colnames(sand)[c(1,length(sand))] = c('r2', 'final_pars')
    sand = sand[colnames(sand) %in% colnames(perf)]
    col_idx = match(colnames(sand), colnames(perf))
    col_idx = col_idx[!is.na(col_idx)]
    perf[1,col_idx] = sand[1,]
    silt = cbind(goof(rfo_pred$SiltObs, rfo_pred$SiltPred, type = 'DSM'), NA)
    colnames(silt)[c(1,length(silt))] = c('r2', 'final_pars')
    silt = silt[colnames(silt) %in% colnames(perf)]
    perf[2,col_idx] = silt[1,]
    clay = cbind(goof(rfo_pred$ClayObs, rfo_pred$ClayPred, type = 'DSM'), NA)
    colnames(clay)[c(1,length(clay))] = c('r2', 'final_pars')
    clay = clay[colnames(clay) %in% colnames(perf)]
    perf[3,col_idx] = clay[1,]

    sandpi = fl_sand$fitted.values[,2]-fl_sand$fitted.values[,1]
    siltpi = fl_silt$fitted.values[,2]-fl_silt$fitted.values[,1]
    claypi = fl_clay$fitted.values[,2]-fl_clay$fitted.values[,1]
    MPI = data.frame(MPI=c(mean(sandpi), mean(siltpi), mean(claypi)),
                     MPI_SD=c(sd(sandpi), sd(siltpi), sd(claypi)))
    perf = cbind(perf, MPI)

    model_details = data.frame(
      Method=c(method_flag, method_flag, method_flag),
      Variable=vars,
      Depth_cm=c(depth_str, depth_str, depth_str)
    )
    perf = cbind(model_details, perf)
    performance_data = rbind(performance_data, perf)

    ## Apply the models to the covariate stacks
    if(spatial_prediction == T){
      cores=detectCores()
      beginCluster(n=as.integer(cores*0.5)) # avoid memory bottleneck, less i/o on disk

      # predict means and inverse transform
      print(paste0('Predicting means and inverse ILR ' ,method_flag,' ', colnames(df)[i],' cm...'))
      qrfpred2 <- clusterR(CovariateStack2, predict, args = list(model=rf_final2, what=mean), progress="text")
      qrfpred1 <- clusterR(CovariateStack1, predict, args = list(model=rf_final1, what=mean), progress="text")
      invilr = clusterR(stack(qrfpred1, qrfpred2), calc, args=list(fun=ilrInv), progress='text') #calc(stack(qrfpred1, qrfpred2), ilrInv) #clusterR(stack(qrfpred1, qrfpred2), ilrInv)
      writeRaster(invilr$layer.1, file=paste0(output_dir, vars[1], '_', gsub(x=depth_str, '-','_'), '.tif'), format='GTiff', overwrite=TRUE)
      writeRaster(invilr$layer.2, file=paste0(output_dir, vars[2], '_', gsub(x=depth_str, '-','_'), '.tif'), format='GTiff', overwrite=TRUE)
      writeRaster(invilr$layer.3, file=paste0(output_dir, vars[3], '_', gsub(x=depth_str, '-','_'), '.tif'), format='GTiff', overwrite=TRUE)
      rm(qrfpred1, qrfpred2, invilr)

      # predict upper and inverse transform
      print(paste0('Predicting uppers and inverse ILR ' ,method_flag,' ', colnames(df)[i],' cm...'))
      upper2 <- clusterR(CovariateStack2, predict, args = list(model=rf_final2, what=0.95), progress="text")
      upper1 <- clusterR(CovariateStack1, predict, args = list(model=rf_final1, what=0.95), progress="text")
      invilr = clusterR(stack(upper1, upper2), calc, args=list(fun=ilrInv), progress='text') #calc(stack(qrfpred1, qrfpred2), ilrInv) #clusterR(stack(qrfpred1, qrfpred2), ilrInv)
      writeRaster(invilr$layer.1, file=paste0(output_dir, vars[1], '_', gsub(x=depth_str, '-','_'), '_upper.tif'), format='GTiff', overwrite=TRUE)
      writeRaster(invilr$layer.2, file=paste0(output_dir, vars[2], '_', gsub(x=depth_str, '-','_'), '_upper.tif'), format='GTiff', overwrite=TRUE)
      writeRaster(invilr$layer.3, file=paste0(output_dir, vars[3], '_', gsub(x=depth_str, '-','_'), '_upper.tif'), format='GTiff', overwrite=TRUE)
      sandupper = invilr$layer.1
      siltupper = invilr$layer.2
      clayupper = invilr$layer.3
      rm(upper1, upper2, invilr)

      # predict lower and inverse transform
      print(paste0('Predicting lowers and inverse ILR for ' ,method_flag,' ', colnames(df)[i],' cm...'))
      lower2 <- clusterR(CovariateStack2, predict, args = list(model=rf_final2, what=0.05), progress="text")
      lower1 <- clusterR(CovariateStack1, predict, args = list(model=rf_final1, what=0.05), progress="text")
      invilr = clusterR(stack(lower1, lower2), calc, args=list(fun=ilrInv), progress='text') #calc(stack(qrfpred1, qrfpred2), ilrInv) #clusterR(stack(qrfpred1, qrfpred2), ilrInv)
      writeRaster(invilr$layer.1, file=paste0(output_dir, vars[1], '_', gsub(x=depth_str, '-','_'), '_lower.tif'), format='GTiff', overwrite=TRUE)
      writeRaster(invilr$layer.2, file=paste0(output_dir, vars[2], '_', gsub(x=depth_str, '-','_'), '_lower.tif'), format='GTiff', overwrite=TRUE)
      writeRaster(invilr$layer.3, file=paste0(output_dir, vars[3], '_', gsub(x=depth_str, '-','_'), '_lower.tif'), format='GTiff', overwrite=TRUE)
      sandlower = invilr$layer.1
      siltlower = invilr$layer.2
      claylower = invilr$layer.3
      rm(lower1, lower2, invilr)

      #sand pred interval and stats
      print(paste0('Calculating predictions intervals for ' ,method_flag,' ', colnames(df)[i],' cm...'))
      sandpi = abs(sandupper - sandlower)
      writeRaster(sandpi,file=paste0(output_dir, vars[1], '_', gsub(x=depth_str, '-','_'), '_PI.tif'), format="GTiff", overwrite=TRUE)
      pi_stats = data.frame(
        MPI_ras = cellStats(sandpi, stat = 'mean', na.rm=TRUE),
        MPI_SD_ras = cellStats(sandpi, stat = 'sd', na.rm=TRUE)
      )
      perf = cbind(performance1, pi_stats)
      perf$Variable = vars[1]
      rm(sandupper, sandlower, sandpi)
      # silt pred interval and stats
      siltpi = abs(siltupper - siltlower)
      writeRaster(siltpi,file=paste0(output_dir, vars[2], '_', gsub(x=depth_str, '-','_'), '_PI.tif'), format="GTiff", overwrite=TRUE)
      pi_stats = data.frame(
        MPI_ras = cellStats(siltpi, stat = 'mean', na.rm=TRUE),
        MPI_SD_ras = cellStats(siltpi, stat = 'sd', na.rm=TRUE)
      )
      perf = cbind(performance2, pi_stats)
      perf$Variable = vars[2]
      rm(siltpi, siltlower, siltupper)
      # clay pred interval and stats
      claypi = abs(clayupper - claylower)
      writeRaster(claypi,file=paste0(output_dir, vars[3], '_', gsub(x=depth_str, '-','_'), '_PI.tif'), format="GTiff", overwrite=TRUE)
      pi_stats = data.frame(
        MPI_ras = cellStats(claypi, stat = 'mean', na.rm=TRUE),
        MPI_SD_ras = cellStats(claypi, stat = 'sd', na.rm=TRUE)
      )
      performance3 = performance2
      performance3[,-c(1,2,3)] = NA
      perf = cbind(performance3, pi_stats)
      perf$Variable = vars[3]
      #performance_data = rbind(performance_data, perf)
      rm(claypi, claylower, clayupper)
      endCluster()
    }

    print(paste0("QRF texture predictions complete for  ", depth_str, " cm using ", method_flag))

    # Clean temporary objects and save an .RData file
    rm(rfo,vimp,rf_final1, rf_final2)

    save.image(paste0(output_dir,'_TEX_',gsub(x=depth_str, '-','_'),"_DATA.RData"))
    write.csv(performance_data,file=paste0(out_dir, "performance", ".csv"),row.names=FALSE) # constantly overwrite to keep up to date
  }# END OF DEPTH
  closeAllConnections()
} # END OF METHOD
write.csv(performance_data,file=paste0(out_dir, "performance", ".csv"),row.names=FALSE)
