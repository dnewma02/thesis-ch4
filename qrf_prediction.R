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
vars = c("CEC", "pH_H2O", "TOC")#, "Sand", "Silt", "Clay")
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

  # Loop through soil variables
  for(v in 1:(length(vars))){
    var <- vars[v]
    # Set Output Directory for method and variable
    if(!dir.exists(paste0(out_dir, method_flag, '/', var))){
      dir.create(paste0(out_dir, method_flag, '/', var), recursive = T)
    }
    output_dir <- paste0(out_dir, method_flag, '/', var,'/')

    csv_Data <- read.csv(file=paste0(data_dir, "data.csv"), header=TRUE, sep=",")

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
    # beautify column headers
    colnames(df) <- gsub(x=colnames(df), "-", "_")
    colnames(df) <- gsub(x=colnames(df), " cm", "")
    colnames(df) <- gsub(x=colnames(df), var, paste0(var,'_'))

    # Bind with site id so we can attached coordinates
    VarData <- cbind(id,df)
    # Prepare the regression matrices by joining harmonized data to covariate values
    colnames(VarData)[1] = 'SiteID'
    Soil_Data = merge(x=SiteData, y=VarData, by='SiteID')

    #Removing rows where covariate values are NA (where sampling may be outside the DEM, oops)
    Soil_Data <- Soil_Data[complete.cases(Soil_Data[,c(which(colnames(Soil_Data) %in% names(CovariateStack)))]),]
    rm(csv_Data, DataSpline, VarData, id)

    ## Set up modelling parameters for caret
    ## Use the caret package for training and validation
    ## Then the quantreg package for final predictions for uncertainty maps
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

    # Loop through the depths
    for(i in 1:(length(depths)-1)){ # each loop is ~45 mins
      #Column number of the soil variable
      col_number <- which(colnames(Soil_Data) %in% colnames(df)[i])
      depth_str = gsub(x=substring(gsub(x=colnames(df)[i], var, "",), 2), "_", "-")

      #Remove NAs in response variable and create training regression matrices
      Train_Data = Soil_Data[!is.na(Soil_Data[,col_number]),]

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
          if(method==9){
            top_lsps = c(top_lsps, vimp$Covariate[1])
          } else {
            top_lsps = c(top_lsps, vimp$Covariate[1:2])
          }
        }
        CovariateStack = stack(paste0(method_dir, top_lsps, '.tif'))
      } # New  raster stack

      # Set up the tuning grid for Random Forest/QRF
      RF_tune <- data.frame(mtry=c(seq(2,nlayers(CovariateStack), by = 2)))

      # Create an equation that will be fed into the train function
      ML_Equation <- as.formula(paste(names(df)[i], "~", paste(names(CovariateStack), collapse="+"),sep=""))

      # Quantile Regression Forest training
      print('Optimizing quantile regression forest...')
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

      ccc_rf <- oss.getCCC(rfo)

      ## Train the final qrf model using quantreg package to be used for predictions
      print("Training final QRF...")
      set.seed(17)
      rf_final <- quantregForest(
        x=Train_Data[,names(CovariateStack)],
        y=Train_Data[,col_number],
        mtry=ccc_rf$ccc_best_tune$mtry
      )
      # assign(paste0('rf_final_', colnames(df)[i]), rf_final)

      # Log model performance
      rfo_pred = rfo$pred[rfo$pred$mtry == ccc_rf$ccc_best_tune$mtry,]
      performance = ccc_rf$ccc_optimal_model[,-1]
      PI = data.frame(predict(rf_final, newdata=Train_Data[,names(CovariateStack)], what=c(0.05, 0.95)))
      names(PI) = c("q5", 'q95')
      PI$PI = PI$q95 - PI$q5
      performance$MPI = mean(PI$PI)
      performance$MPI_SD = sd(PI$PI)
      model_details = data.frame(method_flag, var, depth_str)
      colnames(model_details) = c('Method', "Variable", "Depth_cm")
      performance = cbind(model_details, performance)
      performance_data = rbind(performance_data, performance)

      if(spatial_prediction == T){
        ## Apply the models to the covariate stacks
        print(paste0('Predicting (mean, upper, lower) for ' ,method_flag,' ', colnames(df)[i],' cm...'))
        if(!file.exists(paste0(output_dir, colnames(df)[i],".tif"))){
          cores=detectCores()
          beginCluster(n=as.integer(cores*0.75)) # avoid memory bottleneck, less i/o on disk
          qrfpred <- clusterR(CovariateStack, predict, args = list(model=rf_final, what=mean), progress="text")
          writeRaster(qrfpred, file=paste0(output_dir, colnames(df)[i],".tif"), format="GTiff", overwrite=TRUE)
          rm(qrfpred)
          upper <- clusterR(CovariateStack, predict, args = list(model=rf_final, what=0.95), progress="text")
          writeRaster(upper, file=paste0(output_dir, colnames(df)[i],"_upper.tif"), format="GTiff", overwrite=TRUE)
          lower <- clusterR(CovariateStack, predict, args = list(model=rf_final, what=0.05), progress="text")
          writeRaster(lower, file=paste0(output_dir, colnames(df)[i],"_lower.tif"), format="GTiff", overwrite=TRUE)
          pred_interval <- upper-lower
          rm(upper, lower)
          writeRaster(pred_interval,file=paste0(output_dir, colnames(df)[i],"_PI.tif"), format="GTiff", overwrite=TRUE)
          endCluster()

          pi_stats = data.frame(
            MPI_ras = cellStats(pred_interval, stat = 'mean', na.rm=TRUE),
            MPI_SD_ras = cellStats(pred_interval, stat = 'sd', na.rm=TRUE)
          )
          rm(pred_interval)
          #performance = cbind(performance, pi_stats)
          #performance_data = rbind(performance_data, performance)
        } else {
          pred_interval = stack(paste0(output_dir, colnames(df)[i],"_PI.tif"))
          pi_stats = data.frame(
            MPI_ras = cellStats(pred_interval, stat = 'mean', na.rm=TRUE),
            MPI_SD_ras = cellStats(pred_interval, stat = 'sd', na.rm=TRUE)
          )
          rm(pred_interval)
          #performance = cbind(performance, pi_stats)
          #performance_data = rbind(performance_data, performance)
        }
      }

      print(paste0("QRF prediction complete for  ", colnames(df)[i], " cm, using ", method_flag))

      # Clean temporary objects and save an .RData file
      rm(qrfpred,lower,upper, pred_interval)

      save.image(paste0(output_dir,var,"_",colnames(df)[i],"_DATA.RData"))
      write.csv(performance_data,file=paste0(out_dir, "performance", ".csv"),row.names=FALSE) # constantly overwrite to keep up to date
    } # END OF DEPTH
  } # END OF VAR LOOP
  closeAllConnections()
} # END OF METHOD LOOP
write.csv(performance_data,file=paste0(out_dir, "performance", ".csv"),row.names=FALSE)
