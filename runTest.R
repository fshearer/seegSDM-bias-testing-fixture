library(raster)
library(sp)
library(devtools)
install_github("SEEG-Oxford/seegSDM", ref="0.1-8")
library(seegSDM)
## This is fork of runABRAID from SEEG-Oxford/seegSDM @ 0.1-8
runTest <- function (name,
                     mode, 
                     disease,
                     admin_extract_mode="random",
                     crop_bias=TRUE, # used with mode "bias"
                     filter_bias=TRUE, # used with mode "bias"
                     use_weights=TRUE,
                     use_temporal_covariates=TRUE) {
  # Get file paths
  occurrence_path <- paste0(disease,"_data/occurrences.csv")
  extent_path <- paste0(disease,"_data/extent.tif")
  supplementary_occurrence_path <- paste0(disease,"_data/supplementary_occurrences.csv")
  admin_path <- list(
    "admin0" <- "admins/admin0.tif",
    "admin1" <- "admins/admin1.tif",
    "admin2" <- "admins/admin2.tif",
    "admin3" <- "admins/admin2.tif" # This one wont be used, but is needed for compatablity with older bits of seegSDM
  )
  water_mask <- "admins/waterbodies.tif"
  disease_type <- list(
    "cchf"="virus",
    "chik"="virus",
    "deng"="virus",
    "hat"="parasite",
    "melio"="bacteria",
    "nwcl"="parasite",
    "nwvl"="parasite",
    "owcl"="parasite",
    "owvl"="parasite",
    "scrub"="bacteria"
  )[[disease]]
  
  all_covs <- list(
    "access"= "covariates/access.tif",
    "c10"=list(
      "2001"="covariates/2001.Class10_Grasslands.5km.Percentage.ABRAID.tif",
      "2002"="covariates/2002.Class10_Grasslands.5km.Percentage.ABRAID.tif",
      "2003"="covariates/2003.Class10_Grasslands.5km.Percentage.ABRAID.tif",
      "2004"="covariates/2004.Class10_Grasslands.5km.Percentage.ABRAID.tif",
      "2005"="covariates/2005.Class10_Grasslands.5km.Percentage.ABRAID.tif",
      "2006"="covariates/2006.Class10_Grasslands.5km.Percentage.ABRAID.tif",
      "2007"="covariates/2007.Class10_Grasslands.5km.Percentage.ABRAID.tif",
      "2008"="covariates/2008.Class10_Grasslands.5km.Percentage.ABRAID.tif",
      "2009"="covariates/2009.Class10_Grasslands.5km.Percentage.ABRAID.tif",
      "2010"="covariates/2010.Class10_Grasslands.5km.Percentage.ABRAID.tif",
      "2011"="covariates/2011.Class10_Grasslands.5km.Percentage.ABRAID.tif",
      "2012"="covariates/2012.Class10_Grasslands.5km.Percentage.ABRAID.tif"
    ),
    "c6"=list(
      "2001"="covariates/2001.Class06_Closed_Shrublands.5km.Percentage.ABRAID.tif",
      "2002"="covariates/2002.Class06_Closed_Shrublands.5km.Percentage.ABRAID.tif",
      "2003"="covariates/2003.Class06_Closed_Shrublands.5km.Percentage.ABRAID.tif",
      "2004"="covariates/2004.Class06_Closed_Shrublands.5km.Percentage.ABRAID.tif",
      "2005"="covariates/2005.Class06_Closed_Shrublands.5km.Percentage.ABRAID.tif",
      "2006"="covariates/2006.Class06_Closed_Shrublands.5km.Percentage.ABRAID.tif",
      "2007"="covariates/2007.Class06_Closed_Shrublands.5km.Percentage.ABRAID.tif",
      "2008"="covariates/2008.Class06_Closed_Shrublands.5km.Percentage.ABRAID.tif",
      "2009"="covariates/2009.Class06_Closed_Shrublands.5km.Percentage.ABRAID.tif",
      "2010"="covariates/2010.Class06_Closed_Shrublands.5km.Percentage.ABRAID.tif",
      "2011"="covariates/2011.Class06_Closed_Shrublands.5km.Percentage.ABRAID.tif",
      "2012"="covariates/2012.Class06_Closed_Shrublands.5km.Percentage.ABRAID.tif"
    ),
    "c7"=list(
      "2001"="covariates/2001.Class07_Open_Shrublands.5km.Percentage.ABRAID.tif",
      "2002"="covariates/2002.Class07_Open_Shrublands.5km.Percentage.ABRAID.tif",
      "2003"="covariates/2003.Class07_Open_Shrublands.5km.Percentage.ABRAID.tif",
      "2004"="covariates/2004.Class07_Open_Shrublands.5km.Percentage.ABRAID.tif",
      "2005"="covariates/2005.Class07_Open_Shrublands.5km.Percentage.ABRAID.tif",
      "2006"="covariates/2006.Class07_Open_Shrublands.5km.Percentage.ABRAID.tif",
      "2007"="covariates/2007.Class07_Open_Shrublands.5km.Percentage.ABRAID.tif",
      "2008"="covariates/2008.Class07_Open_Shrublands.5km.Percentage.ABRAID.tif",
      "2009"="covariates/2009.Class07_Open_Shrublands.5km.Percentage.ABRAID.tif",
      "2010"="covariates/2010.Class07_Open_Shrublands.5km.Percentage.ABRAID.tif",
      "2011"="covariates/2011.Class07_Open_Shrublands.5km.Percentage.ABRAID.tif",
      "2012"="covariates/2012.Class07_Open_Shrublands.5km.Percentage.ABRAID.tif"
    ),
    "evi_mean"="covariates/EVI_Fixed_Mean_5km_Mean_DEFLATE.ABRAID_Extent.Gapfilled.tif",
    "evi_sd"="covariates/EVI_Fixed_SD_5km_Mean_DEFLATE.ABRAID_Extent.Gapfilled.tif",
    "gecon"="covariates/gecon.tif",
    "lst_day_mean"="covariates/LST_Day_Mean_5km_Mean_DEFLATE.ABRAID_Extent.Gapfilled.tif",
    "lst_day_sd"="covariates/LST_Day_SD_5km_Mean_DEFLATE.ABRAID_Extent.Gapfilled.tif",
    "lst_night_mean"="covariates/LST_Night_Mean_5km_Mean_DEFLATE.ABRAID_Extent.Gapfilled.tif",
    "lst_night_sd"="covariates/LST_Night_SD_5km_Mean_DEFLATE.ABRAID_Extent.Gapfilled.tif",
    "prec57mn"="covariates/prec57mn.tif",
    "prec57mx"="covariates/prec57mx.tif",
    "tcb_mean"="covariates/TCB_Mean_5km_Mean_DEFLATE.ABRAID_Extent.Gapfilled.tif",
    "tcw_mean"="covariates/TCW_Mean_5km_Mean_DEFLATE.ABRAID_Extent.Gapfilled.tif",
    "tempsuit"="covariates/tempsuit.tif",
    "upr_p"="covariates/upr_p.tif",
    "upr_u"="covariates/upr_u.tif",
    "wd0107mn"="covariates/wd0107mn.tif",
    "wd0107mx"="covariates/wd0107mx.tif",
    "wd0114a0"="covariates/wd0114a0.tif"
  )
  all_discrete <- list(
    "access"=FALSE,
    "c10"=FALSE,
    "c6"=FALSE,
    "c7"=FALSE,
    "evi_mean"=FALSE,
    "evi_sd"=FALSE,
    "gecon"=FALSE,
    "lst_day_mean"=FALSE,
    "lst_day_sd"=FALSE,
    "lst_night_mean"=FALSE,
    "lst_night_sd"=FALSE,
    "prec57mn"=FALSE,
    "prec57mx"=FALSE,
    "tcb_mean"=FALSE,
    "tcw_mean"=FALSE,
    "tempsuit"=FALSE,
    "upr_p"=TRUE,
    "upr_u"=TRUE,
    "wd0107mn"=FALSE,
    "wd0107mx"=FALSE,
    "wd0114a0"=FALSE
  )

  covariate_diseases <- list(
    "cchf"=c("c10", "c6", "c7", "evi_mean", "evi_sd", "lst_day_mean", "lst_day_sd", "lst_night_mean", "lst_night_sd", "tcb_mean", "tcw_mean"),
    "chik"=c(),
    "deng"=c("access", "gecon", "prec57mn", "prec57mx", "tempsuit", "upr_p", "upr_u", "wd0114a0"),
    "hat"=c(),
    "melio"=c(),
    "nwcl"=c("prec57mn", "prec57mx", "upr_p", "wd0107mn", "wd0107mx"),
    "nwvl"=c(),
    "owcl"=c(),
    "owvl"=c(),
    "scrub"=c()
  )
  
  covariate_path <- all_covs[covariate_diseases[[disease]]]
  discrete <- all_discrete[covariate_diseases[[disease]]]
  
  # Functions to assist in the loading of raster data. 
  # This works around the truncation of crs metadata in writen geotiffs.
  abraidCRS <- crs("+init=epsg:4326")
  abraidStack <- function(paths) {
    s <- stack(paths)
    crs(s) <- abraidCRS
    extent(s) <- extent(-180, 180, -60, 85)
    return (s)
  }
  abraidRaster <- function(path) {
    r <- raster(path)
    crs(r) <- abraidCRS
    extent(r) <- extent(-180, 180, -60, 85)
    return (r)
  }
  
  # ~~~~~~~~
  # load data
  
  # occurrence data
  occurrence <- read.csv(occurrence_path, stringsAsFactors = FALSE)
  occurrence <- occurrence2SPDF(occurrence, crs=abraidCRS)
  
  # load the definitve extent raster
  extent <- abraidRaster(extent_path)
  
  # load the admin rasters as a stack
  admin <- abraidStack(admin_path)
  
  if (!use_temporal_covariates) {
    # For our data sets selectLatestCovariates is simple enough (will pick 2012 layer),
    # for more complex covariate sets a better approach may be nessesary
    # load_stack = noop, we just want the strings
    covariate_path <- selectLatestCovariates(covariate_path, load_stack=function(x) { return (x) })
  }
  
  # get the required number of cpus
  nboot <- 64

  # start the cluster
  sfInit(parallel = TRUE, cpus = nboot)
  
  # load seegSDM and dependencies on every cluster
  sfLibrary(seegSDM)

  cat('\nseegSDM loaded on cluster\n\n')

  # prepare absence data
  if (mode == 'bhatt') {
    sub <- function(i, pars) {
      # get the $i^{th}$ row of pars 
      pars[i, ]
    }
    
    # set up range of parameters for use in `extractBhatt`
    # use a similar, but reduced set to that used in Bhatt et al. (2013)
    # for dengue.
    
    # number of pseudo-absences per occurrence
    na <- c(1, 4, 8, 12)
    
    # number of pseudo-presences per occurrence (none)
    np <- c(0, 0, 0, 0)
    
    # maximum distance from occurrence data
    mu <- c(10, 20, 30, 40)
    
    # get all combinations of these
    pars <- expand.grid(na = na,
                        np = np,
                        mu = mu)
    
    # convert this into a list
    par_list <- lapply(1:nrow(pars),
                       sub,
                       pars)
    
    # generate pseudo-data in parallel
    data_list <- sfLapply(par_list,
                          abraidBhatt,
                          occurrence = occurrence,
                          covariates = covariate_path,
                          consensus = extent,
                          admin = admin, 
                          factor = discrete,
                          admin_mode = admin_extract_mode,
                          load_stack = abraidStack)
    cat('extractBhatt done\n\n')
  } else if (substr(mode, 1, 4) == "bias") {
    # Load bias data
    supplementary_occurrence <- read.csv(supplementary_occurrence_path, stringsAsFactors = FALSE)
    supplementary_occurrence <- occurrence2SPDF(supplementary_occurrence, crs=abraidCRS)
    
    presence <- occurrence
    presence <- occurrence2SPDF(cbind(PA=1, presence@data), crs=abraidCRS)
    absence <- supplementary_occurrence
    absence <- occurrence2SPDF(cbind(PA=0, absence@data[, 1:2], Weight=1, absence@data[, 3:6]), crs=abraidCRS)
    
    disease_filters <- list(
      'bacteria'=c(34, 35, 43, 48, 52, 64, 71, 92, 132, 135, 187, 192, 193, 194, 199, 201, 211, 245, 266, 271, 284, 304, 315, 325, 340, 343, 351, 361, 362, 364, 367, 377, 34, 35, 43, 48, 52, 64, 71, 92, 132, 135, 187, 192, 193, 194, 199, 201, 211, 245, 266, 271, 284, 304, 315, 325, 340, 343, 351, 361, 362, 364, 367, 377, 34, 35, 43, 48, 52, 64, 71, 92, 132, 135, 187, 192, 193, 194, 199, 201, 211, 245, 266, 271, 284, 304, 315, 325, 340, 343, 351, 361, 362, 364, 367, 377, 34, 35, 43, 48, 52, 64, 71, 92, 132, 135, 187, 192, 193, 194, 199, 201, 211, 245, 266, 271, 284, 304, 315, 325, 340, 343, 351, 361, 362, 364, 367, 377),
      'fungus'=c(72, 80, 152, 254, 331),
      'parasite'=c(22, 26, 41, 81, 84, 93, 96, 109, 118, 119, 131, 133, 157, 171, 215, 240, 255, 258, 341, 349, 350, 353, 356, 359, 360),
      'prion'=c(78),
      'virus'=c(4, 42, 60, 74, 79, 85, 87, 97, 98, 141, 142, 143, 144, 145, 149, 154, 159, 164, 173, 186, 208, 220, 222, 234, 277, 278, 286, 290, 302, 305, 307, 313, 332, 374, 386, 387, 391, 393)
    )
    if (filter_bias) {
      # Filter by disease type
      if (as.character(disease_type) %in% names(disease_filters)) {
        absence <- absence[absence$Disease %in% disease_filters[[disease_type]], ]
        cat('filtered bias to disease subset\n\n')
      }
    } else {
      if (as.character(disease_type) %in% names(disease_filters)) {
        absence <- absence[absence$Disease %in% unlist(disease_filters, use.names=FALSE), ]
        cat('filtered bias to all classified diseases \n\n')
      }
    }
    
    if (crop_bias) {
      # Filter to consensus
      absence_consensus <- extractBatch(absence, list("consensus"=extent), list("consensus"=TRUE), admin, admin_mode="latlong", load_stack=abraidStack)
      absence_cropped <- !is.na(absence_consensus$consensus) & (absence_consensus$consensus == 100 | absence_consensus$consensus == 50)
      absence <- absence[absence_cropped, ]
      cat('filtered bias to extent\n\n')
    }

    all <- rbind(presence, absence)
    
    # create batches
    batches <- replicate(nboot, subsample(all@data, nrow(all), replace=TRUE), simplify=FALSE)
    batches <- lapply(batches, occurrence2SPDF, crs=abraidCRS)
    cat('batches ready for extract\n\n')

    # Do extractions
    data_list <- sfLapply(batches,
                          extractBatch,
                          covariates = covariate_path,
                          admin = admin, 
                          factor = discrete,
                          admin_mode = admin_extract_mode)
  } else if (mode == "uniform") {
    presence <- occurrence
    presence <- occurrence2SPDF(cbind(PA=1, presence@data), crs=abraidCRS)
    selection_mask <- calc(extent, function (cells) {
      return (ifelse(cells %in% c(100,50), 1, 0))
    })
    absence <- bgSample(selection_mask, n=nrow(presence), prob=crop_bias, replace=TRUE, spatial=FALSE)
    absence <- xy2AbraidSPDF(absence, abraidCRS, 0, 1, sample(presence$Date, nrow(presence), replace=TRUE))
    all <- rbind(presence, absence)
    cat('random bias generated\n\n')
  } else {
    exit(1)
  }
  
  cat('extraction done\n\n')

  if (use_weights) {
    # balance weights
    data_list <- sfLapply(data_list, balanceWeights)
    cat('balance done\n\n')
    
    # run BRT submodels in parallel
    model_list <- sfLapply(data_list,
                           runBRT,
                           wt = 'Weight',
                           gbm.x = names(covariate_path),
                           gbm.y = 'PA',
                           pred.raster = selectLatestCovariates(covariate_path, load_stack=abraidStack),
                           gbm.coords = c('Longitude', 'Latitude'),
                           verbose = TRUE)
    cat('model fitting done (with weights)\n\n')
  } else {
    # balance weights
    data_list <- lapply(data_list, function (extracted_batch) {
      return (extracted_batch[, names(extracted_batch) != 'Weight', drop=FALSE])
    })
    
    # run BRT submodels in parallel
    model_list <- sfLapply(data_list,
                           runBRT,
                           gbm.x = names(covariate_path),
                           gbm.y = 'PA',
                           pred.raster = selectLatestCovariates(covariate_path, load_stack=abraidStack),
                           gbm.coords = c('Longitude', 'Latitude'),
                           verbose = TRUE)
    cat('model fitting done (without weights)\n\n')
  }

  # get cross-validation statistics in parallel
  stat_lis <- sfLapply(model_list,
                       getStats)
  
  cat('statistics extracted\n\n')

  # combine and output results
  
  # make a results directory
  dir.create(paste0(name, '/results'))
  
  # cross-validation statistics (with pairwise-weighted distance sampling)
  stats <- do.call("rbind", stat_lis)
  
  # keep only the relevant statistics
  stats <- stats[, c('auc', 'sens', 'spec', 'pcc', 'kappa',
                     'auc_sd', 'sens_sd', 'spec_sd', 'pcc_sd', 'kappa_sd')]
  
  # write stats to disk
  write.csv(stats,
            paste0(name,'/results/statistics.csv'),
            na = "",
            row.names = FALSE)
  
  # relative influence statistics
  relinf <- getRelInf(model_list)
  
  # append the names to the results
  relinf <- cbind(name = rownames(relinf), relinf)
  
  # output this file
  write.csv(relinf,
            paste0(name,'/results/relative_influence.csv'),
            na = "",
            row.names = FALSE)
  
  # marginal effect curves
  effects <- getEffectPlots(model_list)
  
  # convert the effect curve information into the required format
  
  # keep only the first four columns of each dataframe
  effects <- lapply(effects,
                    function (x) {
                      x <- x[, 1:4]
                      names(x) <- c('x',
                                    'mean',
                                    'lower',
                                    'upper')
                      return(x)
                    })
  
  # paste the name of the covariate in as extra columns
  for(i in 1:length(effects)) {
    
    # get number of evaluation points
    n <- nrow(effects[[i]])
    
    # append name to effect curve
    effects[[i]] <- cbind(name = rep(names(effects)[i], n),
                          effects[[i]])
  }
  
  # combine these into a single dataframe
  effects <- do.call(rbind, effects)
  
  # clean up the row names
  rownames(effects) <- NULL
  
  # save the results
  write.csv(effects,
            paste0(name,'/results/effect_curves.csv'),
            na = "",
            row.names = FALSE)
  
  # get summarized prediction raster layers
  
  # lapply to extract the predictions into a list
  preds <- lapply(model_list,
                  function(x) {x$pred})
  
  # coerce the list into a RasterStack
  preds <- stack(preds)
  
  # summarize predictions
  preds <- combinePreds(preds, parallel=TRUE)
  
  # stop the cluster
  sfStop()
  
  # get the width of the 95% confidence envelope as a metric of uncertainty
  uncertainty <- preds[[4]] - preds[[3]]
  
  # save the mean predicitons and uncerrtainty as rasters
  writeRaster(preds[[1]],
              paste0(name,'/results/mean_prediction'),
              format = 'GTiff',
              NAflag = -9999,
              options = c("COMPRESS=DEFLATE",
                          "ZLEVEL=9"),
              overwrite = TRUE)
  
  writeRaster(uncertainty,
              paste0(name,'/results/prediction_uncertainty'),
              format = 'GTiff',
              NAflag = -9999,
              options = c("COMPRESS=DEFLATE",
                          "ZLEVEL=9"),
              overwrite = TRUE)
  
  ###crop/mask
  water <- abraidRaster(water_mask)
  mean <- preds[[1]]
  mean <- mask(mean, extent, maskvalue=-100, updatevalue=9999)
  uncertainty <- mask(uncertainty, extent, maskvalue=-100, updatevalue=9999)
  mean <- mask(mean, water, inverse=TRUE)
  uncertainty <- mask(uncertainty, water, inverse=TRUE)
  
  writeRaster(mean,
              paste0(name,'/results/mean_prediction_masked'),
              format = 'GTiff',
              NAflag = -9999,
              options = c("COMPRESS=DEFLATE",
                          "ZLEVEL=9"),
              overwrite = TRUE)
  
  writeRaster(uncertainty,
              paste0(name,'/results/prediction_uncertainty_masked'),
              format = 'GTiff',
              NAflag = -9999,
              options = c("COMPRESS=DEFLATE",
                          "ZLEVEL=9"),
              overwrite = TRUE)
  
  cols <- colorRampPalette(c('#91ab84', '#c3d4bb', '#ffffcb', '#cf93ba', '#a44883'))
  
  png(paste0(name,'/results/mean_prediction_masked.png'),
      width = 1656,
      height = 667)
  
  par(mfrow=c(1,1), 
      mar = c(0, 0, 0, 0),
      oma = rep(0, 4))
  plot(mean,
       zlim = c(9999, 9999),
       col = colorRampPalette(c('#eaeaea', '#eaeaea'))(2),
       axes = FALSE,
       box = FALSE,
       legend=FALSE)
  plot(mean,
       zlim = c(0, 1),
       col = cols(1000),
       axes = FALSE,
       box = FALSE,
       legend=FALSE,
       add=TRUE)
  
  dev.off()
  
  png(paste0(name,'/results/prediction_uncertainty_masked.png'),
      width = 1656,
      height = 667)
  
  par(mfrow=c(1,1), 
      mar = c(0, 0, 0, 0),
      oma = rep(0, 4))
  plot(uncertainty,
       zlim = c(9999, 9999),
       col = colorRampPalette(c('#eaeaea', '#eaeaea'))(2),
       axes = FALSE,
       box = FALSE,
       legend=FALSE)
  plot(uncertainty,
       zlim = c(0, 1),
       col = cols(1000),
       axes = FALSE,
       box = FALSE,
       legend=FALSE,
       add=TRUE)
  
  dev.off()
  
  png(paste0(name,'/results/effects.png'),
      width = 2000,
      height = 2500,
      pointsize = 30)
  
  par(mfrow = n2mfrow(length(covariate_path)))
  
  getEffectPlots(model_list, plot = TRUE)
  
  dev.off()
  
  write.csv(presence,
            paste0(name,'/results/presence.csv'),
            na = "",
            row.names = FALSE)
  
  write.csv(absence,
            paste0(name,'/results/absence.csv'),
            na = "",
            row.names = FALSE)
  
  # return an exit code of 0, as in the ABRAID-MP code
  return (0)
}