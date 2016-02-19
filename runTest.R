library(raster)
library(sp)
## This is fork of runABRAID from SEEG-Oxford/seegSDM @ 0.1-8
runTest <- function (mode, 
                     disease_type,
                     occurrence_path,
                     extent_path,
                     supplementary_occurrence_path,
                     admin_path,
                     covariate_path,
                     discrete,
                     water_mask) {

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
  
  # occurrence data
  supplementary_occurrence <- read.csv(supplementary_occurrence_path, stringsAsFactors = FALSE)
  supplementary_occurrence <- occurrence2SPDF(supplementary_occurrence, crs=abraidCRS)
  
  # load the definitve extent raster
  extent <- abraidRaster(extent_path)
  
  # load the admin rasters as a stack
  admin <- abraidStack(admin_path)
  
  # get the required number of cpus
  nboot <- 64

  # start the cluster
  sfInit(parallel = TRUE, cpus = nboot)
  
  # load seegSDM and dependencies on every cluster
  sfClusterCall(load_seegSDM)

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
                          load_stack = abraidStack)
    cat('extractBhatt done\n\n')
  } else if (substr(mode, 1, 4) == "bias") {
    presence <- occurrence
    presence <- occurrence2SPDF(cbind(PA=1, presence@data), crs=abraidCRS)
    absence <- supplementary_occurrence
    absence <- occurrence2SPDF(cbind(PA=0, absence@data[, 1:2], Weight=1, absence@data[, 3:6]), crs=abraidCRS)
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
                          admin_mode="random")
  } else {
    exit(1)
  }
  
  cat('extraction done\n\n')

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
  
  cat('model fitting done\n\n')

  # get cross-validation statistics in parallel
  stat_lis <- sfLapply(model_list,
                       getStats)
  
  cat('statistics extracted\n\n')

  # combine and output results
  
  # make a results directory
  dir.create('results')
  
  # cross-validation statistics (with pairwise-weighted distance sampling)
  stats <- do.call("rbind", stat_lis)
  
  # keep only the relevant statistics
  stats <- stats[, c('auc', 'sens', 'spec', 'pcc', 'kappa',
                     'auc_sd', 'sens_sd', 'spec_sd', 'pcc_sd', 'kappa_sd')]
  
  # write stats to disk
  write.csv(stats,
            'results/statistics.csv',
            na = "",
            row.names = FALSE)
  
  # relative influence statistics
  relinf <- getRelInf(model_list)
  
  # append the names to the results
  relinf <- cbind(name = rownames(relinf), relinf)
  
  # output this file
  write.csv(relinf,
            'results/relative_influence.csv',
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
            'results/effect_curves.csv',
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
              'results/mean_prediction',
              format = 'GTiff',
              NAflag = -9999,
              options = c("COMPRESS=DEFLATE",
                          "ZLEVEL=9"),
              overwrite = TRUE)
  
  writeRaster(uncertainty,
              'results/prediction_uncertainty',
              format = 'GTiff',
              NAflag = -9999,
              options = c("COMPRESS=DEFLATE",
                          "ZLEVEL=9"),
              overwrite = TRUE)
  
  ##TODO add crop/mask
  
  # return an exit code of 0, as in the ABRAID-MP code
  return (0)
}