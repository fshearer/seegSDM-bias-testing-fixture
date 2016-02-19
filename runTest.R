library(raster)
library(sp)
## This is fork of runABRAID from SEEG-Oxford/seegSDM @ 0.1-8
runTest <- function (mode, 
                     disease_type, # used with "filter_bias"
                     occurrence_path,
                     extent_path,
                     supplementary_occurrence_path,
                     admin_path,
                     covariate_path,
                     discrete,
                     water_mask,
                     admin_extract_mode="random",
                     crop_bias=TRUE, # used with mode "bias"
                     filter_bias=TRUE, # used with mode "bias"
                     use_weights=TRUE,
                     use_temporal_covariates=TRUE) {

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
    
    if (filter_bias) {
      # Filter by disease type
      disease_filters <- list(
        'viral'=c(60, 79, 97, 173, 212, 234, 302, 386, 391, 393)
      )
      if (as.character(disease_type) %in% names(disease_filters)) {
        absence <- absence[absence$Disease %in% disease_filters[[disease_type]], ]
        cat('filtered bias to disease subset\n\n')
      }
    } else {
      if (as.character(disease_type) %in% names(disease_filters)) {
        absence <- absence[absence$Disease %in% unlist(disease_filters, use.names=FALSE), ]
        cat('filtered bias to disease subset\n\n')
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