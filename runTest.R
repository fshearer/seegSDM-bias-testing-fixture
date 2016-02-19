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
                     water_mask,
                     verbose = TRUE,
                     max_cpus = 32,
                     parallel_flag = TRUE) {
  
  # Given the locations of: a csv file containing disease occurrence data
  # (`occurrence_path`, a character), a GeoTIFF raster giving the definitive
  # extents of the disease (`extent_path`, a character), a csv file containing 
  # disease occurrence data for other diseases (`supplementary_occurrence_path`,
  # a character), GeoTIFF rasters giving standardised admin units (`admin0_path`,
  # `admin1_path`, `admin2_path`) and GeoTIFF rasters giving the covariates to
  # use (`covariate_path`,  a character vector). Run a predictive model to produce
  # a risk map (and associated outputs) for the disease.
  # The file given by `occurrence_path` must contain the columns 'Longitude',
  # 'Latitude' (giving the coordinates of points), 'Weight' (giving the degree
  # of weighting to assign to each occurrence record), 'Admin' (giving the
  # admin level of the record - e.g. 1, 2 or 3 for polygons or -999 for points),
  # 'GAUL' (the GAUL code corresponding to the admin unit for polygons, or
  # NA for points) and 'Disease' a numeric identifer for the disease of the occurrence.
  # The file given by `supplementary_occurrence_path` must contain the columns 
  # 'Longitude', 'Latitude' (giving the coordinates of points), 'Admin' (giving the
  # admin level of the record - e.g. 1, 2 or 3 for polygons or -999 for points),
  # 'GAUL' (the GAUL code corresponding to the admin unit for polygons, or
  # NA for points) and 'Disease' a numeric identifer for the disease of the occurrence.
  # To treat any covariates as discrete variables, provide a logical vector
  # `discrete` with `TRUE` if the covariate is a discrete variable and `FALSE`
  # otherwise. By default, all covariates are assumed to be continuous.
  # Set the maximum number of CPUs to use with `max_cpus`. At present runABRAID
  # runs 64 bootstrap submodels, so the number of cpus used in the cluster will
  # be set at `min(64, max_cpus)`.
  
  # ~~~~~~~~
  # check inputs are of the correct type and files exist
  abraidCRS <- crs("+init=epsg:4326")
  modes <- c("bhatt", "bias_all", "bias_cropped", "bias_cropped_filtered", "uniform")
  stopifnot(class(mode) == 'character' &&
              is.element(mode, modes))
  
  stopifnot(is.numeric(disease))
  
  stopifnot(class(occurrence_path) == 'character' &&
              file.exists(occurrence_path))
  
  stopifnot(class(extent_path) == 'character' &&
              file.exists(extent_path) && 
              compareCRS(raster(extent_path), abraidCRS))
  
  stopifnot(class(supplementary_occurrence_path) == 'character' &&
              file.exists(supplementary_occurrence_path))
  
  stopifnot(file.exists(admin0_path) && 
              compareCRS(raster(admin0_path), abraidCRS))
  
  stopifnot(file.exists(admin1_path) && 
              compareCRS(raster(admin1_path), abraidCRS))
  
  stopifnot(file.exists(admin2_path) && 
              compareCRS(raster(admin2_path), abraidCRS))
  
  stopifnot(class(unlist(discrete)) == 'logical' &&
              length(discrete == length(covariate_path)))
  
  stopifnot(is.function(load_seegSDM))
  
  stopifnot(is.logical(parallel_flag))
  
  stopifnot(names(discrete) == names(covariate_path))
  
  stopifnot(class(unlist(covariate_path, recursive=TRUE)) == 'character' &&
              all(file.exists(unlist(covariate_path, recursive=TRUE))) &&
              all(sapply(sapply(unlist(covariate_path, recursive=TRUE), raster), compareCRS, abraidCRS)))
  
  # ~~~~~~~~
  # load data
  
  # occurrence data
  occurrence <- read.csv(occurrence_path,
                         stringsAsFactors = FALSE)
  
  # check column names are as expected
  stopifnot(sort(colnames(occurrence)) == sort(c('Admin',
                                                 'Date',     
                                                 'Disease',
                                                 'GAUL',
                                                 'Latitude',
                                                 'Longitude',
                                                 'Weight')))
  
  # convert it to a SpatialPointsDataFrame
  # NOTE: `occurrence` *must* contain columns named 'Latitude' and 'Longitude'
  occurrence <- occurrence2SPDF(occurrence, crs=abraidCRS)
  
  # occurrence data
  supplementary_occurrence <- read.csv(supplementary_occurrence_path,
                                       stringsAsFactors = FALSE)
  
  # check column names are as expected
  stopifnot(sort(colnames(supplementary_occurrence)) == sort(c('Admin',
                                                               'Date',    
                                                               'Disease',
                                                               'GAUL',
                                                               'Latitude',
                                                               'Longitude')))
  # Functions to assist in the loading of raster data. 
  # This works around the truncation of crs metadata in writen geotiffs.
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
  
  # convert it to a SpatialPointsDataFrame
  # NOTE: `occurrence` *must* contain columns named 'Latitude' and 'Longitude'
  supplementary_occurrence <- occurrence2SPDF(supplementary_occurrence, crs=abraidCRS)
  
  # load the definitve extent raster
  extent <- abraidRaster(extent_path)
  
  # load the admin rasters as a stack
  # Note the horrible hack of specifying admin 3 as the provided admin 1.
  # These should be ignored since ABRAID should never contain anything other
  # than levels 0, 1 and 2
  admin <- abraidStack(admin_path)
  
  # get the required number of cpus
  nboot <- 64
  ncpu <- min(nboot,
              max_cpus)
  
  # start the cluster
  sfInit(parallel = parallel_flag,
         cpus = ncpu)
  
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
                         verbose = verbose)
  
  if (verbose) {
    cat('model fitting done\n\n')
  }
  
  # get cross-validation statistics in parallel
  stat_lis <- sfLapply(model_list,
                       getStats)
  
  if (verbose) {
    cat('statistics extracted\n\n')
  }
  
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
  preds <- combinePreds(preds, parallel=parallel_flag)
  
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