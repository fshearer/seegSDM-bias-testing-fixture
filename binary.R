library("raster")
library("sp")

diseases <- c("cchf", "chik", "deng", "hat", "melio", "nwcl", "owcl", "scrub")
variants <- c("1A", "1B", "1C", "2C", "2D", "3B", "4A", "5A", "6A")


tidy <- function(r) {
    crs(r) <- crs("+init=epsg:4326")
    extent(r) <- extent(-180, 180, -60, 85)
    return(r)
}

load <- function(path) {
    r <- raster(path)
    r <- tidy(r)
    return (r)
}
for (disease in diseases){
    for(variant  in variants){
        print(paste(disease, variant))
        mean_path <- paste("all/results/", disease, "_", variant,"/mean_prediction.tif", sep="")
        masked_mean_path <- paste("all/results/", disease, "_", variant,"/mean_prediction_masked_v2.tif", sep="")
        pres_path <- paste("all/", disease, "_data/occurrences.csv", sep="")
        out_path <- paste("all/results/", disease, "_", variant,"/binary_prediction.tif", sep="")
        mean <- load(mean_path)
        pres <- read.csv(pres_path, stringsAsFactors = FALSE)
        pres <- pres[pres$Admin == -999, c('Longitude', 'Latitude')]
        pres <- SpatialPointsDataFrame(pres, pres, coords.nrs  = c(1,2), proj4string = crs("+init=epsg:4326"))
        pres <- extract(mean, pres)
        thresh <- quantile(pres, probs=c(0.25))[[1]]
        print(thresh)
        masked_mean <- load(masked_mean_path)
        binary <- masked_mean >= thresh & masked_mean != 9999
        writeRaster(binary, out_path, "GTiff", NAflag=-9999, options=c("COMPRESS=DEFLATE", "ZLEVEL=9"), overwrite=TRUE)
    }
}

