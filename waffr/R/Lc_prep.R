#' There are different export options in the  USDA, NASS, CropScape and Cropland Data Layers tool.
#' Users can export subsets of the national data set under different regional masks, such as state
#' boundaries, however all inherit the default `WGS 84 / Lon Lat` CRS, reproduced below as OGC WKT:
#'
#' >     GEOGCS["WGS 84",
#' >         DATUM["WGS_1984",
#' >             SPHEROID["WGS 84",6378137,298.257223563,
#' >                 AUTHORITY["EPSG","7030"]],
#' >             AUTHORITY["EPSG","6326"]],
#' >         PRIMEM["Greenwich",0,
#' >             AUTHORITY["EPSG","8901"]],
#' >         UNIT["degree",0.0174532925199433,
#' >             AUTHORITY["EPSG","9122"]],
#' >         AUTHORITY["EPSG","4326"]]

#TODO: remove/ moce these comments. once the paths are sent to funtions as parameters, remove these declarations as well.
##ca_boundary defines the boundary for the area/(s) of analysis (includes counties)
## ca_boundary <- readOGR("input/MAPS/cb_2014_CA_5m.shp", layer="cb_2014_CA_5m")
ca_boundary.path <- "C:/Users/shivn/Box/wfar/input/SHAPEFILES/cnty24k09/cnty24k09_state_poly_s100.shp"
lc.dir <- "C:/Users/shivn/Box/wfar/input/CDL_CA_WGS84"
lc.crs <- CRS("+init=epsg:4326")



tagAnomalous <- function(input.table, column.var, criteria) {
  output.tag <- rep(FALSE, nrow(input.table))
  # TODO: There's probably a more effecient way of doing this
  for (item in column.var) {
    output.tag[input.table[item] != criteria(input.table[item])] <- TRUE
  }
  return(output.tag)
}

rtTagParam <- function(rt.table) {
  rt.table[c("xmin", "xmax", "ymin", "ymax")] <- t(sapply(rt.table[["abs_path"]],
                                                          function(x) as.vector(extent(raster(x)))))
  rt.table["ncell"] <- sapply(rt.table[["abs_path"]],
                              function(x) ncell(raster(x)))
  rt.table[c("xres", "yres")] <- t(sapply(rt.table[["abs_path"]],
                                          function(x) res(raster(x))))
  rt.table["nlayers"] <- sapply(rt.table[["abs_path"]],
                                function(x) ifelse(nlayers(stack(x))>1,nlayers(stack(x)),NA))

  return(rt.table)
}


################ 1. Collect landcover rasters and inspect ####################

#' First we collect all of the landcover rasters in our input directory and organize
#' them into a table with relevant metadata, e.g. resolution and extent. Rasters are
#' tagged as **anomalous** if their extent or resolution is different from the most
#' common values of the group (criteria `mode`).
#'
#' @param lc.dir: input directory containing landcover raster images - for now see example under VICE Lab/RESEARCH/PROJECTS/wfar/input/CDL_CA_WGS84
prepareAnalysis <- function(lc.dir){
  lc.paths <- list.files(path=lc.dir, pattern=".(tif)$", full.names=T, recursive=TRUE)

  lc.table <- data.frame(
    abs_path = lc.paths,
    source = sapply(strsplit(file_path_sans_ext(lc.paths),"/"),'[[',2),
    product_name = sapply(strsplit(file_path_sans_ext(lc.paths),"/"),'[[',4),
    date = as.Date(sprintf("%s-01-01",regmatches(lc.paths,regexpr("(?<=CA_WGS84\\/CDL_)(.{4})(?=_clip)",lc.paths, perl=TRUE))), format="%F"),
    stringsAsFactors = FALSE
  )

  # HACK: Clean up ugly product name

  lc.table[["product_name"]] <- paste0("CDL_", 2007:2016)

  # Tag (x,y)[min,max,res], ncell, and nlayers if applicable
  lc.table <- rtTagParam(lc.table)
  # Tag anomalous parameters according to metadata parameters specified
  lc.table["anoml"] <- tagAnomalous(lc.table, c("xmin", "xmax", "ymin", "ymax", "ncell", "xres", "yres"), modal)
  rm(lc.paths) #do we need this line since I changed this to a function? I do not know if R uses scope the same way as other languages

  return(lc.table)
}


#' 2. Clean anomalous rasters and re-inspect

#' 2.1 Upscale coarse landcover rasters
#' At this point, we visually inspect anom.indicies and perform manual cleaning
#' depending on what the problem is. For CDL, the resolution of older products is 56,
#' so we must upscale to match the later years.
#' **NOTE:** Here, the criterion for cleaning is `anoml == TRUE`, since there are no other anomalous features, but this probably won't be the same for different data sources!

#' DATA CLEANING
#'anom.table = makeAnomtable(lc.table)

#' Upscale with nearest neighbor interpolation (to preserve categorical variable)
#' to match resolution of most recent CDL (56 meters)
#'
#' `resample_cdl` is a function that resamples each raster that you feed into it to match the dimensions of *the last raster in the landcover file table*. It can be run sequentially
#' on one thread (commented out below), or in parallel on multiple threads using `parallel`.
resample_cdl <- function(lc.table, abs_path, source, product_name) {
  dir.create(paste0("output/cleaned_inputs/", source, "_30m/"), recursive = TRUE, showWarnings = FALSE)
  outpath <- paste0("output/cleaned_inputs/", source, "_30m/", product_name, ".tif")
  resample(raster(abs_path), raster(lc.table[["abs_path"]][nrow(lc.table)]),
           method = "ngb", filename = outpath, format = "GTiff",
           prj = FALSE, progress = "text", datatype = "INT1U", overwrite = TRUE)
}

# TODO: Figure out why apply returns extra row, creating NA directory
# Non parallel version
# apply(lc.table[lc.table["xres"] == 56,], 1,
#       function(x)
#         resample_cdl(lc.table, x["abs_path"], x["source"], x["product_name"]))

# Make only as many clusters as necessary, bound by available cores


  cl <- makeCluster(min((detectCores() - 1), sum(lc.table["anoml"] == TRUE)))
  clusterExport(cl, list("lc.table", "resample_cdl"))
  clusterEvalQ(cl,library(raster))
  parRapply(cl, lc.table[lc.table["anoml"] == TRUE,],
         function(x)
           resample_cdl(lc.table, x["abs_path"], x["source"], x["product_name"]))
  stopCluster(cl)


  ######################


  project_cdl <- function(lc.table, abs_path, source, product_name) {
    dir.create(paste0("output/cleaned_inputs/", source, "_30m/"), recursive = TRUE, showWarnings = FALSE)
    outpath <- paste0("output/cleaned_inputs/", source, "_30m/", product_name, ".tif")
    projectRaster(raster(abs_path), crs = CRS("+init=epsg:4326"),
                  method = "ngb", filename = outpath, format = "GTiff",
                  prj = TRUE, progress = "text", datatype = "INT1U", overwrite = TRUE)
  }

  cl <- makeCluster(min((detectCores() - 1), nrow(lc.table)), outfile = "debug.txt")
  clusterExport(cl, list("lc.table", "project_cdl"))
  clusterEvalQ(cl,library(raster))
  parRapply(cl, lc.table,
            function(x)
              project_cdl(lc.table, x["abs_path"], x["source"], x["product_name"]))
  stopCluster(cl)



  ######### 2.2 Import and re-inspect upscaled rasters###########
  # Load cleaned files into file table
  # Load cleaned files into file table
  # TODO: DRY
  lc.paths.cleaned <- list.files(path="output/cleaned_inputs/CDL_CA_30m", pattern=".(tif)$", full.names=T, recursive=TRUE)

  # TODO: remove hardcoded gsub
  lc.table.cleaned <- data.frame(
    abs_path = lc.paths.cleaned,
    source =  gsub('.{4}$', '', sapply(strsplit(file_path_sans_ext(lc.paths.cleaned),"/"),'[[',3)),
    product_name = sapply(strsplit(file_path_sans_ext(lc.paths.cleaned),"/"),'[[',4),
    date = as.Date(sprintf("%s-01-01",regmatches(lc.paths.cleaned,regexpr("(?<=CA_WGS84_30m\\/CDL_)(.*)(?=.tif)",lc.paths.cleaned, perl=TRUE))), format="%F"),
    stringsAsFactors = FALSE
  )

  # Merge cleaned into file table and redo checks
  # TODO: DRY
  # TODO: Add test to ensure that all of "anoml" column == FALSE
  lc.table[lc.table[["date"]] %in% lc.table.cleaned[["date"]],] <- lc.table.cleaned[lc.table.cleaned[["date"]] %in% lc.table[["date"]],]

  # Tag (x,y)[min,max,res], ncell, and nlayers if applicable
  lc.table <- rtTagParam(lc.table)
  # Tag anomalous parameters according to metadata parameters specified
  lc.table["anoml"] <- tagAnomalous(lc.table, c("xmin", "xmax", "ymin", "ymax", "ncell", "xres", "yres"), modal)

  rm(lc.paths.cleaned, lc.table.cleaned)
