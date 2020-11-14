#' Functions for aggregrating SpatialCIMIS data
#'
#' @param etoTable A table containing columns that contain at least the following columns:
#' \code{"abs_path"}, contining the path to the daily/weekly/monthy raster layers
#' \code{yearCounter}, a column containing the indexes that you are going to aggregrate _over_
#' (eg, a column of water years if you want to aggregrate by water year).
#' @param yearCounter The name of the column in \code{prismTable} containing your aggregation index.
#' @param areaMask Raster* object or a Spatial* object containing an area mask (if needed).
#' @param fun The aggregration function used to collapse the different raster layers into one
#' summary layer (eg, \code{sum} or \code{mean}).
#' @return A matrix containing statistics run on all of the values within the summary layer.
#' @examples
#' No examples for now
#'
source("../R_orig/utilities.R")

pkgTest("rgdal", "raster")

tagAnomalous <- function(input.table, column.var, criteria) {
  output.tag <- rep(FALSE, nrow(input.table))
  # TODO: There's probably a more effecient way of doing this
    for (item in column.var) {
    output.tag[input.table[item] != criteria(input.table[item])] <- TRUE
  }
  return(output.tag)
}

