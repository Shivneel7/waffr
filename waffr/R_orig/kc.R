#' Functions for aggregrating NASS CDL data
#'
#' @param lc.Table A table containing columns that contain at least the following columns:
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

pkgTest("readr", "lubridate")

tagAnomalous <- function(input.table, column.var, criteria) {
  output.tag <- rep(FALSE, nrow(input.table))
  # TODO: There's probably a more effecient way of doing this
  for (item in column.var) {
    output.tag[input.table[item] != criteria(input.table[item])] <- TRUE
  }
  return(output.tag)
}

makeAnomtable <- function(parameter.table) {
  return(as.data.frame(sapply(c("xmin", "xmax", "ymin", "ymax", "ncell", "xres", "yres"),
                              function(x) tagAnomalous(parameter.table, x, modal))))
}

tagKcDaily <- function(t, crop) {
  # TODO: Refactor into cases
  # TODO: Use interpolation function to allow for different smoothers
  fKc <- ifelse(t < crop[["t_a.1"]], 0,
  ifelse(t >= crop[["t_a.1"]] & t <= crop[["t_e.1"]],  crop[["1st VALUE"]], 0))

  # Only execute next calc if dual crop parameters exist
  if (!is.na(crop[["t_a.2"]])) {
    fKc2 <- ifelse(t < crop[["t_a.2"]], 0,
          ifelse(t >= crop[["t_a.2"]] & t <= crop[["t_e.2"]],  crop[["2nd VALUE"]], 0))

    # TEST: for overlap between crop 1 and crop 2
    # DRY: This test should occur when CDL.Kc.LUT is imported.
    # TODO: Move to unit tests
    if(any(fKc!=0 & fKc2!=0))
      warning(paste0("oops, you may have double crops with
                     overlapping growing periods for ", crop[["cdl_name"]]))
    # Merge dual-cropping scenarios
    fKc[fKc2!=0] = fKc2[fKc2!=0]
  }
  return(fKc)
}

calculateKcDaily <- function(t, crop) {
  # TODO: Refactor into cases
  # TODO: Use interpolation function to allow for different smoothers
  fKc <- ifelse(t < crop[["t_a.1"]], 0,
  ifelse(t >= crop[["t_a.1"]] & t <= crop[["t_b.1"]],  crop[["Kc AB"]],
  ifelse(t > crop[["t_b.1"]] & t < crop[["t_c.1"]],
         crop[["Kc AB"]] + round(as.numeric(t - crop[["t_b.1"]], units = "days") *
         ((crop[["Kc CD"]] - crop[["Kc AB"]]) / as.numeric((crop[["t_c.1"]] - crop[["t_b.1"]]), units = "days"))),
  ifelse(t >= crop[["t_c.1"]] & t <= crop[["t_d.1"]],  crop[["Kc CD"]],
  ifelse(t > crop[["t_d.1"]] & t < crop[["t_e.1"]],
         crop[["Kc CD"]] + round(as.numeric(t - crop[["t_d.1"]], units = "days") *
         ((crop[["Kc E"]] - crop[["Kc CD"]]) / as.numeric((crop[["t_e.1"]] - crop[["t_d.1"]]), units = "days"))),
  ifelse(t == crop[["t_e.1"]], crop[["Kc E"]], 0))))))

  # Only execute next calc if dual crop parameters exist
  if (!is.na(crop[["t_a.2"]])) {
    fKc2 <- ifelse(t < crop[["t_a.2"]], 0,
    ifelse(t >= crop[["t_a.2"]] & t <= crop[["t_b.2"]],  crop[["2nd Kc AB"]],
    ifelse(t > crop[["t_b.2"]] & t < crop[["t_c.2"]],
           crop[["2nd Kc AB"]] + round(as.numeric(t - crop[["t_b.2"]], units = "days") *
           ((crop[["2nd Kc CD"]] - crop[["2nd Kc AB"]]) / as.numeric((crop[["t_c.2"]] - crop[["t_b.2"]]), units = "days"))),
    ifelse(t >= crop[["t_c.2"]] & t <= crop[["t_d.2"]], crop[["2nd Kc CD"]],
    ifelse(t > crop[["t_d.2"]] & t < crop[["t_e.2"]],
           crop[["2nd Kc CD"]] + round(as.numeric(t - crop[["t_d.2"]], units = "days") *
           ((crop[["2nd Kc E"]] - crop[["2nd Kc CD"]]) / as.numeric((crop[["t_e.2"]] - crop[["t_d.2"]]), units = "days"))),
    ifelse(t == crop[["t_e.2"]], crop[["2nd Kc E"]], 0))))))

    # TEST: for overlap between crop 1 and crop 2
    # DRY: This test should occur when CDL.Kc.LUT is imported.
    # TODO: Move to unit tests
    if(any(fKc!=0 & fKc2!=0))
      warning(paste0("oops, you may have double crops with
                     overlapping growing periods for ", crop[["cdl_name"]]))
    # Merge dual-cropping scenarios
    fKc[fKc2!=0] = fKc2[fKc2!=0]
  }
  return(fKc)
}
