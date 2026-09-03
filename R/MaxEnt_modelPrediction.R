
# Prediction Functions ---------------------------------------------

# Maxnet solution

#' Generate MaxEnt Prediction Using Maxnet
#' @description Projects a fitted maxnet model onto a set of environmental raster layers
#' @param env.vars SpatRaster of environmental predictor layers to predict onto
#' @param SDM.obj ENMeval object containing model results
#' @param model Character string, name of the selected model within SDM.obj (as used to index `eval.models(SDM.obj)`)
#' @param taxon.name Character string, name applied to the output raster layer
#' @param clamp.env Logical, whether to clamp predictions to the range of the model's training data (default TRUE)
#' @return SpatRaster of cloglog suitability predictions, named according to taxon.name
#' @details Extracts the selected model from SDM.obj and projects it onto env.vars using `ENMeval:::enm.maxnet@predict()` with cloglog output
#' @examples
#' \dontrun{
#' predict_maxnet(env.vars, SDM.obj, model = "fc.LQ_rm.1", taxon.name = "species_a")
#' }

predict_maxnet <- function(env.vars, SDM.obj, model, taxon.name, clamp.env = T) {

  if (!requireNamespace("maxnet", quietly = TRUE)) {
    stop("Package 'maxnet' is required for predict_maxnet() but is not installed.")
  }
  # maxnet is a Suggested dependency, so its namespace must be loaded before ENMeval:::enm.maxnet@predict() dispatches
  # https://github.com/jamiemkass/ENMeval/issues/117

  best.model <- ENMeval::eval.models(SDM.obj)[[model]]
  # Extract model for best model

  message("Generating prediction: ", taxon.name)

  future.pred <- ENMeval:::enm.maxnet@predict(best.model, env.vars, list(pred.type ="cloglog", doClamp = clamp.env))
  # Project model onto future climate, allow choice of clamping to be passed up the stack.

  names(future.pred) <- taxon.name
  # Change name from model type to name of taxon

  return(future.pred)

  ## SUPERSEDED
  # if (clamp.env = T) {
  #
  #   pts <- dplyr::bind_rows(SDM.obj@occs[1:2], SDM.obj@bg[1:2])
  #   v <- terra::extract(env.vars, pts, method = 'simple', na.rm = TRUE)
  #   # Extract occurrence and background points to clamp the model
  #   # Model clamping means that environmental variation is clamped to those the model was trained in
  #   # Stops the model extrapolating beyond it's training data. (considered bad practice)
  #
  #   env.vars <- clamp.vars(env.vars, v, left = NULL, right = NULL, categoricals = NULL)
  #   # Clamp future scenario variables by occurrence locations
  #   # derive future scenario predictions using optimal model
  #
  # }
  # ## Is this needed now that clamp option is in following function?
  # Manual method for clamping variables.
  # Made obsolete with doClamp in following function

}

## Maxent.jar solution

#' Generate MaxEnt Prediction Using Maxent.jar
#' @description Projects a fitted maxent.jar model onto a set of environmental raster layers
#' @param env.vars SpatRaster of environmental predictor layers to predict onto
#' @param SDM.obj ENMeval object containing model results
#' @param model Character string, name of the selected model within SDM.obj (as used to index `eval.models(SDM.obj)`)
#' @param taxon.name Character string, name applied to the output raster layer
#' @return SpatRaster of cloglog suitability predictions, named according to taxon.name
#' @details Extracts the selected model from SDM.obj and projects it onto env.vars using `dismo::predict()` with cloglog output format
#' @examples
#' \dontrun{
#' predict_maxentJar(env.vars, SDM.obj, model = "fc.LQ_rm.1", taxon.name = "species_a")
#' }

predict_maxentJar <- function(env.vars, SDM.obj, model, taxon.name) {

  best.model <- ENMeval::eval.models(SDM.obj)[[model]]
  # Extract model for best model

  message("Generating prediction: ", taxon.name)

  future.pred <- dismo::predict(best.model, env.vars, args = c("outputformat=cloglog"))
  # Should be able to getaway with not clamping in maxent.jar as default "should" be to do it
  # Would be nice to confirm 100%, but does appear to be the default from some light googling

  names(future.pred) <- taxon.name
  # Change name from model type to name of taxon

  return(future.pred)

}

#' Calculate Overlap Between Historic and Future Suitability (deprecated)
#' @description Deprecated. Combines a historic thresholded suitability raster with a future thresholded ensemble raster to visualise loss/retention of habitat.
#' @param hist.pred.rast.thresh List of historic thresholded suitability SpatRasters, indexed by species
#' @param future.thresh.ensemble Nested list of future thresholded ensemble SpatRasters, indexed by species, time, then scenario
#' @param species.idx Numeric index identifying the species within `hist.pred.rast.thresh` and `future.thresh.ensemble`
#' @param time.idx Numeric index identifying the time period within `future.thresh.ensemble`
#' @param scenario.idx Numeric index identifying the scenario within `future.thresh.ensemble`

overlap_ensemble_dist.dep <- function(hist.pred.rast.thresh, future.thresh.ensemble, species.idx, time.idx, scenario.idx) {

  # plot(hist.pred.rast.thresh[[species.idx]])
  # plot(future.thresh.ensemble[[species.idx]][[time.idx]][[scenario.idx]])
  # Test plots

  diff.rast <- hist.pred.rast.thresh[[species.idx]] - future.thresh.ensemble[[species.idx]][[time.idx]][[scenario.idx]]
  # plot(diff.rast)
  #  Calculate difference between current and future raster of a particular index

  loss.rast <- terra::ifel(diff.rast == 0, NA, diff.rast * -1)
  # plot(loss.rast)
  # Anywhere that wasn't historically considered habitat, make NA, otherwise, invert the raster.

  combined.rast <- terra::ifel(is.na(loss.rast) & hist.pred.rast.thresh[[species.idx]] == 1, 1, hist.pred.rast.thresh[[species.idx]] + loss.rast)
  # plot(combined.rast)
  # Overlay original habitat to crop out
  # Logic is something like: If no change from historic, rast == 1 based on historic suitability,
  # If anything else, add or subtract the loss

  # test.rast <- ifel(alps.thresh[[11]] == 0, NA, 1)
  # test.rast <- mask(alps.future.thresh.ensemble[[11]][[3]][[2]], test.rast)
  # plot(test.rast)
  # Inital method was to do a simple mask of previous habitat, but wanted to allow areas of increase to appear on map

  return(combined.rast)
  # Returns two object list, a rast and a frequency table
}
# Function for generating a layer of future suitability, while considering historic distribution
# Could be accomplished with a simple mask or crop, but this fills in gaps within the distribution as well.
# Also is designed for ease of looping over and subsetting.
