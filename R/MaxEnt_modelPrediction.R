
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

# Ensemble Futures ---------------------------------------------------------

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


#' Create Ensemble Predictions Across Future Scenarios
#'
#' @description Averages multiple model prediction rasters into a single ensemble raster for each combination of species, time period, and SSP scenario
#' @param rast.list Nested list of prediction rasters structured as `rast.list[[species]][[time.period]][[ssp]]`, where each leaf element is a SpatRaster (or a list coercible to one) containing one layer per model to be ensembled
#' @return A nested list with the same structure as `rast.list`, where each leaf element has been replaced by a single-layer SpatRaster representing the mean ensemble of its inputs
#' @details For each species/time period/SSP combination, layers are summed and divided by the number of layers to produce a simple mean ensemble. The output layer is renamed to `"{species}_{ssp}_{period}_ensemble"`.
#' @note Deprecated in this project — ensembling is now performed at the environment level before model generation for efficiency, rather than post-hoc at the model prediction level. Retained here for reference.
#' @examples
#' \dontrun{
#' # Requires a nested list of prediction rasters (rast.list)
#' ensemble_modelFutures(rast.list = future.pred.list)
#' }

ensemble_modelFutures <- function(rast.list){

  species.list <- names(rast.list)
  time.periods <- names(rast.list[[1]])
  ssps <- names(rast.list[[1]][[1]])
  # Extract names of tiers of prediction rasters
  # Currently fixed to standard format for simplicity

  for(species in species.list) {

    for(period in time.periods) {

      for(ssp in ssps) {

        message("Creating ensembles of ", species, ": ", period, ": ", ssp)

        ## Average each model in output

        rast.subset <- terra::rast(rast.list[[species]][[period]][[ssp]])

        nlyr(rast.subset)

        for (i in seq(nlyr(rast.subset))) {
          if (i == 1) {
            ensemble.subset <- rast.subset[[i]]
          } else {
            ensemble.subset <- ensemble.subset + rast.subset[[i]]
          }
        }
        ensemble <- ensemble.subset / nlyr(rast.subset)
        # Create mean average ensemble of all provided layers

        names(ensemble) <- paste(species, ssp, period, "ensemble", sep = "_")
        # Change layer name to desired string

        plot(ensemble)

        rast.list[[species]][[period]][[ssp]] <- ensemble
        # Overwrite individual models in list structure.

      }
    }
  }

  return(rast.list)

}
# Function for ensembling post future if models were generated at the model level.
# Depreciated in this project as I reran models by ensembling the environments beforehand for efficiency

# alps.future.ensemble <- ensemble_futures(rast.list = alps.future)
# alps.future.thresh.ensemble <- ensemble_futures(rast.list = alps.future.thresh)
# Combine each future raster stack into ensembles

#' Threshold Future Ensemble Prediction Rasters
#'
#' @description Applies a threshold to future ensemble prediction rasters and writes the thresholded rasters to disk, using model selection choices from the validation pipeline
#' @param future.pred.ensemble Nested list of future ensemble prediction rasters structured as `future.pred.ensemble[[species]][[time.period]][[scenario]]`, one species per element of `sdm.outputs`/`modelSelection`
#' @param data.dir Character string specifying the base data directory containing model outputs and `modelSelection.RDS` (default `"data/maxent/"`)
#' @param thresh.method Character string specifying threshold method passed to `write_model_threshold()`: "maxsss" (default), "or.10.threshold", or "or.min.threshold"
#' @return No return value, called for side effects. Writes thresholded prediction rasters to `{data.dir}rastOutputs/thresholdFuturePredictionsRast/`
#' @details Iterates over `future.pred.ensemble`, matching each species to its best-performing model (per `modelSelection.RDS`) and its `SDM.obj` (per `load_results()`), then loops over each time period and scenario within that species' ensemble list, calling `write_model_threshold()` to threshold and write the ensemble raster for each combination.
#' @examples
#' \dontrun{
#' # Requires a nested list of future ensemble prediction rasters (future.pred.ensemble)
#' threshold_ensembleFutures(future.pred.ensemble = future.ensemble.list, data.dir = "data/maxent/")
#' }

threshold_ensembleFutures <- function(future.pred.ensemble, data.dir = "data/maxent/", thresh.method = "maxsss") {

  sdm.outputs <- load_results(data.dir)
  # Load outputs into a list

  modelSelection <- readRDS(paste0(data.dir, "modelSelection.RDS"))
  # Load model selection choices from validation pipeline

  ## TODO Loop over the raster, not the outputs
  # attempting to make this work with rapply rather than a loop

  for (i in seq(length(future.pred.ensemble))) {

    SDM.obj <- sdm.outputs[[i]]
    # Select model iteration

    taxon.name <- SDM.obj@taxon.name
    model.name <- "ensemble"
    # Extract details for filename

    model <- modelSelection |>
      dplyr::filter(name == taxon.name) |>
      dplyr::select(modelSelection) |>
      as.character()
    # Selects best model from model selection table

    for (time.period in names(future.pred.ensemble[[i]])) {

      for (scenario in names(future.pred.ensemble[[i]][[time.period]])) {

        future.pred <- future.pred.ensemble[[i]][[time.period]][[scenario]]
        # Select ensemble raster for this species/time period/scenario combination

        message("Thresholding future ensemble model raster: ", taxon.name, ": ", time.period, ": ", scenario)

        write_model_threshold(output.dir = paste0(data.dir, "rastOutputs/thresholdFuturePredictionsRast/"),
                              SDM.obj = SDM.obj,
                              model = model,
                              thresh.method = thresh.method,
                              future.pred = future.pred,
                              filename = paste(taxon.name, time.period, scenario, model.name, sep = "_"))
        # Write prediction raster to output folder
        # Appears to work with both maxnet and maxent.jar
        # TODO MaxSSS Thresholding not working with Maxnet

      }
    }
  }
}
## TODO move into maxent toolkit?
# TODO Need to refactor this so it loops over future ensembles
# Ideally ensembling might occurr in the main maxent pipeline to reduce steps, but just need a
# quick fix for now.
# Can use top level names as taxon reference.

# test <- rapply(
#   alps.future.ensemble,
#   f = function(r) terra::crop(r, ibra.vect),
#   classes = "SpatRaster",
#   how = "replace"
# )
# Use this bit to apply a single threshold over the list element. rapply is a recursive apply

# Should I be re-thresholding the ensemble model output, rather than doing this?
# Pull the function from maxent toolkit.


#' Threshold Ensemble Prediction Rasters by Model Convergence
#'
#' @description Applies a model convergence threshold to ensemble prediction rasters, converting values to binary (1/0) based on whether they meet or exceed a specified number of models
#' @param rast.ensemble Nested list of ensemble prediction rasters structured as `rast.ensemble[[species]][[time.period]][[scenario]]`
#' @param models.threshold Numeric specifying the minimum number of models that must agree (i.e. raster value >= this threshold) for a cell to be retained as suitable habitat (default `1`). Values are converted to 1 if >= threshold, 0 otherwise
#' @return Nested list of binary (1/0) rasters with the same structure as `rast.ensemble`, where 1 indicates convergence at or above `models.threshold`
#' @details Iterates over species, time periods, and scenarios within the nested raster ensemble structure, applying `ifel(raster >= models.threshold, 1, 0)` to each raster to create binary convergence maps. Printed messages show progress at species and time period levels. Intended for use with ensemble rasters where cell values represent counts of agreeing models (e.g. from `ensemble_modelFutures()` averaged across multiple model runs)
#' @examples
#' \dontrun{
#' # Requires a nested list of ensemble prediction rasters
#' binary.ensemble <- converge_ensemble(rast.ensemble = future.ensemble.list, models.threshold = 3)
#' }
#' @export

converge_ensemble <- function(rast.ensemble, models.threshold = 1){

  message("Prepping ensemble thresholds for divraster...")

  for (i in seq_along(names(rast.ensemble))) {

    taxon.name <- names(rast.ensemble)[[i]]
    print(taxon.name)

    for (y in seq_along(names(rast.ensemble[[i]]))) {

      time.period <- names(rast.ensemble[[i]])[[y]]
      print(time.period)

      for (z in seq_along(names(rast.ensemble[[i]][[y]]))) {

        rast.ensemble[[i]][[y]][[z]] <- ifel(rast.ensemble[[i]][[y]][[z]] >= models.threshold, 1, 0)

      }
    }
  }
  # Loop over each of the lists and convert

  return(rast.ensemble)

}
# Function for setting model convergence threshold that you wish to include


# Beta Diversity  ---------------------------------------------------------

#' Reorganise Nested Rasters by Time Period and Scenario
#'
#' @description Restructures a nested list of prediction rasters from species-first ordering to time-period-first ordering, stacking all species rasters together at each time period/scenario combination
#' @param rast Nested list of rasters structured as `rast[[species]][[time.period]][[scenario]]`, where each element is a single-layer raster
#' @return Nested list of rasters restructured to `rast[[time.period]][[scenario]][[species]]` where each element is a multi-layer raster stack with one layer per species. Layer names correspond to species names from the input
#' @details Reorganises the nesting hierarchy by creating the new structure in two passes: first, copies individual species rasters into the new `[[time.period]][[scenario]][[species]]` structure; second, stacks all species rasters within each time period/scenario combination into a single multi-layer raster using `c()`. This is useful for computing metrics (e.g. species richness) that require all species rasters to be in a single stack
#' @examples
#' \dontrun{
#' # Requires a nested list of species prediction rasters
#' rast.by.period <- cluster_bySpecies(rast = species.predictions)
#' }
#' @export

cluster_bySpecies <- function(rast) {

  new.list <- list()
  # Create shell of new structure

  for (species in names(rast)) {
    for (timeperiod in names(rast[[species]])) {
      for (ssp in names(rast[[species]][[timeperiod]])) {

        # Create structure if it doesn't exist
        if (!timeperiod %in% names(new.list)) {
          new.list[[timeperiod]] <- list()
        }
        if (!ssp %in% names(new.list[[timeperiod]])) {
          new.list[[timeperiod]][[ssp]] <- list()
        }

        # Stack old raster into new list structure
        new.list[[timeperiod]][[ssp]][[species]] <- rast[[species]][[timeperiod]][[ssp]]

      }
    }
  }
  # Reorder raster stack to have species

  for (timeperiod in names(new.list)) {
    for (ssp in names(new.list[[timeperiod]])) {
      for (i in seq_along(names(new.list[[timeperiod]][[ssp]]))) {

        if (i == 1) {
          new.rast <- new.list[[timeperiod]][[ssp]][[i]]
        } else {
          new.rast <- c(new.rast, new.list[[timeperiod]][[ssp]][[i]])
        }
        # Creating and stacking species into a single raster

      }

      new.list[[timeperiod]][[ssp]] <- new.rast
      names(new.list[[timeperiod]][[ssp]]) <- names(rast)
      # Replace species rasts with new stacked rasts

    }
  }
  # Stack the reorganised rasters, reducing one level of the lists.

  return(new.list)

}

#' Calculate Species Richness from Stacked Rasters
#'
#' @description Sums all layers in a multi-layer raster stack to produce a species richness map
#' @param rast.stack Multi-layer raster (SpatRaster) with one layer per species, where cell values represent presence/absence (binary) or probability of presence
#' @return Single-layer SpatRaster with layer name "Species Richness", where each cell value is the sum of all input layers (count of species present or predicted suitable at each cell)
#' @details Iterates over all layers in the raster stack and accumulates their sum into a single output raster. Assumes cell values are numeric (e.g. binary 0/1 or continuous suitability scores). Intended for use with raster stacks produced by `cluster_bySpecies()`, which organises species prediction rasters into a single stack per time period/scenario combination
#' @examples
#' \dontrun{
#' # Requires a multi-layer raster stack with one species per layer
#' richness.rast <- calc_species_richness(rast.stack = species.raster.stack)
#' }
#' @export

calc_species_richness <- function(rast.stack) {

  for (i in seq_along(names(rast.stack))) {
    if (i == 1) {
      richness <- rast.stack[[i]]
    } else {
      richness <- richness + rast.stack[[i]]
    }
  }
  # Add's all rasters in a stack together

  names(richness) <- "Species Richness"

  return(richness)

}
