# Utilities ---------------------------------------------------------------

#' Extract Names at Every Level of a Nested List
#' @description Walks down a nested list structure and extracts the names present at each depth level
#' @param rast.list Nested list object (e.g. of SpatRasters) with names at one or more levels
#' @return Named list with `names` (a list of character vectors, one per depth level) and `depth` (total depth of the list)
#' @details Walks down the first element at each level to determine names and depth. Assumes a symmetric list structure where all branches share the same names and depth.
#' @examples
#' \dontrun{
#' get_nested_names(future.env)
#' }
#' @export

get_nested_names <- function(rast.list) {
  depth_names <- list()
  current <- rast.list
  depth <- 1
  # Initialise trackers for current position and depth

  while (is.list(current) && length(current) > 0) {
    depth_names[[paste0("level_", depth)]] <- names(current)
    current <- current[[1]]
    depth <- depth + 1
  }
  # Walks down the first element at each level to extract names
  # Assumes a symmetric list structure where all branches share the same names and depth

  list(
    names = depth_names,
    depth = length(depth_names)
  )
  # Returns names at each depth level and total depth count

}
# Extracts names at every level of a nested list without hardcoding depth
# Could be useful helper at various points

## Env Data Functions ------------------------------------------------------

#' Create Ensemble Averages from Multiple GCM Climate Rasters
#' @description Averages environmental variables across GCMs within each time period and SSP, and writes the resulting ensemble rasters to file
#' @param rast.list Nested list of SpatRasters structured as time.period > ssp > gcm, each containing the same set of environmental variables
#' @param ensemble.dir Character string specifying output directory for ensemble rasters (default "data/ensembleFutureEnv/")
#' @return The original rast.list, unmodified (ensemble rasters are written to file rather than returned)
#' @details For each time period and SSP, loops over environmental variables and averages the corresponding layer across all GCMs, then writes the resulting ensemble SpatRaster to `ensemble.dir/time.period/ssp/`. Rasters are written to file rather than held in memory to avoid excessive memory use over a full run.
#' @examples
#' \dontrun{
#' ensemble_climateFutures(future.env, ensemble.dir = "data/ensembleFutureEnv/")
#' }
#' @export

ensemble_climateFutures <- function(rast.list,
                                    ensemble.dir = "data/ensembleFutureEnv/") {

  if (!dir.exists(ensemble.dir)) {
    dir.create(ensemble.dir, recursive = TRUE)
  }
  # Create env folder for ensembles

  for (time.period in names(rast.list)) {

    cli::cli_inform("Processing time period: {.val {time.period}}")

    for (ssp in names(rast.list[[time.period]])) {

      cli::cli_inform("Creating ensembles of {.val {time.period}}: {.val {ssp}}")

      rast.subset <- rast.list[[time.period]][[ssp]]
      # Subset the period and ssp to get gcm level data

      gcm.layers <- list()

      for (env.var in names(rast.subset[[1]])) {

        cli::cli_inform("Processing variable: {.val {env.var}}")

        for (gcm in names(rast.subset)) {

          cli::cli_inform("Processing GCM: {.val {gcm}}")
          tictoc::tic()

          if (gcm == names(rast.subset)[1]) {
            var.subset <- rast.subset[[gcm]][[env.var]]
          } else {
            var.subset <- var.subset + rast.subset[[gcm]][[env.var]]
            # Add all rasters together
          }

          tictoc::toc()

        }
        # Loop over to add all rasters together

        var.subset <- var.subset / length(names(rast.subset))
        # Find mean average across all GCMs for this layer

        gcm.layers[[env.var]] <- var.subset
        # Collect ensemble layer for this variable

        rm(var.subset)
        gc()

      }

      gcm.subset <- terra::rast(gcm.layers)
      # Build named SpatRaster from ensemble layers

      # rast.list[[time.period]][[ssp]] <- gcm.subset
      # Replace original model level groupings with ensemble outputs
      ## This blew memory up over a full run, so changing to writing to file.

      cli::cli_inform("Saving ensemble to file: {.val {time.period}}, {.val {ssp}}")

      filename <- paste0(time.period, "_", ssp, "_", "ensemble" , ".tif")
      dir <- paste0(ensemble.dir, "/", time.period, "/", ssp, "/")
      # Create filename and dir

      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
      }
      # Create output directory if it does not exist

      terra::writeRaster(gcm.subset, filename = paste0(dir, filename), overwrite = TRUE)

      rm(gcm.subset)
      gc()

    }
    # End ssp loop
  }
  # End time.period loop
  return(rast.list)
}
## Stacking function for muliple future env's
# TODO Could be parallellised one day, but currently runs overnight so thats fine

#' Read Future Climate Rasters into a Nested List
#' @description Reads future climate rasters from a directory tree into a nested list, with list depth auto-detected
#' @param data.dir Character string specifying the parent directory containing time period subdirectories
#' @return Nested list of SpatRasters structured as time.period > ssp, or time.period > ssp > gcm if GCM-level subdirectories are present
#' @details Directory structure is expected as `data.dir/period/ssp/` (ensemble depth) or `data.dir/period/ssp/gcm/` (GCM depth). Depth is auto-detected per ssp based on whether GCM subdirectories exist.
#' @examples
#' \dontrun{
#' read_futureEnv("data/futureEnv/")
#' }
#' @export

read_futureEnv <- function(data.dir) {
  time.periods <- list.dirs(data.dir, recursive = FALSE) |>
    basename()

  future.env <- list()

  for (period in time.periods) {
    future.env[[period]] <- list()
    ssps <- list.dirs(file.path(data.dir, period), recursive = FALSE) |>
      basename()

    for (ssp in ssps) {
      future.env[[period]][[ssp]] <- list()
      gcms <- list.dirs(file.path(data.dir, period, ssp), recursive = FALSE) |>
        basename()

      if (length(gcms) > 0) {
        cli::cli_inform("Reading {.val {period}} / {.val {ssp}}: GCM-level depth detected ({length(gcms)} GCMs)")
        for (gcm in gcms) {
          tif.file <- list.files(file.path(data.dir, period, ssp, gcm),
                                 pattern = "\\.tif$", full.names = TRUE
          )
          future.env[[period]][[ssp]][[gcm]] <- terra::rast(tif.file)
        }
      } else {
        cli::cli_inform("Reading {.val {period}} / {.val {ssp}}: ensemble depth detected")
        tif.file <- list.files(file.path(data.dir, period, ssp),
                               pattern = "\\.tif$", full.names = TRUE
        )
        future.env[[period]][[ssp]] <- terra::rast(tif.file)
      }
    }
  }

  return(future.env)
}
# Read future climate rasters into a nested list; depth auto-detected.
# If GCM subdirectories exist under each SSP, reads period > ssp > gcm.
# If not, reads period > ssp directly.

#' Reorganise Future Predictions into a Nested List
#' @description Splits a named SpatRaster stack of future predictions into a nested list structured as species > time > ssp > model
#' @param pred.rast.future SpatRaster stack of future predictions, with layer names formatted as "species_time_ssp_model"
#' @return Nested list of SpatRaster layers structured as species > time > ssp > model
#' @details Layer names are split on underscores into their four components; assumes a fixed species_time_ssp_model naming convention
#' @examples
#' \dontrun{
#' stack_futureEnvs(pred.rast.future)
#' }
#' @export

stack_futureEnvs <- function(pred.rast.future) {

  layer.names <- names(pred.rast.future)
  future.list <- list()

  for(name in layer.names) {
    message("Reorganising future prediction: ", name)

    # Split name by underscore to get components
    parts <- strsplit(name, "_")[[1]]

    # Extract components (assuming format: species_time_ssp_model)
    species <- parts[1]
    time <- parts[2]
    ssp <- parts[3]
    model <- parts[4]

    # Create nested structure if it doesn't exist
    if(is.null(future.list[[species]])) future.list[[species]] <- list()
    if(is.null(future.list[[species]][[time]])) future.list[[species]][[time]] <- list()
    if(is.null(future.list[[species]][[time]][[ssp]])) future.list[[species]][[time]][[ssp]] <- list()

    # Add raster layer to appropriate nested location
    future.list[[species]][[time]][[ssp]][[model]] <- pred.rast.future[[name]]
  }

  return(future.list)
}
# Function for stacking outputs of multiple future environmental predictions
# Built on a Species -> Time -> SSP -> Model structure
# TODO RDS does not survive saving and loading, probably because of the way terra references and stores rasters in memory
#' Split Future Prediction Naming into Components
#' @description Splits the taxonName column of an sf polygon object into separate taxonName, time.period, ssp, and model columns
#' @param sf.future.thresh sf polygon object with a taxonName column formatted as "taxonName_time.period_ssp_model"
#' @return sf polygon object with taxonName split into separate taxonName, time.period, ssp, and model columns
#' @details Uses `tidyr::separate()` on the taxonName column, retaining the original column
#' @examples
#' \dontrun{
#' split_futureEnvs_sf(sf.future.thresh)
#' }
#' @export

split_futureEnvs_sf <- function(sf.future.thresh) {

  sf.future.thresh <- sf.future.thresh |>
    tidyr::separate(
      taxonName,
      into = c("taxonName", "time.period", "ssp", "model"),
      sep = "_",
      remove = FALSE
    )

  return(sf.future.thresh)
}
# Function for splitting up sf polygon into model outputs if needed

#' Fill NA Values in a Raster Using a Focal Window
#' @description Fills small gaps of NA cells in a raster by averaging surrounding valid cells with a focal window
#' @param rast SpatRaster containing NA gaps to fill
#' @param focal.width Numeric value specifying the width of the focal window (default 9)
#' @param iter Numeric value specifying the number of times to repeat the focal fill (default 1)
#' @param mask.poly sf or SpatVector polygon to mask the output back to, or NA to skip masking (default NA)
#' @return SpatRaster with NA gaps filled, matching the input layer names
#' @details Applies `terra::focal()` with `na.policy = "only"` so that only NA cells are altered. Focal creates a buffering effect at masked boundaries, so an optional mask.poly can be supplied to crop back to original boundaries.
#' @examples
#' \dontrun{
#' rast_naFill(rast, focal.width = 9, mask.poly = study.area)
#' }
#' @export

rast_naFill <- function(rast, focal.width = 9, iter = 1, mask.poly = NA) {

  for (y in seq(iter)) {
    for (i in seq(length(names(rast)))) {
      layer.name <- names(rast)[i]

      message("Filling NA's in rast: ", names(rast)[[i]])
      output <- terra::focal(rast[[i]], w = focal.width, fun = mean, na.policy = "only", na.rm = T)
      # focal parameters to fill in NA cells with surrounding data.
      # Good for patchy data or to fill in NA cells
      # Originally thought of to stop aggregating rasters causing data loss, but na.rm = T in that function also works

      if (!is.na(mask.poly)) {
        message("Masking output raster...")
        output <- terra::mask(output, terra::vect(mask.poly))
      }
      # Mask back to original boundaries if required
      # Focal creates a buffering effect on masked rasters, this might not always be desired.

      names(output) <- layer.name

      rast[[i]] <- output

      gc()
    }
  }

  return(rast)
}
# Function for filling in NA's in datasets to account for small gaps in various rasters
## TODO, bug here where memory leak occurs if running code from the top. Seems ok if rerunning.
# Added gc() to maybe fix, but not tested

#' Replace NA Values in a Raster with Zero over Valid Cells
#' @description Replaces NA cells in a raster with 0 where a reference layer indicates a valid land cell, leaving NAs outside the reference untouched
#' @param x SpatRaster layer containing NA values to replace, where NA indicates absence (e.g. no snow) rather than missing data
#' @param reference SpatRaster used as a land/valid-cell mask
#' @return SpatRaster with NA cells replaced by 0 where reference indicates a valid cell, otherwise left as NA
#' @details Used for layers such as snow cover, where NA represents no snow rather than missing data, and should be distinguished from true missing/no-data cells outside the reference mask
#' @examples
#' \dontrun{
#' rast_naZero(snow.rast, reference = land.mask)
#' }
#' @export

rast_naZero <- function(x, reference) {
  # Replace NAs in x with 0 where reference indicates valid land cells
  na.mask <- terra::mask(x, reference, inverse = TRUE, updatevalue = 0)
  terra::ifel(!is.na(x), x, na.mask)
}
# Helper: fills NAs in a raster layer with 0 over valid land, using a reference layer as land mask
# Used for snow layers where NA indicates no snow (not missing data)

## SDM Object Functions -------------------------------------------------------

#' Load and Repair ENMeval Model Results
#' @description Loads saved ENMeval model results and repairs raster connections
#' @param data.dir Character string specifying path to directory containing model outputs
#' @return List of ENMeval objects with restored raster predictions
#' @details Reads saved model objects from modelOutput/ and raster files from modelOutputRast/
#'          and recombines them into complete ENMeval objects
#' @examples
#' \dontrun{
#' sdm_outputs <- load_results("path/to/outputs/")
#' }
#' @export

load_results <- function(data.dir) {

  message("Loading ENMeval results...")

  sdm.outputs <- list()

  for (i in seq_along(list.files(paste0(data.dir, "modelOutput/")))) {

    file <-  list.files(paste0(data.dir, "modelOutput/"))[i]
    SDM.obj <- readRDS(paste0(data.dir, "modelOutput/", file))
    # Read SDM obj from file

    if (dir.exists(paste0(data.dir, "modelOutputRast/"))) {
      rast <- list.files(paste0(data.dir, "modelOutputRast/"))[i]
      rast <- terra::rast(paste0(data.dir, "modelOutputRast/", rast))
      SDM.obj@predictions <- rast
      # Re-attach raster saved from model outputs
      ## Rebuilds raster predictions for each SDM.obj with raster from file
    }

    sdm.outputs[[i]] <- SDM.obj

    names(sdm.outputs)[i] <- stringr::str_remove(file, "\\.RDS$")
  }

  return(sdm.outputs)
  # Return

}
# Function for loading results in a list and repairing raster connection

#' Model Selection from ENMeval Output
#' @description Selects optimal model from ENMeval outputs using sequential criteria
#' @param sdm.outputs List of ENMeval objects containing model results, as produced by `load_results()`
#' @param data.dir Character string specifying path to directory to save the model selection RDS output
#' @return data.frame with parameters of selected optimal model
#' @details Uses sequential criteria from Kass 2021:
#'   1. Minimizes omission rate (or.10p.avg)
#'   2. Maximizes CBI validation score (cbi.val.avg)
#'   3. If multiple models tie, randomly select one of the models (statistically equivalent)
#' @references Kass et al. 2021 https://doi.org/10.1111/2041-210X.13628, https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0-vignette.html#select
#' @examples
#' \dontrun{
#' # sdm.outputs is a list of ENMeval objects, as produced by load_results()
#' modelSelection <- best_model(sdm.outputs, data.dir = "path/to/outputs/")
#' }
#' @export

best_model <- function(sdm.outputs, data.dir){

  modelSelection <- data.frame(
    index = as.integer(),
    name = character(),
    modelSelection = character()
  )

  for (i in seq(length(sdm.outputs))) {

    SDM.obj <- sdm.outputs[[i]]
    # Select model iteration

    name <- SDM.obj@taxon.name
    # Extract model output name

    result <- ENMeval::eval.results(SDM.obj)
    # Select the results from the objects

    # bestModel <- result %>%
    #   dplyr::filter(!is.na(or.10p.avg)) %>%
    #   dplyr::filter(or.10p.avg == min(or.10p.avg, na.rm = T)) %>%
    #   dplyr::filter(cbi.val.avg == max(cbi.val.avg))
    # # Select based model based on select criteria
    # Trying a new method to help reduce errors

    bestModel <- result |>
      dplyr::filter(!is.na(or.10p.avg)) |>
      dplyr::filter(or.10p.avg == min(or.10p.avg, na.rm = T))
    # Select based model based on select criteria

    if (nrow(bestModel) > 1) {
      if(all(is.na(bestModel$cbi.val.avg))) {

        message(paste0("Best model contains all NA's in cbi.val.avg for ", name), ", selecting randomly...")
        warning(paste0("Best model contains all NA's in cbi.val.avg for ", name), ", selecting randomly...")

        bestModel <- bestModel[sample(nrow(bestModel), size = 1),]

      } else {

        bestModel <- bestModel |>
          dplyr::filter(cbi.val.avg == max(cbi.val.avg, na.rm = T))

        if (length(unique(bestModel$cbi.val.avg)) == 1) {

          message(paste0("Best model contains identical values in cbi.val.avg for ", name), ", selecting randomly...")
          warning(paste0("Best model contains identical values in cbi.val.avg for ", name), ", selecting randomly...")

          bestModel <- bestModel[sample(nrow(bestModel), size = 1),]

        }

      }
    }
    # Added more complex logic to handle model ties and failures
    # If more than one model is required from or.10p, pick maximum cbi.val.avg, if all NA, randomly select.
    # If cbi.values are identical even after picking max value, select randomly

    modelSelection <- dplyr::bind_rows(modelSelection,
                                       data.frame(index = i,
                                                  name = name,
                                                  modelSelection = as.character(bestModel$tune.args)))

  }

  saveRDS(modelSelection, file = paste0(data.dir, "modelSelection.RDS"))

  return(modelSelection)

}
## Function for selection of best model fit from SDM object
## Currently using Kamie Kass 3 step optimal sequential method

## Raster and Thresholding Functions -----------------------------------------------------

#' Write Best Model Raster to File
#' @description Writes the prediction raster from the best performing model to file
#' @param output.dir Character string specifying output directory path
#' @param SDM.obj ENMeval object containing model results
#' @param model Character string specifying the name of the best model settings within `SDM.obj`, as selected by `best_model()`
#' @return None (writes raster file to disk)
#' @examples
#' \dontrun{
#' write_sdm_rast("path/to/output.dir", SDM.obj, model = "fc.LQ_rm.1")
#' }
#' @export

write_sdm_rast <- function(output.dir, SDM.obj, model) {

  if(dir.exists(output.dir) == F){
    dir.create(output.dir, recursive = T)
  }
  # Create output folder

  taxon.name <- SDM.obj@taxon.name
  # Extract model output name

  pred.rast <- ENMeval::eval.predictions(SDM.obj)[[model]]
  # Plot output prediction raster of best model

  names(pred.rast) <- taxon.name
  # Change name from model type to name of taxon

  terra::writeRaster(x = pred.rast, filename = paste0(output.dir, taxon.name, ".tif"), overwrite = T)
  # Write to new output directory

}
# Writes best model fit raster to file

#' Threshold ENMeval Model Predictions
#' @description Applies threshold to continuous model predictions to create binary suitable/unsuitable maps
#' @param SDM.obj ENMeval object containing model results
#' @param output.dir Character string specifying output directory path
#' @param thresh.method Character string specifying threshold method: "maxsss" (default), "or.10.threshold", or "or.min.threshold"
#' @param model Character string specifying the name of the best model settings within `SDM.obj`, as selected by `best_model()`
#' @param filename Character string specifying an optional output filename override. If NULL (default), a filename is derived from `SDM.obj@taxon.name`
#' @param future.pred Optional SpatRaster of a future prediction to threshold instead of the historic `SDM.obj` prediction. If NULL (default), the historic prediction is thresholded
#' @return None (writes thresholded raster to disk)
#' @details Available thresholding methods:
#'   - maxsss: Maximum training sensitivity plus specificity (Liu et al. 2013)
#'   - or.10.threshold: 10% training presence threshold
#'   - or.min.threshold: Minimum training presence threshold
#' @references Liu et al. 2013 https://doi.org/10.1111/jbi.12058
#' @examples
#' \dontrun{
#' write_model_threshold(SDM.obj, "path/to/output/", thresh.method = "maxsss", model = "fc.LQ_rm.1")
#' }
#' @export

write_model_threshold <- function(output.dir, SDM.obj, thresh.method = "maxsss", model, filename = NULL, future.pred = NULL){

  if(dir.exists(output.dir) == F){
    dir.create(output.dir, recursive = T)
  }
  # Create output folder

  taxon.name <- SDM.obj@taxon.name
  # Extract model output name

  pred.rast <- ENMeval::eval.predictions(SDM.obj)[[model]]
  # Output prediction raster of best model

  names(pred.rast) <- taxon.name
  # Change name from model type to name of taxon

  pred.vals <- terra::extract(pred.rast, SDM.obj@occs[,1:2])[,2]
  # Extract prediction values at each occurrence point
  # Used in or.10 and min threshold

  ## Methods of thresholding
  # https://github.com/jamiemkass/ENMeval/issues/53

  if(thresh.method == "maxsss") {

    if(SDM.obj@algorithm == "maxent.jar") {

      ## Method using MaxSSS Method (Jarrod's Method)
      # Liu et al (2013)
      # https://onlinelibrary.wiley.com/doi/10.1111/jbi.12058

      eval <- ENMeval::eval.models(SDM.obj)[[model]]
      # Extract model evaluation metric list

      maxsss <- as.numeric(eval@results['Maximum.training.sensitivity.plus.specificity.Cloglog.threshold',])
      # Extract the Maxsss threshold from model results
      # Maxent defaults to cloglog outputs unless pred.type is changed to another variable.
      # Therfore using this threshold seems appropriate.
      # Was using .area before which might have been something else. Area under the curve statistic perhaps?

      if (!is.null(future.pred)) {
        pred.rast <- future.pred > maxsss
      } else {
        pred.rast <- pred.rast > maxsss
      }
      # Subset the prediction raster based on on defined maxSSS value from model results
      # Seems most defensible method of thresholding
      # Applys to future prediction if provided

    }

    if (SDM.obj@algorithm == "maxnet") {

      message("TODO Still working out maxnet implimentation of maxsss")
      warning("TODO Still working out maxnet implimentation of maxsss")
      # https://consbiol-unibern.github.io/SDMtune/reference/thresholds.html#details
      # TODO Try thresholds?
    }

  }

  if(thresh.method == "or.10.threshold") {

    ## Method 1: 10 percent presence threshold
    # https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/

    if (length(pred.vals) < 10) {
      p10 <- floor(length(pred.vals) * 0.9)
    } else {
      p10 <- ceiling(length(pred.vals) * 0.9)
    }
    threshold <- rev(sort(pred.vals))[p10]

    if (!is.null(future.pred)) {
      pred.rast <- future.pred > threshold
    } else {
      pred.rast <- pred.rast > threshold
    }
    # This method assumes 10 percent lowest occurrence points aren't suitable.
    # More conservative than or.minimum.threshold

  }

  if(thresh.method == "or.min.threshold") {

    ## Alternative Method: Minimum omission threshold
    # https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold
    # Not very conservative, this represents the absolute minimum possible prediction based on model.

    if (!is.null(future.pred)) {
      pred.rast <- future.pred > min(pred.vals)
    } else {
      pred.rast <- pred.rast > min(pred.vals)
    }
    # Calculated absolute minium threshold

  }

  ## Write to file ##
  if (!is.null(filename)) {
    terra::writeRaster(x = pred.rast, filename = paste0(output.dir, filename, ".tif"), overwrite = T)
    # Write to new output directory as rast
  } else {
    terra::writeRaster(x = pred.rast, filename = paste0(output.dir, taxon.name, ".tif"), overwrite = T)
    # Write to new output directory as rast
  }

}

#' Convert Thresholded Rasters to SF Polygons
#'
#' @description Converts a directory of thresholded binary rasters into a single SF polygon object
#'
#' @param thresh.dir Character string specifying directory containing thresholded rasters
#' @param output.dir Character string specifying directory for output files
#' @param output.name Character string specifying name for output SF file
#' @param simp.factor Numeric simplification factor passed to `sf::st_simplify()` when converting rasters to polygons
#' @details
#' Converts each binary raster in the input directory to polygons and combines them
#' into a single SF object. Only keeps polygons where raster value equals 1
#' (suitable habitat). Raster filenames are used as taxon names in output.
#' @examples
#' \dontrun{
#' convert_thresh_toSf("data/thresholded/", "output/polygons/", "habitat_polygons.shp")
#' }
#' @export

convert_thresh_toSf <- function(thresh.dir, output.dir, output.name, simp.factor) {

  if (length(list.files(thresh.dir)) > 0) {

    if(dir.exists(output.dir) == F){
      dir.create(output.dir, recursive = T)
    }
    # Create output dir if none exists.

    for (i in seq_along(list.files(thresh.dir))) {

      taxon.name <- stringr::str_remove(list.files(thresh.dir)[i], ".tif")
      # Extract taxon name

      message("Converting threshold raster to sf polygon for: ", taxon.name)
      # User feedback

      pred.rast <- terra::rast(paste0(thresh.dir, taxon.name, ".tif"))
      # Load relevant rast

      poly <- terra::as.polygons(pred.rast) |>
        sf::st_as_sf()
      poly.crs <- sf::st_crs(poly)
      # Cast rast to sf and capture crs

      names(poly)[1] <- "taxonName"
      poly <- poly |>
        dplyr::filter(taxonName == 1)
      # Col name correction and drop polygon which == F (0 Value)

      if (nrow(poly) != 0) {
        # Usual operation if polygon is generated from raster

        if (i == 1) {
          sf <- data.frame(
            taxonName = taxon.name,
            geometry = sf::st_geometry(poly)
          ) |>
            sf::st_as_sf(crs = poly.crs)
        } else {
          sf.row <- data.frame(
            taxonName = taxon.name,
            geometry = sf::st_geometry(poly)
          ) |>
            sf::st_as_sf(crs = poly.crs)
          sf <- rbind(sf, sf.row)
        }
        # Either create a new sf dataframe from poly and taxon.name, or add and amend a row.

      } else {
        # Safety catch to write an empty polygon to the row if no geometry is generated

        if (i == 1) {
          sf <- data.frame(
            taxonName = taxon.name,
            geometry = sf::st_geometry(sf::st_multipolygon())
          ) |>
            sf::st_as_sf(crs = poly.crs)
        } else {
          sf.row <- data.frame(
            taxonName = taxon.name,
            geometry = sf::st_geometry(sf::st_multipolygon())
          ) |>
            sf::st_as_sf(crs = poly.crs)
          sf <- rbind(sf, sf.row)
        }
        # Either create a new sf dataframe from poly and taxon.name, or add and amend a row.

      }

    }
    # Loop for extracting polygons from each raster in threshold.dir and making a single sf object

    sf <- sf::st_simplify(x = sf,
                          preserveTopology = T,
                          dTolerance = simp.factor)
    # Simplify geometry to reduce size of output

    sf <- sf |>
      dplyr::mutate(area_km2 = sf::st_area(sf) |>
                      units::set_units("km^2") |>
                      as.numeric())
    # Amend area for potential other use cases

    message("Writing sf polygon to file...")
    sf::st_write(sf, dsn = paste0(output.dir, output.name))

  } else {
    warning("Threshold directory is empty for convert_thresh_toSf()")
  }

}

# TODO Made this from a previous project to help slotting rasters back in if raster.pred == F
# Not sure if it needs to be maintained but it's here just in case

## Graphing Functions ------------------------------------------------------

## TODO Make a graphing pipeline to speed up all standard plot generations

#' Graph Site Map (stub)
#' @description Not yet implemented. Placeholder for plotting a site map of occurrence records over a basemap.

graph_siteMap <- function(){

  basemap <- basemap_ggplot(st_bbox(ibra.alps), map_service = "osm", map_type = "streets_de")

}

#' Graph Historical Climate (stub)
#' @description Not yet implemented. Placeholder for plotting historical climate summaries.

graph_historicalClim <- function(){

}

#' Graph Future Climate (stub)
#' @description Not yet implemented. Placeholder for plotting future climate summaries.

graph_futureClim <- function(){

}


## Metadata Functions ------------------------------------------------------

## Look at ENMeval

## TODO need a framework for capturing metadata of a model run, so it can be provided alongside any raster products
## Look at other report and work and use as template?
