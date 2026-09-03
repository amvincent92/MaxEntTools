# ENMeval Workflow  -----------------------------------------------------------------

# Main default and use of ENMeval::ENMevaluate functions in a convenient pipeline
# Includes all standard processes needed on occurrence data for model generation
# Consideration of paramaterisation should be made for each use case and question being asked

## Notes:

# Built on ENMeval 2.0.5, maxent.jar v3.4.4, dismo 1.3-14
# rJava and Java versions are 1.0-11, and Version 24) respectively.
# If running java not possible, use the "maxnet" algorithm
# This may have downstream issues with model validation and prediction.

# ENMeval will remove occurrence localities that share the same grid cell.
# But often good practice to conduct all desired data cleaning outside of the pipeline

# Can ignore warnings for closing connections when parallel processing

## Known Bugs

# If running prediction (and validation?) pipelines not sequentially, java config becomes and issue
## Injecting it in each function has caused bugs in past, would be good to to maybe get a conditional statement working if time.

# Extract env_config out of functions as it regularly gets annoying if you need to rerun stuff if failed mid function

## Bottlenecks
# Major bottleneck in validation pipeline, writing of response curves in the maxent.jar iteration
# Major bottleneck in making model prediction rasters which could be turned off for model fitting stage

## TODO Development plans
# Predication pipeline is approaching needing big refactor
# Summary plot of all diagnostics?
# Metadata capture
# ROC Curves
# Maxsss for maxnet
# Fix maxnet response curves into a ggplot format.

# Raster pred = F
## Would be nice to train on high def data and not save all different raster outputs
## However, might be a big rework
# Package up at the end? Have set it up reasonably well structured for this.

# Ensembling tools could be brought in from snow gum phd project.
# Reading from a directory of rasters might be a better option than saving all as a single geotiff.
## Getting some issues in reading back super large rasters and it double saves unneccesarly.
# Add denoise function from ACT government work
# Data checks and prep could be its own function pre starting the modelling pipeline

# TODO Could really pull up and paramaterise directory management stuff

## TODO BUG Error in loading the SDM object via RDS. I think it breaks the saved raster inside and
# therefore can't handle re-running the future predictions?
# This might not have ever come up before because you always run in a straight line?
# To confirm.


# Active  -----------------------------------------------------------------

## Data Validation Pipeline ----------------------------------------------------

#' Maxent Data Validation Pipeline
#' @description Pipeline for validating occurrence, environmental, background area, and study area inputs before model generation
#' @param occ.sf.pts sf Points object, occurrence records
#' @param env.vars SpatRaster, environmental variables
#' @param bg.areas sf Polygon, background area polygons, such as ibra bioregions etc
#' @param study.area sf Polygon, study area polygon, must contain all occurrence points
#' @param data.dir Character, path to data directory
#' @param output.dir Character, path to output directory
#' @param subset.col Character, column name for subsetting data by species (default "scientificName")
#' @return None. Called for side effects: prints warnings, and writes an environmental correlation plot to disk
#' @details Checks that occurrence records fall within the study area, checks each species has a minimum number of records via `check_numRecords()`, and writes an environmental correlation plot via `check_envCorr()`
#' @export

maxent_dataValidation_pipe <- function(occ.sf.pts,
                                       env.vars,
                                       bg.areas,
                                       study.area,
                                       data.dir = "data/maxent/",
                                       output.dir = "outputs/maxent/",
                                       subset.col = "scientificName") {

  ## TODO Fill out as testing progresses
  ## Check that occurrence records, env and bg.areas, and study area are all in alignment

  if (dir.exists(paste0(data.dir))) {
    warning("Output directory detected, resuming partially run pipelines can result in errors...")
  }


  if (!all(st_intersects(occ.sf.pts$geometry, study.area, sparse = F))) {
    warning("Occurrence records exist outside of study area of interest")
  }
  # Test if any occurrence records exist outside of defined study aoi
  # Important as any points drawing from NA data may cause errors
  # Also forces correct usage. No occurrence records should be outside study area
  # However, using bg.area and env data is ok as it is cropped and intersected in following two functions for convenience

  check_numRecords(occ.sf.pts, check.pts = 100, subset.col = subset.col)
  # Warns user if any species have less than 100 records

  ### Env Correlations ----------------------------------------

  check_envCorr(output.dir = paste0(output.dir, "modelValidation/"),
                env.vars = env.vars,
                occ.sf.pts = occ.sf.pts)
  # Saves a graph of environmental correlations

  ## TODO Load VIF results and display to user

}
# TODO Extracting these tests from Main pipeline to future proof


## Model Generation Pipeline -----------------------------------------------

#' ENMeval Pipeline Main Function
#' @description Pipeline for standard workflow of generating maxent models from occurrence points, environmental raster stack, background areas and overall study area.
#' @param occ.sf.pts sf Points object, occurrence records
#' @param env.vars SpatRaster, environmental variables
#' @param bg.areas sf Polygon, background area polygons, such as ibra bioregions etc
#' @param study.area sf Polygon, study area polygon, must contain all occurrence points.
#' @param data.dir Character, path to data directory; subfolders (biasLayer, thinnedRecords, bgAreas, bgPoints, modelOutput, modelOutputRast) are created beneath it for intermediate and final model outputs
#' @param seed Numeric, random seed for reproducibility
#' @param java.mem.gb Numeric, Java memory in GB
#' @param bg.buffer.dist Numeric, buffer distance for background points
#' @param bg.pts.num Numeric, number of background points
#' @param subset.col Character, column name for subsetting data
#' @param clampStudyArea Logical, restrict background areas to within area of interest to reduce computation time
#' @param raster.preds Logical, whether to generate and save raster predictions for each model (default TRUE). If FALSE, the modelOutputRast subfolder is not created
#' @param doClamp Logical, whether to clamp predictions to avoid extrapolation beyond training data
#' @param categoricals Character vector, use to tell MaxEnt which of your variables are categorical. Used terra::rasterise and assign field for factor levels. Do not aggregate post as this will break layer for enmeval.
#' @param partitions Character, type of spatial partitioning
#' @param algorithm Character, modeling algorithm
#' @param overlap Logical, whether to calculate niche overlap
#' @param overlapStat Character, statistic for niche overlap
#' @param tune.args List, model tuning arguments
#' @param enable.parallel Logical, enable parallel processing
#' @param numCores Numeric, number of cores for parallel processing
#' @return None. Saves model outputs to disk
#' @export

#
##
###

maxent_modelGeneration_pipe <- function(occ.sf.pts,
                                        env.vars,
                                        bg.areas,
                                        study.area,
                                        data.dir,
                                        # User defined variables
                                        seed = 343,
                                        java.mem.gb = 16,
                                        bg.buffer.dist = NULL,
                                        bg.pts.num = 10000,
                                        subset.col = "scientificName",
                                        clampStudyArea = T,
                                        # Default Pipeline parameters
                                        raster.preds = T,
                                        doClamp = T,
                                        categoricals = NULL,
                                        partitions = "block",
                                        algorithm = "maxent.jar",
                                        overlap = T,
                                        overlapStat = c("D", "I"),
                                        tune.args = list(fc = c("L","LQ","LQH", "LQHP", "LQHPT"),
                                                         rm = c(0.5,1,2,4,5)),
                                        enable.parallel = F,
                                        numCores = 4
                                        # Default maxent function settings
) {

  env_config(algorithm = algorithm,
             seed = seed,
             java.mem.gb = java.mem.gb)
  # Run config function to prepare environment and dependencies

  ## Data Prep ---------------------------------------------------------------
  ### Subset to Study Area ----------------------------------------------------

  if (clampStudyArea == T) {
    bg.areas <- subset_BGpointsByStudyArea(bg.areas = bg.areas, study.area = study.area)
    # Intersects background area data with study area to reduce processing time

    env.vars <- subset_EnvByStudyArea(env.vars = env.vars, study.area = study.area)
    # Crops environmental variables by study area extent to reduce processing time
  } else {
    message("clampStudyArea == F: Environmental variables and background area not subset to by study area...")
  }

  ### Create Bias Layer --------------------------------------------------------------

  if(dir.exists(paste0(data.dir, "biasLayer/")) == F){
    dir.create(paste0(data.dir, "biasLayer/"), recursive = T)
  } else {
    message("SDM bias layer detected...")
  }
  # Output folder for bias layer

  if (file.exists(paste0(data.dir, "biasLayer/biasLayer.rds")) == F) {
    bias.layer <- bias_layer_rast(data.sf = occ.sf.pts,
                                  env.vars = env.vars,
                                  bg.poly = study.area,
                                  output.dir = paste0(data.dir, "biasLayer/"))
  } else {
    message("Bias layer detected, reading from file...")
    bias.layer <- readRDS(paste0(data.dir, "biasLayer/biasLayer.rds"))
  }

  ### Thin Records ------------------------------------------------------------

  if(dir.exists(paste0(data.dir, "thinnedRecords/")) == F){
    dir.create(paste0(data.dir, "thinnedRecords/"), recursive = T)
  } else {
    message("Thinned occurrence records detected...")
  }
  # Output folder for thinned occurrence records

  if (file.exists(paste0(data.dir, "thinnedRecords/thinnedSf_byRast.RDS")) == F) {

    occ.sf.pts <- thin_records_sf_rast(data.sf = occ.sf.pts,
                                       subset.col = subset.col,
                                       ref.rast = env.vars,
                                       output.dir = paste0(data.dir, "thinnedRecords/"))
  } else {
    message("Thinned records by raster detected, reading from file...")
    occ.sf.pts <- readRDS(paste0(data.dir, "thinnedRecords/thinnedSF_byRast.RDS"))
  }
  # If thinned records don't exist, thin them and save output

  check_numRecords(occ.sf.pts, check.pts = 100, subset.col = subset.col)
  # Warns user if any species have less than 100 records after thinning

  ## Background Areas --------------------------------------------------------

  if(dir.exists(paste0(data.dir, "bgAreas/")) == F){
    dir.create(paste0(data.dir, "bgAreas/"), recursive = T)
  } else {
    message("SDM BG areas detected...")
  }
  # Output folder for background areas

  if (file.exists(paste0(data.dir, "bgAreas/bgAreas.RDS")) == F) {
    bg.areas <- define_background_area(bg.areas = bg.areas,
                                       occ.sf.pts = occ.sf.pts,
                                       study.area = study.area,
                                       subset.col = subset.col,
                                       bg.buffer.dist = bg.buffer.dist,
                                       simplify.poly = T,
                                       output.dir = paste0(data.dir, "bgAreas/"))
  } else {
    message("Background areas detected, reading from file...")
    bg.areas <- readRDS(paste0(data.dir, "bgAreas/bgAreas.RDS"))
  }
  # If thinned records don't exist, thin them and save output

  ## Background Points -------------------------------------------------------

  if(dir.exists(paste0(data.dir, "bgPoints/")) == F){
    dir.create(paste0(data.dir, "bgPoints/"), recursive = T)
  } else {
    message("SDM BG points detected...")
  }
  # Output folder for background Points

  if (file.exists(paste0(data.dir, "bgPoints/bgPoints.RDS")) == F) {
    bg.points <- sample_bg_points(bg.areas = bg.areas,
                                  bias.layer = bias.layer,
                                  env.vars = env.vars,
                                  sample.size = bg.pts.num,
                                  output.dir = paste0(data.dir, "bgPoints/"))

    # bg.points.thinned <- thin_records_sf_rast(bg.points,
    #                                           subset.col = "subset_value",
    #                                           ref.rast = env.vars,
    #                                           output.dir = paste0(data.dir, "bgPoints/"))
    ## Could thin bg points to stop occasional overlap between bg points to same rasters.
    ## Needs some alteration to thin records for each of the output sf dataframe
    ## ENMeval does checks for you so not currently essential... Main records thinned above

  } else {
    message("Background points detected, reading from file...")
    bg.points <- readRDS(paste0(data.dir, "bgPoints/bgPoints.RDS"))
  }
  # If thinned records don't exist, thin them and save output

  ## ENM Eval Function -------------------------------------------------------

  if(dir.exists(paste0(data.dir, "modelOutput/")) == F){
    dir.create(paste0(data.dir, "modelOutput/"), recursive = T)
  } else {
    message("SDM Output folder detected...")
  }
  ## Output folder for SDM objects

  if (raster.preds == T) {
    if(dir.exists(paste0(data.dir, "modelOutputRast/")) == F){
      dir.create(paste0(data.dir, "modelOutputRast/"), recursive = T)
    } else {
      message("SDM Raster Output folder detected...")
    }
  } else {
    message("Raster predictions turned off... Skipping directory creation...")
  }
  ## Output folder for SDM raster outputs

  for (i in seq(nrow(unique(as.data.frame(occ.sf.pts)[subset.col])))) {

    sp.name <- unique(as.data.frame(occ.sf.pts)[subset.col])[i,1]

    occ.subset <- occ.sf.pts |>
      dplyr::filter(.data[[subset.col]] == sp.name)
    # Select all records which the equal the subset column equals the unique scientific name subset to the loop,
    occ.subset <- as.data.frame(occ.subset) |>
      dplyr::select(decimalLongitude, decimalLatitude)
    # Subset occurrence records to species of interest
    # Prep for passing to ENMeval

    bg.points.subset <- bg.points |>
      dplyr::filter(.data$subset_value == sp.name) |>
      sf::st_coordinates() |>
      as.data.frame() |>
      dplyr::select(X, Y) |>
      dplyr::rename(decimalLongitude = X,
             decimalLatitude = Y)
    # Subset background points dataframe to species of interest
    # Prep for passing to ENMeval

    ## Debug Test Plots ##

    # plot(occ.subset)
    # plot(bg.points.subset)
    # plot(env.vars)
    #
    # results <- ENMevaluate(
    #   occs = occ.subset,
    #   envs = env.vars,
    #   bg = bg.points.subset,
    #   algorithm = "maxent.jar",
    #   tune.args = list(fc = "L", rm = 1),  # Simple settings
    #   partitions = "block",
    #   parallel = FALSE
    # )

    ## ENM Eval ##
    SDM.obj <- ENMeval::ENMevaluate(occs = occ.subset,
                                    envs = env.vars,
                                    bg = bg.points.subset,
                                    tune.args = tune.args,
                                    partitions = partitions,
                                    taxon.name = sp.name,
                                    overlap = overlap,
                                    overlapStat = overlapStat,
                                    algorithm = algorithm,
                                    numCores = numCores,
                                    raster.preds = raster.preds,
                                    parallel = enable.parallel,
                                    categoricals = categoricals
    )
    # Updated to ENMeval v2.0.5, maxent.jar v3.4.3, predicts v0.1.17
    # Needs available C disk space to write temporary files too.
    # Needs java to be correctly configured
    # Will only run once without resetting R in between

    message("Saving SDM.obj to file...")
    saveRDS(SDM.obj, file = paste0(data.dir, "modelOutput/", sp.name, ".RDS"))
    ## Save Model output

    if (raster.preds == T) {

      message("Saving prediction raster stack to file...")
      terra::writeRaster(SDM.obj@predictions,
                         filename = paste0(data.dir, "modelOutputRast/", sp.name, ".tif"),
                         overwrite = T)
      ## Save raster output
      # Need to do this as an R appears to break raster connections when writing to RDS as a whole.

    } else {
      message("Raster.preds = F, cannot write prediction rasters to file...")
    }



    gc(full = T, reset = T)

  }
  # Loop for generating a model for each species

}

###
##
#

## Model Validation Pipeline -----------------------------------------------------

#
##
###

#' Validation Pipeline for Maxent Models
#' @description Runs validation analyses and generates output files for Maxent models created with maxent_pipe()
#' @param env.vars SpatRaster, environmental variables used in original models
#' @param occ.sf.pts sf Points object, occurrence records used in original models
#' @param data.dir Character, path to directory containing model outputs
#' @param output.dir Character, path to directory for validation outputs
#' @param raster.preds Logical, whether prediction rasters were generated for each model (default TRUE); controls whether prediction rasters are re-written and thresholded
#' @param thresh.method Character, threshold method passed to `write_model_threshold()`: "maxsss" (default), "or.10.threshold", or "or.min.threshold"
#' @return None. Saves validation outputs to disk:
#'   - Environmental correlations
#'   - Prediction rasters
#'   - Variable contribution plots
#'   - Response curves
#'   - Thresholded predictions
#'   - Stacked prediction rasters
#' @export

maxent_modelValidation_pipe <- function(env.vars, occ.sf.pts,
                                        data.dir = "data/maxent/",
                                        output.dir = "outputs/maxent/",
                                        raster.preds = T,
                                        thresh.method = "maxsss") {

  sdm.outputs <- load_results(data.dir)

  # Select Best Model -------------------------------------------------------

  modelSelection <- best_model(sdm.outputs = sdm.outputs,
                               data.dir = data.dir)
  # Selects best fitting model

  # Extract & Thresh Historic Models -------------------------------------------------

  ## Extract prediction raster for each species ##

  if (raster.preds == T) {
    for (i in seq(length(sdm.outputs))) {

      SDM.obj <- sdm.outputs[[i]]
      # Select model iteration

      taxon.name <- SDM.obj@taxon.name
      # Extract model output name

      model <- modelSelection |>
        dplyr::filter(name == taxon.name) |>
        dplyr::select(modelSelection) |>
        as.character()
      # Selects best model from model selection table

      message("Writing and thresholding historic model raster: ", taxon.name)

      write_sdm_rast(output.dir = paste0(data.dir, "modelOutputSelect/predictionRast/"),
                     SDM.obj = SDM.obj,
                     model = model)
      # Write prediction raster to output folder from model outputs

      write_model_threshold(output.dir = paste0(data.dir, "modelOutputSelect/thresholdRast/"),
                            SDM.obj = SDM.obj,
                            model = model,
                            thresh.method = thresh.method)
      # Write prediction raster to output folder
      # Appears to work with both maxnet and maxent.jar
      # TODO MaxSSS Thresholding not working with Maxnet

    }

    ## Stack prediction rasters ##

    if(dir.exists(paste0(output.dir, "outputRasts/")) == F){
      dir.create(paste0(output.dir, "outputRasts/"), recursive = T)
    }
    # Create output raster directory

    message("Writing prediction raster stack...")
    rast_files <- list.files(paste0(data.dir, "modelOutputSelect/predictionRast/"),
                             pattern = "\\.tif$",
                             full.names = TRUE)
    predict.stack <- terra::rast(rast_files)
    terra::writeRaster(predict.stack, filename = paste0(output.dir, "outputRasts/", "predictRastStack.tif"), overwrite = T)
    # Stack all rasters into a singular stack for return to user

    ## Stack prediction thresholded rasters ##

    message("Writing prediction thresholded raster stack...")
    rast_files <- list.files(paste0(data.dir, "modelOutputSelect/thresholdRast/"),
                             pattern = "\\.tif$",
                             full.names = TRUE)
    predict.stack.thresh <- terra::rast(rast_files)
    terra::writeRaster(predict.stack.thresh, filename = paste0(output.dir, "outputRasts/", "predictThresholdRastStack.tif"), overwrite = T)
    # Stack all rasters into a singular stack for return to user

  } else {
    message("Raster predicitions not present in SDM.obj, skipping raster extraction...")
  }

  # Model Validation Figures ------------------------------------------------

  ## Make model validation graphs

  for (i in seq(length(sdm.outputs))) {

    SDM.obj <- sdm.outputs[[i]]
    # Select model iteration
    taxon.name <- SDM.obj@taxon.name
    # Extract model output name
    model <- modelSelection |>
      dplyr::filter(name == taxon.name) |>
      dplyr::select(modelSelection) |>
      as.character()
    # Selects best model from model selection table

    message("Writing validation graphs: ", taxon.name)

    write_var_contrib(SDM.obj = SDM.obj,
                      model = model,
                      output.dir = paste0(output.dir, "modelValidation/varContrib/"),
                      data.dir = paste0(data.dir, "varContribution/")
    )
    # Write varContrib graph to output folder
    # Not available with maxnet models

    write_response_curves(SDM.obj = SDM.obj,
                          model = model,
                          output.dir = paste0(output.dir, "modelValidation/responseCurves/"),
                          data.dir = paste0(data.dir, "responseVariables/")
    )
    # Write response curve graph to output folder

    ## TODO ROC Curves?

  }

  # Metadata ----------------------------------------------------------------

  # Function for generating appropriate metadata file
  # TODO
}

###
##
#

## Future Prediction Pipeline ------------------------------------------------------

## Built on Maxent.jar output, but should handle maxnet ok.

#' Maxent Future Prediction Pipeline
#' @description Pipeline for projecting saved MaxEnt models onto future climate scenarios and writing thresholded prediction rasters to disk
#' @param env.vars.future SpatRaster or nested list of SpatRasters (time.period > ssp, or time.period > ssp > gcm) of future environmental variables
#' @param ensemble.climates Logical, whether future climate rasters have been ensembled across GCMs (default TRUE)
#' @param thresh.method Character, thresholding method passed to `write_model_threshold()` (default "maxsss")
#' @param data.dir Character, path to data directory containing saved model outputs
#' @param output.dir Character, path to output directory
#' @param ensemble.dir Character, path to directory containing ensembled climate rasters (default "data/ensembleFutureEnv")
#' @return None. Writes future prediction and threshold rasters, and stacked prediction rasters, to disk
#' @details Loads saved SDM outputs and model selection choices, then for each taxon projects the selected model (via `predict_maxnet()` or `predict_maxentJar()`) onto the supplied future environment(s), writing individual and stacked prediction and threshold rasters to disk. Handles a single SpatRaster, a time > ssp ensemble list, or a time > ssp > gcm list.
#' @export

maxent_futurePrediction_pipe <- function(env.vars.future,
                                         ensemble.climates = T,
                                         thresh.method = "maxsss",
                                         data.dir = "data/maxent/",
                                         output.dir = "outputs/maxent/",
                                         ensemble.dir = "data/ensembleFutureEnv") {

  sdm.outputs <- load_results(data.dir)
  # Load outputs into a list
  modelSelection <- readRDS(paste0(data.dir, "modelSelection.RDS"))
  # Load model selection choices from validation pipeline

  ## TODO Build each of these processes into their own function
  ## TODO Save raw outputs in data, final stacked product in outputs.

  # Generate Future Predictions  ------------------------------------------------------

    if(dir.exists(paste0(output.dir, "outputRasts/")) == F){
      dir.create(paste0(output.dir, "outputRasts/"), recursive = T)
    }
    # Create parent output folder

    ## Prediction Loop ##

    message("Creating future environmental predictions...")

    for (i in seq(length(sdm.outputs))) {

      SDM.obj <- sdm.outputs[[i]]
      # Select model iteration
      taxon.name <- SDM.obj@taxon.name
      # Extract model output name
      model <- modelSelection |>
        dplyr::filter(name == taxon.name) |>
        dplyr::select(modelSelection) |>
        as.character()
      # Selects best model from model selection table

        if(dir.exists(paste0(data.dir, "futurePredictions/futurePredictionsRast/")) == F){
          dir.create(paste0(data.dir, "futurePredictions/futurePredictionsRast/"), recursive = T)
        }
        if(dir.exists(paste0(data.dir, "futurePredictions/thresholdFuturePredictionsRast/")) == F){
          dir.create(paste0(data.dir, "futurePredictions/thresholdFuturePredictionsRast/"), recursive = T)
        }
        # Create output raster directories

        message("Creating future prediction for: ", taxon.name)

      ## TODO Will need some check to tailor ensemble not ensemble distinction.
      # Worth it to keep the pipeline simple though, ideally no more branches of features

      ## Single Future Env -------------------------------------------------------

      if (inherits(future.env, "SpatRaster")) {
        # If spat raster, likely just a single future

        future.pred <- switch(SDM.obj@algorithm,
                              "maxnet"     = predict_maxnet(env.vars = env.vars.future, SDM.obj = SDM.obj,
                                                            model = model, taxon.name = taxon.name,
                                                            clamp.env = clamp.env),
                              "maxent.jar" = predict_maxentJar(env.vars = env.vars.future, SDM.obj = SDM.obj,
                                                               model = model, taxon.name = taxon.name),
                              stop("Unknown algorithm: ", SDM.obj@algorithm)
        )
        # Takes best model and outputs a prediction onto the future environment

        message("Writing future prediction raster to file: ", taxon.name)

        terra::writeRaster(x = future.pred, filename = paste0(data.dir, "futurePredictions/futurePredictionsRast/", taxon.name, ".tif"), overwrite = T)
        ## Write raster output

        write_model_threshold(output.dir = paste0(data.dir, "futurePredictions/thresholdFuturePredictionsRast/"),
                              SDM.obj = SDM.obj,
                              model = model,
                              thresh.method = thresh.method,
                              future.pred = future.pred)
        ## Write threshold raster output

      }

      ## Multiple Future Env -----------------------------------------------------

      if (is.list(future.env)) {
        # If class is a list, likely a collection of future env stacks

        ## Detect list depth: gcm (time > ssp > gcm) or ensemble (time > ssp) ##

        first.ssp <- future.env[[1]][[1]]

        if (inherits(first.ssp, "SpatRaster")) {
          list.depth <- "ensemble"
        } else if (is.list(first.ssp) && inherits(first.ssp[[1]], "SpatRaster")) {
          list.depth <- "gcm"
        } else {
          cli::cli_abort(
            "Future environment list structure not recognised. \\
       Expected time -> scenario -> SpatRaster or time -> scenario -> gcm -> SpatRaster."
          )
        }

        cli::cli_inform("Future environment list depth detected: {.val {list.depth}}")

        ## Loop for drilling down on listed future env structure ##

        for (time.period in names(future.env)) {
          for (scenario in names(future.env[[time.period]])) {

            if (list.depth == "gcm") {
              for (model.name in names(future.env[[time.period]][[scenario]])) {

                cli::cli_inform("Predicting future climate: {.val {time.period}} / {.val {scenario}} / {.val {model.name}}")

                env.vars.future <- future.env[[time.period]][[scenario]][[model.name]]
                # Extract relevant raster stack

                pred.name <- paste(taxon.name, time.period, scenario, model.name, sep = "_")
                # Filename component shared across predict calls and outputs

                future.pred <- switch(SDM.obj@algorithm,
                                      "maxnet"     = predict_maxnet(env.vars = env.vars.future, SDM.obj = SDM.obj,
                                                                    model = model, taxon.name = pred.name,
                                                                    clamp.env = clamp.env),
                                      "maxent.jar" = predict_maxentJar(env.vars = env.vars.future, SDM.obj = SDM.obj,
                                                                       model = model, taxon.name = pred.name),
                                      cli::cli_abort("Unknown algorithm: {.val {SDM.obj@algorithm}}")
                )
                # Use predict function matching algorithm; filename accounts for multiple future scenarios

                cli::cli_inform("Writing future prediction raster to file...")

                output.path <- paste0(data.dir, "futurePredictions/futurePredictionsRast/",
                                      pred.name, ".tif")
                terra::writeRaster(future.pred, filename = output.path, overwrite = TRUE)
                # Write raster to file with gcm-level name amendments

                write_model_threshold(output.dir = paste0(data.dir, "futurePredictions/thresholdFuturePredictionsRast/"),
                                      SDM.obj = SDM.obj,
                                      model = model,
                                      thresh.method = thresh.method,
                                      future.pred = future.pred,
                                      filename = pred.name)
                # Threshold output model w/ gcm-level filename amendments
              }

            } else {
              # list.depth == "ensemble": no gcm loop, read directly at ssp level

              cli::cli_inform("Predicting future climate: {.val {time.period}} / {.val {scenario}}")

              env.vars.future <- future.env[[time.period]][[scenario]]
              # Extract relevant raster stack

              pred.name <- paste(taxon.name, time.period, scenario, sep = "_")
              # Filename component shared across predict calls and outputs

              future.pred <- switch(SDM.obj@algorithm,
                                    "maxnet"     = predict_maxnet(env.vars = env.vars.future, SDM.obj = SDM.obj,
                                                                  model = model, taxon.name = pred.name,
                                                                  clamp.env = clamp.env),
                                    "maxent.jar" = predict_maxentJar(env.vars = env.vars.future, SDM.obj = SDM.obj,
                                                                     model = model, taxon.name = pred.name),
                                    cli::cli_abort("Unknown algorithm: {.val {SDM.obj@algorithm}}")
              )
              # Use predict function matching algorithm; filename accounts for multiple future scenarios

              cli::cli_inform("Writing future prediction raster to file...")

              output.path <- paste0(data.dir, "futurePredictions/futurePredictionsRast/",
                                    pred.name, ".tif")
              terra::writeRaster(future.pred, filename = output.path, overwrite = TRUE)
              # Write raster to file with ensemble-level name amendments

              write_model_threshold(output.dir = paste0(data.dir, "futurePredictions/thresholdFuturePredictionsRast/"),
                                    SDM.obj = SDM.obj,
                                    model = model,
                                    thresh.method = thresh.method,
                                    future.pred = future.pred,
                                    filename = pred.name)
              # Threshold output model w/ ensemble-level filename amendments
            }
          }
        }
        # Updated loop to account for both types of future.env's, ensemble and gcm.
      }

      ## Depreciated original version ##


      # if (class(future.env) == "list") {
      #   # If class is a list, likely a collection a future env stacks
      #
      #   if(class(future.env[[1]][[1]][[1]]) != "SpatRaster"){
      #     stop("Future environment structure may not match expected list depth/format: time -> scenario -> model")
      #   }
      #   # Stops function if raster list is not exactly 3 levels deep.
      #   # Not necessarily the best test but a stopgap for now.
      #
      #   ## Loop for drilling down on a listed future env structure ##
      #
      #   for (time.period in names(future.env)) {
      #     for (scenario in names(future.env[[time.period]])) {
      #       for (model.name in names(future.env[[time.period]][[scenario]])) {
      #
      #         message(paste("Predicting future climate prediction for:", time.period, scenario, model.name))
      #         # Message user
      #
      #         env.vars.future <- future.env[[time.period]][[scenario]][[model.name]]
      #         # Extract relevant raster stack
      #
      #         pred.name <- paste(taxon.name, time.period, scenario, model.name, sep = "_")
      #         # Filename component shared across predict calls and outputs
      #
      #         future.pred <- switch(SDM.obj@algorithm,
      #                               "maxnet"     = predict_maxnet(env.vars = env.vars.future, SDM.obj = SDM.obj,
      #                                                             model = model, taxon.name = pred.name,
      #                                                             clamp.env = clamp.env),
      #                               "maxent.jar" = predict_maxentJar(env.vars = env.vars.future, SDM.obj = SDM.obj,
      #                                                                model = model, taxon.name = pred.name),
      #                               stop("Unknown algorithm: ", SDM.obj@algorithm)
      #         )
      #         # Use predict function matching algorithm; filename accounts for multiple future scenarios
      #
      #         message("Writing future prediction raster to file...")
      #
      #         output.path <- paste0(data.dir, "futurePredictions/futurePredictionsRast/",
      #                               pred.name, ".tif")
      #         terra::writeRaster(future.pred, filename = output.path, overwrite = TRUE)
      #         # Write raster to file with model amendments
      #
      #         write_model_threshold(output.dir = paste0(data.dir, "futurePredictions/thresholdFuturePredictionsRast/"),
      #                               SDM.obj = SDM.obj,
      #                               model = model,
      #                               thresh.method = thresh.method,
      #                               future.pred = future.pred,
      #                               filename = pred.name)
      #         # Threshold output model w/ filename amendments
      #       }
      #     }
      #   }
      # }
      #
    }

    ## Stack prediction rasters ##


    ## Prediction Rasters

    message("Writing future prediction raster stack...")
    rast_files <- list.files(paste0(data.dir, "futurePredictions/futurePredictionsRast/"),
                             pattern = "\\.tif$",
                             full.names = TRUE)
    future.pred.stack <- terra::rast(rast_files)
    terra::writeRaster(future.pred.stack, filename = paste0(output.dir, "outputRasts/", "predictFutureRastStack.tif"), overwrite = T)
    # Stack all rasters into a singular stack for return to user

    # if (class(future.env) == "list") {
    #   future.pred.stack <- stack_futureEnvs(future.pred.stack)
    #   saveRDS(future.pred.stack, file = paste0(output.dir, "outputRasts/", "predictFutureRastStack.RDS"))
    #   # Create list structure object for ease of use later.
    #   ## Write pred raster stack and list
    # }
    ## TODO This is not working due to terra breaking file links on save

    ## Thresholded Rasters

    message("Writing threshold future prediction raster stack...")
    rast_files <- list.files(paste0(data.dir, "futurePredictions/thresholdFuturePredictionsRast/"),
                             pattern = "\\.tif$",
                             full.names = TRUE)
    future.pred.stack.thresh <- terra::rast(rast_files)
    terra::writeRaster(future.pred.stack.thresh, filename = paste0(output.dir, "outputRasts/", "predictThresholdFutureRastStack.tif"), overwrite = T)
    # Stack all rasters into a singular stack for return to user

    # if (class(future.env) == "list") {
    #   future.pred.stack.thresh <- stack_futureEnvs(future.pred.stack.thresh)
    #   saveRDS(future.pred.stack.thresh, file = paste0(output.dir, "outputRasts/", "predictThresholdFutureRastStack.RDS"))
    # }
    ## TODO This is not working due to terra breaking file links on save
    # Create list structure object for ease of use later.
    ## Write thresholded prediction raster stack and list


    # Placeholder Ensemble Predictions ----------------------------------------

    # TODO When refactoring, need to work in a workflow for creating future ensembles, and also
    # threshold the ensembles in a consistant format.


  }


###
##
#

# TODO NEEDS REFACTOR, attempts too much while interlinked, too locked onto the future + historical comparision model of chapter 1
# Model pipeline starting to become quite unwieldy.
# I have hacked on a solution to using raster.preds = F. However it's not built for long term anymore
# Need to refactor this at some point when time allows
# Main issue is that writing model thresholds requires both current and future environments.
# I would like it to work more independently so stuff will happy in various combinations of env, env.future. envReproject
# Also had issues with thresholding rasters with raster pred == F



## Null Pipeline ----------------------------------------------------------

#
##
###

#' Null Model Pipeline for Maxent Models
#' @description Generates and evaluates null models for Maxent models created with maxent_pipe()
#' @param data.dir Character, path to directory containing model outputs
#' @param output.dir Character, path to directory for null model outputs
#' @param null.iter Numeric, number of null model iterations to run
#' @param parallel Logical, whether to enable parallel processing
#' @param numCores Numeric, number of cores to use for parallel processing
#' @return None. Saves null model outputs to disk:
#'   - Null model RDS files
#'   - Evaluation histograms
#'   - Evaluation violin plots
#' @export

maxent_nullModel_pipe <- function(null.iter = 1000,
                                  data.dir = "data/maxent/",
                                  output.dir = "outputs/maxent/",
                                  # User parameters
                                  parallel = F,
                                  numCores = 4
                                  # Default parallel settings
) {

  sdm.outputs <- load_results(data.dir)
  # Load outputs into a list
  modelSelection <- readRDS(paste0(data.dir, "modelSelection.RDS"))
  # Load model selection choices from validation pipeline

  for (i in seq(length(sdm.outputs))) {

    SDM.obj <- sdm.outputs[[i]]
    # Select model iteration
    taxon.name <- SDM.obj@taxon.name
    # Extract model output name
    model <- modelSelection |>
      dplyr::filter(name == taxon.name) |>
      dplyr::select(modelSelection) |>
      as.character()
    # Selects best model from model selection table

    mod.null <- generate_nullModels(output.dir = paste0(data.dir, "nullModels/"),
                                    SDM.obj = SDM.obj,
                                    model = model,
                                    taxon.name = taxon.name,
                                    null.iter = null.iter,
                                    parallel = parallel,
                                    numCores = numCores,
                                    iter = i)
    # Function for generating a null for each model

    plot_nullModels(output.dir = paste0(output.dir, "nullModels/"),
                    mod.null = mod.null,
                    SDM.obj = SDM.obj,
                    taxon.name = taxon.name)
    # Function for writing eval graphs to file

    gc()
  }

}

###
##
#

# In Development ----------------------------------------------------------

## Results Pipeline --------------------------------------------------------

## TODO Results pipeline function
# Function to spit out any summary graphs that might be a commonly used output?
## Maps of current environments
## Maps of future environments Model + ensemble
## Beta diversity statistics

# IN DEV

#' Maxent Results Summary Pipeline (in development)
#' @description Not yet implemented. Placeholder for a pipeline to generate summary graphs and statistics (site maps, historical/future model maps, threshold overlaps) across a completed modelling run.
#' @param aoi sf Polygon, area of interest to crop summary maps to (distinct from `study.area`, which is the extent modelling was trained on)
#' @param study.area sf Polygon, study area polygon used for model training
#' @param pred.rast SpatRaster, historic model prediction raster(s)
#' @param pred.rast.future SpatRaster or list, future model prediction raster(s)
#' @param env.vars SpatRaster, historic environmental variables
#' @param env.vars.future SpatRaster or list, future environmental variables
#' @param data.dir Character, path to directory containing model outputs
#' @param output.dir Character, path to directory for summary outputs

results_pipe <- function(aoi, study.area,
                         pred.rast,
                         pred.rast.future,
                         env.vars,
                         env.vars.future,
                         data.dir = "data/maxent/",
                         output.dir = "outputs/maxent/"){
  # Crop by AOI -------------------------------------------------------------

  ## crop all maps to study area for mapping interests
  # Distinct from study area (which is where modelling was trained on)

  # Sitemap + Occurrence ----------------------------------------------------------------



  # Historical Model --------------------------------------------------------


  # Future Models -----------------------------------------------------------


  # Threshold Overlaps -------------------------------------------------------



}



## Env Data Pipeline -------------------------------------------------------

## TODO Extract environmental data functions into the toolkit
# Snow gum project and ACT gov modelling the most recent versions

## Reprojection Pipeline ---------------------------------------------------

## TODO Extracting these workflows from Future prediction pipeline to simplify processes
## I haven't really needed them after full run

#' Reproject Maxent Models to High Resolution (in development)
#' @description Incomplete/in development. Intended to reproject saved MaxEnt models onto a higher-resolution environmental raster and convert thresholded outputs to sf polygons. Currently references objects (e.g. `sdm.outputs`, `modelSelection`, `repredict.env`, `thresh.method`, `future.env`, `env.vars`, `env.vars.future`) that are not defined within the function's own arguments, so it will not run standalone in its current form.
#' @param simp.factor Numeric, simplification factor passed to `convert_thresh_toSf()` (default 2000)
#' @param clamp.env Logical, whether to clamp predictions to the range of the model's training data (default TRUE)

maxent_reprojModels <- function(simp.factor = 2000, clamp.env = T) {

  # Env Re-predictions -------------------------------------------------------

  if (repredict.env == T) {

    if(dir.exists(paste0(data.dir, "futurePredictions/envRePredictionsRast/")) == F){
      dir.create(paste0(data.dir, "futurePredictions/envRePredictionsRast/"), recursive = T)
    }
    # Create data folder

    message("Re-predicting hi-res environmental predictions...")

    for (i in seq(length(sdm.outputs))) {

      SDM.obj <- sdm.outputs[[i]]
      # Select model iteration
      taxon.name <- SDM.obj@taxon.name
      # Extract model output name
      model <- modelSelection |>
        dplyr::filter(name == taxon.name) |>
        dplyr::select(modelSelection) |>
        as.character()
      # Selects best model from model selection table

      message("Creating hi-res environmental predictions for: ", taxon.name)

      ## Pick model class

      if (SDM.obj@algorithm == "maxnet") {

        future.pred <- predict_maxnet(env.vars = env.vars,
                                      SDM.obj = SDM.obj,
                                      model = model,
                                      taxon.name = taxon.name,
                                      clamp.env = clamp.env)
        # Takes best model and outputs a prediction onto the future environment

      }

      if (SDM.obj@algorithm == "maxent.jar") {

        future.pred <- predict_maxentJar(env.vars = env.vars,
                                         SDM.obj = SDM.obj,
                                         model = model,
                                         taxon.name = taxon.name)

      }

      terra::writeRaster(x = future.pred, filename = paste0(data.dir, "futurePredictions/envRePredictionsRast/", taxon.name, ".tif"), overwrite = T)
      ## Write output reprojection

      ## Write thresholded model ##

      message("Writing hi-res environmental prediction raster to file: ", taxon.name)

      if (!is.null(env.vars.future)) {
        write_model_threshold(output.dir = paste0(data.dir, "futurePredictions/thresholdEnvRePredictionsRast/"),
                              SDM.obj = SDM.obj,
                              model = model,
                              thresh.method = thresh.method,
                              future.pred = future.pred)
        # Take prediction and write to file
      } else {
        message("Skipping creating model thresholds...")
        message("Disabled if raster.preds == F")
      }
    }

    ## Stack prediction rasters ##

    if(dir.exists(paste0(output.dir, "outputRasts/")) == F){
      dir.create(paste0(output.dir, "outputRasts/"), recursive = T)
    }
    # Create output raster directory

    message("Writing hi-res environmental prediction raster stack...")
    rast_files <- list.files(paste0(data.dir, "futurePredictions/envRePredictionsRast/"),
                             pattern = "\\.tif$",
                             full.names = TRUE)
    pred.stack <- terra::rast(rast_files)
    terra::writeRaster(pred.stack, filename = paste0(output.dir, "outputRasts/", "predictRastStack_hiRes.tif"), overwrite = T)
    # Stack all rasters into a singular stack for return to user

    if (!is.null(env.vars.future)) {
      message("Writing hi-res threshold environmental prediction raster stack...")
      rast_files <- list.files(paste0(data.dir, "futurePredictions/thresholdEnvRePredictionsRast/"),
                               pattern = "\\.tif$",
                               full.names = TRUE)
      pred.stack <- terra::rast(rast_files)
      terra::writeRaster(pred.stack, filename = paste0(output.dir, "outputRasts/", "predictThresholdRastStack_hiRes.tif"), overwrite = T)
      # Stack all rasters into a singular stack for return to user
    } else {
      message("Skipping writing model thresholds...")
      message("Disabled if raster.preds == F")
    }

  }
  ## Reprojects model outputs back onto high resolution rasters
  # Can be used to train model on sensible resolution, then project onto a higher resolution raster
  # I think this is reasonable, as sample size in original model should be large

  # Thresholds to SF Polygons --------------------------------------------------

  if (dir.exists("data/maxent/modelOutputRast/") == T) {

    ## Historical predictions
    message("Converting threshold rasters to sf polygons: ", "Historical Projections")
    convert_thresh_toSf(thresh.dir = paste0(data.dir, "futurePredictions/thresholdRast/"),
                       output.dir = paste0(output.dir, "outputSf/"),
                       output.name = "predictThreshSf.geojson",
                       simp.factor = simp.factor)

  } else {
    message("Raster predicitions present in SDM.obj, skipping threshold conversion to polygon...")
  }

  if (!is.null(env.vars.future)) {

    ## Future climate predictions
    message("Converting threshold rasters to sf polygons: ", "Future Projections")
    convert_thresh_toSf(thresh.dir = paste0(data.dir, "futurePredictions/thresholdFuturePredictionsRast/"),
                       output.dir = paste0(output.dir, "outputSf/"),
                       output.name = "predictThreshFutureSf.geojson",
                       simp.factor = simp.factor)

    if (is.list(future.env)) {
      sf.future.thresh <- sf::st_read(paste0(output.dir, "outputSf/predictThreshFutureSf.geojson"))
      sf.future.thresh <- split_futureEnvs_sf(sf.future.thresh)
      sf::st_write(sf.future.thresh, paste0(output.dir, "outputSf/predictThreshFutureSf.geojson"), delete_dsn = T)
    }
    # If multiple environments were predicted split out taxonName into time.period, ssp, model etc

  } else {
    message("Future environment not found, skipping threshold conversion to polygon...")
  }

  ## Re-projections
  if (!is.null(env.vars.future)) {
    if (repredict.env == T) {
      message("Converting threshold rasters to sf polygons: ", "Historical Re-projections")
      convert_thresh_toSf(thresh.dir = paste0(data.dir, "futurePredictions/thresholdEnvRePredictionsRast/"),
                         output.dir = paste0(output.dir, "outputSf/"),
                         output.name = "predictThreshSf_hiRes.geojson",
                         simp.factor = simp.factor)
    }
    ## Transform threshold prediction rasters into SF objects and writes to file
  } else {
    message("Skipping converting reprojection thresholds...")
    message("Disabled if raster.preds == F")
    message("Hack solution, need to fix this at some point")
  }


}
