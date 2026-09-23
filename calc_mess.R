# MESS Plots --------------------------------------------------------------

#' Calculate a Multivariate Environmental Similarity Surface (MESS)
#'
#' @description Computes a MESS (Elith et al. 2010) comparing a projection
#' surface against the environmental conditions a species' MaxEnt model was
#' trained on (occurrence and background records). Negative values flag areas
#' where the projection surface is environmentally novel relative to training
#' data, i.e. where clamped predictions are most likely to be extrapolation
#' artefacts rather than genuine biological signal.
#'
#' @param SDM.obj ENMevaluation object for the taxon of interest. Training env
#'   values are taken from SDM.obj@occs and SDM.obj@bg.
#' @param future.rast SpatRaster to project the MESS onto (e.g. a future climate
#'   scenario, or env cropped to a different region). Layer names must match
#'   the env variable names used to fit SDM.obj.
#' @param taxon.name Character, name of taxon. Used for messaging and to name
#'   the output layer.
#' @param full Logical, if TRUE returns a SpatRaster with per-variable MESS
#'   layers plus the overall (minimum) MESS layer. Default FALSE returns only
#'   the overall MESS layer.
#' @param agg.fact Numeric, aggregation factor passed to terra::aggregate()
#'   applied to future.rast before the MESS calculation (e.g. 4 reduces
#'   resolution by a factor of 4 in each dimension, i.e. 16x fewer cells).
#'   Default NULL skips aggregation and runs at native resolution.
#'
#' @return SpatRaster of MESS values. Negative values indicate the pixel falls
#'   outside the training range for at least one variable (extrapolation);
#'   positive values indicate the pixel is within the training envelope for
#'   all variables.
#' @export
#'
#' @details Categorical layers (identified via terra::is.factor()) are
#'   excluded, since MESS as defined by Elith et al. (2010) applies to
#'   continuous predictors only. Run this alongside clamped MaxEnt
#'   predictions, not instead of them — clamping stops extrapolation, MESS
#'   tells you where and how much extrapolation would have occurred.
#'
#'   agg.fact trades resolution for speed. MESS is a smooth function of the
#'   underlying predictors, so aggregating (averaging) future.rast before the
#'   calculation should have limited effect on the broad-scale pattern of
#'   novelty, but it will blur fine-scale MESS boundaries and mask small,
#'   locally novel patches (e.g. narrow gullies, ridgelines). Aggregate is
#'   applied with fun = "mean", na.rm = TRUE, only to the continuous
#'   variables retained after categorical layers are dropped.

calc_mess <- function(SDM.obj, future.rast, taxon.name, full = FALSE, agg.fact = NULL) {

  message("Calculating MESS: ", taxon.name)

  cat.vars <- names(future.rast)[terra::is.factor(future.rast)]
  if (length(cat.vars) > 0) {
    message("Excluding categorical variable(s) from MESS: ", paste(cat.vars, collapse = ", "))
    future.rast <- future.rast[[base::setdiff(names(future.rast), cat.vars)]]
  }
  # MESS (Elith et al. 2010) is only defined for continuous predictors

  train.vars <- base::intersect(names(future.rast), names(SDM.obj@occs))
  missing.vars <- base::setdiff(names(future.rast), names(SDM.obj@occs))
  if (length(missing.vars) > 0) {
    message("Variable(s) in future.rast not found in SDM.obj training data, dropping: ", paste(missing.vars, collapse = ", "))
  }
  future.rast <- future.rast[[train.vars]]
  # Restrict projection surface to variables the model was actually trained on
  # base::intersect() preserves future.rast's variable order, kept consistent below

  if (!is.null(agg.fact)) {
    message("Aggregating future.rast by a factor of ", agg.fact, " before MESS calculation...")
    future.rast <- terra::aggregate(future.rast, fact = agg.fact, fun = "mean", na.rm = TRUE)
    # Aggregating only the continuous variables retained above
    # Trades resolution for speed, output MESS raster will match this coarser resolution
  }

  train.vals <- dplyr::bind_rows(SDM.obj@occs, SDM.obj@bg)[, train.vars]
  # Reference environment = occurrence + background records used in training
  # predicts::mess() matches reference columns to raster layers by POSITION, not name
  # so train.vars order must exactly match future.rast's layer order

  mess.rast <- predicts::mess(future.rast, train.vals, full = full)

  if (full == FALSE) {
    names(mess.rast) <- taxon.name
  }

  return(mess.rast)

}

## Function: single MESS map for one taxon x timeperiod x ssp ##

plot_messMap <- function(taxon.name, timeperiod, ssp, agg.fact = NULL, mask.poly = NULL) {

  # Calculates and plots a MESS surface for one taxon under one future
  # timeperiod x ssp combination, saving to
  # outputs/figures/mess/<timeperiod>/
  # mask.poly optionally restricts the input raster to a study area boundary
  # (e.g. ibra.alps, ibra.buffer) before MESS is calculated, accepting either
  # an sf/sfc object or a terra SpatVector. Masking before calc_mess() rather
  # than after reduces the extent/cell count going into aggregation and the
  # MESS calculation itself, rather than computing over the full raster and
  # discarding cells outside mask.poly afterwards.

  future.rast <- future.env[[timeperiod]][[ssp]]

  if (!is.null(mask.poly)) {

    if (inherits(mask.poly, "sf") || inherits(mask.poly, "sfc")) {
      mask.poly <- terra::vect(mask.poly)
    }
    # Coerce sf/sfc input to SpatVector for terra::mask()

    if (!terra::same.crs(future.rast, mask.poly)) {
      cli::cli_abort("mask.poly CRS does not match future.rast CRS. Reproject mask.poly before calling {.fn plot_messMap}.")
    }
    # Flag CRS mismatch rather than silently reprojecting, since reprojecting
    # a vector mask on the fly can mask misaligned inputs upstream

    future.rast <- future.rast |>
      terra::crop(mask.poly) |>
      terra::mask(mask.poly)
    # Crop first to reduce the extent before masking

  }

  mess.rast <- calc_mess(
    SDM.obj = sdm.results[[taxon.name]],
    future.rast = future.rast,
    taxon.name = taxon.name,
    full = FALSE,
    agg.fact = agg.fact
  )

  gg.plot <- ggplot() +
    geom_spatraster(data = mess.rast) +
    scale_fill_viridis_c(na.value = "transparent", name = "MESS") +
    ggtitle(paste0(taxon.name, ", ", timeperiod, ", ", ssp)) +
    theme_minimal()

  taxon.filename <- taxon.name |>
    stringr::str_replace_all("[[:space:]\\(\\)\\.]+", "_") |>
    stringr::str_replace_all("_+$", "")
  # Sanitise taxon name for use in a file path

  mess.dir <- paste0("outputs/mess/", taxon.filename, "/", timeperiod, "/")

  if (dir.exists(mess.dir) == FALSE) {
    dir.create(mess.dir, recursive = TRUE)
  }

  ggsave(filename = paste0(mess.dir, taxon.filename, "_", timeperiod, "_", ssp, ".jpeg"),
         plot = gg.plot, width = 7, height = 5, bg = "white")

  return(gg.plot)

}

## Loop: generate MESS maps for every taxon x timeperiod x ssp combination ##

generate_messMaps <- function(taxon.list = names(sdm.results),
                              timeperiods = names(future.env),
                              ssps = names(future.env[[1]]),
                              agg.fact = 10,
                              mask.poly = NULL) {

  # Generates MESS maps for every combination of taxon.list x timeperiods x
  # ssps, allowing a subset of taxa/time periods/ssps to be requested rather
  # than always plotting the full taxon x 3 x 3 combination set.

  cli::cli_progress_bar("Generating MESS maps",
                        total = length(taxon.list) * length(timeperiods) * length(ssps))

  for (taxon.name in taxon.list) {
    for (timeperiod in timeperiods) {
      for (ssp in ssps) {

        plot_messMap(taxon.name, timeperiod, ssp, agg.fact = agg.fact, mask.poly = mask.poly)

        cli::cli_progress_update()

      }
    }
  }

  cli::cli_progress_done()

}

generate_messMaps(taxon.list = names(alps)[str_detect(string = names(alps), pattern = "pauciflora")],
                  mask.poly = ibra.buffer)
# Generate just Snow gum maps
