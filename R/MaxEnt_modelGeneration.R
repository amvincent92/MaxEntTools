# Model Generation Functions ---------------------------------------------------------------


## Data Checks ---------------------------------------------------------------

#' Check Number of Occurrence Records per Species
#'
#' @description Checks if species have sufficient occurrence records for modeling
#'
#' @param occ.sf.pts sf object containing species occurrence points
#' @param check.pts numeric threshold for minimum number of records (default 250)
#' @param subset.col Character string specifying the column name used to subset by species (default "scientificName")
#'
#' @return No return value, called for side effects. Prints messages about species with insufficient records.
#'
#' @details Iterates through unique species names and checks if the number of occurrence points
#' meets the minimum threshold. Messages are printed for species with insufficient records.
#'
#' @examples
#' \dontrun{
#' # Requires an sf object of occurrence points (occ.sf.pts)
#' check_numRecords(occ.sf.pts, check.pts = 100)
#' }

check_numRecords <- function(occ.sf.pts, check.pts = 250, subset.col = "scientificName") {

  for (i in seq_along(unique(occ.sf.pts[[subset.col]]))) {

    species <- unique(occ.sf.pts[[subset.col]])[i]

    message("Checking number of points in: ", species)

    n.pts <- occ.sf.pts |>
      dplyr::filter(occ.sf.pts[[subset.col]] == species) |>
      nrow()
    # Counts n species in each species

    if (n.pts < check.pts) {
      message(unique(occ.sf.pts[[subset.col]])[i], " has less than ", check.pts, " records: ", n.pts)
      Sys.sleep(1)
    }
    # Return warning if less than 100 records

  }
}

#' Write Environmental Correlation Plot to File
#' @description Creates and saves correlation plot showing relationships between environmental variables
#' @param output.dir Character string specifying output directory path
#' @param env.vars SpatRaster object containing environmental predictor layers
#' @param occ.sf.pts sf points object containing species occurrence locations
#' @return None (writes plot file to disk)
#' @details Extracts environmental values at occurrence points and creates correlation plot showing
#'          Pearson correlations between all pairs of variables using ggcorrplot
#' @examples
#' \dontrun{
#' check_envCorr(output.dir = "path/to/output/", env.vars = env.vars, occ.sf.pts = occ.points)
#' }

check_envCorr <- function(output.dir, env.vars, occ.sf.pts) {

  if(dir.exists(output.dir) == F){
    dir.create(output.dir, recursive = T)
  }
  # Create output folder

  ## Variance Inflation Factor Test ##

  message("Calculating Variance Inflation Factor...")

  if (any(is.factor(env.vars))) {
    message("Skipping any categorical variables in VIF analysis...")
    envs.vif <- usdm::vifstep(env.vars[[which(!is.factor(env.vars))]])
    # envs.rem <- envs.vif@excluded
  } else {
    envs.vif <- usdm::vifstep(env.vars)
    # envs.rem <- envs.vif@excluded
  }

  # saveRDS(envs.vif, file = paste0(output.dir, "vif_test.Rdata"))
  save(envs.vif, file = paste0(output.dir, "vif_test.Rdata"))

  ## Env Correlation Plot ##

  message("Generating environmental correlations figure...")

  if (any(is.factor(env.vars))) {
    env.vars <- env.vars[[which(!is.factor(env.vars))]]
    message("Skipping any categorical variables in correlation plot...")
  }


  # Extract values from rasters at occurrence points
  corr.test <- terra::extract(env.vars, occ.sf.pts, ID = FALSE)

  # Convert to data frame as needed
  corr.test <- as.data.frame(corr.test)

  # Check for NA values
  any(is.na(corr.test))
  corr.test <- na.omit(corr.test)
  # Drops any NA variables from env_extract output
  # Note: NAs can occur when points fall outside raster extent or on no-data pixels

  # Generate correlation plot
  ggcorrplot::ggcorrplot(corr = stats::cor(corr.test),
                         method = "circle",
                         type = "lower",
                         p.mat = ggcorrplot::cor_pmat(corr.test),
                         tl.cex = 16,
                         # insig = "blank"
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 24),
      legend.title = ggplot2::element_text(size = 20),
      legend.text = ggplot2::element_text(size = 16),
      # axis.text = element_text(size = 12),
      # axis.title = element_text(size = 14)
    ) +
    ggplot2::labs(title = "Environmental Variable Correlations")

  ggplot2::ggsave(paste0(output.dir, "envCorrelations.jpeg"), height = 12, width = 12)

}
# Function for testing correlation between environmental variables.


## Crop by Study Area ------------------------------------------------------

#' Subset Occurrence Records by Study Area (stub)
#' @description Not yet implemented. Placeholder for subsetting occurrence records to the study area, if this is ever needed separately from `subset_BGpointsByStudyArea()` and `subset_EnvByStudyArea()`.
#' @param variables Placeholder argument, not yet used

subset_OccByStudyArea.dev <- function(variables) {

  ## TODO IF NEEDED?

}

#' Subset Background Areas by Study Area
#'
#' @description Intersects background area polygons with a defined study area polygon
#'
#' @param bg.areas sf object containing background area polygons
#' @param study.area sf object defining the study area boundary
#'
#' @return sf object containing the intersected background areas
#'
#' @details
#' This function:
#' - Intersects background polygons with study area boundary
#' - Reduces processing time by only keeping relevant background areas
#' - Suppresses intersection warnings
#'
#' @examples
#' \dontrun{
#' bg_subset <- subset_BGpointsByStudyArea(bg_areas, study_area)
#' }

subset_BGpointsByStudyArea <- function(bg.areas, study.area) {

  message("Intersecting background area data by defined study area...")
  bg.areas <- bg.areas[sf::st_intersects(bg.areas, study.area, sparse = F),] |>
    sf::st_intersection(study.area) |>
    suppressWarnings()
  # Intersects all background area polygons to only within defined study area
  # This reduces processing time by not intersecting any polygons well outside the occurrence data

  return(bg.areas)

}

#' Subset Environmental Data by Study Area
#'
#' @description Crops environmental raster data to the extent of a defined study area
#'
#' @param env.vars SpatRaster object containing environmental variables
#' @param study.area sf object defining the study area boundary
#'
#' @return SpatRaster cropped to study area extent
#'
#' @details This function reduces processing time by cropping environmental data
#' to only the relevant study area
#'
#' @examples
#' \dontrun{
#' env_subset <- subset_EnvByStudyArea(env_vars, study_area)
#' }

subset_EnvByStudyArea <- function(env.vars, study.area) {

  message("Cropping environmental data by defined study area...")
  env.vars <- terra::crop(env.vars, terra::ext(terra::vect(study.area$geometry)))
  # Crop environmental variables by defined study area to reduce procesesing time
  # If the study area (and therfore the occurrence records) are targeted to a specific area and data, no point processing whole raster

  return(env.vars)

}

## Thin Records ------------------------------------------------------------

#' Thin SF records (Raster version)
#'
#' Thin spatial points by gridding and selecting one record per cell
#'
#' @param data.sf sf object containing point data
#' @param subset.col Column name for grouping (e.g., species names)
#' @param ref.rast Grid cell resolution extracted from reference raster object
#' @param output.dir Directory path for output file
#' @return Thinned sf object
#' Thin spatial points by gridding and selecting one record per cell

thin_records_sf_rast <- function(data.sf,
                                 subset.col,
                                 ref.rast,
                                 output.dir = NULL) {

  message("Thinning records...")

  data.crs <- as.character(sf::st_crs(data.sf)[1])
  data.sf <- sf::st_transform(data.sf, crs = terra::crs(ref.rast))

  stopifnot(
    inherits(data.sf, "sf"),
    inherits(ref.rast, "SpatRaster"),
    is.character(subset.col),
    subset.col %in% names(data.sf)
  )

  n_initial <- nrow(data.sf)
  message(glue::glue("Original N records: {n_initial}"))

  data_vect <- terra::vect(data.sf)

  thinned.sf <- data.sf |>
    dplyr::mutate(
      cellID = terra::extract(
        ref.rast,
        data_vect,
        cells = TRUE
      )$cell
    ) |>
    dplyr::group_by(.data[[subset.col]], cellID) |>
    dplyr::slice_sample(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(-cellID)

  n_final <- nrow(thinned.sf)
  message(glue::glue("Thinned N records: {n_final}"))

  thinned.sf <- sf::st_transform(thinned.sf, crs = data.crs)

  if (!is.null(output.dir)) {
    saveRDS(thinned.sf, file = paste0(output.dir, "thinnedSf_byRast.rds"))
  }
  # Save output

  return(thinned.sf)
}
# Thin records by a reference raster resolution

#' Thin SF records (Resolution version)
#'
#' Thin spatial points by gridding and selecting one record per cell
#'
#' @param data.sf sf object containing point data
#' @param subset.col Column name for grouping (e.g., species names)
#' @param res Grid cell resolution in degrees (default 0.02)
#' @param output.dir Directory path for output file
#' @return Thinned sf object
#' Thin spatial points by gridding and selecting one record per cell

thin_records_sf_res <- function(data.sf,
                                subset.col,
                                res = 0.02,
                                output.dir = NULL) {

  message("Thinning records...")

  # Transform CRS to a coordinate system
  data.crs <- as.character(sf::st_crs(data.sf)[1])
  data.sf <- sf::st_transform(data.sf, crs = "EPSG:4326")
  ## Resolutions get weird in non-coordinate CRS's

  # Input validation
  stopifnot(
    inherits(data.sf, "sf"),
    is.character(subset.col),
    subset.col %in% names(data.sf)
  )

  # Log initial count
  n_initial <- nrow(data.sf)
  message(glue::glue("Original N records: {n_initial}"))

  # Create base raster grid and assign cell IDs
  # Get extent from input data
  data_vect <- terra::vect(data.sf)
  data_ext <- terra::ext(data_vect)

  # Create raster with proper extent
  base_grid <- terra::rast(
    extent = data_ext,
    resolution = res,
    crs = terra::crs(data.sf)
  )
  base_grid <- terra::setValues(base_grid, 1)
  # base_grid <- terra::init(base_grid, 1)

  # Extract cell IDs and thin points
  thinned.sf <- data.sf |>
    dplyr::mutate(
      cellID = terra::extract(
        base_grid,
        data_vect,  # Use already created vect object
        cells = TRUE
      )$cell
    ) |>
    dplyr::group_by(.data[[subset.col]], cellID) |>
    dplyr::slice_sample(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(-cellID)

  # Log results
  n_final <- nrow(thinned.sf)
  message(glue::glue("Thinned N records: {n_final}"))

  # Transform back to orignal CRS
  thinned.sf <- sf::st_transform(thinned.sf, crs = data.crs)

  if (!is.null(output.dir)) {
    saveRDS(thinned.sf, file = paste0(output.dir, "thinnedSf_byRes.rds"))
  }
  # Save output

  return(thinned.sf)
}
# Thin records by a defined resolution

## A version of this is also maintained in MadR package
# If changing please reconcile!

## Background Areas ---------------------------------------------------------

#' Define Background Area for Species Distribution Modeling
#'
#' @description Creates a background area by intersecting polygons with occurrence points,
#' buffering the points, and creating a masked region within the polygons where the species occurs.
#'
#' @param bg.areas sf object containing background area polygons (e.g., biomes)
#' @param occ.sf.pts sf object containing species occurrence points. Records should be thinned to reduce computation time
#' @param study.area sf object defining the study area boundary
#' @param subset.col Character, column name for subsetting data by species (default "scientificName")
#' @param output.dir Character, path to output directory for saving `bgAreas.rds`, or NULL to skip saving (default NULL)
#' @param bg.buffer.dist numeric buffer distance around occurrence points
#' @param simplify.poly logical whether to simplify the output polygon
#' @param simplify.tolerance numeric simplification tolerance passed to `sf::st_simplify()` when simplify.poly is TRUE (default 1000)
#'
#' @return sf polygon object of the background area
#'

define_background_area <- function(bg.areas,
                                   occ.sf.pts,
                                   study.area,
                                   subset.col = "scientificName",
                                   output.dir = NULL,
                                   bg.buffer.dist = NULL,
                                   simplify.poly = TRUE,
                                   simplify.tolerance = 1000) {

  message("Defining background areas of each point subset...")

  # Input validation
  if(!subset.col %in% names(occ.sf.pts)) {
    stop("subset.col not found in occurrence data")
  }

  # Get unique values to iterate over
  subset.values <- unique(occ.sf.pts[[subset.col]])

  # Initialize output list
  output.list <- list()

  # Iterate over each unique value
  for(i in seq(length(subset.values))) {

    subset <- subset.values[[i]]
    pts.subset <- occ.sf.pts[occ.sf.pts[[subset.col]] == subset,]
    # Subset points

    message(paste("Processing", subset))
    intersecting.polys <- bg.areas |>
      dplyr::filter(sf::st_intersects(geometry, pts.subset, sparse = FALSE) |>
                      apply(1, any))
    # Create logical intersection between points and polygons

    if (!is.null(bg.buffer.dist)) {
      # If buffer option is selected

      buffered.points <- pts.subset |>
        sf::st_transform("+proj=moll") |>
        sf::st_buffer(bg.buffer.dist) |>
        sf::st_union() |>
        sf::st_transform(sf::st_crs(pts.subset))
      # Create buffer around points if value provided

      background.area <- sf::st_intersection(buffered.points, intersecting.polys) |>
        sf::st_union()
      # Intersect this buffer with any intersecting polys

    } else {
      background.area <- intersecting.polys |>
        sf::st_union()
      # If no buffer, just provide intersecting polys
    }

    if(simplify.poly) {
      background.area <- background.area |>
        sf::st_simplify(dTolerance = simplify.tolerance, preserveTopology = TRUE)
    }
    # If simplify.poly is T, simplify the output polygon

    output.list[[i]] <- background.area
    names(output.list[[i]]) <- subset
    # Store result
  }

  output.df <- data.frame(
    # subset_value = names(output.list),
    # This broke when I fixed the loop above?
    geometry = do.call(c, output.list)
  ) |>
    sf::st_as_sf()
  # Convert list to dataframe

  output.df <- output.df |>
    dplyr::mutate(subset_value = subset.values) |>
    dplyr::relocate(subset_value)
  # This replicates above functionality.

  if (!is.null(output.dir)) {
    saveRDS(output.df, file = paste0(output.dir, "bgAreas.rds"))
  }
  # Save output

  return(output.df)
}

## Background Points -------------------------------------------------------
### Bias Layer --------------------------------------------------------------

#' Create a Bias Layer Raster for Species Distribution Modeling
#'
#' @description Creates a bias layer to account for sampling effort across the study area
#' using occurrence record density.
#'
#' @param data.sf sf object containing species occurrence points
#' @param env.vars SpatRaster of environmental variables
#' @param bg.poly sf polygon defining the background area
#' @param method Character string specifying the bias layer method (default "percent"). Currently only "percent" is implemented
#' @param output.dir Character string specifying output directory path to save the bias layer RDS. If NULL (default), the raster is not written to disk
#'
#' @return SpatRaster of sampling bias, with values between 0.1 and 1
#'
#' @details
#' Values are normalized to maximum density and 0.1 is added to ensure minimum sampling.
#' Areas outside the background polygon are masked to NA.
#'
#' @import terra sf
#'

bias_layer_rast <- function(data.sf, env.vars, bg.poly, method = "percent", output.dir = NULL){

  ## Density layer (Grid Cell method) ##
  # Creates bias layer for the sampling error across the entire dataset
  # This seems preferred as subsetting to each species is likely to bias the sample too heavily to the biology
  # Resulting in overfitting of the results.

  ## Data Checks ##

  if (!identical(sf::st_crs(data.sf), sf::st_crs(bg.poly))) {
    stop("CRS mismatch between occurrence data and background polygon")
  }

  ## Sampling Bias Process (Percentage) ##

  if (method == "percent") {

    occ <- terra::vect(data.sf)
    bg <- terra::vect(bg.poly)

    bias.layer <- terra::rasterize(occ, env.vars, fun = "sum")
    # Rasterise occurrence records count for each grid cells
    terra::values(bias.layer)[is.na(terra::values(bias.layer))] <- 0
    # Writes any NA's to any unsampled area within the study area.
    # We want this as 0's should equal valid area with no occurrence, rather than NA area which we want NA (i.e ocean)
    terra::values(bias.layer) <- terra::values(bias.layer) / max(terra::values(bias.layer), na.rm = TRUE) + 0.1
    # Calculate highest sampled value for whole raster
    # Add 0.1 to values to make sure they are sampled sometimes, but weighted to areas of higher sampling
    # Sets a 10 percent sampling threshold for area's without a occurrence point.
    bias.layer <- terra::mask(bias.layer, bg)
    # Mask this layer by the background polygon to re-insert NA's for anything we don't want the model sampling (i.e ocean)

    if (!is.null(output.dir)) {
      saveRDS(bias.layer, file = paste0(output.dir, "biasLayer.rds"))
    }

    return(bias.layer)

  }

  ## Kernal Density Estimate Method ##

  ## TODO Density layer based on Kernel density. Develop if needed. (Jarrod method)

}
# Creates bias sampling layer using raster intersection method.

### Absence Points  -------------------------------------------------------

#' Sample Background Points for Species Distribution Modeling
#'
#' @details
#' Samples background points within each background area polygon using the bias layer as weights.
#' Points are returned as MULTIPOINT geometries grouped by subset_value.
#'
#' @param bg.areas sf object containing background area polygons for each species
#' @param bias.layer SpatRaster with sampling bias weights (generated by bias_layer function)
#' @param env.vars SpatRaster of environmental variables
#' @param sample.size numeric number of background points to sample (default 10000)
#' @param output.dir Character string specifying output directory path to save the sampled points RDS. If NULL (default), the object is not written to disk
#'
#' @return sf object with background points as MULTIPOINT geometries

sample_bg_points <- function(bg.areas, bias.layer, env.vars, sample.size = 10000, output.dir = NULL) {

  message("Sampling background points...")

  bg.points.list <- list()
  # Create list to iterate over

  for(i in seq(nrow(bg.areas))) {

    bg.area <- bg.areas[i,]
    # Grab bg.area of loop

    bias.layer.mask <- terra::mask(bias.layer, terra::vect(bg.area))
    # Uses generated BG area to mask the bias layer for SDM models
    # Works on dataset as a whole. Currently doesn't appear sensible to do on a species by species level as seems biased

    bg.points <- terra::spatSample(
      x = bias.layer.mask,
      size = sample.size,
      method = "weights",
      as.points = TRUE,
      replace = TRUE
    )
    # Sample a pseudo-random points based on a weights method, uses bias layer as the weighted

    coords <- terra::crds(bg.points) |>
      as.data.frame()
    colnames(coords) <- c("X", "Y")
    # Coerce column names from terra version to match SF st_coordinates output (caps vs no caps)

    bg.points.list[[i]] <- data.frame(
      subset_value = bg.area$subset_value,
      coords
    )
    # Create list of output sample points

    message("Sampling: ", as.data.frame(bg.area)[1,1])
  }

  bg.points.df <- do.call(rbind, bg.points.list)
  # Covert list into dataframe using do.call

  bg.points.sf <- sf::st_as_sf(bg.points.df, coords = c("X", "Y"), crs = sf::st_crs(bg.areas))
  # Convert coordinate df to a sf object

  bg.points.sf <- bg.points.sf |>
    dplyr::group_by(subset_value) |>
    dplyr::summarise(geometry = sf::st_union(geometry)) |>
    dplyr::ungroup()
  # Union SF object into multipoint objects to reduce list size

  if (!is.null(output.dir)) {
    saveRDS(bg.points.sf, file = paste0(output.dir, "bgPoints.rds"))
  }

  return(bg.points.sf)
}
# Function for randomly sampling background points using bias layer mask and environmental rasters.
