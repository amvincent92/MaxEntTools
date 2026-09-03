
# Null Model Functions ----------------------------------------------------

#' Generate Null Models for Species Distribution Model Evaluation
#'
#' @description Creates null models to assess significance of SDM predictions by randomizing occurrence data
#'
#' @param output.dir Character string specifying path to output directory
#' @param SDM.obj ENMeval object containing model results
#' @param model Character string specifying the name of the best model settings within `SDM.obj`, as selected by `best_model()`
#' @param taxon.name Character string specifying the species/taxon name, used for output file naming and messages
#' @param iter Character or numeric value used to distinguish this null model run in the output filename (e.g. species iteration index)
#' @param null.iter Numeric value specifying number of null iterations (default 100)
#' @param parallel Logical indicating whether to use parallel processing (default FALSE)
#' @param numCores Numeric value specifying number of cores for parallel processing (default 4)
#'
#' @return ENMnulls object containing null model results
#'
#' @details
#' Uses best model parameters (feature classes and regularization multiplier) to generate
#' null models by randomly sampling points from the background area. The null distribution
#' can be used to assess statistical significance of the actual model.
#'
#' @examples
#' \dontrun{
#' nulls <- generate_nullModels("output/nulls/", sdm_result, null.iter = 100)
#' }

generate_nullModels <- function(output.dir, SDM.obj, model, taxon.name, iter, null.iter = 100, parallel = F, numCores = 4){

  if(dir.exists(output.dir) == F){
    dir.create(output.dir, recursive = T)
  }
  # Create output folder

  best.model <- ENMeval::eval.models(SDM.obj)[[model]]
  # Extract model for best model

  fc <- stringr::str_extract(as.character(model), "(?<=\\.).*?(?=_)")
  rm <- stringr::str_extract(as.character(model), "(?<=rm\\.)\\d+(\\.\\d+)?") |>
    as.numeric()
  # Extract parameters from best model

  message("Simulating null model for: ", taxon.name)
  mod.null <- ENMeval::ENMnulls(SDM.obj, mod.settings = list(fc = fc, rm = rm), no.iter = null.iter,
                       parallel = parallel, numCores = numCores)
  # Generate null based on best fit model

  saveRDS(mod.null, file = paste0(output.dir, taxon.name, "_", iter, ".RDS"))
  # Save output null RDS

  return(mod.null)

}

#' Write Null Model Evaluation Plots to File
#'
#' @description Creates and saves histogram and violin plots comparing actual model
#' performance metrics to null distribution
#'
#' @param output.dir Character string specifying output directory path
#' @param mod.null ENMnulls object containing null model results
#' @param SDM.obj ENMeval object containing model results
#' @param taxon.name Character string specifying the species/taxon name, used for output file naming and messages
#'
#' @return None (writes plot files to disk)
#'
#' @details
#' Generates two types of visualization plots comparing the actual model performance
#' metrics to the null distribution:
#' - Histogram plots showing frequency distribution of null metrics
#' - Violin plots showing density distribution of null metrics
#'
#' Currently evaluates:
#' - 10th percentile training presence threshold (or.10p)
#' - Continuous Boyce Index (cbi.val)
#'
#' @examples
#' \dontrun{
#' plot_nullModels("output/nulls/", null_results, sdm_results)
#' }

plot_nullModels <- function(output.dir, mod.null, SDM.obj, taxon.name) {

  ## Histogram

  if(dir.exists(paste0(output.dir, "evalHistogram/")) == F){
    dir.create(paste0(output.dir, "evalHistogram/"), recursive = T)
  }
  # Create output folder

  try(output.plot <- ENMeval::evalplot.nulls(mod.null,
                                    # stats = c("auc.val", "auc.diff", "cbi.val", "or.mtp", "or.10p"),
                                    stats = c("or.10p", "cbi.val", "auc.val"),
                                    plot.type = "histogram",
                                    return.tbl = F))

  if (exists("output.plot")) {
    ggplot2::ggsave(plot = output.plot, filename = paste0(output.dir, "evalHistogram/", taxon.name, ".jpeg"),
           height = 8, width = 6)
  } else {
    message(paste0("Unable to plot null outputs for: ", taxon.name))
  }

  ## Violin Plots

  if(dir.exists(paste0(output.dir, "evalViolin/")) == F){
    dir.create(paste0(output.dir, "evalViolin/"), recursive = T)
  }
  # Create output folder

  try(output.plot <- ENMeval::evalplot.nulls(mod.null,
                                    # stats = c("auc.val", "auc.diff", "cbi.val", "or.mtp", "or.10p"),
                                    stats = c("or.10p", "cbi.val", "auc.val"),
                                    plot.type = "violin",
                                    return.tbl = F))

  if (exists("output.plot")) {
    ggplot2::ggsave(plot = output.plot, filename = paste0(output.dir, "evalViolin/", taxon.name, ".jpeg"),
           height = 8, width = 6)
  } else {
    message(paste0("Unable to plot null outputs for: ", taxon.name))
  }

  # Save output null summary graphs
  # Return table = T will return raw data if you want to custom figures

}
