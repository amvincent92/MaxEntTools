

# Model Validation Functions --------------------------------------------------------

## Note that validation tools may be different between maxent.jar outputs and maxnet

#' Write Variable Contribution Plot to File
#' @description Creates and saves variable contribution plot for the best performing model
#' @param output.dir Character string specifying output directory path
#' @param SDM.obj ENMeval object containing model results
#' @param model Character string specifying the name of the best model settings within `SDM.obj`, as selected by `best_model()`
#' @param data.dir Character string specifying directory to write the intermediate variable importance CSV to
#' @return None (writes plot file to disk)
#' @details Creates a bar plot showing the percent contribution of each environmental variable to the best model selected by best_model()
#' @examples
#' \dontrun{
#' write_var_contrib(SDM.obj, model = "fc.LQ_rm.1",
#'                    output.dir = "path/to/output/", data.dir = "path/to/data/")
#' }

write_var_contrib <- function(SDM.obj, model, output.dir, data.dir) {

  if(dir.exists(output.dir) == F){
    dir.create(output.dir, recursive = T)
  }
  # Create output folder

  if(dir.exists(data.dir) == F){
    dir.create(data.dir, recursive = T)
  }
  # Create data folder

  taxon.name <- SDM.obj@taxon.name
  # Extract model output name

  if (SDM.obj@algorithm == "maxent.jar") {

    # plot(eval.models(SDM.obj)[[model]])
    # plotting Variable Contribution for the best model
    ## Basic plot of the variable importance

    message(paste0("Writing variable contribution graph to file: ", taxon.name))

    var.df <- SDM.obj@variable.importance[[model]]
    # Extract variable importance from object

    utils::write.csv(var.df, file = paste0(data.dir, taxon.name, ".csv"))
    # Write itermediate data to file

    var.gg <- ggplot2::ggplot(var.df, ggplot2::aes(x = stats::reorder(variable, percent.contribution), y = percent.contribution)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Variable", y = "Percent Contribution",
           title = paste0("Variable Contribution: "),
           subtitle = paste0(taxon.name)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
    ## Method for creating ggplot of var contribution

    ggplot2::ggsave(plot = var.gg, filename = paste0(output.dir, taxon.name, ".jpeg"))
    # Write plot to file

  }
  # Maxent.jar method

  if (SDM.obj@algorithm == "maxnet") {
    message("Variable contribution not availible for maxnet models, skipping...")
  }

}

## Model validation pipeline

#' Write Response Curves to File
#' @description Creates and saves response curves plot for the best performing model
#' @param output.dir Character string specifying output directory path
#' @param SDM.obj ENMeval object containing model results
#' @param model Character string specifying the name of the best model settings within `SDM.obj`, as selected by `best_model()`
#' @param data.dir Character string specifying directory to write the intermediate response curve CSV to
#' @return None (writes plot file to disk)
#' @details Creates response curves showing predicted suitability across the range of each environmental variable, using the best model selected by best_model()
#' @examples
#' \dontrun{
#' write_response_curves(SDM.obj, model = "fc.LQ_rm.1",
#'                        output.dir = "path/to/output/", data.dir = "path/to/data/")
#' }

write_response_curves <- function(SDM.obj, model, output.dir, data.dir) {

  if(dir.exists(output.dir) == F){
    dir.create(output.dir, recursive = T)
  }
  # Create output folder

  if(dir.exists(data.dir) == F){
    dir.create(data.dir, recursive = T)
  }
  # Create data folder

  taxon.name <- SDM.obj@taxon.name
  # Extract model output name

  best.model <- ENMeval::eval.models(SDM.obj)[[model]]
  # Extract model for best model

  ## Method Maxent.jar #

  message(paste0("Writing response curves to file: ", taxon.name))

  if (SDM.obj@algorithm == "maxent.jar") {

    # Extract data without plotting
    response_data <- predicts::partialResponse(best.model, plot = FALSE)
    # TODO This is super slow on large datasets, not sure of any workaround?
    # TODO Save results at least?

    response_df <- do.call(rbind, lapply(names(response_data), function(var) {
      data.frame(
        variable = var,
        value = response_data[[var]][, 1],
        prediction = response_data[[var]][, 2]
      )
    }))
    # response_df <- data.table::rbindlist(
    #   lapply(names(response_data), function(var) {
    #     data.frame(
    #       variable = var,
    #       value = response_data[[var]][, 1],
    #       prediction = response_data[[var]][, 2]
    #     )
    #   }),
    #   use.names = TRUE
    # )
    #
    # response_df <- bind_rows(
    #   map(names(response_data), ~data.frame(
    #     variable = .x,
    #     value = response_data[[.x]][, 1],
    #     prediction = response_data[[.x]][, 2]
    #   ))
    # )
    ## Other Methods suggest by claude

    utils::write.csv(response_df, file = paste0(data.dir, taxon.name, ".csv"))

    # Create custom ggplot
    response.gg <- ggplot2::ggplot(response_df, ggplot2::aes(x = value, y = prediction)) +
      ggplot2::geom_line(color = "black", linewidth = 1) +
      ggplot2::facet_wrap(~variable, scales = "free_x") +
      ggplot2::labs(title = paste0("Response Curves: ", taxon.name),
           subtitle = paste("Model:", model),
           x = "Variable Value",
           y = "Predicted Suitability") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5),
            strip.text = ggplot2::element_text(size = 10))
    ## Ggplot of response curves

    ggplot2::ggsave(plot = response.gg, filename = paste0(output.dir, taxon.name, ".jpeg"),
           height = 8, width = 8)
    # Save to file

    # predicts::partialResponse(e.mxjar@models[[opt.seq$tune.args]],
    #                           var = "bio5")
    # Vignette recommended example

  }

  if (SDM.obj@algorithm == "maxnet") {

    try({
      grDevices::jpeg(filename = paste0(output.dir, taxon.name, ".jpeg"), height = 800, width = 800, quality = 400)
      plot(best.model, type = "cloglog")
    })
    grDevices::dev.off()

    # maxnet::response.plot(best.model,
    #                    v = "bio1",
    #                    type = "cloglog",
    #                    plot = F)
    # Works for one plot, but not combined

    ## TODO Plotting having issues for MaxNet models
    # Could we use maxnet response method and create the grid myself?
    # Could set plot is == F, manuall extract the data, and then create

  }

  ## Method Testing

  # # Response curves of best model
  # predicts::partialResponse(SDM.obj@models[[model]])
  # title(main = "Response Curves for Best Model", outer = TRUE, line = -1)

  # dismo::response(x = eval.models(SDM.obj)[[model]])
  ## Can't get this to work for some reasons
  ## Supposedly the recommend way

}

## TODO
#' Write ROC Curves to File (stub)
#' @description Not yet implemented. Placeholder for writing ROC curves for both maxent.jar and maxnet model algorithms.
#' @param variables Placeholder argument, not yet used

write_roc_curves.dev  <- function(variables) {
  ## TODO Would still like some rocCurves for both model algorithm's if possible
}
## TODO In development
