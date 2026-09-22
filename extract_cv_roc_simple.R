#' Extract Cross-Validated ROC Curves from ENMeval Results
#'
#' Simplified, non-parallel companion to extract_cv_roc_parallel.R. ENMeval
#' never exposes per-fold predictions, so a ROC curve still has to be built
#' by refitting each fold's model here - that part isn't affected by the
#' cv.enm() slicing bug in scripts/superseded/enmeval_github_issue.md, since
#' this function calls ENMeval:::clamp.vars() directly with the correct
#' slice rather than going through the buggy cv.enm(). What this script
#' drops (versus extract_cv_roc_parallel.R) is the per-fold auc.val.avg
#' bookkeeping kept only to compare against ENMeval's own buggy metric, and
#' the foreach/doSNOW parallel dispatch - intended as the leaner version to
#' keep in the pipeline once the upstream bug is fixed and that comparison
#' is no longer needed.
#'
#' @param sdm.results List of ENMevaluation objects, as returned by load_results()
#' @param model.selection Data frame with columns `name` and `modelSelection`,
#'   giving the chosen tune.args per species (as produced in treeMaxentPostHoc.R)
#' @param output.dir Directory to save ROC plots to. Created if it does not
#'   already exist. Defaults to here::here("outputs", "maxent", "modelValidation", "rocCurves")
#' @param save.plots Logical, whether to save a ROC plot per species to output.dir.
#'   Defaults to TRUE
#' @return List with three elements:
#'   - roc: tibble of threshold/TPR/FPR per species and dataset (training/validation)
#'   - auc: tibble of pooled training and cross-validated AUC per species
#'   - plots: named list of ggplot objects, one per species
#' @examples
#' model.selection <- readRDS("data/maxent/modelSelection.RDS")
#' cv.roc <- extract_cv_roc(sdm.results, model.selection)

extract_cv_roc <- function(sdm.results,
                            model.selection,
                            output.dir = here::here("outputs", "maxent", "modelValidation", "rocCurves"),
                            save.plots = TRUE) {

  if (save.plots && !dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  enm.mxjar <- ENMeval:::enm.maxent.jar
  # ENMdetails object holding maxent.jar's own fit/predict/args methods, so
  # refitting here matches what ENMeval itself would have done

  roc.list <- list()
  auc.list <- list()
  plot.list <- list()

  # Species Loop ------------------------------------------------------------

  for (y in names(sdm.results)) {

    SDM.obj <- sdm.results[[y]]

    if (SDM.obj@algorithm != "maxent.jar") {
      cli::cli_alert_warning("Skipping {.val {y}}: algorithm {.val {SDM.obj@algorithm}} not supported")
      next
    }

    model.tune <- model.selection |>
      dplyr::filter(name == y) |>
      dplyr::pull(modelSelection)

    tune.tbl.i <- SDM.obj@tune.settings |>
      dplyr::filter(tune.args == model.tune) |>
      dplyr::select(fc, rm)

    other.settings <- SDM.obj@other.settings
    envs.names <- colnames(SDM.obj@occs)[3:ncol(SDM.obj@occs)]

    occs.env <- SDM.obj@occs[, envs.names]
    bg.env <- SDM.obj@bg[, envs.names]

    n.folds <- nlevels(SDM.obj@occs.grp)

    occ.val.preds <- list()
    bg.val.preds <- list()
    occ.train.preds <- list()
    bg.train.preds <- list()

    ## Per-Fold Cross-Validation ##

    for (k in seq_len(n.folds)) {

      occs.train.z <- occs.env[SDM.obj@occs.grp != k, , drop = FALSE]
      occs.val.z <- occs.env[SDM.obj@occs.grp == k, , drop = FALSE]
      bg.train.z <- bg.env[SDM.obj@bg.grp != k, , drop = FALSE]
      bg.val.z <- bg.env[SDM.obj@bg.grp == k, , drop = FALSE]

      val.z <- ENMeval:::clamp.vars(
        orig.vals = rbind(occs.val.z, bg.val.z),
        ref.vals = rbind(occs.train.z, bg.train.z),
        left = other.settings$clamp.directions$left,
        right = other.settings$clamp.directions$right,
        categoricals = other.settings$categoricals
      )
      occs.val.z <- val.z[seq_len(nrow(occs.val.z)), ]
      bg.val.z <- val.z[(nrow(occs.val.z) + 1):nrow(val.z), ]
      # Re-split using nrow(val.z) as the upper bound - the correct slice,
      # regardless of whether cv.enm()'s own bug has been fixed upstream

      mod.k <- do.call(
        enm.mxjar@fun,
        enm.mxjar@args(occs.train.z, bg.train.z, tune.tbl.i, other.settings)
      )

      occ.val.preds[[k]] <- enm.mxjar@predict(mod.k, occs.val.z, other.settings)
      occ.train.preds[[k]] <- enm.mxjar@predict(mod.k, occs.train.z, other.settings)
      bg.train.preds[[k]] <- enm.mxjar@predict(mod.k, bg.train.z, other.settings)

      if (other.settings$validation.bg == "full") {
        bg.val.preds[[k]] <- enm.mxjar@predict(mod.k, rbind(bg.train.z, bg.val.z), other.settings)
      } else {
        bg.val.preds[[k]] <- enm.mxjar@predict(mod.k, bg.val.z, other.settings)
      }
    }

    ## Pool Predictions and Compute ROC/AUC ##

    eval.val <- predicts::pa_evaluate(p = unlist(occ.val.preds), a = unlist(bg.val.preds))
    eval.train <- predicts::pa_evaluate(p = unlist(occ.train.preds), a = unlist(bg.train.preds))

    roc.i <- dplyr::bind_rows(
      eval.train@tr_stats |>
        dplyr::select(treshold, TPR, FPR) |>
        dplyr::mutate(species = y, dataset = "Training"),
      eval.val@tr_stats |>
        dplyr::select(treshold, TPR, FPR) |>
        dplyr::mutate(species = y, dataset = "Validation")
    )

    gg.plot <- ggplot2::ggplot(roc.i, ggplot2::aes(x = FPR, y = TPR, colour = dataset)) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
      ggplot2::geom_line() +
      ggplot2::scale_colour_viridis_d(name = NULL) +
      ggplot2::labs(
        title = y,
        subtitle = paste0(
          "Training AUC = ", round(eval.train@stats$auc, 3),
          ", Validation AUC = ", round(eval.val@stats$auc, 3)
        ),
        x = "False positive rate",
        y = "True positive rate"
      ) +
      ggplot2::theme_minimal()

    if (save.plots) {
      ggplot2::ggsave(
        filename = file.path(output.dir, paste0(gsub(" ", "_", y), "_roc.jpeg")),
        plot = gg.plot,
        width = 6,
        height = 5
      )
    }

    roc.list[[y]] <- roc.i
    auc.list[[y]] <- tibble::tibble(
      species = y,
      auc.val = eval.val@stats$auc,
      auc.train = eval.train@stats$auc,
      n.folds = n.folds
    )
    plot.list[[y]] <- gg.plot

    cli::cli_alert_success("{.val {y}} done")
  }

  return(list(
    roc = dplyr::bind_rows(roc.list),
    auc = dplyr::bind_rows(auc.list),
    plots = plot.list
  ))
}
# No parallelisation or per-fold error handling is included, for readability -
# add foreach/doSNOW dispatch (see extract_cv_roc_parallel.R) if refitting
# n.folds x n.species models sequentially becomes too slow
