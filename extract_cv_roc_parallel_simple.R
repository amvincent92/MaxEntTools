#' Extract Cross-Validated ROC Curves from ENMeval Results (parallel, no bug workaround)
#'
#' Parallelised version of extract_cv_roc.R, for once the ENMeval::cv.enm()
#' clamped-background slicing bug described in
#' scripts/superseded/enmeval_github_issue.md is fixed upstream. This still
#' has to refit each fold's model, since ENMeval never exposes per-fold
#' predictions - but it drops the per-fold auc.val.avg bookkeeping that
#' extract_cv_roc_parallel.R keeps only to compare against ENMeval's own
#' buggy metric.
#'
#' @param sdm.results List of ENMevaluation objects, as returned by load_results()
#' @param model.selection Data frame with columns `name` and `modelSelection`,
#'   giving the chosen tune.args per species (as produced in treeMaxentPostHoc.R)
#' @param output.dir Directory to save ROC plots to. Created if it does not
#'   already exist. Defaults to here::here("outputs", "maxent", "modelValidation", "rocCurves")
#' @param save.plots Logical, whether to save a ROC plot per species to output.dir.
#'   Defaults to TRUE
#' @param n.cores Number of PSOCK cluster workers to use. Defaults to 6
#' @return List with three elements:
#'   - roc: tibble of threshold/TPR/FPR per species and dataset (training/validation)
#'   - auc: tibble of pooled training and cross-validated AUC per species
#'   - plots: named list of ggplot objects, one per species
#' @examples
#' model.selection <- readRDS("data/maxent/modelSelection.RDS")
#' cv.roc <- extract_cv_roc_parallel_simple(sdm.results, model.selection)

extract_cv_roc_parallel_simple <- function(sdm.results,
                                            model.selection,
                                            output.dir = here::here("outputs", "maxent", "modelValidation", "rocCurves"),
                                            save.plots = TRUE,
                                            n.cores = 6) {

  if (save.plots && !dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  # Per-Species Worker Function -----------------------------------------------

  extract_cv_roc_species <- function(y, SDM.obj, model.tune, output.dir, save.plots) {

    enm.mxjar <- ENMeval:::enm.maxent.jar
    # ENMdetails object holding maxent.jar's own fit/predict/args methods, so
    # refitting here matches what ENMeval itself would have done

    if (SDM.obj@algorithm != "maxent.jar") {
      cli::cli_alert_warning("Skipping {.val {y}}: algorithm {.val {SDM.obj@algorithm}} not supported")
      return(NULL)
    }

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

    auc.i <- tibble::tibble(
      species = y,
      auc.val = eval.val@stats$auc,
      auc.train = eval.train@stats$auc,
      n.folds = n.folds
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

    return(list(roc = roc.i, auc = auc.i, plot = gg.plot))
  }
  # Per-species worker, dispatched to cluster workers via foreach %dopar%
  # below. Takes only this species' own ENMevaluation object and tune.args
  # string (not the full sdm.results list/model.selection table) to keep
  # what's serialised to each worker small

  # Parallel Dispatch -----------------------------------------------------------

  model.tune.vec <- purrr::map_chr(
    names(sdm.results),
    ~ model.selection |>
      dplyr::filter(name == .x) |>
      dplyr::pull(modelSelection)
  )
  # Looked up once, sequentially (cheap), so each worker receives only its
  # own scalar tune.args string

  `%dopar%` <- foreach::`%dopar%`
  # Explicit namespace assignment rather than library(foreach)

  cl <- parallel::makeCluster(n.cores, type = "PSOCK")
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  # doSNOW (rather than doParallel) is registered so .options.snow can drive
  # a live progress bar below - doParallel gives no cross-worker feedback.
  # PSOCK is used since forking (e.g. doMC) isn't available on Windows

  cli::cli_inform("Computing cross-validated ROC curves for {.val {length(sdm.results)}} species across {.val {n.cores}} workers")

  pb <- txtProgressBar(max = length(sdm.results), style = 3)
  progress.fun <- function(n) setTxtProgressBar(pb, n)
  snow.opts <- list(progress = progress.fun)
  on.exit(close(pb), add = TRUE)

  results.list <- foreach::foreach(
    y = names(sdm.results),
    SDM.obj = sdm.results,
    model.tune = model.tune.vec,
    .packages = c("ENMeval", "predicts", "dplyr", "ggplot2", "tibble", "cli"),
    .errorhandling = "pass",
    .options.snow = snow.opts
  ) %dopar% {
    extract_cv_roc_species(y, SDM.obj, model.tune, output.dir, save.plots)
  }
  names(results.list) <- names(sdm.results)
  # Iterates sdm.results, its names, and the pre-computed tune.args together,
  # so each worker only receives the single ENMevaluation object/tune string
  # it needs. .errorhandling = "pass" means one species failing doesn't halt
  # the other workers

  # Combine Results -----------------------------------------------------------

  is.failed <- purrr::map_lgl(results.list, ~ inherits(.x, "error") || inherits(.x, "condition"))

  if (any(is.failed)) {
    cli::cli_alert_warning("{.val {sum(is.failed)}} species failed: {.val {names(results.list)[is.failed]}}")
  }

  is.null.result <- purrr::map_lgl(results.list, is.null)

  results.ok <- purrr::discard(results.list, is.failed | is.null.result)
  # Drops failed species (warned above) and skipped species (e.g. non-
  # maxent.jar, which return NULL) before combining results

  roc.table <- purrr::map(results.ok, "roc") |>
    dplyr::bind_rows()

  auc.table <- purrr::map(results.ok, "auc") |>
    dplyr::bind_rows()

  plot.list <- purrr::map(results.ok, "plot")

  cli::cli_inform("Computed cross-validated ROC curves for {.val {nrow(auc.table)}} species")

  return(list(roc = roc.table, auc = auc.table, plots = plot.list))
}
