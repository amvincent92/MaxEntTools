#' Extract Cross-Validated ROC Curves and AUC from ENMeval Results (parallel)
#' @param sdm.results List of ENMevaluation objects, as returned by load_results()
#' @param model.selection Data frame with columns `name` and `modelSelection`,
#'   giving the chosen tune.args per species (as produced in treeMaxentPostHoc.R)
#' @param output.dir Directory to save ROC plots to. Created if it does not
#'   already exist. Defaults to here::here("outputs", "maxent", "modelValidation", "rocCurves")
#' @param save.plots Logical, whether to save a ROC plot per species to output.dir.
#'   Defaults to TRUE
#' @param n.cores Number of PSOCK cluster workers to use. Defaults to 6
#' @return List with three elements:
#'   - roc: tibble of threshold/TPR/FPR per species and dataset (training/validation),
#'     for plotting ROC curves (pooled across folds)
#'   - auc: tibble of auc.val.avg (mean of per-fold validation AUC, matching
#'     ENMeval's own aggregation) and pooled training AUC per species
#'   - plots: named list of ggplot objects, one per species
#' @examples
#' model.selection <- readRDS("data/maxent/modelSelection.RDS")
#' cv.roc <- extract_cv_roc_parallel(sdm.results, model.selection)
#' 

## AI Generated Function, use with care ## 

# TODO This fixes computation bug that 

extract_cv_roc_parallel <- function(sdm.results,
                                     model.selection,
                                     output.dir = here::here("outputs", "maxent", "modelValidation", "rocCurves"),
                                     save.plots = TRUE,
                                     n.cores = 6) {

  # Version Check -------------------------------------------------------------

  tested.enmeval.version <- "2.0.5"
  installed.enmeval.version <- as.character(utils::packageVersion("ENMeval"))

  if (installed.enmeval.version != tested.enmeval.version) {
    cli::cli_alert_warning(
      "extract_cv_roc_parallel() was tested against ENMeval {.val {tested.enmeval.version}}, but {.val {installed.enmeval.version}} is installed - internal functions this relies on may have changed"
    )
  }
  # Relies on unexported ENMeval internals (clamp.vars(), enm.maxent.jar,
  # ENMevaluation slot structure) that could change between releases

  if (save.plots && !dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  # Per-Species Worker Function -------------------------------------------------

  extract_cv_roc_species <- function(y, SDM.obj, model.tune, output.dir, save.plots) {

    enm.mxjar <- ENMeval:::enm.maxent.jar
    # ENMdetails object holding maxent.jar's own fit/predict/args methods, so
    # refitting here matches what ENMeval itself would have done

    if (SDM.obj@algorithm != "maxent.jar") {
      cli::cli_alert_warning("Skipping {.val {y}}: algorithm {.val {SDM.obj@algorithm}} not supported")
      return(NULL)
    }
    # Only maxent.jar is supported

    if (SDM.obj@doClamp != TRUE) {
      cli::cli_abort("{.val {y}}: extract_cv_roc_parallel() assumes doClamp = TRUE, but this model was fit with doClamp = FALSE")
    }
    # Validation data is always clamped to the training range here,
    # replicating ENMeval's own behaviour

    tune.tbl.i <- SDM.obj@tune.settings |>
      dplyr::filter(tune.args == model.tune) |>
      dplyr::select(fc, rm)
    # Only the fc/rm combination ENMeval selected as best for this species

    other.settings <- SDM.obj@other.settings
    envs.names <- colnames(SDM.obj@occs)[3:ncol(SDM.obj@occs)]

    occs.env <- SDM.obj@occs[, envs.names]
    bg.env <- SDM.obj@bg[, envs.names]

    n.folds <- nlevels(SDM.obj@occs.grp)

    fold.occ.preds <- list()
    fold.bg.preds <- list()
    fold.occ.train.preds <- list()
    fold.bg.train.preds <- list()

    auc.val.k <- rep(NA_real_, n.folds)
    # Per-fold validation AUC, averaged below to match ENMeval's own
    # auc.val.avg aggregation (mean of per-fold AUC), rather than pooling
    # predictions across folds into a single AUC

    fold.error <- FALSE

    ## Per-Fold Cross-Validation ##

    for (k in seq_len(n.folds)) {

      occs.train.z <- occs.env[SDM.obj@occs.grp != k, , drop = FALSE]
      occs.val.z   <- occs.env[SDM.obj@occs.grp == k, , drop = FALSE]
      bg.train.z   <- bg.env[SDM.obj@bg.grp != k, , drop = FALSE]
      bg.val.z     <- bg.env[SDM.obj@bg.grp == k, , drop = FALSE]
      # occs.grp/bg.grp store which fold (1:n.folds) each point belongs to,
      # matching the split ENMeval assigned for this model

      val.z <- ENMeval:::clamp.vars(
        orig.vals = rbind(occs.val.z, bg.val.z),
        ref.vals = rbind(occs.train.z, bg.train.z),
        left = other.settings$clamp.directions$left,
        right = other.settings$clamp.directions$right,
        categoricals = other.settings$categoricals
      )
      occs.val.z <- val.z[seq_len(nrow(occs.val.z)), ]
      bg.val.z <- val.z[(nrow(occs.val.z) + 1):nrow(val.z), ]
      # Clamp held-out data to the training range, then re-split using
      # nrow(val.z) as the upper bound. ENMeval:::cv.enm() instead uses
      # nrow(bg.val.z) (the pre-clamp background count), which is wrong
      # whenever a fold's background count is <= its occurrence count - see
      # scripts/superseded/enmeval_github_issue.md for the full writeup

      mod.k.args <- enm.mxjar@args(occs.train.z, bg.train.z, tune.tbl.i, other.settings)

      mod.k <- tryCatch(
        do.call(enm.mxjar@fun, mod.k.args),
        error = function(cond) {
          cli::cli_alert_warning("Fold {.val {k}} failed for {.val {y}}: {cond$message}")
          NULL
        }
      )
      # Refit maxent.jar from scratch using only this fold's training data

      if (is.null(mod.k)) {
        fold.error <- TRUE
        next
      }

      occs.val.pred <- enm.mxjar@predict(mod.k, occs.val.z, other.settings)

      if (other.settings$validation.bg == "full") {
        bg.pred <- enm.mxjar@predict(mod.k, rbind(bg.train.z, bg.val.z), other.settings)
      } else {
        bg.pred <- enm.mxjar@predict(mod.k, bg.val.z, other.settings)
      }
      # Background prediction set follows validation.bg as set in the
      # original ENMevaluate() call: "full" pools the fold's training
      # background back in as the background sample

      occs.train.pred <- enm.mxjar@predict(mod.k, occs.train.z, other.settings)
      bg.train.pred <- enm.mxjar@predict(mod.k, bg.train.z, other.settings)
      # Also predict on this fold's own training data, giving a comparable
      # training-side curve once pooled across folds below

      eval.k <- predicts::pa_evaluate(p = occs.val.pred, a = bg.pred)
      auc.val.k[k] <- eval.k@stats$auc
      # Matches ENMeval:::tune.validate() exactly - one pa_evaluate() call
      # per fold, using the same validation.bg-dependent background set above

      fold.occ.preds[[k]] <- occs.val.pred
      fold.bg.preds[[k]] <- bg.pred
      fold.occ.train.preds[[k]] <- occs.train.pred
      fold.bg.train.preds[[k]] <- bg.train.pred
    }

    if (fold.error) {
      cli::cli_alert_warning("{.val {y}}: one or more folds failed; pooled curve uses remaining folds")
    }

    ## Pool Predictions and Compute ROC/AUC ##

    occ.preds.pooled <- unlist(fold.occ.preds)
    bg.preds.pooled <- unlist(fold.bg.preds)
    occ.train.preds.pooled <- unlist(fold.occ.train.preds)
    bg.train.preds.pooled <- unlist(fold.bg.train.preds)
    # Pooled across folds to draw one validation/training ROC curve per
    # species. The reported validation AUC below is NOT taken from this
    # pooled evaluation - it is the mean of the per-fold auc.val.k values
    # computed in the loop above, matching ENMeval's own auc.val.avg
    # aggregation. The pooled curve and the averaged AUC are therefore
    # drawn from different aggregations of the same underlying predictions

    eval.cv <- predicts::pa_evaluate(p = occ.preds.pooled, a = bg.preds.pooled)
    eval.train <- predicts::pa_evaluate(p = occ.train.preds.pooled, a = bg.train.preds.pooled)
    # eval.cv is used only for the pooled ROC curve plotted below, not for
    # the reported validation AUC (see auc.i)

    auc.val.avg <- mean(auc.val.k, na.rm = TRUE)
    # Directly comparable to ENMeval's own auc.val.avg column, differing
    # only in the corrected clamped-background slice above

    roc.val.i <- eval.cv@tr_stats |>
      dplyr::select(treshold, TPR, FPR) |>
      dplyr::mutate(species = y, dataset = "Validation")

    roc.train.i <- eval.train@tr_stats |>
      dplyr::select(treshold, TPR, FPR) |>
      dplyr::mutate(species = y, dataset = "Training")

    roc.i <- dplyr::bind_rows(roc.train.i, roc.val.i)

    auc.i <- tibble::tibble(
      species = y,
      auc.val.avg = auc.val.avg,
      auc.train = eval.train@stats$auc,
      n.folds = n.folds,
      n.folds.failed = sum(is.na(auc.val.k))
    )

    gg.plot <- ggplot2::ggplot(roc.i, ggplot2::aes(x = FPR, y = TPR, colour = dataset)) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
      ggplot2::geom_line() +
      ggplot2::scale_colour_viridis_d(name = NULL) +
      ggplot2::labs(
        title = y,
        subtitle = paste0(
          "Training AUC = ", round(eval.train@stats$auc, 3),
          ", Validation AUC (mean of folds) = ", round(auc.val.avg, 3)
        ),
        x = "False positive rate",
        y = "True positive rate"
      ) +
      ggplot2::theme_minimal()
    # Pooled training and validation ROC curves plotted together, for a
    # quick per-species visual check of overfitting/generalisation gap.
    # The validation AUC quoted in the subtitle is the per-fold average,
    # not the AUC of the pooled curve shown

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

# ENMeval only stores summary auc.val.avg per species; it doesn't keep the
# per-fold models or expose TPR/FPR needed to plot a ROC curve, or to see how
# variable cross-validated performance is within a species. This function
# re-does ENMeval's own cross-validation from scratch: for each species,
# refit a maxent.jar model on each fold's training data (using the same fold
# assignments and best tune.args ENMeval already chose), predict on that
# fold's held-out data, and evaluate each fold separately via pa_evaluate(),
# then average the per-fold AUCs - matching ENMeval:::tune.validate()'s own
# aggregation exactly, aside from the corrected clamped-background slice.
# Predictions are also pooled across folds to draw one cross-validated ROC
# curve per species, but that pooled curve's AUC is not what is reported in
# auc.val.avg. The training-data predictions are pooled the same way, giving
# a comparable "training" curve/AUC.
#
# Species are processed in parallel (one per foreach task) since refitting
# n.folds models per species is slow (~40s/species sequentially).
#
# Caveats: only algorithm == "maxent.jar" is supported; doClamp = TRUE is
# assumed and errors otherwise; the pooled validation AUC is conceptually
# different to ENMeval's auc.val.avg (see note at eval.cv above); and the
# "training" curve/AUC is fold-pooled (n.folds - 1 models per point), not a
# single full-data model fit.
