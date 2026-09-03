# Java Config (staging replacement for env_config() in MaxEnt_modelGeneration.R) -----------

# TODO: this is a staging replacement for env_config() in MaxEnt_modelGeneration.R.
# Once validated across a full pipeline run (both algorithm = "maxent.jar" and "maxnet",
# including repeated calls within a single session), move this into
# MaxEnt_modelGeneration.R in place of the original and delete the original.

# Session-level state tracking whether the JVM has already been initialised.
# rJava embeds the JVM via JNI; once a JVM has started inside the R process it
# cannot be safely torn down and reinitialised. Calling detach()/unloadNamespace()
# on rJava after a live JVM exists will crash R. This flag stops env_config() below
# from repeating the unload/reload sequence more than once per session, so it can
# be safely called at the start of every pipeline function without needing a
# fresh R restart in between.
.java.state <- new.env(parent = emptyenv())
.java.state$initialised <- FALSE
.java.state$java.mem.gb <- NA_real_

#' Configure MaxEnt Environment Settings (staging replacement)
#'
#' @description Sets up the R environment for running MaxEnt models by configuring Java memory and loading required packages
#'
#' @param seed Numeric value for random seed to ensure reproducibility
#' @param java.mem.gb Numeric value specifying Java memory allocation in GB (default 8)
#' @param algorithm Character string specifying the MaxEnt algorithm to configure for, either "maxent.jar" or "maxnet". Only "maxent.jar" requires the Java environment reset
#' @param force.reset Logical. If `TRUE`, forces the detach/reload sequence to run again even if the JVM has already been initialised this session. This is unsafe once a JVM is live and should only be used if you have independently confirmed rJava is not currently attached (e.g. after manually detaching it yourself). Default `FALSE`.
#'
#' @return No return value, called for side effects
#'
#' @details This function:
#' - Loads required namespaces (without attaching them to the search path)
#' - Configures Java memory settings, once per session
#' - Sets random seed for reproducibility
#'
#' @note Java must be installed locally to run MaxEnt. See java_setup.R for more details
#' @note This is a staging replacement for the original `env_config()` in
#' MaxEnt_modelGeneration.R. Unlike the original, it uses `requireNamespace()`/`unloadNamespace()`
#' instead of `library()`/`detach()`, since all downstream calls to rJava, dismo and ENMeval in
#' this package are already namespace-qualified (`dismo::`, `ENMeval::`) and do not require these
#' packages to be attached to the search path. A session-level guard (`.java.state`) prevents the
#' unload/reload sequence from running more than once per session, because rJava's embedded JVM
#' cannot be safely reinitialised once started.
#'
#' @note Caveat: `.java.state` is defined at package top-level and is re-created (reset to
#' uninitialised) whenever the package namespace is reloaded, e.g. via `devtools::load_all()`
#' during development. If the JVM was already started in a previous `load_all()` call within the
#' same R session, this guard will not detect that and a further `load_all()` followed by
#' `algorithm = "maxent.jar"` may still crash. This caveat does not apply to normal use of an
#' installed/attached package, only to iterative package development.

env_config <- function(algorithm, seed = 343, java.mem.gb = 8, force.reset = FALSE) {

  ## Java Config

  if (algorithm == "maxent.jar") {

    if (.java.state$initialised && !force.reset) {

      if (!identical(.java.state$java.mem.gb, java.mem.gb)) {
        cli::cli_alert_warning(
          "Java memory was already set to {.val {.java.state$java.mem.gb}}GB earlier this session.
          The JVM cannot be reinitialised, so the requested {.val {java.mem.gb}}GB will have no
          effect. Restart R to change {.arg java.mem.gb}, or call with {.arg force.reset = TRUE}
          if you have confirmed it is safe to do so."
        )
      } else {
        cli::cli_inform("Java environment already initialised this session, skipping reset...")
      }

    } else {

      if (force.reset) {
        cli::cli_alert_warning(
          "Forcing java environment reset ({.arg force.reset = TRUE}). This will crash R if a
          JVM is already running."
        )
      }

      cli::cli_inform("Initalising java environment...")
      # Unload namespaces if loaded. These are not attached via library(), so there is nothing
      # to detach from the search path, only the namespace itself needs unloading.
      if (isNamespaceLoaded("rJava")) unloadNamespace("rJava")
      if (isNamespaceLoaded("dismo")) unloadNamespace("dismo")
      if (isNamespaceLoaded("ENMeval")) unloadNamespace("ENMeval")
      if (isNamespaceLoaded("rJavaEnv")) unloadNamespace("rJavaEnv")

      cli::cli_inform("Setting java memory capacity to {java.mem.gb}gb...")
      options(java.parameters = paste0("-Xmx", java.mem.gb, "g"))
      # Order is very important: options() must be set before rJava's namespace is (re)loaded,
      # or the heap-size setting won't be picked up.

      requireNamespace("rJavaEnv", quietly = TRUE)
      requireNamespace("rJava", quietly = TRUE)
      requireNamespace("dismo", quietly = TRUE)
      requireNamespace("ENMeval", quietly = TRUE)

      .java.state$initialised <- TRUE
      .java.state$java.mem.gb <- java.mem.gb

    }

  }
  ## Java needs to be installed locally to run MaxEnt
  ## Maxnet parameterisation might get around this?

  ## Dependencies
  cli::cli_inform("Loading dependancies...")

  requireNamespace("terra", quietly = TRUE)
  requireNamespace("sf", quietly = TRUE)
  requireNamespace("ENMeval", quietly = TRUE)
  requireNamespace("dismo", quietly = TRUE)

  cli::cli_inform("Setting seed for reproducibility...")
  set.seed(seed)
  # Set seed for reproducibility

}
