## Instructions for setting up Java to get maxent working ##

## TODO Still haven't solved this bug, will revisit later
# Occasionally works but a bit unreliable

# Built on ENMeval 2.0.5, maxent.jar v3.4.4, dismo 1.3-14
# rJava and Java versions are 1.0-11, and Version 24) respectively.

# Marcel had success running Java version 21.0.2

# Java Setup --------------------------------------------------------------

# https://cran.r-project.org/web//packages/rJavaEnv/vignettes/rJavaEnv.html

# install.packages(rJavaEnv)
# library(rJavaEnv)

# .jinit()
# .jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
## Check what java is getting read from system


# rJavaEnv::java_quick_install(version = 21)
# Works to set a unique java environment
# rJavaEnv::java_quick_install(version = 8)
# rJavaEnv::java_check_version_rjava()
## Tools for setting Java versions

# java_env_unset()
# java_clear("project", delete_all = TRUE)

# Update Dismo MaxEnt -----------------------------------------------------

# system.file("java", package = "dismo")

# https://biodiversityinformatics.amnh.org/open_source/maxent/
# Should update the dismo version of maxent.jar at this path to lastest version from

## Java Config ------------------------------------------------------------

# devtools::install_github("jamiemkass/ENMeval")
# Latest version of ENMeval from github.
# install.packages("ENMeval")
# Or CRAN

#' Configure MaxEnt Environment Settings
#'
#' @description Sets up the R environment for running MaxEnt models by configuring Java memory and loading required packages
#'
#' @param seed Numeric value for random seed to ensure reproducibility
#' @param java.mem.gb Numeric value specifying Java memory allocation in GB (default 8)
#'
#' @param algorithm Character string specifying the MaxEnt algorithm to configure for, either "maxent.jar" or "maxnet". Only "maxent.jar" requires the Java environment reset
#'
#' @return No return value, called for side effects
#'
#' @details This function:
#' - Loads required packages
#' - Configures Java memory settings
#' - Sets random seed for reproducibility
#'
#' @note Java must be installed locally to run MaxEnt. See java_setup.R for more details
#' @note This function intentionally uses `library()` to attach packages (rather than `::` or
#' `requireNamespace()`), since options(java.parameters = ...) only takes effect for rJava if set
#' before rJava is attached, and detach()/library() is the mechanism used here to force a clean
#' reload of the Java stack. This is a deliberate exception to the project's no-library()-in-functions
#' convention; see AGENTS.md discussion.

## DEPRECIATED ##
# Attempting new fix in R/env_config()

env_config <- function(algorithm, seed = 343, java.mem.gb = 8){

  ## Java Config

  if (algorithm == "maxent.jar") {

    message("Initalising java environment...")
    # Detach packages if loaded
    if("rJava" %in% (.packages())) detach("package:rJava", unload=TRUE)
    if("dismo" %in% (.packages())) detach("package:dismo", unload=TRUE)
    if("ENMeval" %in% (.packages())) detach("package:ENMeval", unload=TRUE)
    if("rJavaEnv" %in% (.packages())) detach("package:rJavaEnv", unload=TRUE)

    message("Setting java memory capacity to ", java.mem.gb, "gb...")
    options(java.parameters = paste0("-Xmx", java.mem.gb, "g"))
    # Detaching and changing java.parameters makes more memory available in rJava calls
    # ENMevaluate failing at larger raster resolutions when using algorith == "maxent.jar"
    # Order is very important, if rJava or dismo is loaded the options() call won't work.

    # Restart clean
    library(rJavaEnv)
    library(rJava)
    library(dismo)
    library(ENMeval)

  }
  ## Java needs to be installed locally to run MaxEnt
  ## Maxnet parameterisation might get around this?

  ## Dependencies
  message("Loading dependancies...")

  library(terra)
  library(sf)
  library(ENMeval)
  library(dismo)

  message("Setting seed for reproducibility...")
  set.seed(seed)
  # Set seed for reproducibility

}
