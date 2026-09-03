## Project Context

This is a R package to consolidate all written tools and pipelines for a MaxEnt species distribution modelling framework.

## User Expertise

The user is a competent R programmer familiar with tidyverse, sf, terra, ggplot2, and MaxEnt SDM workflows. Do not explain basic R syntax or package concepts unless asked. Focus on solving the problem, not on teaching fundamentals.

## Assistant Behaviour

The following instructions govern how Posit Assistant should interact with this project and the user.

- Don't be overly agreeable, I would like balanced and scientific assessments of arguments and their merits. 
- Prefer R tidyverse syntax over base R when possible.
- Try to keep output coding style to the input style, although make suggestions if you think bad syntax or form is being provided.
- Avoid altering the style of the code or the comments unless expressly directed.
- Avoid adding extra features not specified by the user, only highlight their potential utility.
- When refactoring or altering provided code, please confirm with user before creating a subfunction or helper function unless requested.
- Don't default to editing existing code files, either provide in chat window, or write to a new R file.
- Never suggest or call git commands under any circumstances. Version control should be implemented by user.
- Prefer explicit, readable code over compact or clever solutions. Prioritise clarity — verbose is fine.
- When writing spatial operations (sf, terra), always consider CRS consistency, raster resolution/extent alignment, and memory implications for large rasters. Flag potential issues proactively.
- This is a research project; statistical and methodological choices may appear in publications. Flag assumptions, limitations, or alternative approaches for any inferential or modelling decisions, rather than just providing the most convenient solution.
- Check the active R session variables before suggesting to reload or recompute data that may already be in memory.
- Provide code in the chat window by default. Only write to a new file if the code is a self-contained script or function file, or if explicitly asked.
- Prioritise native posit assist tools over powershell, python, and other external tools without checking and justifying the decision
- Use Australian English in all edits.
- Don't remove browser functions unless asked, just issue a reminder that it's there.
- After long series of activity, please offer pauses to allow user to /compact the token space.

---

## R Syntax & Style Guide

The following style conventions are derived from the project codebase. Apply them consistently when writing or editing any R code for this project.

### Naming Conventions

- **Variables**: dot.case — `aus.poly`, `sp.list`, `pred.rast`, `occ.sf.pts`. Use this by default for all R objects.
- **Functions**: snake_case with an action verb prefix — `rast_naZero()`, `check_envCorr()`, `subset_BGpointsByStudyArea()`. Simple or general-purpose functions may use camelCase when specificity is not needed (e.g. `createEnvStack()`).
- **Column names in data frames**: camelCase — `scientificName`, `decimalLongitude`, `growthForm`. This applies to user-created columns; preserve original casing for columns sourced from external packages or data (ALA, GBIF, etc.).
- **Loop/index variables**: Single letters — `i`, `y`, `z`, `n`.
- **Abbreviations**: Use domain-standard short forms: `sp.` (species), `rast` (raster), `pts` (points), `sf` (simple features), `bg` (background), `env` (environment/environmental), `occ` (occurrence).

### Assignment & Operators

- Always use `<-` for assignment. Never use `=` for assignment or `->`.
- Use `=` only for function arguments and named list elements.
- **Never use `T` or `F` as shorthand.** Always write `TRUE` and `FALSE` in full.

### Spacing & Indentation

- **2-space indentation** throughout. No tabs.
- Spaces around all binary operators (`<-`, `==`, `+`, `~`, etc.).
- When breaking function arguments across lines, align subsequent arguments with the opening parenthesis.

### Pipes

- **Prefer the native pipe `|>`** in all new and edited code. The codebase is being migrated from `%>%` to `|>`; convert `%>%` to `|>` in any code you write or touch unless it would break functionality (e.g. when `.` placeholder behaviour is required by magrittr).
- Place the pipe operator at the **end** of the line, not the start of the next.
- Each verb in a dplyr/tidyr chain goes on its own line, indented 2 spaces.

### Strings

- Always use **double quotes** `"`. Single quotes are not used.

### Packages & Namespacing

- Load packages with `library()` at the top of scripts. Never use `require()`.
- Never call `library()` inside a function body. Use explicit `package::function()` namespace calls instead, or `requireNamespace("pkg", quietly = TRUE)` to test for optional dependencies.
- Group `library()` calls roughly by function (spatial, tidyverse, visualisation).
- Use **explicit namespace calls** for functions where disambiguation matters or clarity is helpful: `dplyr::filter()`, `terra::rast()`, `stringr::str_detect()`, `sf::st_read()`.

### Comments & Section Headers

- Use a four-level hierarchy for section headers, formed with trailing dashes:
  ```r
  # Major Section ===========================================================
  # Section -----------------------------------------------------------------
  ## Subsection -------------------------------------------------------------
  ## Minor block ##
  ```
- Comments go **below** the relevant code or function, describing what it does and the reason for its implementation. Do not place comments inline at the end of a line.
- Comments placed **above** a function are reserved for roxygen2 documentation (`#'`) only.
- Use the `cli` package (not `message()`, `cat()`, or `print()`) for all user-facing output inside functions or loops. See the Functions section for details.
- Mark deferred work with `# TODO` or `## TODO`.

### Functions

- Opening brace on the same line as `function()`: `my_fn <- function(x) {`.
- Include a brief comment at the top of a function body explaining its purpose.
- Provide `return()` explicitly when the return value may not be obvious; otherwise implicit return (last expression) is fine.
- Use `cli` for all messaging and user feedback. Key functions:
  - `cli::cli_inform()` — neutral progress or status messages (replaces `message()`)
  - `cli::cli_alert_success()` — confirm a step completed successfully
  - `cli::cli_alert_warning()` — non-fatal warnings the user should know about
  - `cli::cli_abort()` — fatal errors that stop execution (replaces `stop()`)
  - `cli::cli_progress_bar()` / `cli::cli_progress_update()` / `cli::cli_progress_done()` — progress bars for loops
  - Use inline markup to format values: `{.val {x}}`, `{.field fieldName}`, `{.fn fn_name}`, `{.path path}`

### Control Flow

- Opening brace on the same line: `if (condition) {`.
- `else` on the same line as the closing brace: `} else {`.
- Prefer `dplyr::case_when()` over nested `if/else` for multi-condition assignments inside pipelines.
- Use `for (x in collection)` loops with `cli::cli_progress_bar()` / `cli::cli_progress_update()` to track progress in loops, or `cli::cli_inform()` for milestone messages in nested loops.

### Object Management

- Use `rm()` to remove large intermediate objects when they are no longer needed.
- Follow `rm()` with `gc()` after heavy `terra` raster operations, where external C-level memory may not be released automatically. For regular R objects, `rm()` alone is sufficient.

### File I/O

- Prefer `here::here("data", "file.tif")` over bare relative paths for file I/O. This ensures paths resolve correctly regardless of working directory.
- Check for file/directory existence before reading or writing:
  ```r
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)
  if (file.exists(here::here("data", "file.tif"))) { ... } else { ... }
  ```
- Use `overwrite = TRUE` (not `T`) in `writeRaster()` and similar calls.
- Prefer `saveRDS()` / `readRDS()` for R objects; `terra::writeRaster()` for rasters; `sf::st_write()` for vector data.

### ggplot2

- Use `ggplot()` (not piped from a data frame) and add layers with `+` at the **end** of each line.
- Default theme: `theme_minimal()`.
- Default continuous fill: `scale_fill_viridis_c()`.
- Use `tidyterra::geom_spatraster()` for raster layers.
- Use `ggpubr::ggarrange()` for multi-panel compositions.
- Assign plots to a variable (e.g. `gg.plot <-`) before saving with `ggsave()`.
