# --- Effect of the prevalence filter ---
# How many genera or species does the filter remove per cohort, how much abundance is lost,
# and are CD and controls affected differently? Only prints to the console.
library(tidyverse)
source("R/00_setup.R")

cohorts <- c("franzosa", "hmp2", "hall")

# relative abundances of one cohort, without taxa that occur in no sample
load_cohort <- function(cohort, level) {
  folder <- file.path(dir_processed, cohort)
  abund <- read_csv(file.path(folder, paste0(cohort, "_", level, ".csv")), show_col_types = FALSE)
  meta  <- read_csv(file.path(folder, paste0(cohort, "_meta.csv")), show_col_types = FALSE)
  X <- abund %>% filter(sample_id %in% meta$sample_id) %>% column_to_rownames("sample_id") %>% as.matrix()
  X <- X[rowSums(X) > 0, , drop = FALSE]
  X <- X / rowSums(X)
  list(X = X[, colSums(X > 0) > 0, drop = FALSE], meta = meta[match(rownames(X), meta$sample_id), ])
}

for (level in c("genus", "species")) {
  cat("\n\n#### Level:", level, "\n")
  data <- set_names(map(cohorts, load_cohort, level = level), cohorts)
  
  # what is removed at 5 and 10 %?
  overview <- imap_dfr(data, function(d, key) {
    X <- d$X
    map_dfr(c(0.05, 0.10), function(t) {
      keep <- passes_filter(X, key, level, t)
      lost <- rowSums(X[, !keep, drop = FALSE])
      tibble(cohort = key, thr = t, samples = nrow(X), taxa_all = ncol(X),
             taxa_kept = sum(keep), taxa_removed = sum(!keep),
             zeros_all = mean(X == 0), zeros_after = mean(X[, keep] == 0),
             lost_median = median(lost) * 100, lost_max = max(lost) * 100,
             taxa_per_sample_all = median(rowSums(X > 0)),
             taxa_per_sample_after = median(rowSums(X[, keep] > 0)))
    })
  })
  print(overview %>% mutate(across(where(is.numeric), ~ round(.x, 3))), width = Inf)
  
  # which samples lose the most at 10 %, and through which taxa?
  for (key in cohorts) {
    X <- data[[key]]$X
    meta <- data[[key]]$meta
    prev <- colMeans(X > 0)
    removed <- !passes_filter(X, key, level, 0.10)
    lost <- rowSums(X[, removed, drop = FALSE])
    cat("\n==", key, " samples losing >1%:", sum(lost > 0.01), " >5%:", sum(lost > 0.05), "\n")
    for (i in order(lost, decreasing = TRUE)[1:3]) {
      g <- X[i, removed]
      g <- sort(g[g > 0], decreasing = TRUE)[1:2]
      cat(rownames(X)[i], meta$diagnosis[i], round(lost[i] * 100, 1), "% :",
          paste0(names(g), " ", round(g * 100, 1), "% (prev ", round(prev[names(g)] * 100, 1), "%)", collapse = "; "), "\n")
    }
    # are CD and controls affected differently?
    dx <- meta$diagnosis
    cat("  >5% loss: CD", sum(lost > .05 & dx == "cd"), "/", sum(dx == "cd"),
        " control", sum(lost > .05 & dx == "control"), "/", sum(dx == "control"),
        " | people:", n_distinct(meta$subject_id[lost > .05]), "\n")
  }
}