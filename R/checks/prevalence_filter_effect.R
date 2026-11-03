# Effect of the prevalence filter per cohort
# uses loading + exclusions from R/07_extremes.R (incl. minimum read rule)
suppressMessages({library(tidyverse); library(readxl)})
source("R/00_setup.R")
src <- readLines("R/07_extremes.R")
# reuse loader + exclusions from the extremes script (everything up to prepare())
end <- grep("^# ---- 1. TSS", src) - 1
eval(parse(text = src[1:end]))
res <- map_dfr(c(franzosa="Franzosa", hmp2="HMP2", hall="Hall") %>% imap(~.y) %>% unlist() %>% set_names(c("Franzosa","HMP2","Hall")), function(key) {
  d <- load_cohort(key)
  X <- d$abund %>% filter(sample_id %in% d$meta$sample_id) %>% column_to_rownames("sample_id") %>% as.matrix()
  X <- X[rowSums(X) > 0, , drop = FALSE]; X <- X / rowSums(X)
  X <- X[, colSums(X > 0) > 0, drop = FALSE]
  prev <- colMeans(X > 0)
  out <- map_dfr(c(0.05, 0.10), function(t) {
    keep <- prev >= t
    lost <- rowSums(X[, !keep, drop = FALSE])
    tibble(thr = t, genera_kept = sum(keep), genera_removed = sum(!keep),
           zeros_after = mean(X[, keep] == 0),
           lost_median = median(lost) * 100, lost_max = max(lost) * 100,
           genera_per_sample_after = median(rowSums(X[, keep] > 0)))
  })
  out %>% mutate(cohort = key, samples = nrow(X), genera_all = ncol(X), zeros_all = mean(X == 0),
                 genera_per_sample_all = median(rowSums(X > 0)),
                 removed_prev_max = max(prev[prev < 0.10]) * 100)
})
print(res %>% mutate(across(where(is.numeric), ~ round(.x, 3))), width = 200)
for (key in c("franzosa","hmp2","hall")) {
  d <- load_cohort(key)
  X <- d$abund %>% filter(sample_id %in% d$meta$sample_id) %>% column_to_rownames("sample_id") %>% as.matrix()
  X <- X[rowSums(X) > 0, , drop = FALSE]; X <- X / rowSums(X); X <- X[, colSums(X > 0) > 0]
  prev <- colMeans(X > 0); rm_ <- prev < 0.10
  lost <- rowSums(X[, rm_, drop = FALSE])
  cat("\n==", key, " samples losing >1%:", sum(lost > 0.01), " >5%:", sum(lost > 0.05), "\n")
  top <- order(lost, decreasing = TRUE)[1:3]
  for (i in top) { g <- X[i, rm_]; g <- sort(g[g > 0], decreasing = TRUE)[1:2]
    cat(rownames(X)[i], d$meta$diagnosis[d$meta$sample_id == rownames(X)[i]], round(lost[i]*100,1), "% :", paste0(names(g), " ", round(g*100,1), "% (prev ", round(prev[names(g)]*100,1), "%)", collapse = "; "), "\n") }
}
cat("\n\n#### >5% loss by diagnosis\n")
for (key in c("franzosa","hmp2","hall")) {
  d <- load_cohort(key)
  X <- d$abund %>% filter(sample_id %in% d$meta$sample_id) %>% column_to_rownames("sample_id") %>% as.matrix()
  X <- X[rowSums(X) > 0, , drop = FALSE]; X <- X / rowSums(X); X <- X[, colSums(X > 0) > 0]
  prev <- colMeans(X > 0); lost <- rowSums(X[, prev < 0.10, drop = FALSE])
  dx <- d$meta$diagnosis[match(rownames(X), d$meta$sample_id)]
  cat(key, ": >5% loss CD", sum(lost > .05 & dx == "cd"), "/", sum(dx == "cd"), " control", sum(lost > .05 & dx == "control"), "/", sum(dx == "control"),
      " | people:", n_distinct(d$meta$subject_id[match(rownames(X)[lost > .05], d$meta$sample_id)]), "\n")
}
