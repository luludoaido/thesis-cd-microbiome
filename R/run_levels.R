# --- Run the analyses for genus and species ---
# Runs 05 for every cohort, and then 06 and 07, once per level.
for (level in c("genus", "species")) {
  for (cohort in c("franzosa", "hall", "hmp2")) {
    source("R/05_alpha_diversity.R")
    source("R/05_beta_diversity.R")
  }
  source("R/06_extremes.R")
  source("R/07_robustness.R")
}
