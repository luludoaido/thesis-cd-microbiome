# --- Run the analyses for genus and species ---
# Runs the 05 scripts for every cohort, then 06 and 07, once per level.
# The loop sets "level" and "cohort", which the scripts use.

for (level in c("genus", "species")) {
  for (cohort in c("franzosa", "hall", "hmp2")) {
    source("R/05_alpha_diversity.R")
    source("R/05_beta_diversity.R")
  }
  source("R/06_extremes.R")
  source("R/07_robustness.R")
}
