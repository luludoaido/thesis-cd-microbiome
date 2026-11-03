library(tidyverse)
library(readxl)
source("R/00_setup.R")

# =====================================================================
# HMP2: match clinical metadata (HBI, calprotectin, antibiotics) to the
# metagenome samples of the analysis.
# - Values are in hmp2_metadata_2018-08-20.csv, often in rows of other
#   data types (metabolomics, viromics ...) from the same visit.
# - Matching: 1) directly via External ID, 2) otherwise same visit
#   (person + week_num), 3) otherwise nearest visit of the same person
#   within max_gap weeks.
# Output: results/hmp2_clinical/hmp2_clinical_per_sample.csv + summary
# =====================================================================

path     <- "data"
outdir   <- file.path(dir_results, "hmp2_clinical")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
max_gap  <- 2      # weeks

raw <- read.csv(file.path(path, "raw/hmp2/hmp2_metadata_2018-08-20.csv"), sep = ";", check.names = FALSE)

# samples of the analysis (same exclusions as 06_diversity.R, see R/00_setup.R)
meta <- read.csv(file.path(path, "processed/hmp2/hmp2_meta.csv"), check.names = FALSE) %>%
  transmute(sample_id = `External ID`, subject_id = `Participant ID`,
            diagnosis = tolower(diagnosis), week = week_num)
genus <- read_csv(file.path(path, "processed/hmp2/hmp2_genus_merged.csv"), show_col_types = FALSE)
reads <- raw %>% filter(data_type == "metagenomics") %>% select(sample_id = `External ID`, reads = reads_filtered)
meta  <- meta %>% filter(sample_id %in% genus$sample_id) %>% left_join(reads, by = "sample_id") %>%
  filter(is.na(reads) | reads >= min_reads, !sample_id %in% hmp2_empty_samples)
cat("Samples in the analysis:", nrow(meta), "\n")

# all visits with clinical data (any data type)
abx_col <- "Has the subject used any antibiotics since the last visit?"
clin <- raw %>%
  transmute(sample_id = `External ID`, subject_id = `Participant ID`, week = week_num,
            hbi = suppressWarnings(as.numeric(hbi)),
            calpro = suppressWarnings(as.numeric(fecalcal)),
            abx = na_if(trimws(.data[[abx_col]]), ""))

# match each variable separately so that e.g. HBI and calprotectin may come from different rows
match_var <- function(var) {
  v <- clin %>% filter(!is.na(.data[[var]])) %>% select(sample_id, subject_id, week, value = all_of(var))
  direct <- meta %>% inner_join(v %>% select(sample_id, value) %>% distinct(sample_id, .keep_all = TRUE), by = "sample_id") %>%
    transmute(sample_id, value, how = "same sample", gap = 0)
  rest <- meta %>% filter(!sample_id %in% direct$sample_id) %>%
    inner_join(v %>% filter(!is.na(week)) %>% select(subject_id, week_v = week, value), by = "subject_id",
               relationship = "many-to-many") %>%
    mutate(gap = abs(week - week_v)) %>% filter(gap <= max_gap) %>%
    group_by(sample_id) %>% slice_min(gap, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(sample_id, value, how = ifelse(gap == 0, "same visit", "nearest visit"), gap)
  bind_rows(direct, rest) %>% rename_with(~ paste0(var, "_", .x), -sample_id)
}

out <- meta %>%
  left_join(match_var("hbi"), by = "sample_id") %>%
  left_join(match_var("calpro"), by = "sample_id") %>%
  left_join(match_var("abx"), by = "sample_id") %>%
  rename(hbi = hbi_value, calpro = calpro_value, abx = abx_value) %>%
  mutate(hbi = as.numeric(hbi), calpro = as.numeric(calpro),
         hbi_cat = cut(hbi, c(-Inf, 4, 7, 16, Inf), labels = c("remission (<5)", "mild (5-7)", "moderate (8-16)", "severe (>16)")),
         calpro_cat = cut(calpro, c(-Inf, 50, 250, Inf), labels = c("<50", "50-250", ">250")))
write_csv(out, file.path(outdir, "hmp2_clinical_per_sample.csv"))

# ---- Summary ------------------------------------------------------------
cat("\nCoverage per diagnosis (samples with value):\n")
out %>% group_by(diagnosis) %>%
  summarise(samples = n(), hbi = sum(!is.na(hbi)), calpro = sum(!is.na(calpro)), abx = sum(!is.na(abx)),
            persons = n_distinct(subject_id),
            persons_with_hbi = n_distinct(subject_id[!is.na(hbi)]),
            persons_with_calpro = n_distinct(subject_id[!is.na(calpro)])) %>% print()

cat("\nHow matched:\n")
print(table(hbi = out$hbi_how, useNA = "ifany")); print(table(calpro = out$calpro_how, useNA = "ifany"))

cat("\nHBI (CD samples only):\n"); print(table(out$hbi_cat[out$diagnosis == "cd"], useNA = "ifany"))
cat("\nCalprotectin (ug/g) median per diagnosis:\n")
out %>% filter(!is.na(calpro)) %>% group_by(diagnosis) %>%
  summarise(n = n(), median = median(calpro), q1 = quantile(calpro, .25), q3 = quantile(calpro, .75)) %>% print()
cat("\nCalprotectin categories:\n"); print(table(out$diagnosis, out$calpro_cat, useNA = "ifany"))
cat("\nAntibiotics since last visit:\n"); print(table(out$diagnosis, out$abx, useNA = "ifany"))
cat("\nPersons with at least one antibiotics Yes:\n")
out %>% filter(abx == "Yes") %>% distinct(subject_id, diagnosis) %>% count(diagnosis) %>% print()

# per person: median activity (CD), e.g. for later subgroups
person <- out %>% filter(diagnosis == "cd") %>% group_by(subject_id) %>%
  summarise(n = n(), hbi_median = median(hbi, na.rm = TRUE), calpro_median = median(calpro, na.rm = TRUE),
            any_abx = any(abx == "Yes", na.rm = TRUE))
write_csv(person, file.path(outdir, "hmp2_clinical_cd_per_person.csv"))
cat("\nCD persons: HBI median < 5:", sum(person$hbi_median < 5, na.rm = TRUE), "of", sum(!is.na(person$hbi_median)),
    "| calpro median > 250:", sum(person$calpro_median > 250, na.rm = TRUE), "of", sum(!is.na(person$calpro_median)), "\n")
cat("Spearman HBI vs calprotectin (CD samples): ")
d <- out %>% filter(diagnosis == "cd", !is.na(hbi), !is.na(calpro))
ct <- suppressWarnings(cor.test(d$hbi, d$calpro, method = "spearman")); cat(round(ct$estimate, 2), "p =", signif(ct$p.value, 2), " n =", nrow(d), "\n")
