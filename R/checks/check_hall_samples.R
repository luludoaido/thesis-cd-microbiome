library(tidyverse)
library(vegan)
source("R/00_setup.R")

# =====================================================================
# Check conspicuous Hall samples with the same three criteria as SKST006:
#   1. composition: which genus dominates, is that plausible?
#   2. sequencing depth: number_reads compared with all Hall samples
#   3. other samples of the same person: is the sample a one-off?
# SKST006 is included as reference.
# =====================================================================

path <- "data"
check_ids <- c("p9216_mo2", "SKST037_3_G103026", "SKST006_10_G102994")

meta  <- read.csv(file.path(path, "raw/cmd/HallAB_2017_metadata.csv"))
genus <- read_csv(file.path(path, "processed/hall/hall_genus_merged.csv"), show_col_types = FALSE)
first_bac <- which(colnames(genus) == "Abiotrophia")
X <- genus %>% select(sample_id, all_of(first_bac:ncol(genus))) %>%
  column_to_rownames("sample_id") %>% as.matrix()
raw_sum <- rowSums(X)                 # share of reads assigned to a genus (in %)
X <- X / rowSums(X)                   # TSS, without prevalence filter

per_sample <- tibble(sample_id = rownames(X),
                     raw_sum   = raw_sum,
                     richness  = rowSums(X > 0),
                     shannon   = diversity(X, "shannon"),
                     top_genus = colnames(X)[apply(X, 1, which.max)],
                     top_pct   = apply(X, 1, max) * 100) %>%
  left_join(meta %>% select(sample_id = X, subject_id, disease_subtype, number_reads,
                            days_from_first_collection, calprotectin, HBI), by = "sample_id")

# ---- 1. Composition: the 5 largest genera per checked sample ------------
for (id in check_ids) {
  cat("\n==", id, "==\n")
  top <- sort(X[id, ], decreasing = TRUE)[1:5] * 100
  print(round(top[top > 0], 2))
  # how common is the dominant genus in all other Hall samples? (median %)
  g <- names(top)[1]
  cat("Median", g, "in all other Hall samples:", round(median(X[rownames(X) != id, g]) * 100, 2), "%\n")
}

# ---- 2. Sequencing depth -------------------------------------------------
cat("\nReads of all Hall samples (millions):\n")
print(round(quantile(per_sample$number_reads / 1e6, c(0, .05, .25, .5, .75, 1), na.rm = TRUE), 1))
per_sample %>% filter(sample_id %in% check_ids) %>%
  mutate(reads_million = round(number_reads / 1e6, 1),
         reads_rank_pct = round(ecdf(per_sample$number_reads)(number_reads) * 100)) %>%
  select(sample_id, reads_million, reads_rank_pct, raw_sum, richness, shannon, top_genus, top_pct) %>%
  print()

# ---- 3. Other samples of the same person ---------------------------------
for (id in check_ids) {
  subj <- per_sample$subject_id[per_sample$sample_id == id]
  cat("\n== Person", subj, "==\n")
  per_sample %>% filter(subject_id == subj) %>% arrange(days_from_first_collection) %>%
    transmute(sample_id, day = days_from_first_collection, reads_million = round(number_reads / 1e6, 1),
              richness, shannon = round(shannon, 2), top_genus, top_pct = round(top_pct, 1),
              calprotectin = round(calprotectin), HBI,
              checked = sample_id == id) %>%
    print(n = Inf)
}

# for comparison: distribution of richness and Shannon over all Hall samples
cat("\nRichness all Hall samples:\n"); print(quantile(per_sample$richness, c(0, .05, .5, 1)))
cat("Shannon all Hall samples:\n");  print(round(quantile(per_sample$shannon, c(0, .05, .5, 1)), 2))
