# --- 10 Diet- and CD-associated taxa ---
# Tests a literature list of taxa linked to diet (in healthy people) or to CD for a CD effect
# Franzosa and HMP2 only; exploratory (see methods 12.1).

# Results (per cohort):
# - CD effect of 17 taxa (MaAsLin2, 95 % CI, Holm over the list), direction vs. diet and CD literature
# - forest plots of the diet- and the CD-associated taxa: only taxa with Holm p < 0.05 in at least one
#   cohort, and all 17 taxa (_all); filled points: Holm p < 0.05
# - prevalence and median abundance when present of the 17 taxa, CD vs. control (descriptive)

# Methods: docs/methods.md, section 12

library(tidyverse)
library(Maaslin2)
library(patchwork)
source("R/00_setup.R")

outdir <- file.path(dir_results, "hypotheses")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cohorts  <- c("franzosa", "hmp2")
hmp2_tax <- read_csv(file.path(dir_processed, "hmp2/hmp2_taxonomy.csv"), show_col_types = FALSE)

# --- Literature list (see methods 12.1) ---
# diet: link in diet studies of healthy people; cd: direction reported in CD/IBD; NA = no link
# level genus: all species of the genus summed, as in section 11
lit <- tribble(
  ~taxon,                         ~level,    ~label,            ~diet,          ~cd,
  "Faecalibacterium_prausnitzii", "species", "F. prausnitzii",  "fibre/plant",  "lower",
  "Roseburia",                    "genus",   "Roseburia",       "fibre/plant",  "lower",
  "Eubacterium_rectale",          "species", "E. rectale",      "fibre/plant",  NA,
  "Anaerostipes_hadrus",          "species", "A. hadrus",       "fibre/plant",  NA,
  "Eubacterium_hallii",           "species", "E. hallii",       "fibre/plant",  NA,
  "Coprococcus",                  "genus",   "Coprococcus",     NA,             "lower",
  "Ruminococcus_bromii",          "species", "R. bromii",       "fibre/plant",  NA,
  "Bifidobacterium_adolescentis", "species", "B. adolescentis", "fibre/plant",  "lower",
  "Prevotella",                   "genus",   "Prevotella",      "fibre/plant",  NA,
  "Bacteroides",                  "genus",   "Bacteroides",     "animal/fat",   NA,
  "Alistipes",                    "genus",   "Alistipes",       "animal/fat",   NA,
  "Bilophila_wadsworthia",        "species", "B. wadsworthia",  "animal/fat",   NA,
  "Escherichia_coli",             "species", "E. coli",         NA,             "higher",
  "Ruminococcus_gnavus",          "species", "R. gnavus",       NA,             "higher",
  "Fusobacterium_nucleatum",      "species", "F. nucleatum",    NA,             "higher",
  "Dialister_invisus",            "species", "D. invisus",      NA,             "lower",
  "Akkermansia_muciniphila",      "species", "A. muciniphila",  NA,             "lower") %>%
  # expected in CD if CD patients eat less fibre / more animal food (diet proxy, see methods 0)
  mutate(diet_expected = recode(diet, "fibre/plant" = "lower", "animal/fat" = "higher"))

# --- Load ---
load_cohort <- function(cohort) {
  folder <- file.path(dir_processed, cohort)
  meta <- read_csv(file.path(folder, paste0(cohort, "_meta.csv")), show_col_types = FALSE)
  X <- read_csv(file.path(folder, paste0(cohort, "_species.csv")), show_col_types = FALSE) %>%
    filter(sample_id %in% meta$sample_id) %>%
    arrange(match(sample_id, meta$sample_id)) %>%
    column_to_rownames("sample_id") %>% as.matrix()
  list(meta = meta, X = X / rowSums(X))
}
data <- set_names(map(cohorts, load_cohort), cohorts)

# relative abundance of each listed taxon; no prevalence filter (see methods 11.1)
taxon_values <- function(d, cohort) {
  g <- if (cohort == "hmp2") hmp2_tax$genus[match(colnames(d$X), hmp2_tax$species)] else sub("_.*", "", colnames(d$X))
  stopifnot(lit$taxon[lit$level == "species"] %in% colnames(d$X))
  vals <- map(set_names(lit$taxon), function(t)
    if (lit$level[lit$taxon == t] == "genus") rowSums(d$X[, g %in% t, drop = FALSE]) else d$X[, t])
  as.data.frame(vals, row.names = rownames(d$X))
}

# --- MaAsLin2, same settings as H2b (see methods 11.3, 12.1) ---
run_maaslin <- function(cohort, random) {
  d <- data[[cohort]]
  run_dir <- file.path(outdir, "maaslin2", paste0(cohort, "_cd_diet"))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  m <- d$meta %>% select(sample_id, diagnosis, site, subject_id) %>%
    column_to_rownames("sample_id") %>% as.data.frame()
  Maaslin2(input_data = taxon_values(d, cohort), input_metadata = m, output = run_dir,
           min_abundance = 0, min_prevalence = 0, min_variance = 0,
           normalization = "NONE", transform = "LOG", analysis_method = "LM",
           correction = "holm", standardize = FALSE,
           fixed_effects = c("diagnosis", "site"), random_effects = random,
           reference = "diagnosis,control",
           plot_heatmap = FALSE, plot_scatter = FALSE, cores = 1)
  # Holm over the 17 taxa, diagnosis term only
  read_tsv(file.path(run_dir, "all_results.tsv"), show_col_types = FALSE) %>%
    filter(metadata == "diagnosis") %>%
    transmute(cohort, taxon = feature, samples = N, not_zero = N.not.0, coef, stderr,
              ci_low = coef - qnorm(0.975) * stderr, ci_high = coef + qnorm(0.975) * stderr,
              pval, p_holm = p.adjust(pval, "holm"))
}

res <- bind_rows(run_maaslin("franzosa", NULL), run_maaslin("hmp2", "subject_id")) %>%
  left_join(lit, by = "taxon") %>%
  mutate(observed = if_else(coef < 0, "lower", "higher"),
         matches_diet = if_else(is.na(diet), NA, observed == diet_expected),
         matches_cd   = if_else(is.na(cd),   NA, observed == cd))
print(res %>% select(cohort, label, coef, ci_low, ci_high, p_holm, diet, cd, matches_diet, matches_cd),
      n = Inf, width = Inf)
write_csv(res, file.path(outdir, "species_cd_diet.csv"))

# --- Forest plots: taxa with Holm p < 0.05 in at least one cohort, and all taxa (see methods 12.1) ---
# link: column that groups the plot (diet or cd); expected: column with the expected direction
forest <- function(res, link, expected, title, file, only_sig = TRUE) {
  d <- res %>% filter(!is.na(.data[[link]])) %>%
    group_by(taxon) %>% filter(!only_sig | any(p_holm < 0.05)) %>% ungroup() %>%
    mutate(cohort = recode(cohort, franzosa = "Franzosa", hmp2 = "HMP2"),
           label = paste0(label, if_else(.data[[expected]] == "lower", " (↓)", " (↑)")),
           label = fct_reorder(label, coef, .fun = min),
           facet = factor(.data[[link]], c("fibre/plant", "animal/fat", "lower", "higher"),
                          c("Fibre / plant", "Animal / fat", "Lower in CD", "Higher in CD")),
           significant = if_else(p_holm < 0.05, "Holm p < 0.05", "not significant"))
  p <- ggplot(d, aes(coef, label, colour = cohort)) +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y",
                  position = position_dodge(width = .6)) +
    geom_point(aes(shape = significant), size = 2.6, fill = "white", position = position_dodge(width = .6)) +
    scale_shape_manual(values = c("Holm p < 0.05" = 16, "not significant" = 21), name = NULL) +
    scale_colour_manual(values = c(Franzosa = okabe_ito[["green"]], HMP2 = okabe_ito[["orange"]]), name = NULL) +
    facet_grid(facet ~ ., scales = "free_y", space = "free_y") +
    labs(x = "CD vs. control (log2 fold change, 95 % CI)", y = NULL, title = title,
         caption = "Arrow: expected direction in CD. MaAsLin2; HMP2 with person as random effect.") +
    theme_plot +
    theme(strip.text.y = element_text(angle = 0, hjust = 0))
  ggsave(file.path(outdir, file), p, width = 8, height = if (only_sig) 6 else 1.6 + 0.42 * n_distinct(d$taxon), dpi = 200)
}

forest(res, "diet", "diet_expected", "Diet-associated taxa in CD vs. controls", "species_diet_forest.png")
forest(res, "cd", "cd", "CD-associated taxa in CD vs. controls", "species_cd_forest.png")
forest(res, "diet", "diet_expected", "Diet-associated taxa in CD vs. controls (all)", "species_diet_forest_all.png", only_sig = FALSE)
forest(res, "cd", "cd", "CD-associated taxa in CD vs. controls (all)", "species_cd_forest_all.png", only_sig = FALSE)

# --- How common are the taxa? Prevalence and abundance when present (descriptive, see methods 12.1) ---
# HMP2: all samples, so persons with many samples count more; no test
prev_abund <- imap_dfr(data, function(d, cohort) {
  as_tibble(taxon_values(d, cohort), rownames = "sample_id") %>%
    left_join(d$meta %>% select(sample_id, diagnosis), by = "sample_id") %>%
    pivot_longer(all_of(lit$taxon), names_to = "taxon", values_to = "abund") %>%
    group_by(diagnosis, taxon) %>%
    summarise(prevalence = 100 * mean(abund > 0),
              abundance  = 100 * median(abund[abund > 0]), .groups = "drop") %>%
    mutate(cohort = cohort)
}) %>%
  left_join(lit %>% select(taxon, label), by = "taxon") %>%
  mutate(cohort = recode(cohort, franzosa = "Franzosa", hmp2 = "HMP2"),
         diagnosis = recode(diagnosis, cd = "CD", control = "Control"))
write_csv(prev_abund, file.path(outdir, "species_prevalence_abundance.csv"))

# one panel per cohort and measure; taxa ordered by mean prevalence
taxon_order <- prev_abund %>% group_by(label) %>% summarise(m = mean(prevalence)) %>% arrange(m) %>% pull(label)
prev_panel <- function(coh, measure, xlab, log = FALSE, y_text = TRUE) {
  d <- prev_abund %>% filter(cohort == coh) %>%
    mutate(label = factor(label, taxon_order), value = .data[[measure]])
  p <- ggplot(d, aes(value, label)) +
    geom_line(aes(group = label), colour = "grey75", linewidth = 1) +
    geom_point(aes(colour = diagnosis), size = 2.6) +
    scale_colour_manual(values = c(Control = okabe_ito[["sky_blue"]], CD = okabe_ito[["vermillion"]]), name = NULL) +
    labs(x = xlab, y = NULL, title = paste(coh, if (log) "· abundance" else "· prevalence")) +
    theme_plot + theme(plot.title = element_text(size = 12))
  p <- if (log) p + scale_x_log10(labels = function(x) format(x, drop0trailing = TRUE, scientific = FALSE))
       else p + scale_x_continuous(limits = c(0, 100))
  if (!y_text) p <- p + theme(axis.text.y = element_blank())
  p
}
p_prev <- prev_panel("Franzosa", "prevalence", "% of samples with the taxon") +
  prev_panel("Franzosa", "abundance", "median % when present", log = TRUE, y_text = FALSE) +
  prev_panel("HMP2", "prevalence", "% of samples with the taxon", y_text = FALSE) +
  prev_panel("HMP2", "abundance", "median % when present", log = TRUE, y_text = FALSE) +
  plot_layout(nrow = 1, guides = "collect") +
  plot_annotation(title = "How common are the 17 taxa in CD and controls?",
                  caption = "Prevalence: share of samples with the taxon (HMP2: all 571 samples). Abundance: median relative abundance in the samples that have it (log scale).",
                  theme = theme_plot) &
  theme(legend.position = "top")
ggsave(file.path(outdir, "species_prevalence_abundance.png"), p_prev, width = 13, height = 5.6, dpi = 200)
