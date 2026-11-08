# --- Extreme genera and samples ---
# Which genera are the most and least abundant, and which samples have the lowest and
# highest Shannon diversity? One overview per cohort, after the 10 % prevalence filter.
# Started by run_levels.R, which sets level.
library(tidyverse)
library(vegan)
source("R/00_setup.R")

unit   <- if (level == "genus") "genera" else "species"
outdir <- if (level == "genus") file.path(dir_results, "extremes") else file.path(dir_results, "species", "extremes")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

prev_thr <- 0.10
n_top    <- 10
cohorts  <- c(franzosa = "Franzosa", hmp2 = "HMP2", hall = "Hall")

load_cohort <- function(cohort) {
  folder <- file.path(dir_processed, cohort)
  abund <- read_csv(file.path(folder, paste0(cohort, "_", level, ".csv")), show_col_types = FALSE)
  meta  <- read_csv(file.path(folder, paste0(cohort, "_meta.csv")), show_col_types = FALSE) %>%
    transmute(sample_id, subject_id, diagnosis = factor(diagnosis, levels = c("control", "cd")))
  list(abund = abund, meta = meta)
}

# relative abundances after the prevalence filter, same steps as 04_load_abundance.R
prepare <- function(d, cohort) {
  X <- d$abund %>%
    filter(sample_id %in% d$meta$sample_id) %>%
    column_to_rownames("sample_id") %>% as.matrix()
  X <- X / rowSums(X)
  X <- X[, passes_filter(X, cohort, level, prev_thr), drop = FALSE]
  X <- X[rowSums(X) > 0, , drop = FALSE]
  X / rowSums(X)
}

data <- imap(cohorts, function(label, key) {
  d <- load_cohort(key)
  X <- prepare(d, key)
  list(X = X, meta = d$meta %>% filter(sample_id %in% rownames(X)))
})

# sorts the bars separately within each panel; the ___ part is removed again in the labels
sort_within <- function(x, by, within) reorder(paste(x, within, sep = "___"), by)
remove_suffix <- function(x) sub("___.*$", "", x)

# --- Most and least abundant genera ---
genera <- imap_dfr(data, function(d, key) {
  X <- d$X
  dx <- d$meta$diagnosis[match(rownames(X), d$meta$sample_id)]
  tibble(cohort       = cohorts[[key]],
         genus        = colnames(X),
         mean_pct     = colMeans(X) * 100,
         prevalence   = colMeans(X > 0) * 100,
         mean_cd      = colMeans(X[dx == "cd", , drop = FALSE]) * 100,
         mean_control = colMeans(X[dx == "control", , drop = FALSE]) * 100)
})

top_genera <- bind_rows(
  genera %>% group_by(cohort) %>% slice_max(mean_pct, n = n_top, with_ties = FALSE) %>% mutate(group = "most abundant"),
  genera %>% group_by(cohort) %>% slice_min(mean_pct, n = n_top, with_ties = FALSE) %>% mutate(group = "least abundant")
) %>% ungroup() %>% arrange(cohort, group, desc(mean_pct))

print(top_genera, n = Inf)
write_csv(top_genera, file.path(outdir, paste0(unit, "_top_bottom_overall.csv")))

plot_genera <- function(which, file, log_x = FALSE) {
  d <- top_genera %>% filter(group == which) %>%
    pivot_longer(c(mean_cd, mean_control), names_to = "diagnosis", values_to = "pct") %>%
    mutate(diagnosis = factor(sub("mean_", "", diagnosis), levels = c("control", "cd")),
           cohort = factor(cohort, levels = cohorts),
           genus_lab = sort_within(genus, mean_pct, cohort))
  p <- ggplot(d, aes(pct, genus_lab, fill = diagnosis)) +
    geom_col(position = position_dodge(width = .75), width = .7) +
    facet_wrap(~ cohort, scales = "free_y") +
    scale_y_discrete(labels = remove_suffix) +
    scale_fill_manual(values = col_diag) +
    labs(x = "mean relative abundance (%)", y = NULL, fill = NULL,
         title = sprintf("%d %s %s per cohort (%s %% prevalence filter)", n_top, which, unit, prev_thr * 100)) +
    theme_plot
  if (log_x) p <- p + scale_x_log10()
  ggsave(file.path(outdir, file), p, width = 12, height = 5, dpi = 200)
}
plot_genera("most abundant",  paste0(unit, "_most_abundant_overall.png"))
plot_genera("least abundant", paste0(unit, "_least_abundant_overall.png"), log_x = TRUE)

# --- Most and least abundant genera, separately for controls and CD ---
# Samples are averaged per person first, so that people with many samples do not count more.
genera <- imap_dfr(data, function(d, key) {
  m <- d$meta[match(rownames(d$X), d$meta$sample_id), ]
  map_dfr(c("control", "cd"), function(g) {
    in_group <- m$diagnosis == g
    per_person <- apply(d$X[in_group, , drop = FALSE], 2, function(x) tapply(x, m$subject_id[in_group], mean))
    tibble(cohort     = cohorts[[key]],
           diagnosis  = g,
           genus      = colnames(d$X),
           people     = nrow(per_person),
           mean_pct   = colMeans(per_person) * 100,
           median_pct = apply(per_person, 2, median) * 100,
           prevalence = colMeans(per_person > 0) * 100)
  })
})

# least abundant: only taxa that occur in this group at all
top_genera <- bind_rows(
  genera %>% group_by(cohort, diagnosis) %>% slice_max(mean_pct, n = n_top, with_ties = FALSE) %>% mutate(group = "most abundant"),
  genera %>% filter(mean_pct > 0) %>% group_by(cohort, diagnosis) %>% slice_min(mean_pct, n = n_top, with_ties = FALSE) %>% mutate(group = "least abundant")
) %>% ungroup() %>% arrange(cohort, diagnosis, group, desc(mean_pct))

print(top_genera, n = Inf)
write_csv(top_genera, file.path(outdir, paste0(unit, "_top_bottom.csv")))

plot_genera <- function(which, file, log_x = FALSE) {
  d <- top_genera %>% filter(group == which) %>%
    mutate(diagnosis = factor(diagnosis, levels = c("control", "cd")),
           cohort = factor(cohort, levels = cohorts),
           genus_lab = sort_within(genus, mean_pct, paste(cohort, diagnosis)))
  p <- ggplot(d, aes(mean_pct, genus_lab, fill = diagnosis)) +
    facet_wrap(cohort ~ diagnosis, scales = "free_y", ncol = 2) +
    scale_y_discrete(labels = remove_suffix) +
    scale_fill_manual(values = col_diag) +
    labs(x = "mean relative abundance per person (%)", y = NULL, fill = NULL,
         title = sprintf("%d %s %s per cohort and group (%s %% prevalence filter)", n_top, which, unit, prev_thr * 100)) +
    theme_plot
  # bars have no meaningful start on a log axis, so the log plot uses points
  p <- if (log_x) p + geom_point(shape = 21, size = 3) + scale_x_log10() else p + geom_col(width = .7)
  ggsave(file.path(outdir, file), p, width = 11, height = 11, dpi = 200)
}
plot_genera("most abundant",  paste0(unit, "_most_abundant.png"))
plot_genera("least abundant", paste0(unit, "_least_abundant.png"), log_x = TRUE)

# --- Samples with the lowest and highest Shannon diversity ---
samples <- imap_dfr(data, function(d, key) {
  tibble(cohort    = cohorts[[key]],
         sample_id = rownames(d$X),
         Shannon   = diversity(d$X, index = "shannon"),
         Richness  = rowSums(d$X > 0)) %>%
    left_join(d$meta, by = "sample_id")
})

extremes <- bind_rows(
  samples %>% group_by(cohort) %>% slice_min(Shannon, n = n_top, with_ties = FALSE) %>% mutate(group = "lowest diversity"),
  samples %>% group_by(cohort) %>% slice_max(Shannon, n = n_top, with_ties = FALSE) %>% mutate(group = "highest diversity")
) %>% ungroup() %>% arrange(cohort, group, Shannon)

print(extremes %>% select(cohort, group, sample_id, subject_id, diagnosis, Shannon, Richness), n = Inf)
write_csv(extremes, file.path(outdir, "samples_lowest_highest_shannon.csv"))

# with several samples per person, the extremes often come from the same person
print(extremes %>% group_by(cohort, group) %>%
  summarise(samples = n(), people = n_distinct(subject_id), .groups = "drop"))

d <- extremes %>%
  mutate(cohort = factor(cohort, levels = cohorts),
         group  = factor(group, levels = c("lowest diversity", "highest diversity")),
         label  = ifelse(sample_id == subject_id, sample_id, paste0(sample_id, " (", subject_id, ")")),
         label  = sort_within(label, Shannon, paste(cohort, group)))

plot_samples <- function(d, title, label_size) {
  ggplot(d, aes(Shannon, label, colour = diagnosis)) +
    geom_segment(aes(x = 0, xend = Shannon, yend = label), colour = "grey85", linewidth = 1) +
    geom_point(size = 3) +
    scale_y_discrete(labels = remove_suffix) +
    scale_colour_manual(values = col_diag) +
    labs(x = "Shannon diversity", y = NULL, colour = NULL, title = title) +
    theme_plot + theme(axis.text.y = element_text(size = label_size))
}

p <- plot_samples(d, sprintf("%d samples with the lowest and highest diversity per cohort", n_top), 8) +
  facet_wrap(group ~ cohort, scales = "free_y", ncol = 3)
ggsave(file.path(outdir, "samples_lowest_highest_shannon.png"), p, width = 14, height = 8, dpi = 200)

# the same plots separately, with larger labels
for (g in levels(d$group)) {
  p <- plot_samples(filter(d, group == g), sprintf("%d samples with the %s per cohort", n_top, g), 10) +
    facet_wrap(~ cohort, scales = "free_y", ncol = 3) +
    scale_x_continuous(limits = c(0, ceiling(max(d$Shannon))))
  ggsave(file.path(outdir, paste0("samples_", sub(" diversity", "", g), "_shannon.png")), p, width = 14, height = 5, dpi = 200)
}
