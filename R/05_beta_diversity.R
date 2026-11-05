# --- Beta diversity ---
# Does the compositio of the samples differ between CD and controls?
# Started by run_levels.R, which sets cohort and level.

library(tidyverse)
library(vegan)
library(uwot)    
library(limma)     
library(sva)       
library(MMUPHin)   
source("R/00_setup.R")

source("R/04_load_abundance.R")
alpha <- readRDS(file.path(outdir, "alpha_results.rds"))$alpha

# CLR with pseudocount of half the smallest non-zero value
clr_tf <- function(X) {
  pc <- min(X[X > 0]) / 2
  L <- log(X + pc)
  L - rowMeans(L)
}

#PCA and UMAP coordinates
pca_umap <- function(Z) {
  pca <- prcomp(Z, center = TRUE, scale. = FALSE)
  set.seed(seed)
  um <- umap(Z, n_neighbors = min(15, nrow(Z) - 1), min_dist = 0.3)
  tibble(sample_id = rownames(Z),
         PC1 = pca$x[, 1], PC2 = pca$x[, 2],
         UMAP1 = um[, 1], UMAP2 = um[, 2],
         var1 = 100 * summary(pca)$importance[2, 1],
         var2 = 100 * summary(pca)$importance[2, 2])
}

# Aitchison (Euclidean on CLR) is the main distance
# Bray-Curtis the sensitivity analysis
distances <- function(X, ids) {
  list(aitchison = dist(clr_tf(X)[ids, ]),
       bray      = vegdist(X[ids, ], method = "bray"))
}

# --- PCA and UMAP at 0, 5 and 10 % ---
ord <- imap_dfr(X_thr, function(X, nm) {
  sc <- pca_umap(clr_tf(X)) %>% mutate(threshold = nm)
  # outlier: distance to the median above median + 3 MAD
  d <- sqrt((sc$PC1 - median(sc$PC1))^2 + (sc$PC2 - median(sc$PC2))^2)
  sc$outlier <- d > median(d) + 3 * mad(d)
  sc
}) %>%
  left_join(meta, by = "sample_id") %>%
  left_join(select(few, sample_id, threshold, n_genera, few_genera),
            by = c("sample_id", "threshold"))

print(ord %>% distinct(threshold, var1, var2))

p_pca <- ggplot(ord, aes(PC1, PC2, colour = diagnosis, shape = outlier)) +
  geom_point(alpha = .7) + facet_wrap(~ threshold, scales = "free") +
  theme_plot
save_plot(p_pca, "pca_thresholds.png", 11, 4)

# do samples with few genera lie at the edge?
p_pca_n <- ggplot(ord, aes(PC1, PC2, colour = n_genera)) +
  geom_point(alpha = .8) + scale_colour_viridis_c() + facet_wrap(~ threshold, scales = "free") +
  theme_plot
save_plot(p_pca_n, "pca_n_genera.png", 11, 4)

p_umap <- ggplot(ord, aes(UMAP1, UMAP2, colour = diagnosis, shape = few_genera)) +
  geom_point(alpha = .8) + facet_wrap(~ threshold, scales = "free") +
  theme_plot
save_plot(p_umap, "umap_thresholds.png", 11, 4)

# samples of the same person connected in time order
if (longitudinal) {
  p_traj <- ord %>% filter(threshold == main_thr) %>%
    arrange(subject_id, time) %>%
    ggplot(aes(PC1, PC2, group = subject_id, colour = diagnosis)) +
    geom_path(alpha = .3) + geom_point(alpha = .7) +
    theme_plot
  save_plot(p_traj, "pca_subjects.png", 7, 5)
}

# --- PCoA with Bray-Curtis ---
# Bray-Curtis is not Euclidean, so PCoA instead of PCA
pcoa <- imap_dfr(X_thr, function(X, nm) {
  pc <- cmdscale(vegdist(X, method = "bray"), k = 2, eig = TRUE)
  expl <- 100 * pc$eig / sum(pc$eig[pc$eig > 0])
  tibble(sample_id = rownames(X), threshold = nm,
         PCoA1 = pc$points[, 1], PCoA2 = pc$points[, 2],
         var1 = expl[1], var2 = expl[2])
}) %>% left_join(meta, by = "sample_id")

print(pcoa %>% distinct(threshold, var1, var2))

p_pcoa <- ggplot(pcoa, aes(PCoA1, PCoA2, colour = diagnosis)) +
  geom_point(alpha = .7) + facet_wrap(~ threshold, scales = "free") +
  labs(title = "PCoA, Bray-Curtis") +
  theme_plot
save_plot(p_pcoa, "pcoa_bray_thresholds.png", 11, 4)

# --- PERMANOVA: do CD and controls differ? ---
# first sample of each person only, so that nobody counts more than once
permanova <- imap_dfr(X_thr, function(X, nm) {
  ids <- intersect(rownames(X), meta$sample_id[meta$baseline])
  d <- meta[match(ids, meta$sample_id), ]
  imap_dfr(distances(X, ids), function(D, dn) {
    set.seed(seed)
    a <- adonis2(D ~ diagnosis, data = d, permutations = 999)
    # a different spread of the groups can also make PERMANOVA significant
    set.seed(seed)
    bd <- permutest(betadisper(D, d$diagnosis), permutations = 999)
    tibble(threshold = nm, distance = dn, n = length(ids),
           R2 = a$R2[1], p = a$`Pr(>F)`[1], p_dispersion = bd$tab$`Pr(>F)`[1])
  })
})
print(permanova %>% arrange(distance, threshold))
write_csv(permanova, file.path(outdir, "permanova.csv"))

# without persons on antibiotics (Franzosa only)
if (any(meta$antibiotic == "Yes", na.rm = TRUE)) {
  permanova_abx <- imap_dfr(X_thr, function(X, nm) {
    ids <- intersect(rownames(X), meta$sample_id[meta$baseline & (meta$antibiotic != "Yes" | is.na(meta$antibiotic))])
    d <- meta[match(ids, meta$sample_id), ]
    imap_dfr(distances(X, ids), function(D, dn) {
      set.seed(seed)
      a <- adonis2(D ~ diagnosis, data = d, permutations = 999)
      tibble(threshold = nm, distance = dn, n = length(ids), R2 = a$R2[1], p = a$`Pr(>F)`[1])
    })
  })
  print(permanova_abx %>% arrange(distance, threshold))
  write_csv(permanova_abx, file.path(outdir, "permanova_without_antibiotics.csv"))
}

# diagnosis effect adjusted for centre
base_meta <- meta %>% filter(baseline)
if (n_distinct(base_meta$site) > 1) {
  permanova_site <- imap_dfr(X_thr, function(X, nm) {
    ids <- intersect(rownames(X), base_meta$sample_id)
    d <- meta[match(ids, meta$sample_id), ]
    imap_dfr(distances(X, ids), function(D, dn) {
      set.seed(seed)
      a <- tryCatch(adonis2(D ~ site + diagnosis, data = d, by = "margin", permutations = 999), error = function(e) NULL)
      if (is.null(a)) return(NULL)
      tibble(threshold = nm, distance = dn, R2_diagnosis = a["diagnosis", "R2"], p_diagnosis = a["diagnosis", "Pr(>F)"],
             R2_site = a["site", "R2"], p_site = a["site", "Pr(>F)"])
    })
  })
  print(permanova_site %>% arrange(distance, threshold))
  write_csv(permanova_site, file.path(outdir, "permanova_by_centre.csv"))
}

# --- Batch correction, for figures only ---
# batch = site. The tests above use site as a covariate. They must not run on corrected data,
# because the correction keeps diagnosis in the model and makes p-values too optimistic.
# limma and ComBat work on CLR, MMUPHin on relative abundances (zeros stay zeros).
X_main <- X_thr[[main_thr]]
Z_main <- clr_tf(X_main)
m_main <- meta[match(rownames(Z_main), meta$sample_id), ]

batch_design <- m_main %>% distinct(subject_id, diagnosis, site) %>%
  count(site, diagnosis, name = "persons") %>%
  pivot_wider(names_from = diagnosis, values_from = persons, values_fill = 0)
print(batch_design)
write_csv(batch_design, file.path(outdir, "batch_design.csv"))

# if a batch has only one diagnosis, no correction can separate the two
batch_confounded <- any(batch_design$control == 0 | batch_design$cd == 0)
several_sites <- n_distinct(m_main$site) > 1
if (several_sites && batch_confounded)
  warning("batch and diagnosis are partly confounded (one batch has only one diagnosis) -> correction unreliable")

if (several_sites) {
  design <- model.matrix(~ diagnosis, data = m_main)   # protects the diagnosis effect
  
  # limma and ComBat want samples in columns
  Z_limma <- t(removeBatchEffect(t(Z_main), batch = m_main$site, design = design))
  Z_combat <- t(ComBat(dat = t(Z_main), batch = m_main$site, mod = design))
  
  mm_data <- m_main %>% mutate(site = factor(site)) %>% column_to_rownames("sample_id")
  mmuphin <- adjust_batch(feature_abd = t(X_main), batch = "site", covariates = "diagnosis",
                          data = mm_data, control = list(verbose = FALSE))
  X_mmuphin <- t(mmuphin$feature_abd_adj)
  X_mmuphin <- X_mmuphin / rowSums(X_mmuphin)
  Z_mmuphin <- clr_tf(X_mmuphin)
  
  Z_list <- list(before = Z_main, limma = Z_limma, ComBat = Z_combat, MMUPHin = Z_mmuphin)
  
  ord_batch <- map_dfr(Z_list, pca_umap, .id = "correction") %>%
    mutate(correction = factor(correction, levels = names(Z_list))) %>%
    left_join(meta, by = "sample_id")
  
  p_ub <- ggplot(ord_batch, aes(UMAP1, UMAP2, colour = diagnosis, shape = site)) +
    geom_point(alpha = .7) + facet_wrap(~ correction, scales = "free", ncol = 2) +
    labs(title = paste0("UMAP before/after batch correction (", main_thr, ")"),
         subtitle = if (batch_confounded) "Note: batch and diagnosis partly confounded" else NULL) +
    theme_plot
  save_plot(p_ub, "umap_batch_before_after.png", 11, 9)
  
  p_pb <- ggplot(ord_batch, aes(PC1, PC2, colour = diagnosis, shape = site)) +
    geom_point(alpha = .7) + facet_wrap(~ correction, scales = "free", ncol = 2) +
    labs(title = paste0("PCA before/after batch correction (", main_thr, ")")) +
    theme_plot
  save_plot(p_pb, "pca_batch_before_after.png", 11, 9)
  
  # variation explained by batch before and after, only to describe the correction
  batch_r2 <- map_dfr(Z_list, function(Z) {
    ids <- intersect(rownames(Z), meta$sample_id[meta$baseline])
    d <- meta[match(ids, meta$sample_id), ]
    set.seed(seed)
    a <- adonis2(dist(Z[ids, ]) ~ site, data = d, permutations = 999)
    tibble(R2_batch = a$R2[1], p_batch = a$`Pr(>F)`[1])
  }, .id = "correction")
  print(batch_r2)
  write_csv(batch_r2, file.path(outdir, "batch_r2.csv"))
}

# --- Main figures ---
# ComBat if every batch has both diagnoses, otherwise uncorrected (Hall: Atlanta has only controls)
use_combat <- several_sites && !batch_confounded
Z_fig <- if (use_combat) Z_combat else Z_main
label_fig <- if (use_combat) "ComBat-corrected, batch = site" else "without batch correction"
note_fig <- c(
  if (use_combat && min(as.matrix(batch_design[, c("control", "cd")])) < 3)
    "Note: one batch has fewer than 3 persons of one diagnosis",
  if (!use_combat && several_sites)
    "no correction: batch and diagnosis confounded",
  if (longitudinal) "all time points; dense clusters are usually one person")

ord_fig <- pca_umap(Z_fig) %>% left_join(meta, by = "sample_id")

p_uh <- ggplot(ord_fig, aes(UMAP1, UMAP2, colour = diagnosis, shape = site)) +
  geom_point(alpha = .75, size = 1.8) +
  labs(title = paste0("UMAP ", cohort, " (", main_thr, ", ", label_fig, ")"),
       subtitle = paste(note_fig, collapse = "\n")) +
  theme_plot
save_plot(p_uh, "umap_main.png", 8, 6)

p_ph <- ggplot(ord_fig, aes(PC1, PC2, colour = diagnosis, shape = site)) +
  geom_point(alpha = .75, size = 1.8) +
  labs(x = paste0("PC1 (", round(ord_fig$var1[1], 1), " %)"), y = paste0("PC2 (", round(ord_fig$var2[1], 1), " %)"),
       title = paste0("PCA ", cohort, " (", main_thr, ", ", label_fig, ")"),
       subtitle = paste(note_fig, collapse = "\n")) +
  theme_plot
save_plot(p_ph, "pca_main.png", 8, 6)

# are the clusters driven by the person rather than the diagnosis?
if (longitudinal) {
  p_up <- ord %>% filter(threshold == main_thr) %>%
    ggplot(aes(UMAP1, UMAP2, colour = subject_id, shape = diagnosis)) +
    geom_point(alpha = .8) + guides(colour = "none") +
    labs(title = paste0("UMAP coloured by person (", main_thr, ")")) +
    theme_plot
  save_plot(p_up, "umap_person.png", 7, 5)
}

# --- Samples that stand out ---
outliers <- ord %>%
  group_by(sample_id, subject_id, diagnosis) %>%
  summarise(n_thr_outlier = sum(outlier),
            ever_few_genera = any(few_genera),
            min_n_genera = min(n_genera),
            total_raw = rs[first(sample_id)], .groups = "drop") %>%
  filter(n_thr_outlier > 0) %>%
  arrange(desc(n_thr_outlier))
print(outliers)
write_csv(outliers, file.path(outdir, "outliers.csv"))
saveRDS(list(outliers = outliers), file.path(outdir, "beta_results.rds")) 

# same alpha test without the outliers
# a person whose first sample is an outlier drops out
outlier_test <- alpha %>% filter(baseline) %>%
  mutate(keep = !sample_id %in% outliers$sample_id) %>%
  group_by(threshold) %>%
  summarise(n_with = n(), p_with = wilcox.test(Shannon ~ diagnosis)$p.value,
            n_without = sum(keep),
            p_without = wilcox.test(Shannon[keep] ~ diagnosis[keep])$p.value)
print(outlier_test)