library(tidyverse)
library(readxl)
source("R/00_setup.R")

# =====================================================================
# EDA for the hypotheses (before the analysis plan)
# Part A: detectability of the hypothesis taxa (SQ1, SQ2) in all three cohorts
#         + taxonomy for the Gevers dysbiosis index (family / order)
# Part B: EDA of measured butyrate (SQ3a, Franzosa only)
# Part C: EDA of the functional potential for butyrate (SQ3b, HMP2 + Franzosa)
#
# Important: NO prevalence filter here. The 10 % filter only applies to the multivariate
# analyses (PCA, UMAP, PERMANOVA); otherwise rare but relevant taxa (e.g. Fusobacterium)
# would be dropped for single taxa and the Gevers index.
# Hall is only used for detectability (no butyrate or enzyme data,
# centre/kit confounded with diagnosis).
# =====================================================================

path   <- "data"
outdir <- file.path(dir_results, "hypotheses_eda")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# exclusions from R/00_setup.R, same as 06/07 (incl. samples below min_reads)
exclude_samples <- c(exclude_samples, hmp2_empty_samples, low_read_samples())

# ---- 0. Metadata per cohort --------------------------------------------
meta_hall <- read.csv(file.path(path, "processed/hall/hall_meta.csv")) %>%
  transmute(sample_id = ID, subject_id = subject_id, diagnosis = tolower(diagnosis),
            time = days_from_first_collection,
            site = ifelse(grepl("^SKST", subject_id), "Atlanta", "Boston"))

meta_hmp2 <- read.csv(file.path(path, "processed/hmp2/hmp2_meta.csv"), check.names = FALSE) %>%
  transmute(sample_id = `External ID`, subject_id = `Participant ID`, diagnosis = tolower(diagnosis),
            time = week_num, site = substr(`Participant ID`, 1, 1), age = consent_age)

meta_fr <- read.csv(file.path(path, "processed/franzosa/franzosa_meta.csv")) %>%
  transmute(sample_id = sample_id, subject_id = sample_id, diagnosis = tolower(Diagnosis),
            time = NA_real_, site = sub("[|].*", "", sample_id), age = Age,
            antibiotic = antibiotic, calprotectin = Fecal.Calprotectin)

meta_all <- bind_rows(hall = meta_hall, hmp2 = meta_hmp2, franzosa = meta_fr, .id = "cohort") %>%
  filter(!sample_id %in% exclude_samples) %>%
  mutate(diagnosis = factor(diagnosis, levels = c("control", "cd")))

# baseline = earliest sample per person (as in R/06_diversity.R)
meta_all <- meta_all %>%
  group_by(cohort, subject_id) %>%
  mutate(baseline = row_number(if_else(is.na(time), 0, time)) == 1) %>%
  ungroup()

# ---- A1. Species tables with taxonomy (long) ----------------------------
# goal: one table cohort, sample_id, species, genus, family, order, abundance (share 0-1)

# Hall: species names + taxonomy from curatedMetagenomicData (new NCBI taxonomy)
hall_sp <- read_csv(file.path(path, "processed/hall/hall_species_merged.csv"), show_col_types = FALSE)
first_sp <- which(colnames(hall_sp) == "visit_number") + 1
hall_tax <- read_csv(file.path(path, "raw/cmd/HallAB_2017_taxonomy.csv"),
                       show_col_types = FALSE, name_repair = "minimal")[, -1] %>%
  distinct(species, .keep_all = TRUE)
hall_long <- hall_sp %>%
  select(sample_id, all_of(first_sp:ncol(hall_sp))) %>%
  pivot_longer(-sample_id, names_to = "species", values_to = "abundance") %>%
  left_join(hall_tax %>% select(species, genus, family, order), by = "species")

# HMP2: full lineage in the column name (MetaPhlAn2)
hmp2_sp <- read_csv(file.path(path, "processed/hmp2/hmp2_species_merged.csv"), show_col_types = FALSE)
sp_cols <- grep("s__", colnames(hmp2_sp), value = TRUE)
hmp2_long <- hmp2_sp %>%
  select(sample_id, all_of(sp_cols)) %>%
  pivot_longer(-sample_id, names_to = "lineage", values_to = "abundance") %>%
  mutate(order   = str_match(lineage, "o__([^|]+)")[, 2],
         family  = str_match(lineage, "f__([^|]+)")[, 2],
         genus   = str_match(lineage, "g__([^|]+)")[, 2],
         species = str_replace_all(str_match(lineage, "s__([^|]+)")[, 2], "_", " ")) %>%
  select(-lineage)

# Franzosa: species names only (MetaPhlAn2) -> look up family/order via the HMP2 lineage
fr_sp <- read_csv(file.path(path, "processed/franzosa/franzosa_features_merged.csv"), show_col_types = FALSE)
first_fr <- which(colnames(fr_sp) == "Methanobrevibacter_smithii")
genus_lookup <- hmp2_long %>% distinct(genus, family, order)
# "genera" that are actually higher ranks (e.g. "Lachnospiraceae_bacterium_...") or are missing
manual_lookup <- tribble(
  ~genus,                  ~family,                  ~order,
  "Bacteroidales",         NA,                       "Bacteroidales",
  "Clostridiales",         NA,                       "Clostridiales",
  "Clostridiaceae",        "Clostridiaceae",         "Clostridiales",
  "Lachnospiraceae",       "Lachnospiraceae",        "Clostridiales",
  "Peptostreptococcaceae", "Peptostreptococcaceae",  "Clostridiales",
  "Ruminococcaceae",       "Ruminococcaceae",        "Clostridiales",
  "Subdoligranulum",       "Ruminococcaceae",        "Clostridiales",
  "Erysipelotrichaceae",   "Erysipelotrichaceae",    "Erysipelotrichales",
  "Burkholderiales",       NA,                       "Burkholderiales",
  "Methanosphaera",        "Methanobacteriaceae",    "Methanobacteriales",
  "Propionibacterium",     "Propionibacteriaceae",   "Propionibacteriales")
fr_long <- fr_sp %>%
  select(sample_id, all_of(first_fr:ncol(fr_sp))) %>%
  pivot_longer(-sample_id, names_to = "species", values_to = "abundance") %>%
  mutate(genus = sub("_.*", "", species), species = str_replace_all(species, "_", " ")) %>%
  left_join(bind_rows(manual_lookup, genus_lookup %>% filter(!genus %in% manual_lookup$genus)) %>%
              distinct(genus, .keep_all = TRUE), by = "genus")

tax_long <- bind_rows(hall = hall_long, hmp2 = hmp2_long, franzosa = fr_long, .id = "cohort") %>%
  semi_join(meta_all, by = c("cohort", "sample_id")) %>%
  # harmonise taxonomy: new NCBI names -> old MetaPhlAn2 names
  mutate(order = recode(order, Eubacteriales = "Clostridiales"),
         genus_harm = recode(genus, Phocaeicola = "Bacteroides")) %>%   # Phocaeicola = formerly Bacteroides
  group_by(cohort, sample_id) %>%
  mutate(abundance = abundance / sum(abundance)) %>%                     # TSS at species level
  ungroup()

# how much abundance has no order/family? (should be small)
tax_long %>% group_by(cohort, sample_id) %>%
  summarise(no_order = sum(abundance[is.na(order)]), no_family = sum(abundance[is.na(family)]), .groups = "drop") %>%
  group_by(cohort) %>% summarise(median_no_order = median(no_order), max_no_order = max(no_order),
                                 median_no_family = median(no_family))
# which species have no order? (check of the assignment)
tax_long %>% filter(is.na(order)) %>% group_by(cohort, species) %>%
  summarise(mean_abund = mean(abundance), .groups = "drop") %>% arrange(desc(mean_abund)) %>% print(n = 20)

# ---- A2. Hypothesis taxa per sample -------------------------------------
# Gevers et al. 2014 (Cell Host Microbe 15:382-392), core version from abstract/text:
#   increased in CD (numerator):   families Enterobacteriaceae, Pasteurellaceae, Veillonellaceae, Fusobacteriaceae
#   decreased in CD (denominator): orders Bacteroidales, Clostridiales, Erysipelotrichales
# The complete list is in Table S2A of the paper -> still to be cross-checked.
# Veillonellaceae is its own order in MetaPhlAn2 and NCBI (Veillonellales/Selenomonadales),
# so it is not counted twice.
gevers_up   <- c("Enterobacteriaceae", "Pasteurellaceae", "Veillonellaceae", "Fusobacteriaceae")
gevers_down <- c("Bacteroidales", "Clostridiales", "Erysipelotrichales")

per_sample <- tax_long %>%
  group_by(cohort, sample_id) %>%
  summarise(
    Prevotella         = sum(abundance[genus_harm %in% "Prevotella"]),
    Bacteroides        = sum(abundance[genus_harm %in% "Bacteroides"]),
    F_prausnitzii      = sum(abundance[species %in% "Faecalibacterium prausnitzii"]),
    Roseburia          = sum(abundance[genus_harm %in% "Roseburia"]),
    Enterobacteriaceae = sum(abundance[family %in% "Enterobacteriaceae"]),
    Pasteurellaceae    = sum(abundance[family %in% "Pasteurellaceae"]),
    Veillonellaceae    = sum(abundance[family %in% "Veillonellaceae"]),
    Fusobacteriaceae   = sum(abundance[family %in% "Fusobacteriaceae"]),
    Bacteroidales      = sum(abundance[order %in% "Bacteroidales"]),
    Clostridiales      = sum(abundance[order %in% "Clostridiales"]),
    Erysipelotrichales = sum(abundance[order %in% "Erysipelotrichales"]),
    Gevers_up   = sum(abundance[family %in% gevers_up]),
    Gevers_down = sum(abundance[order %in% gevers_down]),
    .groups = "drop") %>%
  left_join(meta_all, by = c("cohort", "sample_id"))

target_taxa <- c("Prevotella", "Bacteroides", "F_prausnitzii", "Roseburia",
                 "Enterobacteriaceae", "Pasteurellaceae", "Veillonellaceae", "Fusobacteriaceae",
                 "Bacteroidales", "Clostridiales", "Erysipelotrichales", "Gevers_up", "Gevers_down")

# ---- A3. Detectability -------------------------------------------------
# per cohort x diagnosis: detectable in what % of samples / persons, and how much?
detect <- per_sample %>%
  pivot_longer(all_of(target_taxa), names_to = "taxon", values_to = "abundance") %>%
  group_by(cohort, taxon, diagnosis) %>%
  summarise(samples = n(), persons = n_distinct(subject_id),
            prevalence_samples = mean(abundance > 0),
            prevalence_persons = mean(tapply(abundance > 0, subject_id, any)),
            median_pct = 100 * median(abundance),
            mean_pct = 100 * mean(abundance),
            median_pct_if_present = 100 * median(abundance[abundance > 0]),
            .groups = "drop") %>%
  mutate(taxon = factor(taxon, levels = target_taxa), cohort = factor(cohort, levels = c("franzosa", "hmp2", "hall"))) %>%
  arrange(taxon, cohort, diagnosis)
print(detect, n = Inf, width = Inf)
write_csv(detect, file.path(outdir, "detectability_taxa.csv"))

# compact overview: prevalence in % (samples) per taxon x cohort/diagnosis
detect %>%
  transmute(taxon, group = paste(cohort, diagnosis), prevalence = round(100 * prevalence_samples)) %>%
  pivot_wider(names_from = group, values_from = prevalence) %>%
  print(n = Inf, width = Inf)

# distribution of the hypothesis taxa (log10, zeros as separate value at the bottom)
pc_plot <- 1e-6
p_taxa <- per_sample %>%
  pivot_longer(c(Prevotella, Bacteroides, F_prausnitzii, Roseburia), names_to = "taxon", values_to = "abundance") %>%
  mutate(cohort = factor(cohort, levels = c("franzosa", "hmp2", "hall"))) %>%
  ggplot(aes(diagnosis, pmax(abundance, pc_plot), colour = diagnosis)) +
  geom_jitter(width = .2, alpha = .35, size = .8) +
  geom_boxplot(fill = NA, outlier.shape = NA, colour = "black") +
  scale_y_log10() + facet_grid(taxon ~ cohort) +
  labs(y = "relative abundance (log10; zeros at 1e-6)", x = NULL,
       title = "Hypothesis taxa (SQ2), all samples, no prevalence filter") +
  theme(legend.position = "none")
p_taxa
ggsave(file.path(outdir, "hypotheses_taxa_distribution.png"), p_taxa, width = 9, height = 10)

# ---- A4. Prevotella/Bacteroides (SQ2a): how many zeros? ------------------
pb <- per_sample %>%
  group_by(cohort, diagnosis) %>%
  summarise(samples = n(),
            prevotella_null = mean(Prevotella == 0),
            bacteroides_null = mean(Bacteroides == 0),
            both_present = mean(Prevotella > 0 & Bacteroides > 0),
            # Prevotella-dominated (enterotype-like): Prevotella > Bacteroides
            prevotella_dominated = mean(Prevotella > Bacteroides),
            .groups = "drop")
print(pb)
write_csv(pb, file.path(outdir, "prevotella_bacteroides_zeros.csv"))

# distribution of the log ratio (only samples with both taxa) - usually bimodal (enterotypes)
p_pb <- per_sample %>% filter(Prevotella > 0, Bacteroides > 0) %>%
  mutate(cohort = factor(cohort, levels = c("franzosa", "hmp2", "hall"))) %>%
  ggplot(aes(log10(Prevotella / Bacteroides), fill = diagnosis)) +
  geom_histogram(bins = 40, alpha = .6, position = "identity") +
  facet_wrap(~ cohort, scales = "free_y") +
  labs(x = "log10(Prevotella / Bacteroides), only samples with both taxa", y = "samples")
p_pb
ggsave(file.path(outdir, "prevotella_bacteroides_logratio.png"), p_pb, width = 10, height = 4)

# ---- A5. Gevers index: components (SQ1) --------------------------------
# numerator = 0 -> log undefined; how often does that happen?
gevers_check <- per_sample %>%
  group_by(cohort, diagnosis) %>%
  summarise(samples = n(), numerator_zero = mean(Gevers_up == 0), denominator_zero = mean(Gevers_down == 0),
            median_up_pct = 100 * median(Gevers_up), median_down_pct = 100 * median(Gevers_down),
            .groups = "drop")
print(gevers_check)
write_csv(gevers_check, file.path(outdir, "gevers_components.csv"))

p_gev <- per_sample %>%
  pivot_longer(c(Enterobacteriaceae, Pasteurellaceae, Veillonellaceae, Fusobacteriaceae,
                 Bacteroidales, Clostridiales, Erysipelotrichales), names_to = "taxon", values_to = "abundance") %>%
  mutate(direction = ifelse(taxon %in% gevers_up, "increased in CD (numerator)", "decreased in CD (denominator)"),
         cohort = factor(cohort, levels = c("franzosa", "hmp2", "hall"))) %>%
  group_by(cohort, diagnosis, direction, taxon) %>%
  summarise(prevalence = 100 * mean(abundance > 0), .groups = "drop") %>%
  ggplot(aes(taxon, prevalence, fill = diagnosis)) +
  geom_col(position = "dodge") + facet_grid(cohort ~ direction, scales = "free_x", space = "free_x") +
  labs(y = "prevalence (% samples)", x = NULL, title = "Gevers taxa: detectable in how many samples?") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
p_gev
ggsave(file.path(outdir, "gevers_taxa_prevalence.png"), p_gev, width = 10, height = 7)

# ---- B. Measured butyrate (SQ3a, Franzosa only) --------------------------
# source: Franzosa et al. 2019, LC-MS (HILIC-neg). According to the lookup table:
#   HILIC-neg_Cluster_0013 = "butyrate / isobutyrate*" -> butyrate and isobutyrate not separated
#   values = peak intensities (relative), NOT absolute concentrations
read.csv(file.path(path, "processed/franzosa/franzosa_scfa_lookup.csv"))

scfa <- read_csv(file.path(path, "processed/franzosa/franzosa_h3a_merged.csv"), show_col_types = FALSE) %>%
  select(sample_id, butyrate, propionate) %>%
  inner_join(meta_all %>% filter(cohort == "franzosa"), by = "sample_id")

scfa %>% summarise(samples = n(), butyrate_NA = sum(is.na(butyrate)), butyrate_zero = sum(butyrate == 0, na.rm = TRUE),
                   propionate_NA = sum(is.na(propionate)), propionate_zero = sum(propionate == 0, na.rm = TRUE))
scfa %>% group_by(site, diagnosis) %>%
  summarise(n = n(), median_butyrate = median(butyrate, na.rm = TRUE),
            q1 = quantile(butyrate, .25, na.rm = TRUE), q3 = quantile(butyrate, .75, na.rm = TRUE),
            median_propionat = median(propionate, na.rm = TRUE), .groups = "drop")

# right-skewed? -> look at log10
p_but <- scfa %>%
  pivot_longer(c(butyrate, propionate), names_to = "scfa", values_to = "value") %>%
  ggplot(aes(diagnosis, value, colour = diagnosis, shape = site)) +
  geom_jitter(width = .2, alpha = .6) + geom_boxplot(fill = NA, outlier.shape = NA, colour = "black") +
  scale_y_log10() + facet_wrap(~ scfa, scales = "free_y") +
  labs(y = "LC-MS intensity (log10)", x = NULL, title = "Franzosa: measured SCFA (relative)") +
  theme(legend.position = "bottom")
p_but
ggsave(file.path(outdir, "franzosa_scfa.png"), p_but, width = 8, height = 5)

# associations (descriptive): butyrate vs. propionate, calprotectin, antibiotics
cor.test(scfa$butyrate, scfa$propionate, method = "spearman", exact = FALSE)
cor.test(scfa$butyrate, scfa$calprotectin, method = "spearman", exact = FALSE)
scfa %>% group_by(diagnosis, antibiotic) %>% summarise(n = n(), median_butyrate = median(butyrate, na.rm = TRUE), .groups = "drop")

# butyrate vs. hypothesis taxa (preview of SQ3: "accompanied by")
scfa_taxa <- scfa %>% left_join(per_sample %>% filter(cohort == "franzosa") %>%
                                  select(sample_id, F_prausnitzii, Roseburia, Prevotella, Bacteroides), by = "sample_id")
map_dfr(c("F_prausnitzii", "Roseburia", "Prevotella", "Bacteroides"), function(tx) {
  ct <- cor.test(scfa_taxa$butyrate, scfa_taxa[[tx]], method = "spearman", exact = FALSE)
  tibble(taxon = tx, rho = ct$estimate, p = ct$p.value)
})

# ---- C. Functional potential (SQ3b) --------------------------------------
# HMP2 (HUMAnN): 3 butyrate ECs and 3 butanoate pathways
hmp2_ec <- read_csv(file.path(path, "processed/hmp2/hmp2_ec_butyrate_merged.csv"), show_col_types = FALSE)
hmp2_pw <- read_csv(file.path(path, "processed/hmp2/hmp2_pathway_butanoate_merged.csv"), show_col_types = FALSE)
ec_cols <- grep("^[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+:", colnames(hmp2_ec), value = TRUE)
pw_cols <- grep("PWY", colnames(hmp2_pw), value = TRUE)
ec_cols; pw_cols

func_hmp2 <- hmp2_ec %>% select(sample_id, all_of(ec_cols)) %>%
  inner_join(hmp2_pw %>% select(sample_id, all_of(pw_cols)), by = "sample_id") %>%
  pivot_longer(-sample_id, names_to = "feature", values_to = "value") %>%
  inner_join(meta_all %>% filter(cohort == "hmp2"), by = "sample_id")

# do the samples match the abundance data?
n_distinct(func_hmp2$sample_id); sum(meta_all$cohort == "hmp2")

func_hmp2 %>% group_by(feature, diagnosis) %>%
  summarise(samples = n(), share_zero = mean(value == 0), median = median(value),
            q1 = quantile(value, .25), q3 = quantile(value, .75), .groups = "drop") %>%
  print(n = Inf, width = Inf)

# per person (median over time), so that persons with many samples do not dominate
p_func_hmp2 <- func_hmp2 %>%
  group_by(feature, subject_id, diagnosis, site) %>% summarise(value = median(value), .groups = "drop") %>%
  ggplot(aes(diagnosis, value + 1, colour = diagnosis, shape = site)) +
  geom_jitter(width = .2, alpha = .7) + geom_boxplot(fill = NA, outlier.shape = NA, colour = "black") +
  scale_y_log10() + facet_wrap(~ feature, scales = "free_y", ncol = 3,
                               labeller = labeller(feature = label_wrap_gen(35))) +
  labs(y = "abundance + 1 (log10), median per person", x = NULL,
       title = "HMP2: butyrate ECs and butanoate pathways") +
  theme(legend.position = "bottom")
p_func_hmp2
ggsave(file.path(outdir, "hmp2_butyrate_function.png"), p_func_hmp2, width = 10, height = 7)

# are the functions associated with the taxa? (Spearman, all samples, descriptive)
func_taxa <- func_hmp2 %>% select(sample_id, feature, value) %>%
  pivot_wider(names_from = feature, values_from = value) %>%
  inner_join(per_sample %>% filter(cohort == "hmp2") %>% select(sample_id, F_prausnitzii, Roseburia), by = "sample_id")
cor(func_taxa %>% select(-sample_id), method = "spearman") %>% round(2)

# Franzosa: enzyme table (ECs) -> find butyrate-relevant ECs
fr_enz <- read_csv(file.path(path, "processed/franzosa/franzosa_enz_merged.csv"), show_col_types = FALSE)
first_ec <- which(colnames(fr_enz) == "steroids") + 1
enz_cols <- colnames(fr_enz)[first_ec:ncol(fr_enz)]
# check unit: do the ECs per sample sum to 1?
summary(rowSums(fr_enz[, enz_cols]))
# the same 3 ECs as for HMP2 + butyryl-CoA dehydrogenase (1.3.8.1)
butyrate_ec <- c("X2.3.1.19", "X2.7.2.7", "X2.8.3.8", "X1.3.8.1")
grep(paste0("^(", paste(gsub(".", "\\.", butyrate_ec, fixed = TRUE), collapse = "|"), ")\\.\\."), enz_cols, value = TRUE)

fr_ec_cols <- grep(paste0("^(", paste(gsub(".", "\\.", butyrate_ec, fixed = TRUE), collapse = "|"), ")\\.\\."), enz_cols, value = TRUE)
func_fr <- fr_enz %>% select(sample_id, all_of(fr_ec_cols)) %>%
  pivot_longer(-sample_id, names_to = "feature", values_to = "value") %>%
  inner_join(meta_all %>% filter(cohort == "franzosa"), by = "sample_id")
func_fr %>% group_by(feature, diagnosis) %>%
  summarise(samples = n(), share_zero = mean(value == 0), median = median(value), .groups = "drop") %>%
  print(n = Inf, width = Inf)

p_func_fr <- func_fr %>%
  ggplot(aes(diagnosis, value, colour = diagnosis, shape = site)) +
  geom_jitter(width = .2, alpha = .6) + geom_boxplot(fill = NA, outlier.shape = NA, colour = "black") +
  scale_y_log10() + facet_wrap(~ feature, scales = "free_y", labeller = labeller(feature = label_wrap_gen(35))) +
  labs(y = "relative abundance (log10)", x = NULL, title = "Franzosa: butyrate-relevant ECs") +
  theme(legend.position = "bottom")
p_func_fr
ggsave(file.path(outdir, "franzosa_butyrate_ecs.png"), p_func_fr, width = 9, height = 7)

# gene potential vs. measured butyrate (Franzosa only, descriptive)
func_fr %>% select(sample_id, feature, value) %>%
  inner_join(scfa %>% select(sample_id, butyrate), by = "sample_id") %>%
  group_by(feature) %>%
  summarise(rho = cor(value, butyrate, method = "spearman", use = "complete.obs"),
            p = cor.test(value, butyrate, method = "spearman", exact = FALSE)$p.value)

saveRDS(list(per_sample = per_sample, detect = detect, pb = pb, gevers_check = gevers_check,
             scfa = scfa, func_hmp2 = func_hmp2, func_fr = func_fr),
        file.path(outdir, "hypotheses_eda.rds"))
