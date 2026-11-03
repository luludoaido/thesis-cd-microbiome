# From Taxa to Function: Diet-Related Microbiome Signatures in Crohn's Disease

Master's thesis, Luana da Silva do Aido, ACLS (ZHAW), supervisor: Dr. Tugçe Bilgin Sonay.
Target submission 07.01.2027 (official deadline 11.01.2027).

## About the study

### What this study is about

Crohn's disease (CD) is a chronic inflammatory bowel disease. In people with CD the composition of the gut bacteria is often altered, which is called **dysbiosis** (Gevers et al. 2014; Lloyd-Price et al. 2019). Bacteria that ferment dietary fibre into **short-chain fatty acids** (SCFA), above all **butyrate**, are particularly affected. Butyrate is the main energy source of the colonic mucosa (Roediger 1980). The butyrate producer *Faecalibacterium prausnitzii* is anti-inflammatory and less abundant in CD (Sokol et al. 2008).

This thesis compares CD patients with healthy controls in three public datasets: Franzosa et al. 2019, Hall et al. 2017 and HMP2 (Lloyd-Price et al. 2019). It proceeds in three steps:

1. **Is the microbiome altered?** How strongly does the bacterial composition in CD differ from healthy people (dysbiosis)?
2. **Which bacteria are missing?** Are there fewer fibre-associated bacteria and butyrate producers such as *F. prausnitzii* and *Roseburia*?
3. **Is their function missing too?** Do the bacteria carry fewer genes for butyrate production, and is less butyrate measurable in stool?

Hence the title "From Taxa to Function": from the question of *who* is there to the question of *what* the bacteria do.

### Why it matters

- **Frequency:** CD and other inflammatory bowel diseases are now found worldwide, including in countries where they used to be rare (Ng et al. 2017). Diet and the microbiome are considered important environmental factors.
- **Function, not just names:** Many studies only count which bacteria are missing. For the gut, however, what matters is which substances they produce. In CD, the function of the microbiome changes in part more clearly than its composition (Morgan et al. 2012). This thesis therefore looks directly at butyrate genes and measured butyrate.
- **Starting point for diet:** If CD mainly lacks the capacity to use fibre, this is a possible starting point for diet or therapy.
- **Several cohorts:** A finding is more credible if it appears in several independent datasets. The thesis therefore checks whether the differences point in the same direction in all three cohorts.
- **Data availability:** There are hardly any public datasets in which both the microbiome and SCFA were measured in CD. This is itself a result of the thesis and shows where future studies should start.

Diet itself is not recorded in the data. Fibre-associated bacteria and their products serve as an indirect measure (Wu et al. 2011). Propionate is less specific for this than butyrate, because it is also produced from protein and other sources (Reichardt et al. 2014).

### Research question and hypotheses

Main research question: *How do fiber-associated bacteria and their metabolic function differ between CD patients and healthy controls?*

**SQ1 – Dysbiosis.** Does the dataset used reproduce the dysbiosis pattern known from the literature for Crohn's disease, and are the taxa required for this reliably detectable in the dataset at all?

- **H1:** CD patients differ significantly from healthy controls both in the taxon-based dysbiosis index according to Gevers et al. 2014 and in a distance-based dysbiosis measure (Bray-Curtis distance).

**SQ2 – Fibre-associated bacteria.** Does the ratio of fiber-associated, potentially SCFA-producing bacteria differ between CD patients and controls?

- **H2a:** The *Prevotella*/*Bacteroides* ratio is significantly lower in CD patients than in controls.
- **H2b:** The relative abundance of *Faecalibacterium prausnitzii* and *Roseburia* spp. is significantly reduced in CD patients compared to controls.

**SQ3 – Butyrate metabolism.** Are compositional differences in fiber-associated bacteria between CD patients and controls accompanied by corresponding differences in SCFA/butyrate metabolism?

- **H3a:** The measured butyrate concentration is significantly lower in CD patients than in controls.
- **H3b:** The predicted functional potential for SCFA/butyrate biosynthesis pathways is significantly reduced in CD patients compared with controls.

Changes since the original wording (details in [`docs/analysis_log.md`](docs/analysis_log.md)):

- **H3b is the main analysis of SQ3**, H3a is exploratory (D10). Reason: measured butyrate is only available in one cohort, the enzymes in both HMP2 and Franzosa.
- **H3a also covers propionate** (meeting with Tugçe). Butyrate and propionate are relative peak intensities, not concentrations (D05).
- **"predicted" in H3b** is inaccurate: the enzymes and pathways are gene potential measured from metagenomes (HUMAnN), not a prediction from 16S data. Better: "metagenomic functional potential".

### Key terms

| Term | Meaning |
|---|---|
| Shotgun metagenome | All DNA in a stool sample is sequenced. This shows which bacteria are present and which genes they carry. |
| Relative abundance | Share of one bacterium among all bacteria in a sample (sums to 100 %). |
| Species / genus | The analysis runs mainly at genus level. |
| Alpha diversity | Diversity within one sample (richness = number of genera, Shannon = number and evenness). |
| Beta diversity | How different two samples are (e.g. Bray-Curtis distance), tested with PERMANOVA. |
| MDI | Microbial Dysbiosis Index: ratio of bacteria that increase in CD to those that decrease. |
| EC / pathway | EC number = one enzyme. Pathway = a chain of enzymes that together produce, for example, butyrate. |
| Metabolomics | Measurement of metabolic products in stool, here butyrate and propionate (as peak intensity, not concentration). |
| Longitudinal | Hall and HMP2 have several samples per person over time. These are not independent, hence mixed models with person as random effect. |
| Batch effect | Differences that come only from the centre or lab, not from the disease. |

### Cohorts

| Cohort | Samples | Persons (CD / control) | Available data | Processing |
|---|---|---|---|---|
| Franzosa 2019 | 144 | 88 / 56 (1 sample per person) | species, enzymes, metabolites | bioBakery 2 (MetaPhlAn 2, HUMAnN 2) |
| Hall 2017 | 158 | 9 / 12 (several samples per person) | species, genus | bioBakery 3 (MetaPhlAn 3), via curatedMetagenomicData |
| HMP2 (IBDMDB) | 571 | 30 / 14 (several samples per person) | species, genus, enzymes, pathways | bioBakery 3 (MetaPhlAn 3, HUMAnN 3) |

Which samples are included and why is documented in [`docs/decisions.md`](docs/decisions.md).

## Folder structure

```
thesis-cd-microbiome/
├── README.md                 this file
├── thesis-cd-microbiome.Rproj   open in RStudio
├── renv.lock                 list of all R packages with version
├── .Rprofile                 starts renv automatically
├── renv/                     renv settings (packages themselves not in Git)
│
├── data/                     not in Git
│   ├── raw/                  original data, never modify
│   │   ├── cmd/              Hall 2017
│   │   ├── hmp2/             HMP2
│   │   └── franzosa/         Franzosa 2019
│   └── processed/            cleaned data, created by 02
│       ├── franzosa/
│       ├── hall/
│       └── hmp2/
│
├── R/                        all scripts, numbered in run order
│   ├── 00_setup.R … 09_hypotheses_eda.R
│   └── checks/               checks of individual anomalies
│
├── results/                  plots (.png) and tables (.csv), in Git
│   ├── franzosa/  hall/  hmp2/   main analysis per cohort (06)
│   ├── species/              as above, at species level (only created by `06 … species`)
│   ├── hmp2_clinical/        clinical values HMP2 (05)
│   ├── extremes/             most/least abundant genera (07)
│   ├── robustness/           robustness checks (08)
│   └── hypotheses_eda/       EDA for the hypotheses (09)
│
└── docs/                     documentation, see below
```

## Data

The data are not in the repo because they are too large (the HMP2 enzyme table alone is 700 MB). All data are public. Download them and place them in these folders with exactly these file names. The folders already exist (empty).

### Hall 2017 → `data/raw/cmd/`

No manual download. The script `R/01_download_cmd.R` fetches the data from the Bioconductor package curatedMetagenomicData (version 2021-03-31) and saves them as:

- `HallAB_2017_taxa.csv`
- `HallAB_2017_metadata.csv`
- `HallAB_2017_taxonomy.csv`

### HMP2 (IBDMDB) → `data/raw/hmp2/`

Source: [ibdmdb.org](https://ibdmdb.org), Inflammatory Bowel Disease Multi'omics Database. Open the HMP2 data under Downloads.

| File | Content | Where on ibdmdb.org |
|---|---|---|
| `hmp2_metadata_2018-08-20.csv` | metadata for all samples (diagnosis, age, centre, reads, …) | HMP2 → Metadata |
| `taxonomic_profiles_3.tsv` | taxa (MetaPhlAn 3) | HMP2 → Metagenomes → bioBakery 3 products |
| `pathabundances_3.tsv` | pathways (HUMAnN 3) | HMP2 → Metagenomes → bioBakery 3 products |
| `ecs_3.tsv` | enzymes (HUMAnN 3), approx. 700 MB | HMP2 → Metagenomes → bioBakery 3 products |
| `hmp2_dysbiosis_scores.tsv` | dysbiosis scores from Lloyd-Price 2019 (currently not used by the scripts) | HMP2 → Metagenomes |

Note: `hmp2_metadata_2018-08-20.csv` is stored locally with `;` as separator. `R/02_clean_hmp2.R` reads it with `sep = ";"`. A freshly downloaded file uses `,`; in that case this line must be adjusted.

### Franzosa 2019 → `data/raw/franzosa/`

Source: Franzosa et al. 2019, "Gut microbiome structure and metabolic activity in inflammatory bowel disease", *Nature Microbiology* 4, 293–305. On the article page, download the Supplementary Data under Supplementary Information and rename them as follows:

| File | Content |
|---|---|
| `SupplDataset1_metabolite_info.xlsx` | Supplementary Data 1: description of the metabolite clusters |
| `SupplDataset2_metabolites.xlsx` | Supplementary Data 2: metabolite values per sample |
| `SupplDataset4_species.xlsx` | Supplementary Data 4: species per sample |
| `SupplDataset6_enzymes.xlsx` | Supplementary Data 6: enzymes (EC) per sample |

The raw sequences (105–352 GB) are not needed; the tables are enough.

## Setup

Goal: on every computer the scripts run with exactly the same package versions as on mine. **renv** takes care of this.

**How renv works:** Normally R installs packages once for the whole computer, in whatever version is current. renv instead creates a separate package library for this project (in the folder `renv/library/`). Which packages belong in it, in which version, is listed in `renv.lock`. The file `.Rprofile` makes R use this project library automatically when the project is opened.

### One-time: set up the project

You only do these steps once per computer.

1. **Install R 4.6** (from [cran.r-project.org](https://cran.r-project.org)) and RStudio (recommended).
2. **Download the repo.** In the terminal:
   ```bash
   git clone https://github.com/luludoaido/thesis-cd-microbiome.git
   ```
   Alternatively on GitHub: green button "Code" → "Download ZIP" and unzip.
3. **Open the project:** double-click `thesis-cd-microbiome.Rproj`. RStudio opens in the project folder. The console shows a message from renv that packages are missing. This is normal the first time.
4. **Install packages.** In the RStudio console:
   ```r
   renv::restore()
   ```
   renv shows a list of packages and asks whether to install them. Confirm with `y`. The first time this takes 10–30 minutes, because more than 300 packages are installed, including the Bioconductor packages (curatedMetagenomicData, sva, limma, MMUPHin).
   If you get a message that renv is missing: first `install.packages("renv")`, then restart RStudio and run `renv::restore()` again.
5. **Add the data** as described in [Data](#data).
6. **Test:** in the console, `source("R/00_setup.R")`. If no error appears, everything is ready.

### Day to day

- **Open the project:** always via `thesis-cd-microbiome.Rproj`. Then the working directory and packages are set automatically.
- **Run a script:** open it in RStudio and click "Source" at the top right. Or in the terminal from the project folder, e.g. `Rscript R/06_diversity.R hmp2 genus`.
- **Need a new package?** Install it with `install.packages("name")`, then run `renv::snapshot()` so it is added to `renv.lock`, and commit `renv.lock`.

Paths in the code are relative to the project folder, so it runs on any computer. `R/00_setup.R` stops with a message if the working directory is wrong.

## Scripts

### How the scripts fit together

The scripts build on each other: each reads files written by an earlier one.

```
01  download Hall             → data/raw/cmd/
02  clean (3 scripts)         → data/processed/<cohort>/*_meta.csv and *_merged.csv
04–09  analyses               → results/
```

There is no script 03: merging the metadata with the feature tables is now part of the 02 scripts.

**You never need to run `00_setup.R` yourself.** Every script has the line `source("R/00_setup.R")` at the top. `source()` runs another R file. So when you start a script, it first loads packages, paths, exclusion rules, seed and colours from `00_setup.R` automatically, and only then computes.

### When to run which script

| Situation | What to run |
|---|---|
| First time (new computer, fresh clone) | all in order: `01` → `02` (all three) → `04` to `09` |
| Continue working on one analysis | only that script, e.g. only `06` |
| Changed an early script or a rule in `00_setup.R` | that script **and all after it**, otherwise later scripts compute with old data |

### Overview

| Script | Purpose | Reads | Writes |
|---|---|---|---|
| `R/00_setup.R` | paths, exclusion rules, seed, colours, plot theme. Loaded by all scripts. | – | – |
| `R/01_download_cmd.R` | download Hall 2017 from curatedMetagenomicData | internet | `data/raw/cmd/` |
| `R/02_clean_franzosa.R` | Franzosa: reshape supplement tables, metadata without UC, join metadata and features | `data/raw/franzosa/` | `data/processed/franzosa/` |
| `R/02_clean_hall.R` | Hall: sum species to genus, metadata CD + controls, join metadata and features | `data/raw/cmd/` | `data/processed/hall/` |
| `R/02_clean_hmp2.R` | HMP2: metadata (adults, CD + controls), taxa, butyrate pathways and enzymes, join metadata and features | `data/raw/hmp2/` | `data/processed/hmp2/` |
| `R/04_eda_metadata.R` | EDA of the metadata: samples, diagnosis, age, NAs, calprotectin, HBI, medication | `data/processed/` | console/plots only |
| `R/05_hmp2_clinical.R` | HMP2: match HBI, calprotectin and antibiotics to the metagenome samples | `data/` | `results/hmp2_clinical/` |
| `R/06_diversity.R` | main analysis per cohort: sparsity, alpha diversity, mixed models, time trend, PERMANOVA, PCA/UMAP, batch | `data/processed/` | `results/<cohort>/` |
| `R/07_extremes.R` | most and least abundant genera, samples with extreme diversity | `data/processed/` | `results/extremes/` |
| `R/08_robustness.R` | robustness checks: clinical scores, Franzosa outliers, sample size, time stratification, Boston power | `results/` (from 05 and 06) | `results/robustness/` |
| `R/09_hypotheses_eda.R` | EDA for the hypotheses: detectability of taxa, butyrate (measured and functional) | `data/` | `results/hypotheses_eda/` |

`R/06_diversity.R` runs once per cohort. In RStudio, set `cohort` and `level` at the top of the script (lines 22–23) and click "Source", then the next cohort. In the terminal:

```bash
Rscript R/06_diversity.R franzosa genus
```

```bash
Rscript R/06_diversity.R hall genus
```

```bash
Rscript R/06_diversity.R hmp2 genus
```

With `species` instead of `genus` the analysis runs at species level. Results then go to `results/species/`.

Further scripts:

| Script | Purpose |
|---|---|
| `R/checks/check_hall_samples.R` | inspect conspicuous Hall samples (composition, reads, other samples of the same person) |
| `R/checks/prevalence_filter_effect.R` | how many genera the prevalence filter removes per cohort (needs `R/07_extremes.R`) |

Reproducibility: the seed is fixed (`seed <- 42` in `R/00_setup.R`). Two runs of all scripts produce identical files.

## Documentation (`docs/`)

| File | Purpose | When to write |
|---|---|---|
| [`docs/analysis_log.md`](docs/analysis_log.md) | The project diary: timeline, all decisions with number (D01, D02, …) and rationale, one entry per week (section 4), open items and limitations for the thesis. | One entry in section 4 at the end of each week; decisions immediately. |
| [`docs/decisions.md`](docs/decisions.md) | The technical rules in the code: which samples are excluded and why, with numbers before/after and where in the code. Refers to the D numbers in `analysis_log.md`. | As soon as a rule in the code changes. |
| [`docs/analysis_plan.md`](docs/analysis_plan.md) | Analysis plan: which tests with which models, fixed **before** looking at the results. Protects against choosing the method afterwards that gives the nicest result. (currently empty) | Before each new analysis. Add deviations with date and reason. |

In short: `analysis_log.md` says **what and why** was decided, `decisions.md` says **how it is implemented in the code**.

## References

- Franzosa EA et al. (2019). Gut microbiome structure and metabolic activity in inflammatory bowel disease. *Nature Microbiology* 4, 293–305.
- Gevers D et al. (2014). The treatment-naive microbiome in new-onset Crohn's disease. *Cell Host & Microbe* 15(3), 382–392.
- Hall AB et al. (2017). A novel *Ruminococcus gnavus* clade enriched in inflammatory bowel disease patients. *Genome Medicine* 9, 103.
- Lloyd-Price J et al. (2019). Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. *Nature* 569, 655–662.
- Morgan XC et al. (2012). Dysfunction of the intestinal microbiome in inflammatory bowel disease and treatment. *Genome Biology* 13, R79.
- Ng SC et al. (2017). Worldwide incidence and prevalence of inflammatory bowel disease in the 21st century: a systematic review of population-based studies. *The Lancet* 390, 2769–2778.
- Reichardt N et al. (2014). Phylogenetic distribution of three pathways for propionate production within the human gut microbiota. *The ISME Journal* 8, 1323–1335.
- Roediger WEW (1980). Role of anaerobic bacteria in the metabolic welfare of the colonic mucosa in man. *Gut* 21, 793–798.
- Sokol H et al. (2008). *Faecalibacterium prausnitzii* is an anti-inflammatory commensal bacterium identified by gut microbiota analysis of Crohn disease patients. *PNAS* 105(43), 16731–16736.
- Wu GD et al. (2011). Linking long-term dietary patterns with gut microbial enterotypes. *Science* 334, 105–108.
