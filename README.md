# Colorectal–Hippo–Dysbiosis (Genetic Codon)

**Goal:** Integrate bulk RNA-seq from colorectal cancer (TCGA-CRC) and an IBD dysbiosis cohort (GSE235236) to examine how Hippo/YAP signaling in CRC and IBD shares transcriptional programs, and explore dysbiosis-driven Hippo/YAP–TAZ activation as a mechanistic link between dysbiosis and CRC.

**Organization:** Genetic Codon  
**PI:** Sana Noor  
**Sprint:** 2 months  
**Team:** A (Data Integration) | B (Analysis & Reporting)

## Cohorts
- **TCGA-CRC:** bulk RNA-seq (counts/expression) + clinical metadata.
- **GSE235236 (IBD):** colon mucosa expression (UC/CD vs controls).

## Specific Aims
1. Quantify Hippo/YAP activity in CRC tumor vs normal and IBD vs control; define a Dysbiosis–Hippo Response (DHR) module.
2. Compare Hippo/DHR-related differential signals between IBD and CRC.
3. Generate candidate mechanistic links between dysbiosis signatures and Hippo pathway activation for follow-up.

## Notebooks (detailed)

### `notebooks/00_env_and_utils.ipynb`
**Purpose:** Environment bootstrapping and utility functions setup.

**What it does:**
- Mounts Google Drive (if using Colab)
- Installs required packages (pandas, scipy, statsmodels, matplotlib, seaborn, etc.)
- Creates the folder tree for data and results
- Defines shared utility functions for:
  - Expression matrix normalization & QC
  - DEG statistical functions (t-tests, FDR correction)
  - Plotting themes and standard visualizations
  - Common gene set operations

**Output:** Utility module and environment ready for downstream notebooks. **Run this first.**

---

### `notebooks/01A_IBD_intake_qc_deg.ipynb`
**Owner:** Sai  
**Purpose:** IBD cohort (GSE235236) intake, quality control, normalization, and differential expression analysis.

**What it does:**
1. **Data Loading & Cleaning:**
   - Loads GSE235236 TPM matrix and sample metadata
   - Identifies and strips gene identifier column (symbol/Hugo symbol)
   - Detects and removes duplicate genes (by averaging)
   - Imputes missing values and logs normality check

2. **Quality Control:**
   - Plots sum(TPM) per sample to detect outliers/low-quality samples
   - Validates sample metadata alignment with expression data
   - Generates sample count statistics by group (UC, CD, Control)

3. **Data Transformation:**
   - Applies log₂(TPM + 1) transformation for downstream analysis
   - Aligns expression matrix with cleaned metadata (exact sample ID matching)

4. **Differential Expression Analysis:**
   - Compares UC vs Control using Welch's t-tests
   - Compares CD vs Control using Welch's t-tests
   - Corrects p-values for multiple testing (FDR-BH)
   - Calculates log₂FC and mean expression per group

5. **Visualization:**
   - Volcano plots for UC vs Control and CD vs Control
   - Displays fold-change thresholds (±log₂ 1.0) and FDR threshold (0.05)

**Outputs:**
- `data_processed/ibd/ibd_tpm_clean.csv` — QC-cleaned TPM matrix
- `data_processed/ibd/ibd_log2_tpm.csv` — log₂-transformed expression matrix
- `data_processed/ibd/ibd_metadata_clean.csv` — aligned sample metadata
- `results/de/ibd/IBD_DEG_UC_vs_Control.csv` — DEG table (UC vs healthy)
- `results/de/ibd/IBD_DEG_CD_vs_Control.csv` — DEG table (CD vs healthy)
- `results/figures/ibd/volcano_*.png` — volcano plots

---

### `notebooks/01B_TCGA_Analysis.ipynb`
**Purpose:** TCGA-CRC cohort intake, QC, and Hippo/YAP pathway scoring.

**What it does:**
1. **TCGA Data Loading & QC:**
   - Loads TCGA-COAD expression counts and clinical metadata
   - Filters for CRC samples; filters low-abundance genes
   - Normalizes expression (e.g., TMM, CPM, or log-cpm)
   - Plots sequencing depth, gene detection rate for QC

2. **Hippo/YAP Pathway Scoring:**
   - Loads Hippo pathway and YAP target gene signatures from config
   - Computes single-sample pathway activity scores (e.g., mean z-score of signature genes)
   - Stratifies samples into Hippo-High vs Hippo-Low groups

3. **Differential Expression & Validation:**
   - DEG analysis: Hippo-High vs Hippo-Low tumors
   - Validates Hippo/YAP signature genes are significantly enriched in Hippo-High group
   - Performs GSEA or pathway enrichment on the DEG list

4. **Visualization:**
   - Hippo score distributions by tumor grade, stage, mutation status
   - Volcano plot: Hippo-High vs Hippo-Low DEGs
   - Heatmap of top Hippo/YAP target genes

**Outputs:**
- `data_processed/tcga/tcga_log2_cpm.csv` — normalized expression matrix
- `data_processed/tcga/hippo_scores.csv` — per-sample Hippo pathway scores
- `results/de/tcga/TCGA_DEG_HippoHigh_vs_Low.csv` — DEG table
- `results/hippo/hippo_enrichment.csv` — pathway enrichment results
- `results/figures/tcga/volcano_hippo.png` — volcano plot
- `results/figures/tcga/heatmap_hippo_targets.png` — target gene heatmap

---

### `notebooks/02_Overlap between IBD and CRC based on Hippo pathway.ipynb`
**Purpose:** Cross-cohort comparison of Hippo/DHR pathway signals and overlap between IBD and CRC.

**What it does:**
1. **Load & Integrate:**
   - Loads IBD DEGs (UC/CD vs Control)
   - Loads TCGA DEGs (Hippo-High vs Low)
   - Loads Hippo/YAP and Dysbiosis–Hippo Response (DHR) gene signatures

2. **Pathway Activity in Both Cohorts:**
   - Computes Hippo and DHR scores in IBD samples
   - Compares Hippo scores: UC/CD vs Control (should be elevated in disease)
   - Compares Hippo scores: TCGA Hippo-High vs Low tumors

3. **Signature Overlap Analysis:**
   - Identifies genes differentially expressed in *both* IBD and CRC along Hippo pathway
   - Venn diagram: IBD UC DEGs ∩ TCGA Hippo-High DEGs
   - Venn diagram: IBD CD DEGs ∩ TCGA Hippo-High DEGs
   - Computes Jaccard indices and statistical significance of overlaps

4. **Correlation & Association:**
   - Correlates Hippo scores with disease severity (IBD) or tumor characteristics (TCGA)
   - Tests whether high DHR activity in IBD predicts Hippo activation patterns

5. **Visualization:**
   - Scatter plots: Hippo scores (IBD vs TCGA)
   - Venn diagrams: DEG overlaps
   - Heatmap of commonly dysregulated Hippo/DHR genes
   - Bar plots: Hippo score comparisons (with p-values)

**Outputs:**
- `results/hippo/ibd_hippo_scores.csv` — IBD Hippo pathway scores
- `results/hippo/overlap_genes.csv` — shared DEGs across cohorts
- `results/hippo/overlap_statistics.csv` — Jaccard, Fisher exact, overlap p-values
- `results/figures/venn_ibd_tcga_hippo.png` — Venn diagram
- `results/figures/hippo_scores_comparison.png` — score distributions & comparisons
- `results/figures/heatmap_shared_genes.png` — heatmap of overlapping genes

---

### `notebooks/03_Mechansitic Bridge b_w CRC and Dysbiosis.ipynb`
**Purpose:** Exploratory mechanistic analyses linking dysbiosis (microbe signatures) and Hippo pathway activation in CRC.

**What it does:**
1. **Dysbiosis Signature Loading:**
   - Loads microbial dysbiosis signatures (enriched bacterial taxa/functional modules from GSE235236 or external DB)
   - Computes dysbiosis scores in both IBD and TCGA samples

2. **Dysbiosis–Hippo Correlation:**
   - Tests correlation between dysbiosis and Hippo pathway scores across cohorts
   - Stratifies by dysbiosis level and compares Hippo scores (low dysb vs high dysb)

3. **Mechanistic Integration:**
   - Identifies genes co-regulated in dysbiosis signature genes and Hippo targets
   - Explores possible causal mechanisms:
     - Genes elevated in dysbiosis-high tumors AND Hippo-high tumors (convergence)
     - Candidate TFs linking dysbiosis signals to Hippo activation

4. **Pathway Integration Analysis:**
   - Combines IBD pathway perturbations (from DEG analysis) with CRC Hippo activation
   - Tests hypothesis: IBD-dysbiosis activates Hippo pathway → promotes CRC via aberrant Hippo signaling

5. **Visualization:**
   - Scatter plot: dysbiosis score vs Hippo score (all samples, colored by cohort)
   - Box plots: Hippo scores stratified by dysbiosis level
   - Sankey or network plot: dysbiosis → Hippo targets → downstream effects
   - Heatmap: expression of bridging genes (dysbiosis AND Hippo-regulated)

**Outputs:**
- `results/microbe/dysbiosis_scores.csv` — dysbiosis pathway scores (IBD + TCGA)
- `results/microbe/dysbiosis_hippo_correlation.csv` — correlation & p-values
- `results/hippo/mechanistic_genes.csv` — bridging genes (dysbiosis + Hippo)
- `results/figures/dysbiosis_hippo_scatter.png` — correlation scatter
- `results/figures/dysbiosis_stratified_hippo.png` — box plots by dysbiosis level
- `results/figures/mechanistic_bridge_network.png` — pathway integration visualization
- `results/tables/mechanistic_summary.csv` — candidate genes & mechanisms

---

## Getting Started (Colab)

1. **Clone or download this repo** to your Google Drive:
   ```
   GeneticCodonShared/Hippo_Dysbiosis/
   ```

2. **Prepare raw data** in `data_raw/`:
   - Place `GSE235236_TPM_matrix_with_symbols.csv` in `data_raw/ibd/`
   - Place `ibd_metadata.csv` in `data_raw/ibd/`
   - Place TCGA files in `data_raw/tcga/`

3. **Fill manifest files** (`config/manifest_tcga.csv`, `config/manifest_ibd.csv`) with exact paths to your datasets.

4. **Run in recommended order:**
   - `notebooks/00_env_and_utils.ipynb` (environment & utilities)
   - `notebooks/01A_IBD_intake_qc_deg.ipynb` (IBD intake, QC, DEG, Hippo/DHR scores)
   - `notebooks/01B_TCGA_Analysis.ipynb` (TCGA intake, QC, Hippo scoring, DEG)
   - `notebooks/02_Overlap between IBD and CRC based on Hippo pathway.ipynb` (cross-cohort overlap & validation)
   - `notebooks/03_Mechansitic Bridge b_w CRC and Dysbiosis.ipynb` (mechanistic integration)

5. **Monitor outputs** in `data_processed/` and `results/`:
   - Check figures for data quality
   - Verify DEG statistics and significance
   - Review pathway scores for biological plausibility

**Notes:**
- Notebooks are designed to run in **Google Colab**; if running locally, ensure Python 3.8+, required packages, and correct folder paths.
- Estimated runtime: ~30–60 min per notebook on Colab with CPU.
- Save outputs to Drive frequently; use checkpoints between notebooks.

## Acceptance Criteria (Sprint)

✓ **YAP (and DHR) scores:**
  - ↑ in TCGA Hippo-High vs Low tumors (p < 0.05)
  - ↑ in IBD (UC/CD) vs Control (p < 0.05)

✓ **All notebooks:**
  - Run end-to-end without errors
  - Generate expected outputs (DEG tables, figures, scores)
  - Colab "Run all" succeeds

✓ **Cross-cohort validation:**
  - Significant overlap between IBD and CRC Hippo-related DEGs
  - Reproducible mechanistic hypotheses linking dysbiosis and Hippo activation

## Key Analysis Outputs

| Analysis | Input Files | Output Files | Key Metric |
|----------|-------------|--------------|-----------|
| IBD QC & DEG | GSE235236 TPM, metadata | ibd_log2_tpm.csv, DEG tables, volcano plots | FDR < 0.05, \|log₂FC\| > 1 |
| TCGA QC & Hippo | TCGA counts, clinical | tcga_log2_cpm.csv, hippo_scores.csv, DEG table | Hippo-High vs Low p < 0.05 |
| Cross-cohort Overlap | IBD DEGs, TCGA DEGs, gene sets | overlap_genes.csv, Venn diagrams, shared heatmaps | Jaccard index, Fisher p-value |
| Dysbiosis–Hippo Bridge | Dysbiosis signatures, Hippo scores | dysbiosis_scores.csv, mechanistic_genes.csv, networks | Correlation r, bridging gene count |

## Team & Contacts

- **PI:** Sana Noor (sana11100noor@gmail.com | ceo@geneticcodon.com)
- **Data Integration (Team A):** [contact]
- **Analysis & Reporting (Team B):** [contact]

---

**Last updated:** 2026-07-06  
**Status:** In progress (notebooks in development & testing phase)
