# Colorectal–Hippo–Dysbiosis (Genetic Codon)

**Goal:** Integrate bulk RNA-seq (TCGA-COAD) with tumor microbiome (TCMA) and an IBD dysbiosis cohort (GSE235236) to examine how microbial imbalance activates Hippo/YAP signaling and its relationships with microbial signatures and transcriptional programs.

🧠 **Organization:** Genetic Codon  
👩‍🔬 **PI:** Sana Noor  
📅 **Sprint:** 2 months  
👥 **Team:** A (Data Integration) | B (Analysis & Reporting)

## Cohorts
- **TCGA-CRC:** counts + clinical + microbiome (via TCMA).
- **GSE235236 (IBD):** colon mucosa (UC/CD vs controls).

## Specific Aims
1. Quantify Hippo/YAP activity in CRC tumor vs normal and IBD vs control; define a Dysbiosis–Hippo Response (DHR) module.  
2. Associate YAP/DHR with Fusobacteriales/*F. nucleatum* abundance (TCMA) in TCGA-COAD.  
3. Generate candidate mechanistic links between dysbiosis signatures and Hippo pathway activation for follow-up.

## Repo Layout
```
data_raw/ ... raw inputs (not versioned)
data_processed/
results/{qc,de,hippo,microbe,figures,tables}
notebooks/
  ├─ 00_env_and_utils.ipynb
  ├─ 01A_IBD_intake_qc_deg.ipynb
  ├─ 01B_TCGA_Analysis.ipynb
  ├─ 02_Overlap between IBD and CRC based on Hippo pathway.ipynb
  └─ 03_Mechansitic Bridge b_w CRC and Dysbiosis.ipynb
config/ {gene_sets.json, manifest_tcga.csv, manifest_ibd.csv}
docs/
```

## Notebooks (brief)
- `notebooks/00_env_and_utils.ipynb` — Environment bootstrapping and utility functions. Run first (installs packages, imports, helper functions, and common plotting/theme settings).
- `notebooks/01A_IBD_intake_qc_deg.ipynb` — IBD cohort intake, QC, normalisation, and differential expression (DEG) analyses comparing UC/CD vs controls; generates IBD-specific Hippo/DHR scores and QC outputs.
- `notebooks/01B_TCGA_Analysis.ipynb` — TCGA-CRC intake, QC, Hippo/YAP scoring. Primary TCGA analysis notebook focused on transcriptional signatures for Hippo high vs low tumours.
- `notebooks/02_Overlap between IBD and CRC based on Hippo pathway.ipynb` — Cross-cohort comparison: overlap between IBD and CRC differential signals focused on Hippo pathway and DHR module.
- `notebooks/03_Mechansitic Bridge b_w CRC and Dysbiosis.ipynb` — Exploratory mechanistic analyses linking dysbiosis (microbe signatures) and Hippo pathway activation in CRC; includes integrative visualisations and candidate hypotheses.

## Getting Started (Colab)
1. Open `notebooks/00_env_and_utils.ipynb` and run the boot cell to install dependencies and mount your Drive (if using Google Drive).  
2. Fill `config/manifest_tcga.csv` and `config/manifest_ibd.csv` with paths to your datasets (Drive paths or local paths, depending on runtime).  
3. Recommended run order:
   - `notebooks/00_env_and_utils.ipynb` (environment & utils)
   - `notebooks/01A_IBD_intake_qc_deg.ipynb` (IBD intake, QC, DEG, IBD Hippo/DHR)
   - `notebooks/01B_TCGA_Analysis.ipynb` (TCGA intake, QC, Hippo scoring)
   - `notebooks/02_Overlap between IBD and CRC based on Hippo pathway.ipynb` (cross-cohort overlap)
   - `notebooks/03_Mechansitic Bridge b_w CRC and Dysbiosis.ipynb` (mechanistic integration)
4. Save results in `data_processed/` and `results/` as configured in the notebooks.

Notes:
- Notebooks are designed to run in Colab; if running locally, confirm the environment and paths.

## Acceptance Criteria (Sprint)
- YAP (and DHR) scores ↑ in TCGA Hippo High vs Low tumours, and in IBD vs control (p < 0.05).  
- Positive association YAP/DHR ~ Fusobacteriales/*F. nucleatum* (p < 0.05) adjusted for MSI/purity.  
- All notebooks run end-to-end (Colab “Run all” succeeds).
