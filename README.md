# Colorectal–Hippo–Dysbiosis (Genetic Codon)

**Goal:** Integrate bulk RNA-seq (TCGA-COAD) with tumor microbiome (TCMA), mutations (MAF), methylation (450k), and an IBD dysbiosis cohort (GSE235236) to examine how microbial imbalance activates Hippo/YAP–TAZ signaling in human colon.

🧠 **Organization:** Genetic Codon  
👩‍🔬 **PI:** Dr. Sana Noor  
📅 **Sprint:** 2 weeks  
👥 **Team:** A (Data Integration) | B (Analysis & Reporting)

## Cohorts
- **TCGA-COAD:** counts + clinical + MAF + 450k + microbiome (via TCMA).
- **GSE235236 (IBD):** colon mucosa (UC/CD vs controls).

## Specific Aims
1. Quantify Hippo/YAP activity in CRC tumor vs normal and IBD vs control; define a Dysbiosis–Hippo Response (DHR) module.
2. Associate YAP/DHR with Fusobacteriales/*F. nucleatum* abundance (TCMA) in TCGA-COAD.
3. Test mutation and methylation relationships with Hippo core genes; evaluate survival.

## Repo Layout
```
data_raw/ ... raw inputs (not versioned)
data_processed/
results/{qc,de,hippo,microbe,mutation,methylation,survival,figures,tables}
notebooks/ (00..09)
config/ {gene_sets.json, manifest_tcga.csv, manifest_ibd.csv}
docs/
```

## Getting Started (Colab)
1. Open `notebooks/00_env_and_utils.ipynb` → run boot cell.
2. Fill `config/manifest_tcga.csv` and `config/manifest_ibd.csv` with paths in your Drive.
3. Run `01_tcga_intake_qc.ipynb` and `02_ibd_intake_qc.ipynb`.

## Acceptance Criteria (Sprint)
- YAP (and DHR) scores ↑ in TCGA tumor vs normal, and in IBD vs control (p < 0.01).
- Positive association YAP/DHR ~ Fusobacteriales/*F. nucleatum* (p < 0.05) adjusted for MSI/purity.
- ≥1 Hippo core gene shows promoter methylation–expression anti-correlation (|r|>0.3, p<0.05).
- Survival (TCGA): YAP or DHR score significant in multivariable Cox (p < 0.05).
- All notebooks run end-to-end (Colab “Run all” succeeds).