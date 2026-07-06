"""
Generate all 7 paper figures for the Colorectal-Hippo-Dysbiosis manuscript.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
from scipy import stats
import os, warnings
warnings.filterwarnings('ignore')

BASE = "/home/sana/Colorectal_Hippo_Dysbiosis/results_uploaded/Colon_cancer_project"
FIGS = "/home/sana/Colorectal_Hippo_Dysbiosis/docs/figures"
os.makedirs(FIGS, exist_ok=True)

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 11,
    'axes.titlesize': 13, 'axes.labelsize': 12,
    'xtick.labelsize': 10, 'ytick.labelsize': 10,
    'figure.dpi': 150, 'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

YAP_TARGETS = ['CTGF', 'CYR61', 'ANKRD1', 'AREG', 'AMOTL2', 'BIRC5', 'FSTL1']

# ── Load data ──────────────────────────────────────────────────────────────────
meta   = pd.read_csv(f"{BASE}/data_processed/ibd/ibd_metadata_clean.csv")
ibd_l2 = pd.read_csv(f"{BASE}/data_processed/ibd/ibd_log2_tpm.csv", index_col=0)
ibd_l2.index.name = 'gene'

ibd_deg = pd.read_csv(f"{BASE}/results/de/ibd/IBD_DEG_IBD_vs_Control.csv")
ibd_deg['logFC'] = ibd_deg['logFC_IBD_minus_Control']
ibd_deg['padj']  = pd.to_numeric(ibd_deg['padj_fdr'], errors='coerce')
ibd_deg['-log10FDR'] = -np.log10(ibd_deg['padj'].clip(lower=1e-300))

tcga_deg = pd.read_csv(f"{BASE}/results/de/tcga/TCGA_DEG_HippoHigh_vs_HippoLow.csv")
tcga_deg['logFC'] = pd.to_numeric(tcga_deg['logFC_Hippo_high_minus_Hippo_low'], errors='coerce')
tcga_deg['padj']  = pd.to_numeric(tcga_deg['padj_fdr'], errors='coerce')
tcga_deg['-log10FDR'] = -np.log10(tcga_deg['padj'].clip(lower=1e-300))

scores  = pd.read_csv(f"{BASE}/results/hippo/tcga_hippo_yap_scores.csv")
core_df = pd.read_csv(f"{BASE}/results/novelty/core_signature_genes.csv")
core_genes = core_df['core_signature_genes'].tolist()
ovlp = pd.read_csv(f"{BASE}/results/novelty/overlap_IBD_TCGA_same_direction.csv")
ovlp_all = pd.read_csv(f"{BASE}/results/novelty/overlap_IBD_TCGA_all.csv")
enrich  = pd.read_csv(f"{BASE}/results/novelty/core_signature_enrichment.csv")
enrich['Adjusted P-value'] = pd.to_numeric(enrich['Adjusted P-value'], errors='coerce')

# Compute IBD YAP scores
avail_yap = [g for g in YAP_TARGETS if g in ibd_l2.index]
subset_yap = ibd_l2.loc[avail_yap]
z_yap = subset_yap.apply(lambda row: (row - row.mean()) / (row.std() + 1e-8), axis=1)
yap_score = z_yap.mean(axis=0)
yap_df = pd.DataFrame({'sample': yap_score.index, 'yap_score': yap_score.values})
yap_df = yap_df.merge(meta[['index','group']], left_on='sample', right_on='index')

# Compute IBD core signature scores
avail_core = [g for g in core_genes if g in ibd_l2.index]
subset_core = ibd_l2.loc[avail_core]
z_core = subset_core.apply(lambda row: (row - row.mean()) / (row.std() + 1e-8), axis=1)
core_score = z_core.mean(axis=0)
core_score_df = pd.DataFrame({'sample': core_score.index, 'core_score': core_score.values})
combined_df = yap_df.merge(core_score_df, on='sample')

PALETTE = {'UC': '#e05c5c', 'CD': '#f5a623', 'Control': '#4a90d9'}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — IBD volcano + heatmap
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# 1A: Volcano IBD vs HC
ax = axes[0]
sig_up = (ibd_deg['padj'] < 0.05) & (ibd_deg['logFC'] > 0)
sig_dn = (ibd_deg['padj'] < 0.05) & (ibd_deg['logFC'] < 0)
ns = ~(sig_up | sig_dn)

ax.scatter(ibd_deg.loc[ns, 'logFC'], ibd_deg.loc[ns, '-log10FDR'],
           s=4, c='lightgrey', alpha=0.4, rasterized=True)
ax.scatter(ibd_deg.loc[sig_dn, 'logFC'], ibd_deg.loc[sig_dn, '-log10FDR'],
           s=6, c='#4a90d9', alpha=0.6, label=f'Down ({sig_dn.sum()})')
ax.scatter(ibd_deg.loc[sig_up, 'logFC'], ibd_deg.loc[sig_up, '-log10FDR'],
           s=6, c='#e05c5c', alpha=0.6, label=f'Up ({sig_up.sum()})')

# label top genes
top_up = ibd_deg[sig_up].nlargest(5, 'logFC')
top_dn = ibd_deg[sig_dn].nsmallest(5, 'logFC')
for _, r in pd.concat([top_up, top_dn]).iterrows():
    ax.annotate(r['gene'], (r['logFC'], r['-log10FDR']),
                fontsize=7.5, ha='center', va='bottom',
                xytext=(0, 3), textcoords='offset points')

ax.axhline(-np.log10(0.05), color='k', lw=0.8, ls='--', alpha=0.5)
ax.axvline(0, color='k', lw=0.5, alpha=0.3)
ax.set_xlabel('log₂ Fold Change (IBD vs HC)')
ax.set_ylabel('–log₁₀(FDR)')
ax.set_title('A  IBD vs Healthy Controls', fontweight='bold', loc='left')
ax.legend(frameon=False, fontsize=9)

# 1B: Heatmap top 30 DEGs
ax = axes[1]
top_genes = pd.concat([
    ibd_deg[sig_up].nlargest(15, '-log10FDR'),
    ibd_deg[sig_dn].nlargest(15, '-log10FDR')
])['gene'].tolist()
avail_top = [g for g in top_genes if g in ibd_l2.index]

hm_data = ibd_l2.loc[avail_top]
# Z-score across samples
hm_z = hm_data.apply(lambda row: (row - row.mean()) / (row.std() + 1e-8), axis=1)

# Sort samples by group
sample_order = meta.sort_values('group')['index'].tolist()
sample_order = [s for s in sample_order if s in hm_z.columns]
hm_z = hm_z[sample_order]

group_colors = meta.set_index('index').loc[sample_order, 'group'].map(PALETTE)

im = ax.imshow(hm_z.values, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
ax.set_xticks([])
ax.set_yticks(range(len(avail_top)))
ax.set_yticklabels(avail_top, fontsize=7.5)
ax.set_xlabel('Samples (HC → CD → UC)')
ax.set_title('B  Top DEG Expression Heatmap (Z-score)', fontweight='bold', loc='left')

# Color bar on top for groups
n = len(sample_order)
for i, s in enumerate(sample_order):
    g = meta.set_index('index').loc[s, 'group']
    ax.add_patch(plt.Rectangle((i - 0.5, -1.8), 1, 1.2,
                                color=PALETTE[g], clip_on=False, transform=ax.transData))

cb = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
cb.set_label('Z-score', fontsize=9)

patches = [mpatches.Patch(color=v, label=k) for k, v in PALETTE.items()]
ax.legend(handles=patches, loc='lower right', fontsize=8, frameon=False,
          bbox_to_anchor=(1.02, 0))

plt.tight_layout()
plt.savefig(f"{FIGS}/Fig1_IBD_Volcano_Heatmap.png")
plt.close()
print("Fig 1 saved")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — YAP score boxplot IBD vs HC
# ══════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(6, 5))

order = ['Control', 'UC', 'CD']
data_plot = [yap_df[yap_df['group'] == g]['yap_score'].values for g in order]
colors = [PALETTE[g] for g in order]
labels = ['HC\n(n=8)', 'UC\n(n=26)', 'CD\n(n=22)']

bp = ax.boxplot(data_plot, patch_artist=True, widths=0.5,
                medianprops=dict(color='black', lw=2),
                whiskerprops=dict(lw=1.2),
                capprops=dict(lw=1.2),
                flierprops=dict(marker='o', markersize=4, alpha=0.5))
for patch, c in zip(bp['boxes'], colors):
    patch.set_facecolor(c)
    patch.set_alpha(0.7)

# Add individual points
for i, (g, c) in enumerate(zip(order, colors)):
    y = yap_df[yap_df['group'] == g]['yap_score'].values
    x = np.random.normal(i + 1, 0.06, len(y))
    ax.scatter(x, y, s=20, c=c, alpha=0.8, zorder=5, edgecolors='white', lw=0.5)

ibd_sc = yap_df[yap_df['group'].isin(['UC', 'CD'])]['yap_score']
hc_sc  = yap_df[yap_df['group'] == 'Control']['yap_score']
_, mw_p = stats.mannwhitneyu(ibd_sc, hc_sc, alternative='two-sided')
p_str = f"p = {mw_p:.3f}" if mw_p >= 0.001 else f"p = {mw_p:.2e}"

ax.set_xticklabels(labels)
ax.set_ylabel('YAP Target Score (mean Z-score)')
ax.set_title('YAP/Hippo Activity Score Across IBD Cohort', fontweight='bold')
ax.text(0.98, 0.97, f'IBD vs HC: {p_str}\n(Mann-Whitney U)',
        transform=ax.transAxes, ha='right', va='top', fontsize=9,
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='grey', alpha=0.8))
ax.axhline(0, color='k', lw=0.5, ls='--', alpha=0.3)

plt.tight_layout()
plt.savefig(f"{FIGS}/Fig2_YAP_Score_Boxplot.png")
plt.close()
print("Fig 2 saved")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — TCGA Volcano Hippo-High vs Hippo-Low
# ══════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(7, 6))

sig_up = (tcga_deg['padj'] < 0.05) & (tcga_deg['logFC'] > 0)
sig_dn = (tcga_deg['padj'] < 0.05) & (tcga_deg['logFC'] < 0)
ns = ~(sig_up | sig_dn)

ax.scatter(tcga_deg.loc[ns, 'logFC'], tcga_deg.loc[ns, '-log10FDR'],
           s=4, c='lightgrey', alpha=0.3, rasterized=True)
ax.scatter(tcga_deg.loc[sig_dn, 'logFC'], tcga_deg.loc[sig_dn, '-log10FDR'],
           s=6, c='#4a90d9', alpha=0.6, label=f'Down ({sig_dn.sum()})', rasterized=True)
ax.scatter(tcga_deg.loc[sig_up, 'logFC'], tcga_deg.loc[sig_up, '-log10FDR'],
           s=6, c='#e05c5c', alpha=0.6, label=f'Up ({sig_up.sum()})', rasterized=True)

top_up = tcga_deg[sig_up].nlargest(5, 'logFC')
top_dn = tcga_deg[sig_dn].nsmallest(5, 'logFC')
for _, r in pd.concat([top_up, top_dn]).iterrows():
    ax.annotate(r['gene'], (r['logFC'], r['-log10FDR']),
                fontsize=7.5, ha='center', va='bottom',
                xytext=(0, 3), textcoords='offset points')

ax.axhline(-np.log10(0.05), color='k', lw=0.8, ls='--', alpha=0.5)
ax.axvline(0, color='k', lw=0.5, alpha=0.3)
ax.set_xlabel('log₂ Fold Change (Hippo-High vs Hippo-Low)')
ax.set_ylabel('–log₁₀(FDR)')
ax.set_title('TCGA-COAD: Hippo-High vs Hippo-Low DEGs', fontweight='bold')
ax.legend(frameon=False, fontsize=9)

plt.tight_layout()
plt.savefig(f"{FIGS}/Fig3_TCGA_Volcano.png")
plt.close()
print("Fig 3 saved")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Cross-cohort logFC scatter
# ══════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(6, 6))

# Use the overlap_all file
ovlp_all2 = ovlp_all.copy()
ovlp_all2['logFC_IBD'] = pd.to_numeric(ovlp_all2['logFC_IBD_minus_Control'], errors='coerce')
ovlp_all2['logFC_TCGA'] = pd.to_numeric(ovlp_all2['logFC_Hippo_high_minus_Hippo_low'], errors='coerce')
ovlp_all2 = ovlp_all2.dropna(subset=['logFC_IBD','logFC_TCGA'])

same_dir = ovlp_all2['direction_match'] == 'same_direction'
ax.scatter(ovlp_all2.loc[~same_dir, 'logFC_IBD'], ovlp_all2.loc[~same_dir, 'logFC_TCGA'],
           s=30, c='#888888', alpha=0.6, label=f'Discordant (n={(~same_dir).sum()})')
ax.scatter(ovlp_all2.loc[same_dir, 'logFC_IBD'], ovlp_all2.loc[same_dir, 'logFC_TCGA'],
           s=40, c='#e05c5c', alpha=0.8, label=f'Concordant (n={same_dir.sum()})')

# Diagonal
lim = max(abs(ovlp_all2['logFC_IBD']).max(), abs(ovlp_all2['logFC_TCGA']).max()) * 1.1
ax.plot([-lim, lim], [-lim, lim], 'k--', lw=1, alpha=0.5)
ax.axhline(0, color='k', lw=0.5, alpha=0.3)
ax.axvline(0, color='k', lw=0.5, alpha=0.3)

# Highlight top concordant genes
top_conc = ovlp_all2[same_dir].nlargest(5, 'logFC_IBD')
for _, r in top_conc.iterrows():
    ax.annotate(r['gene'], (r['logFC_IBD'], r['logFC_TCGA']),
                fontsize=7.5, xytext=(3, 2), textcoords='offset points')

ax.set_xlabel('log₂FC IBD vs HC')
ax.set_ylabel('log₂FC TCGA Hippo-High vs Hippo-Low')
ax.set_title('Cross-Cohort DEG Concordance\n(IBD vs TCGA Hippo-High)', fontweight='bold')
ax.legend(frameon=False, fontsize=9)

rho, rho_p = stats.spearmanr(ovlp_all2['logFC_IBD'], ovlp_all2['logFC_TCGA'])
ax.text(0.05, 0.95, f'ρ = {rho:.2f}, p = {rho_p:.3e}',
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='grey', alpha=0.8))

plt.tight_layout()
plt.savefig(f"{FIGS}/Fig4_CrossCohort_Scatter.png")
plt.close()
print("Fig 4 saved")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — DHR heatmaps (A: IBD, B: TCGA)
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 7))

# Subset to top 30 core genes by IBD logFC
ibd_core_fc = ibd_deg[ibd_deg['gene'].isin(avail_core)].sort_values('logFC', ascending=False)
top30 = ibd_core_fc['gene'].head(30).tolist()
avail30 = [g for g in top30 if g in ibd_l2.index]

# 5A: IBD
ax = axes[0]
hm = ibd_l2.loc[avail30]
hm_z = hm.apply(lambda row: (row - row.mean()) / (row.std() + 1e-8), axis=1)
sample_order = meta.sort_values('group')['index'].tolist()
sample_order = [s for s in sample_order if s in hm_z.columns]
hm_z = hm_z[sample_order]

im = ax.imshow(hm_z.values, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
ax.set_xticks([])
ax.set_yticks(range(len(avail30)))
ax.set_yticklabels(avail30, fontsize=7.5)
ax.set_title('A  DHR Core Signature — IBD Cohort', fontweight='bold', loc='left')
ax.set_xlabel('Samples (HC → CD → UC)')
plt.colorbar(im, ax=ax, shrink=0.5, label='Z-score')

for i, s in enumerate(sample_order):
    g = meta.set_index('index').loc[s, 'group']
    ax.add_patch(plt.Rectangle((i - 0.5, -1.8), 1, 1.2,
                                color=PALETTE[g], clip_on=False, transform=ax.transData))

# 5B: TCGA
ax = axes[1]
# Load TCGA expression
tcga_l2 = pd.read_csv(f"{BASE}/data_processed/tcga/tcga_log2_rsem.csv", index_col=0)
print(f"TCGA expr shape: {tcga_l2.shape}")
tcga_l2.index = tcga_l2.index.astype(str)

# Stratify
q25 = scores['yap_targets'].quantile(0.25)
q75 = scores['yap_targets'].quantile(0.75)
high_s = scores[scores['yap_targets'] >= q75]['Unnamed: 0'].tolist()
low_s  = scores[scores['yap_targets'] <= q25]['Unnamed: 0'].tolist()

avail_t = [g for g in avail30 if g in tcga_l2.index]
tcga_cols = [c for c in high_s + low_s if c in tcga_l2.columns]
if avail_t and tcga_cols:
    tcga_hm = tcga_l2.loc[avail_t, tcga_cols]
    tcga_hm_z = tcga_hm.apply(lambda row: (row - row.mean()) / (row.std() + 1e-8), axis=1)
    im2 = ax.imshow(tcga_hm_z.values, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
    ax.set_xticks([])
    ax.set_yticks(range(len(avail_t)))
    ax.set_yticklabels(avail_t, fontsize=7.5)
    ax.set_title('B  DHR Core Signature — TCGA Cohort', fontweight='bold', loc='left')
    ax.set_xlabel('Samples (Hippo-Low → Hippo-High)')
    plt.colorbar(im2, ax=ax, shrink=0.5, label='Z-score')
    n_h = sum(1 for c in tcga_cols if c in high_s)
    n_l = sum(1 for c in tcga_cols if c in low_s)
    for i, c in enumerate(tcga_cols):
        color = '#e05c5c' if c in high_s else '#4a90d9'
        ax.add_patch(plt.Rectangle((i - 0.5, -1.8), 1, 1.2,
                                    color=color, clip_on=False, transform=ax.transData))
    patches2 = [mpatches.Patch(color='#e05c5c', label=f'Hippo-High (n={n_h})'),
                mpatches.Patch(color='#4a90d9', label=f'Hippo-Low (n={n_l})')]
    ax.legend(handles=patches2, loc='lower right', fontsize=8, frameon=False)

plt.tight_layout()
plt.savefig(f"{FIGS}/Fig5_DHR_Heatmaps.png")
plt.close()
print("Fig 5 saved")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 6 — Enrichment dotplots
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(16, 6))
libs = ['KEGG_2021_Human', 'GO_Biological_Process_2021', 'MSigDB_Hallmark_2020']
titles = ['A  KEGG 2021', 'B  GO Biological Process', 'C  MSigDB Hallmark']

for ax, lib, title in zip(axes, libs, titles):
    sub = enrich[enrich['Gene_set'] == lib].dropna(subset=['Adjusted P-value'])
    sub = sub.sort_values('Adjusted P-value').head(10)
    sub = sub.sort_values('Combined Score')
    
    y_pos = range(len(sub))
    sc = ax.scatter(sub['Combined Score'], y_pos,
                    s=80, c=-np.log10(sub['Adjusted P-value']+1e-10),
                    cmap='Reds', vmin=0, zorder=5)
    ax.barh(y_pos, sub['Combined Score'], height=0.4, color='lightgrey', alpha=0.5)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels([t[:45] for t in sub['Term']], fontsize=8)
    ax.set_xlabel('Combined Score')
    ax.set_title(title, fontweight='bold', loc='left')
    cb = plt.colorbar(sc, ax=ax, shrink=0.5, pad=0.02)
    cb.set_label('–log₁₀(FDR)', fontsize=8)

plt.tight_layout()
plt.savefig(f"{FIGS}/Fig6_Enrichment.png")
plt.close()
print("Fig 6 saved")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 7 — Archetype schematic + YAP vs DHR scatter
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# 7A: Archetype schematic
ax = axes[0]
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis('off')
ax.set_title('A  DHR Archetype Classification', fontweight='bold', loc='left')

# Central node
central = plt.Circle((5, 7.5), 1.2, color='#2c3e50', zorder=5)
ax.add_patch(central)
ax.text(5, 7.5, 'YAP/TAZ\nActivation', ha='center', va='center',
        color='white', fontsize=9, fontweight='bold', zorder=6)

# Three archtypes
arch_data = [
    (2.0, 3.5, '#e74c3c', 'Immune Evasion\n(DHR-IE)',   ['CD70', 'TNFRSF4', 'TNFRSF6B', 'MIR663A']),
    (5.0, 1.8, '#27ae60', 'ECM/Fibrosis\n(DHR-EF)',     ['SFRP2', 'COL10A1', 'THBS4', 'COMP']),
    (8.0, 3.5, '#8e44ad', 'Metabolic\nImpairment\n(DHR-MI)', ['ADH1C', 'HADHB', 'G0S2', 'SPHK1']),
]

for (x, y, color, label, genes) in arch_data:
    rect = mpatches.FancyBboxPatch((x-1.5, y-1.2), 3, 2.4,
                                    boxstyle='round,pad=0.15',
                                    facecolor=color, edgecolor='white',
                                    alpha=0.85, zorder=4)
    ax.add_patch(rect)
    ax.text(x, y + 0.5, label, ha='center', va='center',
            color='white', fontsize=8.5, fontweight='bold', zorder=5)
    ax.text(x, y - 0.5, '\n'.join(genes[:3]), ha='center', va='top',
            color='white', fontsize=7, zorder=5, style='italic')
    # Arrow from central
    ax.annotate('', xy=(x, y + 1.2), xytext=(5, 6.3),
                arrowprops=dict(arrowstyle='->', color='#555555', lw=2),
                zorder=3)

# Dysbiosis input
ax.text(5, 9.7, 'Gut Dysbiosis', ha='center', va='center',
        fontsize=10, fontweight='bold', color='#c0392b',
        bbox=dict(boxstyle='round,pad=0.4', fc='#fadbd8', ec='#c0392b'))
ax.annotate('', xy=(5, 8.7), xytext=(5, 9.3),
            arrowprops=dict(arrowstyle='->', color='#c0392b', lw=2))

# CRC output
ax.text(5, 0.4, '→ Colorectal Cancer Risk', ha='center', va='center',
        fontsize=9, color='#922b21',
        bbox=dict(boxstyle='round,pad=0.3', fc='#f9ebea', ec='#922b21'))

# 7B: YAP vs DHR scatter
ax = axes[1]
ibd_only = combined_df[combined_df['group'].isin(['UC','CD'])]
colors_sc = [PALETTE[g] for g in ibd_only['group']]
ax.scatter(ibd_only['yap_score'], ibd_only['core_score'],
           c=colors_sc, s=50, alpha=0.8, edgecolors='white', lw=0.5)

rho, rho_p = stats.spearmanr(combined_df['yap_score'], combined_df['core_score'])
m, b = np.polyfit(ibd_only['yap_score'], ibd_only['core_score'], 1)
x_line = np.linspace(ibd_only['yap_score'].min(), ibd_only['yap_score'].max(), 100)
ax.plot(x_line, m * x_line + b, 'k--', lw=1.5, alpha=0.7)

p_str = f"p = {rho_p:.2e}"
ax.text(0.05, 0.95, f'Spearman ρ = {rho:.3f}\n{p_str}',
        transform=ax.transAxes, fontsize=10,
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='grey', alpha=0.9),
        va='top')

patches = [mpatches.Patch(color=PALETTE['UC'], label='UC (n=26)'),
           mpatches.Patch(color=PALETTE['CD'], label='CD (n=22)')]
ax.legend(handles=patches, frameon=False, fontsize=9)
ax.set_xlabel('YAP Target Score')
ax.set_ylabel('DHR Core Signature Score')
ax.set_title('B  YAP Activity vs DHR Score (IBD Cohort)', fontweight='bold', loc='left')

plt.tight_layout()
plt.savefig(f"{FIGS}/Fig7_Archetype_Scatter.png")
plt.close()
print("Fig 7 saved")

print("\n✅ All figures saved to", FIGS)
