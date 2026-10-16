# analysis

Statistics and figures, run on the tables from `pipeline/03_collect_qc.py`. Every script writes a README
(files and columns) and `sessionInfo.txt` next to its outputs.

| Script | Purpose |
|---|---|
| `01_mirna_yield.R` | Small-RNA library yield by group in all samples (Wilcoxon, Fisher); Figure 1A panel and supplementary yield figure |
| `02_mirna_de.R` | Differential miRNA expression (DESeq2 + edgeR QL; viral miRNAs as separate family); Fig 1A volcano. Sets: `primary` (ERCC criterion 1) and `sensitivity` (>=10k human miRNA counts) |
| `03_compare_nog2024.py` | Looks up the 29 miRNAs reported by Nog et al. 2024 in our DE results |
| `04_mir200_shedding.R` | Pre-specified check: miR-200 family vs % human DNA and urinary BKPyV DNA; cell-origin miRNA profile |
| `05_microbiome.R` | Bacterial microbiome (MetaPhlAn 4): mean-RA > 0.1% filter, optional cross-talk filter, alpha (raw + rarefied), Aitchison PERMANOVA, ALDEx2; urinary BK/JC DNA; Fig 1B. Runs: all 22 (supplement), `0.01` cross-talk sensitivity, `0 S18` without the UTI control, and the 11 paired samples (input for Figure 2) |
| `06_sankey.R` | Figure 2: descriptive Sankey genus -> sample -> miRNA (paired samples; no tests). Args: genus table, DE dir, QC table, out dir, MetaPhlAn est-reads table (for phylum colours), optional legend image (`assets/sankey_legend.png`, author-designed) |
| `07_main_figures.R` | Assembles Figure 1 (2x2: A yield, all 22; B volcano, n=11; C alpha, D PCA, all 22) and Figure 2 (A heatmap; B Sankey) from saved panels |
| `08_heatmap.R` | Figure 2A: Spearman correlations, nominal DE miRNAs + viral vs genera (prevalence >= 3/11) and urinary BKPyV DNA; BH-FDR |



## Running

```bash
R=~/Documents/BKPyV_miRNA_microbiome/results
Rscript 01_mirna_yield.R $R/qc/qc_combined_final_mirna.tsv <sample_manifest.tsv> $R/mirna/yield
for s in primary sensitivity; do
  Rscript 02_mirna_de.R $R/qc/mirna_counts_raw.tsv $R/qc/qc_combined_final_mirna.tsv $s $R/mirna/de_$s
  python3 03_compare_nog2024.py $R/mirna/de_$s/de_human_mirna.tsv $R/mirna/de_$s/de_viral_mirna.tsv \
    $R/mirna/comparison_nog2024/nog2024_vs_$s.tsv
done
```

## Environment

R 4.5.1 (macOS arm64). Packages: DESeq2 1.50.2, edgeR 4.8.2, limma 3.66.0, apeglm 1.32.0, vegan 2.7.6, ALDEx2 1.42.0, ggplot2, ggrepel 0.9.8, ggalluvial 0.12.6, patchwork 1.3.2, scales 1.4.0, readr, dplyr, tidyr.
