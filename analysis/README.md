# analysis

Statistics and figures, run on the tables from `pipeline/03_collect_qc.py`. Every script writes a README
(files and columns) and `sessionInfo.txt` next to its outputs.

| Script | Purpose |
|---|---|
| `01_mirna_yield.R` | Small-RNA library yield by group in all samples (Wilcoxon, Fisher); supplementary yield figure |
| `02_mirna_de.R` | Differential miRNA expression (DESeq2 + edgeR QL; viral miRNAs as separate family); Fig 1A volcano. Sets: `primary` (ERCC criterion 1) and `sensitivity` (>=10k human miRNA counts) |
| `03_compare_nog2024.py` | Looks up the 29 miRNAs reported by Nog et al. 2024 in our DE results |

Scripts for Fig 1B (microbiome) and Fig 2 (Sankey) will be added.

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

R 4.5.1 (macOS arm64). Packages: DESeq2 1.50.2, edgeR 4.8.2, limma 3.66.0, apeglm 1.32.0, ANCOMBC 2.12.1,
vegan 2.7.6, ggplot2, ggrepel 0.9.8, ggalluvial 0.12.6, patchwork 1.3.2, scales 1.4.0, readr, dplyr, tidyr.
