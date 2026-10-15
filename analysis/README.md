# analysis

Statistics and figures, run on the tables from `pipeline/03_collect_qc.py`.

| Script | Purpose |
|---|---|
| `01_mirna_yield.R` | Small-RNA library yield by group in all samples (Wilcoxon, Fisher); supplementary yield figure |

Scripts for Fig 1A (volcano), Fig 1B (microbiome) and Fig 2 (Sankey) will be added.

## Environment

R 4.5.1 (macOS arm64). Packages: DESeq2 1.50.2, edgeR 4.8.2, limma 3.66.0, apeglm 1.32.0, ANCOMBC 2.12.1,
vegan 2.7.6, ggplot2, ggrepel 0.9.8, ggalluvial 0.12.6, patchwork 1.3.2, scales 1.4.0, readr, dplyr.
Each script also records `sessionInfo()` next to its outputs.
