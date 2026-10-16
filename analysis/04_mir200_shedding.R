#!/usr/bin/env Rscript
# Pre-specified exploratory check (PROJECT_LOG 2026-10-15): does the miR-200 family signal track
# urinary cell content? Specified before the metagenome data were complete.
#   1. miR-200 family score vs % human DNA in the metagenome (independent extraction + sequencing run)
#   2. miR-200 family score vs urinary BKPyV DNA (BK pairs per million QC pairs)
#   + descriptive cell-origin profile: epithelial/kidney-enriched vs blood/immune-cell miRNAs by group
# Spearman correlations, miRNA analysis set (ERCC criterion 1). No further tests.
#
# Usage: Rscript 04_mir200_shedding.R <mirna_counts_raw.tsv> <qc_combined.tsv> <out_dir>
suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(DESeq2); library(edgeR); library(patchwork) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: 04_mir200_shedding.R <counts.tsv> <qc_combined.tsv> <out_dir>")
out <- args[3]; dir.create(out, showWarnings = FALSE, recursive = TRUE)

MIR200 <- c("hsa-miR-141-3p", "hsa-miR-200a-3p", "hsa-miR-200b-3p", "hsa-miR-200c-3p", "hsa-miR-429")
ORIGIN <- list(
  "Epithelial / kidney-enriched" = c(MIR200, "hsa-miR-192-5p", "hsa-miR-194-5p"),
  "Blood / immune cells"         = c("hsa-miR-451a", "hsa-miR-486-5p", "hsa-miR-16-5p", "hsa-miR-223-3p",
                                     "hsa-miR-142-3p", "hsa-miR-150-5p"))

qc <- read_tsv(args[2], show_col_types = FALSE) |> filter(ercc_c1_pass) |>
  mutate(group = factor(group, levels = c("control", "DNAemia")),
         pct_human_dna = 100 * (mg_pairs_after_fastp - mg_pairs_nonhuman) / mg_pairs_after_fastp,
         bk_per_million = mg_bk_per_million_qc)
cnt <- read_tsv(args[1], show_col_types = FALSE)
m <- as.matrix(cnt[, qc$sample_id]); rownames(m) <- cnt$miRNA; storage.mode(m) <- "integer"
hsa <- m[grepl("^hsa-", rownames(m)), ]
hsa <- hsa[filterByExpr(DGEList(hsa, group = qc$group), group = qc$group), ]   # same filter as 02_mirna_de.R
cd <- data.frame(group = qc$group, row.names = qc$sample_id)
dds <- estimateSizeFactors(DESeqDataSetFromMatrix(hsa, cd, ~ group), type = "poscounts")
v <- assay(varianceStabilizingTransformation(dds, blind = TRUE))

fam <- intersect(MIR200, rownames(v))
qc$mir200_score <- colMeans(v[fam, qc$sample_id, drop = FALSE])

cor_one <- function(y, label) {
  ct <- suppressWarnings(cor.test(qc$mir200_score, y, method = "spearman", exact = TRUE))
  tibble(comparison = label, n = sum(!is.na(y)), rho = unname(ct$estimate), p = ct$p.value)
}
res <- bind_rows(cor_one(qc$pct_human_dna, "miR-200 score vs % human DNA (metagenome)"),
                 cor_one(qc$bk_per_million, "miR-200 score vs urinary BKPyV DNA (pairs per million)"))

prof <- bind_rows(lapply(names(ORIGIN), function(o) {
  g <- intersect(ORIGIN[[o]], rownames(v))
  tibble(sample_id = qc$sample_id, group = qc$group, origin = o, n_mirnas = length(g),
         mirnas = paste(sub("hsa-", "", g), collapse = ", "), score = colMeans(v[g, qc$sample_id, drop = FALSE]))
}))
prof_test <- prof |> group_by(origin, n_mirnas, mirnas) |>
  summarise(median_control = median(score[group == "control"]), median_DNAemia = median(score[group == "DNAemia"]),
            diff_median = median_DNAemia - median_control,
            p_wilcoxon = wilcox.test(score ~ group, exact = TRUE)$p.value, .groups = "drop")

write_tsv(qc |> select(sample_id, group, mir200_score, pct_human_dna, bk_per_million, mi_counts_hsa),
          file.path(out, "mir200_per_sample.tsv"))
write_tsv(res, file.path(out, "mir200_correlations.tsv"))
write_tsv(prof, file.path(out, "cell_origin_per_sample.tsv"))
write_tsv(prof_test, file.path(out, "cell_origin_by_group.tsv"))

cols <- c(control = "#4C78A8", DNAemia = "#E45756")
sc <- function(x, xl, r) ggplot(qc, aes(.data[[x]], mir200_score, colour = group)) +
  geom_point(size = 2.4) + geom_text(aes(label = sample_id), size = 2.6, vjust = -0.9, show.legend = FALSE) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = xl, y = "miR-200 family score (mean VST)",
       subtitle = sprintf("Spearman rho = %.2f, p = %.2g, n = %d", r$rho, r$p, r$n)) +
  theme_bw(base_size = 10)
p1 <- sc("pct_human_dna", "% human DNA in metagenome", res[1, ])
p2 <- sc("bk_per_million", "Urinary BKPyV DNA (pairs per million, log scale)", res[2, ]) + scale_x_log10()
p3 <- ggplot(prof, aes(group, score, colour = group)) +
  geom_boxplot(outlier.shape = NA, colour = "grey40", width = 0.5) +
  geom_point(position = position_jitter(width = 0.1, seed = 1), size = 2) +
  facet_wrap(~origin) + scale_colour_manual(values = cols, guide = "none") +
  labs(x = NULL, y = "Mean VST of marker miRNAs") + theme_bw(base_size = 10)
ggsave(file.path(out, "FigS_mir200_shedding.pdf"), (p1 | p2) / p3 + plot_layout(guides = "collect"), width = 9, height = 7)
ggsave(file.path(out, "FigS_mir200_shedding.png"), (p1 | p2) / p3 + plot_layout(guides = "collect"), width = 9, height = 7, dpi = 300)

writeLines(c(
  "# miR-200 family and urinary cell content (generated by analysis/04_mir200_shedding.R; do not edit by hand)", "",
  "Pre-specified exploratory check (PROJECT_LOG 2026-10-15). Samples: miRNA analysis set (ERCC criterion 1).", "",
  "| File | Contents |", "|---|---|",
  "| `mir200_per_sample.tsv` | miR-200 family score, % human DNA, BKPyV DNA per sample. |",
  "| `mir200_correlations.tsv` | The two pre-specified Spearman correlations (`rho`, exact `p`). |",
  "| `cell_origin_per_sample.tsv` | Mean VST of epithelial/kidney-enriched vs blood/immune-cell marker miRNAs per sample. |",
  "| `cell_origin_by_group.tsv` | Group medians, difference (DNAemia - control) and exact Wilcoxon p, descriptive. |",
  "| `FigS_mir200_shedding.pdf/.png` | Scatter plots of both correlations and the cell-origin profile. |", "",
  paste0("miR-200 family score = mean variance-stabilised expression (DESeq2 VST, blind) of the family members passing ",
         "the expression filter: ", paste(sub("hsa-", "", fam), collapse = ", "), "."),
  "% human DNA = read pairs aligning concordantly to GRCh38 / read pairs after fastp (metagenome).",
  "Note: % human DNA is lowered by bacteriuria and by BKPyV viruria, so it is an imperfect cell-content proxy."
), file.path(out, "README.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(as.data.frame(res)); print(as.data.frame(prof_test |> select(-mirnas)))
print(as.data.frame(qc |> select(sample_id, group, mir200_score, pct_human_dna, bk_per_million) |> arrange(group)))
