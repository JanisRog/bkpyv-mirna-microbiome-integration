#!/usr/bin/env Rscript
# Small-RNA library yield, DNAemia vs control, in ALL samples (no QC exclusion).
# Question: do control urines yield less miRNA? Raw sequencing depth is tested as a negative control
# (it should not differ if the difference is in RNA content, not in sequencing).
#
# Usage: Rscript 01_mirna_yield.R <qc_combined.tsv> <sample_manifest.tsv> <out_dir>
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: 01_mirna_yield.R <qc_combined.tsv> <sample_manifest.tsv> <out_dir>")
qc_f <- args[1]; man_f <- args[2]; out <- args[3]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

ERCC_MIN_ANNOTATED <- 100000   # ERCC criterion 1, the exclusion rule for the DE analysis

qc  <- read_tsv(qc_f, show_col_types = FALSE)
man <- read_tsv(man_f, show_col_types = FALSE) |> select(sample_id, group)
d <- qc |> select(-any_of("group")) |> inner_join(man, by = "sample_id") |>
  mutate(group = factor(group, levels = c("control", "DNAemia")))
stopifnot(nrow(d) == nrow(man))

metrics <- c(
  mi_counts_hsa        = "Human miRNA reads",
  mi_reads_17_30nt     = "Reads 17-30 nt after trimming",
  ercc_reads_annotated = "Reads on annotated human RNA (ERCC criterion 1)",
  mi_reads_raw         = "Raw reads (sequencing depth)"
)

res <- lapply(names(metrics), function(m) {
  x <- d[[m]]
  w <- wilcox.test(log10(x + 1) ~ d$group, exact = TRUE, conf.int = TRUE)
  med <- tapply(x, d$group, median)
  tibble(metric = m, label = metrics[[m]],
         median_control = med[["control"]], median_DNAemia = med[["DNAemia"]],
         fold_median = med[["DNAemia"]] / med[["control"]],
         # Hodges-Lehmann shift on log10 scale -> fold change (DNAemia vs control; sign from factor order)
         hl_fold = 10^(-unname(w$estimate)), hl_fold_lo = 10^(-w$conf.int[2]), hl_fold_hi = 10^(-w$conf.int[1]),
         W = unname(w$statistic), p_wilcoxon = w$p.value)
}) |> bind_rows()

d <- d |> mutate(ercc_c1_pass = ercc_reads_annotated >= ERCC_MIN_ANNOTATED)
tab <- table(d$group, d$ercc_c1_pass)
fis <- fisher.test(tab)
pass <- tibble(group = rownames(tab), pass = tab[, "TRUE"], fail = tab[, "FALSE"]) |>
  mutate(n = pass + fail, pct_pass = round(100 * pass / n, 1))

write_tsv(res, file.path(out, "mirna_yield_by_group.tsv"))
write_tsv(pass, file.path(out, "ercc_c1_pass_by_group.tsv"))
write_tsv(d |> select(sample_id, group, all_of(names(metrics)), ercc_c1_pass) |> arrange(group, sample_id),
          file.path(out, "mirna_yield_per_sample.tsv"))
sink(file.path(out, "mirna_yield_stats.txt"))
cat("Small-RNA yield by group, all samples (n =", nrow(d), ")\n")
cat("Wilcoxon rank-sum (exact) on log10(x+1); HL = Hodges-Lehmann fold change DNAemia/control\n\n")
print(as.data.frame(res |> mutate(across(where(is.numeric), ~ signif(.x, 3)))))
cat("\nERCC C1 pass (>=", ERCC_MIN_ANNOTATED, " annotated reads) by group:\n", sep = "")
print(tab); cat("Fisher exact p =", signif(fis$p.value, 3), "; OR =", signif(fis$estimate, 3),
                "(95% CI", signif(fis$conf.int[1], 3), "-", signif(fis$conf.int[2], 3), ")\n")
sink()

# supplementary figure: per-sample yield, log scale, threshold marked
long <- bind_rows(lapply(names(metrics)[c(1, 3, 4)], function(m)
  tibble(sample_id = d$sample_id, group = d$group, metric = metrics[[m]], value = d[[m]])))
long$metric <- factor(long$metric, levels = metrics[c(1, 3, 4)])
pl <- res |> filter(metric %in% names(metrics)[c(1, 3, 4)]) |>
  mutate(metric = factor(metrics[metric], levels = metrics[c(1, 3, 4)]),
         lab = paste0("p = ", formatC(p_wilcoxon, format = "g", digits = 2)))
p <- ggplot(long, aes(group, value)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey40") +
  geom_point(aes(colour = group), position = position_jitter(width = 0.12, seed = 1), size = 2) +
  geom_hline(data = tibble(metric = factor(metrics[[3]], levels = metrics[c(1, 3, 4)]), y = ERCC_MIN_ANNOTATED),
             aes(yintercept = y), linetype = "dashed", colour = "grey30") +
  geom_text(data = pl, aes(x = 1.5, y = Inf, label = lab), vjust = 1.5, size = 3.3, inherit.aes = FALSE) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_colour_manual(values = c(control = "#4C78A8", DNAemia = "#E45756"), guide = "none") +
  facet_wrap(~metric, scales = "free_y", nrow = 1, labeller = label_wrap_gen(width = 28)) +
  labs(x = NULL, y = "Reads (log scale)") +
  theme_bw(base_size = 11) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(out, "FigS_mirna_yield.pdf"), p, width = 9, height = 3.6)
ggsave(file.path(out, "FigS_mirna_yield.png"), p, width = 9, height = 3.6, dpi = 300)

# main Figure 1A panel: ONE panel. y = human miRNA reads (pre-specified primary yield metric);
# filled = passed ERCC criterion 1 (used in Fig 1B), open = failed; all samples
f1 <- d |> mutate(qc = factor(ifelse(ercc_c1_pass, "Passed QC (used in B)", "Failed QC"),
                              levels = c("Passed QC (used in B)", "Failed QC")))
f1lab <- sprintf("p = %.2g", res$p_wilcoxon[res$metric == "mi_counts_hsa"])
pA <- ggplot(f1, aes(group, mi_counts_hsa)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey40") +
  geom_point(aes(colour = group, shape = qc), position = position_jitter(width = 0.12, seed = 1), size = 2.2, stroke = 0.9) +
  annotate("text", x = 1.5, y = Inf, label = f1lab, vjust = 1.5, size = 3.2) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                expand = expansion(mult = c(0.05, 0.12))) +
  scale_colour_manual(values = c(control = "#4C78A8", DNAemia = "#E45756"), guide = "none") +
  scale_shape_manual(values = c("Passed QC (used in B)" = 16, "Failed QC" = 1), name = NULL) +
  labs(x = NULL, y = "Human miRNA reads per sample (log scale)") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
saveRDS(pA, file.path(out, "Fig1A_yield.rds"))   # for analysis/07_main_figures.R
writeLines(c(
  "# miRNA library yield by group (generated by analysis/01_mirna_yield.R; do not edit by hand)",
  "",
  paste0("Input: `", basename(qc_f), "`, all ", nrow(d), " samples (no QC exclusion). ",
         "Question: do control urines yield less small RNA? Raw sequencing depth is the negative control."),
  "",
  "| File | Contents |", "|---|---|",
  "| `mirna_yield_by_group.tsv` | One row per metric (columns below). |",
  "| `ercc_c1_pass_by_group.tsv` | Samples passing / failing ERCC criterion 1 per group. |",
  "| `mirna_yield_per_sample.tsv` | Per-sample values of the yield metrics and `ercc_c1_pass`. |",
  "| `mirna_yield_stats.txt` | Full test output, incl. Fisher exact test on the pass rates. |",
  "| `FigS_mirna_yield.pdf/.png` | Supplementary figure: per-sample yield, log scale; dashed line = ERCC criterion 1. |",
  "",
  "## Metrics",
  "",
  "| Metric | Meaning |", "|---|---|",
  "| `mi_counts_hsa` | Sum of human mature miRNA counts (**primary yield metric**). |",
  "| `mi_reads_17_30nt` | Reads 17-30 nt after adapter trimming. |",
  "| `ercc_reads_annotated` | Reads on annotated human RNA (ERCC criterion 1 metric; GENCODE + tRNA + piRNA if run). |",
  "| `mi_reads_raw` | Raw sequencing reads: negative control (sequencing depth should not differ). |",
  "",
  "## Columns of mirna_yield_by_group.tsv",
  "",
  "| Column | Meaning |", "|---|---|",
  "| `median_control`, `median_DNAemia` | Group medians. |",
  "| `fold_median` | median_DNAemia / median_control. |",
  "| `hl_fold`, `hl_fold_lo`, `hl_fold_hi` | Hodges-Lehmann fold change DNAemia vs control with 95% CI (from the Wilcoxon test on log10(x+1)). |",
  "| `W`, `p_wilcoxon` | Exact two-sided Wilcoxon rank-sum test on log10(x+1). Not adjusted for multiple metrics; `mi_counts_hsa` is the pre-specified primary metric. |",
  "",
  paste0("ERCC criterion 1 = >=", format(ERCC_MIN_ANNOTATED, big.mark = ","),
         " reads on annotated human RNA (ERCC small RNA-seq QC standard v1.0, exrna.org; Murillo et al. Cell 2019).")
), file.path(out, "README.md"))
cat(readLines(file.path(out, "mirna_yield_stats.txt")), sep = "\n")

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
