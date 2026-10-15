#!/usr/bin/env Rscript
# Differential miRNA expression, DNAemia vs control (Fig 1A).
#   human miRNAs : DESeq2 Wald (primary) + edgeR quasi-likelihood (robustness); BH-FDR within human miRNAs
#   viral miRNAs : bkv-miR-B1-5p (BKPyV-specific), bkv-miR-B1-3p (= jcv-miR-J1-3p), jcv-miR-J1-5p;
#                  tested in the same DESeq2 model but as a SEPARATE family (own BH adjustment),
#                  plus exact Wilcoxon on CPM as a distribution-free check
#   size factors : from human miRNAs only (viral counts must not drive normalisation)
#
# Sample sets:
#   primary     = ERCC criterion 1 pass (>=100k reads on annotated human RNA)
#   sensitivity = >=10,000 human miRNA counts (looser, not guideline-based)
#
# Usage: Rscript 02_mirna_de.R <mirna_counts_raw.tsv> <qc_combined.tsv> <primary|sensitivity> <out_dir>
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
  library(DESeq2); library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: 02_mirna_de.R <counts.tsv> <qc_combined.tsv> <primary|sensitivity> <out_dir>")
counts_f <- args[1]; qc_f <- args[2]; set_name <- args[3]; out <- args[4]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

ERCC_MIN_ANNOTATED <- 100000
SENS_MIN_HSA       <- 10000
FDR                <- 0.05
VIRAL <- c("bkv-miR-B1-5p", "bkv-miR-B1-3p", "jcv-miR-J1-5p")   # jcv-miR-J1-3p is identical to bkv-miR-B1-3p
VIRAL_LABEL <- c("bkv-miR-B1-5p" = "bkv-miR-B1-5p (BKPyV-specific)",
                 "bkv-miR-B1-3p" = "bkv-miR-B1-3p (= jcv-miR-J1-3p)",
                 "jcv-miR-J1-5p" = "jcv-miR-J1-5p (JCPyV-specific)")

# ---- samples ------------------------------------------------------------------
qc <- read_tsv(qc_f, show_col_types = FALSE)
keep <- switch(set_name,
  primary     = qc$ercc_reads_annotated >= ERCC_MIN_ANNOTATED,
  sensitivity = qc$mi_counts_hsa >= SENS_MIN_HSA,
  stop("set must be 'primary' or 'sensitivity'"))
cd <- qc[keep, c("sample_id", "group")] |> as.data.frame()
cd$group <- factor(cd$group, levels = c("control", "DNAemia"))
rownames(cd) <- cd$sample_id

cnt <- read_tsv(counts_f, show_col_types = FALSE)
m <- as.matrix(cnt[, cd$sample_id]); rownames(m) <- cnt$miRNA
storage.mode(m) <- "integer"
hsa <- m[grepl("^hsa-", rownames(m)), , drop = FALSE]
vir <- m[VIRAL, , drop = FALSE]

# ---- filter (edgeR filterByExpr on human miRNAs) ---------------------------------
y0 <- DGEList(hsa, group = cd$group)
keep_f <- filterByExpr(y0, group = cd$group)
hsa_f <- hsa[keep_f, , drop = FALSE]

# ---- DESeq2: human + viral rows, size factors from human only -------------------
dds_h <- DESeqDataSetFromMatrix(hsa_f, cd, design = ~ group)
dds_h <- estimateSizeFactors(dds_h, type = "poscounts")
sf <- sizeFactors(dds_h)
dds <- DESeqDataSetFromMatrix(rbind(hsa_f, vir), cd, design = ~ group)
sizeFactors(dds) <- sf
dds <- DESeq(dds, quiet = TRUE)
coef <- "group_DNAemia_vs_control"
res <- results(dds, name = coef, independentFiltering = FALSE, cooksCutoff = TRUE)
shr <- lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)

de <- tibble(miRNA = rownames(res), baseMean = res$baseMean,
             log2FC = res$log2FoldChange, lfcSE = res$lfcSE,
             log2FC_lo = res$log2FoldChange - qnorm(0.975) * res$lfcSE,
             log2FC_hi = res$log2FoldChange + qnorm(0.975) * res$lfcSE,
             log2FC_shrunk = shr$log2FoldChange, p_deseq2 = res$pvalue)

hum <- de |> filter(grepl("^hsa-", miRNA)) |> mutate(q_deseq2 = p.adjust(p_deseq2, "BH"))
vi  <- de |> filter(miRNA %in% VIRAL)    |> mutate(q_deseq2 = p.adjust(p_deseq2, "BH"))

# ---- edgeR QL (robustness, human only) -----------------------------------------
y <- normLibSizes(DGEList(hsa_f, group = cd$group))
design <- model.matrix(~ group, data = cd)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = "groupDNAemia")
er <- topTags(qlf, n = Inf, sort.by = "none")$table
er <- tibble(miRNA = rownames(er), log2FC_edger = er$logFC, p_edger = er$PValue, q_edger = er$FDR)

hum <- hum |> left_join(er, by = "miRNA") |>
  mutate(direction = ifelse(log2FC > 0, "up in DNAemia", "down in DNAemia"),
         sig_deseq2 = !is.na(q_deseq2) & q_deseq2 < FDR,
         sig_edger  = !is.na(q_edger) & q_edger < FDR,
         robust     = sig_deseq2 & sig_edger & sign(log2FC) == sign(log2FC_edger),
         nominal    = !is.na(p_deseq2) & p_deseq2 < 0.05) |>
  arrange(p_deseq2)

# ---- viral: Wilcoxon on CPM (per million human miRNA counts) --------------------
lib <- colSums(hsa)
vi <- vi |> rowwise() |> mutate(
  label = VIRAL_LABEL[[miRNA]],
  cpm_median_control = median(1e6 * vir[miRNA, cd$group == "control"] / lib[cd$group == "control"]),
  cpm_median_DNAemia = median(1e6 * vir[miRNA, cd$group == "DNAemia"] / lib[cd$group == "DNAemia"]),
  p_wilcoxon_cpm = wilcox.test(1e6 * vir[miRNA, ] / lib ~ cd$group, exact = TRUE)$p.value) |> ungroup()

# ---- outputs -------------------------------------------------------------------
write_tsv(hum, file.path(out, "de_human_mirna.tsv"))
write_tsv(vi,  file.path(out, "de_viral_mirna.tsv"))
write_tsv(tibble(sample_id = cd$sample_id, group = cd$group, size_factor = sf[cd$sample_id],
                 human_mirna_counts = colSums(hsa)), file.path(out, "samples_used.tsv"))
norm <- counts(dds, normalized = TRUE)
write_tsv(as_tibble(norm, rownames = "miRNA"), file.path(out, "normalized_counts_deseq2.tsv"))

summ <- c(
  sprintf("set: %s | samples: %d (DNAemia %d, control %d)", set_name, nrow(cd),
          sum(cd$group == "DNAemia"), sum(cd$group == "control")),
  sprintf("human miRNAs in reference: %d | tested after filterByExpr: %d", nrow(hsa), nrow(hsa_f)),
  sprintf("DESeq2 q<%.2f: %d (up %d, down %d) | edgeR q<%.2f: %d | robust (both, same direction): %d",
          FDR, sum(hum$sig_deseq2), sum(hum$sig_deseq2 & hum$log2FC > 0), sum(hum$sig_deseq2 & hum$log2FC < 0),
          FDR, sum(hum$sig_edger), sum(hum$robust)),
  sprintf("nominal p<0.05 (DESeq2): %d; expected by chance ~%.0f", sum(hum$nominal), 0.05 * nrow(hum)),
  "", "Top 15 human miRNAs (DESeq2):",
  capture.output(print(as.data.frame(hum |> head(15) |>
    transmute(miRNA, baseMean = round(baseMean), log2FC = round(log2FC, 2),
              CI = sprintf("[%.2f, %.2f]", log2FC_lo, log2FC_hi),
              p = signif(p_deseq2, 2), q = signif(q_deseq2, 2), q_edger = signif(q_edger, 2), robust)))),
  "", "Viral miRNAs (separate family):",
  capture.output(print(as.data.frame(vi |> transmute(miRNA, log2FC = round(log2FC, 2),
    p = signif(p_deseq2, 2), q = signif(q_deseq2, 2), cpm_median_control = round(cpm_median_control, 1),
    cpm_median_DNAemia = round(cpm_median_DNAemia, 1), p_wilcoxon_cpm = signif(p_wilcoxon_cpm, 2))))))
writeLines(summ, file.path(out, "de_summary.txt"))

# ---- volcano (Fig 1A) ---------------------------------------------------------------
# Colour: FDR-significant (q < 0.05) if any; otherwise nominal p < 0.05 shown in lighter colours and
# labelled as nominal (exploratory study). Viral miRNAs: triangles, separate test family.
lev <- c("Up, q < 0.05", "Down, q < 0.05", "Up, nominal p < 0.05", "Down, nominal p < 0.05",
         "Not significant", "Viral miRNA (BKPyV, JCPyV)")
pd <- bind_rows(
  hum |> transmute(miRNA, log2FC, p = p_deseq2, kind = case_when(
    sig_deseq2 & log2FC > 0 ~ lev[1], sig_deseq2 & log2FC < 0 ~ lev[2],
    nominal & log2FC > 0 ~ lev[3], nominal & log2FC < 0 ~ lev[4], TRUE ~ lev[5]),
    lab = ifelse(nominal, sub("^hsa-", "", miRNA), NA)),
  vi |> transmute(miRNA, log2FC, p = p_deseq2, kind = lev[6],
                  lab = ifelse(miRNA == "bkv-miR-B1-3p", "bkv-miR-B1-3p / jcv-miR-J1-3p", miRNA))) |>
  filter(!is.na(p)) |> mutate(kind = factor(kind, levels = lev))
present <- lev[lev %in% pd$kind & lev != "Not significant"]
cols <- c("#B2182B", "#2166AC", "#F4A582", "#92C5DE", "grey70", "#7B3294"); names(cols) <- lev
shps <- c(16, 16, 16, 16, 16, 17); names(shps) <- lev
v <- ggplot(pd, aes(log2FC, -log10(p))) +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "grey50") +
  { if (any(hum$sig_deseq2)) geom_hline(yintercept = -log10(max(hum$p_deseq2[hum$sig_deseq2])),
                                         linetype = "dashed", colour = "grey30") } +
  geom_point(data = ~ filter(.x, kind == "Not significant"), colour = "grey75", size = 1.2) +
  geom_point(data = ~ filter(.x, kind != "Not significant"), aes(colour = kind, shape = kind), size = 2.3) +
  geom_text_repel(aes(label = lab), size = 2.8, max.overlaps = 40, min.segment.length = 0, seed = 42, box.padding = 0.45, force = 3,
                  segment.colour = "grey60", na.rm = TRUE) +
  scale_colour_manual(values = cols[present], breaks = present, name = NULL) +
  scale_shape_manual(values = shps[present], breaks = present, name = NULL) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  labs(x = expression(log[2]~fold~change~(DNAemia~vs~control)), y = expression(-log[10]~italic(p)),
       caption = sprintf(paste0("%d DNAemia vs %d control; %d human miRNAs tested (DESeq2); ",
                                "%d with FDR q < 0.05. Dotted line: nominal p = 0.05."),
                         sum(cd$group == "DNAemia"), sum(cd$group == "control"), nrow(hum), sum(hum$sig_deseq2))) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), plot.caption = element_text(size = 8))
ggsave(file.path(out, "Fig1A_volcano.pdf"), v, width = 6.5, height = 5.8)
ggsave(file.path(out, "Fig1A_volcano.png"), v, width = 6.5, height = 5.8, dpi = 300)

# ---- README + session ---------------------------------------------------------------
writeLines(c(
  sprintf("# Differential miRNA expression, %s set (generated by analysis/02_mirna_de.R; do not edit by hand)", set_name), "",
  if (set_name == "primary") "Samples: ERCC criterion 1 pass (>=100,000 reads on annotated human RNA)." else
    "Samples: >=10,000 human miRNA counts (looser sensitivity set, not guideline-based).",
  "Comparison: DNAemia vs control (positive log2FC = higher in DNAemia).", "",
  "| File | Contents |", "|---|---|",
  "| `de_summary.txt` | Numbers of tested/significant miRNAs, top 15, viral results. |",
  "| `de_human_mirna.tsv` | All tested human miRNAs (columns below). |",
  "| `de_viral_mirna.tsv` | Viral miRNAs, tested as a separate family. |",
  "| `samples_used.tsv` | Samples, group, DESeq2 size factor, human miRNA counts. |",
  "| `normalized_counts_deseq2.tsv` | DESeq2-normalised counts (size factors from human miRNAs). |",
  "| `Fig1A_volcano.pdf/.png` | Volcano plot. Dark red/blue = FDR q < 0.05; light red/blue = nominal p < 0.05 only (labelled); triangles = viral miRNAs (BKPyV, JCPyV only; adjusted as a separate family). |", "",
  "## Columns of de_human_mirna.tsv / de_viral_mirna.tsv", "",
  "| Column | Meaning |", "|---|---|",
  "| `baseMean` | Mean DESeq2-normalised count over all samples. |",
  "| `log2FC`, `lfcSE`, `log2FC_lo`, `log2FC_hi` | DESeq2 maximum-likelihood log2 fold change, its standard error and Wald 95% CI. |",
  "| `log2FC_shrunk` | apeglm-shrunken log2FC (less noisy for low counts; for ranking/plotting). |",
  "| `p_deseq2`, `q_deseq2` | DESeq2 Wald p; Benjamini-Hochberg q within the family (human or viral). |",
  "| `log2FC_edger`, `p_edger`, `q_edger` | edgeR quasi-likelihood F-test (robustness check; human only). |",
  "| `sig_deseq2`, `sig_edger` | q < 0.05 in that method. |",
  "| `robust` | Significant in both methods with the same direction. |",
  "| `nominal` | DESeq2 p < 0.05 (unadjusted; to be labelled nominal). |",
  "| `cpm_median_*`, `p_wilcoxon_cpm` | Viral only: median counts per million human miRNA counts per group; exact Wilcoxon test. |", "",
  "Filtering: edgeR `filterByExpr` on human miRNAs. Size factors: DESeq2 `poscounts` on human miRNAs only.",
  "bkv-miR-B1-3p is sequence-identical to jcv-miR-J1-3p and cannot distinguish BKPyV from JCPyV; bkv-miR-B1-5p is BKPyV-specific."
), file.path(out, "README.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat(summ, sep = "\n")
