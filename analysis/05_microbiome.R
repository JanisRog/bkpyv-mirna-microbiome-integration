#!/usr/bin/env Rscript
# Urinary bacterial microbiome, DNAemia vs control (Fig 1B).
#   input       : MetaPhlAn 4 estimated reads per species (Bacteria only), all samples passing metagenome QC
#   abundance   : species kept if mean relative abundance (within Bacteria) > 0.1% across samples
#                 (same filter as Bates et al., J Med Microbiol 2026, on this dataset)
#   cross-talk  : optional (sensitivity): a count is set to 0 if < XTALK x the species' maximum count in any
#                 sample (UNCROSS2 principle, Edgar 2018; index-hopping rate 0.1-2%, Costello et al. 2018)
#   alpha       : observed species richness, Shannon (Wilcoxon, exact)
#   beta        : Aitchison distance (CLR, pseudocount 0.5), PERMANOVA (adonis2, 9999 perm.) + dispersion test
#   DA          : ALDEx2 (128 Monte Carlo CLR instances, Wilcoxon + Welch, BH), genus and species level;
#                 same method/version as Bates et al. 2026 on this dataset
#   viral       : urinary BKPyV / JCPyV DNA (pairs per million QC pairs) by group (Wilcoxon), descriptive
#
# Usage: Rscript 05_microbiome.R <metaphlan_species_est_reads.tsv> <qc_combined.tsv> <out_dir> [xtalk_fraction=0] [exclude=S18,...]
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(vegan)
  library(ALDEx2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: 05_microbiome.R <est_reads.tsv> <qc_combined.tsv> <out_dir> [xtalk] [exclude]")
est_f <- args[1]; qc_f <- args[2]; out <- args[3]
XTALK   <- if (length(args) >= 4) as.numeric(args[4]) else 0
EXCLUDE <- if (length(args) >= 5 && nzchar(args[5])) strsplit(args[5], ",")[[1]] else character(0)
MIN_MEAN_RA <- 0.001
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

# ---- data ----------------------------------------------------------------------
qc <- read_tsv(qc_f, show_col_types = FALSE) |>
  filter(!is.na(mg_est_bacterial_reads), !sample_id %in% EXCLUDE) |>
  mutate(group = factor(group, levels = c("control", "DNAemia")))
est <- read_tsv(est_f, show_col_types = FALSE) |> filter(grepl("^k__Bacteria\\|", clade))
tax <- est$clade
cnt <- round(as.matrix(est[, qc$sample_id])); rownames(cnt) <- tax
cnt <- cnt[rowSums(cnt) > 0, , drop = FALSE]
get_rank <- function(x, r) sub(paste0(".*\\|", r, "__([^|]+).*"), "\\1", x)
n_species_raw <- nrow(cnt)

# ---- cross-talk filter (sensitivity only) ---------------------------------------------
n_xtalk_zeroed <- 0
empty <- character(0)
if (XTALK > 0) {
  mx <- apply(cnt, 1, max)
  low <- cnt > 0 & cnt < XTALK * mx
  n_xtalk_zeroed <- sum(low)
  cnt[low] <- 0
  cnt <- cnt[rowSums(cnt) > 0, , drop = FALSE]
}

# samples left without any bacterial counts (entire profile explained by cross-talk) are reported and dropped
empty <- colnames(cnt)[colSums(cnt) == 0]
if (length(empty)) {
  cnt <- cnt[, colSums(cnt) > 0, drop = FALSE]
  qc  <- qc |> filter(!sample_id %in% empty)
  cnt <- cnt[rowSums(cnt) > 0, , drop = FALSE]
}

# ---- abundance filter --------------------------------------------------------------
ra <- sweep(cnt, 2, pmax(colSums(cnt), 1), "/")
keep <- rowMeans(ra) > MIN_MEAN_RA
cnt <- cnt[keep, , drop = FALSE]
ra  <- sweep(cnt, 2, pmax(colSums(cnt), 1), "/")
sp_name <- get_rank(rownames(cnt), "s")
genus   <- get_rank(rownames(cnt), "g")
gcnt <- rowsum(cnt, genus)
gra  <- sweep(gcnt, 2, pmax(colSums(gcnt), 1), "/")

# ---- alpha diversity -----------------------------------------------------------------
alpha <- tibble(sample_id = qc$sample_id, group = qc$group,
                observed = colSums(cnt > 0), shannon = diversity(t(ra), index = "shannon"),
                bacterial_reads = qc$mg_est_bacterial_reads)
# depth-matched richness: rarefy to the smallest bacterial read count, mean over 100 draws
rare_depth <- min(colSums(cnt))
alpha$observed_rarefied <- rowMeans(sapply(1:100, function(i) { set.seed(i); rowSums(rrarefy(t(cnt), rare_depth) > 0) }))
alpha_test <- bind_rows(lapply(c("observed", "observed_rarefied", "shannon"), function(m) {
  w <- suppressWarnings(wilcox.test(alpha[[m]] ~ alpha$group, exact = TRUE))
  tibble(metric = m, median_control = median(alpha[[m]][alpha$group == "control"]),
         median_DNAemia = median(alpha[[m]][alpha$group == "DNAemia"]), W = unname(w$statistic), p = w$p.value)
}))
depth_rho <- suppressWarnings(cor.test(alpha$observed, alpha$bacterial_reads, method = "spearman"))

# ---- beta diversity (Aitchison) ---------------------------------------------------------
clr <- function(m) { l <- log(m + 0.5); sweep(l, 2, colMeans(l), "-") }
cl <- t(clr(cnt))
d  <- dist(cl)
perm <- adonis2(d ~ group, data = data.frame(group = qc$group), permutations = 9999)
disp <- permutest(betadisper(d, qc$group), permutations = 9999)
pca  <- prcomp(cl, center = TRUE, scale. = FALSE)
ve   <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pcs  <- tibble(sample_id = qc$sample_id, group = qc$group, PC1 = pca$x[, 1], PC2 = pca$x[, 2])

# ---- differential abundance (ALDEx2) -----------------------------------------------------------
run_aldex <- function(m, level) {
  set.seed(1)
  # ALDEx2 orders conditions alphabetically ("DNAemia" < "control"), so diff.btw/effect = control - DNAemia;
  # signs are flipped below so that positive = higher in DNAemia (checked against raw abundances)
  x <- aldex(as.data.frame(m), conditions = as.character(qc$group), mc.samples = 128, test = "t",
             effect = TRUE, denom = "all", verbose = FALSE)
  lv <- sort(unique(as.character(qc$group)))
  sgn <- if (identical(lv, c("DNAemia", "control"))) -1 else 1
  tibble(level = level, taxon = rownames(x), diff_btw = sgn * x$diff.btw, effect = sgn * x$effect,
         p_wilcoxon = x$wi.ep, q_wilcoxon = x$wi.eBH, p_welch = x$we.ep, q_welch = x$we.eBH,
         significant = x$wi.eBH < 0.05) |> arrange(p_wilcoxon)
}
rownames(cnt) <- make.unique(sp_name)
da <- bind_rows(run_aldex(gcnt, "genus"), run_aldex(cnt, "species"))

# ---- urinary polyomavirus DNA ----------------------------------------------------------------
vir <- qc |> transmute(sample_id, group, bk_per_million = mg_bk_per_million_qc,
                       jc_per_million = 1e6 * mg_pairs_jc / mg_pairs_after_fastp)
vir_test <- bind_rows(lapply(c("bk_per_million", "jc_per_million"), function(m) {
  w <- suppressWarnings(wilcox.test(vir[[m]] ~ vir$group, exact = TRUE))
  tibble(metric = m, median_control = median(vir[[m]][vir$group == "control"]),
         median_DNAemia = median(vir[[m]][vir$group == "DNAemia"]), p = w$p.value)
}))

# ---- write tables ----------------------------------------------------------------------------
write_tsv(alpha, file.path(out, "alpha_diversity.tsv")); write_tsv(alpha_test, file.path(out, "alpha_tests.tsv"))
write_tsv(tibble(test = c("PERMANOVA (adonis2) group", "Dispersion (betadisper permutest)"),
                 F = c(perm$F[1], disp$tab$F[1]), R2 = c(perm$R2[1], NA), p = c(perm$`Pr(>F)`[1], disp$tab$`Pr(>F)`[1])),
          file.path(out, "beta_tests.tsv"))
write_tsv(pcs, file.path(out, "pca_aitchison.tsv"))
write_tsv(da, file.path(out, "differential_abundance_aldex2.tsv"))
write_tsv(as_tibble(gra, rownames = "genus"), file.path(out, "genus_relative_abundance.tsv"))
write_tsv(as_tibble(ra |> `rownames<-`(make.unique(sp_name)), rownames = "species"), file.path(out, "species_relative_abundance.tsv"))
write_tsv(vir, file.path(out, "polyomavirus_dna.tsv")); write_tsv(vir_test, file.path(out, "polyomavirus_tests.tsv"))

summ <- c(
  sprintf("samples: %d (DNAemia %d, control %d)%s", nrow(qc), sum(qc$group == "DNAemia"), sum(qc$group == "control"),
          if (length(EXCLUDE)) paste0(" | excluded: ", paste(EXCLUDE, collapse = ",")) else ""),
  sprintf("bacterial species (MetaPhlAn): %d | cross-talk filter: %s (%d counts set to 0) | after mean RA > 0.1%%: %d species, %d genera",
          n_species_raw, if (XTALK > 0) paste0(100 * XTALK, "% of max") else "off", n_xtalk_zeroed, nrow(cnt), nrow(gcnt)),
  if (length(empty)) sprintf("samples with no bacterial counts left after the cross-talk filter (dropped): %s", paste(empty, collapse = ", ")) else "samples dropped for no remaining counts: none",
  "", "Alpha diversity (Wilcoxon):", capture.output(print(as.data.frame(alpha_test))),
  sprintf("observed richness vs bacterial reads: Spearman rho = %.2f, p = %.2g", depth_rho$estimate, depth_rho$p.value),
  sprintf("observed_rarefied = richness after rarefying to %d bacterial reads (mean of 100 draws)", rare_depth),
  "", sprintf("Beta (Aitchison): PERMANOVA R2 = %.3f, F = %.2f, p = %.3g | dispersion p = %.3g",
              perm$R2[1], perm$F[1], perm$`Pr(>F)`[1], disp$tab$`Pr(>F)`[1]),
  "", "ALDEx2 (top 10 per level; diff_btw = median CLR difference DNAemia - control):",
  capture.output(print(as.data.frame(da |> group_by(level) |> slice_head(n = 10) |> ungroup() |>
    transmute(level, taxon, diff_btw = round(diff_btw, 2), effect = round(effect, 2),
              p_wilcoxon = signif(p_wilcoxon, 2), q_wilcoxon = signif(q_wilcoxon, 2), q_welch = signif(q_welch, 2))))),
  sprintf("ALDEx2 BH q < 0.05 (Wilcoxon): genus %d, species %d", sum(da$significant[da$level == "genus"], na.rm = TRUE),
          sum(da$significant[da$level == "species"], na.rm = TRUE)),
  "", "Urinary polyomavirus DNA (pairs per million QC pairs; Wilcoxon):", capture.output(print(as.data.frame(vir_test))))
writeLines(summ, file.path(out, "microbiome_summary.txt"))

# ---- Fig 1B ------------------------------------------------------------------------------------
cols <- c(control = "#4C78A8", DNAemia = "#E45756")
top <- names(sort(rowMeans(gra), decreasing = TRUE))[1:min(12, nrow(gra))]
bar <- as_tibble(gra, rownames = "genus") |> pivot_longer(-genus, names_to = "sample_id", values_to = "ra") |>
  mutate(genus = ifelse(genus %in% top, genus, "Other")) |> group_by(sample_id, genus) |> summarise(ra = sum(ra), .groups = "drop") |>
  left_join(qc |> dplyr::select(sample_id, group), by = "sample_id") |>
  mutate(genus = factor(genus, levels = c(top, "Other")),
         sample_id = factor(sample_id, levels = qc$sample_id[order(qc$group, -gra[top[1], qc$sample_id])]))
pal12 <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
           "#9C755F", "#17BECF", "#8C564B", "#BCBD22")
pal <- c(pal12[seq_along(top)], "grey85"); names(pal) <- c(top, "Other")
pA <- ggplot(bar, aes(sample_id, 100 * ra, fill = genus)) + geom_col(width = 0.9) +
  facet_grid(~group, scales = "free_x", space = "free_x") + scale_fill_manual(values = pal, name = "Genus") +
  labs(x = NULL, y = "Relative abundance (%)") + theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7), strip.background = element_blank(),
        strip.text = element_text(face = "bold"), legend.text = element_text(face = "italic", size = 8),
        legend.key.size = unit(0.35, "cm"))
alab <- c("Richness", "Richness\n(rarefied)", "Shannon")
al <- alpha |> dplyr::select(sample_id, group, observed, observed_rarefied, shannon) |>
  pivot_longer(c(observed, observed_rarefied, shannon)) |>
  mutate(name = factor(alab[match(name, c("observed", "observed_rarefied", "shannon"))], levels = alab))
al_lab <- alpha_test |> mutate(name = factor(alab, levels = alab), lab = sprintf("p = %.2g", p))
pB <- ggplot(al, aes(group, value, colour = group)) + geom_boxplot(outlier.shape = NA, colour = "grey40", width = 0.5) +
  geom_point(position = position_jitter(width = 0.1, seed = 1), size = 1.8) +
  geom_text(data = al_lab, aes(x = 1.5, y = Inf, label = lab), vjust = 1.4, size = 3, inherit.aes = FALSE) +
  facet_wrap(~name, scales = "free_y") + scale_colour_manual(values = cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = "Alpha diversity") + theme_bw(base_size = 10) +
  theme(strip.background = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
pC <- ggplot(pcs, aes(PC1, PC2, colour = group)) + geom_point(size = 2.3) +
  stat_ellipse(level = 0.68, linetype = "dashed", show.legend = FALSE) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2]),
       subtitle = bquote("Aitchison; PERMANOVA" ~ R^2 == .(sprintf("%.2f", perm$R2[1])) * "," ~ italic(p) == .(sprintf("%.2g", perm$`Pr(>F)`[1])))) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
fig <- pA / wrap_plots(pB, pC, widths = c(1.25, 1)) + plot_layout(heights = c(1.1, 1)) + plot_annotation(tag_levels = "i")
saveRDS(list(bars = pA, alpha = pB, pca = pC), file.path(out, "fig_panels.rds"))   # used by 07_main_figures.R
ggsave(file.path(out, "Fig1B_microbiome.pdf"), fig, width = 9, height = 8)
ggsave(file.path(out, "Fig1B_microbiome.png"), fig, width = 9, height = 8, dpi = 300)

# ---- README + session ----------------------------------------------------------------------------
writeLines(c(
  "# Urinary bacterial microbiome (generated by analysis/05_microbiome.R; do not edit by hand)", "",
  sprintf("Samples: %d passing metagenome QC%s. Comparison DNAemia vs control (positive LFC = higher in DNAemia).",
          nrow(qc), if (length(EXCLUDE)) paste0(", excluding ", paste(EXCLUDE, collapse = ", ")) else ""),
  sprintf("Cross-talk filter: %s.", if (XTALK > 0) sprintf("on, counts < %g%% of the species' maximum set to 0 (UNCROSS2 principle, Edgar 2018; hopping rate Costello et al. 2018)", 100 * XTALK) else "off"),
  "Abundance filter: species with mean relative abundance (within Bacteria) > 0.1% (as in Bates et al., J Med Microbiol 2026).", "",
  "| File | Contents |", "|---|---|",
  "| `microbiome_summary.txt` | All test results in one place. |",
  "| `alpha_diversity.tsv`, `alpha_tests.tsv` | Observed species richness (raw and rarefied to the smallest bacterial read count, mean of 100 draws) and Shannon index per sample; exact Wilcoxon tests. |",
  "| `beta_tests.tsv` | PERMANOVA on Aitchison distance (9999 permutations) and test of equal dispersion. |",
  "| `pca_aitchison.tsv` | PCA coordinates of CLR-transformed counts. |",
  "| `differential_abundance_aldex2.tsv` | ALDEx2 per genus and species: `diff_btw` (median CLR difference, DNAemia - control; positive = higher in DNAemia), `effect` (ALDEx2 standardised effect size), Wilcoxon and Welch expected p and BH q (`q_wilcoxon` used for `significant`). |",
  "| `genus_relative_abundance.tsv`, `species_relative_abundance.tsv` | Relative abundance within Bacteria after filtering (fractions). |",
  "| `polyomavirus_dna.tsv`, `polyomavirus_tests.tsv` | Urinary BKPyV and JCPyV DNA, read pairs per million QC pairs; Wilcoxon by group. |",
  "| `Fig1B_microbiome.pdf/.png` | (i) genus composition per sample, (ii) alpha diversity, (iii) Aitchison PCA. |",
  "| `fig_panels.rds` | The three ggplot panels, used by `analysis/07_main_figures.R`. |", "",
  "Counts are MetaPhlAn 4 estimated reads per species (rel_ab_w_read_stats). Richness in low-biomass urine depends on bacterial read depth; see the depth correlation in the summary."
), file.path(out, "README.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat(summ, sep = "\n")
