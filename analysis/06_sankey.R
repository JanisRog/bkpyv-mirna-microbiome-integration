#!/usr/bin/env Rscript
# Fig 2: descriptive three-column Sankey  genus -> sample -> miRNA  (paired samples only).
# Links show shared sample membership, NOT associations; no test is performed here.
#   left  width : genus relative abundance within Bacteria in that sample (%), after the mean-RA > 0.1% filter
#   right width : |log2(normalised count + 1) - control mean| for that sample, shown if >= 1 (2-fold)
#                 colour = above (up) or below (down) the control mean; viral miRNAs in purple
#   miRNAs      : the 8 human miRNAs with the smallest DESeq2 p (Fig 1A) + bkv-miR-B1-5p
#   node order  : genera and miRNAs ordered by the mean position of their samples (fewer crossings)
#   middle      : samples grouped control / DNAemia; dot = urinary BKPyV DNA (metagenome)
# Concept follows the descriptive Sankey of the first submission, simplified (genus level, fewer elements).
#
# Usage: Rscript 06_sankey.R <genus_relative_abundance.tsv> <de_dir> <qc_combined.tsv> <out_dir>
suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: 06_sankey.R <genus_relative_abundance.tsv> <de_dir> <qc_combined.tsv> <out_dir>")
gra_f <- args[1]; de_dir <- args[2]; qc_f <- args[3]; out <- args[4]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
N_GENERA <- 10; N_MIRNA <- 8; MIN_DELTA <- 1; VIRAL <- "bkv-miR-B1-5p"   # MIN_DELTA 1 = 2-fold from control mean

# ---- inputs ----------------------------------------------------------------------------
gra <- read_tsv(gra_f, show_col_types = FALSE)
samples_used <- read_tsv(file.path(de_dir, "samples_used.tsv"), show_col_types = FALSE)
qc <- read_tsv(qc_f, show_col_types = FALSE) |> filter(sample_id %in% samples_used$sample_id) |>
  mutate(group = factor(group, levels = c("control", "DNAemia")))
stopifnot(setequal(qc$sample_id, intersect(names(gra), qc$sample_id)))
ord <- qc |> arrange(group, as.integer(sub("S", "", sample_id))) |> pull(sample_id)

de  <- read_tsv(file.path(de_dir, "de_human_mirna.tsv"), show_col_types = FALSE) |> arrange(p_deseq2)
mir <- c(head(de$miRNA[!is.na(de$p_deseq2)], N_MIRNA), VIRAL)
norm <- read_tsv(file.path(de_dir, "normalized_counts_deseq2.tsv"), show_col_types = FALSE)
lg <- as.matrix(log2(norm[match(mir, norm$miRNA), ord] + 1)); rownames(lg) <- mir
ctrl <- ord[qc$group[match(ord, qc$sample_id)] == "control"]
delta <- lg - rowMeans(lg[, ctrl, drop = FALSE])
w <- abs(delta); w[w < MIN_DELTA] <- 0

g <- as.matrix(gra[, ord]) * 100; rownames(g) <- gra$genus
top <- names(sort(rowMeans(g), decreasing = TRUE))[seq_len(min(N_GENERA, nrow(g)))]
g <- rbind(g[top, , drop = FALSE], Other = colSums(g[setdiff(rownames(g), top), , drop = FALSE]))

# ---- geometry --------------------------------------------------------------------------------
X_G <- 0; X_S <- 1; X_M <- 2; NODE_W <- 0.035
ns <- length(ord)
is_case <- qc$group[match(ord, qc$sample_id)] == "DNAemia"
s_y <- setNames(seq(0.93, 0.04, length.out = ns), ord); s_y[is_case] <- s_y[is_case] - 0.05
# barycentre ordering: place each genus / miRNA at the weighted mean height of its samples
bary <- function(m) apply(m, 1, function(r) if (sum(r) > 0) sum(r * s_y[colnames(m)]) / sum(r) else 0)
g_ord <- c(names(sort(bary(g[rownames(g) != "Other", , drop = FALSE]), decreasing = TRUE)), "Other")
g <- g[g_ord, , drop = FALSE]
h_ord <- setdiff(mir, VIRAL); h_ord <- c(names(sort(bary(w[h_ord, , drop = FALSE]), decreasing = TRUE)), VIRAL)
w <- w[h_ord, , drop = FALSE]; delta <- delta[h_ord, , drop = FALSE]; mir <- h_ord
L_SCALE <- 0.055 / 100                                   # 100 % = sample node height on the left
R_SCALE <- min(0.055 / max(colSums(w)), 0.07 / max(rowSums(w)))
stack <- function(h, top = 0.97, bottom = 0.03) {          # centred, evenly gapped nodes
  gap <- (top - bottom - sum(h)) / max(length(h) - 1, 1)
  y <- top - cumsum(c(0, head(h, -1) + gap)) - h / 2; setNames(y, names(h))
}
g_h <- rowSums(g) * L_SCALE; g_y <- stack(g_h)
m_h <- rowSums(w) * R_SCALE; m_y <- stack(m_h, top = 0.95, bottom = 0.05)

ribbon <- function(x0, x1, y0a, y0b, y1a, y1b, id) {      # smooth band between two vertical segments
  t <- seq(0, 1, length.out = 40); s <- 3 * t^2 - 2 * t^3; x <- x0 + (x1 - x0) * t
  tibble(id = id, x = c(x, rev(x)), y = c(y0a + (y1a - y0a) * s, rev(y0b + (y1b - y0b) * s)))
}
polys <- list(); meta <- list(); k <- 0
# left: genus -> sample (each sample's 100 % stacked top-down in genus order)
g_cur <- g_y + g_h / 2; s_cur <- s_y + 100 * L_SCALE / 2
for (gn in rownames(g)) for (s in ord) {
  v <- g[gn, s]; if (v <= 0) next
  h <- v * L_SCALE; k <- k + 1
  polys[[k]] <- ribbon(X_G + NODE_W / 2, X_S - NODE_W / 2, g_cur[gn], g_cur[gn] - h, s_cur[s], s_cur[s] - h, k)
  meta[[k]] <- tibble(id = k, side = "genus", source = gn, target = s, value = v, fill = gn)
  g_cur[gn] <- g_cur[gn] - h; s_cur[s] <- s_cur[s] - h
}
# right: sample -> miRNA
s_cur <- s_y + colSums(w) * R_SCALE / 2; m_cur <- m_y + m_h / 2
for (s in ord) for (m in mir) {
  v <- w[m, s]; if (v <= 0) next
  h <- v * R_SCALE; k <- k + 1
  dir <- if (delta[m, s] > 0) "up" else "down"
  fill <- if (m == VIRAL) paste0("viral_", dir) else paste0("human_", dir)
  polys[[k]] <- ribbon(X_S + NODE_W / 2, X_M - NODE_W / 2, s_cur[s], s_cur[s] - h, m_cur[m], m_cur[m] - h, k)
  meta[[k]] <- tibble(id = k, side = "miRNA", source = s, target = m, value = v, signed_delta = delta[m, s], fill = fill)
  s_cur[s] <- s_cur[s] - h; m_cur[m] <- m_cur[m] - h
}
P <- bind_rows(polys) |> left_join(bind_rows(meta), by = "id")
links <- bind_rows(meta)

# ---- colours ---------------------------------------------------------------------------------------
pal12 <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
           "#9C755F", "#17BECF", "#8C564B", "#BCBD22")
gcol <- c(setNames(pal12[seq_along(top)], top), Other = "grey80")
dcol <- c(human_up = "#D6604D", human_down = "#4393C3", viral_up = "#7B3294", viral_down = "#C2A5CF")
fills <- c(gcol, dcol)
grp_col <- c(control = "#4C78A8", DNAemia = "#E45756")

nodes <- bind_rows(
  tibble(x = X_G, y = g_y, h = g_h, fill = "grey35"),
  tibble(x = X_S, y = s_y, h = 100 * L_SCALE, fill = grp_col[as.character(qc$group[match(ord, qc$sample_id)])]),
  tibble(x = X_M, y = m_y, h = pmax(m_h, 0.004), fill = ifelse(mir == VIRAL, "#7B3294", "grey35")))
bk <- qc |> transmute(sample_id, bk = log10(mg_bk_per_million_qc + 1)) |> mutate(y = s_y[sample_id])
lab_g <- tibble(x = X_G - NODE_W, y = g_y, lab = names(g_y), it = names(g_y) != "Other")
lab_m <- tibble(x = X_M + NODE_W, y = m_y, lab = sub("^hsa-", "", mir))

p <- ggplot() +
  geom_polygon(data = P, aes(x, y, group = id, fill = fill), alpha = 0.55, colour = NA) +
  geom_rect(data = nodes, aes(xmin = x - NODE_W / 2, xmax = x + NODE_W / 2, ymin = y - h / 2, ymax = y + h / 2),
            fill = nodes$fill, colour = NA) +
  geom_text(data = lab_g, aes(x, y, label = lab, fontface = ifelse(it, "italic", "plain")), hjust = 1, size = 3) +
  geom_text(data = lab_m, aes(x, y, label = lab), hjust = 0, size = 3) +
  geom_label(data = tibble(x = X_S, y = s_y, lab = ord), aes(x, y, label = lab), size = 2.7,
             label.size = 0, label.padding = unit(0.08, "lines"), fill = "white") +
  geom_point(data = bk, aes(x = X_S + NODE_W * 1.6, y = y, colour = bk), size = 2.4) +
  annotate("text", x = c(X_G, X_S, X_M), y = 1.04, label = c("Genus", "Sample", "miRNA"), fontface = "bold", size = 3.6) +
  annotate("text", x = X_S, y = max(s_y[!is_case]) + 0.045, label = "Control",
           colour = grp_col[["control"]], size = 3.2, fontface = "bold") +
  annotate("text", x = X_S, y = max(s_y[is_case]) + 0.045, label = "DNAemia",
           colour = grp_col[["DNAemia"]], size = 3.2, fontface = "bold") +
  # scale bars
  annotate("rect", xmin = 0.25, xmax = 0.29, ymin = -0.06, ymax = -0.06 + 25 * L_SCALE, fill = "grey40") +
  annotate("text", x = 0.31, y = -0.06 + 25 * L_SCALE / 2, label = "25% of bacteria", hjust = 0, size = 2.7) +
  annotate("rect", xmin = 1.25, xmax = 1.29, ymin = -0.06, ymax = -0.06 + 2 * R_SCALE, fill = "grey40") +
  annotate("text", x = 1.31, y = -0.06 + 2 * R_SCALE / 2, label = "4-fold (2 log2 units) from control mean", hjust = 0, size = 2.7) +
  scale_fill_manual(values = fills, breaks = names(dcol),
                    labels = c("Human, above control mean", "Human, below control mean",
                               "bkv-miR-B1-5p, above control mean", "bkv-miR-B1-5p, below control mean"), name = NULL) +
  scale_colour_gradient(low = "grey90", high = "black", name = expression("Urinary BKPyV DNA (log"[10]*" pairs per million)")) +
  coord_cartesian(xlim = c(-0.45, 2.5), ylim = c(-0.08, 1.05), clip = "off") +
  guides(fill = guide_legend(nrow = 2, override.aes = list(alpha = 0.8)), colour = guide_colourbar(barheight = 0.5, barwidth = 6)) +
  theme_void(base_size = 10) + theme(legend.position = "bottom", legend.box = "vertical", plot.margin = margin(10, 10, 5, 10),
                                     plot.background = element_rect(fill = "white", colour = NA))
ggsave(file.path(out, "Fig2_sankey.pdf"), p, width = 8, height = 7.5, bg = "white")
ggsave(file.path(out, "Fig2_sankey.png"), p, width = 8, height = 7.5, dpi = 300, bg = "white")

write_tsv(links, file.path(out, "sankey_links.tsv"))
write_tsv(tibble(miRNA = mir, rank_by_p = match(mir, de$miRNA), log2FC = de$log2FC[match(mir, de$miRNA)],
                 p_deseq2 = de$p_deseq2[match(mir, de$miRNA)]), file.path(out, "sankey_mirnas.tsv"))
stopifnot(all(abs(tapply(links$value[links$side == "genus"], links$target[links$side == "genus"], sum) - 100) < 1e-6))
writeLines(c(
  "# Fig 2: descriptive Sankey (generated by analysis/06_sankey.R; do not edit by hand)", "",
  sprintf("Samples: the %d paired samples (miRNA analysis set; all have metagenomes). No statistical test is performed; links show shared sample membership, not associations.", ns), "",
  "| File | Contents |", "|---|---|",
  "| `Fig2_sankey.pdf/.png` | The figure. |",
  "| `sankey_links.tsv` | Every ribbon: `side` (genus or miRNA), `source`, `target`, `value` (width), `signed_delta` (miRNA side), `fill`. |",
  "| `sankey_mirnas.tsv` | miRNAs shown: the 8 human miRNAs with the smallest DESeq2 p in Fig 1A, plus bkv-miR-B1-5p. |", "",
  "Left ribbon width: genus relative abundance within Bacteria in that sample (%), after the mean relative abundance > 0.1% filter; top 10 genera by mean abundance, the rest pooled as Other. Each sample sums to 100%.",
  sprintf("Right ribbon width: |log2(DESeq2-normalised count + 1) - mean of the controls| for that sample; only values >= %.1f are drawn. Colour: above or below the control mean; viral miRNA in purple.", MIN_DELTA),
  "The two sides use independent scales (see scale bars). Dot next to each sample: urinary BKPyV DNA from the metagenome."
), file.path(out, "README.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("links:", nrow(links), "| genera:", nrow(g), "| miRNAs:", length(mir), "\n")
