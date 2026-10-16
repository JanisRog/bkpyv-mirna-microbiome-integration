#!/usr/bin/env Rscript
# Fig 2: descriptive three-column Sankey  genus -> sample -> miRNA  (paired samples only).
# Links show shared sample membership, NOT associations; no test is performed here.
#   left  width : genus relative abundance within Bacteria in that sample (%), after the mean-RA > 0.1% filter;
#                 genera shown if >= 10% in at least one sample, the rest pooled as Other; colours grouped by phylum
#   right width : |log2(normalised count + 1) - control mean| for that sample, shown if >= 1 (2-fold)
#                 colour = above (up) or below (down) the control mean; viral miRNAs in purple
#   miRNAs      : the 8 human miRNAs with the smallest DESeq2 p (Fig 1A) + bkv-miR-B1-5p and bkv-miR-B1-3p
#                 (3p sequence identical to jcv-miR-J1-3p)
#   node order  : genera and miRNAs ordered by the mean position of their samples (fewer crossings)
#   middle      : one stacked bar per sample = its genus composition (ribbons flow into their segment),
#                 no outline (ribbons run seamlessly into it); grey strip right of it = urinary BKPyV DNA (metagenome)
# Concept follows the descriptive Sankey of the first submission, simplified (genus level, fewer elements).
#
# Usage: Rscript 06_sankey.R <genus_relative_abundance.tsv> <de_dir> <qc_combined.tsv> <out_dir> <metaphlan_species_est_reads.tsv> [legend.png]
#   legend.png (optional): legend image placed below the plot instead of the drawn legend
#   (the figure uses analysis/assets/sankey_legend.png, designed by the first author)
suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 5:6) stop("usage: 06_sankey.R <genus_relative_abundance.tsv> <de_dir> <qc_combined.tsv> <out_dir> <metaphlan_est_reads.tsv> [legend.png]")
gra_f <- args[1]; de_dir <- args[2]; qc_f <- args[3]; out <- args[4]; tax_f <- args[5]
legend_png <- if (length(args) == 6) args[6] else NA
dir.create(out, showWarnings = FALSE, recursive = TRUE)
MIN_GENUS_MAX <- 10; N_MIRNA <- 8; MIN_DELTA <- 1; VIRAL <- c("bkv-miR-B1-5p", "bkv-miR-B1-3p")   # MIN_DELTA 1 = 2-fold

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
top <- rownames(g)[apply(g, 1, max) >= MIN_GENUS_MAX]          # genera reaching >= 10% in any sample
g <- rbind(g[top, , drop = FALSE], Other = colSums(g[setdiff(rownames(g), top), , drop = FALSE]))

# ---- geometry --------------------------------------------------------------------------------
X_G <- 0; X_S <- 1; X_M <- 2; NODE_W <- 0.035; SNODE_W <- 0.08; BK_W <- 0.03   # sample bar and BK strip widths
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
polys <- list(); meta <- list(); segs <- list(); k <- 0
# left: genus -> sample (each sample's 100 % stacked top-down in genus order)
g_cur <- g_y + g_h / 2; s_cur <- s_y + 100 * L_SCALE / 2
for (gn in rownames(g)) for (s in ord) {
  v <- g[gn, s]; if (v <= 0) next
  h <- v * L_SCALE; k <- k + 1
  polys[[k]] <- ribbon(X_G + NODE_W / 2, X_S - SNODE_W / 2, g_cur[gn], g_cur[gn] - h, s_cur[s], s_cur[s] - h, k)
  meta[[k]] <- tibble(id = k, side = "genus", source = gn, target = s, value = v, fill = gn)
  segs[[k]] <- tibble(sample = s, genus = gn, ymax = s_cur[s], ymin = s_cur[s] - h)
  g_cur[gn] <- g_cur[gn] - h; s_cur[s] <- s_cur[s] - h
}
# right: sample -> miRNA
s_cur <- s_y + colSums(w) * R_SCALE / 2; m_cur <- m_y + m_h / 2
for (s in ord) for (m in mir) {
  v <- w[m, s]; if (v <= 0) next
  h <- v * R_SCALE; k <- k + 1
  dir <- if (delta[m, s] > 0) "up" else "down"
  fill <- if (m %in% VIRAL) paste0("viral_", dir) else paste0("human_", dir)
  polys[[k]] <- ribbon(X_S + SNODE_W / 2 + BK_W, X_M - NODE_W / 2, s_cur[s], s_cur[s] - h, m_cur[m], m_cur[m] - h, k)
  meta[[k]] <- tibble(id = k, side = "miRNA", source = s, target = m, value = v, signed_delta = delta[m, s], fill = fill)
  s_cur[s] <- s_cur[s] - h; m_cur[m] <- m_cur[m] - h
}
P <- bind_rows(polys) |> left_join(bind_rows(meta), by = "id")
links <- bind_rows(meta)

# ---- colours ---------------------------------------------------------------------------------------
# genus -> phylum from the MetaPhlAn clade paths; one hue family per phylum, shades within
tx <- read_tsv(tax_f, show_col_types = FALSE, col_select = 1)$clade
phy <- setNames(sub(".*\\|p__([^|]+).*", "\\1", tx), sub(".*\\|g__([^|]+).*", "\\1", tx))
phy <- phy[!duplicated(names(phy))]
fam <- list(Pseudomonadota = c("#08306B", "#9ECAE1"), Bacillota = c("#7F2704", "#FDAE6B"), Firmicutes = c("#7F2704", "#FDAE6B"),
            Actinomycetota = c("#00441B", "#A1D99B"), Bacteroidota = c("#3F007D", "#BCBDDC"),
            Mycoplasmatota = c("#67000D", "#FC9272"))
gp <- ifelse(is.na(phy[top]), "other", phy[top]); names(gp) <- top
gcol <- unlist(lapply(split(top, gp), function(gs) {
  base <- fam[[gp[gs[1]]]]; if (is.null(base)) base <- c("#252525", "#BDBDBD")
  setNames(if (length(gs) == 1) colorRampPalette(base)(3)[2] else colorRampPalette(base)(length(gs) + 1)[seq_along(gs)], gs)
}), use.names = FALSE) |> setNames(unlist(split(top, gp), use.names = FALSE))
gcol <- c(gcol[top], Other = "grey85")
phylum_note <- paste(sprintf("%s: %s", names(split(top, gp)), sapply(split(top, gp), paste, collapse = ", ")), collapse = "; ")
dcol <- c(human_up = "#D6604D", human_down = "#4393C3", viral_up = "#7B3294", viral_down = "#E7298A")
fills <- c(gcol, dcol)
grp_col <- c(control = "#4C78A8", DNAemia = "#E45756")

nodes <- bind_rows(
  tibble(x = X_G, y = g_y, h = g_h, fill = "grey35"),
  tibble(x = X_M, y = m_y, h = pmax(m_h, 0.004), fill = ifelse(mir %in% VIRAL, "#7B3294", "grey35")))
S_H <- 100 * L_SCALE
seg <- bind_rows(segs) |> mutate(xmin = X_S - SNODE_W / 2, xmax = X_S + SNODE_W / 2,
                                  fill = scales::alpha(gcol[genus], 0.55))   # same transparency as the ribbons
bk <- qc |> transmute(sample_id, bk = log10(mg_bk_per_million_qc + 1)) |> mutate(y = s_y[sample_id])
bk_pal <- scales::seq_gradient_pal("grey90", "black")
bk <- bk |> mutate(fill = bk_pal(scales::rescale(bk, from = range(bk))),
                   xmin = X_S + SNODE_W / 2, xmax = X_S + SNODE_W / 2 + BK_W, ymin = y - S_H / 2, ymax = y + S_H / 2)
lab_g <- tibble(x = X_G - NODE_W, y = g_y, lab = names(g_y), it = names(g_y) != "Other")
lab_m <- tibble(x = X_M + NODE_W, y = m_y, lab = ifelse(mir == "bkv-miR-B1-3p", "bkv-miR-B1-3p / jcv-miR-J1-3p", sub("^hsa-", "", mir)))

p <- ggplot() +
  geom_polygon(data = P, aes(x, y, group = id, fill = fill), alpha = 0.55, colour = NA) +
  geom_rect(data = nodes |> mutate(xmin = x - NODE_W / 2, xmax = x + NODE_W / 2, ymin = y - h / 2, ymax = y + h / 2),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = nodes$fill, colour = NA) +
  geom_text(data = lab_g, aes(x, y, label = lab, fontface = ifelse(it, "italic", "plain")), hjust = 1, size = 3) +
  geom_text(data = lab_m, aes(x, y, label = lab), hjust = 0, size = 3) +
  geom_rect(data = seg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = seg$fill, colour = NA) +
  geom_rect(data = bk, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = bk$fill, colour = NA) +
  geom_point(data = bk, aes(x = xmin, y = y, colour = bk), alpha = 0) +          # invisible; drives the colour bar
  geom_label(data = tibble(x = X_S - SNODE_W / 2 - 0.012, y = s_y, lab = ord), aes(x, y, label = lab),
             hjust = 1, size = 2.7, fontface = "bold", fill = scales::alpha("white", 0.85), label.size = 0,
             label.padding = unit(0.06, "lines")) +
  annotate("text", x = X_S + SNODE_W / 2 + BK_W / 2, y = max(s_y) + S_H / 2 + 0.015, label = "BK", size = 2.3) +
  annotate("text", x = c(X_G, X_S, X_M), y = 1.06, label = c("Genus", "Sample", "miRNA"), fontface = "bold", size = 3.6) +
  annotate("text", x = X_S, y = max(s_y[!is_case]) + 0.075, label = "Control",
           colour = "black", size = 3.2, fontface = "bold") +
  annotate("text", x = X_S, y = max(s_y[is_case]) + 0.055, label = "DNAemia",
           colour = "black", size = 3.2, fontface = "bold") +
  # scale bars
  annotate("rect", xmin = 0.25, xmax = 0.29, ymin = -0.06, ymax = -0.06 + 25 * L_SCALE, fill = "grey40") +
  annotate("text", x = 0.31, y = -0.06 + 25 * L_SCALE / 2, label = "= 25% of bacteria", hjust = 0, size = 2.7) +
  annotate("rect", xmin = 1.25, xmax = 1.29, ymin = -0.06, ymax = -0.06 + 2 * R_SCALE, fill = "grey40") +
  annotate("text", x = 1.31, y = -0.06 + 2 * R_SCALE / 2, label = "= 4-fold above/below control mean", hjust = 0, size = 2.7) +
  scale_fill_manual(values = fills, breaks = names(dcol),
                    labels = c("Human miRNA, above control mean", "Human miRNA, below control mean",
                               "Viral miRNA, above control mean", "Viral miRNA, below control mean"), name = NULL) +
  scale_colour_gradient(low = "grey90", high = "black", name = expression(atop("Urinary BKPyV DNA", "(log"[10]*" pairs per million)"))) +
  coord_cartesian(xlim = c(-0.45, 2.85), ylim = c(-0.08, 1.07), clip = "off") +
  theme_void(base_size = 10) + theme(legend.position = "none", plot.margin = margin(10, 10, 5, 10),
                                     plot.background = element_rect(fill = "white", colour = NA))

# ---- custom legend panel below the Sankey (one row) -------------------------------------------------
# 1) mini illustration: ribbons into a stacked bar + brace "Relative abundance of bacterial genera"
# 2) BKPyV DNA colour bar (same palette and range as the strips)   3) miRNA ribbon colours
ex_h <- c(0.14, 0.08, 0.18, 0.12, 0.28); ex_col <- unname(gcol[head(top, length(ex_h))])
ex_top <- 0.92 - cumsum(c(0, head(ex_h, -1) + 0.012)); ex_bar_top <- 0.92 - cumsum(c(0, head(ex_h, -1)))
ex <- bind_rows(lapply(seq_along(ex_h), function(i)
  ribbon(3.0, 3.9, ex_top[i], ex_top[i] - ex_h[i], ex_bar_top[i], ex_bar_top[i] - ex_h[i], i) |> mutate(col = ex_col[i])))
ex_bar <- tibble(ymax = ex_bar_top, ymin = ex_bar_top - ex_h, col = ex_col)
brace <- function(x, y0, y1, d = 0.07) {                      # curly brace, cusp pointing left
  r <- min(0.08, (y1 - y0) / 6); ym <- (y0 + y1) / 2
  arc <- function(cx, cy, a0, a1) { t <- seq(a0, a1, length.out = 15); tibble(x = cx + d * cos(t), y = cy + r * sin(t)) }
  bind_rows(arc(x + d, y1 - r, pi / 2, pi),          # top end -> vertical
            arc(x - d, ym + r, 0, -pi / 2),          # down into the cusp
            arc(x - d, ym - r, pi / 2, 0),           # out of the cusp
            arc(x + d, y0 + r, pi, 3 * pi / 2))      # vertical -> bottom end
}
br <- brace(2.85, min(ex_bar$ymin), max(ex_top))
cb_v <- seq(min(bk$bk), max(bk$bk), length.out = 60)
cb <- tibble(v = cb_v, ymin = 0.1 + 0.8 * (seq_along(cb_v) - 1) / 60, ymax = 0.1 + 0.8 * seq_along(cb_v) / 60,
             fill = bk_pal(scales::rescale(cb_v)))
cb_ticks <- tibble(v = intersect(2:6, floor(min(cb_v)):ceiling(max(cb_v)))) |> filter(v >= min(cb_v), v <= max(cb_v)) |>
  mutate(y = 0.1 + 0.8 * (v - min(cb_v)) / diff(range(cb_v)))
keys <- tibble(lab = c("Human miRNA, above control mean", "Human miRNA, below control mean",
                       "Viral miRNA, above control mean", "Viral miRNA, below control mean"),
               col = scales::alpha(unname(dcol), 0.8), y = c(0.84, 0.62, 0.40, 0.18))
leg <- ggplot() +
  geom_polygon(data = ex, aes(x, y, group = id), fill = scales::alpha(ex$col, 0.55), colour = NA) +
  geom_rect(data = ex_bar, aes(xmin = 3.9, xmax = 3.97, ymin = ymin, ymax = ymax), fill = scales::alpha(ex_bar$col, 0.55)) +
  geom_path(data = br, aes(x, y), linewidth = 0.5) +
  annotate("text", x = 2.68, y = 0.5, label = "Relative abundance of\nbacterial genera", hjust = 1, size = 3.2, lineheight = 0.95) +
  geom_rect(data = cb, aes(xmin = 5.0, xmax = 5.18, ymin = ymin, ymax = ymax), fill = cb$fill) +
  geom_segment(data = cb_ticks, aes(x = 5.18, xend = 5.23, y = y, yend = y), linewidth = 0.3) +
  geom_text(data = cb_ticks, aes(x = 5.27, y = y, label = v), hjust = 0, size = 3) +
  annotate("text", x = 5.09, y = 1.13, label = "Urinary BKPyV DNA", size = 3) +
  annotate("text", x = 5.09, y = 1.0, label = "(log[10]~pairs~per~million)", parse = TRUE, size = 3) +
  geom_rect(data = keys, aes(xmin = 6.1, xmax = 6.38, ymin = y - 0.085, ymax = y + 0.085), fill = keys$col) +
  geom_text(data = keys, aes(x = 6.48, y = y, label = lab), hjust = 0, size = 3.2) +
  coord_cartesian(xlim = c(0, 11), ylim = c(0, 1.2), expand = FALSE, clip = "off") +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
if (!is.na(legend_png)) {                                       # author-designed legend image replaces the drawn one
  img <- png::readPNG(legend_png)
  leg <- wrap_elements(full = grid::rasterGrob(img, interpolate = TRUE, height = unit(0.92, "npc")))   # aspect kept
}
p_full <- p / leg + plot_layout(heights = c(1, if (is.na(legend_png)) 0.22 else 0.26))
saveRDS(p_full, file.path(out, "Fig2_sankey.rds"))   # for analysis/07_main_figures.R
ggsave(file.path(out, "Fig2_sankey.pdf"), p_full, width = 12, height = 8, bg = "white")
ggsave(file.path(out, "Fig2_sankey.png"), p_full, width = 12, height = 8, dpi = 300, bg = "white")

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
  "| `sankey_mirnas.tsv` | miRNAs shown: the 8 human miRNAs with the smallest DESeq2 p in Fig 1A, plus bkv-miR-B1-5p and bkv-miR-B1-3p (= jcv-miR-J1-3p). |", "",
  "Left ribbon width: genus relative abundance within Bacteria in that sample (%), after the mean relative abundance > 0.1% filter; genera reaching >= 10% in at least one sample are shown (colours grouped by phylum), the rest pooled as Other. Each sample sums to 100%.",
  sprintf("Right ribbon width: |log2(DESeq2-normalised count + 1) - mean of the controls| for that sample; only values >= %.1f are drawn. Colour: above or below the control mean; viral miRNA in purple.", MIN_DELTA),
  "The two sides use independent scales (see scale bars). Each sample is drawn as a stacked bar of its genus composition; the grey strip next to it shows urinary BKPyV DNA from the metagenome (log10 read pairs per million; colour bar).",
  paste0("Genus colours by phylum (blue Pseudomonadota, orange Bacillota (Firmicutes), green Actinomycetota, purple Bacteroidota, red Mycoplasmatota): ", phylum_note, ".")
), file.path(out, "README.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("links:", nrow(links), "| genera shown:", length(top), "| miRNAs:", length(mir), "\n"); cat(phylum_note, "\n")
