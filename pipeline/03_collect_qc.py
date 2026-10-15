#!/usr/bin/env python3
"""Merge per-sample outputs into count matrices and one QC table (stdlib only).

Run on the login node after 01/01b/02 finish:
    python3 $SCRIPTS/03_collect_qc.py
Writes to $V2/qc/. No thresholds are applied here; that decision is made on qc_combined.tsv.
"""
import csv, glob, json, os

V2 = os.environ.get("V2", "/shares/amr.imm.uzh/data/projects/jaroge/urine_microbiome_bk/v2")
QC = f"{V2}/qc"
os.makedirs(QC, exist_ok=True)

manifest = list(csv.DictReader(open(f"{V2}/sample_manifest.tsv"), delimiter="\t"))
samples = [r["sample_id"] for r in manifest]
group = {r["sample_id"]: r["group"] for r in manifest}


def load_json(path):
    return json.load(open(path)) if os.path.exists(path) else None


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
    print(f"[INFO] wrote {path} ({len(rows)} rows)")


# ---- miRNA count matrix (hsa + bkv, raw counts) ----
counts = {}
for s in samples:
    f = f"{V2}/mirna/{s}/{s}.counts.tsv"
    if os.path.exists(f):
        counts[s] = {l.split("\t")[0]: int(float(l.split("\t")[1])) for l in open(f)}
mirnas = sorted({m for c in counts.values() for m in c}, key=lambda m: (m.startswith("hsa-"), m))
done = [s for s in samples if s in counts]
write_tsv(f"{QC}/mirna_counts_raw.tsv", ["miRNA"] + done,
          [[m] + [counts[s].get(m, 0) for s in done] for m in mirnas])

# ---- miRTrace RNA-type composition (optional) ----
# mirtrace-stats-rnatype.tsv: rows = RNA type (miRNA, rRNA, tRNA, artifact, unknown), columns = samples
rnatype = {}
rt_files = glob.glob(f"{V2}/mirna/mirtrace/report/mirtrace-stats-rnatype.tsv")
if rt_files:
    rt_rows = list(csv.reader(open(rt_files[0]), delimiter="\t"))
    rt_samples = rt_rows[0][1:]
    for r in rt_rows[1:]:
        for smp, v in zip(rt_samples, r[1:]):
            rnatype.setdefault(smp, {})[r[0].lower()] = float(v)

# ---- MetaPhlAn species matrices ----
rel, est = {}, {}
for s in samples:
    f = f"{V2}/metagenome/{s}/{s}.metaphlan.txt"
    if not os.path.exists(f):
        continue
    rel[s], est[s] = {}, {}
    for l in open(f):
        if l.startswith("#"):
            continue
        x = l.rstrip("\n").split("\t")
        if x[0].split("|")[-1].startswith("s__"):
            rel[s][x[0]] = float(x[2])
            est[s][x[0]] = int(float(x[4])) if len(x) > 4 else 0
species = sorted({k for d in rel.values() for k in d})
mg_done = [s for s in samples if s in rel]
write_tsv(f"{QC}/metaphlan_species_relab.tsv", ["clade"] + mg_done,
          [[k] + [rel[s].get(k, 0.0) for s in mg_done] for k in species])
write_tsv(f"{QC}/metaphlan_species_est_reads.tsv", ["clade"] + mg_done,
          [[k] + [est[s].get(k, 0) for s in mg_done] for k in species])

# ---- combined QC table ----
# Column definitions: the single source for the table header and for qc/README.md.
COLUMNS = {
    "sample_id": "Sample ID (S<number>; samples 6, 9, 10 do not exist).",
    "group": "DNAemia = BKPyV-DNAemia >4 log10 copies/mL with immunosuppression reduction ('treatment' in the metadata); control = no BKPyV-DNAemia.",
    # small RNA (01_mirna)
    "mi_reads_raw": "Small RNA: raw reads in R1 (50 bp).",
    "mi_pct_with_adapter": "Small RNA: % of raw reads containing the TruSeq small-RNA 3' adapter (TGGAATTCTCGGGTGCCAAGG).",
    "mi_reads_17_30nt": "Small RNA: reads kept after trimming (adapter found, 17-30 nt, Q20). Most lost reads are adapter dimers (insert <17 nt).",
    "mi_counts_hsa": "Small RNA: sum of human mature miRNA counts (miRDeep2 quantifier, miRBase 23.1). Used as library size for viral CPM.",
    "mi_pct_reads_hsa_mirna": "Small RNA: mi_counts_hsa as % of mi_reads_17_30nt.",
    "mi_n_hsa_ge10": "Small RNA: number of human miRNAs with >=10 counts (library complexity; descriptive only).",
    "mi_bkv_5p": "Small RNA: bkv-miR-B1-5p counts. BKPyV-specific (differs from jcv-miR-J1-5p at 3/22 nt).",
    "mi_bkv_jcv_3p_shared": "Small RNA: bkv-miR-B1-3p counts. Sequence IDENTICAL to jcv-miR-J1-3p, so it cannot distinguish BKPyV from JCPyV.",
    "mi_jcv_5p": "Small RNA: jcv-miR-J1-5p counts (JCPyV-specific).",
    "mi_bkv_5p_cpm": "Small RNA: bkv-miR-B1-5p per million human miRNA counts.",
    # ERCC (01c_ercc_qc)
    "ercc_reads_rrna": "ERCC: trimmed reads aligning to human rRNA (45S NR_046235.3, 5S NR_023363.1); removed before genome alignment.",
    "ercc_reads_genome": "ERCC: remaining reads aligned to GRCh38 (bowtie2, primary alignments).",
    "ercc_reads_pirna": "ERCC: genome-aligned reads outside GENCODE/tRNA annotation that match a piRBase v3.0 human gold piRNA.",
    "ercc_reads_annotated": "ERCC: reads on annotated human RNA = GENCODE v50 genes (miRNA, snoRNA, snRNA, lncRNA, mRNA, ...) + GENCODE tRNAs + piRNA. Ribosomal features excluded.",
    "ercc_annotated_fraction_without_pirna": "ERCC: annotated fraction without piRNA (GENCODE + tRNA only); kept to show the effect of adding piRNA.",
    "ercc_annotated_fraction": "ERCC: ercc_reads_annotated / (ercc_reads_genome - reads on ribosomal features).",
    "ercc_c1_pass": "ERCC criterion 1 (ERCC small RNA-seq QC standard v1.0, exrna.org; Murillo et al. Cell 2019): >=100,000 reads on annotated human RNA. Used as the sample exclusion rule for differential expression.",
    "ercc_c2_pass": "ERCC criterion 2: annotated fraction > 0.5 (i.e. most genome-mapped, non-ribosomal reads fall in annotated RNA). Reported per sample.",
    "ercc_pass": "Both ERCC criteria met.",
    # miRTrace (01b)
    "mt_pct_mirna": "miRTrace: % of miRTrace QC-passed reads classified as miRNA.",
    "mt_pct_rrna": "miRTrace: % rRNA.",
    "mt_pct_trna": "miRTrace: % tRNA.",
    "mt_pct_artifact": "miRTrace: % artifacts (mostly adapter dimers / low complexity).",
    "mt_pct_unknown": "miRTrace: % unclassified (not miRNA/rRNA/tRNA/artifact; includes snoRNA, Y RNA, mRNA fragments, piRNA).",
    # metagenome (02_metagenome)
    "mg_pairs_raw": "Metagenome: raw read pairs.",
    "mg_pairs_after_fastp": "Metagenome: read pairs after fastp (Q20, >=50 bp).",
    "mg_pairs_nonhuman": "Metagenome: read pairs not concordantly aligned to GRCh38 (bowtie2 --very-sensitive).",
    "mg_pct_nonhuman": "Metagenome: mg_pairs_nonhuman as % of mg_pairs_after_fastp (= % non-human DNA).",
    "mg_pairs_bk": "Metagenome: non-human pairs with >=1 mate aligned to BKPyV (NC_001538.1).",
    "mg_pairs_jc": "Metagenome: non-human pairs with >=1 mate aligned to JCPyV (NC_001699.1).",
    "mg_bk_per_million_qc": "Metagenome: mg_pairs_bk per million mg_pairs_after_fastp (urinary BKPyV DNA load proxy).",
    "mg_pairs_clean": "Metagenome: pairs after removing human and polyomavirus reads; input to MetaPhlAn.",
    "mg_est_bacterial_reads": "Metagenome: MetaPhlAn 4.2.6 estimated reads assigned to Bacteria (rel_ab_w_read_stats).",
    "mg_n_species": "Metagenome: number of species detected by MetaPhlAn.",
}
cols = list(COLUMNS)
rows = []
for s in samples:
    mi = load_json(f"{V2}/mirna/{s}/{s}.run_manifest.json") or {}
    mg = load_json(f"{V2}/metagenome/{s}/{s}.run_manifest.json") or {}
    er = load_json(f"{V2}/mirna/{s}/{s}.ercc.json") or {}
    rt = rnatype.get(s, {})

    def pct(a, b):
        return round(100 * a / b, 3) if a is not None and b else ""

    def rt_pct(key):
        tot = sum(rt.values())
        return round(100 * rt[key] / tot, 2) if key in rt and tot else ""

    lib = mi.get("mirna_counts_sum_hsa") or 0   # CPM relative to human miRNA counts
    rows.append([
        s, group[s],
        mi.get("reads_raw", ""), pct(mi.get("reads_with_adapter"), mi.get("reads_raw")),
        mi.get("reads_17_30nt", ""), mi.get("mirna_counts_sum_hsa", ""),
        pct(mi.get("mirna_counts_sum_hsa"), mi.get("reads_17_30nt")),
        mi.get("n_hsa_mirnas_ge10", ""), mi.get("bkv_miR_B1_5p", ""), mi.get("bkv_jcv_miR_3p", ""),
        mi.get("jcv_miR_J1_5p", ""), round(1e6 * mi["bkv_miR_B1_5p"] / lib, 1) if lib else "",
        er.get("reads_rrna", ""), er.get("reads_genome_mapped", ""), er.get("reads_pirna", ""), er.get("reads_annotated_rna", ""),
        er.get("annotated_fraction_without_pirna", ""), er.get("annotated_fraction", ""), er.get("ercc_c1_pass", ""), er.get("ercc_c2_pass", ""), er.get("ercc_pass", ""),
        rt_pct("mirna"), rt_pct("rrna"), rt_pct("trna"), rt_pct("artifact"), rt_pct("unknown"),
        mg.get("pairs_raw", ""), mg.get("pairs_after_fastp", ""), mg.get("pairs_nonhuman", ""),
        pct(mg.get("pairs_nonhuman"), mg.get("pairs_after_fastp")),
        mg.get("pairs_bk", ""), mg.get("pairs_jc", ""),
        round(1e6 * mg["pairs_bk"] / mg["pairs_after_fastp"], 2) if mg.get("pairs_after_fastp") else "",
        mg.get("pairs_clean", ""), mg.get("metaphlan_est_bacterial_reads", ""), mg.get("metaphlan_n_species", ""),
    ])
write_tsv(f"{QC}/qc_combined.tsv", cols, rows)

# ---- README for the qc folder (generated, so it always matches the tables) ----
with open(f"{QC}/README.md", "w") as fh:
    fh.write("# QC tables (generated by pipeline/03_collect_qc.py; do not edit by hand)\n\n")
    fh.write("## Files\n\n| File | Contents |\n|---|---|\n")
    fh.write("| `qc_combined.tsv` | One row per sample; all QC metrics (columns below). Empty cells = step not run/finished. |\n")
    fh.write("| `mirna_counts_raw.tsv` | Raw miRNA counts, miRNA x sample (human + BKPyV + JCPyV; miRBase 23.1). Mature miRNAs on several identical precursors: max over precursors. |\n")
    fh.write("| `metaphlan_species_relab.tsv` | MetaPhlAn species relative abundance (%), clade x sample. |\n")
    fh.write("| `metaphlan_species_est_reads.tsv` | MetaPhlAn estimated reads per species, clade x sample. |\n\n")
    fh.write("## Columns of qc_combined.tsv\n\nPrefixes: `mi_` small RNA, `ercc_` ERCC QC, `mt_` miRTrace, `mg_` metagenome.\n\n")
    fh.write("| Column | Meaning |\n|---|---|\n")
    for c, d in COLUMNS.items():
        fh.write(f"| `{c}` | {d} |\n")
print(f"[INFO] wrote {QC}/README.md")
missing = [s for s in samples if s not in counts or s not in rel]
if missing:
    print(f"[WARN] incomplete samples: {', '.join(missing)}")
