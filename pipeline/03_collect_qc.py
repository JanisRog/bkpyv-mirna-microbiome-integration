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
cols = ["sample_id", "group",
        "mi_reads_raw", "mi_pct_with_adapter", "mi_reads_17_30nt", "mi_counts_hsa", "mi_pct_reads_hsa_mirna",
        "mi_n_hsa_ge10", "mi_bkv_5p", "mi_bkv_jcv_3p_shared", "mi_jcv_5p", "mi_bkv_5p_cpm",
        "ercc_reads_rrna", "ercc_reads_genome", "ercc_reads_annotated", "ercc_annotated_fraction",
        "ercc_c1_pass", "ercc_c2_pass", "ercc_pass",
        "mt_pct_mirna", "mt_pct_rrna", "mt_pct_trna", "mt_pct_artifact", "mt_pct_unknown",
        "mg_pairs_raw", "mg_pairs_after_fastp", "mg_pairs_nonhuman", "mg_pct_nonhuman",
        "mg_pairs_bk", "mg_pairs_jc", "mg_bk_per_million_qc", "mg_pairs_clean",
        "mg_est_bacterial_reads", "mg_n_species"]
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
        er.get("reads_rrna", ""), er.get("reads_genome_mapped", ""), er.get("reads_annotated_rna", ""),
        er.get("annotated_fraction", ""), er.get("ercc_c1_pass", ""), er.get("ercc_c2_pass", ""), er.get("ercc_pass", ""),
        rt_pct("mirna"), rt_pct("rrna"), rt_pct("trna"), rt_pct("artifact"), rt_pct("unknown"),
        mg.get("pairs_raw", ""), mg.get("pairs_after_fastp", ""), mg.get("pairs_nonhuman", ""),
        pct(mg.get("pairs_nonhuman"), mg.get("pairs_after_fastp")),
        mg.get("pairs_bk", ""), mg.get("pairs_jc", ""),
        round(1e6 * mg["pairs_bk"] / mg["pairs_after_fastp"], 2) if mg.get("pairs_after_fastp") else "",
        mg.get("pairs_clean", ""), mg.get("metaphlan_est_bacterial_reads", ""), mg.get("metaphlan_n_species", ""),
    ])
write_tsv(f"{QC}/qc_combined.tsv", cols, rows)
missing = [s for s in samples if s not in counts or s not in rel]
if missing:
    print(f"[WARN] incomplete samples: {', '.join(missing)}")
