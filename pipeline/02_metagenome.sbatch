#!/usr/bin/env bash
#SBATCH --job-name=bk_v2_mg
#SBATCH --partition=lowprio
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --array=1-22%6
#SBATCH --output=/shares/amr.imm.uzh/data/projects/jaroge/urine_microbiome_bk/v2/logs/%x_%A_%a.out
#SBATCH --error=/shares/amr.imm.uzh/data/projects/jaroge/urine_microbiome_bk/v2/logs/%x_%A_%a.err

# Per sample: fastp -> host removal (T2T-CHM13) -> BK/JC read counts + removal -> MetaPhlAn4.
# Steps whose output exists are skipped, so a timed-out task can simply be resubmitted.
set -euo pipefail
source "${SCRIPTS:-/shares/amr.imm.uzh/data/projects/jaroge/urine_microbiome_bk/v2/scripts}/config.sh"
source "$CONTAINERS_ENV"
load_apptainer
T=${SLURM_CPUS_PER_TASK:-16}

SID=$(awk -F'\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR==i+1{print $1}' "$MANIFEST")
R1=$(manifest_get "$SID" mg_r1); R2=$(manifest_get "$SID" mg_r2)
OUT=$OUT_MG/$SID; TMP=$TMP_ROOT/mg_$SID
mkdir -p "$OUT" "$TMP"
echo "[INFO] $SID  R1=$R1"

# 1. read QC (trimmed reads only kept in tmp; ~10 GB per sample)
if [[ ! -s "$OUT/$SID.fastp.json" ]]; then
  run fastp fastp -i "$R1" -I "$R2" -o "$TMP/qc_1.fq.gz" -O "$TMP/qc_2.fq.gz" \
    --detect_adapter_for_pe -q "$MG_QUAL" -l "$MG_MIN_LEN" -w "$T" \
    -j "$OUT/$SID.fastp.json" -h "$OUT/$SID.fastp.html" 2> "$OUT/$SID.fastp.log"
fi

# 2. host removal: keep pairs not concordantly aligned to human
if [[ ! -s "$OUT/$SID.nonhuman_1.fq.gz" ]]; then
  run bowtie2 bowtie2 -p "$T" --very-sensitive -x "$HOST_INDEX" -1 "$TMP/qc_1.fq.gz" -2 "$TMP/qc_2.fq.gz" \
    --un-conc-gz "$TMP/nonhuman_%.fq.gz" -S /dev/null 2> "$OUT/$SID.bowtie2_host.log"
  mv "$TMP/nonhuman_1.fq.gz" "$OUT/$SID.nonhuman_1.fq.gz"
  mv "$TMP/nonhuman_2.fq.gz" "$OUT/$SID.nonhuman_2.fq.gz"
  rm -f "$TMP/qc_1.fq.gz" "$TMP/qc_2.fq.gz"
fi

# 3. polyomavirus: count read pairs with >=1 primary mapped mate per genome, remove them
if [[ ! -s "$OUT/$SID.virus_counts.tsv" ]]; then
  run bowtie2 bowtie2 -p "$T" --very-sensitive --no-unal -x "$REF/polyomavirus/bk_jc" \
    -1 "$OUT/$SID.nonhuman_1.fq.gz" -2 "$OUT/$SID.nonhuman_2.fq.gz" \
    --un-conc-gz "$TMP/clean_%.fq.gz" 2> "$OUT/$SID.bowtie2_virus.log" \
  | awk -F'\t' '!/^@/ && int($2/4)%2==0 && int($2/256)%2==0 && int($2/2048)%2==0 {seen[$3"\t"$1]=1}
                END{for(k in seen){split(k,a,"\t"); n[a[1]]++} for(r in n) print r"\t"n[r]}' \
  > "$TMP/virus_counts.tsv"
  mv "$TMP/clean_1.fq.gz" "$OUT/$SID.clean_1.fq.gz"
  mv "$TMP/clean_2.fq.gz" "$OUT/$SID.clean_2.fq.gz"
  mv "$TMP/virus_counts.tsv" "$OUT/$SID.virus_counts.tsv"
  [[ -s "$OUT/$SID.virus_counts.tsv" ]] || printf "none\t0\n" > "$OUT/$SID.virus_counts.tsv"
fi

# 4. bacterial profile
if [[ ! -s "$OUT/$SID.metaphlan.txt" ]]; then
  run metaphlan metaphlan -1 "$OUT/$SID.clean_1.fq.gz" -2 "$OUT/$SID.clean_2.fq.gz" --input_type fastq \
    --db_dir "$MPA_DB" -x "$MPA_INDEX" --offline --nproc "$T" -t rel_ab_w_read_stats \
    --mapout "$OUT/$SID.metaphlan_map.bz2" --tmp_dir "$TMP" --sample_id "$SID" \
    -o "$TMP/metaphlan.txt" 2> "$OUT/$SID.metaphlan.log"
  mv "$TMP/metaphlan.txt" "$OUT/$SID.metaphlan.txt"
fi

# 5. per-sample run manifest
count_pairs() { echo $(( $(zcat "$1" | wc -l) / 4 )); }
NONHUMAN=$(count_pairs "$OUT/$SID.nonhuman_1.fq.gz")
CLEAN=$(count_pairs "$OUT/$SID.clean_1.fq.gz")
python3 - "$SID" "$OUT" "$NONHUMAN" "$CLEAN" <<PY
import json, sys
sid, out, nonhuman, clean = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
fj = json.load(open(f"{out}/{sid}.fastp.json"))["summary"]
virus = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(f"{out}/{sid}.virus_counts.tsv")}
mpa_reads, n_species, bact_reads = None, 0, 0
for l in open(f"{out}/{sid}.metaphlan.txt"):
    if l.startswith("#") and "reads processed" in l: mpa_reads = int(l[1:].split()[0])
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t")
    if f[0].split("|")[-1].startswith("s__"): n_species += 1
    if f[0] == "k__Bacteria" and len(f) > 4: bact_reads = int(float(f[4]))
m = {
  "sample_id": sid, "input_r1": "$R1", "input_r2": "$R2",
  "slurm_job": "${SLURM_ARRAY_JOB_ID:-}_${SLURM_ARRAY_TASK_ID:-}",
  "container_fastp": "$SIF_FASTP", "container_bowtie2": "$SIF_BOWTIE2", "container_metaphlan": "$SIF_METAPHLAN",
  "host_index": "$HOST_INDEX", "virus_accessions": "$VIRUS_ACCS", "metaphlan_index": "$MPA_INDEX",
  "fastp_qual": $MG_QUAL, "fastp_min_len": $MG_MIN_LEN, "bowtie2_mode": "--very-sensitive",
  "pairs_raw": fj["before_filtering"]["total_reads"] // 2,
  "pairs_after_fastp": fj["after_filtering"]["total_reads"] // 2,
  "pairs_nonhuman": nonhuman,
  "pairs_bk": virus.get("NC_001538.1", 0), "pairs_jc": virus.get("NC_001699.1", 0),
  "pairs_clean": clean,
  "metaphlan_reads_processed": mpa_reads,
  "metaphlan_est_bacterial_reads": bact_reads,
  "metaphlan_n_species": n_species,
}
json.dump(m, open(f"{out}/{sid}.run_manifest.json", "w"), indent=1)
print(json.dumps(m, indent=1))
PY
rm -rf "$TMP"
echo "[DONE] $SID"
