#!/usr/bin/env bash
#SBATCH --job-name=bk_v2_mirna
#SBATCH --partition=lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-22
#SBATCH --output=/shares/amr.imm.uzh/data/projects/jaroge/urine_microbiome_bk/v2/logs/%x_%A_%a.out
#SBATCH --error=/shares/amr.imm.uzh/data/projects/jaroge/urine_microbiome_bk/v2/logs/%x_%A_%a.err

# Per sample: cutadapt (R1 only) -> collapse -> miRDeep2 quantifier on ONE hsa+bkv+jcv reference.
set -euo pipefail
source "${SCRIPTS:-/shares/amr.imm.uzh/data/projects/jaroge/urine_microbiome_bk/v2/scripts}/config.sh"
source "$CONTAINERS_ENV"
load_apptainer
THREADS=${SLURM_CPUS_PER_TASK:-4}

SID=$(awk -F'\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR==i+1{print $1}' "$MANIFEST")
R1=$(manifest_get "$SID" mirna_r1)
MIRREF=$REF/mirbase_$MIRBASE_VERSION
OUT=$OUT_MIRNA/$SID
mkdir -p "$OUT"
echo "[INFO] $SID  R1=$R1"

# 1. trimming: reads without adapter cannot be 17-30 nt small RNAs at 150 bp -> discarded
run cutadapt cutadapt -j "$THREADS" -a "$ADAPTER_3P" --discard-untrimmed \
  -q "$MIRNA_QUAL" -m "$MIRNA_MIN_LEN" -M "$MIRNA_MAX_LEN" \
  --json "$OUT/$SID.cutadapt.json" -o "$OUT/$SID.trim.fastq.gz" "$R1" > "$OUT/$SID.cutadapt.log"

# 2. collapse to miRDeep2 read format >seq_<i>_x<count>
zcat "$OUT/$SID.trim.fastq.gz" | awk 'NR%4==2{c[$0]++} END{i=0; for(s in c){i++; print ">seq_"i"_x"c[s]; print s}}' > "$OUT/$SID.collapsed.fa"

# 3. quantify (no -t species filter: the reference is already hsa+bkv+jcv only)
WORK=$OUT/quantifier
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
run mirdeep2 quantifier.pl -p "$MIRREF/hairpin.hsa_bkv.fa" -m "$MIRREF/mature.hsa_bkv.fa" \
  -r "$OUT/$SID.collapsed.fa" -y "$SID" -d -j > "$OUT/$SID.quantifier.log" 2>&1
# mature miRNAs on several identical precursors get the same reads on each -> take max, not sum
awk -F'\t' 'NR>1{m=$1; v=$2+0; if(!(m in c) || v>c[m]) c[m]=v} END{for(m in c) print m"\t"c[m]}' \
  "$WORK/miRNAs_expressed_all_samples_$SID.csv" | sort > "$OUT/$SID.counts.tsv"
cd "$OUT"; rm -rf "$WORK" "$OUT/$SID.collapsed.fa"

# 4. per-sample run manifest
python3 - "$SID" "$R1" "$OUT" <<PY
import json, sys
sid, r1, out = sys.argv[1:]
cj = json.load(open(f"{out}/{sid}.cutadapt.json"))
counts = {l.split("\t")[0]: int(float(l.split("\t")[1])) for l in open(f"{out}/{sid}.counts.tsv")}
m = {
  "sample_id": sid, "input_r1": r1, "slurm_job": "${SLURM_ARRAY_JOB_ID:-}_${SLURM_ARRAY_TASK_ID:-}",
  "container_cutadapt": "$SIF_CUTADAPT", "container_mirdeep2": "$SIF_MIRDEEP2",
  "mirbase_version": "$MIRBASE_VERSION", "reference_species": "$MIRNA_SPECIES_PREFIXES",
  "adapter_3p": "$ADAPTER_3P", "min_len": $MIRNA_MIN_LEN, "max_len": $MIRNA_MAX_LEN, "qual": $MIRNA_QUAL,
  "discard_untrimmed": True, "multi_precursor_rule": "max",
  "reads_raw": cj["read_counts"]["input"],
  "reads_with_adapter": cj["read_counts"]["read1_with_adapter"],
  "reads_17_30nt": cj["read_counts"]["output"],
  "mirna_counts_sum_hsa": sum(v for k, v in counts.items() if k.startswith("hsa-")),
  "bkv_miR_B1_5p": counts.get("bkv-miR-B1-5p", 0),
  "bkv_jcv_miR_3p": counts.get("bkv-miR-B1-3p", 0),   # identical to jcv-miR-J1-3p, not virus-specific
  "jcv_miR_J1_5p": counts.get("jcv-miR-J1-5p", 0),
  "n_hsa_mirnas_ge10": sum(1 for k, v in counts.items() if k.startswith("hsa-") and v >= 10),
}
json.dump(m, open(f"{out}/{sid}.run_manifest.json", "w"), indent=1)
print(json.dumps(m, indent=1))
PY
echo "[DONE] $SID"
