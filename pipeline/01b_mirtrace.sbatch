#!/usr/bin/env bash
#SBATCH --job-name=bk_v2_mirtrace
#SBATCH --partition=lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/shares/amr.imm.uzh/data/projects/jaroge/urine_microbiome_bk/v2/logs/%x_%j.out
#SBATCH --error=/shares/amr.imm.uzh/data/projects/jaroge/urine_microbiome_bk/v2/logs/%x_%j.err

# miRTrace QC on raw R1 of all samples (library composition: miRNA / rRNA / tRNA / artifacts).
set -euo pipefail
source "${SCRIPTS:-/shares/amr.imm.uzh/data/projects/jaroge/urine_microbiome_bk/v2/scripts}/config.sh"
source "$CONTAINERS_ENV"
load_apptainer
OUT=$OUT_MIRNA/mirtrace
mkdir -p "$OUT"
# config: path,name,adapter
awk -F'\t' -v a="$ADAPTER_3P" 'NR>1{print $4","$1","a}' "$MANIFEST" > "$OUT/mirtrace_config.csv"
export JAVA_TOOL_OPTIONS="-Xmx20g"
run mirtrace mirtrace qc --species hsa --config "$OUT/mirtrace_config.csv" \
  --output-dir "$OUT/report" --num-threads "${SLURM_CPUS_PER_TASK:-8}" --force
echo "[DONE] mirtrace -> $OUT/report"
