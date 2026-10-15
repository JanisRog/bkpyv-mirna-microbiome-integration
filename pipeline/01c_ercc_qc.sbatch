#!/usr/bin/env bash
#SBATCH --job-name=bk_v2_ercc
#SBATCH --partition=lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=10:00:00
#SBATCH --output=/shares/amr.imm.uzh/data/projects/jaroge/urine_microbiome_bk/v2/logs/%x_%j.out
#SBATCH --error=/shares/amr.imm.uzh/data/projects/jaroge/urine_microbiome_bk/v2/logs/%x_%j.err

# ERCC small-RNA QC standard v1.0 on the trimmed 17-30 nt reads from 01_mirna:
#   rRNA removal (45S + 5S) -> GRCh38 alignment -> overlap with annotated RNA
#   (GENCODE all genes incl. miRNA/snoRNA/snRNA/lncRNA/mRNA + GENCODE tRNA set).
#   criterion 1: annotated reads >= ERCC_MIN_ANNOTATED
#   criterion 2: annotated / genome-mapped non-ribosomal reads > ERCC_MIN_FRACTION
# Deviation from exceRpt: piRNA and circRNA annotations are not included (reported in Methods).
# Small data, so all samples run sequentially in one job; finished samples are skipped.
set -euo pipefail
source "${SCRIPTS:-/shares/amr.imm.uzh/data/projects/jaroge/urine_microbiome_bk/v2/scripts}/config.sh"
source "$CONTAINERS_ENV"
load_apptainer
T=${SLURM_CPUS_PER_TASK:-8}

# ---- reference (built once) -------------------------------------------------
ER=$REF/ercc_gencode$GENCODE_VERSION
mkdir -p "$ER"
if [[ ! -s "$ER/rrna.1.bt2" ]]; then
  curl -fsSL --retry 5 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${RRNA_ACCS// /,}&rettype=fasta" -o "$ER/rrna.fa"
  [[ $(grep -c '^>' "$ER/rrna.fa") -eq $(wc -w <<< "$RRNA_ACCS") ]] || { echo "[ERROR] rRNA download incomplete"; exit 1; }
  run bowtie2 bowtie2-build -q "$ER/rrna.fa" "$ER/rrna"
fi
if [[ ! -s "$ER/annotated.bed" ]]; then
  for f in primary_assembly.annotation tRNAs; do
    [[ -s "$ER/$f.gtf.gz" ]] || curl -fsSL --retry 5 "$GENCODE_URL/gencode.v$GENCODE_VERSION.$f.gtf.gz" -o "$ER/$f.gtf.gz"
  done
  # ribosomal features (excluded from both numerator and denominator) vs. all other annotated RNA
  zcat "$ER/primary_assembly.annotation.gtf.gz" | awk -F'\t' -v OFS='\t' '
    $3=="exon"{match($9,/gene_type "[^"]+"/); gt=substr($9,RSTART+11,RLENGTH-12);
      out=(gt=="rRNA"||gt=="rRNA_pseudogene"||gt=="Mt_rRNA") ? "ribo" : "annot"; print $1,$4-1,$5 > ("'"$ER"'/" out ".unsorted.bed")}'
  zcat "$ER/tRNAs.gtf.gz" | awk -F'\t' -v OFS='\t' '!/^#/ && $3=="tRNA"{print $1,$4-1,$5}' >> "$ER/annot.unsorted.bed"
  for b in ribo annot; do
    sort -k1,1 -k2,2n "$ER/$b.unsorted.bed" > "$ER/$b.sorted.bed"
    run bedtools bedtools merge -i "$ER/$b.sorted.bed" > "$ER/$b.bed"
    rm -f "$ER/$b.unsorted.bed" "$ER/$b.sorted.bed"
  done
  mv "$ER/annot.bed" "$ER/annotated.bed"
fi
# chromosome naming must match the bowtie2 index
run bowtie2 bowtie2-inspect -n "$HOST_INDEX" | awk '{print $1}' | sort -u > "$ER/index_chroms.txt"
shared=$(cut -f1 "$ER/annotated.bed" | sort -u | comm -12 - "$ER/index_chroms.txt")
grep -qx chr1 <<< "$shared" || { echo "[ERROR] chromosome names of GENCODE and $HOST_INDEX do not match"; exit 1; }

# ---- per sample ---------------------------------------------------------------
tail -n +2 "$MANIFEST" | cut -f1 | while read -r SID; do
  OUT=$OUT_MIRNA/$SID
  IN=$OUT/$SID.trim.fastq.gz
  [[ -s "$OUT/$SID.ercc.json" ]] && { echo "[SKIP] $SID"; continue; }
  [[ -s "$IN" ]] || { echo "[WARN] $SID: no trimmed reads, run 01_mirna first"; continue; }
  W=$TMP_ROOT/ercc_$SID; mkdir -p "$W"

  run bowtie2 bowtie2 -p "$T" $SMALLRNA_BT2_OPTS -x "$ER/rrna" -U "$IN" --un-gz "$W/nonrrna.fq.gz" -S /dev/null 2> "$OUT/$SID.ercc_rrna.log"
  run bowtie2 bowtie2 -p "$T" $SMALLRNA_BT2_OPTS -x "$HOST_INDEX" -U "$W/nonrrna.fq.gz" 2> "$OUT/$SID.ercc_genome.log" \
    | run samtools samtools view -b -F 260 -o "$W/genome.bam" -
  count_reads() { run samtools samtools view -c "$1"; }
  n_in=$(( $(zcat "$IN" | wc -l) / 4 ))
  n_nonrrna=$(( $(zcat "$W/nonrrna.fq.gz" | wc -l) / 4 ))
  n_genome=$(count_reads "$W/genome.bam")
  run bedtools bedtools intersect -u -a "$W/genome.bam" -b "$ER/ribo.bed" > "$W/ribo.bam"
  run bedtools bedtools intersect -v -a "$W/genome.bam" -b "$ER/ribo.bed" > "$W/nonribo.bam"
  run bedtools bedtools intersect -u -a "$W/nonribo.bam" -b "$ER/annotated.bed" > "$W/annot.bam"
  n_mtribo=$(count_reads "$W/ribo.bam")
  n_annot=$(count_reads "$W/annot.bam")

  python3 - "$SID" "$OUT" "$n_in" "$n_nonrrna" "$n_genome" "$n_mtribo" "$n_annot" <<PY
import json, sys
sid, out = sys.argv[1], sys.argv[2]
n_in, n_nonrrna, n_genome, n_mtribo, n_annot = map(int, sys.argv[3:])
denom = n_genome - n_mtribo
frac = n_annot / denom if denom else 0.0
m = {
  "sample_id": sid, "standard": "ERCC small RNA-seq QC v1.0 (exrna.org)",
  "gencode_version": $GENCODE_VERSION, "rrna_accessions": "$RRNA_ACCS", "host_index": "$HOST_INDEX",
  "bowtie2_opts": "$SMALLRNA_BT2_OPTS", "annotation_excludes": "piRNA, circRNA",
  "reads_17_30nt": n_in, "reads_rrna": n_in - n_nonrrna,
  "reads_genome_mapped": n_genome, "reads_mt_or_annotated_ribosomal": n_mtribo,
  "reads_annotated_rna": n_annot, "annotated_fraction": round(frac, 4),
  "ercc_c1_pass": n_annot >= $ERCC_MIN_ANNOTATED, "ercc_c2_pass": frac > $ERCC_MIN_FRACTION,
}
m["ercc_pass"] = m["ercc_c1_pass"] and m["ercc_c2_pass"]
json.dump(m, open(f"{out}/{sid}.ercc.json", "w"), indent=1)
print(json.dumps(m))
PY
  rm -rf "$W"
done
echo "[DONE] ERCC QC"
