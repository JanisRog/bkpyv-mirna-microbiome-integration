# Urinary miRNA and microbiome in kidney transplant recipients with BKPyV-DNAemia

Code for an independent reprocessing and integration of urinary small-RNA sequencing and shotgun
metagenomics from 22 kidney transplant recipients (11 with BK polyomavirus DNAemia, 11 controls).
It covers differential expression of human and BKPyV miRNAs, the urinary bacterial microbiome, and
microbiome–miRNA associations.

> Manuscript in preparation (short communication, *Transplant Infectious Disease*). Citation will be added on publication.

## Repository layout

| Path | Contents |
|---|---|
| `pipeline/` | Read processing on a SLURM cluster with Apptainer containers (QC, quantification, profiling) |
| `analysis/` | Statistics and figures from the processed tables |

## Pipeline

All tools run from pinned BioContainers images, listed in `pipeline/config.sh`.

| Step | Script | What it does |
|---|---|---|
| 0 | `00_setup.sbatch` | Pulls containers, builds the miRBase 23.1 *hsa* + *bkv* reference and the BKPyV/JCPyV bowtie2 index, and writes `sample_manifest.tsv` from the sample sheet and the ENA run table |
| 1 | `01_mirna.sbatch` | Small RNA (R1 only): cutadapt (TruSeq small-RNA adapter, untrimmed reads discarded, 17–30 nt, Q20), then miRDeep2 `quantifier.pl` against one combined human + BKPyV reference |
| 1b | `01b_mirtrace.sbatch` | miRTrace library-composition QC |
| 2 | `02_metagenome.sbatch` | Metagenome: fastp, human read removal (bowtie2, GRCh38 no-alt), BKPyV/JCPyV read counting and removal, MetaPhlAn 4.2.6 (`mpa_vJan25_CHOCOPhlAnSGB_202503`) |
| 3 | `03_collect_qc.py` | Merges outputs into count matrices and a per-sample QC table |

Every job writes a `run_manifest.json` next to its outputs, recording inputs, parameters and container images.

### Running

1. Edit the paths in `pipeline/config.sh`, including input folders and the shared database locations.
2. Provide a sample sheet `sample,group`, with `group` set to `treatment` or `control`.
3. Submit the jobs:

```bash
S=/path/to/pipeline
setup=$(sbatch --parsable --export=ALL,SCRIPTS=$S $S/00_setup.sbatch)
sbatch --dependency=afterok:$setup --export=ALL,SCRIPTS=$S $S/01_mirna.sbatch
sbatch --dependency=afterok:$setup --export=ALL,SCRIPTS=$S $S/01b_mirtrace.sbatch
sbatch --dependency=afterok:$setup --export=ALL,SCRIPTS=$S $S/02_metagenome.sbatch
# when finished
V2=/path/to/output python3 $S/03_collect_qc.py
```

The SLURM headers (partition, log paths) are set for the University of Zurich ScienceCluster; adjust them for other systems.

## Data availability

- Shotgun metagenomes: NCBI BioProject [PRJNA1252894](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1252894).
- Small RNA sequencing: available from the corresponding author upon reasonable request.

No patient-level clinical data are included in this repository.

## License

MIT, see [LICENSE](LICENSE).
