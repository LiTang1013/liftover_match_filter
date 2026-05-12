# Modular mtDNA liftover + codon-match + tRNA-match pipeline

This package refactors the previous integrated liftover/codon script into a config-driven pipeline.

## Pipeline logic

```text
pairwise alignment
  -> alignment-derived species->human position map
  -> liftover-only raw VCF/COV
  -> optional codon-match annotation
  -> optional tRNA-match annotation
  -> optional final INFO-based filter
```

The raw lifted VCF is intentionally not codon-filtered or tRNA-filtered. It stores original species coordinates and mapping fields in `MTLIFT_*` INFO tags. Downstream codon and tRNA steps add `MTCODON_*` and `MTTRNA_*` INFO tags but do not drop records.

## Files

```text
mt_pipeline.py              # self-contained Python pipeline
config.example.ini          # all paths and switches
run_mt_pipeline_config.sh   # SLURM array wrapper
```

## Quick start

Edit paths in `config.example.ini`, then run one sample:

```bash
python3 mt_pipeline.py run-sample --config config.example.ini --sample panTro6
```

Run as a SLURM array:

```bash
sbatch --array=0-158 run_mt_pipeline_config.sh
```

or override config path:

```bash
CONFIG=/path/to/config.ini sbatch --array=0-158 run_mt_pipeline_config.sh
```

List samples discovered from `sample_list_file` or `primate_dir`:

```bash
python3 mt_pipeline.py list-samples --config config.example.ini
```

## tRNAscan-SE setup

tRNA coordinates must be in the same coordinate space as the VCF side being annotated:

* species tRNAscan output: same coordinate system as the original species VCF positions stored in `MTLIFT_ORIG_POS`;
* human tRNAscan output: same coordinate system as the final human VCF POS.

Run tRNAscan-SE for human:

```bash
python3 mt_pipeline.py run-trnascan --config config.example.ini --sample human
python3 mt_pipeline.py build-trna-index --config config.example.ini --sample human
```

Run tRNAscan-SE for one species:

```bash
python3 mt_pipeline.py run-trnascan --config config.example.ini --sample panTro6
python3 mt_pipeline.py build-trna-index --config config.example.ini --sample panTro6
```

If `run_trnascan_if_missing = 1`, `run-sample` will try to run tRNAscan-SE and build the species index when missing.

## Important output INFO fields

### Liftover-only raw VCF

```text
MTLIFT_ORIG_CHROM
MTLIFT_ORIG_POS
MTLIFT_ORIG_REF
MTLIFT_ORIG_ALT
MTLIFT_ORIG_ROT_POS
MTLIFT_HUMAN_ROT_POS
MTLIFT_HUMAN_POS
MTLIFT_USABLE
MTLIFT_REFALT_SWAPPED
```

### Codon annotation

Strict codon mode is controlled by `strict_phase_match = 1` in `config.ini`.
When it is enabled, `MTCODON_STATUS=PASS` requires all of the following:

```text
1. same mitochondrial gene after gene-name harmonization
2. same codon_pos_in_triplet, i.e. source and human are both codon base 1, 2, or 3
3. primate/source REF codon or ALT codon equals the human codon
```

Output INFO fields:

```text
MTCODON_STATUS=PASS|SKIPPED_NONCODING|NO_HUMAN_CODON|GENE_MISMATCH|PHASE_MISMATCH|MISMATCH|MISSING_COORD
MTCODON_MATCH=ref_codon|alt_codon|NA
MTCODON_STRICT_PHASE=1|0
MTCODON_GENE_MATCH=yes|no|NA
MTCODON_PHASE_MATCH=yes|no|NA
MTCODON_PRIMATE_GENE
MTCODON_PRIMATE_CODON
MTCODON_PRIMATE_PHASE
MTCODON_HUMAN_GENE
MTCODON_HUMAN_CODON
MTCODON_HUMAN_PHASE
```

### tRNA annotation

```text
MTTRNA_STATUS=OK|NO_SPECIES_TRNA|NO_HUMAN_TRNA|NO_SPECIES_OR_HUMAN_TRNA|MISSING_SPECIES_COORD
MTTRNA_S_CLASS / MTTRNA_H_CLASS
MTTRNA_REGION_MATCH
MTTRNA_S_PAIR_TYPE / MTTRNA_H_PAIR_TYPE
MTTRNA_PAIR_TYPE_MATCH
MTTRNA_S_PAIR_STATE / MTTRNA_H_PAIR_STATE
MTTRNA_PAIR_STATE_MATCH
MTTRNA_S_PAIR_POS
MTTRNA_H_PAIR_POS
MTTRNA_S_PAIR_LIFTED_HPOS
MTTRNA_PAIR_POS_MATCH
```

## Single-step commands

Build or rebuild only the tRNA index:

```bash
python3 mt_pipeline.py build-trna-index --config config.example.ini --sample panTro6
```

Annotate one already lifted VCF with codon-match:

```bash
python3 mt_pipeline.py annotate-codon \
  --config config.example.ini \
  --sample panTro6 \
  --input sample.lifted.raw.vcf \
  --output sample.lifted.codon.vcf
```

Annotate one lifted VCF with tRNA-match:

```bash
python3 mt_pipeline.py annotate-trna \
  --config config.example.ini \
  --sample panTro6 \
  --input sample.lifted.codon.vcf \
  --output sample.lifted.codon.trna.vcf
```

Filter after annotation:

```bash
python3 mt_pipeline.py filter-vcf \
  --config config.example.ini \
  --mode codon_or_trna \
  --input sample.lifted.codon.trna.vcf \
  --output sample.final.vcf
```

Supported filter modes:

```text
none
codon_pass
trna_region_match
trna_pair_state_match
trna_pair_pos_match
codon_or_trna
trna_loose_match
trna_strict_match
region_policy
```

`region_policy` applies region-aware final filtering:
1) coding variants require `MTCODON_STATUS=PASS`;
2) tRNA variants require `MTTRNA_STRICT_MATCH=yes`;
3) other noncoding variants (including control region) are kept.

`MTTRNA_STRICT_MATCH` logic:
- loop-loop: `MTTRNA_REGION_MATCH=yes` AND `MTTRNA_ELEMENT_MATCH=yes`
- stem-stem: `MTTRNA_REGION_MATCH=yes` AND `MTTRNA_ELEMENT_MATCH=yes` AND `MTTRNA_PAIR_STATE_MATCH=yes` AND `MTTRNA_PAIR_POS_MATCH=yes`
- mixed/other cases: `no`

## Notes

The current liftover stage is SNV-oriented because it harmonizes REF/ALT to the human reference and recodes genotypes when REF/ALT are swapped. Non-SNV or multi-allelic records are written to the `.dropped.vcf` output.

The tRNA paired-position check uses the same alignment-derived `posmap.tsv.gz` used by liftover. When `species_trna_coord_space=original`, source paired genomic positions are rotated before liftover; when `species_trna_coord_space=rotated`, paired positions are used as-is to avoid double rotation.

For rotated species tRNA indexes (for example `chrom=Allenopithecus_nigroviridis`, `pos=<rotated position>`), set:

```ini
[settings]
species_trna_coord_space = rotated
species_trna_lookup_ignore_chrom = 1
human_trna_lookup_ignore_chrom = 0
```

`species_trna_lookup_ignore_chrom=1` makes species tRNA lookup key on position only, and `species_trna_coord_space=rotated` switches species lookup to `INFO/MTLIFT_ORIG_ROT_POS`.
