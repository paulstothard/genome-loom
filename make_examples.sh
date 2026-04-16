#!/usr/bin/env bash
set -euo pipefail

EXAMPLES_DIR="examples"
DATA_DIR="${EXAMPLES_DIR}/data"
OUT_DIR="${EXAMPLES_DIR}/output"
ENV_NAME="genome-loom"
PYTHON=(conda run -n "${ENV_NAME}" python)

mkdir -p "${DATA_DIR}" "${OUT_DIR}"

"${PYTHON[@]}" scripts/make_demo_data.py --out-dir "${DATA_DIR}"

rm -rf "${OUT_DIR}/light" "${OUT_DIR}/dark"
rm -rf "${OUT_DIR}/ecoli_complete_reference"
rm -rf "${OUT_DIR}/fragmented_reference"

"${PYTHON[@]}" genome_loom.py \
  --reference "${DATA_DIR}/reference.fasta" \
  --comparisons \
  "${DATA_DIR}/comparison_alpha.fasta" \
  "${DATA_DIR}/comparison_beta.fasta" \
  "${DATA_DIR}/comparison_gamma.fasta" \
  --outdir "${OUT_DIR}/light" \
  --summary-output "${OUT_DIR}/light/genome-loom.summary.json" \
  --views overview reference-pairs all-pairs neighbor \
  --format png \
  --theme light \
  --width 12 \
  --height 8 \
  --dpi 180 \
  --min-contig-length 1000 \
  --max-contigs 0 \
  --min-block-length 800 \
  --threads 2 \
  --title "Genome Loom Demo"

"${PYTHON[@]}" genome_loom.py \
  --reference "${DATA_DIR}/reference.fasta" \
  --comparisons \
  "${DATA_DIR}/comparison_alpha.fasta" \
  "${DATA_DIR}/comparison_beta.fasta" \
  "${DATA_DIR}/comparison_gamma.fasta" \
  --outdir "${OUT_DIR}/dark" \
  --summary-output "${OUT_DIR}/dark/genome-loom.summary.json" \
  --views overview reference-pairs neighbor \
  --format png \
  --theme dark \
  --width 12 \
  --height 8 \
  --dpi 180 \
  --min-contig-length 1000 \
  --max-contigs 0 \
  --min-block-length 800 \
  --threads 2 \
  --title "Genome Loom Demo"

"${PYTHON[@]}" genome_loom.py \
  --reference "${EXAMPLES_DIR}/case_studies/ecoli_complete_reference/reference.fasta" \
  --comparisons \
  "${EXAMPLES_DIR}/case_studies/ecoli_complete_reference/comparisons/GCF_003073835_1.fasta" \
  "${EXAMPLES_DIR}/case_studies/ecoli_complete_reference/comparisons/GCF_002854065_1.fasta" \
  "${EXAMPLES_DIR}/case_studies/ecoli_complete_reference/comparisons/GCF_001900355_1.fasta" \
  --outdir "${OUT_DIR}/ecoli_complete_reference" \
  --summary-output "${OUT_DIR}/ecoli_complete_reference/genome-loom.summary.json" \
  --views overview reference-pairs all-pairs neighbor \
  --format png \
  --theme light \
  --width 13 \
  --height 8 \
  --dpi 180 \
  --min-contig-length 1000 \
  --max-contigs 6 \
  --min-block-length 5000 \
  --threads 2 \
  --title "Real E. coli Example"

"${PYTHON[@]}" genome_loom.py \
  --reference "${EXAMPLES_DIR}/case_studies/fragmented_reference/reference.fasta" \
  --comparisons \
  "${EXAMPLES_DIR}/case_studies/fragmented_reference/comparisons/GCF_002589795_1.fasta" \
  "${EXAMPLES_DIR}/case_studies/fragmented_reference/comparisons/GCF_002011945_1.fasta" \
  "${EXAMPLES_DIR}/case_studies/fragmented_reference/comparisons/GCF_002854065_1.fasta" \
  --outdir "${OUT_DIR}/fragmented_reference" \
  --summary-output "${OUT_DIR}/fragmented_reference/genome-loom.summary.json" \
  --views overview reference-pairs all-pairs neighbor \
  --format png \
  --theme light \
  --width 13 \
  --height 8 \
  --dpi 180 \
  --min-contig-length 1000 \
  --max-contigs 6 \
  --min-block-length 5000 \
  --threads 2 \
  --title "Most-Contig Local Example"

echo "Examples written to ${OUT_DIR}/light, ${OUT_DIR}/dark, ${OUT_DIR}/ecoli_complete_reference, and ${OUT_DIR}/fragmented_reference"
