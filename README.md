# genome-loom

genome-loom creates comparative genome ribbon plots that show how a reference genome aligns with multiple comparison genomes across overview, pairwise, and neighbor-chain views. The reference remains at the top of the stack, comparison genomes remain in input order, and the figure set changes only which ribbon layer is visible. `genome-loom` preserves genome order, contig order, and contig orientation.

| Overview | Reference-Pairs | All-Pairs | Neighbor |
| :-: | :-: | :-: | :-: |
| [![overview](examples/output/light/overview/reference-vs-all.png)](examples/output/light/overview/reference-vs-all.png) | [![reference pair](examples/output/light/reference_pairs/reference-vs-comparison_alpha.png)](examples/output/light/reference_pairs/reference-vs-comparison_alpha.png) | [![all pairs](examples/output/light/all_pairs/comparison_alpha-vs-comparison_beta.png)](examples/output/light/all_pairs/comparison_alpha-vs-comparison_beta.png) | [![neighbor](examples/output/light/neighbor/neighbor-chain.png)](examples/output/light/neighbor/neighbor-chain.png) |

- `overview`: one full-stack figure showing ribbons from the reference to all comparison genomes at once.
- `reference-pairs`: one two-row figure per reference-to-comparison relationship, showing the reference and one selected comparison genome.
- `all-pairs`: one two-row figure per genome pair, showing only the two genomes involved in that selected comparison.
- `neighbor`: one full-stack figure showing ribbons only between adjacent rows in the supplied genome order.
- Comparison contigs are painted from their relationship to the reference, so the same reference-based color context can carry across the full figure set, including pairwise views where the reference row is not shown.
- Ribbon color follows the upper genome in each join: when the reference is the upper genome, ribbons use the reference contig color in that aligned region; when a non-reference genome is the upper genome, ribbons use that genome's reference-based painted color in the aligned region, or fall back to the neutral comparison color if no reference-based paint is present there.

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/paulstothard/genome-loom.git
cd genome-loom
```

### 2. Create the conda environment

```bash
conda env create -f environment.yml
```

This installs Python, matplotlib, Biopython, minimap2, and the other required packages into a self-contained `genome-loom` environment.

### 3. Verify

```bash
conda run -n genome-loom minimap2 --version
conda run -n genome-loom python genome_loom.py --help
```

## Quick Start

```bash
conda activate genome-loom

python genome_loom.py \
  --reference ref.fasta \
  --comparisons genome_a.fasta genome_b.fasta genome_c.fasta \
  --outdir results \
  --views overview reference-pairs neighbor \
  --theme light
```

This writes:

```text
results/
  genome-loom.summary.json
  overview/
    reference-vs-all.png
  reference_pairs/
    reference-vs-genome_a.png
    reference-vs-genome_b.png
    reference-vs-genome_c.png
  neighbor/
    neighbor-chain.png
```

`--output one-plot.png` is still available as a compatibility shortcut for a single overview image, but `--outdir` is the preferred workflow.

## Figure Views

- `overview`: full stack with reference-to-all ribbons.
- `reference-pairs`: one two-row figure per reference-to-comparison relationship.
- `all-pairs`: one two-row figure per selected genome pair.
- `neighbor`: one full-stack neighbor-chain figure with ribbons between adjacent rows.

The two pairwise view families are intentionally easier to read: each image shows only the two genomes involved in that selected comparison, which creates more vertical separation between rows and makes ribbon geometry easier to interpret than in the full stacked views.

## Main Options

`genome_loom.py` is the main single-entry script. It behaves like a predictable, non-interactive batch command for local use or server-side pipelines.

```bash
# Full figure set
python genome_loom.py \
  --reference my_reference.fasta \
  --comparisons my_genomes/ \
  --outdir results

# SVG output, dark theme, all figure families
python genome_loom.py \
  --reference ref.fasta \
  --comparisons set_a/ set_b/ outgroup.fasta \
  --outdir results \
  --views overview reference-pairs all-pairs neighbor \
  --theme dark \
  --format svg \
  --width 14 --height 9 --dpi 300 \
  --title "Genome Loom Example"

# Override display labels used in row labels and titles
python genome_loom.py \
  --reference ref.fasta \
  --comparisons comparison_a.fasta comparison_b.fasta \
  --outdir results \
  --display-names \
    "ref=NC_000913.3" \
    "comparison_a=E. coli O157:H7" \
    "comparison_b=E. coli CFT073"

# Change the role-label prefix shown on the top row
python genome_loom.py \
  --reference ref.fasta \
  --comparisons comparison_a.fasta \
  --outdir results \
  --reference-role-label assembly

# Focus the figure on selected top/reference contigs
python genome_loom.py \
  --reference assembly.fasta \
  --reference-contigs chromosome plasmid_A \
  --comparisons comparisons/ \
  --outdir focused-results

# Server-style run with explicit summary and work directory
python genome_loom.py \
  --reference ref.fasta \
  --comparisons comparisons/ \
  --outdir results \
  --summary-output results/results.json \
  --work-dir results/work \
  --threads 8 \
  --force

# JSON-config driven run
python genome_loom.py --config run.json --outdir results --force
```

Display labels are used exactly as supplied. genome-loom does not rewrite
underscores, dots, or other accession punctuation, so accession labels such as
`GCF_001900355.1` can be shown without renaming input files.

Key options are summarized below; run `python genome_loom.py --help` for the full reference.

| Option | Required | Default | Description |
| --- | --- | --- | --- |
| `--reference` | Yes, unless `--check` or `--version` | — | Reference genome FASTA. |
| `--comparisons` | Yes, unless `--check` or `--version` | — | Comparison FASTAs: files and/or directories scanned one level deep for `.fa`, `.fasta`, `.fna`, or `.fas` files. |
| `--outdir` | Required unless `--output`, `--check`, or `--version` is used | — | Output directory for figure families and summary JSON. |
| `--output` | Required unless `--outdir`, `--check`, or `--version` is used | — | Compatibility shortcut for one overview image. |
| `--summary-output` | No | `OUTDIR/genome-loom.summary.json`, or `OUTPUT` with `.json` suffix | Machine-readable JSON summary path. |
| `--display-names` | No | — | Optional `KEY=LABEL` overrides keyed by FASTA stem or filename for figure labels and titles. |
| `--reference-role-label` | No | `reference` | Optional role-label prefix for the top row; use `none` to omit the prefix. |
| `--reference-contigs` | No | — | Keep only the named contigs from the top/reference genome. Aliases: `--assembly-contigs`, `--query-contigs`. |
| `--views` | No | `overview reference-pairs neighbor` | Figure families to generate when `--outdir` is used. |
| `--format` | No | `png` | Output format: `png`, `pdf`, or `svg`. |
| `--theme` | No | `light` | Figure theme: `light` or `dark`. |
| `--width` / `--height` | No | `12` / `8` | Figure dimensions in inches. |
| `--dpi` | No | `300` | Output resolution. |
| `--title` | No | — | Figure title override. |
| `--min-contig-length` | No | `1000` | Discard contigs shorter than this before rendering. |
| `--max-contigs` | No | `0` | Cap visible contig blocks per genome; `0` keeps all contigs, and `1` is rejected. |
| `--reference-segments` | No | `0` | Split the visible reference into this many equal-length colored segments; `0` keeps standard contig-based coloring. |
| `--min-block-length` | No | `500` | Discard short alignment blocks. |
| `--minimap-preset` | No | `asm5` | minimap2 preset: `asm5`, `asm10`, or `asm20`. |
| `--min-mapq` | No | `0` | Discard low-confidence alignment blocks. |
| `--threads` | No | `1` | Thread count passed to minimap2. |
| `--tmpdir` | No | System temp directory | Parent directory for an auto-created temporary work directory. |
| `--work-dir` | No | Auto-created temporary directory | Explicit intermediate directory for prepared FASTAs and cached files. |
| `--keep-temp` | No | Off | Keep an auto-created temporary work directory after the run. |
| `--force` | No | Off | Allow overwriting existing outputs and reusing a non-empty `--work-dir`. |
| `--config` | No | — | JSON config file. Keys map to CLI options by replacing `_` with `-`; command-line arguments override config values. |
| `--check` | No | Off | Check wrapper and external-tool availability, then exit. |
| `--version` | No | — | Print version, then exit. |

## Output Layout

Typical output when `--outdir` is used:

```text
results/
├── genome-loom.summary.json
├── overview/
│   └── reference-vs-all.png
├── reference_pairs/
│   └── reference-vs-comparison.png
├── all_pairs/
│   └── comparison_a-vs-comparison_b.png
└── neighbor/
    └── neighbor-chain.png
```

The exact figure directories depend on `--views`. When `--output` is used instead of `--outdir`, the wrapper writes one overview image and a sibling summary JSON path derived from the output filename unless `--summary-output` is supplied.

## Server And Batch Use

`genome_loom.py` is the preferred entry point for an analysis server. It requires explicit inputs, never prompts interactively, exits nonzero on failure, and writes a machine-readable JSON summary for the caller.

Example server-style command:

```bash
python genome_loom.py \
  --reference /job/reference.fasta \
  --comparisons /job/comparisons \
  --outdir /job/results \
  --summary-output /job/results/results.json \
  --work-dir /job/results/work \
  --threads 8 \
  --force
```

Practical recommendations for predictable server-side operation:

- Pass absolute paths or paths inside the job directory.
- Run `python genome_loom.py --check` in the target conda environment before accepting jobs.
- Use `--summary-output` with a known filename such as `results.json`.
- Use `--work-dir` if prepared FASTAs should be retained for debugging or caching.
- Use `--tmpdir` when temporary files should live on fast local scratch storage.
- Use `--force` when re-running into an existing results directory or work directory.

If `--work-dir` is omitted, `genome-loom` creates a temporary working directory. That directory is deleted after the run unless `--keep-temp` is used. If `--keep-temp` is supplied without `--work-dir`, the auto-created temporary work directory is kept and reported in the summary JSON.

The JSON summary includes:

- `status`, `tool`, `version`, and timestamp
- input paths and output paths
- figure-generation settings
- filtering and alignment settings
- genome order and per-genome contig summaries
- one record per generated figure
- one record per computed pairwise alignment
- warnings, such as likely crowding in full-stack views

On failure, the summary still attempts to record the error type, message, and the output paths already known to the wrapper.

If `--reference-contigs` is used, only those named contigs from the top/reference genome are prepared, aligned, summarized, and rendered. This is useful when the top genome is a large assembly but the figure should concentrate on one chromosome, plasmid, or a small hand-picked set of contigs.

### JSON Config Files

`--config` accepts a JSON object that uses option names as keys. For example:

```json
{
  "reference": "ref.fasta",
  "comparisons": ["comparison_a.fasta", "comparison_b.fasta"],
  "outdir": "results",
  "views": ["overview", "reference-pairs", "neighbor"],
  "theme": "light",
  "reference_segments": 10,
  "max_contigs": 12,
  "threads": 4
}
```

Command-line arguments override values loaded from the config file, so a caller can keep a stable base config and customize only a few fields per job.

## Alignment Settings

`genome-loom` currently uses `minimap2` for all alignments. The main alignment mode option is the minimap2 preset:

- `--minimap-preset asm5`: best starting point for very close assemblies or strain-level comparisons. In the figures, this usually gives the cleanest long ribbons and the least spurious low-similarity clutter, but it can drop genuinely homologous regions once divergence becomes moderate.
- `--minimap-preset asm10`: a middle ground for somewhat more diverged comparisons. In the figures, this often restores ribbons that `asm5` missed, but the recovered alignments may be shorter, more broken up, or more repetitive-looking than in a very close comparison.
- `--minimap-preset asm20`: most permissive of the three standard assembly presets and best for more divergent exploratory work. In the figures, this can recover coarse synteny and shared regions across larger divergence, but ribbons are more likely to be fragmented, partial, or visually busy, so interpretation should be more cautious.

Other alignment-related options:

- `--min-block-length`: discard short alignment blocks after minimap2 runs.
- `--min-mapq`: discard low-confidence blocks by mapping quality.
- `--threads`: controls the minimap2 thread count.

Practical rule of thumb:

- Use `asm5` first for close same-species or strain comparisons when you want the cleanest and most conservative ribbon structure.
- Move to `asm10` when expected homologous regions start disappearing, or when the figure looks too sparse even though the genomes should still share substantial structure.
- Move to `asm20` for clearly more divergent genomes when the goal is to recover broad shared structure rather than only the strongest high-confidence blocks.
- If a more permissive preset adds many short, broken, or noisy ribbons, raise `--min-block-length` and/or `--min-mapq` to keep the figure readable.
- If a permissive preset changes the figure dramatically, treat the result as a broader structural sketch rather than a crisp one-to-one alignment view.

## Color Propagation and Figure Meaning

For `overview`, `reference-pairs`, `neighbor`, and any `all-pairs` figure that includes the reference, ribbons use reference-flow colors. A reference contig color follows the aligned DNA into comparison rows, and comparison contig bars are painted where direct reference-alignment evidence indicates reference-based sequence context is present.

That propagated comparison-contig coloring is global across the figure set, so it remains visible even when the current image is showing only one selected pair or a neighbor-chain ribbon layer. In `neighbor`, those colors can continue to percolate downward through adjacent comparisons.

By default, the figure row label for the top genome reads `reference | <name>` so the viewer can see both the role and the actual identifier. `--reference-role-label` can override that prefix, and `--reference-role-label none` omits it entirely. In pairwise figures that do not show the reference row, the legend note still names the reference genome that supplied the propagated color context.

For `all-pairs` images that do not include the reference, ribbons are colored locally from the upper genome in that selected pair because the ribbon itself no longer contains a direct reference thread. The comparison contig bars are still painted from direct reference evidence when available. The JSON render metadata records ribbon coloring as either `reference-flow` or `subject-local`, and comparison contig coloring as `reference-based` when global reference coloring is active.

Reference colors are assigned by contig size, not FASTA order. The largest reference contig receives the first palette color, the next largest receives the second, and so on. Once the palette is exhausted, the remaining smaller contigs share the fallback color. That fallback color still flows through matches, but it no longer distinguishes which individual small contig contributed a given fallback thread.

### Segment-Based Reference Coloring

By default, reference-flow coloring is contig-based: each visible reference contig receives one color, and that color is propagated through alignments. Use `--reference-segments N` to switch to segment-based coloring instead. In this mode, the whole visible reference is divided into `N` equal-length colored intervals before reference-flow coloring is computed.

Segmenting works for both single-contig and multicontig references. For a multicontig reference, the segment boundaries are computed across the concatenated visible reference length, using the displayed contig order. A segment can therefore end on one contig and continue on the next; genome-loom records and renders that as separate per-contig intervals with the same segment label and color.

For example, this command renders a neighbor-chain figure where a single complete reference chromosome is divided into ten colored segments:

```bash
python genome_loom.py \
  --reference examples/case_studies/ecoli_complete_reference/reference.fasta \
  --comparisons examples/case_studies/ecoli_complete_reference/comparisons/ \
  --outdir results-segmented \
  --views neighbor \
  --theme dark \
  --reference-segments 10 \
  --minimap-preset asm10 \
  --threads 4 \
  --force
```

The segment colors behave like reference contig colors everywhere else in the figure set:

- The top reference row is painted as `segment 1`, `segment 2`, and so on from left to right.
- On multicontig references, segment labels continue across contig boundaries in displayed reference order.
- Alignment blocks that cross segment boundaries are split at those boundaries before colors are propagated.
- Ribbons use the color of the reference segment that contributed each aligned interval.
- Comparison contig bars receive propagated segment colors where direct reference-alignment evidence supports them.
- In `neighbor`, segment colors can continue to move through adjacent rows, just like contig colors do.
- The legend switches from contig names to segment labels so the color key matches the segmented reference row.

This example renders a fragmented, multicontig reference using twelve equal-length reference segments instead of contig colors:

```bash
python genome_loom.py \
  --reference examples/case_studies/fragmented_reference/reference.fasta \
  --comparisons examples/case_studies/fragmented_reference/comparisons/ \
  --outdir fragmented-segments \
  --views overview neighbor \
  --theme dark \
  --reference-segments 12 \
  --max-contigs 0 \
  --minimap-preset asm10 \
  --threads 4 \
  --force
```

If you want to segment only one contig from a multicontig reference, first select it explicitly:

```bash
python genome_loom.py \
  --reference assembly.fasta \
  --reference-contigs chromosome \
  --comparisons comparisons/ \
  --outdir chromosome-segments \
  --views overview neighbor \
  --reference-segments 12
```

`--reference-segments 0` disables segmentation. `--reference-segments 1` is rejected because it would be visually equivalent to ordinary single-contig coloring. Segment counts greater than the built-in palette size are allowed; colors cycle after the available distinct palette colors have been used.

## Contig Capping and Fragmented Assemblies

Use `--max-contigs` to keep highly fragmented assemblies readable.

- `0` keeps all contigs.
- `N` keeps the largest `N - 1` contigs as separate visible blocks.
- The remaining smaller contigs are merged into one trailing `remaining_contigs_N` block.

The kept contigs preserve their original relative FASTA order. Smaller contigs that are not kept as separate blocks are moved into the trailing merged block. That merged block is aligned and can carry propagated reference color, but it always uses the fallback color because individual small-contig contributions are no longer distinguishable within that bin.

## Practical Recommendations

`genome-loom` is designed for microbial-scale genome comparisons, including small plasmids through typical bacterial chromosomes. Larger assemblies can work, but figure density and rendering time increase as genome size, genome count, and contig count all rise together.

### Genome Size

- Best fit: small plasmids through typical bacterial genomes, roughly tens of kb to about 10 Mb.
- Still practical with care: larger bacterial assemblies and some small eukaryotic scaffolds up to roughly 100 Mb.
- Use caution beyond that: figures can become sparse, slow to render, and harder to interpret unless heavily filtered.

### Total Genome Count Per Figure

For full-stack views such as `overview` and `neighbor`, a good rule of thumb is:

```text
recommended total genomes <= figure height in inches / 0.6
```

So for the default `12 x 8` inch figure, aim for roughly 13 genomes or fewer in the full stacked views. In practice:

- Best readability: about 4 to 10 total genomes.
- Usually still workable: about 11 to 13 total genomes, especially with short names and contig capping.
- Above that: increase figure height, reduce the genome set, or lean more heavily on the pairwise view families.

`genome-loom` prints a warning to `stderr` and records it in the summary JSON when the selected full-stack views are likely to look crowded at the chosen height.

### Contigs Per Genome

- Best readability: keep visible contig blocks in the single digits to low teens.
- For fragmented bacterial assemblies, `--max-contigs 6` to `--max-contigs 24` is often a good range.
- If many genomes are shown at once, lean toward lower `--max-contigs` values.

### Labels and Long Names

Genome row labels shrink as figures get more crowded, and contig-legend items wrap onto additional rows. Long legend labels are truncated with an ellipsis when needed. If row labels or contig legends start to feel dense, the best fixes are:

- increase figure height
- reduce the genome set for full-stack views
- lower `--max-contigs`
- rely on pairwise views for detailed interpretation

## Examples

Generate deterministic synthetic genomes plus two real-data example sets:

```bash
bash rebuild_example_outputs.sh
```

This writes:

```text
examples/data/
examples/case_studies/ecoli_complete_reference/
examples/case_studies/fragmented_reference/
examples/output/light/
examples/output/dark/
examples/output/ecoli_complete_reference/
examples/output/fragmented_reference/
```

The real-data examples are copied from local `genome-artistry` example datasets:

- `ecoli_complete_reference`: one complete E. coli reference plus three comparisons.
- `fragmented_reference`: the most fragmented local FASTA found in the neighboring `genome-artistry` example/montage cache plus three nearby comparisons.

These example sets are useful for checking how contig capping, propagated reference color, and ribbon readability behave on both simple and fragmented assemblies.

## Dependencies

- **minimap2** 2.30 — sequence alignment for genome ribbon figures. [GitHub](https://github.com/lh3/minimap2) · [Li 2018](https://doi.org/10.1093/bioinformatics/bty191)
