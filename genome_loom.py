#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import json
import re
import shutil
import sys
import tempfile
import textwrap
from dataclasses import replace
from pathlib import Path

from scripts.align import PairwiseAlignment, run_minimap2
from scripts.fasta import (
    Genome,
    cap_contig_blocks,
    read_fasta,
    select_contigs,
    write_fasta,
)
from scripts.summary import write_summary


VERSION = "0.5.0"
VIEW_CHOICES = ("overview", "reference-pairs", "all-pairs", "neighbor")
MINIMAP2_PRESETS = ("asm5", "asm10", "asm20")
GENOME_ORDER_CHOICES = ("input", "reference-similarity")
REFERENCE_PALETTE_CHOICES = ("categorical", "continuous")
FULL_STACK_VIEWS = {"overview", "neighbor"}


def _collect_comparisons(values: list[str]) -> list[Path]:
    paths: list[Path] = []
    exts = {".fa", ".fasta", ".fna", ".fas"}
    for value in values:
        p = Path(value)
        if p.is_dir():
            paths.extend(sorted(x for x in p.iterdir() if x.suffix.lower() in exts))
        else:
            paths.append(p)
    seen: set[Path] = set()
    unique: list[Path] = []
    for p in paths:
        rp = p.resolve()
        if rp not in seen:
            seen.add(rp)
            unique.append(rp)
    return unique


def _summary_for_genome(genome: Genome) -> dict:
    return {
        "name": genome.name,
        "display_name": genome.display_name or genome.name,
        "path": str(genome.path),
        "length_bp": genome.length,
        "contig_count": len(genome.contigs),
        "remainder_contig_count": sum(
            c.source_count for c in genome.contigs if c.is_remainder
        ),
        "contigs": [
            {
                "name": c.name,
                "length_bp": c.length,
                "source_count": c.source_count,
                "is_remainder": c.is_remainder,
            }
            for c in genome.contigs
        ],
    }


def _figure_name(genome: Genome) -> str:
    return genome.display_name or genome.name


def _safe_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    return cleaned.strip("._-") or "genome"


def _alignment_summary(aln: PairwiseAlignment) -> dict:
    return {
        "subject": aln.reference.name,
        "comparison": aln.comparison.name,
        "ani": round(aln.ani, 6),
        "aligned_bases": aln.aligned_bases,
        "matching_bases": aln.matching_bases,
        "block_count": len(aln.blocks),
    }


def _overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def _parse_display_name_pairs(values: list[str] | None) -> dict[str, str]:
    overrides: dict[str, str] = {}
    if not values:
        return overrides
    for value in values:
        if "=" not in value:
            raise argparse.ArgumentTypeError(
                f"Malformed --display-names entry {value!r}; expected KEY=LABEL"
            )
        key, label = value.split("=", 1)
        key = key.strip()
        label = label.strip()
        if not key or not label:
            raise argparse.ArgumentTypeError(
                f"Malformed --display-names entry {value!r}; expected KEY=LABEL"
            )
        overrides[key] = label
    return overrides


def _parse_contig_name_filters(values: list[str] | None) -> list[str]:
    if not values:
        return []
    names: list[str] = []
    for value in values:
        for item in value.split(","):
            name = item.strip()
            if name and name not in names:
                names.append(name)
    if not names:
        raise argparse.ArgumentTypeError(
            "contig selection must include at least one contig name"
        )
    return names


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return parsed


def _nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be zero or greater")
    return parsed


def _positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive number")
    return parsed


def _check_tool(name: str) -> dict:
    resolved = shutil.which(name)
    return {
        "requested": name,
        "resolved": resolved,
        "available": resolved is not None,
    }


def _write_check() -> int:
    payload = {
        "tool": "genome-loom",
        "version": VERSION,
        "checks": {
            "minimap2": _check_tool("minimap2"),
        },
    }
    print(json.dumps(payload, indent=2))
    return 0


def _config_args_from_file(path: Path) -> list[str]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError("Config file must contain a JSON object")
    args: list[str] = []
    for key, value in payload.items():
        option = f"--{key.replace('_', '-')}"
        if value is None:
            continue
        if isinstance(value, bool):
            if value:
                args.append(option)
            continue
        if isinstance(value, list):
            if not value:
                continue
            args.append(option)
            args.extend(str(item) for item in value)
            continue
        args.extend([option, str(value)])
    return args


def _parse_args(
    parser: argparse.ArgumentParser, argv: list[str] | None
) -> argparse.Namespace:
    raw = list(argv) if argv is not None else sys.argv[1:]
    bootstrap = argparse.ArgumentParser(add_help=False)
    bootstrap.add_argument("--config", type=Path)
    bootstrap_args, _unknown = bootstrap.parse_known_args(raw)
    if bootstrap_args.config is None:
        return parser.parse_args(raw)
    config_path = bootstrap_args.config.expanduser().resolve()
    if not config_path.exists():
        parser.error(f"--config does not exist: {config_path}")
    try:
        config_args = _config_args_from_file(config_path)
    except Exception as exc:
        parser.error(f"Could not parse --config {config_path}: {exc}")
    filtered_raw: list[str] = []
    skip_next = False
    for token in raw:
        if skip_next:
            skip_next = False
            continue
        if token == "--config":
            skip_next = True
            continue
        filtered_raw.append(token)
    return parser.parse_args(config_args + filtered_raw)


def _recommended_full_stack_limit(height_in: float) -> int:
    return max(4, int(height_in / 0.6))


def _reference_segment_intervals(
    *,
    reference: Genome,
    segment_count: int,
    colors: list[str],
) -> list:
    """Split the visible reference into equal-length colored intervals."""
    intervals = []
    total_length = sum(contig.length for contig in reference.contigs)
    if total_length <= 0:
        return intervals

    contig_offsets: list[tuple[str, int, int]] = []
    offset = 0
    for contig in reference.contigs:
        contig_offsets.append((contig.name, offset, offset + contig.length))
        offset += contig.length

    for index in range(segment_count):
        segment_start = round(index * total_length / segment_count)
        segment_end = round((index + 1) * total_length / segment_count)
        if segment_end <= segment_start:
            continue
        for contig_name, contig_start, contig_end in contig_offsets:
            start = max(segment_start, contig_start)
            end = min(segment_end, contig_end)
            if end <= start:
                continue
            intervals.append(
                {
                    "genome": reference.name,
                    "contig": contig_name,
                    "start": start - contig_start,
                    "end": end - contig_start,
                    "color": colors[index % len(colors)],
                    "origin": f"segment {index + 1}",
                }
            )
    return intervals


def _write_failure_summary(
    summary_path: Path,
    *,
    version: str,
    outdir: Path | None,
    output: Path | None,
    work_dir: Path | None,
    exc: Exception,
) -> None:
    payload = {
        "status": "error",
        "tool": "genome-loom",
        "version": version,
        "outputs": {
            "outdir": str(outdir) if outdir is not None else None,
            "output": str(output) if output is not None else None,
            "summary_json": str(summary_path),
            "work_dir": str(work_dir) if work_dir is not None else None,
        },
        "error": {
            "type": exc.__class__.__name__,
            "message": str(exc),
            "step": getattr(exc, "step", None),
            "returncode": getattr(exc, "returncode", None),
        },
    }
    write_summary(summary_path, payload)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genome_loom.py",
        formatter_class=lambda prog: argparse.RawDescriptionHelpFormatter(
            prog, max_help_position=34, width=92
        ),
        description=textwrap.dedent(
            """\
            Render genome ribbon figure sets from local FASTA files.

            The wrapper runs three steps as one non-interactive batch command:
              1. validate and prepare input FASTAs
              2. align genomes with minimap2
              3. render one or more figure views plus a JSON summary
            """
        ),
        epilog=textwrap.dedent(
            """\
            Server/batch tip:
              Use absolute paths inside a job directory, set --summary-output to a
              predictable filename such as results.json, and pass --work-dir when
              intermediate FASTAs should be retained for debugging or caching.

            JSON config files are supported via --config and use option names as keys,
            for example: {"reference": "...", "comparisons": ["a.fasta", "b.fasta"]}.
            Command-line arguments override values loaded from --config.
            """
        ),
    )
    parser.add_argument(
        "--config",
        type=Path,
        help="Optional JSON config file; command-line arguments override config values.",
    )
    parser.add_argument("--reference", type=Path, help="Reference genome FASTA.")
    parser.add_argument(
        "--reference-contigs",
        "--assembly-contigs",
        "--query-contigs",
        nargs="+",
        metavar="CONTIG",
        help=(
            "Optional contig names to keep from the top/reference genome only. "
            "Names can be supplied as separate values or comma-separated."
        ),
    )
    parser.add_argument(
        "--comparisons",
        nargs="+",
        help="Comparison FASTA files and/or directories of FASTA files.",
    )
    parser.add_argument(
        "--display-names",
        nargs="+",
        metavar="KEY=LABEL",
        help=(
            "Optional display-name overrides keyed by FASTA stem or filename, "
            "for example reference=K-12_MG1655 or NC_000913.3.fasta=K-12 MG1655."
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        help="Directory for generated figure sets and summary JSON.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Compatibility shortcut: write a single overview plot to this path.",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        help="Machine-readable JSON summary path. Defaults next to the output figure(s).",
    )
    parser.add_argument(
        "--views",
        nargs="+",
        choices=VIEW_CHOICES,
        default=["overview", "reference-pairs", "neighbor"],
        help="Figure sets to generate when --outdir is used.",
    )
    parser.add_argument(
        "--genome-order",
        choices=GENOME_ORDER_CHOICES,
        default="input",
        help=(
            "Row order for comparison genomes. Use input to preserve the supplied "
            "comparison order, or reference-similarity to place genomes with the "
            "highest direct ANI to the reference closest to the reference row."
        ),
    )
    parser.add_argument(
        "--format", choices=["png", "pdf", "svg"], default="png", help="Output format."
    )
    parser.add_argument(
        "--theme", choices=["light", "dark"], default="light", help="Figure theme."
    )
    parser.add_argument(
        "--width", type=_positive_float, default=12.0, help="Figure width in inches."
    )
    parser.add_argument(
        "--height", type=_positive_float, default=8.0, help="Figure height in inches."
    )
    parser.add_argument(
        "--dpi", type=_positive_int, default=300, help="Output resolution."
    )
    parser.add_argument(
        "--min-contig-length",
        type=_nonnegative_int,
        default=1000,
        help="Ignore contigs shorter than this many base pairs before rendering.",
    )
    parser.add_argument(
        "--max-contigs",
        type=_nonnegative_int,
        default=0,
        help=(
            "Maximum contig blocks to show per genome. Use 0 to keep all contigs. "
            "When capped, the largest max-1 contigs stay separate and the rest "
            "are merged into one trailing remaining_contigs_N block."
        ),
    )
    parser.add_argument(
        "--reference-segments",
        type=_nonnegative_int,
        default=0,
        help=(
            "Split the visible reference into this many colored segments. "
            "Use 0 to color by reference contig instead."
        ),
    )
    parser.add_argument(
        "--reference-palette",
        choices=REFERENCE_PALETTE_CHOICES,
        default="categorical",
        help=(
            "Color palette for reference contigs or reference segments. "
            "categorical uses distinct colors; continuous uses a theme-specific "
            "ordered gradient so neighboring contigs or segments have related colors."
        ),
    )
    parser.add_argument(
        "--min-block-length",
        type=_nonnegative_int,
        default=500,
        help="Ignore alignment blocks shorter than this many base pairs.",
    )
    parser.add_argument(
        "--minimap-preset",
        choices=MINIMAP2_PRESETS,
        default="asm5",
        help=(
            "minimap2 assembly-alignment preset. "
            "asm5 is the strictest and best for close genomes, "
            "asm10 is more tolerant, and asm20 is the loosest."
        ),
    )
    parser.add_argument(
        "--min-mapq",
        type=_nonnegative_int,
        default=0,
        help="Minimum minimap2 mapping quality to keep an alignment block.",
    )
    parser.add_argument(
        "--threads",
        type=_positive_int,
        default=1,
        help="Thread count passed to minimap2.",
    )
    parser.add_argument(
        "--tmpdir",
        type=Path,
        help="Parent directory for an auto-created temporary work directory.",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        help="Explicit intermediate-work directory for prepared FASTAs and cached files.",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep an auto-created temporary work directory instead of deleting it.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow overwriting existing outputs and reusing a non-empty --work-dir.",
    )
    parser.add_argument("--title", help="Optional figure title override.")
    parser.add_argument(
        "--reference-role-label",
        default="reference",
        help=(
            "Optional role-label prefix for the top reference row. "
            "Use 'none' to omit the prefix."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument(
        "--check",
        action="store_true",
        help="Check wrapper and external-tool availability, then exit.",
    )
    return parser


def _render_one(
    *,
    subject: Genome,
    comparisons: list[Genome],
    alignments: dict[str, PairwiseAlignment],
    context_genomes: list[Genome],
    ribbon_layers: list,
    ribbon_segments: list | None = None,
    color_intervals: dict | None = None,
    full_color_genomes: set[str] | None = None,
    legend_genome: Genome | None = None,
    output: Path,
    width: float,
    height: float,
    dpi: int,
    title: str | None,
    theme: str,
    actual_reference: Genome | None = None,
    reference_role_label: str | None = "reference",
    color_palette_name: str = "categorical",
) -> dict:
    from scripts.render import render_loom

    return render_loom(
        reference=subject,
        comparisons=comparisons,
        alignments=alignments,
        output=output,
        width=width,
        height=height,
        dpi=dpi,
        title=title,
        theme_name=theme,
        context_genomes=context_genomes,
        ribbon_layers=ribbon_layers,
        ribbon_segments=ribbon_segments,
        color_intervals=color_intervals,
        full_color_genomes=full_color_genomes,
        legend_genome=legend_genome,
        actual_reference=actual_reference,
        reference_role_label=reference_role_label,
        color_palette_name=color_palette_name,
    )


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = _parse_args(parser, argv)

    if args.check:
        return _write_check()
    if args.reference is None:
        parser.error("--reference is required")
    if not args.comparisons:
        parser.error("--comparisons is required")
    if args.outdir is None and args.output is None:
        parser.error(
            "Provide --outdir for figure sets or --output for one overview plot"
        )
    if args.max_contigs == 1 or args.max_contigs < 0:
        parser.error("--max-contigs must be 0 or at least 2")
    if args.reference_segments == 1:
        parser.error("--reference-segments must be 0 or at least 2")
    try:
        display_name_overrides = _parse_display_name_pairs(args.display_names)
    except argparse.ArgumentTypeError as exc:
        parser.error(str(exc))
    try:
        selected_reference_contigs = _parse_contig_name_filters(args.reference_contigs)
    except argparse.ArgumentTypeError as exc:
        parser.error(str(exc))
    reference_role_label = args.reference_role_label.strip()
    if reference_role_label.lower() == "none":
        reference_role_label = ""

    reference_path = args.reference.resolve()
    if not reference_path.exists():
        parser.error(f"--reference does not exist: {reference_path}")

    comparison_paths = _collect_comparisons(args.comparisons)
    missing = [p for p in comparison_paths if not p.exists()]
    if missing:
        parser.error(f"Comparison FASTA does not exist: {missing[0]}")
    if not comparison_paths:
        parser.error("--comparisons did not resolve to any FASTA files")

    outdir = args.outdir or args.output.parent
    summary_path = args.summary_output or (
        outdir / "genome-loom.summary.json"
        if args.outdir
        else args.output.with_suffix(".json")
    )
    output_path = args.output.resolve() if args.output is not None else None
    outdir = outdir.resolve()
    summary_path = summary_path.resolve()

    if (
        args.outdir is not None
        and outdir.exists()
        and any(outdir.iterdir())
        and not args.force
    ):
        parser.error(
            f"--outdir already exists and is not empty: {outdir} (use --force to overwrite)"
        )
    if output_path is not None and output_path.exists() and not args.force:
        parser.error(
            f"--output already exists: {output_path} (use --force to overwrite)"
        )
    if summary_path.exists() and not args.force:
        parser.error(
            f"--summary-output already exists: {summary_path} (use --force to overwrite)"
        )
    outdir.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    tmp_parent = args.tmpdir.expanduser().resolve() if args.tmpdir is not None else None
    if tmp_parent is not None:
        tmp_parent.mkdir(parents=True, exist_ok=True)

    temp_ctx: tempfile.TemporaryDirectory[str] | None = None
    work_dir: Path | None = None
    if args.work_dir is not None:
        work_dir = args.work_dir.expanduser().resolve()
        if work_dir.exists() and any(work_dir.iterdir()) and not args.force:
            parser.error(
                f"--work-dir already exists and is not empty: {work_dir} "
                "(use --force to reuse it)"
            )
        work_dir.mkdir(parents=True, exist_ok=True)
    elif args.keep_temp:
        work_dir = Path(tempfile.mkdtemp(prefix="genome-loom-", dir=tmp_parent))
    else:
        temp_ctx = tempfile.TemporaryDirectory(prefix="genome-loom-", dir=tmp_parent)
        work_dir = Path(temp_ctx.name)

    warnings: list[str] = []

    generated: list[dict] = []
    alignment_cache: dict[tuple[str, str], PairwiseAlignment] = {}
    genomes_by_path: dict[Path, Genome] = {}
    alignment_paths_by_path: dict[Path, Path] = {}

    try:

        def prepare_genome(path: Path, index: int) -> Genome:
            genome = read_fasta(path, min_contig_length=args.min_contig_length)
            if path == reference_path and selected_reference_contigs:
                genome = select_contigs(genome, selected_reference_contigs)
            genome = cap_contig_blocks(genome, max_contigs=args.max_contigs)
            display_name = display_name_overrides.get(
                path.stem
            ) or display_name_overrides.get(path.name)
            if display_name:
                genome = replace(genome, display_name=display_name)
            prepared = work_dir / f"{index:03d}_{_safe_name(genome.name)}.fasta"
            write_fasta(
                prepared,
                [(contig.name, contig.sequence) for contig in genome.contigs],
            )
            alignment_paths_by_path[path] = prepared
            return genome

        reference = prepare_genome(reference_path, 0)
        genomes_by_path[reference_path] = reference
        for index, path in enumerate(comparison_paths, start=1):
            genomes_by_path[path] = prepare_genome(path, index)

        applied_reference_segments = (
            args.reference_segments if args.reference_segments else None
        )

        def get_alignment(
            subject_path: Path, comparison_path: Path
        ) -> PairwiseAlignment:
            key = (str(subject_path), str(comparison_path))
            if key not in alignment_cache:
                aln = run_minimap2(
                    alignment_paths_by_path[subject_path],
                    alignment_paths_by_path[comparison_path],
                    threads=args.threads,
                    preset=args.minimap_preset,
                    min_block_length=args.min_block_length,
                    min_mapq=args.min_mapq,
                )
                alignment_cache[key] = replace(
                    aln,
                    reference=subject_path,
                    comparison=comparison_path,
                )
            return alignment_cache[key]

        reference_similarity_scores: dict[Path, float] = {}
        if args.genome_order == "reference-similarity":
            for path in comparison_paths:
                reference_similarity_scores[path] = get_alignment(
                    reference_path, path
                ).ani
            ordered_comparison_paths = [
                path
                for _index, path in sorted(
                    enumerate(comparison_paths),
                    key=lambda item: (
                        -reference_similarity_scores[item[1]],
                        item[0],
                    ),
                )
            ]
        else:
            ordered_comparison_paths = comparison_paths[:]

        ordered_paths = [reference_path, *ordered_comparison_paths]
        ordered_genomes = [genomes_by_path[p] for p in ordered_paths]
        comparison_genomes = [genomes_by_path[p] for p in ordered_comparison_paths]
        recommended_limit = _recommended_full_stack_limit(args.height)
        if (
            set(args.views) & FULL_STACK_VIEWS
            and len(ordered_genomes) > recommended_limit
        ):
            warnings.append(
                "Full-stack views may look crowded with "
                f"{len(ordered_genomes)} genomes at {args.width:g}x{args.height:g} in. "
                f"Recommended total for overview/neighbor at this height is about "
                f"{recommended_limit} genomes or fewer; consider a taller figure, fewer "
                "genomes, or pairwise-only output."
            )
            print(f"Warning: {warnings[-1]}", file=sys.stderr)

        from scripts.render import RibbonLayer
        from scripts.render import (
            ColorInterval,
            THEMES,
            RibbonSegment,
            color_palette_for_theme,
        )

        def make_layer(subject_path: Path, comparison_path: Path) -> RibbonLayer:
            return RibbonLayer(
                genomes_by_path[subject_path],
                genomes_by_path[comparison_path],
                get_alignment(subject_path, comparison_path),
            )

        def ranked_contigs_by_size(genome: Genome):
            return sorted(
                (
                    item
                    for item in enumerate(genome.contigs)
                    if not item[1].is_remainder
                ),
                key=lambda item: (-item[1].length, item[0]),
            )

        def reference_color_contigs(genome: Genome):
            if args.reference_palette == "continuous":
                return [
                    item
                    for item in enumerate(genome.contigs)
                    if not item[1].is_remainder
                ]
            return ranked_contigs_by_size(genome)

        theme_fallback = THEMES[args.theme].fallback_contig
        reference_palette = color_palette_for_theme(args.reference_palette, args.theme)
        reference_colors = {}
        for rank, (_original_index, contig) in enumerate(
            reference_color_contigs(reference)
        ):
            reference_colors[contig.name] = (
                reference_palette[rank]
                if rank < len(reference_palette)
                else theme_fallback
            )
        for contig in reference.contigs:
            if contig.is_remainder:
                reference_colors[contig.name] = theme_fallback

        def initial_reference_intervals() -> dict[str, dict[str, list[ColorInterval]]]:
            if applied_reference_segments:
                segmented: dict[str, dict[str, list[ColorInterval]]] = {
                    reference.name: {}
                }
                for interval in _reference_segment_intervals(
                    reference=reference,
                    segment_count=applied_reference_segments,
                    colors=reference_palette,
                ):
                    color_interval = ColorInterval(**interval)
                    segmented[reference.name].setdefault(
                        color_interval.contig, []
                    ).append(color_interval)
                return segmented
            return {
                reference.name: {
                    contig.name: [
                        ColorInterval(
                            genome=reference.name,
                            contig=contig.name,
                            start=0,
                            end=contig.length,
                            color=reference_colors[contig.name],
                            origin=contig.name,
                        )
                    ]
                    for contig in reference.contigs
                }
            }

        def dominant_thread(
            intervals: dict[str, dict[str, list[ColorInterval]]],
            genome_name: str,
            contig_name: str,
            start: int,
            end: int,
        ) -> ColorInterval | None:
            candidates = intervals.get(genome_name, {}).get(contig_name, [])
            best: tuple[int, ColorInterval] | None = None
            for interval in candidates:
                overlap = _overlap_len(start, end, interval.start, interval.end)
                if overlap <= 0:
                    continue
                if best is None or overlap > best[0]:
                    best = (overlap, interval)
            return best[1] if best else None

        def thread_slices(
            intervals: dict[str, dict[str, list[ColorInterval]]],
            genome_name: str,
            contig_name: str,
            start: int,
            end: int,
        ) -> list[tuple[int, int, ColorInterval]]:
            candidates = [
                interval
                for interval in intervals.get(genome_name, {}).get(contig_name, [])
                if _overlap_len(start, end, interval.start, interval.end) > 0
            ]
            if not candidates:
                return []
            boundaries = {start, end}
            for interval in candidates:
                boundaries.add(max(start, interval.start))
                boundaries.add(min(end, interval.end))
            ordered = sorted(boundaries)
            slices: list[tuple[int, int, ColorInterval]] = []
            for slice_start, slice_end in zip(ordered, ordered[1:]):
                if slice_end <= slice_start:
                    continue
                thread = dominant_thread(
                    intervals,
                    genome_name,
                    contig_name,
                    slice_start,
                    slice_end,
                )
                if thread is not None:
                    slices.append((slice_start, slice_end, thread))
            return slices

        def projected_comparison_span(block, subject_start: int, subject_end: int):
            subject_len = max(1, block.reference_end - block.reference_start)
            comparison_len = block.comparison_end - block.comparison_start
            if block.strand == "-":
                comparison_start = block.comparison_start + round(
                    (block.reference_end - subject_end)
                    / subject_len
                    * comparison_len
                )
                comparison_end = block.comparison_start + round(
                    (block.reference_end - subject_start)
                    / subject_len
                    * comparison_len
                )
            else:
                comparison_start = block.comparison_start + round(
                    (subject_start - block.reference_start)
                    / subject_len
                    * comparison_len
                )
                comparison_end = block.comparison_start + round(
                    (subject_end - block.reference_start)
                    / subject_len
                    * comparison_len
                )
            return sorted((comparison_start, comparison_end))

        def add_interval(
            intervals: dict[str, dict[str, list[ColorInterval]]],
            genome_name: str,
            contig_name: str,
            start: int,
            end: int,
            color: str,
            origin: str,
        ) -> None:
            if end <= start:
                return
            intervals.setdefault(genome_name, {}).setdefault(contig_name, []).append(
                ColorInterval(
                    genome=genome_name,
                    contig=contig_name,
                    start=start,
                    end=end,
                    color=color,
                    origin=origin,
                )
            )

        def build_reference_flow(
            path_pairs: list[tuple[Path, Path]],
            *,
            base_intervals: dict[str, dict[str, list[ColorInterval]]] | None = None,
        ) -> tuple[list[RibbonSegment], dict[str, dict[str, list[ColorInterval]]]]:
            intervals = (
                {
                    genome_name: {
                        contig_name: list(contig_intervals)
                        for contig_name, contig_intervals in contigs.items()
                    }
                    for genome_name, contigs in base_intervals.items()
                }
                if base_intervals is not None
                else initial_reference_intervals()
            )
            segments: list[RibbonSegment] = []
            for subject_path, comparison_path in path_pairs:
                subject = genomes_by_path[subject_path]
                comparison = genomes_by_path[comparison_path]
                aln = get_alignment(subject_path, comparison_path)
                for block in aln.blocks:
                    slices = thread_slices(
                        intervals,
                        subject.name,
                        block.reference_contig,
                        block.reference_start,
                        block.reference_end,
                    )
                    for subject_start, subject_end, thread in slices:
                        comparison_start, comparison_end = projected_comparison_span(
                            block, subject_start, subject_end
                        )
                        if comparison_end <= comparison_start:
                            continue
                        segments.append(
                            RibbonSegment(
                                subject=subject.name,
                                comparison=comparison.name,
                                subject_contig=block.reference_contig,
                                subject_start=subject_start,
                                subject_end=subject_end,
                                comparison_contig=block.comparison_contig,
                                comparison_start=comparison_start,
                                comparison_end=comparison_end,
                                strand=block.strand,
                                color=thread.color,
                                origin=thread.origin,
                            )
                        )
                        add_interval(
                            intervals,
                            comparison.name,
                            block.comparison_contig,
                            comparison_start,
                            comparison_end,
                            thread.color,
                            thread.origin,
                        )
            return segments, intervals

        _global_reference_flow: (
            tuple[list[RibbonSegment], dict[str, dict[str, list[ColorInterval]]]] | None
        ) = None

        def global_reference_flow() -> (
            tuple[list[RibbonSegment], dict[str, dict[str, list[ColorInterval]]]]
        ):
            nonlocal _global_reference_flow
            if _global_reference_flow is None:
                _global_reference_flow = build_reference_flow(
                    [(reference_path, path) for path in ordered_comparison_paths]
                )
            return _global_reference_flow

        def record(
            view: str, output: Path, meta: dict, subject: Genome, comps: list[Genome]
        ) -> None:
            generated.append(
                {
                    "view": view,
                    "path": str(output),
                    "subject": subject.name,
                    "comparisons": [g.name for g in comps],
                    "render": meta,
                }
            )

        if args.output and not args.outdir:
            segments, intervals = global_reference_flow()
            meta = _render_one(
                subject=reference,
                comparisons=comparison_genomes,
                alignments={},
                context_genomes=ordered_genomes,
                ribbon_layers=[],
                ribbon_segments=segments,
                color_intervals=intervals,
                full_color_genomes={reference.name},
                legend_genome=reference,
                output=args.output,
                width=args.width,
                height=args.height,
                dpi=args.dpi,
                title=args.title or f"{_figure_name(reference)} vs all comparisons",
                theme=args.theme,
                actual_reference=reference,
                reference_role_label=reference_role_label or None,
                color_palette_name=args.reference_palette,
            )
            record("overview", args.output, meta, reference, comparison_genomes)
        else:
            selected = set(args.views)
            ext = args.format

            if "overview" in selected:
                view_dir = outdir / "overview"
                output = view_dir / f"{_safe_name(reference.name)}-vs-all.{ext}"
                segments, intervals = global_reference_flow()
                meta = _render_one(
                    subject=reference,
                    comparisons=comparison_genomes,
                    alignments={},
                    context_genomes=ordered_genomes,
                    ribbon_layers=[],
                    ribbon_segments=segments,
                    color_intervals=intervals,
                    full_color_genomes={reference.name},
                    legend_genome=reference,
                    output=output,
                    width=args.width,
                    height=args.height,
                    dpi=args.dpi,
                    title=args.title
                    or f"{_figure_name(reference)} vs all comparisons",
                    theme=args.theme,
                    actual_reference=reference,
                    reference_role_label=reference_role_label or None,
                    color_palette_name=args.reference_palette,
                )
                record("overview", output, meta, reference, comparison_genomes)

            if "reference-pairs" in selected:
                view_dir = outdir / "reference_pairs"
                for genome, path in zip(comparison_genomes, ordered_comparison_paths):
                    output = (
                        view_dir
                        / f"{_safe_name(reference.name)}-vs-{_safe_name(genome.name)}.{ext}"
                    )
                    segments, _selected_intervals = build_reference_flow(
                        [(reference_path, path)]
                    )
                    _global_segments, intervals = global_reference_flow()
                    meta = _render_one(
                        subject=reference,
                        comparisons=[genome],
                        alignments={},
                        context_genomes=[reference, genome],
                        ribbon_layers=[],
                        ribbon_segments=segments,
                        color_intervals=intervals,
                        full_color_genomes={reference.name},
                        legend_genome=reference,
                        output=output,
                        width=args.width,
                        height=args.height,
                        dpi=args.dpi,
                        title=f"{_figure_name(reference)} vs {_figure_name(genome)}",
                        theme=args.theme,
                        actual_reference=reference,
                        reference_role_label=reference_role_label or None,
                        color_palette_name=args.reference_palette,
                    )
                    record("reference-pairs", output, meta, reference, [genome])

            if "all-pairs" in selected:
                view_dir = outdir / "all_pairs"
                for a_path, b_path in itertools.combinations(ordered_paths, 2):
                    subject = genomes_by_path[a_path]
                    comp = genomes_by_path[b_path]
                    output = (
                        view_dir
                        / f"{_safe_name(subject.name)}-vs-{_safe_name(comp.name)}.{ext}"
                    )
                    if a_path == reference_path:
                        segments, _selected_intervals = build_reference_flow(
                            [(a_path, b_path)]
                        )
                        _global_segments, intervals = global_reference_flow()
                        meta = _render_one(
                            subject=subject,
                            comparisons=[comp],
                            alignments={},
                            context_genomes=[subject, comp],
                            ribbon_layers=[],
                            ribbon_segments=segments,
                            color_intervals=intervals,
                            full_color_genomes={reference.name},
                            legend_genome=reference,
                            output=output,
                            width=args.width,
                            height=args.height,
                            dpi=args.dpi,
                            title=f"{_figure_name(subject)} vs {_figure_name(comp)}",
                            theme=args.theme,
                            actual_reference=reference,
                            reference_role_label=reference_role_label or None,
                            color_palette_name=args.reference_palette,
                        )
                    else:
                        layers = [make_layer(a_path, b_path)]
                        _global_segments, intervals = global_reference_flow()
                        meta = _render_one(
                            subject=subject,
                            comparisons=[comp],
                            alignments={},
                            context_genomes=[subject, comp],
                            ribbon_layers=layers,
                            color_intervals=intervals,
                            full_color_genomes={reference.name},
                            legend_genome=reference,
                            output=output,
                            width=args.width,
                            height=args.height,
                            dpi=args.dpi,
                            title=f"{_figure_name(subject)} vs {_figure_name(comp)}",
                            theme=args.theme,
                            actual_reference=reference,
                            reference_role_label=reference_role_label or None,
                            color_palette_name=args.reference_palette,
                        )
                    record("all-pairs", output, meta, subject, [comp])

            if "neighbor" in selected:
                view_dir = outdir / "neighbor"
                output = view_dir / f"neighbor-chain.{ext}"
                neighbor_segments, intervals = build_reference_flow(
                    list(zip(ordered_paths, ordered_paths[1:]))
                )
                meta = _render_one(
                    subject=reference,
                    comparisons=comparison_genomes,
                    alignments={},
                    context_genomes=ordered_genomes,
                    ribbon_layers=[],
                    ribbon_segments=neighbor_segments,
                    color_intervals=intervals,
                    full_color_genomes={reference.name},
                    legend_genome=reference,
                    output=output,
                    width=args.width,
                    height=args.height,
                    dpi=args.dpi,
                    title="Neighbor chain",
                    theme=args.theme,
                    actual_reference=reference,
                    reference_role_label=reference_role_label or None,
                    color_palette_name=args.reference_palette,
                )
                record("neighbor", output, meta, reference, ordered_genomes[1:])

        write_summary(
            summary_path,
            {
                "version": VERSION,
                "inputs": {
                    "reference": str(reference_path),
                    "comparisons": [str(p) for p in comparison_paths],
                    "reference_contigs": selected_reference_contigs or None,
                    "reference_segments": args.reference_segments or None,
                    "reference_palette": args.reference_palette,
                },
                "outputs": {
                    "outdir": str(outdir),
                    "output": str(output_path) if output_path is not None else None,
                    "summary_json": str(summary_path),
                    "work_dir": (
                        str(work_dir) if args.keep_temp or args.work_dir else None
                    ),
                    "figures": generated,
                },
                "settings": {
                    "views": args.views,
                    "format": args.format,
                    "theme": args.theme,
                    "genome_order": args.genome_order,
                    "reference_palette": args.reference_palette,
                    "reference_role_label": reference_role_label or None,
                    "reference_segments": applied_reference_segments,
                    "width": args.width,
                    "height": args.height,
                    "dpi": args.dpi,
                },
                "filters": {
                    "reference_contigs": selected_reference_contigs or None,
                    "reference_segments": applied_reference_segments,
                    "min_contig_length": args.min_contig_length,
                    "max_contigs": args.max_contigs,
                    "min_block_length": args.min_block_length,
                    "min_mapq": args.min_mapq,
                },
                "runtime": {
                    "threads": args.threads,
                    "config": (
                        str(args.config.resolve()) if args.config is not None else None
                    ),
                    "tmpdir": str(tmp_parent) if tmp_parent is not None else None,
                    "keep_temp": args.keep_temp,
                    "force": args.force,
                },
                "alignment": {
                    "tool": "minimap2",
                    "preset": args.minimap_preset,
                },
                "coloring_explanation": {
                    "comparison_contigs": (
                        "Comparison contigs are painted from their relationship to the true project reference."
                    ),
                    "ribbons": (
                        "When the true reference is the upper genome in a join, ribbons use the reference contig color in that aligned region. "
                        "When a non-reference genome is the upper genome, ribbons use that genome's reference-based painted color in the aligned region, "
                        "and fall back to the neutral comparison color when no reference-based paint is present there."
                    ),
                },
                "ordering": {
                    "method": args.genome_order,
                    "input_comparison_order": [
                        genomes_by_path[p].name for p in comparison_paths
                    ],
                    "comparison_order": [
                        genomes_by_path[p].name for p in ordered_comparison_paths
                    ],
                    "reference_similarity": [
                        {
                            "name": genomes_by_path[p].name,
                            "ani": round(reference_similarity_scores[p], 6),
                        }
                        for p in ordered_comparison_paths
                        if p in reference_similarity_scores
                    ],
                    "genome_order": [g.name for g in ordered_genomes],
                },
                "genomes": [_summary_for_genome(g) for g in ordered_genomes],
                "alignments": {
                    f"{Path(subject).name}::{Path(comp).name}": _alignment_summary(aln)
                    for (subject, comp), aln in alignment_cache.items()
                },
                "warnings": warnings,
            },
        )
    except Exception as exc:
        try:
            _write_failure_summary(
                summary_path,
                version=VERSION,
                outdir=outdir,
                output=output_path,
                work_dir=work_dir,
                exc=exc,
            )
        except Exception:
            pass
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    finally:
        if temp_ctx is not None:
            temp_ctx.cleanup()

    print(f"Generated {len(generated)} figure(s) in {outdir}", file=sys.stderr)
    print(f"Summary written: {summary_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
