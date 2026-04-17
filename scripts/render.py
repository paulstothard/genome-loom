from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
from matplotlib.path import Path as MplPath

from .align import PairwiseAlignment
from .fasta import Genome


DEFAULT_COLORS = [
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#C9B000",
    "#56B4E9",
    "#E69F00",
    "#6A3D9A",
    "#A6761D",
    "#1B9E77",
    "#E7298A",
    "#66A61E",
    "#7570B3",
    "#A6CEE3",
    "#B15928",
    "#FB9A99",
    "#8DD3C7",
]
FALLBACK_CONTIG_COLOR = "#b8b8b8"
ENABLE_RIBBON_OUTLINES = False


@dataclass(frozen=True)
class Theme:
    name: str
    background: str
    text: str
    comparison_fill: str
    comparison_edge: str
    subject_edge: str
    fallback_contig: str
    scale: str
    legend_edge: str
    ribbon_alpha_forward: float
    ribbon_alpha_reverse: float


THEMES = {
    "light": Theme(
        name="light",
        background="#ffffff",
        text="#111111",
        comparison_fill="#c9c9c9",
        comparison_edge="#777777",
        subject_edge="#555555",
        fallback_contig=FALLBACK_CONTIG_COLOR,
        scale="#000000",
        legend_edge="#999999",
        ribbon_alpha_forward=0.22,
        ribbon_alpha_reverse=0.16,
    ),
    "dark": Theme(
        name="dark",
        background="#101214",
        text="#f1efe8",
        comparison_fill="#4d5458",
        comparison_edge="#9aa0a3",
        subject_edge="#d8d2c4",
        fallback_contig="#72777b",
        scale="#f1efe8",
        legend_edge="#8c9296",
        ribbon_alpha_forward=0.30,
        ribbon_alpha_reverse=0.22,
    ),
}


@dataclass(frozen=True)
class ContigLayout:
    x0: float
    x1: float
    length: int


@dataclass(frozen=True)
class RibbonLayer:
    subject: Genome
    comparison: Genome
    alignment: PairwiseAlignment


@dataclass(frozen=True)
class ColorInterval:
    genome: str
    contig: str
    start: int
    end: int
    color: str
    origin: str


@dataclass(frozen=True)
class RibbonSegment:
    subject: str
    comparison: str
    subject_contig: str
    subject_start: int
    subject_end: int
    comparison_contig: str
    comparison_start: int
    comparison_end: int
    strand: str
    color: str
    origin: str


def _nice_scale_length(total_bp: int) -> int:
    target = total_bp / 5.0
    if target <= 0:
        return 1
    scale = 10 ** math.floor(math.log10(target))
    for m in [1, 2, 5, 10]:
        candidate = int(m * scale)
        if candidate >= target:
            return candidate
    return int(10 * scale)


def _format_bp(bp: int) -> str:
    if bp < 10_000:
        return f"{bp:,} bp"
    if bp < 1_000_000:
        return f"{bp / 1_000:g} kb"
    return f"{bp / 1_000_000:g} Mb"


def _shorten(text: str, max_chars: int) -> str:
    if len(text) <= max_chars:
        return text
    return text[: max(1, max_chars - 1)] + "…"


def _display_name(genome: Genome) -> str:
    return (genome.display_name or genome.name).replace("_", " ")


def _row_label(
    genome: Genome,
    reference: Genome,
    reference_role_label: str | None = "reference",
) -> str:
    label = _display_name(genome)
    if genome.name == reference.name and reference_role_label:
        return f"{reference_role_label} | {label}"
    return label


def _text_width_frac(text: str, fontsize_pts: float, width_in: float) -> float:
    """Estimate text width in axis-fraction units.

    This is intentionally conservative; it lets layout avoid obvious collisions
    without needing a renderer pass just to measure text extents.
    """
    return len(text) * fontsize_pts * 0.56 / 72.0 / width_in


def _wrap_legend_entries(
    entries: list[tuple[str, str, str]],
    *,
    plot_width: float,
    width_in: float,
    fontsize_pts: float,
    sw: float,
) -> tuple[list[list[tuple[str, str, str, float]]], float]:
    items: list[tuple[str, str, str, float]] = []
    item_gap = 0.026
    for raw_label, display_label, color in entries:
        item_w = sw + 0.008 + _text_width_frac(display_label, fontsize_pts, width_in)
        item_w = min(item_w, plot_width * 0.70)
        items.append((raw_label, display_label, color, item_w))

    if not items:
        return [[]], plot_width

    cell_w = max(item_w for *_rest, item_w in items)
    cols = max(1, int((plot_width + item_gap) / (cell_w + item_gap)))
    cols = min(cols, len(items))
    rows = [items[i : i + cols] for i in range(0, len(items), cols)]
    return rows, cell_w


def _ranked_contigs_by_size(genome: Genome):
    return sorted(
        (item for item in enumerate(genome.contigs) if not item[1].is_remainder),
        key=lambda item: (-item[1].length, item[0]),
    )


def _contig_colors_by_size(genome: Genome, fallback_color: str) -> dict[str, str]:
    colors: dict[str, str] = {}
    for rank, (_original_index, contig) in enumerate(_ranked_contigs_by_size(genome)):
        colors[contig.name] = (
            DEFAULT_COLORS[rank] if rank < len(DEFAULT_COLORS) else fallback_color
        )
    for contig in genome.contigs:
        if contig.is_remainder:
            colors[contig.name] = fallback_color
    return colors


def _has_intervals(
    intervals: dict[str, dict[str, list[ColorInterval]]], genome_name: str
) -> bool:
    return any(intervals.get(genome_name, {}).values())


def _dominant_interval_color(
    intervals: dict[str, dict[str, list[ColorInterval]]],
    genome_name: str,
    contig_name: str,
    start: int,
    end: int,
) -> str | None:
    candidates = intervals.get(genome_name, {}).get(contig_name, [])
    best: tuple[int, str] | None = None
    for interval in candidates:
        overlap = max(0, min(end, interval.end) - max(start, interval.start))
        if overlap <= 0:
            continue
        if best is None or overlap > best[0]:
            best = (overlap, interval.color)
    return best[1] if best else None


def _mix_colors(color_a: str, color_b: str, frac_b: float) -> str:
    frac_b = min(max(frac_b, 0.0), 1.0)
    frac_a = 1.0 - frac_b
    a = mcolors.to_rgb(color_a)
    b = mcolors.to_rgb(color_b)
    return mcolors.to_hex(
        tuple(frac_a * av + frac_b * bv for av, bv in zip(a, b)),
        keep_alpha=False,
    )


def _ribbon_sort_key(
    top_left: float, top_right: float, bottom_left: float, bottom_right: float
):
    top_span = abs(top_right - top_left)
    bottom_span = abs(bottom_right - bottom_left)
    avg_span = (top_span + bottom_span) / 2
    center_shift = abs(
        ((top_left + top_right) / 2) - ((bottom_left + bottom_right) / 2)
    )
    return (
        -avg_span,
        -center_shift,
        min(top_left, top_right, bottom_left, bottom_right),
    )


def _layout_genome(
    genome: Genome,
    x_left: float,
    x_width: float,
    max_span: int,
    gap_bp: int,
) -> dict[str, ContigLayout]:
    x_per_bp = x_width / max(max_span, 1)
    offset = 0
    layout: dict[str, ContigLayout] = {}
    for contig in genome.contigs:
        x0 = x_left + offset * x_per_bp
        x1 = x0 + contig.length * x_per_bp
        layout[contig.name] = ContigLayout(x0=x0, x1=x1, length=contig.length)
        offset += contig.length + gap_bp
    return layout


def _pos_to_x(layout: ContigLayout, pos: int) -> float:
    frac = min(max(pos / max(layout.length, 1), 0.0), 1.0)
    return layout.x0 + frac * (layout.x1 - layout.x0)


def _ribbon_path(
    x1: float,
    x2: float,
    y_top: float,
    x3: float,
    x4: float,
    y_bottom: float,
    *,
    reverse: bool = False,
) -> MplPath:
    top_left, top_right = sorted((x1, x2))
    bottom_left, bottom_right = sorted((x3, x4))
    if not reverse:
        verts = [
            (top_left, y_top),
            (bottom_left, y_bottom),
            (bottom_right, y_bottom),
            (top_right, y_top),
            (top_left, y_top),
        ]
        codes = [
            MplPath.MOVETO,
            MplPath.LINETO,
            MplPath.LINETO,
            MplPath.LINETO,
            MplPath.CLOSEPOLY,
        ]
        return MplPath(verts, codes)

    y_mid = (y_top + y_bottom) / 2
    top_center = (top_left + top_right) / 2
    bottom_center = (bottom_left + bottom_right) / 2
    pinch_x = (top_center + bottom_center) / 2
    min_width = min(top_right - top_left, bottom_right - bottom_left)
    waist_half = max(0.0004, min(0.006, min_width * 0.06))
    mid_left = pinch_x - waist_half
    mid_right = pinch_x + waist_half
    verts = [
        (top_left, y_top),
        (mid_left, y_mid),
        (bottom_right, y_bottom),
        (bottom_left, y_bottom),
        (mid_right, y_mid),
        (top_right, y_top),
        (top_left, y_top),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.LINETO,
        MplPath.LINETO,
        MplPath.LINETO,
        MplPath.LINETO,
        MplPath.LINETO,
        MplPath.CLOSEPOLY,
    ]
    return MplPath(verts, codes)


def render_loom(
    *,
    reference: Genome,
    comparisons: list[Genome],
    alignments: dict[str, PairwiseAlignment],
    output: Path,
    width: float,
    height: float,
    dpi: int,
    title: str | None = None,
    theme_name: str = "light",
    context_genomes: list[Genome] | None = None,
    ribbon_layers: list[RibbonLayer] | None = None,
    ribbon_segments: list[RibbonSegment] | None = None,
    color_intervals: dict[str, dict[str, list[ColorInterval]]] | None = None,
    full_color_genomes: set[str] | None = None,
    legend_genome: Genome | None = None,
    actual_reference: Genome | None = None,
    reference_role_label: str | None = "reference",
) -> dict:
    """Render a stacked subject-to-comparison ribbon plot."""
    theme = THEMES.get(theme_name)
    if theme is None:
        raise ValueError(
            f"Unknown theme '{theme_name}'. Choose from: {', '.join(THEMES)}"
        )

    genomes = (
        context_genomes if context_genomes is not None else [reference, *comparisons]
    )
    row_count = len(genomes)
    if row_count < 2:
        raise ValueError("At least one comparison genome is required")
    if ribbon_layers is None:
        ribbon_layers = [
            RibbonLayer(reference, comp, alignments[comp.name])
            for comp in comparisons
            if comp.name in alignments
        ]
    if full_color_genomes is None:
        full_color_genomes = {layer.subject.name for layer in ribbon_layers}
    if color_intervals is None:
        color_intervals = {}
    if actual_reference is None:
        actual_reference = reference

    max_len = max(g.length for g in genomes)
    gap_bp = max(1_000, int(max_len * 0.008))
    max_span = max(
        sum(contig.length for contig in genome.contigs)
        + gap_bp * max(0, len(genome.contigs) - 1)
        for genome in genomes
    )

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
    fig.patch.set_facecolor(theme.background)
    ax.set_facecolor(theme.background)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    row_labels = {
        genome.name: _row_label(
            genome,
            actual_reference,
            reference_role_label=reference_role_label,
        )
        for genome in genomes
    }
    longest_genome_label = max(len(label) for label in row_labels.values())
    label_size = max(7.0, min(16.0, height * 0.95, 82.0 / row_count))
    font_size = max(6.5, min(13.0, height * 0.82))
    left_label = max(
        0.16,
        min(
            0.30,
            0.035 + _text_width_frac("M" * longest_genome_label, label_size, width),
        ),
    )
    right_pad = 0.05
    plot_left = left_label
    plot_right = 1.0 - right_pad
    plot_width = plot_right - plot_left

    subjects = []
    seen_subjects: set[str] = set()
    for layer in ribbon_layers:
        if layer.subject.name not in seen_subjects:
            subjects.append(layer.subject)
            seen_subjects.add(layer.subject.name)
    for genome in genomes:
        if genome.name in full_color_genomes and genome.name not in seen_subjects:
            subjects.append(genome)
            seen_subjects.add(genome.name)
    if not subjects:
        subjects = [reference]

    subject_colors: dict[str, dict[str, str]] = {}
    for subject in subjects:
        subject_colors[subject.name] = _contig_colors_by_size(
            subject, theme.fallback_contig
        )

    legend_subject = legend_genome or subjects[0]
    if legend_genome is None and color_intervals:
        for subject in subjects:
            if subject.name in full_color_genomes:
                legend_subject = subject
                break
    legend_colors = subject_colors.get(legend_subject.name)
    if legend_colors is None:
        legend_colors = _contig_colors_by_size(legend_subject, theme.fallback_contig)
    ranked_legend_contigs = _ranked_contigs_by_size(legend_subject)
    distinct_color_count = min(len(ranked_legend_contigs), len(DEFAULT_COLORS))
    legend_entries = [
        (contig.name, _shorten(contig.name, 24), legend_colors[contig.name])
        for _original_index, contig in ranked_legend_contigs[:distinct_color_count]
    ]
    overflow_count = max(0, len(ranked_legend_contigs) - distinct_color_count)
    remainder_count = sum(
        contig.source_count for contig in legend_subject.contigs if contig.is_remainder
    )
    fallback_count = overflow_count + remainder_count
    if fallback_count:
        legend_entries.append(
            (
                "other_contigs",
                f"other contigs ({fallback_count})",
                theme.fallback_contig,
            )
        )

    sw = max(0.010, min(0.018, 0.18 / width))
    legend_rows, legend_cell_w = _wrap_legend_entries(
        legend_entries,
        plot_width=plot_width,
        width_in=width,
        fontsize_pts=font_size,
        sw=sw,
    )
    legend_line_h = max(sw * 1.55, (font_size / 72.0 / height) * 1.55)
    legend_height = len(legend_rows) * legend_line_h
    bottom_pad = 0.035
    legend_y0 = bottom_pad
    scale_label_h = (font_size / 72.0 / height) * 1.25
    scale_band_h = scale_label_h + 0.045
    note_text = (
        "Reference-based colors propagated from "
        f"{_shorten(_display_name(actual_reference), 28)} to comparison contigs using sequence alignments."
    )
    note_h = (font_size / 72.0 / height) * 1.25 if color_intervals else 0.0
    note_gap = 0.010 if color_intervals else 0.0
    note_y = legend_y0 + legend_height + note_gap
    scale_y = note_y + note_h + scale_band_h * 0.72
    bottom = min(0.42, scale_y + 0.055)
    top = 0.88 if title else 0.93
    if bottom >= top - 0.18:
        bottom = top - 0.18
    row_area = top - bottom
    row_gap = row_area / max(row_count - 1, 1)
    block_h = max(0.006, min(0.026, row_gap * 0.20))

    if title:
        ax.text(
            0.5,
            0.965,
            title,
            ha="center",
            va="top",
            fontsize=max(12.0, height * 1.15),
            fontweight="bold",
            color=theme.text,
        )

    layouts = {
        genome.name: _layout_genome(genome, plot_left, plot_width, max_span, gap_bp)
        for genome in genomes
    }
    row_y = {genome.name: top - i * row_gap for i, genome in enumerate(genomes)}
    scale_bp = _nice_scale_length(max_len)
    scale_w = (scale_bp / max(max_span, 1)) * plot_width
    ribbon_edge_forward = mcolors.to_rgba(
        _mix_colors(theme.text, theme.background, 0.40), alpha=0.55
    )
    ribbon_edge_reverse = mcolors.to_rgba(
        _mix_colors("#8b1e3f", theme.background, 0.20), alpha=0.70
    )

    # Ribbons first, so genome blocks and colored landing zones sit above them.
    if ribbon_segments is not None:
        ribbon_items = []
        for segment in ribbon_segments:
            ref_layout = layouts[segment.subject].get(segment.subject_contig)
            comp_layout = layouts[segment.comparison].get(segment.comparison_contig)
            if ref_layout is None or comp_layout is None:
                continue
            rx1 = _pos_to_x(ref_layout, segment.subject_start)
            rx2 = _pos_to_x(ref_layout, segment.subject_end)
            cx1 = _pos_to_x(comp_layout, segment.comparison_start)
            cx2 = _pos_to_x(comp_layout, segment.comparison_end)
            if segment.strand == "-":
                cx1, cx2 = cx2, cx1
            ribbon_items.append(
                (
                    _ribbon_sort_key(rx1, rx2, cx1, cx2),
                    segment,
                    rx1,
                    rx2,
                    cx1,
                    cx2,
                )
            )
        for _sort_key, segment, rx1, rx2, cx1, cx2 in sorted(
            ribbon_items, key=lambda item: item[0]
        ):
            path = _ribbon_path(
                rx1,
                rx2,
                row_y[segment.subject] - block_h / 2,
                cx1,
                cx2,
                row_y[segment.comparison] + block_h / 2,
                reverse=segment.strand == "-",
            )
            ax.add_patch(
                patches.PathPatch(
                    path,
                    facecolor=mcolors.to_rgba(
                        segment.color,
                        alpha=theme.ribbon_alpha_forward,
                    ),
                    edgecolor=(
                        (
                            ribbon_edge_forward
                            if segment.strand == "+"
                            else ribbon_edge_reverse
                        )
                        if ENABLE_RIBBON_OUTLINES
                        else "none"
                    ),
                    linewidth=0.72 if ENABLE_RIBBON_OUTLINES else 0.0,
                    zorder=1,
                )
            )
    else:
        for layer in ribbon_layers:
            subject = layer.subject
            comp = layer.comparison
            colors = subject_colors.get(subject.name, {})
            ribbon_items = []
            for block in layer.alignment.blocks:
                ref_layout = layouts[subject.name].get(block.reference_contig)
                comp_layout = layouts[comp.name].get(block.comparison_contig)
                if ref_layout is None or comp_layout is None:
                    continue
                rx1 = _pos_to_x(ref_layout, block.reference_start)
                rx2 = _pos_to_x(ref_layout, block.reference_end)
                cx1 = _pos_to_x(comp_layout, block.comparison_start)
                cx2 = _pos_to_x(comp_layout, block.comparison_end)
                if block.strand == "-":
                    cx1, cx2 = cx2, cx1
                color = _dominant_interval_color(
                    color_intervals,
                    subject.name,
                    block.reference_contig,
                    block.reference_start,
                    block.reference_end,
                )
                if color is None:
                    if subject.name in full_color_genomes:
                        color = colors.get(
                            block.reference_contig, theme.fallback_contig
                        )
                    else:
                        color = theme.comparison_fill
                ribbon_items.append(
                    (
                        _ribbon_sort_key(rx1, rx2, cx1, cx2),
                        block,
                        color,
                        rx1,
                        rx2,
                        cx1,
                        cx2,
                    )
                )
            for _sort_key, block, color, rx1, rx2, cx1, cx2 in sorted(
                ribbon_items, key=lambda item: item[0]
            ):
                path = _ribbon_path(
                    rx1,
                    rx2,
                    row_y[subject.name] - block_h / 2,
                    cx1,
                    cx2,
                    row_y[comp.name] + block_h / 2,
                    reverse=block.strand == "-",
                )
                ax.add_patch(
                    patches.PathPatch(
                        path,
                        facecolor=mcolors.to_rgba(
                            color,
                            alpha=theme.ribbon_alpha_forward,
                        ),
                        edgecolor=(
                            (
                                ribbon_edge_forward
                                if block.strand == "+"
                                else ribbon_edge_reverse
                            )
                            if ENABLE_RIBBON_OUTLINES
                            else "none"
                        ),
                        linewidth=0.72 if ENABLE_RIBBON_OUTLINES else 0.0,
                        zorder=1,
                    )
                )

    # Genome rows.
    for i, genome in enumerate(genomes):
        y = row_y[genome.name]
        row_label = row_labels[genome.name]
        max_label_chars = max(
            14, int((left_label - 0.04) * width * 72 / (label_size * 0.56))
        )
        ax.text(
            plot_left - 0.018,
            y,
            _shorten(row_label, max_label_chars),
            ha="right",
            va="center",
            fontsize=label_size,
            fontstyle="italic" if i > 0 else "normal",
            color=theme.text,
        )
        for contig in genome.contigs:
            layout = layouts[genome.name][contig.name]
            genome_colors = subject_colors.get(genome.name)
            is_full_color = (
                genome.name in full_color_genomes and genome_colors is not None
            )
            has_propagated_color = _has_intervals(color_intervals, genome.name)
            color = (
                genome_colors.get(contig.name, theme.fallback_contig)
                if is_full_color
                else theme.comparison_fill
            )
            edge = (
                theme.subject_edge
                if is_full_color or has_propagated_color
                else theme.comparison_edge
            )
            ax.add_patch(
                patches.Rectangle(
                    (layout.x0, y - block_h / 2),
                    layout.x1 - layout.x0,
                    block_h,
                    facecolor=color,
                    edgecolor=edge,
                    linewidth=0.45,
                    zorder=3,
                )
            )
            for interval in color_intervals.get(genome.name, {}).get(contig.name, []):
                x0 = _pos_to_x(layout, interval.start)
                x1 = _pos_to_x(layout, interval.end)
                if x1 <= x0:
                    continue
                ax.add_patch(
                    patches.Rectangle(
                        (x0, y - block_h / 2),
                        x1 - x0,
                        block_h,
                        facecolor=interval.color,
                        edgecolor="none",
                        alpha=0.94 if not is_full_color else 0.30,
                        zorder=4,
                    )
                )

    # Scale bar.
    sx0 = plot_left
    sy = scale_y
    ax.plot([sx0, sx0 + scale_w], [sy, sy], color=theme.scale, linewidth=2.0)
    ax.text(
        sx0 + scale_w / 2,
        sy + 0.018,
        _format_bp(scale_bp),
        ha="center",
        va="bottom",
        fontsize=font_size,
        color=theme.text,
    )

    # Contig legend.
    for row_idx, row in enumerate(legend_rows):
        x = plot_left
        y = legend_y0 + (len(legend_rows) - 1 - row_idx) * legend_line_h
        for _raw_label, display_label, color, item_w in row:
            ax.add_patch(
                patches.Rectangle(
                    (x, y),
                    sw,
                    sw,
                    facecolor=color,
                    edgecolor=theme.legend_edge,
                    linewidth=0.4,
                )
            )
            ax.text(
                x + sw + 0.006,
                y + sw / 2,
                display_label,
                ha="left",
                va="center",
                fontsize=font_size,
                color=theme.text,
            )
            x += legend_cell_w + 0.026

    if color_intervals:
        ax.text(
            plot_left,
            note_y,
            note_text,
            ha="left",
            va="bottom",
            fontsize=max(6.0, font_size * 0.86),
            color=theme.text,
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return {
        "subject_contig_colors": subject_colors,
        "legend_contig_colors": legend_colors,
        "legend_subject": legend_subject.name,
        "legend_note": note_text if color_intervals else None,
        "colored_subjects": [s.name for s in subjects],
        "distinct_subject_contig_colors": distinct_color_count,
        "fallback_subject_contigs": fallback_count,
        "legend_rows": len(legend_rows),
        "scale_bar_bp": scale_bp,
        "row_order": [g.name for g in genomes],
        "theme": theme.name,
        "thread_mode": (
            "reference-flow" if ribbon_segments is not None else "subject-local"
        ),
        "color_assignment": "largest-contigs-first",
        "comparison_contig_coloring": (
            "reference-based" if color_intervals else "neutral"
        ),
        "coloring_explanation": {
            "comparison_contigs": (
                "Comparison contigs are painted from their relationship to the true project reference."
                if color_intervals
                else "Comparison contigs use the neutral comparison color."
            ),
            "ribbons": (
                "Ribbons use the upper genome's reference-based painted color in the aligned region, and fall back to the neutral comparison color when that upper region has no reference-based paint."
                if ribbon_segments is None
                else "Ribbons use reference-flow coloring propagated from the true project reference."
            ),
        },
    }
