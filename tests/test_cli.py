from __future__ import annotations

import json
from pathlib import Path

from genome_loom import _parse_args, build_parser, main
from scripts.align import AlignmentBlock, PairwiseAlignment


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{sequence}\n")


def test_reference_contigs_can_be_loaded_from_config(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fasta"
    comparison = tmp_path / "comparison.fasta"
    write_fasta(reference, [("chr", "ACGT"), ("plasmid", "TGCA")])
    write_fasta(comparison, [("comp1", "ACGT")])

    config = tmp_path / "run.json"
    config.write_text(
        json.dumps(
            {
                "reference": str(reference),
                "reference_contigs": ["plasmid"],
                "comparisons": [str(comparison)],
                "outdir": str(tmp_path / "results"),
                "views": ["overview"],
            }
        ),
        encoding="utf-8",
    )

    args = _parse_args(build_parser(), ["--config", str(config)])
    assert args.reference == reference
    assert args.reference_contigs == ["plasmid"]


def test_reference_segments_can_be_loaded_from_config(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fasta"
    comparison = tmp_path / "comparison.fasta"
    write_fasta(reference, [("chr", "ACGT")])
    write_fasta(comparison, [("comp1", "ACGT")])

    config = tmp_path / "run.json"
    config.write_text(
        json.dumps(
            {
                "reference": str(reference),
                "reference_segments": 4,
                "comparisons": [str(comparison)],
                "outdir": str(tmp_path / "results"),
                "views": ["overview"],
            }
        ),
        encoding="utf-8",
    )

    args = _parse_args(build_parser(), ["--config", str(config)])
    assert args.reference == reference
    assert args.reference_segments == 4


def test_reference_segments_are_passed_to_renderer(
    monkeypatch, tmp_path: Path
) -> None:
    reference = tmp_path / "reference.fasta"
    comparison = tmp_path / "comparison.fasta"
    write_fasta(reference, [("chr", "ACGTACGT")])
    write_fasta(comparison, [("comp1", "ACGTACGT")])

    monkeypatch.setattr(
        "genome_loom.run_minimap2",
        lambda reference_fasta, comparison_fasta, **kwargs: PairwiseAlignment(
            reference=reference_fasta,
            comparison=comparison_fasta,
            blocks=[],
            aligned_bases=8,
            matching_bases=8,
            ani=99.0,
        ),
    )

    captured = {}

    def fake_render(**kwargs):
        captured["color_intervals"] = kwargs["color_intervals"]
        return {"figure_path": str(kwargs["output"]), "view": "overview"}

    monkeypatch.setattr("scripts.render.render_loom", fake_render)

    outdir = tmp_path / "results"
    exit_code = main(
        [
            "--reference",
            str(reference),
            "--reference-segments",
            "4",
            "--comparisons",
            str(comparison),
            "--outdir",
            str(outdir),
            "--views",
            "overview",
            "--min-contig-length",
            "0",
        ]
    )
    assert exit_code == 0

    intervals = captured["color_intervals"]["reference"]["chr"]
    assert len(intervals) == 4
    assert [interval.origin for interval in intervals] == [
        "segment 1",
        "segment 2",
        "segment 3",
        "segment 4",
    ]

    summary = json.loads((outdir / "genome-loom.summary.json").read_text())
    assert summary["inputs"]["reference_segments"] == 4
    assert summary["settings"]["reference_segments"] == 4
    assert summary["filters"]["reference_segments"] == 4


def test_reference_segments_can_span_multiple_contigs(
    monkeypatch, tmp_path: Path
) -> None:
    reference = tmp_path / "reference.fasta"
    comparison = tmp_path / "comparison.fasta"
    write_fasta(reference, [("chrA", "A" * 6), ("chrB", "C" * 4)])
    write_fasta(comparison, [("comp1", "A" * 10)])

    monkeypatch.setattr(
        "genome_loom.run_minimap2",
        lambda reference_fasta, comparison_fasta, **kwargs: PairwiseAlignment(
            reference=reference_fasta,
            comparison=comparison_fasta,
            blocks=[],
            aligned_bases=0,
            matching_bases=0,
            ani=0.0,
        ),
    )

    captured = {}

    def fake_render(**kwargs):
        captured["color_intervals"] = kwargs["color_intervals"]
        return {"figure_path": str(kwargs["output"]), "view": "overview"}

    monkeypatch.setattr("scripts.render.render_loom", fake_render)

    outdir = tmp_path / "results"
    exit_code = main(
        [
            "--reference",
            str(reference),
            "--reference-segments",
            "4",
            "--comparisons",
            str(comparison),
            "--outdir",
            str(outdir),
            "--views",
            "overview",
            "--min-contig-length",
            "0",
        ]
    )
    assert exit_code == 0

    intervals = captured["color_intervals"]["reference"]
    assert [(i.start, i.end, i.origin) for i in intervals["chrA"]] == [
        (0, 2, "segment 1"),
        (2, 5, "segment 2"),
        (5, 6, "segment 3"),
    ]
    assert [(i.start, i.end, i.origin) for i in intervals["chrB"]] == [
        (0, 2, "segment 3"),
        (2, 4, "segment 4"),
    ]

    summary = json.loads((outdir / "genome-loom.summary.json").read_text())
    assert summary["genomes"][0]["contig_count"] == 2
    assert summary["settings"]["reference_segments"] == 4


def test_reference_segments_split_long_alignment_blocks(
    monkeypatch, tmp_path: Path
) -> None:
    reference = tmp_path / "reference.fasta"
    comparison = tmp_path / "comparison.fasta"
    write_fasta(reference, [("chr", "ACGTACGT")])
    write_fasta(comparison, [("comp1", "ACGTACGT")])

    monkeypatch.setattr(
        "genome_loom.run_minimap2",
        lambda reference_fasta, comparison_fasta, **kwargs: PairwiseAlignment(
            reference=reference_fasta,
            comparison=comparison_fasta,
            blocks=[
                AlignmentBlock(
                    reference_contig="chr",
                    reference_start=0,
                    reference_end=8,
                    comparison_contig="comp1",
                    comparison_start=0,
                    comparison_end=8,
                    strand="+",
                    matches=8,
                    block_length=8,
                    mapq=60,
                )
            ],
            aligned_bases=8,
            matching_bases=8,
            ani=99.0,
        ),
    )

    captured = {}

    def fake_render(**kwargs):
        captured["color_intervals"] = kwargs["color_intervals"]
        captured["ribbon_segments"] = kwargs["ribbon_segments"]
        return {"figure_path": str(kwargs["output"]), "view": "overview"}

    monkeypatch.setattr("scripts.render.render_loom", fake_render)

    outdir = tmp_path / "results"
    exit_code = main(
        [
            "--reference",
            str(reference),
            "--reference-segments",
            "4",
            "--comparisons",
            str(comparison),
            "--outdir",
            str(outdir),
            "--views",
            "overview",
            "--min-contig-length",
            "0",
        ]
    )
    assert exit_code == 0

    comparison_intervals = captured["color_intervals"]["comparison"]["comp1"]
    assert [(i.start, i.end, i.origin) for i in comparison_intervals] == [
        (0, 2, "segment 1"),
        (2, 4, "segment 2"),
        (4, 6, "segment 3"),
        (6, 8, "segment 4"),
    ]
    assert [
        (s.subject_start, s.subject_end, s.comparison_start, s.comparison_end, s.origin)
        for s in captured["ribbon_segments"]
    ] == [
        (0, 2, 0, 2, "segment 1"),
        (2, 4, 2, 4, "segment 2"),
        (4, 6, 4, 6, "segment 3"),
        (6, 8, 6, 8, "segment 4"),
    ]


def test_reference_contigs_filter_is_recorded_in_summary(
    monkeypatch, tmp_path: Path
) -> None:
    reference = tmp_path / "reference.fasta"
    comparison = tmp_path / "comparison.fasta"
    write_fasta(reference, [("chr", "ACGTACGT"), ("plasmid", "TGCATGCA")])
    write_fasta(comparison, [("comp1", "ACGTACGT")])

    monkeypatch.setattr(
        "genome_loom.run_minimap2",
        lambda reference_fasta, comparison_fasta, **kwargs: PairwiseAlignment(
            reference=reference_fasta,
            comparison=comparison_fasta,
            blocks=[],
            aligned_bases=8,
            matching_bases=8,
            ani=99.0,
        ),
    )
    monkeypatch.setattr(
        "scripts.render.render_loom",
        lambda **kwargs: {"figure_path": str(kwargs["output"]), "view": "overview"},
    )

    outdir = tmp_path / "results"
    exit_code = main(
        [
            "--reference",
            str(reference),
            "--reference-contigs",
            "plasmid",
            "--comparisons",
            str(comparison),
            "--outdir",
            str(outdir),
            "--views",
            "overview",
            "--min-contig-length",
            "0",
        ]
    )
    assert exit_code == 0

    summary = json.loads((outdir / "genome-loom.summary.json").read_text())
    assert summary["inputs"]["reference_contigs"] == ["plasmid"]
    assert summary["filters"]["reference_contigs"] == ["plasmid"]
    assert summary["genomes"][0]["contig_count"] == 1
    assert summary["genomes"][0]["contigs"][0]["name"] == "plasmid"


def test_reference_contigs_must_exist(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fasta"
    comparison = tmp_path / "comparison.fasta"
    write_fasta(reference, [("chr", "ACGT")])
    write_fasta(comparison, [("comp1", "ACGT")])

    outdir = tmp_path / "results"
    exit_code = main(
        [
            "--reference",
            str(reference),
            "--reference-contigs",
            "missing_contig",
            "--comparisons",
            str(comparison),
            "--outdir",
            str(outdir),
            "--views",
            "overview",
            "--min-contig-length",
            "0",
        ]
    )
    assert exit_code == 1

    summary = json.loads((outdir / "genome-loom.summary.json").read_text())
    assert summary["status"] == "error"
    assert "missing_contig" in summary["error"]["message"]
