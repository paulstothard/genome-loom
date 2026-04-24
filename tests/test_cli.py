from __future__ import annotations

import json
from pathlib import Path

from genome_loom import _parse_args, build_parser, main
from scripts.align import PairwiseAlignment


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
