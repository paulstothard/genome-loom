from __future__ import annotations

from pathlib import Path

from scripts.fasta import Contig, Genome
from scripts.render import ColorInterval, _display_name, _legend_entries_from_intervals


def test_display_name_preserves_accession_punctuation() -> None:
    genome = Genome(
        name="GCF_001900355_1",
        path=Path("GCF_001900355_1.fasta"),
        contigs=[],
        display_name="GCF_001900355.1",
    )

    assert _display_name(genome) == "GCF_001900355.1"


def test_display_name_preserves_underscores_without_override() -> None:
    genome = Genome(
        name="GCF_001900355_1",
        path=Path("GCF_001900355_1.fasta"),
        contigs=[],
    )

    assert _display_name(genome) == "GCF_001900355_1"


def test_interval_legend_preserves_accession_punctuation() -> None:
    genome = Genome(
        name="reference",
        path=Path("reference.fasta"),
        contigs=[Contig(name="NZ_CP010133.1", sequence="A" * 10)],
    )
    intervals = {
        "reference": {
            "NZ_CP010133.1": [
                ColorInterval(
                    genome="reference",
                    contig="NZ_CP010133.1",
                    start=0,
                    end=10,
                    color="#0072B2",
                    origin="NZ_CP010133.1",
                )
            ]
        }
    }

    assert _legend_entries_from_intervals(genome, intervals) == [
        ("NZ_CP010133.1", "NZ_CP010133.1", "#0072B2")
    ]
