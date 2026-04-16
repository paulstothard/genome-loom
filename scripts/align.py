from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class AlignmentBlock:
    reference_contig: str
    reference_start: int
    reference_end: int
    comparison_contig: str
    comparison_start: int
    comparison_end: int
    strand: str
    matches: int
    block_length: int
    mapq: int

    @property
    def identity(self) -> float:
        return self.matches / self.block_length if self.block_length else 0.0


@dataclass(frozen=True)
class PairwiseAlignment:
    reference: Path
    comparison: Path
    blocks: list[AlignmentBlock]
    aligned_bases: int
    matching_bases: int
    ani: float


def run_minimap2(
    reference_fasta: Path,
    comparison_fasta: Path,
    *,
    threads: int = 1,
    preset: str = "asm5",
    min_block_length: int = 500,
    min_mapq: int = 0,
) -> PairwiseAlignment:
    """Align comparison to reference and parse minimap2 PAF from stdout."""
    cmd = [
        "minimap2",
        "-x",
        preset,
        "-t",
        str(max(1, threads)),
        str(reference_fasta),
        str(comparison_fasta),
    ]
    proc = subprocess.run(cmd, check=False, text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"minimap2 failed for {comparison_fasta.name} vs {reference_fasta.name}: "
            f"{proc.stderr.strip()}"
        )

    blocks: list[AlignmentBlock] = []
    aligned = 0
    matches = 0
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        q_name = fields[0]
        q_start = int(fields[2])
        q_end = int(fields[3])
        strand = fields[4]
        t_name = fields[5]
        t_start = int(fields[7])
        t_end = int(fields[8])
        n_match = int(fields[9])
        block_len = int(fields[10])
        mapq = int(fields[11])
        if block_len < min_block_length or mapq < min_mapq:
            continue
        blocks.append(
            AlignmentBlock(
                reference_contig=t_name,
                reference_start=t_start,
                reference_end=t_end,
                comparison_contig=q_name,
                comparison_start=q_start,
                comparison_end=q_end,
                strand=strand,
                matches=n_match,
                block_length=block_len,
                mapq=mapq,
            )
        )
        aligned += block_len
        matches += n_match

    ani = matches / aligned if aligned else 0.0
    return PairwiseAlignment(
        reference=reference_fasta,
        comparison=comparison_fasta,
        blocks=blocks,
        aligned_bases=aligned,
        matching_bases=matches,
        ani=ani,
    )
