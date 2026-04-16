from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Contig:
    name: str
    sequence: str
    source_count: int = 1
    is_remainder: bool = False

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class Genome:
    name: str
    path: Path
    contigs: list[Contig]

    @property
    def length(self) -> int:
        return sum(c.length for c in self.contigs)


def read_fasta(path: Path, *, min_contig_length: int = 0) -> Genome:
    """Read a FASTA file, filtering short contigs."""
    contigs: list[Contig] = []
    name: str | None = None
    parts: list[str] = []

    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seq = "".join(parts).upper()
                    if len(seq) >= min_contig_length:
                        contigs.append(Contig(name=name, sequence=seq))
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)

    if name is not None:
        seq = "".join(parts).upper()
        if len(seq) >= min_contig_length:
            contigs.append(Contig(name=name, sequence=seq))

    if not contigs:
        raise ValueError(f"No contigs >= {min_contig_length} bp found in {path}")

    return Genome(name=path.stem, path=path, contigs=contigs)


def cap_contig_blocks(genome: Genome, max_contigs: int = 0, gap_size: int = 100) -> Genome:
    """Collapse smaller contigs so a genome never renders more than max_contigs blocks."""
    if max_contigs <= 0 or len(genome.contigs) <= max_contigs:
        return genome
    if max_contigs < 2:
        raise ValueError("--max-contigs must be 0 or at least 2")

    retain_count = max_contigs - 1
    ranked = sorted(
        enumerate(genome.contigs),
        key=lambda item: (-item[1].length, item[0]),
    )
    retained_indexes = {index for index, _contig in ranked[:retain_count]}
    retained = [
        contig
        for index, contig in enumerate(genome.contigs)
        if index in retained_indexes
    ]
    remainder = [
        contig
        for index, contig in enumerate(genome.contigs)
        if index not in retained_indexes
    ]
    remainder_sequence = ("N" * gap_size).join(contig.sequence for contig in remainder)
    remainder_contig = Contig(
        name=f"remaining_contigs_{len(remainder)}",
        sequence=remainder_sequence,
        source_count=len(remainder),
        is_remainder=True,
    )
    return Genome(name=genome.name, path=genome.path, contigs=[*retained, remainder_contig])


def write_fasta(path: Path, records: list[tuple[str, str]], width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), width):
                handle.write(seq[i : i + width] + "\n")


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]
