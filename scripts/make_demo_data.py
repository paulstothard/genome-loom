from __future__ import annotations

import argparse
import random
from pathlib import Path

from fasta import reverse_complement, write_fasta


def mutate(seq: str, rate: float, rng: random.Random) -> str:
    bases = "ACGT"
    chars = list(seq)
    for i, base in enumerate(chars):
        if rng.random() < rate:
            choices = [b for b in bases if b != base]
            chars[i] = rng.choice(choices)
    return "".join(chars)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("examples/data"))
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    lengths = [42000, 36000, 31000, 26000, 22000, 17000]
    ref_records = []
    ref: dict[str, str] = {}
    for i, length in enumerate(lengths, 1):
        name = f"contig_{i}"
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        ref[name] = seq
        ref_records.append((name, seq))

    write_fasta(args.out_dir / "reference.fasta", ref_records)

    recipes = {
        "comparison_alpha.fasta": [
            ("contig_1", "+"),
            ("contig_2", "+"),
            ("contig_3", "-"),
            ("contig_4", "+"),
            ("contig_5", "+"),
            ("contig_6", "-"),
        ],
        "comparison_beta.fasta": [
            ("contig_2", "+"),
            ("contig_1", "+"),
            ("contig_4", "-"),
            ("contig_3", "+"),
            ("contig_6", "+"),
            ("contig_5", "-"),
        ],
        "comparison_gamma.fasta": [
            ("contig_3", "-"),
            ("contig_4", "-"),
            ("contig_1", "+"),
            ("contig_2", "+"),
            ("contig_5", "+"),
        ],
    }

    for filename, recipe in recipes.items():
        records = []
        for idx, (source, strand) in enumerate(recipe, 1):
            seq = ref[source]
            if strand == "-":
                seq = reverse_complement(seq)
            seq = mutate(seq, 0.012 + idx * 0.0015, rng)
            records.append((f"{Path(filename).stem}_{idx}_{source}_{strand}", seq))
        write_fasta(args.out_dir / filename, records)


if __name__ == "__main__":
    main()
