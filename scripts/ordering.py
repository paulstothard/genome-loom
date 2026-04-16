from __future__ import annotations

from pathlib import Path

from .align import run_minimap2


def order_by_input(comparisons: list[Path]) -> list[Path]:
    return comparisons


def order_by_greedy_ani(
    reference: Path,
    comparisons: list[Path],
    *,
    threads: int = 1,
    preset: str = "asm5",
    min_block_length: int = 500,
) -> tuple[list[Path], dict[str, dict[str, float]]]:
    """Greedily order genomes by ANI, starting with the best match to reference."""
    if len(comparisons) <= 1:
        return comparisons, {}

    scores: dict[str, dict[str, float]] = {}

    def score(a: Path, b: Path) -> float:
        scores.setdefault(a.name, {})
        if b.name not in scores[a.name]:
            aln = run_minimap2(
                a,
                b,
                threads=threads,
                preset=preset,
                min_block_length=min_block_length,
            )
            scores[a.name][b.name] = aln.ani
            scores.setdefault(b.name, {})[a.name] = aln.ani
        return scores[a.name][b.name]

    remaining = comparisons[:]
    first = max(remaining, key=lambda p: score(reference, p))
    ordered = [first]
    remaining.remove(first)

    while remaining:
        prev = ordered[-1]
        nxt = max(remaining, key=lambda p: score(prev, p))
        ordered.append(nxt)
        remaining.remove(nxt)

    return ordered, scores
