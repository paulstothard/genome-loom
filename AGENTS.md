# AGENTS.md

Guidance for coding agents working in this repository.

## Project Role

`genome-loom` is a non-interactive Python CLI tool that creates comparative genome ribbon plots from one reference genome and multiple comparison genomes. It aligns genomes with `minimap2`, renders overview/pairwise/neighbor figures, preserves input order and orientation, and writes a machine-readable summary for server callers.

Treat the CLI, view names, figure paths, summary JSON, and alignment/filtering semantics as public server contracts.

## Entry Points

- Main command: `python genome_loom.py`
- Alignment logic: `scripts/align.py`
- FASTA handling: `scripts/fasta.py`
- Rendering: `scripts/render.py`
- Summary helpers: `scripts/summary.py`
- Tests: `tests/`

## Coding Style

- Keep the tool non-interactive and deterministic for batch use.
- Preserve reference/top-row semantics and comparison input order unless an option explicitly changes the visible subset.
- Keep view names stable: `overview`, `reference-pairs`, `all-pairs`, `neighbor`.
- Keep rendering logic separated from alignment and FASTA preparation so failures are easier to isolate.
- Prefer explicit summary fields over requiring server code to infer behavior from image filenames.
- Avoid changing visual defaults casually; generated images are the user-facing product.

## Modification Approach

- Read `README.md` before changing CLI behavior, figure layout, output paths, summary shape, or alignment filtering.
- If adding options, update config parsing, README option table, summary recording, and tests. README option tables should keep the normalized columns `Option`, `Required`, `Default`, and `Description`, with defaults copied from the parser or documented as path-derived values.
- If changing minimap2 invocation or block filtering, add tests around block parsing/filtering and run at least one example render.
- If changing rendering, inspect PNG/SVG outputs manually and verify all requested view families still write summary records.
- If this repo is used as a submodule elsewhere, coordinate API-affecting changes with `ragtag-scaffolding`.

## Testing

Fast local checks:

```bash
conda run -n genome-loom pytest -q
conda run -n genome-loom python genome_loom.py --check
conda run -n genome-loom python genome_loom.py --help
```

Example/full checks:

```bash
conda run -n genome-loom bash rebuild_example_outputs.sh
```

Run fast tests for parser, config, summary, and mocked alignment/rendering edits. Run example rebuilds when changing rendering, output layout, minimap2 settings, or README image examples. Do not run any tests after changes that affect only README files, comments, or whitespace.

## Release Process

Follow `deploy/README.md` for the authoritative release procedure. The `deploy/` path is a symlink to the external deployment-notes area, but release notes are still referenced as `deploy/release-notes-vX.Y.Z.md` from this repo.

Release order, every time:

1. Pick `X.Y.Z` using semantic versioning.
2. Update `genome_loom.py` so `VERSION = "X.Y.Z"` matches the intended tag.
3. Add release notes at `deploy/release-notes-vX.Y.Z.md`.
4. Run the release checks from `deploy/README.md`, including `bash rebuild_example_outputs.sh`, `python -m py_compile genome_loom.py scripts/*.py`, `conda run -n genome-loom python genome_loom.py --help`, and `python genome_loom.py --version`.
5. Check `git status` and review the diff. The release commit must include the version bump, release notes, and any regenerated tracked outputs or docs.
6. Commit the release prep changes on `main`, for example `git add ...` then `git commit -m "Prepare release vX.Y.Z"`.
7. Push `main` with `git push origin main`.
8. Create the annotated tag on the committed release-prep commit: `git tag -a vX.Y.Z -m "genome-loom vX.Y.Z"`.
9. Push the tag with `git push origin vX.Y.Z`.
10. Create the GitHub release with `gh release create vX.Y.Z --title "genome-loom vX.Y.Z" --notes-file "deploy/release-notes-vX.Y.Z.md"` or the GitHub Releases UI.

Do not tag until `python genome_loom.py --version` prints `X.Y.Z`. Tags should point at the release-prep commit, not at an earlier commit and not at uncommitted local changes.

If `ragtag-scaffolding` uses this release through its submodule, update and test that submodule pointer separately.

## Failed Server Job Reproduction

When debugging a failed server job locally:

1. Copy the complete job directory to `/tmp/genome-loom-job-123`.
2. Preserve reference/comparison FASTAs, config JSON, stdout/stderr, summary JSON, generated figures, and retained work directory.
3. Re-run with explicit work and result paths:

```bash
conda run -n genome-loom python genome_loom.py --config /tmp/genome-loom-job-123/run.json --outdir /tmp/genome-loom-job-123/repro-results --work-dir /tmp/genome-loom-job-123/repro-work --force
```

If rendering failed after alignments completed, reuse the retained work directory when possible. If alignment failed, simplify to one comparison FASTA and the same minimap2 preset/mapq thresholds before editing code.

## Server Contract

- Exit nonzero on failure and try to write summary JSON when output paths are known.
- Keep generated figure records in the summary for each requested view.
- Preserve `--summary-output`, `--work-dir`, `--tmpdir`, `--keep-temp`, and `--force` behavior for reproducible server runs.
- Keep top-row role labeling controlled by `--reference-role-label`.
- Do not make directory scanning recursive unless explicitly requested; current comparison directory handling is one level deep.
