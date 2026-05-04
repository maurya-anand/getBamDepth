# Changelog

All notable changes to this project will be documented in this file.

## [1.1.0] - 2026-04-05

### Added

- Parallelization support via `--threads` to process BED regions concurrently using `Parallel::ForkManager`.
- Optional `--output FILE` flag to write results to a file instead of stdout.
- `make install` target to create the `getbamdepth` Conda environment with `samtools`, `perl-parallel-forkmanager`, and `perl-sys-cpu`.
- `make tests` target with coverage for depth file, BAM, and CRAM inputs.

### Changed

- `get_depth()` now accepts an in-memory hash reference instead of reading from file on each call.
- Depth data is loaded once into memory in the parent process before forking.
- Makefile now detects and prefers `mamba` over `conda` automatically.
- Test targets now write outputs via `getBamDepth --output` instead of shell redirection.
- Installation now uses Conda/Mamba-managed Perl dependencies only.
- Default thread selection now uses detected CPU count minus two, with a minimum of one thread.
- The configured thread count is now passed to `samtools depth`.

## [1.0.0] - 2023-07-04

### Added

- Initial release
- Adapted semantic versioning for Zenodo.
