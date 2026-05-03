# Changelog

All notable changes to this project will be documented in this file.

## [1.1.0] - 2026-03-05

### Added

- Parallelization support via `--threads` flag (default: 4) to process BED regions concurrently using `Parallel::ForkManager`.
- Optional `--output FILE` flag to write results to a file instead of stdout.
- `make install` target to automate installation of Perl dependencies and conda environment (getbamdepth) with `samtools`, and `perl-parallel-forkmanager`.
- `make tests` target with test coverage for depth file, BAM, and CRAM inputs.

### Changed

- `get_depth()` now accepts in-memory hash reference instead of reading from file on each call.
- Depth file pre-loaded once into memory in parent process before forking (eliminates repeated disk I/O).
- Makefile now detects and prefers mamba over conda automatically.

## [1.0.0] - 2023-07-04

### Added

- Initial release
- Adapted semantic versioning for Zenodo.
