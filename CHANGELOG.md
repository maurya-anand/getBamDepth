# Changelog

All notable changes to this project will be documented in this file.

## [2.0.0] - 2026-06-01

## Added

- Official [documentation] (https://getbamdepth.readthedocs.io).

### Changed

- Reworked coverage calculation to stream depth data instead of loading the full depth file into memory.
- Removed Perl-side region parallelization and `Parallel::ForkManager` usage.
- `--threads` now controls only the thread count passed to `samtools depth` for `--bam` input.
- Default thread count now uses detected CPU count minus two, with a minimum of one thread.
- `samtools depth` now runs with `-a` and uses the user-provided BED file directly via `-b`.
- `avg_depth` is now calculated across all bases in each BED interval, including zero-depth positions.
- `--depth` input is now streamed directly and is expected to be a position-sorted tab-delimited file with columns: `chrom`, `position`, `depth`.
- BED input expectations are now documented explicitly: 0-based, end-exclusive, tab-delimited, with column 4 used as the region name when present.
- Runtime logging now reports effective parameters, BED/depth input loading, active BED intervals, and final output location via plain stderr messages.
- Input-source logging now distinguishes precomputed depth input from BAM/CRAM depth generation via `samtools depth`.
- Test targets now write outputs via `getBamDepth --output` instead of shell redirection.
- Installation now uses Conda/Mamba-managed dependencies only.

### Fixed

- Eliminated false test failures caused by `mamba run` stdout capture behavior.
- Updated expected example output to match full-region average depth calculation.

### Removed

- `Parallel::ForkManager` dependency from the runtime and Conda environment setup.
- `Log::Log4perl` dependency from the runtime and Conda environment setup.

### Notes

- This release changes the meaning of `avg_depth` and updates expected output values accordingly.
- Users relying on the previous covered-positions-only average depth behavior should review downstream assumptions before upgrading.

## [1.1.0] - 2026-04-05

### Added

- Optional `--output FILE` flag to write results to a file instead of stdout.
- `make install` target to create the `getbamdepth` Conda environment with `samtools` and `perl-sys-cpu`.
- `make tests` target with coverage for depth file, BAM, and CRAM inputs.

### Changed

- Makefile now detects and prefers `mamba` over `conda` automatically.
- Test targets now write outputs via `getBamDepth --output` instead of shell redirection.
- Installation now uses Conda/Mamba-managed Perl dependencies only.
- Default thread selection now uses detected CPU count minus two, with a minimum of one thread.
- The configured thread count is now passed to `samtools depth`.
- `samtools depth` now uses the user BED file directly to reduce the calculation time.
- Average depth now includes zero-depth positions for full-region accuracy.
- Coverage calculation now streams depth data instead of loading the full depth file into memory.

## [1.0.0] - 2023-07-04

### Added

- Initial release
- Adapted semantic versioning for Zenodo.
