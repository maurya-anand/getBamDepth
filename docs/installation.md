# Installation

## Installation methods

1. Recommended: automatic installation with the Makefile.
2. Alternative (optional): manual install section at the end of this page.

## Which path should I use?

- Use the Makefile path if you want the easiest setup and built-in test/uninstall commands.
- Use the alternative manual path only if you prefer running from the repository without `make install`.

## Recommended: Automatic install (Makefile)

### Requirements

- `mamba` or `conda` available in `PATH`.

`make install` creates `bd-env` and installs `perl` and `samtools` automatically.

### Install

`make install` uses `mamba` when available and falls back to `conda`.
It creates a conda environment named `bd-env` with `samtools` and `perl`, then installs `getBamDepth` into that environment's `bin` directory.

```bash
make install
```

Activate the environment:

```bash
conda activate bd-env
```

Verify installation:

```bash
getBamDepth --help
```

## Maintenance for Makefile install

Remove the `bd-env` environment:

```bash
make uninstall
```

Run quick smoke tests (3 tests):

```bash
make test
```

Run smoke tests plus validation tests:

```bash
make test-full
```

## Optional: Manual install

Manual install runs the script from this project directory and does not run `make install`.

### Requirements for manual install

- Perl.
- `samtools` only if you use `--bam` with BAM/CRAM input.

Run with a depth file:

```bash
./getBamDepth --bed BED_FILE --depth DEPTH_FILE
```

Run with an alignment file:

```bash
./getBamDepth --bed BED_FILE --bam ALIGNMENT_FILE
```
