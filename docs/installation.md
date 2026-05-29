# Installation

getBamDepth can be installed with the project `Makefile`. It can also be run directly as a Perl script.

## Makefile install

The `make install` target uses `mamba` if it is available. If `mamba` is not available, it uses `conda`.

The target creates a conda environment named `bd-env`. It installs `samtools` and `perl` from `conda-forge` and `bioconda`. It also copies `getBamDepth` into the environment `bin` directory.

```bash
make install
```

Activate the environment before running the installed command:

```bash
conda activate bd-env
```

Show the help text by running the command with no options:

```bash
getBamDepth
```

## Manual use

Manual use needs Perl. BAM, SAM, and CRAM input also needs samtools.

Required Perl modules are part of the Perl core used by this script:

- `strict`
- `warnings`
- `Getopt::Long`
- `File::Basename`
- `POSIX`

Run the script from the project directory:

```bash
./getBamDepth --bed BED_FILE --depth DEPTH_FILE
```

or:

```bash
./getBamDepth --bed BED_FILE --bam ALIGNMENT_FILE
```

## Uninstall

The `make uninstall` target removes the `bd-env` conda environment.

```bash
make uninstall
```

The same action can be run directly with conda:

```bash
conda env remove -n bd-env -y
```

## Test targets

The `make test` target runs three smoke tests. It uses the `bd-env` conda environment.

```bash
make test
```

The `make test-full` target runs the smoke tests and the validation test suite.

```bash
make test-full
```
