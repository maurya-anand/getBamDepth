# Installation

## Makefile install

`make install` uses `mamba` when available. Otherwise it uses `conda`.

It creates a conda environment named `bd-env` with `samtools` and `perl`. It copies `getBamDepth` into the environment `bin` directory.

```bash
make install
```

Activate the environment:

```bash
conda activate bd-env
```

Print the usage text:

```bash
getBamDepth
```

## Manual use

Manual use needs Perl. BAM and CRAM input also needs samtools.

Run from the project directory:

```bash
./getBamDepth --bed BED_FILE --depth DEPTH_FILE
```

or:

```bash
./getBamDepth --bed BED_FILE --bam ALIGNMENT_FILE
```

## Uninstall

`make uninstall` removes the `bd-env` conda environment.

```bash
make uninstall
```

## Tests

`make test` runs three smoke tests.

```bash
make test
```

`make test-full` runs the smoke tests and the validation test suite.

```bash
make test-full
```
