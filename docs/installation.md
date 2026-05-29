# Installation

There are two ways to install getBamDepth. The first uses conda or mamba and sets up everything for you. The second is manual, for people who already have Perl and samtools.

## Install with conda or mamba

This is the easy way. It creates a conda environment called `bd-env` with samtools and Perl, then installs the tool into that environment.

```bash
make install
```

Activate the environment before you run the tool:

```bash
conda activate bd-env
```

Check that it works:

```bash
getBamDepth --help
```

## Manual installation

If you do not want to use conda or mamba, install the dependencies yourself.

You need:

- **Perl 5.10 or later.** It already includes the modules the tool uses: strict, warnings, Getopt::Long, File::Basename, and POSIX.
- **samtools.** Needed only for BAM, SAM, and CRAM input. You can skip it if you only use `--depth`.

Then run the script straight from the project directory:

```text
./getBamDepth --bed BED_FILE [--bam BAM_FILE | --depth DEPTH_FILE]
```

## Uninstall

To remove the conda environment:

```bash
make uninstall
```

Or do it by hand:

```bash
conda env remove -n bd-env
```

## Test the install

Run the included tests to confirm everything is set up:

```bash
make test
```

For a fuller run that also checks the input validation:

```bash
make test-full
```
