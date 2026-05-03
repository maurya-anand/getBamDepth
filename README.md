# getBamDepth

[![publish](https://img.shields.io/github/actions/workflow/status/maurya-anand/getBamDepth/release.yml)](https://github.com/maurya-anand/getBamDepth/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13356789.svg)](https://doi.org/10.5281/zenodo.13356789)

This tool calculates the average depth of coverage for regions specified in a BED file. The depth of coverage can be calculated from a BAM/SAM/CRAM file using samtools or from a pre-calculated depth file. The script also calculates the number of bases in each region that meet certain depth thresholds.

Supports parallel processing via configurable thread count for faster analysis on large datasets.

## Requirements

- Perl 5.10 or later
- mamba or conda (for managing dependencies and samtools)
- Parallel::ForkManager (installed in conda environment)
- samtools (installed in conda environment, required for --bam option)

## Installation

### Quick Start

```bash
make install
```

This installs:
- Perl dependencies (Getopt::Long, File::Basename) to ~/perl5
- Conda environment (getbamdepth) with samtools and perl-parallel-forkmanager

Then activate the environment:
```bash
conda activate getbamdepth
```

### Manual Installation

If you prefer not to use make:

1. Install Perl dependencies:
```bash
cpanm --local-lib=~/perl5 Getopt::Long File::Basename Parallel::ForkManager
```

2. Install samtools (required for --bam option):
```bash
mamba create -y -n getbamdepth -c bioconda samtools perl-parallel-forkmanager
mamba activate getbamdepth
```

> [!NOTE]
> If you only use the `--depth` option with pre-calculated depth files, samtools is not required.

## Testing

Run the test suite to verify the installation:

```bash
make tests
```

This runs three test cases covering depth file input, BAM/CRAM input, and custom thresholds.

## Example Usage

If using the conda environment:
```bash
mamba run -n getbamdepth ./getBamDepth --bed BED_FILE [--bam BAM_FILE | --depth DEPTH_FILE] [--thresholds THRESHOLDS] [--threads INT] [--output FILE]
```

Or activate the environment first:
```bash
conda activate getbamdepth
./getBamDepth --bed BED_FILE [--bam BAM_FILE | --depth DEPTH_FILE] [--thresholds THRESHOLDS] [--threads INT] [--output FILE]
```

### Examples

Depth file input with default thresholds:
```bash
./getBamDepth --bed example/example-targets.bed --depth example/sample.depth
```

CRAM file input:
```bash
./getBamDepth --bed example/example-targets.bed --bam example/sample.cram
```

BAM file with custom thresholds and parallel processing:
```bash
./getBamDepth --bed example/example-targets.bed --bam example/sample.bam --thresholds 5,10 --threads 4
```

Use 2 threads for lower resource usage:
```bash
./getBamDepth --bed example/example-targets.bed --depth example/sample.depth --threads 2
```

Write output to a file instead of stdout:
```bash
./getBamDepth --bed example/example-targets.bed --depth example/sample.depth --output results.txt
```

## Inputs

- `--bed BED_FILE` (required): Path to the BED file (0-based coordinates). Example: `example/example-targets.bed`
  |      |       |       |      |   |   |
  |------|-------|-------|------|---|---|
  | chr1 | 631032| 636027| Gene1| . | + |
  | chrM | 5922  | 6115  | Gene2| . | + |
  
- `--bam BAM_FILE`: Path to BAM, SAM, or CRAM file. Must be indexed (use `samtools index`). Either this or `--depth` is required.
- `--depth DEPTH_FILE`: Path to pre-calculated depth file from `samtools depth`. Either this or `--bam` is required.
- `--thresholds THRESHOLDS`: Comma-separated list of depth thresholds (default: `10,50`). Example: `--thresholds 5,10,20`
- `--threads INT`: Number of parallel threads for region processing (default: 4). Higher values speed up processing but use more CPU. Example: `--threads 8`
- `--output FILE`: Write output to a file instead of stdout. If not provided, output is printed to stdout. Example: `--output results.txt`

## Output

The script outputs a tab-delimited table to standard output. Example: `example/sample.coverage.out.txt`
  | ID     | chrom | start  | end    | total_bases | region | avg_depth | 10x | 50x | 10x(%) | 50x(%) |
  |--------|-------|--------|--------|-------------|--------|-----------|-----|-----|--------|--------|
  | sample | chr1  | 631033 | 636027 | 4995        | Gene1  | 187.22    | 928 | 473 | 18.58  | 9.47   |
  | sample | chrM  | 5923   | 6115   | 193         | Gene2  | 614.41    | 193 | 148 | 100.00 | 76.68  |

The columns of the table are:

- `sample`: The name of the sample, derived from the BAM or depth file name.
- `chrom`: The chromosome of the region.
- `start`: The start position of the region.
- `end`: The end position of the region.
- `total_bases`: The total number of bases in the region. This is calculated from the 0-based bed file (`end - start + 1`).
- `region`: The name of the region.
- `avg_depth`: This is the average depth of coverage for the region. It gives you an idea of how many times each base in the region was sequenced on average. A higher average depth means that the region was sequenced more times, which generally leads to more reliable results.
- `10x`, `50x`, ..., `100x`: These columns represent the count of bases in the region that have been sequenced at least a certain number of times. For instance, the `10x` column indicates the number of bases that have been sequenced at least 10 times.
- `10x(%)`, `50x(%)`, ..., `100x(%)`: These columns represent the percentage of bases in the region that have been sequenced at least a certain number of times. For example, the `10x(%)` column shows the percentage of bases that have been sequenced at least 10 times.

## Reference

- Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H, **Twelve years of SAMtools and BCFtools**, GigaScience (2021) 10(2) giab008 [[33590861]](https://pubmed.ncbi.nlm.nih.gov/33590861)
