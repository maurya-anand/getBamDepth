# getBamDepth

[![publish](https://img.shields.io/github/actions/workflow/status/maurya-anand/getBamDepth/release.yml)](https://github.com/maurya-anand/getBamDepth/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13356789.svg)](https://doi.org/10.5281/zenodo.13356789)

This tool calculates the average depth of coverage for regions specified in a BED file. The depth of coverage can be calculated from a BAM/SAM/CRAM file using samtools or from a pre-calculated depth file. The script also calculates the number of bases in each region that meet certain depth thresholds.

Supports configurable threading for `samtools depth` to improve performance on large datasets.

## Requirements

- Perl 5.10 or later
- mamba or conda (for managing dependencies and samtools)
- samtools (installed in conda environment, required for --bam option)

## Installation

### Quick Start

```bash
make install
```

This installs:
- Conda environment (`bd-env`) with `samtools` and `perl`

Then activate the environment:

```bash
conda activate bd-env
```

Now you can call `getBamDepth` directly from anywhere (when the environment is active):

### Uninstall

To remove the conda environment:

```bash
make uninstall
```

Or manually:
```bash
conda env remove -n bd-env
```

### Manual Installation

If you prefer not to use conda/mamba, install the required dependencies yourself:

**Required:**
- **Perl 5.10+** — includes all necessary modules: strict, warnings, Getopt::Long, File::Basename, POSIX
- **samtools** — for BAM/CRAM/SAM processing (optional if using `--depth` only)

Then run the script directly from the project directory:

```bash
./getBamDepth --bed BED_FILE [--bam BAM_FILE | --depth DEPTH_FILE]
```

> [!NOTE]
> If you only use the `--depth` option with pre-calculated depth files, samtools is not required.

## Testing

Run the test suite to verify the installation:

```bash
make test
```

## Example Usage

**Option 1: Activate the environment first**

```bash
conda activate bd-env
getBamDepth --bed BED_FILE [--bam BAM_FILE | --depth DEPTH_FILE] [--thresholds THRESHOLDS] [--threads INT] [--output FILE]
```

**Option 2: Run directly without activating**

```bash
mamba run -n bd-env getBamDepth --bed BED_FILE [--bam BAM_FILE | --depth DEPTH_FILE] [--thresholds THRESHOLDS] [--threads INT] [--output FILE]
```

### Examples

Depth file input with default thresholds:

```bash
getBamDepth --bed example/example-targets.bed --depth example/sample.depth
```

CRAM file input:

```bash
getBamDepth --bed example/example-targets.bed --bam example/sample.cram
```

BAM file with custom thresholds:

```bash
getBamDepth --bed example/example-targets.bed --bam example/sample.bam --thresholds 5,10 --threads 4
```

Use 2 threads for lower resource usage (with BAM/CRAM input):

```bash
getBamDepth --bed example/example-targets.bed --bam example/sample.bam --threads 2
```

Write output to a file instead of stdout:

```bash
getBamDepth --bed example/example-targets.bed --depth example/sample.depth --output results.txt
```

## Inputs

- `--bed BED_FILE` (required): Path to the BED file (0-based coordinates). Example: `example/example-targets.bed`
  |      |       |       |      |   |   |
  |------|-------|-------|------|---|---|
  | chr1 | 631032| 636027| Gene1| . | + |
  | chrM | 5922  | 6115  | Gene2| . | + |
  Expected format:
  - Tab-delimited BED file
  - 0-based start, end-exclusive end
  - At least 3 columns: `chrom`, `start`, `end`
  - Column 4 is used as the region name when present
  - The file should be sorted by chromosome and start position for best results
  
- `--bam BAM_FILE`: Path to BAM, SAM, or CRAM file. Must be indexed (use `samtools index`). Either this or `--depth` is required.
- `--depth DEPTH_FILE`: Path to a pre-calculated depth file. Either this or `--bam` is required.
  Expected format:
  - Tab-delimited with 3 columns: `chrom`, `position`, `depth`
  - `position` must be 1-based
  - The file should be sorted by chromosome and position
  - For best accuracy, generate it from the same BED file with:
    ```bash
    samtools depth -a -b targets.bed sample.bam > sample.depth
    ```
- `--thresholds THRESHOLDS`: Comma-separated list of depth thresholds (default: `10,50`). Example: `--thresholds 5,10,20`
- `--threads INT`: Number of threads passed to `samtools depth` (only used with `--bam` input). If not provided, defaults to: detected CPU count minus 2 (minimum 1). If CPU detection fails, defaults to 1. Example: `--threads 8`
- `--output FILE`: Write output to a file instead of stdout. If not provided, output is printed to stdout. Example: `--output results.txt`

Runtime behavior:
- Tab-delimited results are written to stdout or to `--output FILE`
- Status, progress, and error messages are written to the terminal (or `/dev/tty` when available)
- All log messages include timestamps and severity levels (INFO, WARN, ERROR)

## Validation and Error Handling

The script validates all input files and parameters:

**BED file validation:**
- At least 3 columns required (chr, start, end)
- Start and end coordinates must be non-negative integers
- End coordinate must be >= start coordinate
- Empty regions (start == end) are handled without errors

**Depth file validation:**
- Exactly 3 columns required (chr, position, depth)
- Position must be a positive integer
- Depth must be a non-negative integer

**BAM/CRAM file validation:**
- BAM/CRAM files must be **indexed** before use
- Supported index formats: `.bai` (BAM), `.csi` (BAM), `.crai` (CRAM)
- Create an index with: `samtools index sample.bam` or `samtools index sample.cram`
- Missing index will cause an error: `BAM index not found (run: samtools index ...)`
- Without an index, samtools would perform a slow linear scan of the entire file instead of seeking to target regions

**Threshold validation:**
- All thresholds must be positive numbers (integers or floats)
- Thresholds are automatically sorted in ascending order
- Invalid thresholds are rejected with clear error messages

## Output

The script outputs a tab-delimited table to standard output. Example: `example/sample.coverage.out.txt`

| ID     | chrom | start  | end    | total_bases | region | avg_depth | 10x | 50x | 10x(%) | 50x(%) |
|--------|-------|--------|--------|-------------|--------|-----------|-----|-----|--------|--------|
| sample | chr1  | 631033 | 636027 | 4995        | Gene1  | 61.40     | 928 | 473 | 18.58  | 9.47   |
| sample | chrM  | 5923   | 6115   | 193         | Gene2  | 614.41    | 193 | 148 | 100.00 | 76.68  |

The columns of the table are:

- `sample`: The name of the sample, derived from the BAM or depth file name.
- `chrom`: The chromosome of the region.
- `start`: The start position of the region.
- `end`: The end position of the region.
- `total_bases`: The total number of bases in the region. This is calculated from the 0-based BED interval as `end - start`.
- `region`: The name of the region.
- `avg_depth`: This is the average depth of coverage for the region. It gives you an idea of how many times each base in the region was sequenced on average. A higher average depth means that the region was sequenced more times, which generally leads to more reliable results.
- `10x`, `50x`, ..., `100x`: These columns represent the count of bases in the region that have been sequenced at least a certain number of times. For instance, the `10x` column indicates the number of bases that have been sequenced at least 10 times.
- `10x(%)`, `50x(%)`, ..., `100x(%)`: These columns represent the percentage of bases in the region that have been sequenced at least a certain number of times. For example, the `10x(%)` column shows the percentage of bases that have been sequenced at least 10 times.

## Reference

- Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H, **Twelve years of SAMtools and BCFtools**, GigaScience (2021) 10(2) giab008 [[33590861]](https://pubmed.ncbi.nlm.nih.gov/33590861)
