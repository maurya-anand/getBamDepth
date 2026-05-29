# Usage

## Command

```text
getBamDepth --bed BED_FILE [--bam BAM_FILE | --depth DEPTH_FILE] [--thresholds THRESHOLDS] [--threads INT] [--output FILE]
```

`--bed` is required for coverage calculation.

Use exactly one depth source:

- `--bam` for BAM or CRAM input.
- `--depth` for a depth file.

## Parameters

```bash
getBamDepth --help
```

```text
Usage: getBamDepth --bed BED_FILE [--bam BAM_FILE | --depth DEPTH_FILE] [--thresholds THRESHOLDS] [--threads INT] [--output FILE]
       getBamDepth --help

Parameters:
  --help, -h
      Print this help text and exit.
  --bed BED_FILE
      Required BED file with 0-based intervals.
      Expected columns: chrom, start, end, region_name.
  --bam BAM_FILE
      BAM/CRAM input. The script runs samtools depth -a -b BED_FILE.
  --depth DEPTH_FILE
      Precomputed depth file input. The script streams this file directly.
      Expected format: chrom, 1-based position, depth.
      The file should be position-sorted and should ideally be generated with:
      samtools depth -a -b BED_FILE BAM_FILE
  --thresholds THRESHOLDS
      Comma-separated depth thresholds. Default: 10,50
  --threads INT
      Threads passed to samtools depth for --bam input only.
      Default: detected CPU count minus 2, minimum 1.
  --output FILE
      Write tab-delimited output to a file instead of stdout.
      Status and progress messages are written to the terminal when available.
```

## Examples

Depth file input:

```bash
getBamDepth \
    --bed example/example-targets.bed \
    --depth example/sample.depth
```

CRAM input:

```bash
getBamDepth \
    --bed example/example-targets.bed \
    --bam example/sample.cram
```

BAM input with custom thresholds:

```bash
getBamDepth \
    --bed example/example-targets.bed \
    --bam example/sample.bam \
    --thresholds 5,10
```

BAM input with custom threads:

```bash
getBamDepth \
    --bed example/example-targets.bed \
    --bam example/sample.bam \
    --threads 4
```

Write output to a file:

```bash
getBamDepth \
    --bed example/example-targets.bed \
    --depth example/sample.depth \
    --output results.txt
```

## Messages

The result table goes to standard output unless `--output` is set.

Messages use this format:

```text
INFO  2026/05/29 13:00:00 > message text
```

When standard output is a terminal, messages are written to `/dev/tty` if it is available.
