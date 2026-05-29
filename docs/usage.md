# Usage

## Command form

```text
getBamDepth --bed BED_FILE [--bam ALIGNMENT_FILE | --depth DEPTH_FILE] [--thresholds THRESHOLDS] [--threads INT] [--output FILE]
```

`--bed` is required.

Exactly one depth source is required:

- `--bam` for BAM, SAM, or CRAM input.
- `--depth` for a precomputed depth file.

`--bam` and `--depth` cannot be used together.

## Help text

Running `getBamDepth` with no options prints the help text.

```bash
getBamDepth
```

Invalid options also print the help text after the unknown option message.

## Examples

Depth file input with default thresholds:

```bash
getBamDepth \
    --bed example/example-targets.bed \
    --depth example/sample.depth
```

CRAM input with default thresholds:

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

BAM input with custom thresholds and threads:

```bash
getBamDepth \
    --bed example/example-targets.bed \
    --bam example/sample.bam \
    --thresholds 5,10 \
    --threads 4
```

Write the result table to a file:

```bash
getBamDepth \
    --bed example/example-targets.bed \
    --depth example/sample.depth \
    --thresholds 5,10 \
    --output results.txt
```

## Messages

The result table goes to standard output unless `--output` is set.

Status, warning, and error messages have this form:

```text
INFO  2026/05/29 13:00:00 > message text
```

When standard output is a terminal, messages are written to `/dev/tty` if it is available. When output is redirected and `--output` is not set, messages may share standard output with the result table.

Use `--output FILE` when a separate result file is needed.
