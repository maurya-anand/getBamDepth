# Validation and errors

getBamDepth validates input files and options before and during processing. If a check fails, the script writes an `ERROR` message and exits with status `1`.

## Required options

- `--bed` is required.
- One depth source is required.
- The depth source must be `--bam` or `--depth`.
- `--bam` and `--depth` cannot be used together.

## File checks

- The BED file must exist.
- The BAM, SAM, or CRAM file must exist when `--bam` is used.
- The depth file must exist when `--depth` is used.
- An index file must exist when `--bam` is used.

Accepted index suffixes:

- `.bai`
- `.csi`
- `.crai`

Missing alignment indexes produce this error form:

```text
BAM index not found (run: samtools index ALIGNMENT_FILE)
```

The message says `BAM index` for BAM, SAM, and CRAM input.

## BED checks

- Each non-empty BED line must have at least three columns.
- `start` must be a non-negative integer.
- `end` must be a non-negative integer.
- `end` must be greater than or equal to `start`.
- Empty regions are allowed when `start` equals `end`.

## Depth file checks

- Each non-empty depth line must have exactly three columns.
- Position must be a positive integer.
- Depth must be a non-negative integer.

## Threshold checks

- Each threshold must be a positive number.
- Zero is rejected.
- Negative values are rejected.
- Non-numeric values are rejected.
- Decimal values are accepted.
- Thresholds are sorted from low to high.

## Logging

Messages use this format:

```text
LEVEL YYYY/MM/DD HH:MM:SS > message text
```

Message levels are:

- `INFO`
- `WARN`
- `ERROR`

The result table is written to standard output or to the file set by `--output`.

When standard output is a terminal, status messages are written to `/dev/tty` if it is available. When standard output is redirected and `--output` is not set, status messages may be written to standard output with the result table.
