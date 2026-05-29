# Validation and errors

getBamDepth exits with status `1` when validation fails after option parsing.

## Required input

- `--bed` is required.
- Either `--bam` or `--depth` is required.
- `--bam` and `--depth` cannot be used together.

## File checks

- The BED file must exist.
- The `--bam` file must exist when `--bam` is used.
- The `--depth` file must exist when `--depth` is used.
- The `--bam` file must have an index file.

Accepted index suffixes:

- `.bai`
- `.csi`
- `.crai`

## BED checks

- Each non-empty line must have at least three columns.
- `start` must be a non-negative integer.
- `end` must be a non-negative integer.
- `end` must be greater than or equal to `start`.

## Depth checks

- Each non-empty line must have exactly three columns.
- Position must be a positive integer.
- Depth must be a non-negative integer.

## Threshold checks

- Thresholds must be positive numbers.
- Decimal values are accepted.
- Thresholds are sorted from low to high.

## Messages

Messages use this format:

```text
LEVEL YYYY/MM/DD HH:MM:SS > message text
```

Message levels are `INFO`, `WARN`, and `ERROR`.
