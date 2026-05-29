# Inputs

## `--bed BED_FILE`

Required.

BED file rules:

- Tab-delimited.
- Empty lines are ignored.
- At least three columns: `chrom`, `start`, `end`.
- Column 4 is used as the region name when present.
- Missing region names are reported as `unknown`.
- Coordinates are 0-based.
- `start` and `end` must be non-negative integers.
- `end` must be greater than or equal to `start`.

Example:

```tsv
chr1	631032	636027	Gene1	.	+
chrM	5922	6115	Gene2	.	+
```

The output start is `BED start + 1`.

## `--bam BAM_FILE`

Optional.

BAM, SAM, or CRAM input. Use this option instead of `--depth`.

The input file must exist. One of these index files must also exist:

- `BAM_FILE.bai`
- `BAM_FILE.csi`
- `BAM_FILE.crai`

The script runs:

```text
samtools depth -a -@ THREADS -b BED_FILE BAM_FILE
```

## `--depth DEPTH_FILE`

Optional.

Precomputed depth input. Use this option instead of `--bam`.

Depth file rules:

- Tab-delimited.
- Empty lines are ignored.
- Exactly three columns: `chrom`, `position`, `depth`.
- Position is 1-based.
- Position must be a positive integer.
- Depth must be a non-negative integer.

Recommended depth command:

```bash
samtools depth -a -b targets.bed sample.bam > sample.depth
```

## `--thresholds THRESHOLDS`

Optional.

Comma-separated depth thresholds. Default: `10,50`.

Threshold rules:

- Each threshold must be a positive number.
- Decimal values are accepted.
- Values are sorted from low to high.

Example:

```bash
--thresholds 5,10,20
```

## `--threads INT`

Optional.

Threads passed to `samtools depth` with `--bam` input.

Default: detected CPU count minus 2, minimum 1.

`--threads` is not used with `--depth` input.

## `--output FILE`

Optional.

Writes the result table to a file. Without `--output`, the result table is written to standard output.
