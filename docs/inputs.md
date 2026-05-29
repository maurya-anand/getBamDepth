# Inputs

## `--bed BED_FILE`

Required.

The BED file defines the regions to measure.

Rules:

- The file must be tab-delimited.
- Empty lines are ignored.
- At least three columns are required: `chrom`, `start`, and `end`.
- The fourth column is used as the region name when present.
- If the fourth column is missing, the region name is `unknown`.
- BED coordinates are 0-based.
- The BED start is included.
- The BED end is not included.
- `start` must be a non-negative integer.
- `end` must be a non-negative integer.
- `end` must be greater than or equal to `start`.

Example:

```tsv
chr1	631032	636027	Gene1	.	+
chrM	5922	6115	Gene2	.	+
```

The script converts BED start positions to 1-based positions in the output. For example, BED start `631032` is reported as output start `631033`.

## `--bam ALIGNMENT_FILE`

Optional.

Path to a BAM, SAM, or CRAM file. Use `--bam` or `--depth`, not both.

The file must exist. An index file must also exist. The script accepts these index names:

- `ALIGNMENT_FILE.bai`
- `ALIGNMENT_FILE.csi`
- `ALIGNMENT_FILE.crai`

For BAM input, an index can be made with:

```bash
samtools index sample.bam
```

For CRAM input, an index can be made with:

```bash
samtools index sample.cram
```

With `--bam`, getBamDepth runs this command internally:

```text
samtools depth -a -@ THREADS -b BED_FILE ALIGNMENT_FILE
```

`-a` keeps positions with zero depth in the depth stream.

## `--depth DEPTH_FILE`

Optional.

Path to a precomputed depth file. Use `--depth` or `--bam`, not both.

Rules:

- The file must be tab-delimited.
- Empty lines are ignored.
- Each non-empty line must have exactly three columns.
- Column 1 is chromosome.
- Column 2 is 1-based position.
- Column 3 is depth.
- Position must be a positive integer.
- Depth must be a non-negative integer.
- The file should be sorted by chromosome and position.

Recommended command:

```bash
samtools depth -a -b targets.bed sample.bam > sample.depth
```

Use the same BED file for depth creation and getBamDepth. This keeps zero-depth positions in the target regions.

## `--thresholds THRESHOLDS`

Optional.

Comma-separated depth thresholds. The default is `10,50`.

Rules:

- Each threshold must be a positive number.
- Whole numbers are accepted.
- Decimal values are accepted.
- Thresholds are sorted from low to high before output columns are written.

Example:

```bash
--thresholds 5,10,20
```

For each threshold, the output includes one count column and one percentage column.

## `--threads INT`

Optional.

Number of threads passed to `samtools depth` when `--bam` is used.

Default value:

- Detected CPU count minus 2.
- Minimum value is 1.
- If CPU count cannot be detected, the value is 1.

`--threads` has no effect with `--depth` input.

Example:

```bash
--threads 8
```

## `--output FILE`

Optional.

Writes the result table to a file.

If `--output` is not set, the result table is written to standard output.

Example:

```bash
--output results.txt
```
