# Algorithm and assumptions

## Coordinate model

- BED input uses 0-based half-open intervals: `[start, end)`.
- Output `start` is `BED start + 1`.
- Output `end` is `BED end`.
- `total_bases = BED end - BED start`.

## Per-region calculations

- `avg_depth = sum(depth at each position in region) / total_bases`
- `<threshold>x = count(positions with depth >= threshold)`
- `<threshold>x(%) = 100 * <threshold>x / total_bases`

## Assumptions

- One BED file is required.
- Exactly one depth source is required: `--bam` or `--depth`.
- BAM/CRAM mode runs `samtools depth -a -b BED_FILE`.
- Depth records must be tab-delimited: `chrom`, `1-based position`, `depth`.
- Depth file input should be sorted by chromosome and position.

## Worked example

BED:

```tsv
chr1    0   3   A
chr1    1   4   B
```

Depth:

```tsv
chr1    1   5
chr1    2   10
chr1    3   20
chr1    4   0
```

With `--thresholds 10,20`, output is:

```tsv
ID  chrom   start   end total_bases region  avg_depth   10x 20x 10x(%)  20x(%)
t   chr1    1   3   3   A   11.67   2   1   66.67   33.33
t   chr1    2   4   3   B   10.00   2   1   66.67   33.33
```
