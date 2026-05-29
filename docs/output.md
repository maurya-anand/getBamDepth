# Output

getBamDepth writes a tab-delimited table.

## Example

Example output from `example/sample.coverage.out.txt`:

| ID | chrom | start | end | total_bases | region | avg_depth | 10x | 50x | 10x(%) | 50x(%) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sample | chr1 | 631033 | 636027 | 4995 | Gene1 | 61.40 | 928 | 473 | 18.58 | 9.47 |
| sample | chrM | 5923 | 6115 | 193 | Gene2 | 614.41 | 193 | 148 | 100.00 | 76.68 |

## Columns

| Column | Meaning |
| --- | --- |
| `ID` | Input file name without the final extension. |
| `chrom` | Chromosome from the BED file. |
| `start` | BED start plus 1. |
| `end` | BED end. |
| `total_bases` | BED end minus BED start. |
| `region` | BED column 4, or `unknown`. |
| `avg_depth` | Sum of depth values divided by `total_bases`. |
| `<threshold>x` | Bases with depth greater than or equal to the threshold. |
| `<threshold>x(%)` | Percent of bases with depth greater than or equal to the threshold. |

Threshold columns depend on `--thresholds`.

## Method and assumptions

Coordinate mapping:

- BED input interval is `[start, end)` in 0-based coordinates.
- Output start is `start + 1`.
- Output end is `end`.
- Region length is `total_bases = end - start`.

Per-region calculations:

- `avg_depth = sum(depth at each position in region) / total_bases`
- `<threshold>x = count(positions with depth >= threshold)`
- `<threshold>x(%) = 100 * <threshold>x / total_bases`

Input assumptions:

- One BED file is required.
- Exactly one depth source is required: `--bam` or `--depth`.
- Depth rows are `chrom`, `1-based position`, `depth`.
- BAM/CRAM mode uses `samtools depth -a -b BED_FILE`.
- Depth file input should be position-sorted.

## Worked example

Input BED:

```tsv
chr1    0   3   A
chr1    1   4   B
```

Input depth:

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

## Empty regions

If BED start equals BED end, `total_bases` is `0`.

For empty regions, `avg_depth`, threshold counts, and threshold percentages are `0`.
