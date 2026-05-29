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

## Empty regions

If BED start equals BED end, `total_bases` is `0`.

For empty regions, `avg_depth`, threshold counts, and threshold percentages are `0`.
