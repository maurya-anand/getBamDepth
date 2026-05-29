# Output

getBamDepth writes a tab-delimited table.

By default, the table is written to standard output. With `--output FILE`, the table is written to that file.

## Example

Example output from `example/sample.coverage.out.txt`:

| ID | chrom | start | end | total_bases | region | avg_depth | 10x | 50x | 10x(%) | 50x(%) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sample | chr1 | 631033 | 636027 | 4995 | Gene1 | 61.40 | 928 | 473 | 18.58 | 9.47 |
| sample | chrM | 5923 | 6115 | 193 | Gene2 | 614.41 | 193 | 148 | 100.00 | 76.68 |

## Columns

| Column | Meaning |
| --- | --- |
| `ID` | Sample name from the input file name. The final file extension is removed. |
| `chrom` | Chromosome from the BED file. |
| `start` | 1-based start position. The script adds 1 to the BED start. |
| `end` | End position from the BED file. |
| `total_bases` | Region size from the BED interval. This is `end - start` from the original BED coordinates. |
| `region` | Region name from BED column 4. If column 4 is missing, this is `unknown`. |
| `avg_depth` | Mean depth across the full region. |
| `10x`, `50x`, ... | Number of bases with depth greater than or equal to the threshold. |
| `10x(%)`, `50x(%)`, ... | Percentage of bases with depth greater than or equal to the threshold. |

Threshold column names depend on `--thresholds`. The default thresholds are `10` and `50`.

## Average depth

Average depth is calculated as:

```text
sum of depth values in the region / total_bases
```

Positions with zero depth are included when they are present in the depth stream. For BAM, SAM, and CRAM input, the script runs `samtools depth -a`, so zero-depth positions in the requested BED regions are included.

For `--depth` input, the depth file should include zero-depth positions. Use `samtools depth -a -b BED_FILE` when making the file.

## Empty regions

A BED row with `start` equal to `end` has `total_bases` equal to `0`.

For an empty region:

- `avg_depth` is `0`.
- Threshold counts are `0`.
- Threshold percentages are `0`.
