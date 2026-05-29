# Output

The tool prints a tab-delimited table. By default it prints on the terminal (`stdout`). Use `--output` to write the output to a file.

Here is the example output from `example/sample.coverage.out.txt`:

| ID     | chrom | start  | end    | total_bases | region | avg_depth | 10x | 50x | 10x(%) | 50x(%) |
|--------|-------|--------|--------|-------------|--------|-----------|-----|-----|--------|--------|
| sample | chr1  | 631033 | 636027 | 4995        | Gene1  | 61.40     | 928 | 473 | 18.58  | 9.47   |
| sample | chrM  | 5923   | 6115   | 193         | Gene2  | 614.41    | 193 | 148 | 100.00 | 76.68  |

## What each column means

- **ID**: The sample name. It comes from the name of the BAM or depth file.
- **chrom**: The chromosome of the region.
- **start**: The start position of the region.
- **end**: The end position of the region.
- **total_bases**: The number of bases in the region. It is the 0-based BED interval, so `end - start`.
- **region**: The name of the region. It comes from the fourth column of the BED file.
- **avg_depth**: The average depth across the region. It tells you how many times each base was sequenced on average. A higher number means the region was sequenced more times, which usually means more reliable results. Zero-depth positions are included in this average.
- **10x, 50x, ...**: For each threshold, the count of bases sequenced at least that many times. The `10x` column is the number of bases sequenced at least 10 times. The column names follow the thresholds you set.
- **10x(%), 50x(%), ...**: The same counts as a percentage of the region. The `10x(%)` column is the share of bases sequenced at least 10 times.

## How avg_depth is counted

The average is the total depth across the whole region divided by `total_bases`. Positions with zero depth still count toward the region size, so they pull the average down. This gives you the true average across the full region, not just the covered part.
