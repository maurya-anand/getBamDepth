# getBamDepth

This Perl script calculates the average depth of coverage for regions specified in a BED file. The depth of coverage can be calculated from a BAM file using `samtools` or from a pre-calculated depth file. The script also calculates the number of bases in each region that meet certain depth thresholds.

## Usage

```bash
perl getBamDepth.pl --bed BED_FILE [--bam BAM_FILE | --depth DEPTH_FILE] [--thresholds THRESHOLDS]
```

```bash
perl getBamDepth.pl --bed regions.bed --bam sample.bam --thresholds 10,30,50
```

## Options

- `--bed BED_FILE`: Path to the BED file (0 based). This is a required option.
- `--bam BAM_FILE`: Path to the BAM file. Either this option or `--depth` must be provided.
- `--depth DEPTH_FILE`: Path to the pre-calculated depth file. Either this option or `--bam` must be provided.
- `--thresholds THRESHOLDS`: Comma-separated list of depth thresholds. The script will calculate the number of bases in each region that meet each threshold. Default thresholds are 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.

## Output

The script outputs a tab-delimited table to standard output. The columns of the table are:

- `sample`: The name of the sample, derived from the BAM or depth file name.
- `chrom`: The chromosome of the region.
- `start`: The start position of the region.
- `end`: The end position of the region.
- `region`: The name of the region.
- `avg_depth`: The average depth of coverage for the region.
- `10x`, `20x`, ..., `100x`: The number of bases in the region that meet each depth threshold.

## Dependencies

- Perl 5
- samtools (if using the `--bam` option)
