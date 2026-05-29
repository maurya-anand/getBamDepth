# Inputs

This page explains every input file and option.

## --bed BED_FILE (required)

The BED file lists the regions you want to measure.

- It is tab-delimited.
- Coordinates are 0-based. The start is included and the end is not.
- It needs at least three columns: `chrom`, `start`, and `end`.
- If a fourth column is present, it is used as the region name.

Example, from `example/example-targets.bed`:

```tsv
chr1  631032  636027    Gene1   .   +
chrM    5922    6115    Gene2   .   +
```

## --bam BAM_FILE

Path to a BAM, SAM, or CRAM file. Use this or `--depth`, not both.

The file must be indexed. Make an index with:

```bash
samtools index sample.bam
```

The tool automatically runs `samtools depth -a -b BED_FILE` on this file to get the depth at each position.

## --depth DEPTH_FILE

Path to a depth file you made earlier. Use this or `--bam`, not both.

The format is:

- Tab-delimited with three columns: `chrom`, `position`, and `depth`.
- The `position` is 1-based.
- The file should be sorted by chromosome and then by position.

For the best match with your regions, make the depth file from the same BED file:

```bash
samtools depth -a -b targets.bed sample.bam > sample.depth
```

## --thresholds THRESHOLDS

A comma-separated list of depth thresholds. The default is `10,50`.

For each threshold, the output gets two extra columns: the number of bases that reach that depth, and the percentage. The thresholds are sorted in ascending order for you.

Example:

```bash
--thresholds 5,10,20
```

## --threads INT

The number of threads passed to `samtools depth`. This only matters for `--bam` input. It has no effect when you use `--depth`.

If you do not set it, the tool uses the detected CPU count minus two, with a minimum of one. If it cannot detect the CPU count, it uses one.

Example:

```bash
--threads 8
```

## --output FILE

Write the result table to a file instead of the screen. If you do not set it, the table goes to standard output.

Example:

```bash
--output results.txt
```
