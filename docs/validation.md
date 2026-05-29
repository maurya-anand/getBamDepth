# Validation and errors

The tool checks your input files and options before it runs. If something is wrong, it stops and prints a clear error message. This page lists the checks.

## BED file checks

- The file must have at least three columns: chrom, start, and end.
- Start and end must be non-negative whole numbers.
- The end must be greater than or equal to the start.
- Empty regions, where the start equals the end, are allowed and cause no error.

## Depth file checks

- Each line must have exactly three columns: chrom, position, and depth.
- The position must be a positive whole number.
- The depth must be a non-negative whole number.

## BAM and CRAM checks

- The file must be indexed before use.
- Supported index types are `.bai` and `.csi` for BAM, and `.crai` for CRAM.
- Make an index with `samtools index sample.bam` or `samtools index sample.cram`.
- A missing index gives this error: `BAM index not found (run: samtools index ...)`.

The index matters for speed. Without it, samtools has to scan the whole file instead of jumping straight to your regions.

## Threshold checks

- Each threshold must be a positive number. Whole numbers and decimals both work.
- Thresholds are sorted in ascending order for you.
- A bad threshold is rejected with a clear message.

## Messages and logging

- The result table goes to standard output, or to the file you set with `--output`.
- Status, progress, and error messages go to the terminal, or to `/dev/tty` when it is available.
- Every message has a timestamp and a level: INFO, WARN, or ERROR.
