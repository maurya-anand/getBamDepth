# Usage

Activate the conda environment before you run any command:

```bash
conda activate bd-env
```

## Your first command

This reads the example regions and the example depth file, then prints a coverage table:

```bash
getBamDepth --bed example/example-targets.bed --depth example/sample.depth
```

## Command shape

```text
getBamDepth --bed BED_FILE [--bam BAM_FILE | --depth DEPTH_FILE] [--thresholds THRESHOLDS] [--threads INT] [--output FILE]
```

You always need `--bed`. You also need exactly one source of depth, either `--bam` or `--depth`. You cannot pass both.

See [Inputs](inputs.md) for the full list of options.

## Common examples

Depth file input with the default thresholds:

```bash
getBamDepth --bed example/example-targets.bed --depth example/sample.depth
```

CRAM file input:

```bash
getBamDepth --bed example/example-targets.bed --bam example/sample.cram
```

BAM file with custom thresholds and four threads:

```bash
getBamDepth --bed example/example-targets.bed --bam example/sample.bam --thresholds 5,10 --threads 4
```

Use fewer threads when you want lower resource use:

```bash
getBamDepth --bed example/example-targets.bed --bam example/sample.bam --threads 2
```

Write the result to a file instead of the screen:

```bash
getBamDepth --bed example/example-targets.bed --depth example/sample.depth --output results.txt
```

## Where messages go

The result table goes to standard output, or to the file you set with `--output`.

Status, progress, and error messages go to the terminal. Each message has a timestamp and a level: INFO, WARN, or ERROR. These messages stay separate from the result table, so you can redirect the table to a file without mixing the two.
