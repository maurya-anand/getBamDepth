# getBamDepth

[![publish](https://img.shields.io/github/actions/workflow/status/maurya-anand/getBamDepth/release.yml)](https://github.com/maurya-anand/getBamDepth/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13356789.svg)](https://doi.org/10.5281/zenodo.13356789)

This tool calculates the average depth of coverage for regions specified in a BED file. The depth of coverage can be calculated from a BAM/SAM/CRAM file using `samtools` or from a pre-calculated depth file. The script also calculates the number of bases in each region that meet certain depth thresholds.

> [!IMPORTANT]  
> Please make sure to install the dependencies listed below.

## Requirements

- Perl5
- [samtools](https://www.htslib.org/) (if using the `--bam` option)

## Example Usage

> [!NOTE]  
> Before running the script, make sure it is executable.

```bash
chmod a+x ./getBamDepth
```

```bash
./getBamDepth --bed BED_FILE [--bam BAM_FILE | --depth DEPTH_FILE] [--thresholds THRESHOLDS]
```

```bash
./getBamDepth --bed example/example-targets.bed --depth example/sample.depth
```

```bash
./getBamDepth --bed example/example-targets.bed --bam example/sample.cram
```

```bash
./getBamDepth --bed example/example-targets.bed --bam example/sample.bam --thresholds 5,10
```

## Inputs

- `--bed BED_FILE`: This is a mandatory parameter. It specifies the path to the BED file, which uses a 0-based index. Example: `example/example-targets.bed`
  |      |       |       |      |   |   |
  |------|-------|-------|------|---|---|
  | chr1 | 631032| 636027| Gene1| . | + |
  | chrM | 5922  | 6115  | Gene2| . | + |
  
- `--bam BAM_FILE`: This parameter is used to provide the path to the input file, which can be a BAM, SAM, or CRAM file. Please ensure that the input file is indexed before using it. Use `samtools index -b BAM_FILE` command to create an index. You must provide either this parameter or the `--depth` parameter.
- `--depth DEPTH_FILE`: This parameter is used to specify the path to a depth file that has been pre-calculated using the `samtools depth` command. You must provide either this parameter or the `--bam` parameter.
- `--thresholds THRESHOLDS`: This parameter is used to specify a list of depth thresholds, separated by commas. The script will then calculate the number of bases in each region that meet each of these thresholds. If not provided, the default thresholds used are 10 and 50 (`--thresholds 10,50`).

## Output

The script outputs a tab-delimited table to standard output. Example: `example/sample.coverage.out.txt`
  | ID     | chrom | start  | end    | total_bases | region | avg_depth | 10x | 50x | 10x(%) | 50x(%) |
  |--------|-------|--------|--------|-------------|--------|-----------|-----|-----|--------|--------|
  | sample | chr1  | 631033 | 636027 | 4995        | Gene1  | 187.22    | 928 | 473 | 18.58  | 9.47   |
  | sample | chrM  | 5923   | 6115   | 193         | Gene2  | 614.41    | 193 | 148 | 100.00 | 76.68  |

The columns of the table are:

- `sample`: The name of the sample, derived from the BAM or depth file name.
- `chrom`: The chromosome of the region.
- `start`: The start position of the region.
- `end`: The end position of the region.
- `total_bases`: The total number of bases in the region. This is calculated from the 0-based bed file (`end - start + 1`).
- `region`: The name of the region.
- `avg_depth`: This is the average depth of coverage for the region. It gives you an idea of how many times each base in the region was sequenced on average. A higher average depth means that the region was sequenced more times, which generally leads to more reliable results.
- `10x`, `50x`, ..., `100x`: These columns represent the count of bases in the region that have been sequenced at least a certain number of times. For instance, the `10x` column indicates the number of bases that have been sequenced at least 10 times.
- `10x(%)`, `50x(%)`, ..., `100x(%)`: These columns represent the percentage of bases in the region that have been sequenced at least a certain number of times. For example, the `10x(%)` column shows the percentage of bases that have been sequenced at least 10 times.

## Refernce

- Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H, **Twelve years of SAMtools and BCFtools**, GigaScience (2021) 10(2) giab008 [[33590861]](https://pubmed.ncbi.nlm.nih.gov/33590861)
