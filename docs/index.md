# getBamDepth

getBamDepth calculates the average depth of coverage for regions in a BED file.

You can give it a BAM, SAM, or CRAM file, and it will run samtools to get the depth. You can also give it a depth file that you already made, and it will read that instead. For each region, it also counts how many bases reach the depth thresholds you choose.

## What it does

- Reads a list of regions from a BED file.
- Gets the depth at each position, either from an alignment file or a depth file.
- Reports the average depth for each region.
- Counts and reports the bases that meet your depth thresholds, as both a number and a percentage.

## When to use it

Use getBamDepth when you want a quick coverage summary for a set of target regions. This is common in panel sequencing, exome work, and any case where you need to check how well your regions were covered.

## Requirements

- Perl 5.10 or later.
- samtools, needed only for the `--bam` option.
- mamba or conda, if you want the easy install.

If you only use the `--depth` option with a depth file you made earlier, you do not need samtools.

## Next steps

- [Installation](installation.md) shows how to set up the tool.
- [Usage](usage.md) shows how to run your first command.
- [Inputs](inputs.md) explains each file and option.
- [Output](output.md) explains the result table.
- [Validation and errors](validation.md) explains the checks and error messages.

## Citation

If you use this tool, please cite samtools:

Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. GigaScience (2021) 10(2) giab008. [PubMed 33590861](https://pubmed.ncbi.nlm.nih.gov/33590861)
