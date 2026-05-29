# getBamDepth

[![publish](https://img.shields.io/github/actions/workflow/status/maurya-anand/getBamDepth/release.yml)](https://github.com/maurya-anand/getBamDepth/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13356789.svg)](https://doi.org/10.5281/zenodo.13356789)

getBamDepth calculates depth of coverage for regions in a BED file.

The depth can come from a BAM, SAM, or CRAM file. In that mode, getBamDepth runs `samtools depth`. The depth can also come from a depth file made earlier.

For each BED region, getBamDepth reports average depth. It also reports how many bases reach each depth threshold.

## Main features

- Reads target regions from a BED file.
- Reads depth from an alignment file or from a depth file.
- Reports one output row for each BED region.
- Reports average depth for each region.
- Reports threshold counts and threshold percentages.
- Supports custom depth thresholds.

## Requirements

- Perl 5.10 or later.
- samtools for BAM, SAM, and CRAM input.
- mamba or conda for the `make install` target.

samtools is not needed when only `--depth` input is used.

## Documentation pages

- [Installation](installation.md)
- [Usage](usage.md)
- [Inputs](inputs.md)
- [Output](output.md)
- [Validation and errors](validation.md)

## Citation

If getBamDepth is used in published work, cite the Zenodo record:

Anand Maurya. (2024). maurya-anand/getBamDepth: v1.0.0 (v1.0.0). Zenodo. [https://doi.org/10.5281/zenodo.13356789](https://doi.org/10.5281/zenodo.13356789)

Also cite samtools when BAM, SAM, or CRAM input is used:

Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. GigaScience (2021) 10(2) giab008. [PubMed 33590861](https://pubmed.ncbi.nlm.nih.gov/33590861)
