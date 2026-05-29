# getBamDepth

[![publish](https://img.shields.io/github/actions/workflow/status/maurya-anand/getBamDepth/release.yml)](https://github.com/maurya-anand/getBamDepth/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13356789.svg)](https://doi.org/10.5281/zenodo.13356789)

getBamDepth calculates depth of coverage for regions in a BED file.

Input requires:

- One BED file with the regions to measure.
- One depth source.

The depth source must be one of these:

- A BAM, SAM, or CRAM file.
- A depth file from `samtools depth`.

When the depth source is BAM, SAM, or CRAM, getBamDepth runs `samtools depth`.

For each BED region, getBamDepth reports average depth and depth-threshold counts.

## Requirements

- Perl.
- samtools for BAM, SAM, and CRAM input.
- mamba or conda for `make install`, `make uninstall`, and test targets.

samtools is not needed when `--depth` input is used.

## Pages

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
