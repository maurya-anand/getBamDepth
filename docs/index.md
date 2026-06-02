# getBamDepth

[![publish](https://img.shields.io/github/actions/workflow/status/maurya-anand/getBamDepth/release.yml)](https://github.com/maurya-anand/getBamDepth/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20460474.svg)](https://doi.org/10.5281/zenodo.20460474)

getBamDepth calculates depth of coverage for regions in a BED file.

Input requires:

- One BED file with the regions to measure.
- One depth source.

The depth source must be one of these:

- A BAM or CRAM file.
- A depth file from `samtools depth`.

When the depth source is BAM or CRAM, getBamDepth runs `samtools depth`.

For each BED region, getBamDepth reports average depth and depth-threshold counts.

## Citation

If getBamDepth is used in published work, cite the Zenodo record:

Anand Maurya. (2026). maurya-anand/getBamDepth: v2.0.0 (v2.0.0). Zenodo. [https://doi.org/10.5281/zenodo.20460474](https://doi.org/10.5281/zenodo.20460474)

Please also cite samtools:

Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. GigaScience (2021) 10(2) giab008. [PubMed 33590861](https://pubmed.ncbi.nlm.nih.gov/33590861)
