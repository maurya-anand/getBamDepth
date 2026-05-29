# getBamDepth

[![publish](https://img.shields.io/github/actions/workflow/status/maurya-anand/getBamDepth/release.yml)](https://github.com/maurya-anand/getBamDepth/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13356789.svg)](https://doi.org/10.5281/zenodo.13356789)
[![docs](https://img.shields.io/badge/docs-readthedocs-blue)](https://getbamdepth.readthedocs.io)

This tool calculates the average depth of coverage for regions in a BED file. The depth can come from a BAM, SAM, or CRAM file using samtools, or from a pre-calculated depth file. For each region, it also counts the bases that meet the depth thresholds you choose.

Full documentation is available at [getbamdepth.readthedocs.io](https://getbamdepth.readthedocs.io).

## Requirements

- Perl 5.10 or later
- samtools (needed only for the `--bam` option)
- mamba or conda (for the easy install)

If you only use the `--depth` option with a pre-calculated depth file, samtools is not required.

## Quick Start

Install the conda environment (`bd-env`) with samtools and Perl:

```bash
make install
```

Activate the environment:

```bash
conda activate bd-env
```

Run a first command:

```bash
getBamDepth --bed example/example-targets.bed --depth example/sample.depth
```

For installation options, inputs, output format, and validation rules, see the [documentation](https://getbamdepth.readthedocs.io).

## Testing

```bash
make test
```

## Citation

If you use this tool, please cite it through its Zenodo record:

Anand Maurya. (2024). maurya-anand/getBamDepth: v1.0.0 (v1.0.0). Zenodo. [https://doi.org/10.5281/zenodo.13356789](https://doi.org/10.5281/zenodo.13356789)

Please also cite samtools:

Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H, **Twelve years of SAMtools and BCFtools**, GigaScience (2021) 10(2) giab008 [[33590861]](https://pubmed.ncbi.nlm.nih.gov/33590861)
