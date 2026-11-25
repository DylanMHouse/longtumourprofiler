# nf-core/longtumourprofiler: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of nf-core/longtumourprofiler, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Clonal deconvolution with PyClone-VI (`--pyclonevi`) from the consensus somatic SNVs and the ASCAT copy number and purity. Replaces the unused nf-core `pyclonevi` module.
- Clone trees with PhyClone (`--phyclone`) over the PyClone-VI clusters: MAP tree and tables, topology report and a MAP tree plot. The plot script is adapted from PhyClone-Workflow and stays under GPL-3.0 (`assets/licenses/GPL-3.0.txt`).

### `Fixed`

- Spectre, CNVpytor, SpecImmune and modkit pileup now read the LongPhase haplotagged BAM. Sniffles moves out of the Spectre subworkflow (now `BAM_CALL_CNV_SPECTRE`) so its SVs can still be co-phased. With `--phasing false` every downstream step falls back to the alignment BAM, and ClairS phases the tumour itself ([@DylanMHouse](https://github.com/DylanMHouse))

### `Dependencies`

### `Deprecated`
