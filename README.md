# bigbio/quantmsandes

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/bigbio/quantmsandes)
[![GitHub Actions CI Status](https://github.com/bigbio/quantmsandes/actions/workflows/nf-test.yml/badge.svg)](https://github.com/bigbio/quantmsandes/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/bigbio/quantmsandes/actions/workflows/linting.yml/badge.svg)](https://github.com/bigbio/quantmsandes/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/bigbio/quantmsandes)

## Introduction

**bigbio/quantmsandes** is an exploratory pipeline for data-dependent acquisition (DDA)
proteomics, built on the [andes](https://github.com/bigbio/andes) search engine and following
[nf-core](https://nf-co.re/) guidelines. Use it to try new algorithms and frontier approaches in
proteomics, such as glycoproteomics, proteoform identification and quantification, and agentic
analysis, before they are mature enough for production.

> [!IMPORTANT]
> **The production pipeline is [quantms](https://github.com/bigbio/quantms).** quantmsandes is
> for exploring new methods: its defaults, parameters and outputs may change between releases.
> For routine analyses, use quantms.

## The quantms ecosystem

quantmsandes is one of three pipelines in the quantms ecosystem. Each is independent; pick the
one that fits the use case.

| Pipeline                                               | Use it for                                                                               | Engines                                           |
| ------------------------------------------------------ | ---------------------------------------------------------------------------------------- | ------------------------------------------------- |
| [quantms](https://github.com/bigbio/quantms)           | Production analysis of DDA label-free, TMT/iTRAQ and DIA data                            | OpenMS-based workflow with several search engines |
| [quantmsdiann](https://github.com/bigbio/quantmsdiann) | DIA analysis with DIA-NN                                                                 | DIA-NN                                            |
| **quantmsandes**                                       | Exploratory DDA analysis: new algorithms, glycoproteomics, proteoforms, agentic analysis | andes                                             |

All three read [SDRF](https://github.com/bigbio/proteomics-metadata-standard) sample metadata and
report quality control with [pmultiqc](https://github.com/bigbio/pmultiqc).

## What quantmsandes is for

- **New algorithms in identification and quantification.** andes is a data-driven search
  engine with its own scoring models, co-isolated (chimeric) peptide search, PTM refinement,
  and label-free quantification with match-between-runs.
- **Glycoproteomics.** Intact N-glycopeptide search against species-specific glycan databases.
- **Proteoforms.** Searching sequence variants, isoforms, mature chains and annotated
  modification sites from UniProt, PEFF or annotation files, with knockoff-based error
  estimates for variant identifications, and a per-protein-family proteoform table with
  quantities.
- **Agentic analysis.** andes exposes many advanced parameters. quantmsandes passes them
  through, so an AI agent such as Claude can choose and tune them for a given dataset and
  question, and compare the outcomes.

## Status

> [!WARNING]
> quantmsandes is in early development. The repository currently holds the pipeline template
> and the [design](docs/design.md); the workflow is not runnable yet.

Planned workflow for DDA label-free and TMT/iTRAQ data:

1. **Input validation:** SDRF parsing and validation with
   [sdrf-pipelines](https://github.com/bigbio/sdrf-pipelines).
2. **Search:** andes builds its fragment index once, then searches each run (Thermo `.raw`,
   Bruker `.d` and mzML are read directly).
3. **PSM rescoring:** one [Percolator](https://github.com/percolator/percolator) run over all
   runs.
4. **Quantification:** andes label-free quantification with match-between-runs, or reporter
   ion quantification for TMT/iTRAQ.
5. **Protein inference:** [OpenMS](https://www.openms.de/) protein inference.
6. **Protein quantification:** [mokume](https://github.com/bigbio/mokume).
7. **Quality control:** [pmultiqc](https://github.com/bigbio/pmultiqc).

The containers come from [quantms-containers](https://github.com/bigbio/quantms-containers).

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set up Nextflow.

Once the workflow is available, a run will look like:

```bash
nextflow run bigbio/quantmsandes \
   -profile <docker/singularity/.../institute> \
   --input experiment.sdrf.tsv \
   --database proteins.fasta \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

bigbio/quantmsandes was originally written by Yasset Perez-Riverol.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use bigbio/quantmsandes for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
