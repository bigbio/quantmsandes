# quantmsandes design

Status: draft, 2026-10-19.

## 1. Summary

quantmsandes is an independent pipeline of the quantms ecosystem for exploratory DDA
proteomics with the andes search engine. andes does identification (with predicted
retention-time and fragment features computed in the same search), quantification and
protein inference; Percolator rescores PSMs once over all runs; mokume quantifies proteins.
pmultiqc reports quality control. The pipeline uses no OpenMS tools.

quantms stays the production pipeline and quantmsdiann the DIA-NN pipeline. quantmsandes
shares their input conventions (SDRF) and reporting (pmultiqc), not their code paths.

Where a parameter means the same as in quantms, quantmsandes uses the quantms name, so users
moving between the pipelines do not learn a second name for the same setting. Parameters
that only andes has are prefixed `andes_`; those of the mokume step are prefixed `protein_`.

## 2. Scope

First version:

- DDA label-free (LFQ), with match-between-runs (MBR).
- DDA TMT and iTRAQ.
- PSM rescoring with Percolator, run once over all runs, on andes's features plus predicted
  retention-time (and later fragment-intensity) features that andes computes during the
  search.
- Protein inference with andes.
- Protein quantification with mokume.
- pmultiqc report.

Exploratory modes that pass andes's own options through (section 5.2):

- Glycoproteomics (`--glyco`).
- Proteoforms: UniProt variants, isoforms and chains, annotation files, knockoff FDR for
  variants, and the proteoform table.
- Chimeric search and PTM refinement.

Not in the first version: DIA (use quantmsdiann), MSstats export, site localisation.

## 3. Process graph

### 3.1 LFQ

```
INPUT_CHECK            SDRF parsing and validation                            [1 task]
ANDES_INDEX            build the fragment index once                          [1 task]
ANDES_SEARCH           per run: search, write PIN + quant tables              [1 task per run]
                       andes --spectrum RUN --database DB --mods MODS
                             --output-pin RUN.pin --quant-dir QDIR
                             --predicted-features rt ...
PERCOLATOR             concatenate all PINs, run Percolator once,             [1 task]
                       write target and decoy PSM results
ANDES_QUANT_LIBRARY    andes quant library --quant-dir QDIR                   [1 task]
                             --psms targets.tsv --decoy-psms decoys.tsv
                             --mbr-library-fdr F
ANDES_QUANT_EXTRACT    per run: andes quant extract --quant-dir QDIR --run S  [1 task per run]
ANDES_QUANT_FINALIZE   andes quant finalize --quant-dir QDIR -o BASE          [1 task]
                             --output-parquet OUT
ANDES_PROTEINS         andes protein inference on the PSMs that pass the      [1 task]
                       Percolator q-value cut-off; writes the protein groups
                       and their protein q-values, and the groups into the
                       feature parquet (pg_accessions, anchor_protein, unique,
                       pg_global_qvalue)
MOKUME_PROTEINS        mokume quantify features2proteins                      [1 task]
                             --parquet OUT/quantms.feature.parquet --sdrf SDRF ...
PMULTIQC
```

Notes:

- `--quant-dir` alone captures what the quant subcommands need; the per-run search does not
  need `--lfq`.
- `andes search --mbr` runs library, extract and finalize in one process. The pipeline
  uses the three subcommands instead, so the extract step runs in parallel across runs.
- One Percolator run over all runs is intended: the MBR library is built from one set of
  q-values, so false library entries do not grow with the number of runs. This needs every
  SpecId to be unique across runs and every per-run PIN to have the same header. andes
  does both from the release after 0.3.0 (the SpecId carries the input's file stem; the PIN
  charge columns follow `--charge` only).
- `quant library` joins the merged Percolator result on (run, SpecId). Two inputs with the
  same file stem collide, so the pipeline must reject duplicate stems at input check.
- `ANDES_INDEX` is its own task so every search task reuses one index.
- Predicted features are computed inside each `ANDES_SEARCH` task, not in a task of their
  own: the spectra are already in memory, so nothing is read twice and no intermediate file
  is staged. The feature columns depend on the flag only, never on the data, so every
  run's PIN has the same header (a row the model cannot encode gets a neutral value and a
  flag). Retention time is calibrated per run before the PINs are merged.

### 3.2 TMT / iTRAQ

```
INPUT_CHECK            [1 task]
ANDES_INDEX            [1 task]
ANDES_SEARCH           per run: andes ... --output-pin RUN.pin --tmt PLEX     [1 task per run]
                             --tmt-level 2|3 --output-parquet RUN/
PERCOLATOR             concatenate all PINs, run once                         [1 task]
FILTER_FEATURES        keep feature rows whose PSM passes the Percolator      [1 task]
                       q-value cut-off; concatenate runs
ANDES_PROTEINS         as for LFQ                                             [1 task]
MOKUME_PROTEINS        --quant-method abd                                     [1 task]
PMULTIQC
```

`andes search --tmt` writes the reporter intensities into `quantms.feature.parquet`:
corrected values in `intensities` (labels `TMT126` … or `ITRAQ114` …) and raw values in
`additional_intensities`.

## 4. Inputs

- An SDRF (`.sdrf.tsv`) with the data files, labels, enzyme and modifications.
- A protein database: FASTA, UniProt XML (for proteoform modes) or PEFF.
- Spectra as Thermo `.raw`, Bruker `.d` or mzML; andes reads all three directly, so there
  is no conversion step.

## 5. Parameters

### 5.1 Parameters with quantms names

Used with the same meaning as in quantms: the input/output, SDRF validation, institutional
and generic options, `database`, `enable_pmultiqc` and `contaminant_string`. The protein
inference parameters keep their quantms names where andes's inference has the same setting
(for example `protein_level_fdr_cutoff`, `use_shared_peptides`,
`min_peptides_per_protein`, `picked_fdr`); the set is fixed with the andes inference
(8.3).

Translated into andes flags:

| Parameter                                                                              | andes flag                                                       |
| -------------------------------------------------------------------------------------- | ---------------------------------------------------------------- |
| `precursor_mass_tolerance`, `precursor_mass_tolerance_unit`                            | `--precursor-tol 20ppm`                                          |
| fixed and variable mods (SDRF), `variable_mods`                                        | mods file passed with `--mods`                                   |
| enzyme (SDRF)                                                                          | `--enzyme`                                                       |
| `num_enzyme_termini`                                                                   | `--enzyme-specificity`                                           |
| `allowed_missed_cleavages`                                                             | `--max-missed-cleavages`                                         |
| `min_precursor_charge`, `max_precursor_charge`                                         | `--charge MIN..MAX`                                              |
| `min_peptide_length`, `max_peptide_length`                                             | `--min-length`, `--max-length`                                   |
| `isotope_error_range`                                                                  | `--isotope-error`                                                |
| `ms2_fragment_method`                                                                  | `--fragmentation`                                                |
| `protocol`                                                                             | `--protocol`                                                     |
| `num_hits`                                                                             | `--top-n`                                                        |
| `max_mods`                                                                             | `--max-mods`                                                     |
| `min_peaks`                                                                            | `--min-peaks`                                                    |
| `add_decoys` = true                                                                    | `--decoy-strategy reverse` (andes generates decoys)              |
| `add_decoys` = false, `decoy_string`, `decoy_string_position`                          | `--decoy-strategy none` and `--decoy-prefix` or `--decoy-suffix` |
| `train_FDR`, `test_FDR`, `subset_max_train`, `klammer`, `description_correct_features` | Percolator arguments                                             |
| `run_fdr_cutoff`                                                                       | PSM q-value cut-off applied after Percolator                     |
| SDRF label                                                                             | `--quant-dir` (LFQ), or `--tmt PLEX`                             |
| `reporter_mass_shift`                                                                  | `--tmt-tol`                                                      |
| `min_precursor_purity`                                                                 | `--tmt-min-purity`                                               |
| `isotope_correction`, `plex_corr_matrix_file`                                          | `--tmt-correction FILE`                                          |

There is no fragment tolerance parameter: andes chooses the fragment tolerance per scoring
model, and forcing a value breaks the models. There is no instrument parameter either:
andes detects the instrument from the spectra.

andes hides some of the flags above (`--charge`, `--top-n`, `--max-missed-cleavages` and
others) from its `--help`, but accepts them.

### 5.2 andes parameters

Defaults follow andes's own defaults unless stated.

| Parameter                  | andes flag                                   | Default                                 |
| -------------------------- | -------------------------------------------- | --------------------------------------- |
| `andes_score`              | `--score`                                    | `auto`                                  |
| `andes_chimeric`           | `--chimeric`                                 | `false`                                 |
| `andes_refine`             | `--refine`                                   | `false`                                 |
| `andes_precursor_cal`      | `--precursor-cal`                            | `off`                                   |
| `andes_precursor_mono`     | `--precursor-mono`                           | `auto`                                  |
| `andes_uniprot_variants`   | `--uniprot-variants`                         | `false`                                 |
| `andes_uniprot_isoforms`   | `--uniprot-isoforms`                         | `false`                                 |
| `andes_uniprot_chains`     | `--uniprot-chains`                           | `false`                                 |
| `andes_uniprot_evidence`   | `--uniprot-evidence`                         | `experimental`                          |
| `andes_annotations`        | `--annotations`                              | not set                                 |
| `andes_knockoffs`          | `--knockoffs`                                | `4`                                     |
| `andes_proteoform_fdr`     | `--proteoform-fdr`                           | not set (andes: max(`--fdr`, 0.05))     |
| `andes_glyco`              | `--glyco`                                    | `false`                                 |
| `andes_glyco_species`      | `--glyco-species`                            | not set                                 |
| `andes_glyco_glycan_db`    | `--glyco-glycan-gdb`                         | not set                                 |
| `andes_glyco_min_core_y`   | `--glyco-min-core-y`                         | `0`                                     |
| `andes_glyco_hcd_pair`     | `--glyco-hcd-pair`                           | `true`                                  |
| `andes_tmt_level`          | `--tmt-level`                                | `2`                                     |
| `andes_lfq_tol`            | `--lfq-tol`                                  | not set (andes default)                 |
| `andes_lfq_rt_window`      | `--lfq-rt-window`                            | `60`                                    |
| `andes_lfq_gaussian_width` | `--lfq-gaussian-width`                       | not set                                 |
| `andes_lfq_min_cosine`     | `--lfq-min-cosine`                           | `0.7`                                   |
| `andes_lfq_feature_fdr`    | `--lfq-feature-fdr`                          | `0.05`                                  |
| `andes_mbr`                | runs the library, extract and finalize steps | `true` (andes's own default is `false`) |
| `andes_mbr_library_fdr`    | `--mbr-library-fdr`                          | `0.01`                                  |
| `andes_quant_feature_fdr`  | `--quant-feature-fdr`                        | `off`                                   |
| `andes_predicted_features` | `--predicted-features`                       | `none` until measured; then `rt`        |
| `andes_extra_args`         | extra flags passed to `andes search`         | `''`                                    |

`andes_glyco` requires `andes_glyco_species` or `andes_glyco_glycan_db`; the schema should
enforce this.

`andes_extra_args` is how the exploratory and agentic modes reach the rest of andes's
options: an agent such as Claude can read andes's `--help` and documentation, choose the
options for a dataset and a question, and pass them here without a pipeline release.

### 5.3 Protein quantification parameters (mokume)

| Parameter                      | mokume flag                                     | Default for LFQ                  | Default for TMT |
| ------------------------------ | ----------------------------------------------- | -------------------------------- | --------------- |
| `protein_quant_method`         | `--quant-method`                                | `maxlfq`                         | `abd`           |
| `protein_min_unique_peptides`  | `--min-unique`                                  | `2`                              | `2`             |
| `protein_min_peptide_length`   | `--min-aa`                                      | `7`                              | `7`             |
| `protein_remove_contaminants`  | `--remove-contaminants` / `--keep-contaminants` | `true`                           | `true`          |
| `protein_run_normalization`    | `--run-normalization`                           | not set (mokume: `median`)       | not set         |
| `protein_sample_normalization` | `--sample-normalization`                        | not set (mokume: `globalmedian`) | not set         |
| `protein_irs`                  | `--irs`                                         | –                                | `false`         |
| `protein_irs_reference_regex`  | `--irs-reference-regex`                         | –                                | mokume default  |
| `protein_impute`               | `--impute`                                      | `false`                          | `false`         |
| `protein_extra_args`           | other mokume flags                              | `''`                             | `''`            |

The LFQ and TMT defaults that differ are set in a config file per experiment type, not as
separate parameters.

Normalisation is left to mokume's defaults on purpose. The feature parquet carries andes's
raw quantity; andes's own normalised value (`quantity_normalized`) is only written to
`mbr.tsv`. See 8.1.

## 6. Outputs

| Output                                | Source                                   |
| ------------------------------------- | ---------------------------------------- |
| `quantms.feature.parquet`             | andes (finalize for LFQ, search for TMT) |
| `psms.parquet` and Percolator results | andes, Percolator                        |
| `mbr.tsv`, `mbr_validation.tsv`       | andes quant finalize (LFQ with MBR)      |
| `proteoforms.tsv`                     | andes (proteoform modes)                 |
| protein table                         | mokume                                   |
| pmultiqc report                       | pmultiqc                                 |

## 7. Containers

- **andes and Percolator:** one image, `ghcr.io/bigbio/andes:<version>` plus a Singularity
  image `oras://ghcr.io/bigbio/andes-sif:<version>`, built by quantms-containers
  (`andes-<version>/Dockerfile`). It installs the andes Linux release archive, which reads
  Thermo `.raw` through a bundled .NET 8 runtime and timsTOF `.d` natively, and Percolator
  3.09.
  - Base image `ubuntu:24.04`: the andes release binary needs glibc 2.39.
  - `procps` is installed; Nextflow's task wrapper needs `ps`.
  - Percolator's `.deb` does not declare `libgomp1` and `libboost-filesystem1.83.0`; the
    image installs them, and the build fails if either binary has a missing library.
  - No `ENTRYPOINT` and no non-root `USER`, so Nextflow can run its own command.
  - Tested with the 0.3.0 archive under Docker and under Singularity with `--no-home
--containall`: a Thermo `.raw` search and Percolator on its output ran in both, with
    the same results.
- The predicted-feature models ship with andes: the native retention-time model is pure
  Rust, with its weights in the release next to the scoring models; fragment-intensity
  models need the ONNX runtime in the image.
- **mokume:** the released `mokume` package image.
- **sdrf-pipelines and pmultiqc:** the images quantms uses.

## 8. Open questions

1. **Normalisation.** andes's run normalisation (`quantity_normalized`, a level per run
   plus an RT-local term) improved MBR CVs in benchmarks, but it is not in the feature
   parquet. Options: andes writes it to the parquet and mokume reads it with its
   normalisation off; or keep mokume's median normalisation. Decide by comparing CVs and
   fold-change accuracy, broken down by abundance.
2. **Transferred values.** mokume ignores `evidence_type`, so MBR-transferred quantities
   weigh the same as identified ones. Decide whether that is acceptable or mokume needs an
   option.
3. **andes protein inference.** To be built in andes: the grouping rule (parsimony, razor
   peptides), the protein q-value (picked protein FDR), and the interface (a subcommand over
   the merged Percolator result). It writes the groups into the feature parquet, so mokume's
   `--min-unique` counts peptides unique to a group, not to one accession.
4. **Variant FDR with one external Percolator.** andes builds knockoff variants and applies
   the class and knockoff FDR only when it rescores itself (`--rescore`). With one
   Percolator over all runs, variant PSMs get neither. Options: an andes subcommand that
   applies the class and knockoff FDR to the merged Percolator result, or per-run
   `--rescore` in the proteoform mode.
5. **TMT filtering.** Confirm whether `quantms.feature.parquet` from `andes search --tmt`
   is filtered by FDR, and how its rows join to the merged Percolator result. The
   quantification gate (purity, PEP, S/N) applies to `tmt.tsv` and the parquet alike from
   0.4.0.
6. **MBR off.** Confirm the LFQ path when `andes_mbr = false`.
7. **Shared index.** Confirm how a prebuilt fragment index is passed to `andes search`.
8. **pmultiqc input.** pmultiqc looks for `*.psm.parquet`, which does not match andes's
   `psms.parquet`. Confirm what pmultiqc needs from the andes outputs.
9. **Modification notation.** andes writes stacked modifications (for example
   Carbamidomethyl and a TMT tag on an N-terminal Cys) by mass, without UniMod accessions,
   and places N-terminal modifications on the first residue. A follow-up andes change will
   write them in ProForma terminal notation with an accession per part.
10. **iTRAQ labels.** andes writes `ITRAQ114`-style labels; mokume matches labels to the
    SDRF exactly after lowercasing. Check the SDRF label spelling matches.
11. **Predicted features.** Measure the search time and the yield of
    `--predicted-features rt` (then `rt,ms2`) on Astral, TMT, UPS1 and timsTOF before
    making `rt` the default. Fragment-intensity models: MS2PIP first (Apache-2.0); others
    only once their weight licence is clear.
12. **Decoys in mokume.** mokume drops decoy rows through its contaminant patterns. Confirm
    the pipeline passes a pattern that matches andes decoys (`XXX_`, or the database's own
    decoy string).

## 9. Testing

- `conf/test.config`: a small LFQ dataset with an SDRF.
- `conf/test_tmt.config`: a small TMT dataset.
- nf-test for each module, checking the andes command line built from the parameters (mods
  file, tolerances, decoy flags).
- Benchmarks against quantms on the same datasets: PSMs at 1% FDR, protein counts, and
  quantification accuracy and CV broken down by abundance.

## 10. Implementation order

1. andes container (andes and Percolator, Docker and Singularity).
2. Input check from SDRF, and the parameter schema.
3. LFQ: index, search, Percolator, library, extract, finalize.
4. andes protein inference, mokume and pmultiqc.
5. TMT.
6. Exploratory modes: glyco, proteoforms, chimeric and refinement.
7. Predicted features in the search (retention time, then fragment intensities).
8. Resolve open questions 8.1–8.4 and 8.11 with benchmark runs before fixing defaults.
