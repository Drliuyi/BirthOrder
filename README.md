# Birth Order

Analysis source for *Relative birth order shows divergent associations with ischaemic heart disease and type 2 diabetes*.

**Code only.** No UK Biobank records, individual identifiers, measurements, dates of individual events, sample lists, Bootstrap draw files, fitted models, results, figures or manuscripts are distributed. The four CSV files in `reference/` are variable/phenotype dictionaries, not observations. Access to UK Biobank data must be obtained independently.

## Organisation

The shallow structure follows [migration](https://github.com/zhaoran7/migration): analysis and rendering are separate, with filenames tied to the manuscript rather than a chain of dated runs.

| Entry | Scope |
|---|---|
| `analysis/Figure2*` | PheWAS, relative/absolute/first–last/two-child definitions, LASI and published Swedish trend methods |
| `analysis/Figure3*` | Adjustment, full-follow-up incident Cox, subgroup analysis and standardized probability curves |
| `analysis/Figure4*` | Single mediation, the 121-candidate screen, paired disease differences, representative selection and refinement |
| `analysis/Figure5*` | Joint/domain mediation, BODY conditional donors, metabolite composition and final protein Y estimation |
| `analysis/Figure6*` | Bulk expression, single-cell QC, embedding and pseudobulk analysis |
| `analysis/Table1.R` | Baseline descriptions |
| `analysis/Supplementary*` | Risk-factor/omics/CMR associations, shape and HbA1c/penalty sensitivities |
| `figures/` | Quantitative renderers and original panel-rendering function libraries |
| `reference/` | Variable definitions, frozen model membership and code index only |

`reference/modules.json` lists the extracted method modules and callable functions. Separate single-screen, refinement and joint programs are not combined into one estimator. Obsolete production queues and remote deployment/monitoring scripts are not included.

## Huang-lab paths

```text
DIR0/                         D:/ (native Windows), /mnt/d (WSL)
  scripts/BirthOrder/          this repository; code only
  scripts/0f/                 existing shared helpers; left unchanged
  data/ukb/phe/Rdata/          all.rds, PheWAS.rds, met.rds, prot.rds
  analysis/ukb/birth_order/    private derived inputs, frozen results and new outputs
```

`config.R` and `paths.sh` use `DIR0`, `PHEDIR`, `SCRIPT_ROOT`, `HELPER_DIR` and `UKB_OUT`. All paths may be overridden. Do not put inputs or outputs inside this repository. The starting point is the prepared baseline/PheWAS objects and omics matrices; dated outcome/risk-set and additional covariate inputs remain necessary where they are absent from those objects.

```sh
source paths.sh
bash run.sh check
Rscript --vanilla test.R "$PWD"
bash run.sh methods Figure5
```

The checks use syntax and synthetic inputs only. They do not fit UK Biobank models or start Bootstrap jobs. No packages are installed automatically.

## Execution scope

The independent facade supports baseline descriptions, PheWAS and adjustment/incident analyses after the private inputs have been prepared:

```sh
bash run.sh analyse Table1 "$BIRTH_RESULTS/new_table1"
bash run.sh analyse Figure2 "$BIRTH_RESULTS/new_phewas"
bash run.sh analyse Figure3 "$BIRTH_RESULTS/new_adjustment_cox"
bash run.sh draw Figure5 "$BIRTH_SOURCE/Figure5.xlsx" "$BIRTH_RESULTS/new_figure5"
```

Outputs must be new directories. These commands are not run as part of publishing the repository. Fresh-input facade results have not been certified against every frozen historical analysis.

The detailed method modules preserve the original fitted-model/g-computation implementations as explicit function libraries. Packaging removed source-time job dispatch and moved literal private paths and input locks into required `BIRTH_PRIVATE_*` settings. Each R module also retains its script body in `run_in_private_context(context)`, which is never called automatically. An authorized, correctly prepared private context and matching software are required; do not source unrelated estimator modules into a shared global environment. Python modules similarly require their private stage context; extracting functions is not a claim that missing historical inputs have been reconstructed.

**This is not yet a certified one-command pipeline from the four raw input objects to all final results.** The remaining integration is the fresh mediation input/A/Y adapter, dated risk-set preparation and final author-controlled composite-figure assembly. No `all` command silently runs a subset and claims completion. Figure 1 artwork and manual manuscript composites are not analysis code and are not uploaded. Generic `draw` outputs are quantitative panels, not a promise of pixel-identical author layouts.

## Statistical separation

The Cox models use full available follow-up; the mediation estimands use standardized cumulative incidence at 11 years. Single mediation and joint/domain mediation retain separate approved implementations. Joint models use 50 inner Monte Carlo draws and fixed Bootstrap numbers 1–100, with B0 separate. The selected set remains 11 non-omics indicators, 29 metabolites and four proteins. The original 121 → 64 → 44 → 18 selection and its signed PM threshold are retained. No discarded joint-model branch is dispatched by this repository.

The protein methods include centre-stratified cause-specific baselines, the approved ethnicity coding and native `survival::coxph` ridge implementation. Compatibility of inputs, model versions and random streams must be checked in the private runtime; successful draws from different statistical versions must not be mixed. Statistical P/q and confidence intervals are not recomputed merely to format figures.

## Publication safeguards

```sh
python3 -B privacy_check.py
git config core.hooksPath hooks
```

The checker inspects an explicit, hash-bound file allowlist, text content, dictionary schemas and staged/committed blobs. It rejects unreviewed files, binary data, private directories, common credential patterns and literal individual identifiers. GitHub Actions repeats the committed-tree check. These safeguards supplement, not replace, manual source review. Never add private inputs or generated outputs, even in a private repository.

Required software depends on the stage: R with survival, splines, nnet, digest, jsonlite and data.table; plotting additionally uses ggplot2, readxl, svglite and related packages. Python plotting/expression stages use numpy, pandas, scipy, matplotlib and the original single-cell stack. Exact environments must be supplied privately for historical reproduction. Code is shared privately; no open-source licence is asserted on behalf of the authors or third-party package developers.
