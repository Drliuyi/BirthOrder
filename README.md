# Birth Order

Analysis and plotting code for *Relative birth order shows opposing associations with ischaemic heart disease and type 2 diabetes*.

**Code only.** No UK Biobank records, individual identifiers, measurements, dates of individual events, sample lists, Bootstrap draw files, fitted models, results, figures or manuscripts are distributed. The four CSV files in `reference/` are variable/phenotype dictionaries, not observations. Access to UK Biobank data must be obtained independently.

## Organisation

Analysis and plotting are separate. Filenames follow the main and supplementary figures.

| Entry | Scope |
|---|---|
| `analysis/Figure2*` | PheWAS, relative/absolute/first–last/two-child definitions, LASI and published Swedish trend methods |
| `analysis/Figure3*` | Single mediation, the 121-candidate screen, paired disease differences, representative selection and refinement |
| `analysis/Figure4*` | Joint/domain mediation, BODY conditional donors, metabolite composition and final protein Y estimation |
| `analysis/Figure5*` | Bulk expression, single-cell QC, embedding and pseudobulk analysis |
| `analysis/Table1.R` | Baseline descriptions |
| `analysis/SupplementaryFigure2*` | Adjustment, full-follow-up incident Cox, subgroups and standardized probability curves |
| `analysis/Supplementary*` | Risk-factor/omics/CMR associations, shape and HbA1c/penalty sensitivities |
| `figures/` | Quantitative renderers and original panel-rendering function libraries |
| `reference/` | Variable definitions, frozen model membership and code index only |

`reference/modules.json` lists the method libraries and their functions. Single-mediator and joint-mediator analyses use separate estimators.

## Data locations

```text
DIR0/
  scripts/BirthOrder/          analysis code
  scripts/0f/                 shared helper functions
  data/ukb/phe/Rdata/          all.rds, PheWAS.rds, met.rds, prot.rds
  analysis/ukb/birth_order/    private derived inputs, frozen results and new outputs
```

Set paths in `config.R` or with the environment variables in `paths.sh`. The directory layout above is an example; all locations are configurable. Keep data and outputs outside this repository. Analyses start from prepared baseline/PheWAS objects and omics matrices. Incident analyses also require event dates, risk sets and the relevant covariates.

```sh
source paths.sh
bash run.sh check
Rscript --vanilla test.R "$PWD"
bash run.sh methods Figure4
```

The checks use syntax and synthetic inputs only. They do not fit UK Biobank models or start Bootstrap jobs. No packages are installed automatically.

## Usage

The command-line entry point supports baseline descriptions, PheWAS and adjustment/incident analyses:

```sh
bash run.sh analyse Table1 "$BIRTH_RESULTS/new_table1"
bash run.sh analyse Figure2 "$BIRTH_RESULTS/new_phewas"
bash run.sh analyse SupplementaryFigure2 "$BIRTH_RESULTS/new_adjustment_cox"
bash run.sh draw Figure4 "$BIRTH_SOURCE/Figure4.xlsx" "$BIRTH_RESULTS/new_figure4"
```

Use a new output directory for each run. Plotting commands read existing estimates; they do not refit statistical models.

Set `BIRTH_PYTHON` to a Python environment with the plotting dependencies installed. The R plotting entry points use this interpreter; otherwise they use `python3`.

The mediation and expression analyses are function libraries. R libraries can be loaded using `analysis/private_context.R`; Python libraries require their corresponding input context and dependencies. Inputs and integrity checks use `BIRTH_PRIVATE_*` settings. Load each estimator in its own environment. The `draw` commands produce quantitative panels rather than assembled manuscript pages.

The repository does not yet provide a single command that regenerates every manuscript result from the four input objects. Mediation input/model adapters, dated risk-set preparation and final composite-figure assembly still require stage-specific inputs. Figure 1 illustrations are supplied separately.

## Figure panels

Standalone renderers use explicit input paths and refuse existing outputs. Required columns are described at the top of each script.

| Panel | Renderer | Inputs |
|---|---|---|
| Figure 3d | `figures/Figure3_contrasts.py` | Saved six-contrast estimates and a-reference directions |
| Figure 4b | `figures/Figure4_heatmap.py` | Saved domain estimates/q values and joint estimates/intervals/proportions |
| Figure 5a | `figures/Figure5_bulk.R` | Expression values and saved bracket annotations |
| Figure 5b | `figures/Figure5_UMAP.py` | Coordinates, cell labels and reference geometry |
| Figure 5c | `figures/Figure5_features.py` | Coordinates and normalized gene expression; one shared 0–5 colour scale |
| Figure 5d | `figures/Figure5_donors.py` | Donor-level expression and saved significance labels; `--layout horizontal` for side-by-side cell types |

```sh
"${BIRTH_PYTHON:-python3}" "$SCRIPT_ROOT/BirthOrder/figures/Figure3_contrasts.py" --values "$BIRTH_SOURCE/contrasts.csv" --a-directions "$BIRTH_SOURCE/directions.csv" --out "$BIRTH_RESULTS/new_figure3d"
"${BIRTH_PYTHON:-python3}" "$SCRIPT_ROOT/BirthOrder/figures/Figure4_heatmap.py" --heatmap "$BIRTH_SOURCE/domains.csv" --joint "$BIRTH_SOURCE/joint.csv" --out "$BIRTH_RESULTS/new_figure4b"
```

Adjust the script location if the repository is stored elsewhere. Keep inputs and outputs outside the repository. Figure 3d's ribbon colour is a derived contrast-direction code, not a fitted Y-model coefficient. Figure 4b uses the saved model-by-effect BH families, including DEP in each correction family; it does not calculate a new whole-figure correction. Final multi-panel assembly and illustrations are separate from the quantitative renderers.

## Statistical separation

The Cox models use full available follow-up; mediation effects are standardized cumulative incidence contrasts at 11 years. Joint analyses use 50 inner Monte Carlo draws and 100 bootstrap replicates, with the original-sample estimate separate. The selected set contains 11 non-omics indicators, 29 metabolites and four proteins. Candidate screening and representative selection follow the 121 → 64 → 44 → 18 sequence implemented in `analysis/mediation.R` and the Figure 3 libraries.

Protein models include centre-stratified cause-specific baselines and the `survival::coxph` ridge implementation. Preserve matching inputs, model specifications and random streams when combining estimates. Plotting uses saved P values, BH-adjusted q values and confidence intervals.

## Publication safeguards

```sh
python3 -B privacy_check.py
git config core.hooksPath hooks
```

The checker inspects an explicit, hash-bound file allowlist, text content, dictionary schemas and staged/committed blobs. It rejects unreviewed files, binary data, private directories, common credential patterns and literal individual identifiers. GitHub Actions repeats the committed-tree check. These safeguards supplement, not replace, manual source review. Never add private inputs or generated outputs, even in a private repository.

Required software depends on the stage: R with survival, splines, nnet, digest, jsonlite and data.table; plotting additionally uses ggplot2, readxl, svglite and related packages. Python plotting/expression stages use numpy, pandas, scipy, matplotlib and the single-cell packages imported by each module. No packages are installed automatically. Third-party software retains its own licence terms.
