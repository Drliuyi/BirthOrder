# Change paths here or set the corresponding environment variables.
birth_config <- function() {
  env <- function(name, default) Sys.getenv(name, unset = default)
  drive <- env("DIR0", if (.Platform$OS.type == "windows") "D:/" else "/mnt/d")
  scripts <- env("SCRIPT_ROOT", file.path(drive, "scripts"))
  phenotype <- env("PHEDIR", file.path(drive, "data", "ukb", "phe"))
  output <- env("UKB_OUT", file.path(drive, "analysis", "ukb"))
  frozen <- env("BIRTH_FROZEN", file.path(output, "birth_order", "frozen"))
  list(paths = list(
    phenotype = phenotype, scripts = scripts,
    helpers = env("HELPER_DIR", file.path(scripts, "0f")),
    all = env("ALL_RDS", file.path(phenotype, "Rdata", "all.rds")),
    phewas = env("PHEWAS_RDS", file.path(phenotype, "Rdata", "PheWAS.rds")),
    metabolites = env("METABOLITE_RDS", file.path(phenotype, "Rdata", "met.rds")),
    proteins = env("PROTEIN_RDS", file.path(phenotype, "Rdata", "prot.rds")),
    withdrawals = env("BIRTH_WITHDRAWALS", file.path(phenotype, "common", "withdrawal.csv")),
    covariates = env("BIRTH_COVARIATES", file.path(phenotype, "Rdata", "birth_covariates.rds")),
    incident = env("BIRTH_INCIDENT", file.path(phenotype, "Rdata", "birth_incident.rds")),
    frozen = frozen, results = env("BIRTH_RESULTS", file.path(output, "birth_order", "results"))),
    workbook = env("BIRTH_WORKBOOK", file.path(frozen, "border_supplementary_table.xlsx")),
    source = env("BIRTH_SOURCE", file.path(output, "birth_order", "source")),
    primary_covariates = c("age", "sex", "ethnicity", "center", "assessment_year"),
    disease_codes = c(DM = "250", T2D = "250.2", IHD = "411", Angina = "411.3", CAS = "411.4"),
    primary_sibship = 2:8, cox_ties = "efron", risk_years = 11, mc_draws = 50L,
    bootstrap_ids = 1:100,
    joint_models = c("S44_BASE11", "S44_NMR_REF11", "S44_NMR_PLUS40", "S44_PROT_REF11", "S44_PROT_PLUS15"),
    colours = list(positive = "#B65F61", negative = "#5589AA", exposure = "#1B6866",
      population = c("11 Risk Factors" = "#575757", "+29 Metabolites" = "#A96D31", "+4 Proteins" = "#345F8C"),
      domain = c(SMOKING = "#C29AA1", BODY_COMPOSITION = "#A5C0BE", CLINICAL_LIPIDS = "#B5C0CF",
        GLYCAEMIA = "#D9B88A", NMR_ALL = "#8BAD9E", CANDIDATE_PROTEINS = "#A894B7")))
}
