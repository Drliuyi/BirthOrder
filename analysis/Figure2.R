# Figure 2: disease-history PheWAS. No incident-outcome substitution.
analyse_Figure2 <- function(config, code_root) {
  d <- read_baseline(config, primary_only = FALSE)
  disease <- read_diseases(config)
  dictionary <- read_table(file.path(code_root, "reference", "diseases.csv"))
  restrictions <- read_table(file.path(code_root, "reference", "sex_restrictions.csv"))
  list(phewas = phewas(d, disease, dictionary, restrictions))
}
