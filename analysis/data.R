clean_count <- function(x) {
  x <- numeric_value(x)
  x[!is.finite(x) | x < 0 | x != floor(x)] <- NA_real_
  x
}
birth_order <- function(d) {
  needed <- c("eid", "cnt_older", "cnt_brother", "cnt_sister", "adopted")
  assert(all(needed %in% names(d)), paste("Missing baseline fields:", paste(setdiff(needed, names(d)), collapse = ", ")))
  assert(!anyNA(d$eid) && !anyDuplicated(as.character(d$eid)), "Missing or duplicate baseline IDs")
  older <- clean_count(d$cnt_older); brothers <- clean_count(d$cnt_brother); sisters <- clean_count(d$cnt_sister)
  adopted <- numeric_value(d$adopted); adopted[!adopted %in% c(0, 1)] <- NA_real_
  d$S <- brothers + sisters + 1; d$rank <- older + 1
  valid <- !is.na(older) & !is.na(brothers) & !is.na(sisters) & older <= brothers + sisters & (is.na(adopted) | adopted != 1)
  d$birth_relative <- ifelse(valid & d$S >= 2, older / (d$S - 1), NA_real_)
  d$birth_rank <- ifelse(valid & d$S >= 2, d$rank, NA_real_)
  extreme <- valid & d$S >= 2 & (d$rank == 1 | d$rank == d$S)
  d$birth_first_last <- ifelse(extreme, as.numeric(d$rank == d$S), NA_real_)
  d$birth_two_child <- ifelse(valid & d$S == 2, older, NA_real_)
  d$birth_group <- factor(ifelse(valid & d$S >= 2, ifelse(d$rank == 1, "First", ifelse(d$rank == d$S, "Last", "Middle")), NA_character_), levels = c("First", "Middle", "Last"))
  d
}
primary_covariates <- function(d) {
  need <- c("age", "sex", "ethnicity", "center", "date_attend")
  assert(all(need %in% names(d)), paste("Missing primary covariates:", paste(setdiff(need, names(d)), collapse = ", ")))
  d$age <- numeric_value(d$age); d$age[!is.finite(d$age)] <- NA_real_
  sex <- numeric_value(d$sex); sex[!sex %in% 0:1] <- NA_real_
  d$sex <- factor(sex, levels = 0:1)
  e <- numeric_value(d$ethnicity); ethnicity <- rep(NA_character_, nrow(d))
  categories <- list(White = c(1,1001:1003), Mixed = c(2,2001:2004), `Asian or Asian British` = c(3,3001:3004),
    `Black or Black British` = c(4,4001:4003), Chinese = 5, `Other ethnic group` = 6)
  for (name in names(categories)) ethnicity[e %in% categories[[name]]] <- name
  d$ethnicity <- factor(ethnicity, levels = names(categories))
  center <- as.character(d$center); center[is.na(d$center) | !nzchar(center) | center %in% c("-1","-3")] <- NA_character_
  d$center <- factor(center)
  d$assessment_year <- factor(format(as.Date(d$date_attend), "%Y"))
  d
}
read_baseline <- function(config, primary_only = TRUE) {
  p <- config$paths$all; assert(file.exists(p), paste("Missing all.rds:", p))
  d <- primary_covariates(birth_order(readRDS(p)))
  if (primary_only) d <- d[is.finite(d$birth_relative) & d$S %in% config$primary_sibship, , drop = FALSE]
  d
}
read_diseases <- function(config) {
  p <- config$paths$phewas; assert(file.exists(p), paste("Missing PheWAS.rds:", p))
  x <- readRDS(p); assert(is.data.frame(x) && "id" %in% names(x), "PheWAS.rds must have id and disease columns")
  assert(!anyNA(x$id) && !anyDuplicated(as.character(x$id)) && !anyDuplicated(names(x)), "Invalid PheWAS IDs/columns")
  names(x)[names(x) == "id"] <- "eid"
  y <- x[setdiff(names(x), "eid")]
  assert(all(vapply(y, function(v) all(is.na(v) | v %in% 0:1), logical(1))), "Disease states must remain 0/1/NA")
  x
}

adjustment_covariates <- function(d, config) {
  needed <- c("education", "household_income", "Townsend", "smoking_status", "alcohol_status", "BMI", "moderate_activity_days", "walking_days")
  # These are the actual M1/M2/M3 coding rules, not a new covariate selection.
  if (!"education" %in% names(d) && "edu" %in% names(d)) {
    hierarchy <- c("1"="Degree", "5"="NVQ_HND_HNC", "6"="Professional", "2"="A_AS", "3"="O_GCSE", "4"="CSE", "-7"="None")
    labels <- vapply(as.character(d$edu), function(s) {
      if (is.na(s) || !nzchar(trimws(s))) return(NA_character_)
      s <- gsub('[\"\047]', '', chartr('[]', '  ', s))
      x <- unique(suppressWarnings(as.numeric(strsplit(trimws(s), '[|,;[:space:]]+', perl=TRUE)[[1]])))
      if (length(x)==1L && !is.na(x) && x==7) x <- -7
      if (!length(x) || anyNA(x) || any(!x %in% c(1:6,-7)) || (-7 %in% x && length(x)>1L)) return(NA_character_)
      unname(hierarchy[which(names(hierarchy) %in% as.character(x))[1]])
    }, character(1))
    d$education <- factor(labels, levels=c("None","CSE","O_GCSE","A_AS","Professional","NVQ_HND_HNC","Degree"))
  }
  categorical <- function(x, levels) { x<-numeric_value(x);x[!x %in% levels]<-NA_real_;factor(x,levels=levels) }
  if (!"household_income" %in% names(d) && "inc" %in% names(d)) d$household_income<-categorical(d$inc,1:5)
  if (!"Townsend" %in% names(d) && "tdi" %in% names(d)) {d$Townsend<-numeric_value(d$tdi);d$Townsend[!is.finite(d$Townsend)]<-NA_real_}
  if (!"smoking_status" %in% names(d) && "smoke_status" %in% names(d)) d$smoking_status<-categorical(d$smoke_status,0:2)
  if (!"BMI" %in% names(d) && "bmi" %in% names(d)) d$BMI<-numeric_value(d$bmi)
  for (pair in list(c("moderate_activity_days","p884_i0"),c("walking_days","p864_i0"))) {
    if (!pair[1] %in% names(d) && pair[2] %in% names(d)) {
      x<-numeric_value(d[[pair[2]]]);x[!is.finite(x)|!x %in% 0:7]<-NA_real_;d[[pair[1]]]<-x
    }
  }
  missing <- setdiff(needed, names(d))
  if (length(missing)) {
    assert(file.exists(config$paths$covariates), paste("Required covariates absent from all.rds; supply the frozen baseline supplement:",paste(missing,collapse=", ")))
    supplement <- readRDS(config$paths$covariates)
    assert(all(c("eid",missing) %in% names(supplement)), "Baseline supplement lacks required encoded fields")
    assert(all(as.character(d$eid) %in% as.character(supplement$eid)), "Supplement ID coverage differs from original baseline roster")
    d <- join_by_id(d, supplement[c("eid",missing)])
  }
  if ("alcohol_status" %in% names(d) && !is.factor(d$alcohol_status)) d$alcohol_status<-categorical(d$alcohol_status,0:2)
  assert(all(needed %in% names(d)), "Incomplete adjustment variables")
  expected_levels <- list(education=c("None","CSE","O_GCSE","A_AS","Professional","NVQ_HND_HNC","Degree"),
    household_income=as.character(1:5),smoking_status=as.character(0:2),alcohol_status=as.character(0:2))
  for (name in names(expected_levels)) {
    x<-as.character(d[[name]]);allowed<-expected_levels[[name]]
    assert(all(is.na(x)|x %in% allowed),paste("Unknown frozen categorical coding:",name))
    d[[name]]<-factor(x,levels=allowed)
  }
  assert(all(is.na(d$BMI)|(is.finite(d$BMI)&d$BMI>0)),"BMI requires input review; no silent trimming")
  d
}
