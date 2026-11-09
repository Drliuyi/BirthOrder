association_sample <- function(d, y, exposure, sex_rule = "both", cap = Inf, full = NULL) {
  if (is.null(full)) {
    assert(length(y) == nrow(d) && all(is.na(y) | y %in% 0:1), "Invalid disease states")
    d$y <- y
    variables <- c(exposure, "age", "sex", "ethnicity", "center", "assessment_year", if (exposure != "birth_two_child") "S")
    keep <- is.finite(d[[exposure]]) & complete.cases(d[variables]) & !is.na(d$y)
    if (sex_rule != "both") keep <- keep & as.character(d$sex) == if (sex_rule == "female") "0" else "1"
    d <- droplevels(d[keep, , drop = FALSE])
    terms <- c(exposure, "age", "sex", "ethnicity", "center", "assessment_year", if (exposure != "birth_two_child") "factor(S)")
    constants <- vapply(setdiff(terms, exposure), function(x) length(unique(d[[if (x == "factor(S)") "S" else x]])) < 2, logical(1))
    formula <- stats::reformulate(setdiff(terms, names(constants)[constants]), response = "y")
  } else { d <- full$data; formula <- full$formula }
  d <- droplevels(d[is.finite(d$S) & d$S <= cap, , drop = FALSE])
  list(data = d, formula = formula)
}

logistic_fit <- function(sample, exposure) {
  d <- sample$data; f <- sample$formula
  row <- data.frame(N = nrow(d), cases = sum(d$y == 1), controls = sum(d$y == 0),
    beta = NA_real_, SE = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_,
    status = "NOT_ESTIMABLE", warning = "", reason = "", formula = paste(deparse(f), collapse = " "))
  if (row$cases < 20 || row$controls < 20 || length(unique(d[[exposure]])) < 2) return(row)
  warnings <- character()
  result <- tryCatch(withCallingHandlers({
    X <- stats::model.matrix(f, d)
    assert(all(is.finite(X)) && qr(X)$rank == ncol(X), "Rank-deficient design; no silent variable deletion")
    fit <- stats::glm(f, data = d, family = stats::binomial(), na.action = stats::na.fail)
    cf <- stats::coef(summary(fit))
    assert(isTRUE(fit$converged) && all(is.finite(stats::coef(fit))), "Nonconverged/nonfinite Logistic model")
    assert(exposure %in% rownames(cf) && all(is.finite(cf[exposure, ])) && cf[exposure, 2] > 0, "Unidentifiable exposure")
    list(coefficient = cf[exposure, ], fit = fit)
  }, warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }), error = identity)
  row$warning <- paste(unique(warnings), collapse = "; ")
  if (inherits(result, "error")) { row$reason <- conditionMessage(result); return(row) }
  b <- unname(result$coefficient[1]); se <- unname(result$coefficient[2])
  row$beta <- b; row$SE <- se; row$P <- unname(result$coefficient[4])
  row$OR <- exp(b); row$CI_low <- exp(b - stats::qnorm(.975) * se); row$CI_high <- exp(b + stats::qnorm(.975) * se)
  row$status <- if (length(warnings)) "OK_WARNING_DESCRIPTIVE" else "OK_DESCRIPTIVE"
  attr(row, "fit") <- result$fit
  row
}

phewas <- function(d, disease, dictionary, restrictions) {
  columns <- setdiff(names(disease), "eid")
  assert(length(columns) == 1582 && !anyDuplicated(columns), "Frozen PheWAS requires all 1582 disease positions")
  ix <- match(as.character(d$eid), as.character(disease$eid)); present <- !is.na(ix)
  d <- d[present, , drop = FALSE]; disease <- disease[ix[present], , drop = FALSE]
  mi <- match(columns, as.character(dictionary$phecode)); ri <- match(columns, as.character(restrictions$phecode))
  assert(!anyNA(mi) && !anyNA(ri), "Unmapped PheCode metadata")
  sex <- rep("both", length(columns))
  sex[tolower(as.character(restrictions$female_only[ri])) == "true"] <- "female"
  sex[tolower(as.character(restrictions$male_only[ri])) == "true"] <- "male"
  exposures <- c("birth_two_child", "birth_first_last", "birth_rank", "birth_relative")
  inherited_M <- c(birth_first_last = 1362L, birth_rank = 1419L, birth_relative = 1419L)
  result <- list()
  for (exposure in exposures) {
    # Pre-fit family is determined on the original all-family population.
    planned <- vapply(seq_along(columns), function(i) {
      z <- association_sample(d, disease[[columns[i]]], exposure, sex[i])$data
      sum(z$y == 1) >= 20 && sum(z$y == 0) >= 20 && length(unique(z[[exposure]])) > 1
    }, logical(1))
    M <- sum(planned)
    if (exposure != "birth_two_child") assert(M == inherited_M[[exposure]], paste("Frozen test family changed:", exposure, M))
    for (population in c("All families", if (exposure != "birth_two_child") "Sibship 2-8")) {
      rows <- lapply(seq_along(columns), function(i) {
        full <- association_sample(d, disease[[columns[i]]], exposure, sex[i])
        sample <- if (population == "All families") full else association_sample(d, NULL, exposure, cap = 8, full = full)
        row <- logistic_fit(sample, exposure); attr(row, "fit") <- NULL
        cbind(data.frame(population = population, exposure = exposure, phecode = columns[i],
          description = dictionary$description[mi[i]], category = dictionary$group[mi[i]], planned = planned[i]), row)
      })
      block <- do.call(rbind, rows); block$BH_q <- bh_fixed(block$P, M); block$family_n <- M
      block$Bonferroni_P <- pmin(1, block$P * M)
      result[[paste(population, exposure)]] <- block
    }
  }
  do.call(rbind, result)
}

adjustment_models <- function(d, codes) {
  m1 <- c("birth_relative", "age", "sex", "ethnicity", "center", "assessment_year", "factor(S)")
  m2 <- c(m1, "education", "household_income", "Townsend")
  m3 <- c(m2, "smoking_status", "alcohol_status", "BMI", "moderate_activity_days", "walking_days")
  forms <- lapply(list(M1 = m1, M2 = m2, M3 = m3), stats::reformulate, response = "y")
  assert(all(unique(unlist(lapply(forms, all.vars))) %in% c(names(d), "y")), "Missing encoded adjustment covariates")
  rows <- curves <- list()
  for (outcome in names(codes)) {
    a <- d; a$y <- a[[codes[[outcome]]]]
    common <- complete.cases(a[all.vars(forms$M3)])
    for (label in c("M1_original", "M1", "M2", "M3")) {
      model <- if (label == "M1_original") "M1" else label
      keep <- if (label == "M1_original") complete.cases(a[all.vars(forms$M1)]) else common
      z <- logistic_fit(list(data = droplevels(a[keep, , drop = FALSE]), formula = forms[[model]]), "birth_relative")
      fit <- attr(z, "fit"); attr(z, "fit") <- NULL
      rows[[paste(outcome, label)]] <- cbind(outcome = outcome, model = label, z)
      if (label == "M3" && !is.null(fit)) curves[[outcome]] <- standardized_curve(fit, a[keep, , drop = FALSE], outcome)
    }
  }
  results <- do.call(rbind, rows); results$BH_q <- NA_real_
  common <- results$model != "M1_original"; results$BH_q[common] <- bh_fixed(results$P[common], 15L)
  list(models = results, curves = do.call(rbind, curves))
}

standardized_curve <- function(fit, d, outcome) {
  covariance <- stats::vcov(fit)
  predictions <- lapply(seq(0, 1, .05), function(b) {
    x <- d; x$birth_relative <- b
    X <- stats::model.matrix(stats::delete.response(stats::terms(fit)), x)
    p <- stats::plogis(as.numeric(X %*% stats::coef(fit)))
    gradient <- colMeans(X * (p * (1 - p)))
    list(mean = mean(p), gradient = gradient)
  })
  do.call(rbind, lapply(seq_along(predictions), function(i) {
    delta <- predictions[[i]]$mean - predictions[[1]]$mean
    gradient <- predictions[[i]]$gradient - predictions[[1]]$gradient
    se <- sqrt(as.numeric(t(gradient) %*% covariance %*% gradient))
    data.frame(outcome = outcome, birth_relative = (i - 1) * .05, difference_pp = 100 * delta,
      CI_low_pp = 100 * (delta - qnorm(.975) * se), CI_high_pp = 100 * (delta + qnorm(.975) * se))
  }))
}

incident_cox <- function(baseline, risks, config) {
  require_packages("survival")
  assert(is.data.frame(risks) && all(c("eid", "outcome", "entry", "exit", "event", "eligible") %in% names(risks)),
    "Incident input must contain authoritative dated risk sets, not PheWAS disease-history states")
  rows <- lapply(names(config$disease_codes), function(outcome) {
    r <- risks[risks$outcome == outcome & risks$eligible %in% TRUE, , drop = FALSE]
    assert(!anyDuplicated(as.character(r$eid)), "Duplicate disease risk-set IDs")
    i <- match(as.character(r$eid), as.character(baseline$eid)); r <- r[!is.na(i), , drop = FALSE]
    d <- baseline[i[!is.na(i)], , drop = FALSE]
    d$entry <- r$entry; d$exit <- r$exit; d$event <- r$event
    assert(all(is.finite(d$entry) & is.finite(d$exit) & d$exit > d$entry) && all(d$event %in% 0:1), "Invalid full-follow-up risk set")
    fit <- survival::coxph(survival::Surv(entry, exit, event) ~ birth_relative + age + sex + ethnicity + center + assessment_year + factor(S),
      data = d, ties = config$cox_ties, na.action = stats::na.omit)
    cf <- summary(fit)$coefficients["birth_relative", ]; ci <- summary(fit)$conf.int["birth_relative", ]
    data.frame(outcome = outcome, N = fit$n, events = fit$nevent, HR = ci[1], CI_low = ci[3], CI_high = ci[4],
      P = cf[ncol(summary(fit)$coefficients)], followup_cap = Inf)
  })
  result <- do.call(rbind, rows); result$BH_q <- bh_fixed(result$P, 5L); result
}
