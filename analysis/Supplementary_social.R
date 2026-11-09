# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

bi_ses_derive <- function(d) {
    finite <- is.finite(d$tdi)
    cuts <- if (any(finite)) 
        as.numeric(stats::quantile(d$tdi[finite], c(0.2, 0.4, 0.6, 0.8), type = 7, names = FALSE))
    else rep(NA_real_, 4)
    valid <- all(is.finite(cuts)) && !anyDuplicated(cuts)
    d$tdi_quintile <- if (valid) 
        cut(d$tdi, breaks = c(-Inf, cuts, Inf), labels = paste0("Q", 1:5), right = TRUE)
    else factor(rep(NA_character_, nrow(d)), levels = paste0("Q", 1:5))
    attr(d, "tdi_cuts") <- cuts
    attr(d, "tdi_cut_status") <- if (valid) 
        "OK"
    else "UNAVAILABLE_NONUNIQUE_QUANTILES"
    d
}

bi_ses_extension <- function(ctx, d, exposure, extra, code, label) {
    ff <- bi_formula(exposure, extra = extra)
    unavailable <- bi_c_outcome_unavailable(d, exposure)
    if (!is.null(unavailable)) {
        z <- unavailable[rep(1L, 2L), , drop = FALSE]
        z$model <- c("M0_same_sample", label)
        z$formula <- c(paste(deparse(bi_formula(exposure)), collapse = " "), paste(deparse(ff), collapse = " "))
        z$sample_set <- paste0("complete_case_", label)
        z$extra_covariates <- paste(extra, collapse = ";")
        z$source_N <- nrow(d)
        z$union_complete_N <- z$excluded_union_complete_N <- NA_integer_
        return(z)
    }
    same <- bi_sample(d, ff, exposure)$data
    a <- bi_fit(ctx, same, exposure)
    b <- bi_fit(ctx, same, exposure, extra = extra)
    if (!identical(as.character(a$sample_hash), as.character(b$sample_hash))) 
        stop("SES same-sample hash invariant failed")
    a$model <- "M0_same_sample"
    b$model <- label
    z <- bi_bind(list(a, b))
    z$sample_set <- paste0("complete_case_", label)
    z$extra_covariates <- paste(extra, collapse = ";")
    z$source_N <- nrow(d)
    z$union_complete_N <- nrow(same)
    z$excluded_union_complete_N <- nrow(d) - nrow(same)
    z
}

bi_stage_ses <- function(ctx) {
    cfg <- bi_module_config(ctx, "ses")
    adult <- c("bmi", "smoking", "alcohol")
    if (!is.null(cfg$adult_covariates)) 
        adult <- as.character(unlist(cfg$adult_covariates))
    modifiers <- c("tdi", "tdi_quintile", "income", "income_group", "education", "university")
    extensions <- list(M1a = "tdi", M1b = "income", M1c = "education", M1d = c("tdi", "income", "education"), M2 = adult)
    estimates <- interactions <- models <- missing <- list()
    k <- j <- m <- h <- 0L
    for (code in bi_c_codes(ctx)) {
        d <- bi_ses_derive(bi_data(ctx, code))
        for (exposure in bi_c_exposures()) {
            for (modifier in modifiers) {
                key <- paste("ses_interaction", exposure, code, modifier, sep = "__")
                z <- bi_c_job(ctx, key, exposure, function() bi_c_outcome_unavailable(d, exposure) %||% bi_interaction(ctx, 
                  d, exposure, modifier))
                z$modifier <- modifier
                z$inference_scope <- "Post-hoc exploratory modifier-specific BH; not formal whole-PheWAS q"
                j <- j + 1L
                interactions[[j]] <- bi_c_annotate(z, code, "adult_SES_global_interaction", modifier, key, paste0("ses_exploratory_", 
                  modifier))
                if (modifier != "tdi") 
                  for (level in bi_c_modifier_levels(d, modifier)) {
                    skey <- paste("ses_stratum", exposure, code, modifier, level, sep = "__")
                    z <- bi_c_job(ctx, skey, exposure, function() {
                      unavailable <- bi_c_outcome_unavailable(d, exposure)
                      if (!is.null(unavailable)) 
                        return(unavailable)
                      keep <- !is.na(d[[modifier]]) & as.character(d[[modifier]]) == level
                      bi_fit(ctx, d[keep, , drop = FALSE], exposure)
                    })
                    k <- k + 1L
                    estimates[[k]] <- bi_c_annotate(z, code, paste0("adult_SES_", modifier), level, skey)
                  }
            }
            for (label in names(extensions)) {
                key <- paste("ses_extension", exposure, code, label, sep = "__")
                z <- bi_c_job(ctx, key, exposure, function() bi_ses_extension(ctx, d, exposure, extensions[[label]], 
                  code, label))
                if (!"model" %in% names(z)) {
                  z <- z[rep(1L, 2L), , drop = FALSE]
                  z$model <- c("M0_same_sample", label)
                }
                z <- bi_c_annotate(z, code, "adult_SES_same_sample_extension", label, paste(key, z$model, sep = "__"))
                m <- m + 1L
                models[[m]] <- z
            }
            for (field in unique(c("tdi", "income", "education", adult))) {
                h <- h + 1L
                eligible <- !is.na(d[[exposure]])
                present <- field %in% names(d)
                missing[[h]] <- data.frame(trait = code, exposure = exposure, field = field, N = sum(eligible), 
                  missing_N = if (present) 
                    sum(eligible & is.na(d[[field]]))
                  else sum(eligible), status = if (present && any(!is.na(d[[field]]))) 
                    "AVAILABLE"
                  else "UNAVAILABLE_FIELD")
            }
        }
        cuts <- data.frame(trait = code, quantile = c(0.2, 0.4, 0.6, 0.8), value = attr(d, "tdi_cuts"), status = attr(d, 
            "tdi_cut_status"), reference = "Entire frozen baseline cohort before outcome/exposure complete-case filtering; type7; ties are never split")
        bi_save(ctx, paste0("raw/ses/tdi_quantile_definition_", gsub("[.]", "_", code), ".csv"), cuts)
    }
    estimates <- bi_bind(estimates)
    interactions <- bi_adjust(ctx, bi_bind(interactions))
    models <- bi_bind(models)
    bi_save(ctx, "raw/ses/stratified_results.csv", estimates)
    bi_save(ctx, "raw/ses/interactions.csv", interactions)
    bi_save(ctx, "raw/ses/same_sample_models.csv", models)
    bi_save(ctx, "raw/ses/field_missingness.csv", bi_bind(missing))
    list(status = if (bi_c_complete(estimates) && bi_c_complete(interactions) && bi_c_complete(models)) "PASS" else "PARTIAL", 
        reason = "Adult SES analyses and sample-locked M0/M1/M2 completed; no mediation or childhood-poverty interpretation")
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        bi_ses_derive <- function(d) {
            finite <- is.finite(d$tdi)
            cuts <- if (any(finite)) as.numeric(stats::quantile(d$tdi[finite], c(0.2, 0.4, 0.6, 0.8), type = 7, 
                names = FALSE)) else rep(NA_real_, 4)
            valid <- all(is.finite(cuts)) && !anyDuplicated(cuts)
            d$tdi_quintile <- if (valid) cut(d$tdi, breaks = c(-Inf, cuts, Inf), labels = paste0("Q", 1:5), right = TRUE) else factor(rep(NA_character_, 
                nrow(d)), levels = paste0("Q", 1:5))
            attr(d, "tdi_cuts") <- cuts
            attr(d, "tdi_cut_status") <- if (valid) "OK" else "UNAVAILABLE_NONUNIQUE_QUANTILES"
            d
        }
        bi_ses_extension <- function(ctx, d, exposure, extra, code, label) {
            ff <- bi_formula(exposure, extra = extra)
            unavailable <- bi_c_outcome_unavailable(d, exposure)
            if (!is.null(unavailable)) {
                z <- unavailable[rep(1L, 2L), , drop = FALSE]
                z$model <- c("M0_same_sample", label)
                z$formula <- c(paste(deparse(bi_formula(exposure)), collapse = " "), paste(deparse(ff), collapse = " "))
                z$sample_set <- paste0("complete_case_", label)
                z$extra_covariates <- paste(extra, collapse = ";")
                z$source_N <- nrow(d)
                z$union_complete_N <- z$excluded_union_complete_N <- NA_integer_
                return(z)
            }
            same <- bi_sample(d, ff, exposure)$data
            a <- bi_fit(ctx, same, exposure)
            b <- bi_fit(ctx, same, exposure, extra = extra)
            if (!identical(as.character(a$sample_hash), as.character(b$sample_hash))) stop("SES same-sample hash invariant failed")
            a$model <- "M0_same_sample"
            b$model <- label
            z <- bi_bind(list(a, b))
            z$sample_set <- paste0("complete_case_", label)
            z$extra_covariates <- paste(extra, collapse = ";")
            z$source_N <- nrow(d)
            z$union_complete_N <- nrow(same)
            z$excluded_union_complete_N <- nrow(d) - nrow(same)
            z
        }
        bi_stage_ses <- function(ctx) {
            cfg <- bi_module_config(ctx, "ses")
            adult <- c("bmi", "smoking", "alcohol")
            if (!is.null(cfg$adult_covariates)) adult <- as.character(unlist(cfg$adult_covariates))
            modifiers <- c("tdi", "tdi_quintile", "income", "income_group", "education", "university")
            extensions <- list(M1a = "tdi", M1b = "income", M1c = "education", M1d = c("tdi", "income", "education"), 
                M2 = adult)
            estimates <- interactions <- models <- missing <- list()
            k <- j <- m <- h <- 0L
            for (code in bi_c_codes(ctx)) {
                d <- bi_ses_derive(bi_data(ctx, code))
                for (exposure in bi_c_exposures()) {
                  for (modifier in modifiers) {
                    key <- paste("ses_interaction", exposure, code, modifier, sep = "__")
                    z <- bi_c_job(ctx, key, exposure, function() bi_c_outcome_unavailable(d, exposure) %||% bi_interaction(ctx, 
                      d, exposure, modifier))
                    z$modifier <- modifier
                    z$inference_scope <- "Post-hoc exploratory modifier-specific BH; not formal whole-PheWAS q"
                    j <- j + 1L
                    interactions[[j]] <- bi_c_annotate(z, code, "adult_SES_global_interaction", modifier, key, 
                      paste0("ses_exploratory_", modifier))
                    if (modifier != "tdi") for (level in bi_c_modifier_levels(d, modifier)) {
                      skey <- paste("ses_stratum", exposure, code, modifier, level, sep = "__")
                      z <- bi_c_job(ctx, skey, exposure, function() {
                        unavailable <- bi_c_outcome_unavailable(d, exposure)
                        if (!is.null(unavailable)) return(unavailable)
                        keep <- !is.na(d[[modifier]]) & as.character(d[[modifier]]) == level
                        bi_fit(ctx, d[keep, , drop = FALSE], exposure)
                      })
                      k <- k + 1L
                      estimates[[k]] <- bi_c_annotate(z, code, paste0("adult_SES_", modifier), level, skey)
                    }
                  }
                  for (label in names(extensions)) {
                    key <- paste("ses_extension", exposure, code, label, sep = "__")
                    z <- bi_c_job(ctx, key, exposure, function() bi_ses_extension(ctx, d, exposure, extensions[[label]], 
                      code, label))
                    if (!"model" %in% names(z)) {
                      z <- z[rep(1L, 2L), , drop = FALSE]
                      z$model <- c("M0_same_sample", label)
                    }
                    z <- bi_c_annotate(z, code, "adult_SES_same_sample_extension", label, paste(key, z$model, sep = "__"))
                    m <- m + 1L
                    models[[m]] <- z
                  }
                  for (field in unique(c("tdi", "income", "education", adult))) {
                    h <- h + 1L
                    eligible <- !is.na(d[[exposure]])
                    present <- field %in% names(d)
                    missing[[h]] <- data.frame(trait = code, exposure = exposure, field = field, N = sum(eligible), 
                      missing_N = if (present) sum(eligible & is.na(d[[field]])) else sum(eligible), status = if (present && 
                        any(!is.na(d[[field]]))) "AVAILABLE" else "UNAVAILABLE_FIELD")
                  }
                }
                cuts <- data.frame(trait = code, quantile = c(0.2, 0.4, 0.6, 0.8), value = attr(d, "tdi_cuts"), 
                  status = attr(d, "tdi_cut_status"), reference = "Entire frozen baseline cohort before outcome/exposure complete-case filtering; type7; ties are never split")
                bi_save(ctx, paste0("raw/ses/tdi_quantile_definition_", gsub("[.]", "_", code), ".csv"), cuts)
            }
            estimates <- bi_bind(estimates)
            interactions <- bi_adjust(ctx, bi_bind(interactions))
            models <- bi_bind(models)
            bi_save(ctx, "raw/ses/stratified_results.csv", estimates)
            bi_save(ctx, "raw/ses/interactions.csv", interactions)
            bi_save(ctx, "raw/ses/same_sample_models.csv", models)
            bi_save(ctx, "raw/ses/field_missingness.csv", bi_bind(missing))
            list(status = if (bi_c_complete(estimates) && bi_c_complete(interactions) && bi_c_complete(models)) "PASS" else "PARTIAL", 
                reason = "Adult SES analyses and sample-locked M0/M1/M2 completed; no mediation or childhood-poverty interpretation")
        }
    }), envir = context)
}

