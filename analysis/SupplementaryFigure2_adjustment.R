# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)

num <- function(x) suppressWarnings(as.numeric(as.character(x)))

hash <- function(x) digest(x, algo = "sha256")

filehash <- function(x) digest(file = x, algo = "sha256")

idhash <- function(x) digest(enc2utf8(paste0(paste(sort(as.character(x), method = "radix"), collapse = "\n"), "\n")), 
    algo = "sha256", serialize = FALSE)

csv <- function(x, p) {
    assert(!file.exists(p), paste("REFUSE_EXISTING_OUTPUT", p))
    write.csv(x, p, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}

js <- function(x, p) {
    assert(!file.exists(p), paste("REFUSE_EXISTING_OUTPUT", p))
    write_json(x, p, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null")
}

log_event <- function(event, ...) {
    x <- c(list(time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), event = event), list(...))
    line <- toJSON(x, auto_unbox = TRUE, na = "null", null = "null")
    cat(line, "\n", file = file.path(OUT, "logs/progress.jsonl"), append = TRUE)
    cat(line, "\n")
    flush.console()
}

read_signed <- function(p) {
    z <- readRDS(p)
    assert(is.list(z) && all(c("value", "signature") %in% names(z)), paste("SIGNED_CACHE_REQUIRED", p))
    if ("value_hash" %in% names(z)) 
        assert(identical(z$value_hash, hash(z$value)), "CACHE_VALUE_HASH_MISMATCH")
    else {
        h <- z$hash
        z$hash <- NULL
        assert(identical(h, hash(z)), "CACHE_INTERNAL_HASH_MISMATCH")
    }
    z
}

safe_input <- function(p) {
    assert(file.exists(p), paste("MISSING_INPUT", p))
    p <- normalizePath(p, mustWork = TRUE)
    assert(startsWith(p, normalizePath(INPUT_ROOT)) || identical(p, normalizePath(WITHDRAWALS)), "UNEXPECTED_INPUT_ROOT")
    p
}

needed <- c("eid", "border_relative", "birth_rank", "S", ".exposure_eligible", "age", "sex", "ethnicity", "center", 
    "assessment_year", "date_attend", "edu", "inc", "tdi", "smoke_status", "alcohol_status", "bmi", "cnt_brother", 
    "cnt_sister", "cnt_older", "adopted_num")

codes <- c(DM = "250", T2D = "250.2", IHD = "411", Angina = "411.3", CAS = "411.4")

baseline <- c("border_relative", "age", "sex", "ethnicity", "center", "assessment_year", "S")

edu_tokens <- function(s) {
    if (is.na(s) || !nzchar(trimws(s)) || toupper(trimws(s)) %in% c("NA", "NULL", "[]")) 
        return(numeric())
    s <- chartr("[]", "  ", s)
    s <- gsub("[\"']", "", s, perl = TRUE)
    unique(suppressWarnings(as.numeric(strsplit(trimws(s), "[|,;[:space:]]+", perl = TRUE)[[1]])))
}

education <- function(x) {
    h <- c(`1` = "Degree", `5` = "NVQ_HND_HNC", `6` = "Professional", `2` = "A_AS", `3` = "O_GCSE", `4` = "CSE", 
        `-7` = "None")
    z <- vapply(as.character(x), function(s) {
        v <- edu_tokens(s)
        if (length(v) == 1L && !is.na(v) && v == 7) 
            v <- -7
        if (!length(v) || anyNA(v) || any(!v %in% c(1:6, -7)) || (-7 %in% v && length(v) > 1L)) 
            return(NA_character_)
        unname(h[which(names(h) %in% as.character(v))[1L]])
    }, character(1))
    factor(z, levels = c("None", "CSE", "O_GCSE", "A_AS", "Professional", "NVQ_HND_HNC", "Degree"))
}

status_factor <- function(x, allowed) {
    x <- num(x)
    x[!x %in% allowed] <- NA_real_
    factor(x, levels = allowed)
}

days <- function(x) {
    z <- num(x)
    z[!is.finite(z) | !z %in% 0:7] <- NA_real_
    z
}

extras <- c("education", "household_income", "Townsend", "smoking_status", "alcohol_status", "BMI", "moderate_activity_days", 
    "walking_days")

complete <- function(dd, fields) {
    ok <- complete.cases(dd[, fields, drop = FALSE])
    for (v in fields) if (is.numeric(dd[[v]])) 
        ok <- ok & is.finite(dd[[v]])
    ok
}

fit_one <- function(dd, disease, code, label, model) {
    tm <- proc.time()[["elapsed"]]
    dd <- droplevels(dd)
    f <- forms[[model]]
    w <- character()
    r <- data.frame(disease = disease, phecode = code, sample_set = if (label == "M1_original_sample") 
        "original_sample"
    else "common_sample", model = label, N = nrow(dd), cases = sum(dd$y == 1), controls = sum(dd$y == 0), beta = NA_real_, 
        SE = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_, q_BH15 = NA_real_, estimator = "glm_binomial_logit_ML_model_based_Wald", 
        status = "NOT_ESTIMABLE", warning = "", error = "", converged = NA, iterations = NA_integer_, rank = NA_integer_, 
        columns = NA_integer_, B_variance = var(dd$border_relative), sample_hash = idhash(dd$eid), formula = paste(deparse(f), 
            collapse = " "), reference_levels = "", elapsed_seconds = NA_real_, stringsAsFactors = FALSE)
    tryCatch({
        assert(nrow(dd) > 0 && r$cases >= 20 && r$controls >= 20, "INSUFFICIENT_CASE_CONTROL_SUPPORT")
        assert(is.finite(r$B_variance) && r$B_variance > 0, "B_HAS_NO_VARIATION")
        fac <- names(dd)[vapply(dd, is.factor, logical(1))]
        r$reference_levels <- paste(c(vapply(fac, function(k) paste0(k, "=", levels(dd[[k]])[1]), character(1)), 
            paste0("S=", min(dd$S))), collapse = ";")
        assert(identical(levels(dd$sex)[1], "0") && identical(levels(dd$ethnicity)[1], "White") && identical(levels(dd$center)[1], 
            "c11011") && identical(levels(dd$assessment_year)[1], "2009") && min(dd$S) == 2, "ORIGINAL_M1_REFERENCE_CHANGED")
        x <- model.matrix(f, dd)
        r$columns <- ncol(x)
        r$rank <- qr(x)$rank
        assert(all(is.finite(x)) && r$rank == r$columns, "RANK_DEFICIENT_NO_COLUMN_DELETION")
        rm(x)
        fit <- withCallingHandlers(glm(f, data = dd, family = binomial(link = "logit"), na.action = na.fail, control = glm.control(maxit = 25), 
            model = FALSE, y = FALSE), warning = function(e) {
            w <<- c(w, conditionMessage(e))
            invokeRestart("muffleWarning")
        })
        r$converged <- isTRUE(fit$converged)
        r$iterations <- fit$iter
        cf <- coef(summary(fit))
        assert(r$converged, "GLM_NONCONVERGED")
        assert(all(is.finite(coef(fit))), "NONFINITE_COEFFICIENT")
        assert("border_relative" %in% rownames(cf) && all(is.finite(cf["border_relative", ])) && cf["border_relative", 
            2] > 0, "TARGET_NOT_IDENTIFIABLE")
        assert(nobs(fit) == nrow(dd), "FIT_CHANGED_FIXED_SAMPLE")
        r$beta <- unname(cf["border_relative", 1])
        r$SE <- unname(cf["border_relative", 2])
        r$P <- unname(cf["border_relative", 4])
        r$OR <- exp(r$beta)
        r$CI_low <- exp(r$beta - qnorm(0.975) * r$SE)
        r$CI_high <- exp(r$beta + qnorm(0.975) * r$SE)
        r$status <- if (length(w)) 
            "OK_WARNING_DESCRIPTIVE_PENDING"
        else "OK_DESCRIPTIVE_PENDING"
        if (any(!is.finite(c(r$OR, r$CI_low, r$CI_high)))) 
            r$status <- "WALD_OVERFLOW_DESCRIPTIVE_PENDING"
        coefout <- data.frame(disease = disease, model = label, term = rownames(cf), beta = cf[, 1], SE = cf[, 
            2], P = cf[, 4], row.names = NULL)
        csv(coefout, file.path(OUT, "models", paste0(disease, "_", label, "_coefficients.csv")))
    }, error = function(e) {
        r$error <<- conditionMessage(e)
        r$beta <<- r$SE <<- r$OR <<- r$CI_low <<- r$CI_high <<- r$P <<- NA_real_
    })
    r$warning <- paste(unique(w), collapse = "; ")
    r$elapsed_seconds <- proc.time()[["elapsed"]] - tm
    csv(r, file.path(OUT, "models", paste0(disease, "_", label, ".csv")))
    log_event("MODEL_FINISHED", disease = disease, model = label, N = r$N, cases = r$cases, OR = r$OR, status = r$status, 
        seconds = r$elapsed_seconds)
    gc(FALSE)
    r
}

reference_map <- c(sex = "0", ethnicity = "White", center = "c11011", assessment_year = "2009", education = "None", 
    household_income = "1", smoking_status = "0", alcohol_status = "0")

lev <- list()

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.poly"))
        suppressPackageStartupMessages({
            library(data.table)
            library(jsonlite)
            library(digest)
        })
        args <- commandArgs(TRUE)
        stopifnot(length(args) == 1L)
        OUT <- normalizePath(args[1], mustWork = TRUE)
        INPUT_ROOT <- Sys.getenv("BIRTH_INPUT_ROOT", birth_private_setting("PRIVATE_BIRTH"))
        SOURCE_ROOT <- Sys.getenv("BIRTH_SOURCE_ROOT", file.path(INPUT_ROOT, "border_primary_S8_relative_legacy_locked_v3_run_20261020_004726"))
        SOCIAL <- Sys.getenv("BIRTH_SOCIAL_INPUT", file.path(INPUT_ROOT, "border_socio_behavior_mediation_S8_legacy_phewas_run_20261020_015611_308966/private/social_fields.csv"))
        WITHDRAWALS <- Sys.getenv("BIRTH_WITHDRAWALS", birth_private_setting("PRIVATE_WITHDRAWAL_CSV"))
        HISTORICAL <- Sys.getenv("BIRTH_HISTORICAL_RESULTS", file.path(INPUT_ROOT, "border_sibship_S8_phewas_20261016_v1/S8_phewas_all_results.csv"))
        stopifnot(Sys.info()[["sysname"]] == "Linux", !file.exists(file.path(OUT, "RUN_COMPLETE.json")))
        for (f in c("models", "logs", "private")) dir.create(file.path(OUT, f), showWarnings = FALSE)
        assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
        num <- function(x) suppressWarnings(as.numeric(as.character(x)))
        hash <- function(x) digest(x, algo = "sha256")
        filehash <- function(x) digest(file = x, algo = "sha256")
        idhash <- function(x) digest(enc2utf8(paste0(paste(sort(as.character(x), method = "radix"), collapse = "\n"), 
            "\n")), algo = "sha256", serialize = FALSE)
        csv <- function(x, p) {
            assert(!file.exists(p), paste("REFUSE_EXISTING_OUTPUT", p))
            write.csv(x, p, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
        }
        js <- function(x, p) {
            assert(!file.exists(p), paste("REFUSE_EXISTING_OUTPUT", p))
            write_json(x, p, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null")
        }
        log_event <- function(event, ...) {
            x <- c(list(time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), event = event), list(...))
            line <- toJSON(x, auto_unbox = TRUE, na = "null", null = "null")
            cat(line, "\n", file = file.path(OUT, "logs/progress.jsonl"), append = TRUE)
            cat(line, "\n")
            flush.console()
        }
        read_signed <- function(p) {
            z <- readRDS(p)
            assert(is.list(z) && all(c("value", "signature") %in% names(z)), paste("SIGNED_CACHE_REQUIRED", p))
            if ("value_hash" %in% names(z)) assert(identical(z$value_hash, hash(z$value)), "CACHE_VALUE_HASH_MISMATCH") else {
                h <- z$hash
                z$hash <- NULL
                assert(identical(h, hash(z)), "CACHE_INTERNAL_HASH_MISMATCH")
            }
            z
        }
        safe_input <- function(p) {
            assert(file.exists(p), paste("MISSING_INPUT", p))
            p <- normalizePath(p, mustWork = TRUE)
            assert(startsWith(p, normalizePath(INPUT_ROOT)) || identical(p, normalizePath(WITHDRAWALS)), "UNEXPECTED_INPUT_ROOT")
            p
        }
        started <- Sys.time()
        t0 <- proc.time()[["elapsed"]]
        paths <- c(cohort = file.path(SOURCE_ROOT, "private/checkpoints/cohort.rds"), legacy_states = file.path(SOURCE_ROOT, 
            "private/checkpoints/legacy_states.rds"), social = SOCIAL, withdrawals = WITHDRAWALS, historical = HISTORICAL)
        invisible(lapply(paths, safe_input))
        st0 <- file.info(paths)[, c("size", "mtime")]
        log_event("READ_FROZEN_INPUTS", input_root = INPUT_ROOT, source_root = SOURCE_ROOT, output_root = OUT)
        assert(identical(filehash(HISTORICAL), birth_private_setting("PRIVATE_HASH_3")), "HISTORICAL_FIGURE_RESULT_CHANGED")
        co <- read_signed(paths["cohort"])
        ph <- read_signed(paths["legacy_states"])
        assert(identical(co$signature, ph$signature), "COHORT_STATE_SIGNATURE_DIFFERS")
        source_signature <- co$signature
        d <- co$value$data
        ph <- ph$value
        assert(!anyNA(d$eid) && !anyDuplicated(d$eid) && identical(as.character(d$eid), as.character(ph$ids)), 
            "SOURCE_ID_ALIGNMENT")
        needed <- c("eid", "border_relative", "birth_rank", "S", ".exposure_eligible", "age", "sex", "ethnicity", 
            "center", "assessment_year", "date_attend", "edu", "inc", "tdi", "smoke_status", "alcohol_status", 
            "bmi", "cnt_brother", "cnt_sister", "cnt_older", "adopted_num")
        assert(all(needed %in% names(d)), "MISSING_REQUIRED_COHORT_COLUMNS")
        d <- d[, needed]
        rm(co)
        gc(FALSE)
        codes <- c(DM = "250", T2D = "250.2", IHD = "411", Angina = "411.3", CAS = "411.4")
        assert(all(codes %in% names(ph$states)), "MISSING_FROZEN_OUTCOME")
        assert(all(vapply(ph$states[codes], function(x) all(is.na(x) | x %in% 0:1), logical(1))), "INVALID_OUTCOME_STATE")
        assert(all(vapply(ph$states[codes], function(x) all(is.na(x[!ph$present])), logical(1))), "ABSENT_OLD_ID_BECAME_CONTROL")
        for (code in codes) d[[paste0("Y_", code)]] <- ph$states[[code]]
        d$old_matrix_present <- ph$present
        rm(ph)
        gc(FALSE)
        w <- fread(WITHDRAWALS, header = FALSE, colClasses = "character", data.table = FALSE)[[1]]
        w <- setdiff(w, c("eid", "participant_id"))
        assert(all(grepl("^[0-9]+$", w)), "WITHDRAWAL_SCHEMA")
        new_withdrawals <- sum(d$eid %in% w)
        d <- d[!d$eid %in% w, , drop = FALSE]
        eligible <- d$.exposure_eligible & is.finite(d$border_relative) & is.finite(d$S) & d$S >= 2 & d$S <= 8 & 
            d$border_relative >= 0 & d$border_relative <= 1
        assert(!anyNA(eligible), "EXPOSURE_ELIGIBILITY_UNKNOWN")
        assert(all(abs(d$border_relative[eligible] - (d$birth_rank[eligible] - 1)/(d$S[eligible] - 1)) < 1e-12), 
            "B_IDENTITY")
        assert(all(d$S[eligible] == d$cnt_brother[eligible] + d$cnt_sister[eligible] + 1) && !any(d$adopted_num[eligible] == 
            1, na.rm = TRUE), "BIRTH_ELIGIBILITY_CHANGED")
        d <- d[eligible, , drop = FALSE]
        baseline <- c("border_relative", "age", "sex", "ethnicity", "center", "assessment_year", "S")
        assert(all(vapply(d[c("sex", "ethnicity", "center", "assessment_year")], is.factor, logical(1))), "BASE_COVARIATE_FACTORS_CHANGED")
        assert(identical(levels(d$sex), c("0", "1")) && identical(levels(d$ethnicity), c("White", "Mixed", "Asian or Asian British", 
            "Black or Black British", "Chinese", "Other ethnic group")), "BASE_FACTOR_ENCODING_CHANGED")
        edu_tokens <- function(s) {
            if (is.na(s) || !nzchar(trimws(s)) || toupper(trimws(s)) %in% c("NA", "NULL", "[]")) return(numeric())
            s <- chartr("[]", "  ", s)
            s <- gsub("[\"']", "", s, perl = TRUE)
            unique(suppressWarnings(as.numeric(strsplit(trimws(s), "[|,;[:space:]]+", perl = TRUE)[[1]])))
        }
        education <- function(x) {
            h <- c(`1` = "Degree", `5` = "NVQ_HND_HNC", `6` = "Professional", `2` = "A_AS", `3` = "O_GCSE", `4` = "CSE", 
                `-7` = "None")
            z <- vapply(as.character(x), function(s) {
                v <- edu_tokens(s)
                if (length(v) == 1L && !is.na(v) && v == 7) v <- -7
                if (!length(v) || anyNA(v) || any(!v %in% c(1:6, -7)) || (-7 %in% v && length(v) > 1L)) return(NA_character_)
                unname(h[which(names(h) %in% as.character(v))[1L]])
            }, character(1))
            factor(z, levels = c("None", "CSE", "O_GCSE", "A_AS", "Professional", "NVQ_HND_HNC", "Degree"))
        }
        d$education <- education(d$edu)
        status_factor <- function(x, allowed) {
            x <- num(x)
            x[!x %in% allowed] <- NA_real_
            factor(x, levels = allowed)
        }
        d$household_income <- status_factor(d$inc, 1:5)
        d$Townsend <- num(d$tdi)
        d$Townsend[!is.finite(d$Townsend)] <- NA_real_
        d$smoking_status <- status_factor(d$smoke_status, 0:2)
        d$alcohol_status <- status_factor(d$alcohol_status, 0:2)
        d$BMI <- num(d$bmi)
        assert(all(is.na(d$BMI) | (is.finite(d$BMI) & d$BMI > 0)), "BMI_REQUIRES_INPUT_CODING_REVIEW_NO_AUTOMATIC_TRIMMING")
        social <- fread(SOCIAL, select = c("eid", "p53_i0", "p884_i0", "p864_i0"), colClasses = "character", data.table = FALSE, 
            na.strings = c("", "NA"), showProgress = FALSE)
        assert(!anyNA(social$eid) && !anyDuplicated(social$eid), "SOCIAL_ID_INVALID")
        j <- match(d$eid, social$eid)
        assert(!anyNA(j), "SOCIAL_SOURCE_ROSTER_GAP")
        dates <- suppressWarnings(as.Date(social$p53_i0[j]))
        known_date <- !is.na(d$date_attend)
        assert(all(!is.na(dates[known_date]) & dates[known_date] == d$date_attend[known_date]), "ACTIVITY_BASELINE_INSTANCE_MISMATCH")
        days <- function(x) {
            z <- num(x)
            z[!is.finite(z) | !z %in% 0:7] <- NA_real_
            z
        }
        d$moderate_activity_days <- days(social$p884_i0[j])
        d$walking_days <- days(social$p864_i0[j])
        rm(social)
        gc(FALSE)
        extras <- c("education", "household_income", "Townsend", "smoking_status", "alcohol_status", "BMI", "moderate_activity_days", 
            "walking_days")
        complete <- function(dd, fields) {
            ok <- complete.cases(dd[, fields, drop = FALSE])
            for (v in fields) if (is.numeric(dd[[v]])) ok <- ok & is.finite(dd[[v]])
            ok
        }
        forms <- list(M1 = y ~ border_relative + age + sex + ethnicity + center + assessment_year + factor(S), 
            M2 = y ~ border_relative + age + sex + ethnicity + center + assessment_year + factor(S) + education + 
                household_income + Townsend, M3 = y ~ border_relative + age + sex + ethnicity + center + assessment_year + 
                factor(S) + education + household_income + Townsend + smoking_status + alcohol_status + BMI + moderate_activity_days + 
                walking_days)
        historical <- read.csv(HISTORICAL, check.names = FALSE, colClasses = "character")
        historical <- historical[historical$branch == "border_relative" & historical$phecode %in% codes, , drop = FALSE]
        assert(nrow(historical) == 5 && !anyDuplicated(historical$phecode), "HISTORICAL_FIVE_ROWS")
        csv(historical, file.path(OUT, "HISTORICAL_PHEWAS_REFERENCE.csv"))
        inputmeta <- data.frame(role = names(paths), path = unname(paths), bytes = st0$size, mtime = as.character(st0$mtime), 
            sha256 = vapply(paths, filehash, character(1)))
        csv(inputmeta, file.path(OUT, "INPUT_MANIFEST.csv"))
        js(list(source = "FROZEN_OLD_PHEWAS_AS_CONFIRMED_BY_USER_NOT_E1_OR_REPAIRED_V2", withdrawals_in_already_filtered_roster = new_withdrawals, 
            eligible_roster = nrow(d), old_matrix_absent = sum(!d$old_matrix_present), B_variance = var(d$border_relative), 
            code_signature = source_signature), file.path(OUT, "COHORT_AUDIT.json"))
        results <- flows <- states <- levels_out <- list()
        fit_one <- function(dd, disease, code, label, model) {
            tm <- proc.time()[["elapsed"]]
            dd <- droplevels(dd)
            f <- forms[[model]]
            w <- character()
            r <- data.frame(disease = disease, phecode = code, sample_set = if (label == "M1_original_sample") "original_sample" else "common_sample", 
                model = label, N = nrow(dd), cases = sum(dd$y == 1), controls = sum(dd$y == 0), beta = NA_real_, 
                SE = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_, q_BH15 = NA_real_, 
                estimator = "glm_binomial_logit_ML_model_based_Wald", status = "NOT_ESTIMABLE", warning = "", error = "", 
                converged = NA, iterations = NA_integer_, rank = NA_integer_, columns = NA_integer_, B_variance = var(dd$border_relative), 
                sample_hash = idhash(dd$eid), formula = paste(deparse(f), collapse = " "), reference_levels = "", 
                elapsed_seconds = NA_real_, stringsAsFactors = FALSE)
            tryCatch({
                assert(nrow(dd) > 0 && r$cases >= 20 && r$controls >= 20, "INSUFFICIENT_CASE_CONTROL_SUPPORT")
                assert(is.finite(r$B_variance) && r$B_variance > 0, "B_HAS_NO_VARIATION")
                fac <- names(dd)[vapply(dd, is.factor, logical(1))]
                r$reference_levels <- paste(c(vapply(fac, function(k) paste0(k, "=", levels(dd[[k]])[1]), character(1)), 
                  paste0("S=", min(dd$S))), collapse = ";")
                assert(identical(levels(dd$sex)[1], "0") && identical(levels(dd$ethnicity)[1], "White") && identical(levels(dd$center)[1], 
                  "c11011") && identical(levels(dd$assessment_year)[1], "2009") && min(dd$S) == 2, "ORIGINAL_M1_REFERENCE_CHANGED")
                x <- model.matrix(f, dd)
                r$columns <- ncol(x)
                r$rank <- qr(x)$rank
                assert(all(is.finite(x)) && r$rank == r$columns, "RANK_DEFICIENT_NO_COLUMN_DELETION")
                rm(x)
                fit <- withCallingHandlers(glm(f, data = dd, family = binomial(link = "logit"), na.action = na.fail, 
                  control = glm.control(maxit = 25), model = FALSE, y = FALSE), warning = function(e) {
                  w <<- c(w, conditionMessage(e))
                  invokeRestart("muffleWarning")
                })
                r$converged <- isTRUE(fit$converged)
                r$iterations <- fit$iter
                cf <- coef(summary(fit))
                assert(r$converged, "GLM_NONCONVERGED")
                assert(all(is.finite(coef(fit))), "NONFINITE_COEFFICIENT")
                assert("border_relative" %in% rownames(cf) && all(is.finite(cf["border_relative", ])) && cf["border_relative", 
                  2] > 0, "TARGET_NOT_IDENTIFIABLE")
                assert(nobs(fit) == nrow(dd), "FIT_CHANGED_FIXED_SAMPLE")
                r$beta <- unname(cf["border_relative", 1])
                r$SE <- unname(cf["border_relative", 2])
                r$P <- unname(cf["border_relative", 4])
                r$OR <- exp(r$beta)
                r$CI_low <- exp(r$beta - qnorm(0.975) * r$SE)
                r$CI_high <- exp(r$beta + qnorm(0.975) * r$SE)
                r$status <- if (length(w)) "OK_WARNING_DESCRIPTIVE_PENDING" else "OK_DESCRIPTIVE_PENDING"
                if (any(!is.finite(c(r$OR, r$CI_low, r$CI_high)))) r$status <- "WALD_OVERFLOW_DESCRIPTIVE_PENDING"
                coefout <- data.frame(disease = disease, model = label, term = rownames(cf), beta = cf[, 1], SE = cf[, 
                  2], P = cf[, 4], row.names = NULL)
                csv(coefout, file.path(OUT, "models", paste0(disease, "_", label, "_coefficients.csv")))
            }, error = function(e) {
                r$error <<- conditionMessage(e)
                r$beta <<- r$SE <<- r$OR <<- r$CI_low <<- r$CI_high <<- r$P <<- NA_real_
            })
            r$warning <- paste(unique(w), collapse = "; ")
            r$elapsed_seconds <- proc.time()[["elapsed"]] - tm
            csv(r, file.path(OUT, "models", paste0(disease, "_", label, ".csv")))
            log_event("MODEL_FINISHED", disease = disease, model = label, N = r$N, cases = r$cases, OR = r$OR, 
                status = r$status, seconds = r$elapsed_seconds)
            gc(FALSE)
            r
        }
        for (disease in names(codes)) {
            code <- codes[[disease]]
            dd <- d
            dd$y <- dd[[paste0("Y_", code)]]
            original <- complete(dd, c("y", baseline))
            common <- original & complete(dd, extras)
            states[[disease]] <- data.frame(disease = disease, phecode = code, eligible_roster = nrow(dd), old_matrix_present = sum(dd$old_matrix_present), 
                outcome_known = sum(!is.na(dd$y)), outcome_cases = sum(dd$y == 1, na.rm = TRUE), outcome_controls = sum(dd$y == 
                  0, na.rm = TRUE), unknown_outcome = sum(is.na(dd$y)), unknown_due_absent_old_matrix = sum(!dd$old_matrix_present), 
                known_but_M1_covariates_incomplete = sum(!is.na(dd$y) & !original))
            flow <- data.frame(disease = disease, phecode = code, original_N = sum(original), original_cases = sum(dd$y[original] == 
                1), original_controls = sum(dd$y[original] == 0), common_N = sum(common), common_cases = sum(dd$y[common] == 
                1), common_controls = sum(dd$y[common] == 0), lost_N = sum(original) - sum(common), lost_percent = 100 * 
                (sum(original) - sum(common))/sum(original), original_sample_hash = idhash(dd$eid[original]), common_sample_hash = idhash(dd$eid[common]))
            for (v in extras) {
                miss <- !complete(dd, v)
                flow[[paste0("missing_", v)]] <- sum(original & miss)
                flow[[paste0("missing_", v, "_percent")]] <- 100 * sum(original & miss)/sum(original)
            }
            hist <- historical[historical$phecode == code, , drop = FALSE]
            flow$historical_figure_N <- num(hist$N)
            flow$current_minus_historical_N <- sum(original) - num(hist$N)
            flow$historical_sample_hash_equal <- identical(flow$original_sample_hash, hist$sample_hash)
            flows[[disease]] <- flow
            for (label in c("M1_original_sample", "M1_common", "M2_common", "M3_common")) {
                model <- sub("_.*", "", label)
                sel <- if (label == "M1_original_sample") original else common
                ans <- fit_one(dd[sel, c("eid", "y", baseline, extras), drop = FALSE], disease, code, label, model)
                results[[paste(disease, label)]] <- ans
            }
        }
        res <- as.data.frame(rbindlist(results))
        flow <- as.data.frame(rbindlist(flows))
        state <- as.data.frame(rbindlist(states))
        assert(nrow(res) == 20 && !anyDuplicated(res[c("disease", "model")]), "TWENTY_RESULT_SLOTS_REQUIRED")
        for (g in names(codes)) {
            z <- res[res$disease == g & res$sample_set == "common_sample", ]
            assert(nrow(z) == 3 && length(unique(z$sample_hash)) == 1 && length(unique(z$N)) == 1 && length(unique(z$cases)) == 
                1 && length(unique(z$controls)) == 1, "COMMON_SAMPLE_NOT_IDENTICAL")
        }
        family <- which(res$sample_set == "common_sample")
        assert(length(family) == 15, "BH_FAMILY_CHANGED")
        ok <- family[is.finite(res$P[family]) & res$P[family] >= 0 & res$P[family] <= 1]
        if (length(ok)) res$q_BH15[ok] <- p.adjust(res$P[ok], method = "BH", n = 15)
        csv(res, file.path(OUT, "MODEL_RESULTS_LONG.csv"))
        csv(flow, file.path(OUT, "SAMPLE_FLOW.csv"))
        csv(state, file.path(OUT, "OUTCOME_STATE_AUDIT.csv"))
        reference_map <- c(sex = "0", ethnicity = "White", center = "c11011", assessment_year = "2009", education = "None", 
            household_income = "1", smoking_status = "0", alcohol_status = "0")
        lev <- list()
        for (v in c("sex", "ethnicity", "center", "assessment_year", extras)) {
            x <- d[[v]]
            if (is.factor(x)) for (l in levels(x)) lev[[paste(v, l)]] <- data.frame(variable = v, level = l, model_reference = l == 
                reference_map[[v]], N_in_eligible_roster = sum(x == l, na.rm = TRUE), missing_in_eligible_roster = sum(is.na(x)))
        }
        csv(as.data.frame(rbindlist(lev)), file.path(OUT, "COVARIATE_ENCODING_COUNTS.csv"))
        st1 <- file.info(paths)[, c("size", "mtime")]
        assert(identical(st0, st1), "SOURCE_CHANGED_DURING_ANALYSIS")
        js(list(started = as.character(started), finished = as.character(Sys.time()), elapsed_seconds = proc.time()[["elapsed"]] - 
            t0, source_authority = "FROZEN_OLD_PHEWAS_CONFIRMED_BY_USER", models_processed = nrow(res), finite_B_estimates = sum(is.finite(res$beta) & 
            is.finite(res$SE) & is.finite(res$P)), status_counts = as.list(table(res$status)), fixed_BH_family = 15, 
            common_sample_identity = "PASS_5_OF_5", historical_results_overwritten = FALSE, participant_data_exported = FALSE, 
            R_version = R.version.string, packages = lapply(c("data.table", "jsonlite", "digest"), function(p) list(package = p, 
                version = as.character(packageVersion(p))))), file.path(OUT, "RUN_COMPLETE.json"))
        writeLines(capture.output(sessionInfo()), file.path(OUT, "SESSION_INFO.txt"))
        log_event("RUN_COMPLETE", models = 20, finite = sum(is.finite(res$P)), elapsed_seconds = proc.time()[["elapsed"]] - 
            t0)
    }), envir = context)
}

