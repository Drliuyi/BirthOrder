# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

s8_branches <- c("border1n", "border_rank", "border_relative")

s8_M <- c(border1n = 1362L, border_rank = 1419L, border_relative = 1419L)

s8_make_sample <- function(full, row = NULL) {
    d <- full$data
    keep <- !is.na(d$S) & d$S <= 8
    before <- d
    d <- droplevels(d[keep, , drop = FALSE])
    removed <- character()
    for (k in names(d)) if (is.factor(before[[k]])) {
        lost <- setdiff(levels(before[[k]]), levels(d[[k]]))
        if (length(lost)) 
            removed <- c(removed, paste0(k, ":", lost))
        if (nrow(d) && length(levels(before[[k]])) > 1L && length(levels(d[[k]])) >= 1L) 
            review_assert(identical(levels(before[[k]])[1], levels(d[[k]])[1]), paste("REFERENCE_CHANGED", k))
    }
    f <- full$formula
    terms <- attr(terms(f), "term.labels")
    target <- terms[1L]
    constants <- character()
    if (nrow(d)) 
        for (term in setdiff(terms, target)) {
            column <- if (term == "factor(S)") 
                "S"
            else term
            if (length(unique(d[[column]])) <= 1L) 
                constants <- c(constants, term)
        }
    if (nrow(d) && length(unique(d$S)) > 1 && "factor(S)" %in% terms) 
        review_assert(min(d$S) == min(before$S), "FAMILY_SIZE_REFERENCE_CHANGED")
    list(data = d, formula = f, removed_empty_levels = paste(removed, collapse = ";"), removed_terms = full$removed_terms, 
        new_constants = paste(constants, collapse = ";"), counts = list(N = nrow(d), cases = sum(d$y == 1), noncases = sum(d$y == 
            0)))
}

s8_saved_plan <- function(row) {
    value <- if (row$branch == "border_relative") 
        row$planned_eligible
    else row$planned
    review_assert(value %in% c("TRUE", "FALSE"), "SAVED_PLAN_FLAG_MISSING")
    value == "TRUE"
}

s8_fit_one <- function(row, dat, ph) {
    started <- proc.time()[["elapsed"]]
    b <- row$branch
    p <- row$phecode
    r <- as.list(row[c("branch", "phecode", "description", "category", "groupnum", "color", "sex_rule", "metadata_verified")])
    r$key <- review_key(b, p)
    r$registered_M <- unname(s8_M[b])
    r$planned_eligible <- s8_saved_plan(row)
    for (k in c("N", "cases", "noncases")) r[[paste0("full_", k)]] <- as.integer(row[[k]])
    for (k in c("beta", "se", "p")) r[[paste0("full_", k)]] <- as.numeric(row[[paste0("saved_raw_", k)]])
    for (k in c("OR", "CI_low", "CI_high")) r[[paste0("full_", k)]] <- as.numeric(row[[paste0("saved_", k)]])
    r$full_q <- as.numeric(row$recovered_wald_BH_q)
    r$full_lp_state <- row$separation_assessment
    r$N <- r$cases <- r$noncases <- NA_integer_
    r$beta <- r$se <- r$OR <- r$CI_low <- r$CI_high <- r$p <- NA_real_
    r$raw_target_beta <- r$raw_target_se <- r$raw_target_p <- NA_real_
    r$converged <- NA
    r$matrix_rank <- r$matrix_columns <- NA_integer_
    r$warning <- ""
    r$warning_count <- 0L
    r$status <- "NOT_ESTIMABLE"
    r$error <- ""
    r$family_match <- FALSE
    r$full_sample_hash_match <- FALSE
    r$lp_status <- "NOT_RUN_FULL_LP_NOT_AUTHORIZED"
    r$separation_detected <- NA
    r$target_infinite <- NA
    r$inference_status <- "DESCRIPTIVE_WALD_PENDING_REVIEW"
    r$formal_p <- r$formal_q <- NA_real_
    answer <- tryCatch({
        full <- linear_sample(dat, ph[[p]], b, row$sex_rule)
        for (k in c("N", "cases", "noncases")) review_assert(full$counts[[k]] == as.integer(row[[k]]), paste("FULL_SAMPLE_COUNT_MISMATCH", 
            k))
        h <- review_hash_ids(full$data$eid)
        savedhash <- if (b == "border_relative") 
            row$sample_hash_sha256
        else row$sample_id_digest
        if (!is.na(savedhash) && nzchar(savedhash)) {
            review_assert(identical(if (b == "border_relative") 
                h$win
            else h$mac, savedhash), "FULL_SAMPLE_HASH_MISMATCH")
            r$full_sample_hash_match <- TRUE
        }
        review_assert(identical(review_normalize_formula(review_formula_text(full$formula)), review_normalize_formula(row$actual_formula)), 
            "FULL_FORMULA_MISMATCH")
        planned <- full$counts$cases >= 20L && full$counts$noncases >= 20L && (b == "border_relative" || length(unique(full$data[[b]])) >= 
            2L)
        review_assert(identical(planned, r$planned_eligible), "REGISTERED_TEST_FAMILY_MISMATCH")
        r$family_match <- TRUE
        s <- s8_make_sample(full, row)
        d <- s$data
        for (k in c("N", "cases", "noncases")) r[[k]] <- s$counts[[k]]
        r$sample_hash <- review_hash_ids(d$eid)$win
        r$actual_formula <- review_formula_text(s$formula)
        r$removed_terms <- paste(s$removed_terms, collapse = ";")
        r$removed_empty_levels <- s$removed_empty_levels
        r$effect_unit <- unname(linear_units[b])
        r$S_factor_levels <- paste(sort(unique(d$S)), collapse = ";")
        fac <- names(d)[vapply(d, is.factor, logical(1))]
        r$reference_levels <- paste(vapply(fac, function(k) paste0(k, "=", levels(d[[k]])[1]), character(1)), collapse = ";")
        r$reference_levels <- paste0(r$reference_levels, ";S=", if (nrow(d)) 
            min(d$S)
        else NA)
        branch_ok <- attr(dat, "s8_branch_ok")[[b]]
        review_assert(!is.null(branch_ok), "MISSING_BRANCH_GATE")
        threshold <- branch_ok && r$cases >= 20L && r$noncases >= 20L && length(unique(d[[b]])) >= 2L
        r$s8_threshold_eligible <- threshold
        review_assert(!threshold || r$planned_eligible, "S8_ELIGIBLE_OUTSIDE_ORIGINAL_FAMILY")
        if (row$metadata_verified != "TRUE") {
            r$status <- "PHECODE_METADATA_UNVERIFIED"
            return(r)
        }
        if (!threshold) {
            r$status <- "LOW_CASE_OR_CONTROL_COUNT_OR_CONSTANT_TARGET"
            return(r)
        }
        if (nzchar(s$new_constants)) {
            r$status <- "NEW_CONSTANT_COVARIATE_NOT_REMOVED"
            r$error <- s$new_constants
            return(r)
        }
        mm <- model.matrix(s$formula, d)
        r$matrix_rank <- qr(mm)$rank
        r$matrix_columns <- ncol(mm)
        if (any(!is.finite(mm)) || r$matrix_rank != r$matrix_columns) {
            r$status <- "RANK_DEFICIENT"
            r$error <- "No informative column deletion"
            return(r)
        }
        review_assert(b %in% colnames(mm), "TARGET_COLUMN_MISSING")
        r$design_hash <- digest::digest(list(x = mm, y = d$y), algo = "sha256")
        rm(mm)
        got <- review_capture(glm(s$formula, data = d, family = binomial(link = "logit")))
        r$warning <- paste(got$warnings, collapse = "; ")
        r$warning_count <- length(got$warnings)
        if (inherits(got$value, "error")) {
            r$status <- "FIT_ERROR"
            r$error <- conditionMessage(got$value)
            return(r)
        }
        fit <- got$value
        r$converged <- fit$converged
        r$iterations <- fit$iter
        cf <- summary(fit)$coefficients
        if (b %in% rownames(cf)) {
            r$raw_target_beta <- unname(cf[b, "Estimate"])
            r$raw_target_se <- unname(cf[b, "Std. Error"])
            r$raw_target_p <- unname(cf[b, "Pr(>|z|)"])
        }
        if (!isTRUE(fit$converged)) {
            r$status <- "NONCONVERGED"
            return(r)
        }
        if (!b %in% rownames(cf) || !all(is.finite(cf[b, ])) || cf[b, "Std. Error"] <= 0) {
            r$status <- "TARGET_NOT_IDENTIFIABLE"
            return(r)
        }
        if (!all(is.finite(coef(fit)))) {
            r$status <- "NONFINITE_COEFFICIENT"
            return(r)
        }
        review_assert(nobs(fit) == r$N, "FIT_SAMPLE_CHANGED")
        r$beta <- r$raw_target_beta
        r$se <- r$raw_target_se
        r$p <- r$raw_target_p
        ci <- exp(c(r$beta, r$beta - qnorm(0.975) * r$se, r$beta + qnorm(0.975) * r$se))
        r$OR <- ci[1]
        r$CI_low <- ci[2]
        r$CI_high <- ci[3]
        r$delta_logOR <- r$beta - r$full_beta
        r$status <- if (length(got$warnings)) 
            "OK_WARNING_DESCRIPTIVE_PENDING"
        else "OK_DESCRIPTIVE_PENDING"
        if (any(!is.finite(ci))) 
            r$status <- "WALD_OVERFLOW_DESCRIPTIVE_PENDING"
        r$p_method <- "two-sided Wald z"
        r$ci_method <- "95% normal Wald; model-based SE"
        r
    }, error = function(e) {
        r$status <<- "FAILED"
        r$error <<- conditionMessage(e)
        r
    })
    answer$elapsed_seconds <- proc.time()[["elapsed"]] - started
    answer$peak_rss_mb <- review_peak_mb()
    answer$finished <- review_stamp()
    answer
}

s8_main <- function() {
    args <- commandArgs(TRUE)
    if (length(args) != 3L) 
        stop("Usage: 10.sibship_sensitivity.R MODULE MODE OUT")
    module <- normalizePath(args[1], mustWork = TRUE)
    mode <- args[2]
    out <- args[3]
    if (!mode %in% c("sibship_check", "sibship_focus", "sibship_report")) 
        stop("Unknown stage")
    extra <- Sys.getenv("REVIEW_R_LIB", birth_private_setting("PRIVATE_LIB"))
    if (dir.exists(extra)) 
        .libPaths(unique(c(.libPaths(), extra)))
    for (p in c("digest", "jsonlite", "dplyr", "ggplot2", "ggrepel", "patchwork", "scales")) if (!requireNamespace(p, 
        quietly = TRUE)) 
        stop("Missing dependency ", p)
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        RhpcBLASctl::blas_set_num_threads(1L)
        RhpcBLASctl::omp_set_num_threads(1L)
    }
    source(file.path(module, "frozen", "common.R"), local = .GlobalEnv)
    source(file.path(module, "frozen", "helper.R"), local = .GlobalEnv)
    source(file.path(module, "f", "11.sibship_report.R"), local = .GlobalEnv)
    cfg <- jsonlite::read_json(file.path(module, "config.json"), simplifyVector = TRUE)
    cfg$out <- out
    cfg$module <- module
    cfg$huang_plotting <- file.path(module, "frozen", "huang_plotting.R")
    for (k in c("all_rds", "phewas_rds", "summary_csv")) cfg[[k]] <- Sys.getenv(paste0("S8_", toupper(k)), cfg[[k]])
    for (d in c("raw", "fig", "reports", "audit", "logs", "checkpoints", "branches", "code_snapshot")) dir.create(file.path(out, 
        d), recursive = TRUE, showWarnings = FALSE)
    review_event(out, "START_S8_FULL_PHEWAS", mode = mode, workers = cfg$workers, total_rows = 4746L)
    prov <- do.call(rbind, lapply(names(cfg$expected_sha), function(k) {
        review_assert(file.exists(cfg[[k]]), paste("Missing input", cfg[[k]]))
        review_event(out, "VERIFY_INPUT_SHA", role = k, path = cfg[[k]])
        hash <- review_sha(cfg[[k]])
        review_assert(hash == cfg$expected_sha[[k]], paste("INPUT_SHA_MISMATCH", k))
        data.frame(role = k, path = cfg[[k]], sha256 = hash)
    }))
    files <- sort(list.files(module, recursive = TRUE, full.names = TRUE))
    files <- files[!file.info(files)$isdir]
    cm <- data.frame(path = substring(files, nchar(module) + 2L), sha256 = vapply(files, review_sha, character(1)))
    sig <- digest::digest(list(config = cfg, inputs = prov, code = cm, R = R.version$version.string, package_versions = vapply(c("stats", 
        "digest", "dplyr"), function(p) as.character(packageVersion(p)), character(1))), algo = "sha256")
    owner <- file.path(out, "audit", "ownership.rds")
    if (file.exists(owner)) 
        review_assert(identical(readRDS(owner)$signature, sig), "OUTPUT_SIGNATURE_MISMATCH")
    else {
        review_rds(list(signature = sig, config = cfg, created = review_stamp()), owner)
        review_csv(prov, file.path(out, "audit", "input_provenance.csv"))
        review_csv(cm, file.path(out, "audit", "code_manifest.csv"))
        capture.output(sessionInfo(), file = file.path(out, "audit", "sessionInfo.txt"))
        for (i in seq_along(files)) {
            dst <- file.path(out, "code_snapshot", cm$path[i])
            dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
            review_assert(file.copy(files[i], dst, overwrite = FALSE), "SNAPSHOT_COPY_FAILED")
        }
    }
    if (file.exists(file.path(out, "complete.json")) && file.exists(file.path(out, "audit", "sha256_manifest.txt"))) {
        dd <- jsonlite::read_json(file.path(out, "complete.json"), simplifyVector = TRUE)
        review_assert(identical(dd$signature, sig) && dd$rows == 4746L, "COMPLETE_SIGNATURE_MISMATCH")
        ml <- readLines(file.path(out, "audit", "sha256_manifest.txt"))
        mp <- substring(ml, 66L)
        mh <- substring(ml, 1L, 64L)
        review_assert(all(file.exists(file.path(out, mp))) && all(vapply(file.path(out, mp), review_sha, character(1)) == 
            mh), "COMPLETED_OUTPUT_HASH_MISMATCH")
        review_event(out, "VERIFIED_ALREADY_COMPLETE", rows = 4746L)
        return(invisible(NULL))
    }
    if (mode != "sibship_report") {
        td <- file.path(out, "audit", paste0("synthetic_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
        log <- file.path(out, "logs", paste0(basename(td), ".log"))
        rc <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(file.path(module, "tests", "test_s8.R")), 
            shQuote(module), shQuote(td)), stdout = log, stderr = log)
        review_assert(identical(rc, 0L), "S8_SYNTHETIC_TEST_FAILED")
    }
    if (mode == "sibship_check") {
        review_event(out, "CHECK_PASSED")
        return(invisible(NULL))
    }
    full <- read.csv(cfg$summary_csv, colClasses = "character", check.names = FALSE, fileEncoding = "UTF-8")
    registry <- full[full$branch %in% s8_branches, , drop = FALSE]
    for (b in s8_branches) {
        a <- registry[registry$branch == b, ]
        review_assert(nrow(a) == 1582L && !anyDuplicated(a$phecode), "SOURCE_REGISTRY_MISMATCH")
        flags <- vapply(seq_len(nrow(a)), function(i) s8_saved_plan(a[i, ]), logical(1))
        review_assert(sum(flags) == s8_M[b], "FIXED_DENOMINATOR_MISMATCH")
    }
    keys <- review_key(registry$branch, registry$phecode)
    review_assert(!anyDuplicated(keys), "DUPLICATE_KEY")
    records <- list()
    cp <- function(k) file.path(out, "checkpoints", paste0(k, ".rds"))
    for (k in keys) if (file.exists(cp(k))) {
        r <- review_cache_read(cp(k), sig)
        review_assert(identical(r$key, k) && identical(review_key(r$branch, r$phecode), k), "CHECKPOINT_IDENTITY_MISMATCH")
        records[[k]] <- r
    }
    if (mode != "sibship_report") {
        review_event(out, "READ_BASELINE")
        base <- readRDS(cfg$all_rds)
        need <- c("eid", "age", "sex", "ethnicity", "center", "date_attend", "cnt_older", "cnt_brother", "cnt_sister", 
            "adopted")
        review_assert(all(need %in% names(base)), "MISSING_FIELDS")
        dat <- linear_covars(make_birth_order(as.data.frame(base[need])))
        rm(base)
        gc()
        dat <- dat[, c("eid", linear_exposures, linear_covariates, "S")]
        review_event(out, "READ_ALL_1582_PHENOTYPES")
        ph <- readRDS(cfg$phewas_rds)
        review_assert(!anyNA(ph$id) && !anyDuplicated(ph$id), "PHENOTYPE_ID_ERROR")
        codes <- setdiff(names(ph), "id")
        review_assert(length(codes) == 1582L && !anyDuplicated(codes), "PHENOTYPE_SET_ERROR")
        for (b in s8_branches) review_assert(setequal(codes, registry$phecode[registry$branch == b]), "PHENOTYPE_FAMILY_SET_MISMATCH")
        ix <- match(as.character(dat$eid), as.character(ph$id))
        keep <- !is.na(ix)
        dat <- dat[keep, , drop = FALSE]
        ph <- ph[ix[keep], codes, drop = FALSE]
        gc()
        review_assert(all(vapply(ph, function(y) all(is.na(y) | y %in% c(0, 1)), logical(1))), "OUTCOME_STATE_INVALID")
        distributions <- flows <- list()
        for (b in linear_exposures) {
            fs <- linear_sample(dat, rep(0, nrow(dat)), b, "both")
            s <- s8_make_sample(fs)
            t <- as.data.frame(table(fs$data$S), stringsAsFactors = FALSE)
            names(t) <- c("S", "N")
            t$S <- as.integer(t$S)
            t$branch <- b
            distributions[[b]] <- t
            flows[[b]] <- data.frame(branch = b, originalN = nrow(fs$data), retainedN = nrow(s$data), removedN = nrow(fs$data) - 
                nrow(s$data), retention_pct = 100 * nrow(s$data)/nrow(fs$data), border12_ref_only = b == "border12")
        }
        distribution <- do.call(rbind, distributions)
        flow <- do.call(rbind, flows)
        attr(dat, "s8_branch_ok") <- setNames(as.list(flow$retainedN >= 500L), flow$branch)
        if (!file.exists(file.path(out, "raw", "sibship_S_distribution.csv"))) 
            review_csv(distribution, file.path(out, "raw", "sibship_S_distribution.csv"))
        if (!file.exists(file.path(out, "raw", "sibship_S8_sample_flow.csv"))) 
            review_csv(flow, file.path(out, "raw", "sibship_S8_sample_flow.csv"))
        if (!file.exists(file.path(out, "audit", "source_registry.csv"))) 
            review_csv(registry, file.path(out, "audit", "source_registry.csv"))
        start <- as.numeric(Sys.time())
        prior_runtime <- sum(vapply(records, function(r) if (is.null(r$elapsed_seconds)) 0 else r$elapsed_seconds, 
            numeric(1)), na.rm = TRUE)/cfg$workers
        jobs <- list()
        pending <- keys[!keys %in% names(records)]
        nexti <- 1L
        fatal <- ""
        failed_record <- function(key, status, error, elapsed) {
            row <- registry[match(key, keys), , drop = FALSE]
            r <- as.list(row)
            r$key <- key
            r$full_N <- as.integer(row$N)
            r$full_cases <- as.integer(row$cases)
            r$full_noncases <- as.integer(row$noncases)
            r$N <- r$cases <- r$noncases <- NA_integer_
            for (k in c("beta", "se", "OR", "CI_low", "CI_high", "p", "raw_target_beta", "raw_target_se", "raw_target_p")) r[[k]] <- NA_real_
            r$status <- status
            r$error <- error
            r$elapsed_seconds <- elapsed
            r$family_match <- NA
            r$planned_eligible <- s8_saved_plan(row)
            r$registered_M <- unname(s8_M[row$branch])
            r
        }
        publish_branch <- function(b) {
            bk <- keys[registry$branch == b]
            dest <- file.path(out, "branches", b)
            if (!all(bk %in% names(records)) || file.exists(file.path(dest, "complete.json"))) 
                return(invisible(NULL))
            zz <- review_bind(records[bk])
            zz$p <- as.numeric(zz$p)
            review_assert(!any(zz$family_match == "FALSE", na.rm = TRUE), "BRANCH_FAMILY_MISMATCH")
            valid <- is.finite(zz$p) & zz$p >= 0 & zz$p <= 1
            review_assert(all(zz$planned_eligible[valid] == "TRUE"), "P_OUTSIDE_FAMILY")
            zz$descriptive_BH_q <- zz$descriptive_bonf <- NA_real_
            zz$descriptive_BH_q[valid] <- p.adjust(zz$p[valid], "BH", n = s8_M[b])
            zz$descriptive_bonf[valid] <- p.adjust(zz$p[valid], "bonferroni", n = s8_M[b])
            dir.create(dest, recursive = TRUE, showWarnings = FALSE)
            if (!file.exists(file.path(dest, "results.csv"))) 
                review_csv(zz, file.path(dest, "results.csv"))
            if (!file.exists(file.path(dest, "results.rds"))) 
                review_rds(zz, file.path(dest, "results.rds"))
            jsonlite::write_json(list(status = "COMPLETED_DESCRIPTIVE_PENDING", branch = b, rows = 1582L, M = unname(s8_M[b]), 
                finite_P = sum(valid), signature = sig, finished = review_stamp()), file.path(dest, "complete.json"), 
                pretty = TRUE, auto_unbox = TRUE)
            review_event(out, "BRANCH_COMPLETE", branch = b, rows = 1582L, finite_P = sum(valid), M = unname(s8_M[b]))
        }
        while (nexti <= length(pending) || length(jobs)) {
            while (nexti <= length(pending) && length(jobs) < cfg$workers && !nzchar(fatal)) {
                k <- pending[nexti]
                row <- registry[match(k, keys), , drop = FALSE]
                job <- parallel::mcparallel(s8_fit_one(row, dat, ph), mc.set.seed = FALSE, silent = TRUE)
                jobs[[as.character(job$pid)]] <- list(job = job, key = k, start = as.numeric(Sys.time()))
                review_event(out, "MODEL_START", key = k, pid_worker = job$pid, completed = length(records), total = 4746L)
                nexti <- nexti + 1L
            }
            for (pid in names(jobs)) {
                j <- jobs[[pid]]
                v <- parallel::mccollect(j$job, wait = FALSE)
                elapsed <- as.numeric(Sys.time()) - j$start
                if (is.null(v) && (elapsed > cfg$per_model_seconds || as.numeric(Sys.time()) - start + prior_runtime > 
                  cfg$total_seconds)) {
                  tools::pskill(j$job$pid, 15L)
                  Sys.sleep(1)
                  v <- parallel::mccollect(j$job, wait = FALSE)
                  if (is.null(v)) {
                    tools::pskill(j$job$pid, 9L)
                    parallel::mccollect(j$job, wait = TRUE)
                  }
                  v <- list(failed_record(j$key, "TIMEOUT", "Owned S8 worker exceeded resource limit", elapsed))
                }
                if (!is.null(v)) {
                  r <- v[[1L]]
                  if (!is.list(r) || inherits(r, "try-error")) {
                    r <- failed_record(j$key, "WORKER_FAILED", "Owned child returned invalid result", elapsed)
                    fatal <- "WORKER_FAILURE"
                  }
                  review_assert(identical(r$key, j$key) && identical(review_key(r$branch, r$phecode), j$key), "WORKER_IDENTITY_MISMATCH")
                  if (is.null(r$elapsed_seconds)) 
                    r$elapsed_seconds <- elapsed
                  review_save_checkpoint(r, cp(j$key), sig)
                  records[[j$key]] <- r
                  jobs[[pid]] <- NULL
                  review_event(out, "MODEL_END", key = j$key, status = r$status, error = r$error, completed = length(records), 
                    total = 4746L, seconds = r$elapsed_seconds)
                  if (grepl("MISMATCH|INVALID|cannot allocate", r$error)) 
                    fatal <- paste("SYSTEMIC_CHECK_FAILED", j$key, r$error)
                  if (!nzchar(fatal)) 
                    for (b in s8_branches) publish_branch(b)
                }
            }
            if (as.numeric(Sys.time()) - start + prior_runtime > cfg$total_seconds) 
                fatal <- "TOTAL_BUDGET_EXHAUSTED"
            if (nzchar(fatal) && !length(jobs)) 
                break
            Sys.sleep(0.2)
        }
        if (nzchar(fatal)) {
            review_event(out, "STOPPED", reason = fatal, completed = length(records))
            stop(fatal)
        }
        for (b in s8_branches) publish_branch(b)
        rm(dat, ph)
        gc()
    }
    else {
        flow <- read.csv(file.path(out, "raw", "sibship_S8_sample_flow.csv"))
        distribution <- read.csv(file.path(out, "raw", "sibship_S_distribution.csv"))
    }
    review_assert(length(records) == 4746L, "INCOMPLETE_SCAN_NO_FINAL_CORRECTION")
    z <- review_bind(records[keys])
    review_assert(identical(as.character(z$key), keys) && identical(review_key(z$branch, z$phecode), keys), "RESULT_IDENTITY_MISMATCH")
    numeric <- c("N", "cases", "noncases", "beta", "se", "OR", "CI_low", "CI_high", "p", "raw_target_beta", "raw_target_se", 
        "raw_target_p", "registered_M", "full_N", "full_cases", "full_noncases", "full_beta", "full_se", "full_OR", 
        "full_CI_low", "full_CI_high", "full_p", "full_q", "delta_logOR", "elapsed_seconds", "warning_count", "matrix_rank", 
        "matrix_columns")
    for (k in numeric) if (k %in% names(z)) 
        z[[k]] <- as.numeric(z[[k]])
    review_assert(!any(z$family_match == "FALSE", na.rm = TRUE), "FAMILY_NOT_CONFIRMED_NO_CORRECTION")
    z$descriptive_BH_q <- z$descriptive_bonf <- NA_real_
    for (b in s8_branches) {
        idx <- which(z$branch == b)
        valid <- idx[is.finite(z$p[idx]) & z$p[idx] >= 0 & z$p[idx] <= 1]
        review_assert(all(z$planned_eligible[valid] == "TRUE") && length(valid) <= s8_M[b], "P_OUTSIDE_PLANNED_FAMILY")
        z$descriptive_BH_q[valid] <- p.adjust(z$p[valid], method = "BH", n = s8_M[b])
        z$descriptive_bonf[valid] <- p.adjust(z$p[valid], method = "bonferroni", n = s8_M[b])
        bd <- file.path(out, "branches", b)
        dir.create(bd, recursive = TRUE, showWarnings = FALSE)
        if (!file.exists(file.path(bd, "results.csv"))) 
            review_csv(z[idx, ], file.path(bd, "results.csv"))
    }
    for (i in seq_len(nrow(prov))) review_assert(review_sha(prov$path[i]) == prov$sha256[i], "SOURCE_CHANGED_DURING_RUN")
    bundle_path <- file.path(out, "audit", "report_bundle.rds")
    if (file.exists(bundle_path)) 
        bundle <- review_cache_read(bundle_path, sig)
    else {
        stage <- file.path(out, "report_attempts", paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
        cfg_report <- cfg
        cfg_report$out <- stage
        rr <- s8_report(cfg_report, z, flow[flow$branch %in% s8_branches, , drop = FALSE], distribution, full)
        bundle <- data.frame(path = basename(rr$files), source = rr$files, sha256 = vapply(rr$files, review_sha, 
            character(1)))
        review_save_checkpoint(bundle, bundle_path, sig)
    }
    for (i in seq_len(nrow(bundle))) {
        review_assert(review_sha(bundle$source[i]) == bundle$sha256[i], "REPORT_BUNDLE_HASH_MISMATCH")
        dst <- file.path(out, bundle$path[i])
        if (file.exists(dst)) 
            review_assert(review_sha(dst) == bundle$sha256[i], "REPORT_DESTINATION_CONFLICT")
        else review_assert(file.copy(bundle$source[i], dst, overwrite = FALSE), "REPORT_PUBLISH_FAILED")
    }
    done <- list(status = "COMPLETED_DESCRIPTIVE_PENDING_REVIEW", rows = nrow(z), branches = s8_branches, fixed_M = as.list(s8_M), 
        fitted_finite_P = sum(is.finite(z$p)), full_LP_run = FALSE, full_Firth_run = FALSE, signature = sig, finished = review_stamp())
    if (!file.exists(file.path(out, "complete.json"))) 
        jsonlite::write_json(done, file.path(out, "complete.json"), auto_unbox = TRUE, pretty = TRUE)
    ff <- sort(list.files(out, recursive = TRUE, full.names = TRUE))
    ff <- ff[!file.info(ff)$isdir & !grepl("/logs/|run.lock|sha256_manifest", ff)]
    if (!file.exists(file.path(out, "audit", "sha256_manifest.txt"))) 
        writeLines(paste(vapply(ff, review_sha, character(1)), substring(ff, nchar(out) + 2)), file.path(out, "audit", 
            "sha256_manifest.txt"))
    review_event(out, "FINISHED_S8", rows = nrow(z), valid_P = sum(is.finite(z$p)), status = done$status)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        s8_branches <- c("border1n", "border_rank", "border_relative")
        s8_M <- c(border1n = 1362L, border_rank = 1419L, border_relative = 1419L)
        s8_make_sample <- function(full, row = NULL) {
            d <- full$data
            keep <- !is.na(d$S) & d$S <= 8
            before <- d
            d <- droplevels(d[keep, , drop = FALSE])
            removed <- character()
            for (k in names(d)) if (is.factor(before[[k]])) {
                lost <- setdiff(levels(before[[k]]), levels(d[[k]]))
                if (length(lost)) removed <- c(removed, paste0(k, ":", lost))
                if (nrow(d) && length(levels(before[[k]])) > 1L && length(levels(d[[k]])) >= 1L) review_assert(identical(levels(before[[k]])[1], 
                  levels(d[[k]])[1]), paste("REFERENCE_CHANGED", k))
            }
            f <- full$formula
            terms <- attr(terms(f), "term.labels")
            target <- terms[1L]
            constants <- character()
            if (nrow(d)) for (term in setdiff(terms, target)) {
                column <- if (term == "factor(S)") "S" else term
                if (length(unique(d[[column]])) <= 1L) constants <- c(constants, term)
            }
            if (nrow(d) && length(unique(d$S)) > 1 && "factor(S)" %in% terms) review_assert(min(d$S) == min(before$S), 
                "FAMILY_SIZE_REFERENCE_CHANGED")
            list(data = d, formula = f, removed_empty_levels = paste(removed, collapse = ";"), removed_terms = full$removed_terms, 
                new_constants = paste(constants, collapse = ";"), counts = list(N = nrow(d), cases = sum(d$y == 
                  1), noncases = sum(d$y == 0)))
        }
        s8_saved_plan <- function(row) {
            value <- if (row$branch == "border_relative") row$planned_eligible else row$planned
            review_assert(value %in% c("TRUE", "FALSE"), "SAVED_PLAN_FLAG_MISSING")
            value == "TRUE"
        }
        s8_fit_one <- function(row, dat, ph) {
            started <- proc.time()[["elapsed"]]
            b <- row$branch
            p <- row$phecode
            r <- as.list(row[c("branch", "phecode", "description", "category", "groupnum", "color", "sex_rule", 
                "metadata_verified")])
            r$key <- review_key(b, p)
            r$registered_M <- unname(s8_M[b])
            r$planned_eligible <- s8_saved_plan(row)
            for (k in c("N", "cases", "noncases")) r[[paste0("full_", k)]] <- as.integer(row[[k]])
            for (k in c("beta", "se", "p")) r[[paste0("full_", k)]] <- as.numeric(row[[paste0("saved_raw_", k)]])
            for (k in c("OR", "CI_low", "CI_high")) r[[paste0("full_", k)]] <- as.numeric(row[[paste0("saved_", 
                k)]])
            r$full_q <- as.numeric(row$recovered_wald_BH_q)
            r$full_lp_state <- row$separation_assessment
            r$N <- r$cases <- r$noncases <- NA_integer_
            r$beta <- r$se <- r$OR <- r$CI_low <- r$CI_high <- r$p <- NA_real_
            r$raw_target_beta <- r$raw_target_se <- r$raw_target_p <- NA_real_
            r$converged <- NA
            r$matrix_rank <- r$matrix_columns <- NA_integer_
            r$warning <- ""
            r$warning_count <- 0L
            r$status <- "NOT_ESTIMABLE"
            r$error <- ""
            r$family_match <- FALSE
            r$full_sample_hash_match <- FALSE
            r$lp_status <- "NOT_RUN_FULL_LP_NOT_AUTHORIZED"
            r$separation_detected <- NA
            r$target_infinite <- NA
            r$inference_status <- "DESCRIPTIVE_WALD_PENDING_REVIEW"
            r$formal_p <- r$formal_q <- NA_real_
            answer <- tryCatch({
                full <- linear_sample(dat, ph[[p]], b, row$sex_rule)
                for (k in c("N", "cases", "noncases")) review_assert(full$counts[[k]] == as.integer(row[[k]]), 
                  paste("FULL_SAMPLE_COUNT_MISMATCH", k))
                h <- review_hash_ids(full$data$eid)
                savedhash <- if (b == "border_relative") row$sample_hash_sha256 else row$sample_id_digest
                if (!is.na(savedhash) && nzchar(savedhash)) {
                  review_assert(identical(if (b == "border_relative") h$win else h$mac, savedhash), "FULL_SAMPLE_HASH_MISMATCH")
                  r$full_sample_hash_match <- TRUE
                }
                review_assert(identical(review_normalize_formula(review_formula_text(full$formula)), review_normalize_formula(row$actual_formula)), 
                  "FULL_FORMULA_MISMATCH")
                planned <- full$counts$cases >= 20L && full$counts$noncases >= 20L && (b == "border_relative" || 
                  length(unique(full$data[[b]])) >= 2L)
                review_assert(identical(planned, r$planned_eligible), "REGISTERED_TEST_FAMILY_MISMATCH")
                r$family_match <- TRUE
                s <- s8_make_sample(full, row)
                d <- s$data
                for (k in c("N", "cases", "noncases")) r[[k]] <- s$counts[[k]]
                r$sample_hash <- review_hash_ids(d$eid)$win
                r$actual_formula <- review_formula_text(s$formula)
                r$removed_terms <- paste(s$removed_terms, collapse = ";")
                r$removed_empty_levels <- s$removed_empty_levels
                r$effect_unit <- unname(linear_units[b])
                r$S_factor_levels <- paste(sort(unique(d$S)), collapse = ";")
                fac <- names(d)[vapply(d, is.factor, logical(1))]
                r$reference_levels <- paste(vapply(fac, function(k) paste0(k, "=", levels(d[[k]])[1]), character(1)), 
                  collapse = ";")
                r$reference_levels <- paste0(r$reference_levels, ";S=", if (nrow(d)) min(d$S) else NA)
                branch_ok <- attr(dat, "s8_branch_ok")[[b]]
                review_assert(!is.null(branch_ok), "MISSING_BRANCH_GATE")
                threshold <- branch_ok && r$cases >= 20L && r$noncases >= 20L && length(unique(d[[b]])) >= 2L
                r$s8_threshold_eligible <- threshold
                review_assert(!threshold || r$planned_eligible, "S8_ELIGIBLE_OUTSIDE_ORIGINAL_FAMILY")
                if (row$metadata_verified != "TRUE") {
                  r$status <- "PHECODE_METADATA_UNVERIFIED"
                  return(r)
                }
                if (!threshold) {
                  r$status <- "LOW_CASE_OR_CONTROL_COUNT_OR_CONSTANT_TARGET"
                  return(r)
                }
                if (nzchar(s$new_constants)) {
                  r$status <- "NEW_CONSTANT_COVARIATE_NOT_REMOVED"
                  r$error <- s$new_constants
                  return(r)
                }
                mm <- model.matrix(s$formula, d)
                r$matrix_rank <- qr(mm)$rank
                r$matrix_columns <- ncol(mm)
                if (any(!is.finite(mm)) || r$matrix_rank != r$matrix_columns) {
                  r$status <- "RANK_DEFICIENT"
                  r$error <- "No informative column deletion"
                  return(r)
                }
                review_assert(b %in% colnames(mm), "TARGET_COLUMN_MISSING")
                r$design_hash <- digest::digest(list(x = mm, y = d$y), algo = "sha256")
                rm(mm)
                got <- review_capture(glm(s$formula, data = d, family = binomial(link = "logit")))
                r$warning <- paste(got$warnings, collapse = "; ")
                r$warning_count <- length(got$warnings)
                if (inherits(got$value, "error")) {
                  r$status <- "FIT_ERROR"
                  r$error <- conditionMessage(got$value)
                  return(r)
                }
                fit <- got$value
                r$converged <- fit$converged
                r$iterations <- fit$iter
                cf <- summary(fit)$coefficients
                if (b %in% rownames(cf)) {
                  r$raw_target_beta <- unname(cf[b, "Estimate"])
                  r$raw_target_se <- unname(cf[b, "Std. Error"])
                  r$raw_target_p <- unname(cf[b, "Pr(>|z|)"])
                }
                if (!isTRUE(fit$converged)) {
                  r$status <- "NONCONVERGED"
                  return(r)
                }
                if (!b %in% rownames(cf) || !all(is.finite(cf[b, ])) || cf[b, "Std. Error"] <= 0) {
                  r$status <- "TARGET_NOT_IDENTIFIABLE"
                  return(r)
                }
                if (!all(is.finite(coef(fit)))) {
                  r$status <- "NONFINITE_COEFFICIENT"
                  return(r)
                }
                review_assert(nobs(fit) == r$N, "FIT_SAMPLE_CHANGED")
                r$beta <- r$raw_target_beta
                r$se <- r$raw_target_se
                r$p <- r$raw_target_p
                ci <- exp(c(r$beta, r$beta - qnorm(0.975) * r$se, r$beta + qnorm(0.975) * r$se))
                r$OR <- ci[1]
                r$CI_low <- ci[2]
                r$CI_high <- ci[3]
                r$delta_logOR <- r$beta - r$full_beta
                r$status <- if (length(got$warnings)) "OK_WARNING_DESCRIPTIVE_PENDING" else "OK_DESCRIPTIVE_PENDING"
                if (any(!is.finite(ci))) r$status <- "WALD_OVERFLOW_DESCRIPTIVE_PENDING"
                r$p_method <- "two-sided Wald z"
                r$ci_method <- "95% normal Wald; model-based SE"
                r
            }, error = function(e) {
                r$status <<- "FAILED"
                r$error <<- conditionMessage(e)
                r
            })
            answer$elapsed_seconds <- proc.time()[["elapsed"]] - started
            answer$peak_rss_mb <- review_peak_mb()
            answer$finished <- review_stamp()
            answer
        }
        s8_main <- function() {
            args <- commandArgs(TRUE)
            if (length(args) != 3L) stop("Usage: 10.sibship_sensitivity.R MODULE MODE OUT")
            module <- normalizePath(args[1], mustWork = TRUE)
            mode <- args[2]
            out <- args[3]
            if (!mode %in% c("sibship_check", "sibship_focus", "sibship_report")) stop("Unknown stage")
            extra <- Sys.getenv("REVIEW_R_LIB", birth_private_setting("PRIVATE_LIB"))
            if (dir.exists(extra)) .libPaths(unique(c(.libPaths(), extra)))
            for (p in c("digest", "jsonlite", "dplyr", "ggplot2", "ggrepel", "patchwork", "scales")) if (!requireNamespace(p, 
                quietly = TRUE)) stop("Missing dependency ", p)
            if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
                RhpcBLASctl::blas_set_num_threads(1L)
                RhpcBLASctl::omp_set_num_threads(1L)
            }
            source(file.path(module, "frozen", "common.R"), local = .GlobalEnv)
            source(file.path(module, "frozen", "helper.R"), local = .GlobalEnv)
            source(file.path(module, "f", "11.sibship_report.R"), local = .GlobalEnv)
            cfg <- jsonlite::read_json(file.path(module, "config.json"), simplifyVector = TRUE)
            cfg$out <- out
            cfg$module <- module
            cfg$huang_plotting <- file.path(module, "frozen", "huang_plotting.R")
            for (k in c("all_rds", "phewas_rds", "summary_csv")) cfg[[k]] <- Sys.getenv(paste0("S8_", toupper(k)), 
                cfg[[k]])
            for (d in c("raw", "fig", "reports", "audit", "logs", "checkpoints", "branches", "code_snapshot")) dir.create(file.path(out, 
                d), recursive = TRUE, showWarnings = FALSE)
            review_event(out, "START_S8_FULL_PHEWAS", mode = mode, workers = cfg$workers, total_rows = 4746L)
            prov <- do.call(rbind, lapply(names(cfg$expected_sha), function(k) {
                review_assert(file.exists(cfg[[k]]), paste("Missing input", cfg[[k]]))
                review_event(out, "VERIFY_INPUT_SHA", role = k, path = cfg[[k]])
                hash <- review_sha(cfg[[k]])
                review_assert(hash == cfg$expected_sha[[k]], paste("INPUT_SHA_MISMATCH", k))
                data.frame(role = k, path = cfg[[k]], sha256 = hash)
            }))
            files <- sort(list.files(module, recursive = TRUE, full.names = TRUE))
            files <- files[!file.info(files)$isdir]
            cm <- data.frame(path = substring(files, nchar(module) + 2L), sha256 = vapply(files, review_sha, character(1)))
            sig <- digest::digest(list(config = cfg, inputs = prov, code = cm, R = R.version$version.string, package_versions = vapply(c("stats", 
                "digest", "dplyr"), function(p) as.character(packageVersion(p)), character(1))), algo = "sha256")
            owner <- file.path(out, "audit", "ownership.rds")
            if (file.exists(owner)) review_assert(identical(readRDS(owner)$signature, sig), "OUTPUT_SIGNATURE_MISMATCH") else {
                review_rds(list(signature = sig, config = cfg, created = review_stamp()), owner)
                review_csv(prov, file.path(out, "audit", "input_provenance.csv"))
                review_csv(cm, file.path(out, "audit", "code_manifest.csv"))
                capture.output(sessionInfo(), file = file.path(out, "audit", "sessionInfo.txt"))
                for (i in seq_along(files)) {
                  dst <- file.path(out, "code_snapshot", cm$path[i])
                  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
                  review_assert(file.copy(files[i], dst, overwrite = FALSE), "SNAPSHOT_COPY_FAILED")
                }
            }
            if (file.exists(file.path(out, "complete.json")) && file.exists(file.path(out, "audit", "sha256_manifest.txt"))) {
                dd <- jsonlite::read_json(file.path(out, "complete.json"), simplifyVector = TRUE)
                review_assert(identical(dd$signature, sig) && dd$rows == 4746L, "COMPLETE_SIGNATURE_MISMATCH")
                ml <- readLines(file.path(out, "audit", "sha256_manifest.txt"))
                mp <- substring(ml, 66L)
                mh <- substring(ml, 1L, 64L)
                review_assert(all(file.exists(file.path(out, mp))) && all(vapply(file.path(out, mp), review_sha, 
                  character(1)) == mh), "COMPLETED_OUTPUT_HASH_MISMATCH")
                review_event(out, "VERIFIED_ALREADY_COMPLETE", rows = 4746L)
                return(invisible(NULL))
            }
            if (mode != "sibship_report") {
                td <- file.path(out, "audit", paste0("synthetic_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
                log <- file.path(out, "logs", paste0(basename(td), ".log"))
                rc <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(file.path(module, "tests", "test_s8.R")), 
                  shQuote(module), shQuote(td)), stdout = log, stderr = log)
                review_assert(identical(rc, 0L), "S8_SYNTHETIC_TEST_FAILED")
            }
            if (mode == "sibship_check") {
                review_event(out, "CHECK_PASSED")
                return(invisible(NULL))
            }
            full <- read.csv(cfg$summary_csv, colClasses = "character", check.names = FALSE, fileEncoding = "UTF-8")
            registry <- full[full$branch %in% s8_branches, , drop = FALSE]
            for (b in s8_branches) {
                a <- registry[registry$branch == b, ]
                review_assert(nrow(a) == 1582L && !anyDuplicated(a$phecode), "SOURCE_REGISTRY_MISMATCH")
                flags <- vapply(seq_len(nrow(a)), function(i) s8_saved_plan(a[i, ]), logical(1))
                review_assert(sum(flags) == s8_M[b], "FIXED_DENOMINATOR_MISMATCH")
            }
            keys <- review_key(registry$branch, registry$phecode)
            review_assert(!anyDuplicated(keys), "DUPLICATE_KEY")
            records <- list()
            cp <- function(k) file.path(out, "checkpoints", paste0(k, ".rds"))
            for (k in keys) if (file.exists(cp(k))) {
                r <- review_cache_read(cp(k), sig)
                review_assert(identical(r$key, k) && identical(review_key(r$branch, r$phecode), k), "CHECKPOINT_IDENTITY_MISMATCH")
                records[[k]] <- r
            }
            if (mode != "sibship_report") {
                review_event(out, "READ_BASELINE")
                base <- readRDS(cfg$all_rds)
                need <- c("eid", "age", "sex", "ethnicity", "center", "date_attend", "cnt_older", "cnt_brother", 
                  "cnt_sister", "adopted")
                review_assert(all(need %in% names(base)), "MISSING_FIELDS")
                dat <- linear_covars(make_birth_order(as.data.frame(base[need])))
                rm(base)
                gc()
                dat <- dat[, c("eid", linear_exposures, linear_covariates, "S")]
                review_event(out, "READ_ALL_1582_PHENOTYPES")
                ph <- readRDS(cfg$phewas_rds)
                review_assert(!anyNA(ph$id) && !anyDuplicated(ph$id), "PHENOTYPE_ID_ERROR")
                codes <- setdiff(names(ph), "id")
                review_assert(length(codes) == 1582L && !anyDuplicated(codes), "PHENOTYPE_SET_ERROR")
                for (b in s8_branches) review_assert(setequal(codes, registry$phecode[registry$branch == b]), "PHENOTYPE_FAMILY_SET_MISMATCH")
                ix <- match(as.character(dat$eid), as.character(ph$id))
                keep <- !is.na(ix)
                dat <- dat[keep, , drop = FALSE]
                ph <- ph[ix[keep], codes, drop = FALSE]
                gc()
                review_assert(all(vapply(ph, function(y) all(is.na(y) | y %in% c(0, 1)), logical(1))), "OUTCOME_STATE_INVALID")
                distributions <- flows <- list()
                for (b in linear_exposures) {
                  fs <- linear_sample(dat, rep(0, nrow(dat)), b, "both")
                  s <- s8_make_sample(fs)
                  t <- as.data.frame(table(fs$data$S), stringsAsFactors = FALSE)
                  names(t) <- c("S", "N")
                  t$S <- as.integer(t$S)
                  t$branch <- b
                  distributions[[b]] <- t
                  flows[[b]] <- data.frame(branch = b, originalN = nrow(fs$data), retainedN = nrow(s$data), removedN = nrow(fs$data) - 
                    nrow(s$data), retention_pct = 100 * nrow(s$data)/nrow(fs$data), border12_ref_only = b == "border12")
                }
                distribution <- do.call(rbind, distributions)
                flow <- do.call(rbind, flows)
                attr(dat, "s8_branch_ok") <- setNames(as.list(flow$retainedN >= 500L), flow$branch)
                if (!file.exists(file.path(out, "raw", "sibship_S_distribution.csv"))) review_csv(distribution, 
                  file.path(out, "raw", "sibship_S_distribution.csv"))
                if (!file.exists(file.path(out, "raw", "sibship_S8_sample_flow.csv"))) review_csv(flow, file.path(out, 
                  "raw", "sibship_S8_sample_flow.csv"))
                if (!file.exists(file.path(out, "audit", "source_registry.csv"))) review_csv(registry, file.path(out, 
                  "audit", "source_registry.csv"))
                start <- as.numeric(Sys.time())
                prior_runtime <- sum(vapply(records, function(r) if (is.null(r$elapsed_seconds)) 0 else r$elapsed_seconds, 
                  numeric(1)), na.rm = TRUE)/cfg$workers
                jobs <- list()
                pending <- keys[!keys %in% names(records)]
                nexti <- 1L
                fatal <- ""
                failed_record <- function(key, status, error, elapsed) {
                  row <- registry[match(key, keys), , drop = FALSE]
                  r <- as.list(row)
                  r$key <- key
                  r$full_N <- as.integer(row$N)
                  r$full_cases <- as.integer(row$cases)
                  r$full_noncases <- as.integer(row$noncases)
                  r$N <- r$cases <- r$noncases <- NA_integer_
                  for (k in c("beta", "se", "OR", "CI_low", "CI_high", "p", "raw_target_beta", "raw_target_se", 
                    "raw_target_p")) r[[k]] <- NA_real_
                  r$status <- status
                  r$error <- error
                  r$elapsed_seconds <- elapsed
                  r$family_match <- NA
                  r$planned_eligible <- s8_saved_plan(row)
                  r$registered_M <- unname(s8_M[row$branch])
                  r
                }
                publish_branch <- function(b) {
                  bk <- keys[registry$branch == b]
                  dest <- file.path(out, "branches", b)
                  if (!all(bk %in% names(records)) || file.exists(file.path(dest, "complete.json"))) return(invisible(NULL))
                  zz <- review_bind(records[bk])
                  zz$p <- as.numeric(zz$p)
                  review_assert(!any(zz$family_match == "FALSE", na.rm = TRUE), "BRANCH_FAMILY_MISMATCH")
                  valid <- is.finite(zz$p) & zz$p >= 0 & zz$p <= 1
                  review_assert(all(zz$planned_eligible[valid] == "TRUE"), "P_OUTSIDE_FAMILY")
                  zz$descriptive_BH_q <- zz$descriptive_bonf <- NA_real_
                  zz$descriptive_BH_q[valid] <- p.adjust(zz$p[valid], "BH", n = s8_M[b])
                  zz$descriptive_bonf[valid] <- p.adjust(zz$p[valid], "bonferroni", n = s8_M[b])
                  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
                  if (!file.exists(file.path(dest, "results.csv"))) review_csv(zz, file.path(dest, "results.csv"))
                  if (!file.exists(file.path(dest, "results.rds"))) review_rds(zz, file.path(dest, "results.rds"))
                  jsonlite::write_json(list(status = "COMPLETED_DESCRIPTIVE_PENDING", branch = b, rows = 1582L, 
                    M = unname(s8_M[b]), finite_P = sum(valid), signature = sig, finished = review_stamp()), file.path(dest, 
                    "complete.json"), pretty = TRUE, auto_unbox = TRUE)
                  review_event(out, "BRANCH_COMPLETE", branch = b, rows = 1582L, finite_P = sum(valid), M = unname(s8_M[b]))
                }
                while (nexti <= length(pending) || length(jobs)) {
                  while (nexti <= length(pending) && length(jobs) < cfg$workers && !nzchar(fatal)) {
                    k <- pending[nexti]
                    row <- registry[match(k, keys), , drop = FALSE]
                    job <- parallel::mcparallel(s8_fit_one(row, dat, ph), mc.set.seed = FALSE, silent = TRUE)
                    jobs[[as.character(job$pid)]] <- list(job = job, key = k, start = as.numeric(Sys.time()))
                    review_event(out, "MODEL_START", key = k, pid_worker = job$pid, completed = length(records), 
                      total = 4746L)
                    nexti <- nexti + 1L
                  }
                  for (pid in names(jobs)) {
                    j <- jobs[[pid]]
                    v <- parallel::mccollect(j$job, wait = FALSE)
                    elapsed <- as.numeric(Sys.time()) - j$start
                    if (is.null(v) && (elapsed > cfg$per_model_seconds || as.numeric(Sys.time()) - start + prior_runtime > 
                      cfg$total_seconds)) {
                      tools::pskill(j$job$pid, 15L)
                      Sys.sleep(1)
                      v <- parallel::mccollect(j$job, wait = FALSE)
                      if (is.null(v)) {
                        tools::pskill(j$job$pid, 9L)
                        parallel::mccollect(j$job, wait = TRUE)
                      }
                      v <- list(failed_record(j$key, "TIMEOUT", "Owned S8 worker exceeded resource limit", elapsed))
                    }
                    if (!is.null(v)) {
                      r <- v[[1L]]
                      if (!is.list(r) || inherits(r, "try-error")) {
                        r <- failed_record(j$key, "WORKER_FAILED", "Owned child returned invalid result", elapsed)
                        fatal <- "WORKER_FAILURE"
                      }
                      review_assert(identical(r$key, j$key) && identical(review_key(r$branch, r$phecode), j$key), 
                        "WORKER_IDENTITY_MISMATCH")
                      if (is.null(r$elapsed_seconds)) r$elapsed_seconds <- elapsed
                      review_save_checkpoint(r, cp(j$key), sig)
                      records[[j$key]] <- r
                      jobs[[pid]] <- NULL
                      review_event(out, "MODEL_END", key = j$key, status = r$status, error = r$error, completed = length(records), 
                        total = 4746L, seconds = r$elapsed_seconds)
                      if (grepl("MISMATCH|INVALID|cannot allocate", r$error)) fatal <- paste("SYSTEMIC_CHECK_FAILED", 
                        j$key, r$error)
                      if (!nzchar(fatal)) for (b in s8_branches) publish_branch(b)
                    }
                  }
                  if (as.numeric(Sys.time()) - start + prior_runtime > cfg$total_seconds) fatal <- "TOTAL_BUDGET_EXHAUSTED"
                  if (nzchar(fatal) && !length(jobs)) break
                  Sys.sleep(0.2)
                }
                if (nzchar(fatal)) {
                  review_event(out, "STOPPED", reason = fatal, completed = length(records))
                  stop(fatal)
                }
                for (b in s8_branches) publish_branch(b)
                rm(dat, ph)
                gc()
            } else {
                flow <- read.csv(file.path(out, "raw", "sibship_S8_sample_flow.csv"))
                distribution <- read.csv(file.path(out, "raw", "sibship_S_distribution.csv"))
            }
            review_assert(length(records) == 4746L, "INCOMPLETE_SCAN_NO_FINAL_CORRECTION")
            z <- review_bind(records[keys])
            review_assert(identical(as.character(z$key), keys) && identical(review_key(z$branch, z$phecode), keys), 
                "RESULT_IDENTITY_MISMATCH")
            numeric <- c("N", "cases", "noncases", "beta", "se", "OR", "CI_low", "CI_high", "p", "raw_target_beta", 
                "raw_target_se", "raw_target_p", "registered_M", "full_N", "full_cases", "full_noncases", "full_beta", 
                "full_se", "full_OR", "full_CI_low", "full_CI_high", "full_p", "full_q", "delta_logOR", "elapsed_seconds", 
                "warning_count", "matrix_rank", "matrix_columns")
            for (k in numeric) if (k %in% names(z)) z[[k]] <- as.numeric(z[[k]])
            review_assert(!any(z$family_match == "FALSE", na.rm = TRUE), "FAMILY_NOT_CONFIRMED_NO_CORRECTION")
            z$descriptive_BH_q <- z$descriptive_bonf <- NA_real_
            for (b in s8_branches) {
                idx <- which(z$branch == b)
                valid <- idx[is.finite(z$p[idx]) & z$p[idx] >= 0 & z$p[idx] <= 1]
                review_assert(all(z$planned_eligible[valid] == "TRUE") && length(valid) <= s8_M[b], "P_OUTSIDE_PLANNED_FAMILY")
                z$descriptive_BH_q[valid] <- p.adjust(z$p[valid], method = "BH", n = s8_M[b])
                z$descriptive_bonf[valid] <- p.adjust(z$p[valid], method = "bonferroni", n = s8_M[b])
                bd <- file.path(out, "branches", b)
                dir.create(bd, recursive = TRUE, showWarnings = FALSE)
                if (!file.exists(file.path(bd, "results.csv"))) review_csv(z[idx, ], file.path(bd, "results.csv"))
            }
            for (i in seq_len(nrow(prov))) review_assert(review_sha(prov$path[i]) == prov$sha256[i], "SOURCE_CHANGED_DURING_RUN")
            bundle_path <- file.path(out, "audit", "report_bundle.rds")
            if (file.exists(bundle_path)) bundle <- review_cache_read(bundle_path, sig) else {
                stage <- file.path(out, "report_attempts", paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
                cfg_report <- cfg
                cfg_report$out <- stage
                rr <- s8_report(cfg_report, z, flow[flow$branch %in% s8_branches, , drop = FALSE], distribution, 
                  full)
                bundle <- data.frame(path = basename(rr$files), source = rr$files, sha256 = vapply(rr$files, review_sha, 
                  character(1)))
                review_save_checkpoint(bundle, bundle_path, sig)
            }
            for (i in seq_len(nrow(bundle))) {
                review_assert(review_sha(bundle$source[i]) == bundle$sha256[i], "REPORT_BUNDLE_HASH_MISMATCH")
                dst <- file.path(out, bundle$path[i])
                if (file.exists(dst)) review_assert(review_sha(dst) == bundle$sha256[i], "REPORT_DESTINATION_CONFLICT") else review_assert(file.copy(bundle$source[i], 
                  dst, overwrite = FALSE), "REPORT_PUBLISH_FAILED")
            }
            done <- list(status = "COMPLETED_DESCRIPTIVE_PENDING_REVIEW", rows = nrow(z), branches = s8_branches, 
                fixed_M = as.list(s8_M), fitted_finite_P = sum(is.finite(z$p)), full_LP_run = FALSE, full_Firth_run = FALSE, 
                signature = sig, finished = review_stamp())
            if (!file.exists(file.path(out, "complete.json"))) jsonlite::write_json(done, file.path(out, "complete.json"), 
                auto_unbox = TRUE, pretty = TRUE)
            ff <- sort(list.files(out, recursive = TRUE, full.names = TRUE))
            ff <- ff[!file.info(ff)$isdir & !grepl("/logs/|run.lock|sha256_manifest", ff)]
            if (!file.exists(file.path(out, "audit", "sha256_manifest.txt"))) writeLines(paste(vapply(ff, review_sha, 
                character(1)), substring(ff, nchar(out) + 2)), file.path(out, "audit", "sha256_manifest.txt"))
            review_event(out, "FINISHED_S8", rows = nrow(z), valid_P = sum(is.finite(z$p)), status = done$status)
        }
        if (sys.nframe() == 0L) s8_main()
    }), envir = context)
}

