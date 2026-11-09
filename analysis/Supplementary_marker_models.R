# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

v2_feature_dictionary <- function(root) {
    groups <- c(clinical = "clinical", blood_pressure = "clinical", nmr = "nmr", olink = "olink", cmr = "imaging")
    data.table::rbindlist(lapply(names(groups), function(k) {
        z <- data.table::fread(file.path(root, "config", paste0(k, "_features.csv")), colClasses = "character")
        z[, `:=`(module = k, family = groups[[k]], mediator_id = paste(k, feature_id, sep = "__"), indexed = FALSE)]
        if (k == "cmr") {
            zz <- data.table::copy(z[tolower(bsa_index) == "true"])
            zz[, `:=`(mediator_id = paste0(mediator_id, "_bsa"), indexed = TRUE, label = paste0(label, " / BSA"))]
            z <- data.table::rbindlist(list(z, zz), fill = TRUE)
        }
        z
    }), fill = TRUE)
}

v2_feature_context <- function(ctx) {
    e <- new.env(parent = globalenv())
    e$root <- ctx$root
    e$out <- ctx$cfg$sources$phase2
    e$feature_cache <- new.env()
    e$num <- bi_num
    e$assert <- v2_assert
    e$cp_load <- function(path, strict = TRUE) v2_unwrap(path)
    e$read_json <- jsonlite::read_json
    e$dict <- function() v2_feature_dictionary(ctx$root)
    sys.source(file.path(ctx$root, "vendor/features.R"), envir = e)
    e$core <- v2_unwrap(file.path(e$out, "cache/core.rds"))
    co <- v2_read(ctx, "cohort")
    d <- co$data[co$data$.exposure_eligible, , drop = FALSE]
    ix <- match(d$eid, e$core$data$eid)
    v2_assert(!anyNA(ix), "FROZEN_FEATURE_CACHE_MISSING_NEW_ELIGIBLE_IDS; projection must be rebuilt, not silently filtered")
    for (k in c("S", "border_relative", "age", "sex", "ethnicity", "center", "assessment_year")) v2_assert(isTRUE(all.equal(as.character(d[[k]]), 
        as.character(e$core$data[[k]][ix]))), paste("FEATURE_COHORT_MISMATCH", k))
    proof <- jsonlite::fromJSON(file.path(e$out, "cache/projected/complete.json"))
    for (i in seq_len(nrow(proof$files))) {
        p <- file.path(e$out, "cache/projected", paste0(proof$files$module[i], ".csv"))
        v2_assert(identical(v2_filehash(p), proof$files$sha256[i]), paste("PROJECTED_CACHE_HASH_MISMATCH", p))
    }
    source <- data.table::fread(file.path(e$out, "INPUT_SHA256.csv"), data.table = FALSE)
    st <- file.info(source$path)
    v2_assert(all(st$size == source$bytes) & all(abs(as.numeric(st$mtime) - source$mtime) < 1e-04), "SOURCE_CHANGED_SINCE_PROJECTION")
    v2_csv(ctx, "INDEPENDENT_FEATURE_CACHE_REUSE.csv", data.frame(path = source$path, sha256 = source$sha256, verified = "source size/mtime plus authenticated projected file content; no historical disease selection", 
        model_results_reused = FALSE))
    e
}

v2_medication_flags <- function(ctx) {
    v2_cp(ctx, "baseline_medication", function() {
        co <- v2_read(ctx, "cohort")
        p <- co$manifest$path[co$manifest$role == "baseline"]
        raw <- data.table::fread(p, data.table = FALSE, colClasses = c(eid = "character"))
        contract <- list(verified = TRUE, baseline_only_verified = TRUE, arrays_verified = TRUE, medication_dictionary_verified = TRUE, 
            evidence = "Same authenticated baseline extract and conservative frozen diabetes helper as October16 incident adapter", 
            algorithm = "Eastwood2016_Fig2_S1_conservative_baseline_v1")
        dm <- ev_diabetes_baseline(raw, co$data$eid, contract)
        v2_assert(isTRUE(dm$verified), "MEDICATION_DEFINITION_GATE_FAILED")
        f <- dm$flags
        known <- f$raw_linked & f$medication_complete & f$touchscreen_insulin_known & !f$medication_uncertain
        used <- f$insulin_TS | f$insulin_NI | f$metformin_NI | f$nonmetformin_oral_NI
        touchscreen <- function(x, allowed, target) vapply(as.character(x), function(t) {
            if (is.na(t) || !nzchar(t)) 
                return(NA_integer_)
            v <- suppressWarnings(as.numeric(strsplit(trimws(t), "[|,;[:space:]]+", perl = TRUE)[[1]]))
            if (anyNA(v) || any(!v %in% c(allowed, -7)) || (-7 %in% v && length(v) > 1)) 
                return(NA_integer_)
            as.integer(target %in% v)
        }, integer(1))
        j <- match(co$data$eid, as.character(raw$eid))
        bp <- rep(NA_integer_, nrow(co$data))
        if (all(c("p6153_i0", "p6177_i0") %in% names(raw))) {
            female <- touchscreen(raw$p6153_i0, 1:5, 2)
            male <- touchscreen(raw$p6177_i0, 1:3, 2)
            bp <- ifelse(as.character(co$data$sex) == "0", female[j], male[j])
        }
        data.frame(eid = co$data$eid, DM_algorithm = dm$status, glucose_drug_free = ifelse(known, as.integer(!used), 
            NA_integer_), bp_drug_free = ifelse(is.na(bp), NA_integer_, as.integer(bp == 0)))
    })
}

v2_continuous_one <- function(ctx, e, row, sensitivity = "primary") {
    id <- row$mediator_id
    key <- paste("continuous", id, sensitivity, sep = "_")
    module <- row$module
    fam <- if (module == "clinical") 
        if (row$feature_id %in% c("glucose", "hba1c")) 
            "GLYCAEMIA_2"
        else "BIOCHEMISTRY_7"
    else if (module == "blood_pressure") 
        if (row$feature_id %in% c("SBP", "DBP", "sbp", "dbp")) 
            "BP_MAIN_2"
        else "BP_DERIVED_2"
    else if (module == "cmr") 
        if (isTRUE(row$indexed)) 
            "CMR_INDEXED_7"
        else "CMR_PRIMARY_82"
    else toupper(paste0(module, "_", if (module == "nmr") 
        251
    else 2923))
    if (sensitivity != "primary") 
        fam <- paste(fam, sensitivity, sep = "_")
    v2_cp(ctx, paste0("continuous_detail_", key), function() {
        r <- v2_row(module, key, fam, feature = id)
        r$label <- row$label
        r$unit <- row$unit
        r$analysis <- sensitivity
        tryCatch({
            z <- e$load_feature(e$core, row)
            d <- v2_s8(z$data)
            co <- v2_read(ctx, "cohort")
            d <- d[d$eid %in% co$data$eid, , drop = FALSE]
            if (sensitivity == "no_prescan_IHD_T2D_351") {
                inc <- v2_read(ctx, "incident")
                j <- match(d$eid, inc$data$eid)
                keep <- !is.na(j) & !is.na(d$measurement_date)
                for (k in c("411", "250.2", "351")) {
                  b <- inc$baseline[[k]][j]
                  ev <- inc$events[[k]][j]
                  keep <- keep & !is.na(b) & b == 0 & (is.na(ev) | ev > as.Date(d$measurement_date))
                }
                keep[is.na(keep)] <- FALSE
                d <- d[keep, , drop = FALSE]
            }
            else if (sensitivity != "primary") {
                flags <- v2_medication_flags(ctx)
                j <- match(d$eid, flags$eid)
                if (sensitivity == "no_glucose_drug") 
                  keep <- flags$glucose_drug_free[j] == 1
                if (sensitivity == "no_BP_drug") 
                  keep <- flags$bp_drug_free[j] == 1
                if (sensitivity == "no_baseline_DM") 
                  keep <- flags$DM_algorithm[j] == "no_diabetes"
                keep[is.na(keep)] <- FALSE
                d <- d[keep, , drop = FALSE]
            }
            d$y <- d$marker
            s <- bi_sample(d, bi_formula("border_relative", z$extra), "border_relative")
            sd_y <- sd(s$data$y)
            v2_assert(is.finite(sd_y) && sd_y > 0, "CONSTANT_CONTINUOUS_OUTCOME")
            a <- v2_model(ctx, s$data, module, key, fam, feature = id, extra = z$extra, type = "gaussian")
            orig <- a
            if (row$transform == "log") {
                raw_data <- s$data
                raw_data$y <- raw_data$marker_original
                orig <- v2_model(ctx, raw_data, module, paste0(key, "_original_scale"), fam, feature = id, extra = z$extra, 
                  type = "gaussian")
                v2_assert(orig$N == a$N, "ORIGINAL_SCALE_SAMPLE_CHANGED")
            }
            a$raw_beta <- orig$beta
            a$raw_CI_low <- orig$CI_low
            a$raw_CI_high <- orig$CI_high
            a$raw_scale_p <- orig$p
            a$raw_scale_status <- orig$status
            a$sd_outcome <- sd_y
            a$standardized_beta <- a$beta/sd_y
            a$standardized_CI_low <- a$CI_low/sd_y
            a$standardized_CI_high <- a$CI_high/sd_y
            a$label <- row$label
            a$unit <- row$unit
            a$analysis <- sensitivity
            a$source_qc <- z$source_qc
            a$transform <- row$transform
            a$scale <- if (row$transform == "log") 
                "LOG_OUTCOME_BETA_WITH_RAW_UNIT_AND_SD_COMPANIONS"
            else "RAW_BETA_0_to_1_WITH_SD_COMPANION"
            r <- a
        }, error = function(err) {
            r$reason <<- conditionMessage(err)
            r$status <<- if (grepl("DEPENDENCY_NOT_READY|SOURCE_|CACHE_|FIELD_|MISSING", conditionMessage(err))) 
                "BLOCKED_INPUT"
            else "EXECUTION_OR_MODEL_ERROR"
        })
        r
    })
}

v2_clinical_omics <- function(ctx) {
    e <- v2_feature_context(ctx)
    dd <- v2_feature_dictionary(ctx$root)
    jobs <- lapply(seq_len(nrow(dd)), function(i) list(row = as.list(dd[i]), sensitivity = "primary"))
    for (i in which(dd$feature_id %in% c("glucose", "hba1c") & dd$module == "clinical")) for (s in c("no_baseline_DM", 
        "no_glucose_drug")) jobs[[length(jobs) + 1]] <- list(row = as.list(dd[i]), sensitivity = s)
    for (i in which(dd$module == "blood_pressure")) jobs[[length(jobs) + 1]] <- list(row = as.list(dd[i]), sensitivity = "no_BP_drug")
    family_for <- function(j) {
        r <- j$row
        m <- r$module
        f <- if (m == "clinical") 
            if (r$feature_id %in% c("glucose", "hba1c")) 
                "GLYCAEMIA_2"
            else "BIOCHEMISTRY_7"
        else if (m == "blood_pressure") 
            if (r$feature_id %in% c("SBP", "DBP", "sbp", "dbp")) 
                "BP_MAIN_2"
            else "BP_DERIVED_2"
        else if (m == "cmr") 
            if (isTRUE(r$indexed)) 
                "CMR_INDEXED_7"
            else "CMR_PRIMARY_82"
        else toupper(paste0(m, "_", if (m == "nmr") 
            251
        else 2923))
        if (j$sensitivity != "primary") 
            f <- paste(f, j$sensitivity, sep = "_")
        f
    }
    reg <- bi_bind(lapply(jobs, function(j) data.frame(test_id = paste("continuous", j$row$mediator_id, j$sensitivity, 
        sep = "_"), family_id = family_for(j), in_family = TRUE)))
    v2_register(ctx, "continuous", reg)
    v2_medication_flags(ctx)
    chunks <- split(seq_along(jobs), ceiling(seq_along(jobs)/25L))
    out <- v2_parallel(ctx, chunks, function(ii) bi_bind(lapply(ii, function(i) v2_continuous_one(ctx, e, jobs[[i]]$row, 
        jobs[[i]]$sensitivity))))
    tab <- v2_adjust(bi_bind(out), reg)
    v2_csv(ctx, "CLINICAL_OMICS_RESULTS.csv", tab)
    v2_cp(ctx, "results_continuous", function() tab)
    list(status = "PROCESSED", positions = nrow(tab), valid_p = sum(is.finite(tab$p)))
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        v2_feature_dictionary <- function(root) {
            groups <- c(clinical = "clinical", blood_pressure = "clinical", nmr = "nmr", olink = "olink", cmr = "imaging")
            data.table::rbindlist(lapply(names(groups), function(k) {
                z <- data.table::fread(file.path(root, "config", paste0(k, "_features.csv")), colClasses = "character")
                z[, `:=`(module = k, family = groups[[k]], mediator_id = paste(k, feature_id, sep = "__"), indexed = FALSE)]
                if (k == "cmr") {
                  zz <- data.table::copy(z[tolower(bsa_index) == "true"])
                  zz[, `:=`(mediator_id = paste0(mediator_id, "_bsa"), indexed = TRUE, label = paste0(label, " / BSA"))]
                  z <- data.table::rbindlist(list(z, zz), fill = TRUE)
                }
                z
            }), fill = TRUE)
        }
        v2_feature_context <- function(ctx) {
            e <- new.env(parent = globalenv())
            e$root <- ctx$root
            e$out <- ctx$cfg$sources$phase2
            e$feature_cache <- new.env()
            e$num <- bi_num
            e$assert <- v2_assert
            e$cp_load <- function(path, strict = TRUE) v2_unwrap(path)
            e$read_json <- jsonlite::read_json
            e$dict <- function() v2_feature_dictionary(ctx$root)
            sys.source(file.path(ctx$root, "vendor/features.R"), envir = e)
            e$core <- v2_unwrap(file.path(e$out, "cache/core.rds"))
            co <- v2_read(ctx, "cohort")
            d <- co$data[co$data$.exposure_eligible, , drop = FALSE]
            ix <- match(d$eid, e$core$data$eid)
            v2_assert(!anyNA(ix), "FROZEN_FEATURE_CACHE_MISSING_NEW_ELIGIBLE_IDS; projection must be rebuilt, not silently filtered")
            for (k in c("S", "border_relative", "age", "sex", "ethnicity", "center", "assessment_year")) v2_assert(isTRUE(all.equal(as.character(d[[k]]), 
                as.character(e$core$data[[k]][ix]))), paste("FEATURE_COHORT_MISMATCH", k))
            proof <- jsonlite::fromJSON(file.path(e$out, "cache/projected/complete.json"))
            for (i in seq_len(nrow(proof$files))) {
                p <- file.path(e$out, "cache/projected", paste0(proof$files$module[i], ".csv"))
                v2_assert(identical(v2_filehash(p), proof$files$sha256[i]), paste("PROJECTED_CACHE_HASH_MISMATCH", 
                  p))
            }
            source <- data.table::fread(file.path(e$out, "INPUT_SHA256.csv"), data.table = FALSE)
            st <- file.info(source$path)
            v2_assert(all(st$size == source$bytes) & all(abs(as.numeric(st$mtime) - source$mtime) < 1e-04), "SOURCE_CHANGED_SINCE_PROJECTION")
            v2_csv(ctx, "INDEPENDENT_FEATURE_CACHE_REUSE.csv", data.frame(path = source$path, sha256 = source$sha256, 
                verified = "source size/mtime plus authenticated projected file content; no historical disease selection", 
                model_results_reused = FALSE))
            e
        }
        v2_medication_flags <- function(ctx) {
            v2_cp(ctx, "baseline_medication", function() {
                co <- v2_read(ctx, "cohort")
                p <- co$manifest$path[co$manifest$role == "baseline"]
                raw <- data.table::fread(p, data.table = FALSE, colClasses = c(eid = "character"))
                contract <- list(verified = TRUE, baseline_only_verified = TRUE, arrays_verified = TRUE, medication_dictionary_verified = TRUE, 
                  evidence = "Same authenticated baseline extract and conservative frozen diabetes helper as October16 incident adapter", 
                  algorithm = "Eastwood2016_Fig2_S1_conservative_baseline_v1")
                dm <- ev_diabetes_baseline(raw, co$data$eid, contract)
                v2_assert(isTRUE(dm$verified), "MEDICATION_DEFINITION_GATE_FAILED")
                f <- dm$flags
                known <- f$raw_linked & f$medication_complete & f$touchscreen_insulin_known & !f$medication_uncertain
                used <- f$insulin_TS | f$insulin_NI | f$metformin_NI | f$nonmetformin_oral_NI
                touchscreen <- function(x, allowed, target) vapply(as.character(x), function(t) {
                  if (is.na(t) || !nzchar(t)) return(NA_integer_)
                  v <- suppressWarnings(as.numeric(strsplit(trimws(t), "[|,;[:space:]]+", perl = TRUE)[[1]]))
                  if (anyNA(v) || any(!v %in% c(allowed, -7)) || (-7 %in% v && length(v) > 1)) return(NA_integer_)
                  as.integer(target %in% v)
                }, integer(1))
                j <- match(co$data$eid, as.character(raw$eid))
                bp <- rep(NA_integer_, nrow(co$data))
                if (all(c("p6153_i0", "p6177_i0") %in% names(raw))) {
                  female <- touchscreen(raw$p6153_i0, 1:5, 2)
                  male <- touchscreen(raw$p6177_i0, 1:3, 2)
                  bp <- ifelse(as.character(co$data$sex) == "0", female[j], male[j])
                }
                data.frame(eid = co$data$eid, DM_algorithm = dm$status, glucose_drug_free = ifelse(known, as.integer(!used), 
                  NA_integer_), bp_drug_free = ifelse(is.na(bp), NA_integer_, as.integer(bp == 0)))
            })
        }
        v2_continuous_one <- function(ctx, e, row, sensitivity = "primary") {
            id <- row$mediator_id
            key <- paste("continuous", id, sensitivity, sep = "_")
            module <- row$module
            fam <- if (module == "clinical") if (row$feature_id %in% c("glucose", "hba1c")) "GLYCAEMIA_2" else "BIOCHEMISTRY_7" else if (module == 
                "blood_pressure") if (row$feature_id %in% c("SBP", "DBP", "sbp", "dbp")) "BP_MAIN_2" else "BP_DERIVED_2" else if (module == 
                "cmr") if (isTRUE(row$indexed)) "CMR_INDEXED_7" else "CMR_PRIMARY_82" else toupper(paste0(module, 
                "_", if (module == "nmr") 251 else 2923))
            if (sensitivity != "primary") fam <- paste(fam, sensitivity, sep = "_")
            v2_cp(ctx, paste0("continuous_detail_", key), function() {
                r <- v2_row(module, key, fam, feature = id)
                r$label <- row$label
                r$unit <- row$unit
                r$analysis <- sensitivity
                tryCatch({
                  z <- e$load_feature(e$core, row)
                  d <- v2_s8(z$data)
                  co <- v2_read(ctx, "cohort")
                  d <- d[d$eid %in% co$data$eid, , drop = FALSE]
                  if (sensitivity == "no_prescan_IHD_T2D_351") {
                    inc <- v2_read(ctx, "incident")
                    j <- match(d$eid, inc$data$eid)
                    keep <- !is.na(j) & !is.na(d$measurement_date)
                    for (k in c("411", "250.2", "351")) {
                      b <- inc$baseline[[k]][j]
                      ev <- inc$events[[k]][j]
                      keep <- keep & !is.na(b) & b == 0 & (is.na(ev) | ev > as.Date(d$measurement_date))
                    }
                    keep[is.na(keep)] <- FALSE
                    d <- d[keep, , drop = FALSE]
                  } else if (sensitivity != "primary") {
                    flags <- v2_medication_flags(ctx)
                    j <- match(d$eid, flags$eid)
                    if (sensitivity == "no_glucose_drug") keep <- flags$glucose_drug_free[j] == 1
                    if (sensitivity == "no_BP_drug") keep <- flags$bp_drug_free[j] == 1
                    if (sensitivity == "no_baseline_DM") keep <- flags$DM_algorithm[j] == "no_diabetes"
                    keep[is.na(keep)] <- FALSE
                    d <- d[keep, , drop = FALSE]
                  }
                  d$y <- d$marker
                  s <- bi_sample(d, bi_formula("border_relative", z$extra), "border_relative")
                  sd_y <- sd(s$data$y)
                  v2_assert(is.finite(sd_y) && sd_y > 0, "CONSTANT_CONTINUOUS_OUTCOME")
                  a <- v2_model(ctx, s$data, module, key, fam, feature = id, extra = z$extra, type = "gaussian")
                  orig <- a
                  if (row$transform == "log") {
                    raw_data <- s$data
                    raw_data$y <- raw_data$marker_original
                    orig <- v2_model(ctx, raw_data, module, paste0(key, "_original_scale"), fam, feature = id, 
                      extra = z$extra, type = "gaussian")
                    v2_assert(orig$N == a$N, "ORIGINAL_SCALE_SAMPLE_CHANGED")
                  }
                  a$raw_beta <- orig$beta
                  a$raw_CI_low <- orig$CI_low
                  a$raw_CI_high <- orig$CI_high
                  a$raw_scale_p <- orig$p
                  a$raw_scale_status <- orig$status
                  a$sd_outcome <- sd_y
                  a$standardized_beta <- a$beta/sd_y
                  a$standardized_CI_low <- a$CI_low/sd_y
                  a$standardized_CI_high <- a$CI_high/sd_y
                  a$label <- row$label
                  a$unit <- row$unit
                  a$analysis <- sensitivity
                  a$source_qc <- z$source_qc
                  a$transform <- row$transform
                  a$scale <- if (row$transform == "log") "LOG_OUTCOME_BETA_WITH_RAW_UNIT_AND_SD_COMPANIONS" else "RAW_BETA_0_to_1_WITH_SD_COMPANION"
                  r <- a
                }, error = function(err) {
                  r$reason <<- conditionMessage(err)
                  r$status <<- if (grepl("DEPENDENCY_NOT_READY|SOURCE_|CACHE_|FIELD_|MISSING", conditionMessage(err))) "BLOCKED_INPUT" else "EXECUTION_OR_MODEL_ERROR"
                })
                r
            })
        }
        v2_clinical_omics <- function(ctx) {
            e <- v2_feature_context(ctx)
            dd <- v2_feature_dictionary(ctx$root)
            jobs <- lapply(seq_len(nrow(dd)), function(i) list(row = as.list(dd[i]), sensitivity = "primary"))
            for (i in which(dd$feature_id %in% c("glucose", "hba1c") & dd$module == "clinical")) for (s in c("no_baseline_DM", 
                "no_glucose_drug")) jobs[[length(jobs) + 1]] <- list(row = as.list(dd[i]), sensitivity = s)
            for (i in which(dd$module == "blood_pressure")) jobs[[length(jobs) + 1]] <- list(row = as.list(dd[i]), 
                sensitivity = "no_BP_drug")
            family_for <- function(j) {
                r <- j$row
                m <- r$module
                f <- if (m == "clinical") if (r$feature_id %in% c("glucose", "hba1c")) "GLYCAEMIA_2" else "BIOCHEMISTRY_7" else if (m == 
                  "blood_pressure") if (r$feature_id %in% c("SBP", "DBP", "sbp", "dbp")) "BP_MAIN_2" else "BP_DERIVED_2" else if (m == 
                  "cmr") if (isTRUE(r$indexed)) "CMR_INDEXED_7" else "CMR_PRIMARY_82" else toupper(paste0(m, "_", 
                  if (m == "nmr") 251 else 2923))
                if (j$sensitivity != "primary") f <- paste(f, j$sensitivity, sep = "_")
                f
            }
            reg <- bi_bind(lapply(jobs, function(j) data.frame(test_id = paste("continuous", j$row$mediator_id, 
                j$sensitivity, sep = "_"), family_id = family_for(j), in_family = TRUE)))
            v2_register(ctx, "continuous", reg)
            v2_medication_flags(ctx)
            chunks <- split(seq_along(jobs), ceiling(seq_along(jobs)/25L))
            out <- v2_parallel(ctx, chunks, function(ii) bi_bind(lapply(ii, function(i) v2_continuous_one(ctx, 
                e, jobs[[i]]$row, jobs[[i]]$sensitivity))))
            tab <- v2_adjust(bi_bind(out), reg)
            v2_csv(ctx, "CLINICAL_OMICS_RESULTS.csv", tab)
            v2_cp(ctx, "results_continuous", function() tab)
            list(status = "PROCESSED", positions = nrow(tab), valid_p = sum(is.finite(tab$p)))
        }
    }), envir = context)
}

