# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

bi_c_exposures <- function() c("border12", "border1n", "border_rank", "border_relative")

bi_c_codes <- function(ctx = NULL) {
    fixed <- c("411", "411.3", "411.4", "351", "250", "250.2")
    if (is.null(ctx)) 
        return(fixed)
    configured <- as.character(unlist(ctx$cfg$clinical$core))
    if (!identical(configured, fixed)) 
        stop("FIXED_6_PHECODE_REGISTRY_MISMATCH")
    configured
}

bi_c_failure <- function(exposure, status, reason, formula = NA_character_) {
    data.frame(exposure = exposure, N = NA_integer_, cases = NA_integer_, noncases = NA_integer_, beta = NA_real_, 
        se = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_, p_raw = NA_real_, status = status, 
        warning = "", error = reason, formula = formula, rank = NA_integer_, columns = NA_integer_, sample_hash = NA_character_, 
        stringsAsFactors = FALSE)
}

bi_c_outcome_unavailable <- function(d, exposure) {
    if ("y" %in% names(d) && (!nrow(d) || any(!is.na(d$y)))) 
        return(NULL)
    bi_c_failure(exposure, "DATA_UNAVAILABLE", "REGISTERED_OUTCOME_COLUMN_MISSING_OR_ALL_UNKNOWN", paste(deparse(bi_formula(exposure)), 
        collapse = " "))
}

bi_c_job <- function(ctx, key, exposure, fn) {
    safe_key <- paste0(gsub("[^A-Za-z0-9_.-]", "_", key), "_", substr(bi_hash(key), 1L, 12L))
    bi_checkpoint(ctx, safe_key, function() {
        z <- tryCatch(fn(), error = function(e) bi_c_failure(exposure, "MODEL_ERROR", conditionMessage(e)))
        z$data_snapshot_id <- ctx$signature %||% NA_character_
        z$model_version <- ctx$cfg$version %||% "border_integrated_repair_20261016_v2"
        z$effect_unit <- unname(linear_units[exposure])
        z$variable_encoding <- "Frozen M0: linear age; factor sex/ethnicity/center/assessment_year; categorical S except border12; see explicit formula/group registry"
        z
    })
}

bi_c_annotate <- function(z, code, layer, level, test_id, family_id = NA_character_) {
    z$trait <- code
    z$analysis_layer <- layer
    z$stratum <- as.character(level)
    z$test_id <- test_id
    z$family_id <- family_id
    z$interpretation <- "Hospital ever-recorded disease OR; observational, not incident HR or causal risk"
    z
}

bi_c_complete <- function(z) {
    !nrow(z) || all(z$status %in% c("OK", "PASS", "SUCCESS", "NOT_APPLICABLE_FIXED_S"))
}

bi_c_derive <- function(d, definitions = NULL) {
    if (is.null(definitions)) 
        definitions <- list(age_group_breaks = c(-Inf, 50, 60, Inf), age_group_labels = c("<50", "50-59", ">=60"), 
            bmi_group_breaks = c(-Inf, 18.5, 25, 30, Inf), bmi_group_labels = c("<18.5", "18.5-<25", "25-<30", 
                ">=30"))
    d$S_group <- factor(ifelse(d$S %in% 2:4, as.character(d$S), ifelse(d$S >= 5 & d$S <= 8, "5-8", NA_character_)), 
        levels = c("2", "3", "4", "5-8"))
    d$birth_cohort <- cut(d$birth_year, breaks = c(1937, 1949, 1959, 1971), labels = c("1938-49", "1950-59", "1960-71"), 
        right = TRUE)
    d$age_group <- cut(d$age, breaks = as.numeric(unlist(definitions$age_group_breaks)), labels = as.character(unlist(definitions$age_group_labels)), 
        right = FALSE)
    if ("bmi" %in% names(d)) 
        d$bmi_group <- cut(d$bmi, breaks = as.numeric(unlist(definitions$bmi_group_breaks)), labels = as.character(unlist(definitions$bmi_group_labels)), 
            right = FALSE)
    d
}

bi_c_validate_settings <- function(ctx) {
    cfg <- bi_module_config(ctx, "clinical")
    definition <- ctx$cfg$clinical
    allowed <- c("age_group", "ethnicity", "center", "assessment_year", "bmi_group", "smoking", "alcohol", "family_composition", 
        "sibling_sex", "cnt_brother", "cnt_sister")
    if (is.null(cfg$secondary_modifiers)) 
        stop("LOCKED_SECONDARY_MODIFIER_REGISTRY_REQUIRED")
    modifiers <- as.character(unlist(cfg$secondary_modifiers))
    if (anyDuplicated(modifiers) || any(!modifiers %in% allowed)) 
        stop("INVALID_SECONDARY_MODIFIER_REGISTRY")
    fixed <- list(size_levels = 2:8, cumulative_maxima = 3:8, size_groups = c(2, 2, 3, 3, 4, 4, 5, 8), birth_cohorts = c(1938, 
        1949, 1950, 1959, 1960, 1971))
    for (k in names(fixed)) if (!identical(as.numeric(unlist(definition[[k]])), as.numeric(fixed[[k]]))) 
        stop("FROZEN_CLINICAL_GROUP_CONTRACT_MISMATCH: ", k)
    for (k in c("age", "bmi")) {
        breaks <- as.numeric(unlist(definition[[paste0(k, "_group_breaks")]]))
        labels <- as.character(unlist(definition[[paste0(k, "_group_labels")]]))
        if (length(breaks) < 3L || length(labels) != length(breaks) - 1L || anyNA(breaks) || any(diff(breaks) <= 
            0) || anyDuplicated(labels)) 
            stop("INVALID_LOCKED_EXPLORATORY_GROUPS: ", k)
    }
    modifiers
}

bi_c_interaction_registry <- function(ctx = NULL, modifiers = c("S_group", "sex", "birth_cohort")) {
    z <- expand.grid(trait = bi_c_codes(ctx), exposure = bi_c_exposures(), modifier = modifiers, stringsAsFactors = FALSE)
    z$applicable <- !(z$modifier == "S_group" & z$exposure == "border12")
    z$family_id <- ifelse(z$applicable, paste0("clinical_exploratory_", z$modifier), NA_character_)
    z$test_id <- paste("interaction", z$exposure, z$trait, z$modifier, sep = "__")
    z$inference_scope <- "Post-hoc exploratory modifier-specific BH; not formal whole-PheWAS q or independent replication"
    z
}

bi_c_fixed_strata <- function(d) {
    z <- list()
    for (s in 2:8) z[[paste0("exact_S_", s)]] <- list(layer = "exact_S", level = as.character(s), keep = !is.na(d$S) & 
        d$S == s)
    for (g in c("2", "3", "4", "5-8")) z[[paste0("group_S_", g)]] <- list(layer = "S_group", level = g, keep = !is.na(d$S_group) & 
        d$S_group == g)
    for (s in 3:8) z[[paste0("cumulative_S_", s)]] <- list(layer = "cumulative_S", level = paste0("S<=", s), keep = !is.na(d$S) & 
        d$S <= s)
    z[["cumulative_full"]] <- list(layer = "cumulative_S", level = "unrestricted", keep = rep(TRUE, nrow(d)))
    for (g in c("0", "1")) z[[paste0("sex_", g)]] <- list(layer = "sex", level = g, keep = !is.na(d$sex) & as.character(d$sex) == 
        g)
    for (g in c("1938-49", "1950-59", "1960-71")) z[[paste0("cohort_", g)]] <- list(layer = "birth_cohort", level = g, 
        keep = !is.na(d$birth_cohort) & d$birth_cohort == g)
    z
}

bi_c_modifier_levels <- function(d, modifier) {
    if (!modifier %in% names(d)) 
        return("unavailable")
    x <- d[[modifier]]
    if (is.factor(x)) 
        return(if (length(levels(x))) levels(x) else "unavailable")
    z <- sort(unique(as.character(x[!is.na(x)])))
    if (!length(z)) 
        "unavailable"
    else z
}

bi_c_saved_id_hash <- function(ids, format) {
    ids <- as.character(ids)
    if (anyNA(ids) || anyDuplicated(ids)) 
        stop("REUSE_INVALID_SAMPLE_IDS")
    if (format == "serialized_sorted_ids") 
        return(digest::digest(sort(ids), algo = "sha256"))
    txt <- if (length(ids)) 
        paste0(paste(sort(ids, method = "radix"), collapse = "\n"), "\n")
    else ""
    digest::digest(enc2utf8(txt), algo = "sha256", serialize = FALSE)
}

bi_c_reuse_source <- function(ctx, role) {
    if (is.null(ctx$clinical_reuse_sources)) 
        ctx$clinical_reuse_sources <- list()
    if (!is.null(ctx$clinical_reuse_sources[[role]])) 
        return(ctx$clinical_reuse_sources[[role]])
    source <- tryCatch({
        pinned <- c(all_rds = birth_private_setting("PRIVATE_HASH_1"), phewas_rds = birth_private_setting("PRIVATE_HASH_2"))
        if (!isTRUE(ctx$clinical_reuse_inputs_verified)) {
            for (key in names(pinned)) {
                if (!identical(ctx$cfg$input_sha256[[key]], unname(pinned[[key]]))) 
                  stop("REUSE_SOURCE_INPUT_PIN_MISMATCH: ", key)
                path <- ctx$cfg$inputs[[key]]
                if (is.null(path) || !file.exists(path) || !identical(bi_hash(path, TRUE), unname(pinned[[key]]))) 
                  stop("REUSE_CURRENT_INPUT_BYTES_UNVERIFIED: ", key)
            }
            ctx$clinical_reuse_inputs_verified <- TRUE
        }
        path <- ctx$cfg$inputs[[role]]
        expected <- ctx$cfg$input_sha256[[role]]
        if (is.null(path) || !file.exists(path) || is.null(expected)) 
            stop("REUSE_SOURCE_FILE_OR_HASH_MISSING")
        sha <- bi_hash(path, TRUE)
        if (!identical(sha, expected)) 
            stop("REUSE_SOURCE_FILE_HASH_MISMATCH")
        list(table = utils::read.csv(path, colClasses = "character", check.names = FALSE), path = path, sha256 = sha, 
            error = "")
    }, error = function(e) list(table = NULL, path = ctx$cfg$inputs[[role]] %||% NA_character_, sha256 = NA_character_, 
        error = conditionMessage(e)))
    ctx$clinical_reuse_sources[[role]] <- source
    source
}

bi_c_authenticate_saved <- function(d, exposure, code, role, source) {
    if (!is.data.frame(source$table)) 
        stop(source$error)
    table <- source$table
    if (!all(c("branch", "phecode", "actual_formula", "N", "cases", "noncases", "status") %in% names(table))) 
        stop("REUSE_SOURCE_SCHEMA_UNSUPPORTED")
    row <- table[as.character(table$branch) == exposure & as.character(table$phecode) == code, , drop = FALSE]
    if (nrow(row) != 1L) 
        stop("REUSE_SOURCE_KEY_NOT_UNIQUE_OR_MISSING")
    s <- bi_sample(d, bi_formula(exposure), exposure)
    for (k in c("N", "cases", "noncases")) if (!identical(as.numeric(s$counts[[k]]), as.numeric(row[[k]]))) 
        stop("REUSE_SAMPLE_COUNT_MISMATCH: ", k)
    normalize <- function(x) gsub("[[:space:]]+", "", x)
    formula <- paste(deparse(s$formula), collapse = " ")
    if (!identical(normalize(formula), normalize(row$actual_formula))) 
        stop("REUSE_FORMULA_MISMATCH")
    field <- if (role == "s8_results") 
        "sample_hash"
    else if (exposure == "border_relative") 
        "sample_hash_sha256"
    else "sample_id_digest"
    format <- if (field == "sample_id_digest") 
        "serialized_sorted_ids"
    else "sorted_newline_utf8"
    hash <- bi_c_saved_id_hash(s$data$eid, format)
    if (!field %in% names(row) || is.na(row[[field]]) || !identical(hash, as.character(row[[field]]))) 
        stop("REUSE_SAMPLE_HASH_MISMATCH: ", field)
    if (role == "s8_results") {
        x <- stats::model.matrix(s$formula, s$data)
        design <- digest::digest(list(x = x, y = s$data$y), algo = "sha256")
        if (!"design_hash" %in% names(row) || is.na(row$design_hash) || !identical(design, row$design_hash)) 
            stop("REUSE_S8_DESIGN_HASH_MISMATCH")
    }
    get <- function(name, default = NA_character_) if (name %in% names(row)) 
        row[[name]]
    else default
    num <- function(name) suppressWarnings(as.numeric(get(name)))
    result <- bi_c_failure(exposure, as.character(row$status), "", formula)
    for (k in c("N", "cases", "noncases")) result[[k]] <- as.numeric(s$counts[[k]])
    for (k in c("beta", "se", "OR", "CI_low", "CI_high", "p")) result[[k]] <- num(k)
    result$p_raw <- num("raw_target_p")
    result$raw_beta <- num("raw_target_beta")
    result$raw_se <- num("raw_target_se")
    result$rank <- num("matrix_rank")
    result$columns <- num("matrix_columns")
    result$warning <- get("warning", "")
    result$error <- get("error", get("reason", ""))
    result$sample_hash <- bi_hash(sort(as.character(s$data$eid), method = "radix"))
    result$removed_terms <- paste(s$removed_terms, collapse = ";")
    result$source_status <- get("status")
    result$source_inference <- get("inference_status", get("interpretation"))
    result$source_separation_assessment <- get("separation_assessment", get("lp_status"))
    result$inference <- "REUSED_SAVED_EVIDENCE_ORIGINAL_INFERENCE_LIMITATIONS_RETAINED"
    separation <- any(grepl("SEPARATION_REVIEW_REQUIRED|LP_NUISANCE_SEPARATION|LP_TARGET_SEPARATION|LP_COMPLETE_SEPARATION", 
        c(get("status"), get("lp_status"), get("separation_assessment"))), na.rm = TRUE)
    if (separation) {
        result$status <- "REUSED_SOURCE_REVIEW_REQUIRED"
        for (k in c("beta", "se", "OR", "CI_low", "CI_high", "p")) result[[k]] <- NA_real_
    }
    result$reuse_decision <- "REUSED_CERTIFIED_NO_FIT"
    result$reuse_reason <- "Source/input SHA256, model key, counts, formula and original sample hash match; S8 design hash also checked"
    result$reuse_source_role <- role
    result$reuse_source_path <- source$path
    result$reuse_source_sha256 <- source$sha256
    result$reuse_source_row_sha256 <- bi_hash(row)
    result$reuse_source_hash_field <- field
    result$reuse_source_sample_hash <- hash
    result
}

bi_c_fit_or_reuse <- function(ctx, d, exposure, code, spec) {
    unavailable <- bi_c_outcome_unavailable(d, exposure)
    if (!is.null(unavailable)) {
        unavailable$reuse_decision <- "OUTCOME_UNAVAILABLE_NO_FIT"
        return(unavailable)
    }
    candidate <- spec$layer == "cumulative_S" && spec$level %in% c("S<=8", "unrestricted")
    if (exposure == "border12" && spec$layer %in% c("exact_S", "S_group", "cumulative_S")) 
        candidate <- spec$layer == "cumulative_S" || spec$level == "2"
    if (!candidate) {
        z <- bi_fit(ctx, d, exposure)
        z$reuse_decision <- "NEW_ESTIMAND_REQUIRES_FIT"
        return(z)
    }
    role <- if (exposure != "border12" && spec$level == "S<=8") 
        "s8_results"
    else "main_results"
    source <- bi_c_reuse_source(ctx, role)
    tryCatch(bi_c_authenticate_saved(d, exposure, code, role, source), error = function(e) {
        z <- bi_c_failure(exposure, "REUSE_UNVERIFIED_NOT_REFITTED", conditionMessage(e), paste(deparse(bi_formula(exposure)), 
            collapse = " "))
        z$reuse_decision <- "NOT_CERTIFIED_NO_AUTOMATIC_REFIT"
        z$reuse_reason <- conditionMessage(e)
        z$reuse_source_role <- role
        z$reuse_source_path <- source$path
        z$reuse_source_sha256 <- source$sha256
        z
    })
}

bi_c_restrictions <- function(ctx, d) {
    fields <- bi_module_config(ctx, "clinical")$restriction_fields
    inspect <- function(key, field_id) {
        spec <- fields[[key]]
        column <- spec$column
        ok <- isTRUE(spec$verified) && identical(as.character(spec$field_id), field_id) && identical(as.integer(spec$instance), 
            0L) && identical(as.character(spec$coding), "100349") && length(column) == 1L && !is.na(column) && 
            column %in% names(d) && length(spec$provenance) == 1L && !is.na(spec$provenance) && nzchar(spec$provenance)
        if (!isTRUE(ok)) 
            return(list(available = FALSE, keep = rep(FALSE, nrow(d)), column = column %||% NA_character_, reason = paste0("MISSING_OR_UNVERIFIED_BASELINE_FIELD_", 
                field_id)))
        raw <- as.character(d[[column]])
        value <- suppressWarnings(as.numeric(raw))
        invalid <- !is.na(raw) & (is.na(value) | !value %in% c(0, 1, -1, -3))
        if (any(invalid)) 
            return(list(available = FALSE, keep = rep(FALSE, nrow(d)), column = column, reason = paste0("INVALID_CODING_100349_FIELD_", 
                field_id)))
        list(available = TRUE, keep = !is.na(value) & value == 0, column = column, reason = "Known No (0) only; Yes (1), unknown (-1), prefer-not (-3), and missing excluded")
    }
    a <- inspect("adopted", "1767")
    m <- inspect("multiple_birth", "1777")
    combined <- list(available = a$available && m$available, keep = a$keep & m$keep, column = paste(a$column, m$column, 
        sep = ";"), reason = paste(a$reason, m$reason, sep = " | "))
    list(known_not_adopted = a, known_not_multiple_birth = m, known_not_adopted_and_not_multiple_birth = combined)
}

bi_c_restriction_fit <- function(ctx, d, exposure, code, label, spec) {
    baseline <- bi_sample(d, bi_formula(exposure), exposure)
    unavailable <- bi_c_outcome_unavailable(d, exposure)
    if (!is.null(unavailable)) {
        z <- unavailable
        retained_N <- NA_integer_
    }
    else if (!spec$available) {
        z <- bi_c_failure(exposure, "SKIPPED_UNAVAILABLE", spec$reason)
        retained_N <- NA_integer_
    }
    else {
        restricted <- d[spec$keep, , drop = FALSE]
        ss <- bi_sample(restricted, bi_formula(exposure), exposure)
        unchanged <- identical(sort(as.character(baseline$data$eid)), sort(as.character(ss$data$eid)))
        z <- if (unchanged) 
            bi_c_fit_or_reuse(ctx, restricted, exposure, code, list(layer = "cumulative_S", level = "unrestricted"))
        else bi_fit(ctx, restricted, exposure)
        retained_N <- nrow(ss$data)
    }
    z$restriction <- label
    z$restriction_field <- spec$column
    z$restriction_definition <- spec$reason
    z$base_complete_N <- nrow(baseline$data)
    z$restricted_complete_N <- retained_N
    z$excluded_N <- nrow(baseline$data) - retained_N
    z$inference_scope <- "Post-hoc restriction sensitivity; original main eligibility unchanged; no causal claim"
    z
}

bi_stage_clinical <- function(ctx) {
    secondary <- bi_c_validate_settings(ctx)
    definition <- ctx$cfg$clinical
    groups <- data.frame(variable = c("S_exact", "S_group", "cumulative_S", "birth_cohort", "age_group", "bmi_group", 
        "secondary_modifiers"), definition = c("2;3;4;5;6;7;8", "2;3;4;5-8", "<=3;<=4;<=5;<=6;<=7;<=8;unrestricted", 
        "1938-49;1950-59;1960-71", paste(unlist(definition$age_group_breaks), collapse = ";"), paste(unlist(definition$bmi_group_breaks), 
            collapse = ";"), paste(secondary, collapse = ";")), encoding = c("Exact S", "S>8 excluded from S_group estimand", 
        "Nested not independent comparisons", "Inclusive year boundaries", "Left closed right open", "Left closed right open", 
        "Fixed configuration registry"))
    bi_save(ctx, "raw/clinical/group_definitions.csv", groups)
    estimates <- interactions <- flows <- restrictions <- list()
    k <- j <- h <- r <- 0L
    for (code in bi_c_codes(ctx)) {
        d <- bi_c_derive(bi_data(ctx, code), definition)
        h <- h + 1L
        flows[[h]] <- data.frame(trait = code, metric = c("cohort_rows", "birth_cohort_missing_or_outside_1938_1971", 
            "S_outside_2_8"), N = c(nrow(d), sum(is.na(d$birth_cohort)), sum(is.na(d$S_group))))
        for (exposure in bi_c_exposures()) {
            strata <- bi_c_fixed_strata(d)
            for (id in names(strata)) {
                spec <- strata[[id]]
                key <- paste(exposure, code, id, sep = "__")
                z <- bi_c_job(ctx, key, exposure, function() bi_c_fit_or_reuse(ctx, d[spec$keep, , drop = FALSE], 
                  exposure, code, spec))
                k <- k + 1L
                estimates[[k]] <- bi_c_annotate(z, code, spec$layer, spec$level, key)
            }
            for (modifier in c("S_group", "sex", "birth_cohort", secondary)) {
                key <- paste("interaction", exposure, code, modifier, sep = "__")
                family <- paste0("clinical_exploratory_", modifier)
                z <- bi_c_job(ctx, key, exposure, function() {
                  if (exposure == "border12" && modifier == "S_group") 
                    return(bi_c_failure(exposure, "NOT_APPLICABLE_FIXED_S", "border12 has fixed S=2; cross-S interaction is not an estimand"))
                  unavailable <- bi_c_outcome_unavailable(d, exposure)
                  if (!is.null(unavailable)) 
                    return(unavailable)
                  dd <- if (modifier == "sibling_sex") 
                    d[!is.na(d$S) & d$S == 2, , drop = FALSE]
                  else d
                  bi_interaction(ctx, dd, exposure, modifier)
                })
                if (exposure == "border12" && modifier == "S_group") 
                  family <- NA_character_
                z$modifier <- modifier
                z$inference_scope <- "Post-hoc exploratory modifier-specific BH; not formal whole-PheWAS q"
                j <- j + 1L
                interactions[[j]] <- bi_c_annotate(z, code, "global_interaction", modifier, key, family)
                if (modifier %in% secondary && !modifier %in% c("cnt_brother", "cnt_sister")) {
                  for (level in bi_c_modifier_levels(d, modifier)) {
                    skey <- paste("secondary_stratum", exposure, code, modifier, level, sep = "__")
                    zz <- bi_c_job(ctx, skey, exposure, function() {
                      if (!modifier %in% names(d)) 
                        return(bi_c_failure(exposure, "UNAVAILABLE_FIELD", modifier))
                      unavailable <- bi_c_outcome_unavailable(d, exposure)
                      if (!is.null(unavailable)) 
                        return(unavailable)
                      keep <- !is.na(d[[modifier]]) & as.character(d[[modifier]]) == level
                      if (modifier == "sibling_sex") 
                        keep <- keep & !is.na(d$S) & d$S == 2
                      bi_fit(ctx, d[keep, , drop = FALSE], exposure)
                    })
                    k <- k + 1L
                    estimates[[k]] <- bi_c_annotate(zz, code, paste0("exploratory_", modifier), level, skey)
                  }
                }
            }
            restriction_specs <- bi_c_restrictions(ctx, d)
            for (label in names(restriction_specs)) {
                key <- paste("restriction", exposure, code, label, sep = "__")
                z <- bi_c_job(ctx, key, exposure, function() bi_c_restriction_fit(ctx, d, exposure, code, label, 
                  restriction_specs[[label]]))
                r <- r + 1L
                restrictions[[r]] <- bi_c_annotate(z, code, "adoption_multiple_birth_restriction", label, key)
            }
        }
    }
    estimates <- bi_bind(estimates)
    interactions <- bi_bind(interactions)
    planned <- !is.na(interactions$family_id)
    adjusted <- bi_adjust(ctx, interactions[planned, , drop = FALSE])
    interactions <- bi_bind(list(adjusted, interactions[!planned, , drop = FALSE]))
    registry <- bi_c_interaction_registry(ctx, c("S_group", "sex", "birth_cohort", secondary))
    if (nrow(interactions) != nrow(registry) || !setequal(interactions$test_id, registry$test_id)) 
        stop("6_PHECODE_INTERACTION_REGISTRY_INVARIANT_FAILED")
    restrictions <- bi_bind(restrictions)
    bi_save(ctx, "raw/clinical/stratified_results.csv", estimates)
    reuse <- estimates[!is.na(estimates$reuse_decision) & estimates$reuse_decision != "NEW_ESTIMAND_REQUIRES_FIT", 
        , drop = FALSE]
    bi_save(ctx, "raw/clinical/existing_result_reuse.csv", reuse)
    bi_save(ctx, "raw/clinical/interactions.csv", interactions)
    bi_save(ctx, "raw/clinical/cohort_coverage.csv", bi_bind(flows))
    bi_save(ctx, "raw/clinical/exploratory_interaction_registry.csv", registry)
    bi_save(ctx, "raw/clinical/adoption_multiple_birth_sensitivity.csv", restrictions)
    list(status = if (bi_c_complete(estimates) && bi_c_complete(interactions) && bi_c_complete(restrictions)) "PASS" else "PARTIAL", 
        reason = "All six registered outcomes x four exposures (24 combinations); modifier-specific post-hoc exploratory families; unsupported/restriction-unavailable rows retained")
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        bi_c_exposures <- function() c("border12", "border1n", "border_rank", "border_relative")
        bi_c_codes <- function(ctx = NULL) {
            fixed <- c("411", "411.3", "411.4", "351", "250", "250.2")
            if (is.null(ctx)) return(fixed)
            configured <- as.character(unlist(ctx$cfg$clinical$core))
            if (!identical(configured, fixed)) stop("FIXED_6_PHECODE_REGISTRY_MISMATCH")
            configured
        }
        bi_c_failure <- function(exposure, status, reason, formula = NA_character_) {
            data.frame(exposure = exposure, N = NA_integer_, cases = NA_integer_, noncases = NA_integer_, beta = NA_real_, 
                se = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_, p_raw = NA_real_, 
                status = status, warning = "", error = reason, formula = formula, rank = NA_integer_, columns = NA_integer_, 
                sample_hash = NA_character_, stringsAsFactors = FALSE)
        }
        bi_c_outcome_unavailable <- function(d, exposure) {
            if ("y" %in% names(d) && (!nrow(d) || any(!is.na(d$y)))) return(NULL)
            bi_c_failure(exposure, "DATA_UNAVAILABLE", "REGISTERED_OUTCOME_COLUMN_MISSING_OR_ALL_UNKNOWN", paste(deparse(bi_formula(exposure)), 
                collapse = " "))
        }
        bi_c_job <- function(ctx, key, exposure, fn) {
            safe_key <- paste0(gsub("[^A-Za-z0-9_.-]", "_", key), "_", substr(bi_hash(key), 1L, 12L))
            bi_checkpoint(ctx, safe_key, function() {
                z <- tryCatch(fn(), error = function(e) bi_c_failure(exposure, "MODEL_ERROR", conditionMessage(e)))
                z$data_snapshot_id <- ctx$signature %||% NA_character_
                z$model_version <- ctx$cfg$version %||% "border_integrated_repair_20261016_v2"
                z$effect_unit <- unname(linear_units[exposure])
                z$variable_encoding <- "Frozen M0: linear age; factor sex/ethnicity/center/assessment_year; categorical S except border12; see explicit formula/group registry"
                z
            })
        }
        bi_c_annotate <- function(z, code, layer, level, test_id, family_id = NA_character_) {
            z$trait <- code
            z$analysis_layer <- layer
            z$stratum <- as.character(level)
            z$test_id <- test_id
            z$family_id <- family_id
            z$interpretation <- "Hospital ever-recorded disease OR; observational, not incident HR or causal risk"
            z
        }
        bi_c_complete <- function(z) {
            !nrow(z) || all(z$status %in% c("OK", "PASS", "SUCCESS", "NOT_APPLICABLE_FIXED_S"))
        }
        bi_c_derive <- function(d, definitions = NULL) {
            if (is.null(definitions)) definitions <- list(age_group_breaks = c(-Inf, 50, 60, Inf), age_group_labels = c("<50", 
                "50-59", ">=60"), bmi_group_breaks = c(-Inf, 18.5, 25, 30, Inf), bmi_group_labels = c("<18.5", 
                "18.5-<25", "25-<30", ">=30"))
            d$S_group <- factor(ifelse(d$S %in% 2:4, as.character(d$S), ifelse(d$S >= 5 & d$S <= 8, "5-8", NA_character_)), 
                levels = c("2", "3", "4", "5-8"))
            d$birth_cohort <- cut(d$birth_year, breaks = c(1937, 1949, 1959, 1971), labels = c("1938-49", "1950-59", 
                "1960-71"), right = TRUE)
            d$age_group <- cut(d$age, breaks = as.numeric(unlist(definitions$age_group_breaks)), labels = as.character(unlist(definitions$age_group_labels)), 
                right = FALSE)
            if ("bmi" %in% names(d)) d$bmi_group <- cut(d$bmi, breaks = as.numeric(unlist(definitions$bmi_group_breaks)), 
                labels = as.character(unlist(definitions$bmi_group_labels)), right = FALSE)
            d
        }
        bi_c_validate_settings <- function(ctx) {
            cfg <- bi_module_config(ctx, "clinical")
            definition <- ctx$cfg$clinical
            allowed <- c("age_group", "ethnicity", "center", "assessment_year", "bmi_group", "smoking", "alcohol", 
                "family_composition", "sibling_sex", "cnt_brother", "cnt_sister")
            if (is.null(cfg$secondary_modifiers)) stop("LOCKED_SECONDARY_MODIFIER_REGISTRY_REQUIRED")
            modifiers <- as.character(unlist(cfg$secondary_modifiers))
            if (anyDuplicated(modifiers) || any(!modifiers %in% allowed)) stop("INVALID_SECONDARY_MODIFIER_REGISTRY")
            fixed <- list(size_levels = 2:8, cumulative_maxima = 3:8, size_groups = c(2, 2, 3, 3, 4, 4, 5, 8), 
                birth_cohorts = c(1938, 1949, 1950, 1959, 1960, 1971))
            for (k in names(fixed)) if (!identical(as.numeric(unlist(definition[[k]])), as.numeric(fixed[[k]]))) stop("FROZEN_CLINICAL_GROUP_CONTRACT_MISMATCH: ", 
                k)
            for (k in c("age", "bmi")) {
                breaks <- as.numeric(unlist(definition[[paste0(k, "_group_breaks")]]))
                labels <- as.character(unlist(definition[[paste0(k, "_group_labels")]]))
                if (length(breaks) < 3L || length(labels) != length(breaks) - 1L || anyNA(breaks) || any(diff(breaks) <= 
                  0) || anyDuplicated(labels)) stop("INVALID_LOCKED_EXPLORATORY_GROUPS: ", k)
            }
            modifiers
        }
        bi_c_interaction_registry <- function(ctx = NULL, modifiers = c("S_group", "sex", "birth_cohort")) {
            z <- expand.grid(trait = bi_c_codes(ctx), exposure = bi_c_exposures(), modifier = modifiers, stringsAsFactors = FALSE)
            z$applicable <- !(z$modifier == "S_group" & z$exposure == "border12")
            z$family_id <- ifelse(z$applicable, paste0("clinical_exploratory_", z$modifier), NA_character_)
            z$test_id <- paste("interaction", z$exposure, z$trait, z$modifier, sep = "__")
            z$inference_scope <- "Post-hoc exploratory modifier-specific BH; not formal whole-PheWAS q or independent replication"
            z
        }
        bi_c_fixed_strata <- function(d) {
            z <- list()
            for (s in 2:8) z[[paste0("exact_S_", s)]] <- list(layer = "exact_S", level = as.character(s), keep = !is.na(d$S) & 
                d$S == s)
            for (g in c("2", "3", "4", "5-8")) z[[paste0("group_S_", g)]] <- list(layer = "S_group", level = g, 
                keep = !is.na(d$S_group) & d$S_group == g)
            for (s in 3:8) z[[paste0("cumulative_S_", s)]] <- list(layer = "cumulative_S", level = paste0("S<=", 
                s), keep = !is.na(d$S) & d$S <= s)
            z[["cumulative_full"]] <- list(layer = "cumulative_S", level = "unrestricted", keep = rep(TRUE, nrow(d)))
            for (g in c("0", "1")) z[[paste0("sex_", g)]] <- list(layer = "sex", level = g, keep = !is.na(d$sex) & 
                as.character(d$sex) == g)
            for (g in c("1938-49", "1950-59", "1960-71")) z[[paste0("cohort_", g)]] <- list(layer = "birth_cohort", 
                level = g, keep = !is.na(d$birth_cohort) & d$birth_cohort == g)
            z
        }
        bi_c_modifier_levels <- function(d, modifier) {
            if (!modifier %in% names(d)) return("unavailable")
            x <- d[[modifier]]
            if (is.factor(x)) return(if (length(levels(x))) levels(x) else "unavailable")
            z <- sort(unique(as.character(x[!is.na(x)])))
            if (!length(z)) "unavailable" else z
        }
        bi_c_saved_id_hash <- function(ids, format) {
            ids <- as.character(ids)
            if (anyNA(ids) || anyDuplicated(ids)) stop("REUSE_INVALID_SAMPLE_IDS")
            if (format == "serialized_sorted_ids") return(digest::digest(sort(ids), algo = "sha256"))
            txt <- if (length(ids)) paste0(paste(sort(ids, method = "radix"), collapse = "\n"), "\n") else ""
            digest::digest(enc2utf8(txt), algo = "sha256", serialize = FALSE)
        }
        bi_c_reuse_source <- function(ctx, role) {
            if (is.null(ctx$clinical_reuse_sources)) ctx$clinical_reuse_sources <- list()
            if (!is.null(ctx$clinical_reuse_sources[[role]])) return(ctx$clinical_reuse_sources[[role]])
            source <- tryCatch({
                pinned <- c(all_rds = birth_private_setting("PRIVATE_HASH_3"), phewas_rds = birth_private_setting("PRIVATE_HASH_4"))
                if (!isTRUE(ctx$clinical_reuse_inputs_verified)) {
                  for (key in names(pinned)) {
                    if (!identical(ctx$cfg$input_sha256[[key]], unname(pinned[[key]]))) stop("REUSE_SOURCE_INPUT_PIN_MISMATCH: ", 
                      key)
                    path <- ctx$cfg$inputs[[key]]
                    if (is.null(path) || !file.exists(path) || !identical(bi_hash(path, TRUE), unname(pinned[[key]]))) stop("REUSE_CURRENT_INPUT_BYTES_UNVERIFIED: ", 
                      key)
                  }
                  ctx$clinical_reuse_inputs_verified <- TRUE
                }
                path <- ctx$cfg$inputs[[role]]
                expected <- ctx$cfg$input_sha256[[role]]
                if (is.null(path) || !file.exists(path) || is.null(expected)) stop("REUSE_SOURCE_FILE_OR_HASH_MISSING")
                sha <- bi_hash(path, TRUE)
                if (!identical(sha, expected)) stop("REUSE_SOURCE_FILE_HASH_MISMATCH")
                list(table = utils::read.csv(path, colClasses = "character", check.names = FALSE), path = path, 
                  sha256 = sha, error = "")
            }, error = function(e) list(table = NULL, path = ctx$cfg$inputs[[role]] %||% NA_character_, sha256 = NA_character_, 
                error = conditionMessage(e)))
            ctx$clinical_reuse_sources[[role]] <- source
            source
        }
        bi_c_authenticate_saved <- function(d, exposure, code, role, source) {
            if (!is.data.frame(source$table)) stop(source$error)
            table <- source$table
            if (!all(c("branch", "phecode", "actual_formula", "N", "cases", "noncases", "status") %in% names(table))) stop("REUSE_SOURCE_SCHEMA_UNSUPPORTED")
            row <- table[as.character(table$branch) == exposure & as.character(table$phecode) == code, , drop = FALSE]
            if (nrow(row) != 1L) stop("REUSE_SOURCE_KEY_NOT_UNIQUE_OR_MISSING")
            s <- bi_sample(d, bi_formula(exposure), exposure)
            for (k in c("N", "cases", "noncases")) if (!identical(as.numeric(s$counts[[k]]), as.numeric(row[[k]]))) stop("REUSE_SAMPLE_COUNT_MISMATCH: ", 
                k)
            normalize <- function(x) gsub("[[:space:]]+", "", x)
            formula <- paste(deparse(s$formula), collapse = " ")
            if (!identical(normalize(formula), normalize(row$actual_formula))) stop("REUSE_FORMULA_MISMATCH")
            field <- if (role == "s8_results") "sample_hash" else if (exposure == "border_relative") "sample_hash_sha256" else "sample_id_digest"
            format <- if (field == "sample_id_digest") "serialized_sorted_ids" else "sorted_newline_utf8"
            hash <- bi_c_saved_id_hash(s$data$eid, format)
            if (!field %in% names(row) || is.na(row[[field]]) || !identical(hash, as.character(row[[field]]))) stop("REUSE_SAMPLE_HASH_MISMATCH: ", 
                field)
            if (role == "s8_results") {
                x <- stats::model.matrix(s$formula, s$data)
                design <- digest::digest(list(x = x, y = s$data$y), algo = "sha256")
                if (!"design_hash" %in% names(row) || is.na(row$design_hash) || !identical(design, row$design_hash)) stop("REUSE_S8_DESIGN_HASH_MISMATCH")
            }
            get <- function(name, default = NA_character_) if (name %in% names(row)) row[[name]] else default
            num <- function(name) suppressWarnings(as.numeric(get(name)))
            result <- bi_c_failure(exposure, as.character(row$status), "", formula)
            for (k in c("N", "cases", "noncases")) result[[k]] <- as.numeric(s$counts[[k]])
            for (k in c("beta", "se", "OR", "CI_low", "CI_high", "p")) result[[k]] <- num(k)
            result$p_raw <- num("raw_target_p")
            result$raw_beta <- num("raw_target_beta")
            result$raw_se <- num("raw_target_se")
            result$rank <- num("matrix_rank")
            result$columns <- num("matrix_columns")
            result$warning <- get("warning", "")
            result$error <- get("error", get("reason", ""))
            result$sample_hash <- bi_hash(sort(as.character(s$data$eid), method = "radix"))
            result$removed_terms <- paste(s$removed_terms, collapse = ";")
            result$source_status <- get("status")
            result$source_inference <- get("inference_status", get("interpretation"))
            result$source_separation_assessment <- get("separation_assessment", get("lp_status"))
            result$inference <- "REUSED_SAVED_EVIDENCE_ORIGINAL_INFERENCE_LIMITATIONS_RETAINED"
            separation <- any(grepl("SEPARATION_REVIEW_REQUIRED|LP_NUISANCE_SEPARATION|LP_TARGET_SEPARATION|LP_COMPLETE_SEPARATION", 
                c(get("status"), get("lp_status"), get("separation_assessment"))), na.rm = TRUE)
            if (separation) {
                result$status <- "REUSED_SOURCE_REVIEW_REQUIRED"
                for (k in c("beta", "se", "OR", "CI_low", "CI_high", "p")) result[[k]] <- NA_real_
            }
            result$reuse_decision <- "REUSED_CERTIFIED_NO_FIT"
            result$reuse_reason <- "Source/input SHA256, model key, counts, formula and original sample hash match; S8 design hash also checked"
            result$reuse_source_role <- role
            result$reuse_source_path <- source$path
            result$reuse_source_sha256 <- source$sha256
            result$reuse_source_row_sha256 <- bi_hash(row)
            result$reuse_source_hash_field <- field
            result$reuse_source_sample_hash <- hash
            result
        }
        bi_c_fit_or_reuse <- function(ctx, d, exposure, code, spec) {
            unavailable <- bi_c_outcome_unavailable(d, exposure)
            if (!is.null(unavailable)) {
                unavailable$reuse_decision <- "OUTCOME_UNAVAILABLE_NO_FIT"
                return(unavailable)
            }
            candidate <- spec$layer == "cumulative_S" && spec$level %in% c("S<=8", "unrestricted")
            if (exposure == "border12" && spec$layer %in% c("exact_S", "S_group", "cumulative_S")) candidate <- spec$layer == 
                "cumulative_S" || spec$level == "2"
            if (!candidate) {
                z <- bi_fit(ctx, d, exposure)
                z$reuse_decision <- "NEW_ESTIMAND_REQUIRES_FIT"
                return(z)
            }
            role <- if (exposure != "border12" && spec$level == "S<=8") "s8_results" else "main_results"
            source <- bi_c_reuse_source(ctx, role)
            tryCatch(bi_c_authenticate_saved(d, exposure, code, role, source), error = function(e) {
                z <- bi_c_failure(exposure, "REUSE_UNVERIFIED_NOT_REFITTED", conditionMessage(e), paste(deparse(bi_formula(exposure)), 
                  collapse = " "))
                z$reuse_decision <- "NOT_CERTIFIED_NO_AUTOMATIC_REFIT"
                z$reuse_reason <- conditionMessage(e)
                z$reuse_source_role <- role
                z$reuse_source_path <- source$path
                z$reuse_source_sha256 <- source$sha256
                z
            })
        }
        bi_c_restrictions <- function(ctx, d) {
            fields <- bi_module_config(ctx, "clinical")$restriction_fields
            inspect <- function(key, field_id) {
                spec <- fields[[key]]
                column <- spec$column
                ok <- isTRUE(spec$verified) && identical(as.character(spec$field_id), field_id) && identical(as.integer(spec$instance), 
                  0L) && identical(as.character(spec$coding), "100349") && length(column) == 1L && !is.na(column) && 
                  column %in% names(d) && length(spec$provenance) == 1L && !is.na(spec$provenance) && nzchar(spec$provenance)
                if (!isTRUE(ok)) return(list(available = FALSE, keep = rep(FALSE, nrow(d)), column = column %||% 
                  NA_character_, reason = paste0("MISSING_OR_UNVERIFIED_BASELINE_FIELD_", field_id)))
                raw <- as.character(d[[column]])
                value <- suppressWarnings(as.numeric(raw))
                invalid <- !is.na(raw) & (is.na(value) | !value %in% c(0, 1, -1, -3))
                if (any(invalid)) return(list(available = FALSE, keep = rep(FALSE, nrow(d)), column = column, reason = paste0("INVALID_CODING_100349_FIELD_", 
                  field_id)))
                list(available = TRUE, keep = !is.na(value) & value == 0, column = column, reason = "Known No (0) only; Yes (1), unknown (-1), prefer-not (-3), and missing excluded")
            }
            a <- inspect("adopted", "1767")
            m <- inspect("multiple_birth", "1777")
            combined <- list(available = a$available && m$available, keep = a$keep & m$keep, column = paste(a$column, 
                m$column, sep = ";"), reason = paste(a$reason, m$reason, sep = " | "))
            list(known_not_adopted = a, known_not_multiple_birth = m, known_not_adopted_and_not_multiple_birth = combined)
        }
        bi_c_restriction_fit <- function(ctx, d, exposure, code, label, spec) {
            baseline <- bi_sample(d, bi_formula(exposure), exposure)
            unavailable <- bi_c_outcome_unavailable(d, exposure)
            if (!is.null(unavailable)) {
                z <- unavailable
                retained_N <- NA_integer_
            } else if (!spec$available) {
                z <- bi_c_failure(exposure, "SKIPPED_UNAVAILABLE", spec$reason)
                retained_N <- NA_integer_
            } else {
                restricted <- d[spec$keep, , drop = FALSE]
                ss <- bi_sample(restricted, bi_formula(exposure), exposure)
                unchanged <- identical(sort(as.character(baseline$data$eid)), sort(as.character(ss$data$eid)))
                z <- if (unchanged) bi_c_fit_or_reuse(ctx, restricted, exposure, code, list(layer = "cumulative_S", 
                  level = "unrestricted")) else bi_fit(ctx, restricted, exposure)
                retained_N <- nrow(ss$data)
            }
            z$restriction <- label
            z$restriction_field <- spec$column
            z$restriction_definition <- spec$reason
            z$base_complete_N <- nrow(baseline$data)
            z$restricted_complete_N <- retained_N
            z$excluded_N <- nrow(baseline$data) - retained_N
            z$inference_scope <- "Post-hoc restriction sensitivity; original main eligibility unchanged; no causal claim"
            z
        }
        bi_stage_clinical <- function(ctx) {
            secondary <- bi_c_validate_settings(ctx)
            definition <- ctx$cfg$clinical
            groups <- data.frame(variable = c("S_exact", "S_group", "cumulative_S", "birth_cohort", "age_group", 
                "bmi_group", "secondary_modifiers"), definition = c("2;3;4;5;6;7;8", "2;3;4;5-8", "<=3;<=4;<=5;<=6;<=7;<=8;unrestricted", 
                "1938-49;1950-59;1960-71", paste(unlist(definition$age_group_breaks), collapse = ";"), paste(unlist(definition$bmi_group_breaks), 
                  collapse = ";"), paste(secondary, collapse = ";")), encoding = c("Exact S", "S>8 excluded from S_group estimand", 
                "Nested not independent comparisons", "Inclusive year boundaries", "Left closed right open", "Left closed right open", 
                "Fixed configuration registry"))
            bi_save(ctx, "raw/clinical/group_definitions.csv", groups)
            estimates <- interactions <- flows <- restrictions <- list()
            k <- j <- h <- r <- 0L
            for (code in bi_c_codes(ctx)) {
                d <- bi_c_derive(bi_data(ctx, code), definition)
                h <- h + 1L
                flows[[h]] <- data.frame(trait = code, metric = c("cohort_rows", "birth_cohort_missing_or_outside_1938_1971", 
                  "S_outside_2_8"), N = c(nrow(d), sum(is.na(d$birth_cohort)), sum(is.na(d$S_group))))
                for (exposure in bi_c_exposures()) {
                  strata <- bi_c_fixed_strata(d)
                  for (id in names(strata)) {
                    spec <- strata[[id]]
                    key <- paste(exposure, code, id, sep = "__")
                    z <- bi_c_job(ctx, key, exposure, function() bi_c_fit_or_reuse(ctx, d[spec$keep, , drop = FALSE], 
                      exposure, code, spec))
                    k <- k + 1L
                    estimates[[k]] <- bi_c_annotate(z, code, spec$layer, spec$level, key)
                  }
                  for (modifier in c("S_group", "sex", "birth_cohort", secondary)) {
                    key <- paste("interaction", exposure, code, modifier, sep = "__")
                    family <- paste0("clinical_exploratory_", modifier)
                    z <- bi_c_job(ctx, key, exposure, function() {
                      if (exposure == "border12" && modifier == "S_group") return(bi_c_failure(exposure, "NOT_APPLICABLE_FIXED_S", 
                        "border12 has fixed S=2; cross-S interaction is not an estimand"))
                      unavailable <- bi_c_outcome_unavailable(d, exposure)
                      if (!is.null(unavailable)) return(unavailable)
                      dd <- if (modifier == "sibling_sex") d[!is.na(d$S) & d$S == 2, , drop = FALSE] else d
                      bi_interaction(ctx, dd, exposure, modifier)
                    })
                    if (exposure == "border12" && modifier == "S_group") family <- NA_character_
                    z$modifier <- modifier
                    z$inference_scope <- "Post-hoc exploratory modifier-specific BH; not formal whole-PheWAS q"
                    j <- j + 1L
                    interactions[[j]] <- bi_c_annotate(z, code, "global_interaction", modifier, key, family)
                    if (modifier %in% secondary && !modifier %in% c("cnt_brother", "cnt_sister")) {
                      for (level in bi_c_modifier_levels(d, modifier)) {
                        skey <- paste("secondary_stratum", exposure, code, modifier, level, sep = "__")
                        zz <- bi_c_job(ctx, skey, exposure, function() {
                          if (!modifier %in% names(d)) return(bi_c_failure(exposure, "UNAVAILABLE_FIELD", modifier))
                          unavailable <- bi_c_outcome_unavailable(d, exposure)
                          if (!is.null(unavailable)) return(unavailable)
                          keep <- !is.na(d[[modifier]]) & as.character(d[[modifier]]) == level
                          if (modifier == "sibling_sex") keep <- keep & !is.na(d$S) & d$S == 2
                          bi_fit(ctx, d[keep, , drop = FALSE], exposure)
                        })
                        k <- k + 1L
                        estimates[[k]] <- bi_c_annotate(zz, code, paste0("exploratory_", modifier), level, skey)
                      }
                    }
                  }
                  restriction_specs <- bi_c_restrictions(ctx, d)
                  for (label in names(restriction_specs)) {
                    key <- paste("restriction", exposure, code, label, sep = "__")
                    z <- bi_c_job(ctx, key, exposure, function() bi_c_restriction_fit(ctx, d, exposure, code, label, 
                      restriction_specs[[label]]))
                    r <- r + 1L
                    restrictions[[r]] <- bi_c_annotate(z, code, "adoption_multiple_birth_restriction", label, key)
                  }
                }
            }
            estimates <- bi_bind(estimates)
            interactions <- bi_bind(interactions)
            planned <- !is.na(interactions$family_id)
            adjusted <- bi_adjust(ctx, interactions[planned, , drop = FALSE])
            interactions <- bi_bind(list(adjusted, interactions[!planned, , drop = FALSE]))
            registry <- bi_c_interaction_registry(ctx, c("S_group", "sex", "birth_cohort", secondary))
            if (nrow(interactions) != nrow(registry) || !setequal(interactions$test_id, registry$test_id)) stop("6_PHECODE_INTERACTION_REGISTRY_INVARIANT_FAILED")
            restrictions <- bi_bind(restrictions)
            bi_save(ctx, "raw/clinical/stratified_results.csv", estimates)
            reuse <- estimates[!is.na(estimates$reuse_decision) & estimates$reuse_decision != "NEW_ESTIMAND_REQUIRES_FIT", 
                , drop = FALSE]
            bi_save(ctx, "raw/clinical/existing_result_reuse.csv", reuse)
            bi_save(ctx, "raw/clinical/interactions.csv", interactions)
            bi_save(ctx, "raw/clinical/cohort_coverage.csv", bi_bind(flows))
            bi_save(ctx, "raw/clinical/exploratory_interaction_registry.csv", registry)
            bi_save(ctx, "raw/clinical/adoption_multiple_birth_sensitivity.csv", restrictions)
            list(status = if (bi_c_complete(estimates) && bi_c_complete(interactions) && bi_c_complete(restrictions)) "PASS" else "PARTIAL", 
                reason = "All six registered outcomes x four exposures (24 combinations); modifier-specific post-hoc exploratory families; unsupported/restriction-unavailable rows retained")
        }
    }), envir = context)
}

