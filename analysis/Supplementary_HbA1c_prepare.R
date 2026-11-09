# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

gs_assert <- function(ok, reason) if (!isTRUE(ok)) stop(reason, call. = FALSE)

gs_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)

gs_sha <- function(path) digest::digest(file = path, algo = "sha256")

gs_vec <- function(x) as.character(unlist(x, use.names = FALSE))

gs_utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

gs_new_json <- function(x, path) {
    gs_assert(!file.exists(path), "PRESERVE_EXISTING_JSON")
    tmp <- paste0(path, ".pending_", Sys.getpid())
    gs_assert(!file.exists(tmp), "PRESERVE_EXISTING_PENDING_JSON")
    jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
    gs_assert(!file.exists(path) && file.rename(tmp, path), "ADDITIVE_JSON_PUBLICATION_FAILED")
    invisible(list(path = path, sha256 = gs_sha(path)))
}

gs_new_rds <- function(x, path) {
    gs_assert(!file.exists(path), "PRESERVE_EXISTING_PRIVATE_RDS")
    tmp <- paste0(path, ".pending_", Sys.getpid())
    gs_assert(!file.exists(tmp), "PRESERVE_EXISTING_PENDING_RDS")
    saveRDS(x, tmp, version = 2L, compress = FALSE)
    gs_assert(!file.exists(path) && file.rename(tmp, path), "ADDITIVE_RDS_PUBLICATION_FAILED")
    invisible(list(path = path, sha256 = gs_sha(path), bytes = unname(file.info(path)$size)))
}

gs_norm <- function(path, mustWork = TRUE) {
    gs_assert(is.character(path) && length(path) == 1L && !is.na(path) && nzchar(path) && !grepl("(^|/)[.][.](/|$)", 
        path) && startsWith(path, "/"), "ABSOLUTE_EXACT_LINUX_PATH_REQUIRED")
    normalizePath(path, winslash = "/", mustWork = mustWork)
}

gs_inside <- function(path, root) identical(path, root) || startsWith(path, paste0(root, "/"))

gs_check_config <- function(cfg) {
    gs_assert(identical(cfg$scope, "BASE11_B0_ONLY") && identical(gs_vec(cfg$bootstrap_ids), "0") && length(cfg$b) == 
        1L && is.numeric(cfg$b) && !is.na(cfg$b) && cfg$b == 0, "B0_ONLY_NO_BOOTSTRAP_DISPATCH")
    gs_assert(isTRUE(cfg$authorization$real_B0) && identical(cfg$authorization$bootstrap, FALSE), "EXPLICIT_REAL_B0_ONLY_AUTHORIZATION_REQUIRED")
    gs_assert(length(cfg$seed) == 1L && is.numeric(cfg$seed) && !is.na(cfg$seed) && cfg$seed == 20261103L, "EXACT_ORIGINAL_BASE_B0_SEED_20261103_REQUIRED")
    for (key in c("n_mc", "uniform_dim", "block_size", "eval_year")) {
        expected <- switch(key, n_mc = 50L, uniform_dim = 11L, block_size = 512L, eval_year = 11L)
        if (!is.null(cfg[[key]])) 
            gs_assert(identical(as.numeric(cfg[[key]]), as.numeric(expected)), "FROZEN_INTEGRATION_CONTRACT_CHANGED")
    }
    cfg
}

gs_resolve_config <- function(cfg) {
    cfg <- gs_check_config(cfg)
    keys <- c("input_root", "clinical_path", "source_root", "y_root", "state_root", "output_root", "adapter_root")
    for (key in keys) {
        override <- Sys.getenv(paste0("GS44_", toupper(key)))
        if (nzchar(override)) 
            cfg[[key]] <- override
        gs_assert(!is.null(cfg[[key]]), paste0("MISSING_CONFIG_ROOT_", toupper(key)))
        cfg[[key]] <- gs_norm(cfg[[key]], mustWork = key != "state_root")
    }
    if (!file.exists(cfg$state_root)) 
        cfg$state_root <- file.path(gs_norm(dirname(cfg$state_root)), basename(cfg$state_root))
    namespace <- "border_glycaemia_sensitivity_BASE11_B0_20261106_v1"
    gs_assert(grepl(paste0("/", namespace, "/"), cfg$output_root, fixed = TRUE) && basename(dirname(cfg$state_root)) == 
        namespace && basename(cfg$state_root) == "state", "ISOLATED_NEW_TASK_OUTPUT_AND_STATE_REQUIRED")
    reads <- unlist(cfg[c("input_root", "clinical_path", "source_root", "y_root", "adapter_root")])
    gs_assert(!any(vapply(reads, function(p) gs_inside(cfg$output_root, p) || gs_inside(cfg$state_root, p) || gs_inside(p, 
        cfg$output_root) || gs_inside(p, cfg$state_root), FALSE)), "NEW_OUTPUTS_MUST_NOT_OVERLAP_READONLY_INPUT_OR_SOURCE_ROOTS")
    gs_assert(!gs_inside(cfg$output_root, cfg$state_root) && !gs_inside(cfg$state_root, cfg$output_root), "SEPARATE_NEW_STATE_AND_OUTPUT_TREES_REQUIRED")
    cfg
}

gs_input_pins <- function() list(input = birth_private_setting("PRIVATE_HASH_1"), bootstrap = birth_private_setting("PRIVATE_HASH_2"), 
    weights = birth_private_setting("PRIVATE_HASH_3"), clinical = birth_private_setting("PRIVATE_HASH_4"))

gs_endpoint_counts <- function(r) list(N = nrow(r), events_full = sum(r$event == 1), events_11y = sum(r$event == 
    1 & r$exit <= 11), deaths_full = sum(r$death_event == 1), deaths_11y = sum(r$death_event == 1 & r$exit <= 11), 
    events_after11 = sum(r$event == 1 & r$exit > 11), followup_max = max(r$exit))

gs_raw_hba1c_mask <- function(value, reportability) {
    gs_assert(length(value) == length(reportability), "RAW_HBA1C_QUALITY_LENGTH_MISMATCH")
    h <- suppressWarnings(as.numeric(value))
    q <- suppressWarnings(as.numeric(reportability))
    available <- is.finite(h) & h > 0 & !is.na(q) & q == 1
    exclude <- available & h >= 48
    list(value = h, available = available, exclude = exclude, retain = !exclude)
}

gs_prepare_B0 <- function(cfg, private, receipts, progress = function(...) NULL) {
    gs_check_config(cfg)
    paths <- list(input = file.path(cfg$input_root, "LINE_INPUTS.private.rds"), bootstrap = file.path(cfg$input_root, 
        "draws", "bootstrap_0000.private.rds"), weights = file.path(cfg$input_root, "draws", "weights_0000.private.rds"), 
        clinical = cfg$clinical_path)
    pins <- gs_input_pins()
    for (key in names(paths)) {
        gs_assert(file.exists(paths[[key]]), paste0("MISSING_INPUT_PATH_", toupper(key)))
        if (!is.null(cfg[[paste0(key, "_sha256")]])) 
            gs_assert(identical(cfg[[paste0(key, "_sha256")]], pins[[key]]), "IMMUTABLE_INPUT_PIN_OVERRIDE_REJECTED")
        gs_assert(identical(gs_sha(paths[[key]]), pins[[key]]), paste0("FROZEN_INPUT_HASH_CHANGED_", toupper(key)))
    }
    progress("INPUT_PINS_VERIFIED")
    original <- readRDS(paths$input)
    old_index <- readRDS(paths$bootstrap)
    old_ws <- readRDS(paths$weights)
    ids <- c("smoking", "bmi", "weight", "fat_free_mass", "fat_mass", "hip", "waist", "clinical__hdl", "clinical__tg", 
        "clinical__glucose", "clinical__hba1c")
    base <- c("B", "age", "sex", "ethnicity", "center", "assessment_year", "S")
    order <- c("SMOKING", "BODY_COMPOSITION", "CLINICAL_LIPIDS", "GLYCAEMIA")
    cols <- vapply(original$members, `[[`, "", "column")
    gs_assert(identical(original$schema, "SELECTED44_LINE_FROZEN_INPUTS_V1") && identical(original$line, "BASE") && 
        nrow(original$data) == 123350L && identical(names(original$members), ids) && identical(gs_vec(original$base), 
        base) && original$uniform_dim == 11L && identical(names(original$domains), order) && identical(gs_vec(original$model_specifications$S44_BASE11$allocation_order), 
        order) && identical(gs_vec(original$model_specifications$S44_BASE11$members), ids) && !isTRUE(original$model_specifications$S44_BASE11$total_only), 
        "FROZEN_BASE11_SCIENTIFIC_CONTRACT_CHANGED")
    gs_assert(!anyDuplicated(original$data$eid) && !anyNA(original$data$cluster_id) && identical(gs_hash(as.character(original$data$eid)), 
        original$target_hash) && identical(original$target_hash, birth_private_setting("PRIVATE_HASH_5")) && identical(gs_hash(original$selection_functions), 
        original$selection_source_hash), "ORIGINAL_TARGET_OR_SELECTION_SOURCE_BINDING_CHANGED")
    gs_assert(old_index$bootstrap_id == 0L && identical(as.integer(old_index$indices), seq_len(nrow(original$data))) && 
        identical(old_index$target_hash, original$target_hash) && identical(old_index$cluster_membership_hash, 
        original$cluster_membership_hash) && is.logical(old_ws$complete) && !anyNA(old_ws$complete) && identical(old_ws$complete, 
        complete.cases(original$data[, unique(c(base, cols)), drop = FALSE])) && sum(old_ws$complete) == 95314L && 
        identical(gs_hash(as.character(original$data$eid[old_ws$complete])), birth_private_setting("PRIVATE_HASH_6")), 
        "ORIGINAL_B0_IDENTITY_INDEX_AND_COMPLETE_MASK_CHANGED")
    gs_assert(identical(old_ws$pair_binding$target_hash, original$target_hash) && identical(old_ws$pair_binding$actual_index_sha256, 
        pins$bootstrap) && identical(gs_hash(as.numeric(old_ws$w)), old_ws$pair_binding$selection_weight_hash) && 
        identical(gs_hash(as.integer(old_index$indices)), old_ws$pair_binding$sampling_hash), "ORIGINAL_B0_WEIGHT_INDEX_PROVENANCE_CHANGED")
    header <- names(utils::read.csv(paths$clinical, nrows = 0L, check.names = FALSE))
    raw_cols <- c("eid", "p30750_i0", "p30756_i0")
    gs_assert(all(raw_cols %in% header), "RAW_BASELINE_HBA1C_AND_REPORTABILITY_COLUMNS_MISSING")
    classes <- rep("NULL", length(header))
    classes[match(raw_cols, header)] <- "character"
    raw <- utils::read.csv(paths$clinical, colClasses = classes, check.names = FALSE)
    gs_assert(!anyDuplicated(raw$eid), "DUPLICATE_RAW_CLINICAL_IDENTITY")
    j <- match(as.character(original$data$eid), raw$eid)
    gs_assert(!anyNA(j), "RAW_CLINICAL_TARGET_ALIGNMENT_FAILED")
    mask <- gs_raw_hba1c_mask(raw$p30750_i0[j], raw$p30756_i0[j])
    h <- mask$value
    available <- mask$available
    exclude <- mask$exclude
    retain <- mask$retain
    observed <- original$data$v5__clinical__hba1c
    gs_assert(identical(original$members$clinical__hba1c$column, "v5__clinical__hba1c") && identical(is.na(observed), 
        !available) && all(observed[available] == h[available]) && sum(available) == 113685L && sum(!available) == 
        9665L && sum(exclude) == 719L && sum(old_ws$complete & exclude) == 575L, "OBSERVED_RAW_HBA1C_EXCLUSION_PROOF_CHANGED")
    target_rows <- which(retain)
    z <- original
    z$data <- original$data[target_rows, , drop = FALSE]
    gs_assert(identical(names(original$risks), c("IHD", "T2D")), "PAIRED_ORIGINAL_ENDPOINTS_REQUIRED")
    z$risks <- lapply(original$risks, function(r) {
        gs_assert(identical(as.character(r$eid), as.character(original$data$eid)) && all(c("entry", "exit", "event", 
            "death_event") %in% names(r)) && all(is.finite(r$exit) & r$exit > 0 & r$entry == 0 & r$event %in% 0:1 & 
            r$death_event %in% 0:1 & r$event + r$death_event <= 1), "ORIGINAL_FULLFOLLOWUP_RISK_INVALID")
        r[target_rows, , drop = FALSE]
    })
    z$schema <- "GLYCAEMIA_SENSITIVITY_BASE11_B0_INPUTS_V1_20261106"
    z$target_hash <- gs_hash(as.character(z$data$eid))
    z$cluster_membership_hash <- gs_hash(paste(z$data$eid, z$data$cluster_id, sep = ":"))
    gs_assert(nrow(z$data) == 122631L && identical(z$target_hash, birth_private_setting("PRIVATE_HASH_7")) && length(unique(z$data$cluster_id)) == 
        112038L, "FILTERED_TARGET_HASH_OR_COUNT_CHANGED")
    aux <- gs_vec(z$selection_auxiliary_columns)
    required_aux <- unique(unlist(regmatches(z$selection_functions$b3_missing, gregexpr("aux_[A-Za-z0-9_]+", z$selection_functions$b3_missing))))
    gs_assert(length(required_aux) > 0L && all(required_aux %in% aux) && all(aux %in% names(z$data)) && !any(aux %in% 
        c(base, cols)), "FROZEN_SELECTION_AUXILIARIES_CHANGED")
    se <- new.env(parent = globalenv())
    for (name in names(z$selection_functions)) se[[name]] <- eval(parse(text = z$selection_functions[[name]]), 
        envir = se)
    progress("REFITTING_SELECTION_ON_FILTERED_TARGET", list(N_target = nrow(z$data)))
    ws <- se$b3_missing(z$data, cols, base)
    gs_assert(identical(ws$complete, complete.cases(z$data[, unique(c(base, cols)), drop = FALSE])) && identical(ws$complete, 
        old_ws$complete[target_rows]) && sum(ws$complete) == 94739L && length(ws$w) == sum(ws$complete) && all(is.finite(ws$w) & 
        ws$w > 0) && !is.null(ws$fit), "FRESH_TARGET_SELECTION_MODEL_OR_COMPLETE_MASK_INVALID")
    selection_fit <- ws$fit
    gs_assert(isTRUE(selection_fit$converged) && all(is.finite(stats::coef(selection_fit))) && stats::nobs(selection_fit) == 
        nrow(z$data), "SELECTION_MODEL_MUST_USE_FULL_FILTERED_TARGET")
    ws$design_weights <- rep(1, length(ws$w))
    ws$selection_weights <- as.numeric(ws$w)
    ws$selection_coefficients <- stats::coef(selection_fit)
    ws$selection_formula <- paste(deparse(stats::formula(selection_fit)), collapse = " ")
    ws$fit <- NULL
    d <- z$data[ws$complete, , drop = FALSE]
    risks <- lapply(z$risks, function(r) r[ws$complete, , drop = FALSE])
    sample_hash <- gs_hash(as.character(d$eid))
    gs_assert(identical(sample_hash, birth_private_setting("PRIVATE_HASH_8")) && all(d$v5__clinical__hba1c < 48) && 
        sum(risks$IHD$event) == 6291L && sum(risks$T2D$event) == 3005L && sum(risks$IHD$death_event) == 4551L && 
        sum(risks$T2D$death_event) == 4877L && sum(risks$IHD$event == 1 & risks$IHD$exit > 11) == 936L && sum(risks$T2D$event == 
        1 & risks$T2D$exit > 11) == 548L, "FILTERED_FIT_HASH_OR_FULLFOLLOWUP_ENDPOINT_COUNTS_CHANGED")
    gs_assert(!identical(as.numeric(ws$w), as.numeric(old_ws$w[retain[old_ws$complete]])), "OLD_IPW_SUBSETTING_INSTEAD_OF_NEW_TARGET_REFIT_REJECTED")
    clusters <- unique(as.character(z$data$cluster_id))
    index <- list(bootstrap_id = 0L, indices = seq_len(nrow(z$data)), cluster_multiplicity = rep.int(1L, length(clusters)), 
        target_hash = z$target_hash, cluster_membership_hash = z$cluster_membership_hash, bootstrap_seed = 20261023L, 
        RNG_kind = c("L'Ecuyer-CMRG", "Inversion", "Rejection"), bootstrap_sampling_performed = FALSE)
    z$sensitivity <- list(rule = "OBSERVED_RAW_BASELINE_HBA1C_GE48_REPORTABLE_EXCLUDED", original_target_hash = original$target_hash, 
        original_input_sha256 = pins$input, original_B0_index_sha256 = pins$bootstrap, original_B0_weights_sha256 = pins$weights, 
        raw_clinical_sha256 = pins$clinical, raw_value_field = "p30750_i0", raw_quality_field = "p30756_i0", unknown_HbA1c_retained_in_target = TRUE, 
        selection_refitted_on_filtered_target = TRUE, original_values_auxiliaries_and_fullfollowup_preserved = TRUE, 
        bootstrap_ids = 0L)
    files <- list(input = gs_new_rds(z, file.path(private, "LINE_INPUTS.private.rds")), index = gs_new_rds(index, 
        file.path(private, "bootstrap_0000.private.rds")), selection_fit = gs_new_rds(selection_fit, file.path(private, 
        "SELECTION_FIT.private.rds")))
    ws$pair_binding <- list(target_hash = z$target_hash, sample_hash = sample_hash, design_weight_hash = gs_hash(ws$design_weights), 
        selection_weight_hash = gs_hash(ws$selection_weights), sampling_hash = gs_hash(as.integer(index$indices)), 
        actual_index_sha256 = files$index$sha256)
    files$weights <- gs_new_rds(ws, file.path(private, "weights_0000.private.rds"))
    frozen <- readRDS(files$input$path)
    frozen_ws <- readRDS(files$weights$path)
    old_fit_ids <- as.character(original$data$eid[old_ws$complete])
    new_fit_ids <- as.character(d$eid)
    removed_rows <- old_ws$complete & !retain
    retained_fit_rows <- which(old_ws$complete & retain)
    columns_identical <- vapply(names(original$data), function(column) identical(frozen$data[[column]], original$data[[column]][target_rows]), 
        logical(1))
    risk_rows_identical <- vapply(names(original$risks), function(endpoint) identical(frozen$risks[[endpoint]], 
        original$risks[[endpoint]][target_rows, , drop = FALSE]) && identical(risks[[endpoint]], original$risks[[endpoint]][retained_fit_rows, 
        , drop = FALSE]), logical(1))
    definition_fields <- c("members", "base", "domains", "model_specifications")
    definition_hashes <- setNames(lapply(definition_fields, function(field) list(original = gs_hash(original[[field]]), 
        frozen = gs_hash(frozen[[field]]), identical = identical(original[[field]], frozen[[field]]))), definition_fields)
    qa <- list(schema = "GLYCAEMIA_BASE11_B0_ACTUAL_SAMPLE_RISK_QA_V1", passed = TRUE, original_fit_N = length(old_fit_ids), 
        retained_fit_N = length(new_fit_ids), removed_fit_N = length(setdiff(old_fit_ids, new_fit_ids)), added_fit_N = length(setdiff(new_fit_ids, 
            old_fit_ids)), removed_non_high_N = sum(removed_rows & !exclude), filtered_target_N = nrow(frozen$data), 
        retained_fit_exact_ordered_original_subset = identical(new_fit_ids, as.character(original$data$eid[retained_fit_rows])), 
        all_retained_target_columns_identical = all(columns_identical), retained_target_column_checks_N = length(columns_identical), 
        changed_retained_target_column_N = sum(!columns_identical), retained_entire_risk_rows_identical = as.list(risk_rows_identical), 
        entry_exit_event_death_and_other_risk_fields_unchanged = all(risk_rows_identical), scientific_definition_hashes = definition_hashes, 
        roundtrip = list(input_data_hash_matches = identical(gs_hash(frozen$data), gs_hash(z$data)), input_target_risk_hash_matches = identical(gs_hash(frozen$risks), 
            gs_hash(z$risks)), fit_sample_hash_matches = identical(gs_hash(as.character(frozen$data$eid[frozen_ws$complete])), 
            sample_hash), fit_risk_hash_matches = identical(gs_hash(lapply(frozen$risks, function(r) r[frozen_ws$complete, 
            , drop = FALSE])), gs_hash(risks)), weight_hash_matches = identical(gs_hash(as.numeric(frozen_ws$w)), 
            gs_hash(as.numeric(ws$w))), complete_selection_payload_hash_matches = identical(gs_hash(frozen_ws), 
            gs_hash(ws))), target_data_hash = gs_hash(frozen$data), target_risk_hash = gs_hash(frozen$risks), fit_data_hash = gs_hash(d), 
        fit_risk_hash = gs_hash(risks), fit_sample_hash = sample_hash, new_weight_hash = gs_hash(as.numeric(frozen_ws$w)), 
        participant_rows_in_reports = 0L)
    gs_assert(qa$original_fit_N == 95314L && qa$retained_fit_N == 94739L && qa$removed_fit_N == 575L && qa$added_fit_N == 
        0L && qa$removed_non_high_N == 0L && qa$filtered_target_N == 122631L && isTRUE(qa$retained_fit_exact_ordered_original_subset) && 
        isTRUE(qa$all_retained_target_columns_identical) && isTRUE(qa$entry_exit_event_death_and_other_risk_fields_unchanged) && 
        all(vapply(definition_hashes, function(x) isTRUE(x$identical) && identical(x$original, x$frozen), logical(1))) && 
        all(unlist(qa$roundtrip, use.names = FALSE)), "ACTUAL_SAMPLE_RISK_DIFF_OR_PRIVATE_FREEZE_ROUNDTRIP_FAILED")
    qa_file <- gs_new_json(qa, file.path(receipts, "SAMPLE_RISK_QA.json"))
    manifest <- list(schema = "GLYCAEMIA_SENSITIVITY_BASE11_B0_PREPARED_V1", complete = TRUE, utc = gs_utc(), b = 0L, 
        original_readonly_inputs = Map(function(path, sha256) list(path = path, sha256 = sha256), paths, pins), 
        target_hash = z$target_hash, fit_sample_hash = sample_hash, pair_binding = ws$pair_binding, data_hash = gs_hash(z$data), 
        observed_fit_data_hash = gs_hash(d), risk_hash = gs_hash(risks), original_selected_definition_hash = gs_hash(original[c("members", 
            "base", "domains", "model_specifications")]), filtered_selected_definition_hash = gs_hash(z[c("members", 
            "base", "domains", "model_specifications")]), selection_source_hash = z$selection_source_hash, selection_coefficient_hash = gs_hash(ws$selection_coefficients), 
        weight_hash = gs_hash(as.numeric(ws$w)), selection_formula = ws$selection_formula, N_original_target = 123350L, 
        N_excluded = 719L, N_unknown_retained = 9665L, N_target = nrow(z$data), N_original_fit = 95314L, N_fit = nrow(d), 
        N_standardization = nrow(d), N_target_families = length(clusters), endpoints = lapply(risks, gs_endpoint_counts), 
        IPW = list(min = min(ws$w), median = stats::median(ws$w), max = max(ws$w), sum = sum(ws$w), ESS = sum(ws$w)^2/sum(ws$w^2), 
            selection_fit_target_N = stats::nobs(selection_fit)), full_followup_unchanged = TRUE, report_horizon_years = 11L, 
        old_weight_reuse = FALSE, bootstrap_sampling_performed = FALSE, fresh_selection_refits = 1L, participant_rows_in_reports = 0L, 
        actual_sample_risk_qa = qa_file, files = files)
    gs_assert(identical(manifest$original_selected_definition_hash, manifest$filtered_selected_definition_hash), 
        "ROW_FILTER_ONLY_SELECTED_SCIENTIFIC_DEFINITION_CHANGED")
    gs_new_json(manifest, file.path(receipts, "PREPARATION_MANIFEST.json"))
    progress("PREPARATION_COMPLETE", list(N_target = nrow(z$data), N_fit = nrow(d)))
    list(z = z, d = d, w = as.numeric(ws$w), risks = risks, selection = ws, index = index, manifest = manifest, 
        files = files)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        gs_assert <- function(ok, reason) if (!isTRUE(ok)) stop(reason, call. = FALSE)
        gs_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)
        gs_sha <- function(path) digest::digest(file = path, algo = "sha256")
        gs_vec <- function(x) as.character(unlist(x, use.names = FALSE))
        gs_utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
        gs_new_json <- function(x, path) {
            gs_assert(!file.exists(path), "PRESERVE_EXISTING_JSON")
            tmp <- paste0(path, ".pending_", Sys.getpid())
            gs_assert(!file.exists(tmp), "PRESERVE_EXISTING_PENDING_JSON")
            jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
            gs_assert(!file.exists(path) && file.rename(tmp, path), "ADDITIVE_JSON_PUBLICATION_FAILED")
            invisible(list(path = path, sha256 = gs_sha(path)))
        }
        gs_new_rds <- function(x, path) {
            gs_assert(!file.exists(path), "PRESERVE_EXISTING_PRIVATE_RDS")
            tmp <- paste0(path, ".pending_", Sys.getpid())
            gs_assert(!file.exists(tmp), "PRESERVE_EXISTING_PENDING_RDS")
            saveRDS(x, tmp, version = 2L, compress = FALSE)
            gs_assert(!file.exists(path) && file.rename(tmp, path), "ADDITIVE_RDS_PUBLICATION_FAILED")
            invisible(list(path = path, sha256 = gs_sha(path), bytes = unname(file.info(path)$size)))
        }
        gs_norm <- function(path, mustWork = TRUE) {
            gs_assert(is.character(path) && length(path) == 1L && !is.na(path) && nzchar(path) && !grepl("(^|/)[.][.](/|$)", 
                path) && startsWith(path, "/"), "ABSOLUTE_EXACT_LINUX_PATH_REQUIRED")
            normalizePath(path, winslash = "/", mustWork = mustWork)
        }
        gs_inside <- function(path, root) identical(path, root) || startsWith(path, paste0(root, "/"))
        gs_check_config <- function(cfg) {
            gs_assert(identical(cfg$scope, "BASE11_B0_ONLY") && identical(gs_vec(cfg$bootstrap_ids), "0") && length(cfg$b) == 
                1L && is.numeric(cfg$b) && !is.na(cfg$b) && cfg$b == 0, "B0_ONLY_NO_BOOTSTRAP_DISPATCH")
            gs_assert(isTRUE(cfg$authorization$real_B0) && identical(cfg$authorization$bootstrap, FALSE), "EXPLICIT_REAL_B0_ONLY_AUTHORIZATION_REQUIRED")
            gs_assert(length(cfg$seed) == 1L && is.numeric(cfg$seed) && !is.na(cfg$seed) && cfg$seed == 20261103L, 
                "EXACT_ORIGINAL_BASE_B0_SEED_20261103_REQUIRED")
            for (key in c("n_mc", "uniform_dim", "block_size", "eval_year")) {
                expected <- switch(key, n_mc = 50L, uniform_dim = 11L, block_size = 512L, eval_year = 11L)
                if (!is.null(cfg[[key]])) gs_assert(identical(as.numeric(cfg[[key]]), as.numeric(expected)), "FROZEN_INTEGRATION_CONTRACT_CHANGED")
            }
            cfg
        }
        gs_resolve_config <- function(cfg) {
            cfg <- gs_check_config(cfg)
            keys <- c("input_root", "clinical_path", "source_root", "y_root", "state_root", "output_root", "adapter_root")
            for (key in keys) {
                override <- Sys.getenv(paste0("GS44_", toupper(key)))
                if (nzchar(override)) cfg[[key]] <- override
                gs_assert(!is.null(cfg[[key]]), paste0("MISSING_CONFIG_ROOT_", toupper(key)))
                cfg[[key]] <- gs_norm(cfg[[key]], mustWork = key != "state_root")
            }
            if (!file.exists(cfg$state_root)) cfg$state_root <- file.path(gs_norm(dirname(cfg$state_root)), basename(cfg$state_root))
            namespace <- "border_glycaemia_sensitivity_BASE11_B0_20261106_v1"
            gs_assert(grepl(paste0("/", namespace, "/"), cfg$output_root, fixed = TRUE) && basename(dirname(cfg$state_root)) == 
                namespace && basename(cfg$state_root) == "state", "ISOLATED_NEW_TASK_OUTPUT_AND_STATE_REQUIRED")
            reads <- unlist(cfg[c("input_root", "clinical_path", "source_root", "y_root", "adapter_root")])
            gs_assert(!any(vapply(reads, function(p) gs_inside(cfg$output_root, p) || gs_inside(cfg$state_root, 
                p) || gs_inside(p, cfg$output_root) || gs_inside(p, cfg$state_root), FALSE)), "NEW_OUTPUTS_MUST_NOT_OVERLAP_READONLY_INPUT_OR_SOURCE_ROOTS")
            gs_assert(!gs_inside(cfg$output_root, cfg$state_root) && !gs_inside(cfg$state_root, cfg$output_root), 
                "SEPARATE_NEW_STATE_AND_OUTPUT_TREES_REQUIRED")
            cfg
        }
        gs_input_pins <- function() list(input = birth_private_setting("PRIVATE_HASH_9"), bootstrap = birth_private_setting("PRIVATE_HASH_10"), 
            weights = birth_private_setting("PRIVATE_HASH_11"), clinical = birth_private_setting("PRIVATE_HASH_12"))
        gs_endpoint_counts <- function(r) list(N = nrow(r), events_full = sum(r$event == 1), events_11y = sum(r$event == 
            1 & r$exit <= 11), deaths_full = sum(r$death_event == 1), deaths_11y = sum(r$death_event == 1 & r$exit <= 
            11), events_after11 = sum(r$event == 1 & r$exit > 11), followup_max = max(r$exit))
        gs_raw_hba1c_mask <- function(value, reportability) {
            gs_assert(length(value) == length(reportability), "RAW_HBA1C_QUALITY_LENGTH_MISMATCH")
            h <- suppressWarnings(as.numeric(value))
            q <- suppressWarnings(as.numeric(reportability))
            available <- is.finite(h) & h > 0 & !is.na(q) & q == 1
            exclude <- available & h >= 48
            list(value = h, available = available, exclude = exclude, retain = !exclude)
        }
        gs_prepare_B0 <- function(cfg, private, receipts, progress = function(...) NULL) {
            gs_check_config(cfg)
            paths <- list(input = file.path(cfg$input_root, "LINE_INPUTS.private.rds"), bootstrap = file.path(cfg$input_root, 
                "draws", "bootstrap_0000.private.rds"), weights = file.path(cfg$input_root, "draws", "weights_0000.private.rds"), 
                clinical = cfg$clinical_path)
            pins <- gs_input_pins()
            for (key in names(paths)) {
                gs_assert(file.exists(paths[[key]]), paste0("MISSING_INPUT_PATH_", toupper(key)))
                if (!is.null(cfg[[paste0(key, "_sha256")]])) gs_assert(identical(cfg[[paste0(key, "_sha256")]], 
                  pins[[key]]), "IMMUTABLE_INPUT_PIN_OVERRIDE_REJECTED")
                gs_assert(identical(gs_sha(paths[[key]]), pins[[key]]), paste0("FROZEN_INPUT_HASH_CHANGED_", toupper(key)))
            }
            progress("INPUT_PINS_VERIFIED")
            original <- readRDS(paths$input)
            old_index <- readRDS(paths$bootstrap)
            old_ws <- readRDS(paths$weights)
            ids <- c("smoking", "bmi", "weight", "fat_free_mass", "fat_mass", "hip", "waist", "clinical__hdl", 
                "clinical__tg", "clinical__glucose", "clinical__hba1c")
            base <- c("B", "age", "sex", "ethnicity", "center", "assessment_year", "S")
            order <- c("SMOKING", "BODY_COMPOSITION", "CLINICAL_LIPIDS", "GLYCAEMIA")
            cols <- vapply(original$members, `[[`, "", "column")
            gs_assert(identical(original$schema, "SELECTED44_LINE_FROZEN_INPUTS_V1") && identical(original$line, 
                "BASE") && nrow(original$data) == 123350L && identical(names(original$members), ids) && identical(gs_vec(original$base), 
                base) && original$uniform_dim == 11L && identical(names(original$domains), order) && identical(gs_vec(original$model_specifications$S44_BASE11$allocation_order), 
                order) && identical(gs_vec(original$model_specifications$S44_BASE11$members), ids) && !isTRUE(original$model_specifications$S44_BASE11$total_only), 
                "FROZEN_BASE11_SCIENTIFIC_CONTRACT_CHANGED")
            gs_assert(!anyDuplicated(original$data$eid) && !anyNA(original$data$cluster_id) && identical(gs_hash(as.character(original$data$eid)), 
                original$target_hash) && identical(original$target_hash, birth_private_setting("PRIVATE_HASH_13")) && 
                identical(gs_hash(original$selection_functions), original$selection_source_hash), "ORIGINAL_TARGET_OR_SELECTION_SOURCE_BINDING_CHANGED")
            gs_assert(old_index$bootstrap_id == 0L && identical(as.integer(old_index$indices), seq_len(nrow(original$data))) && 
                identical(old_index$target_hash, original$target_hash) && identical(old_index$cluster_membership_hash, 
                original$cluster_membership_hash) && is.logical(old_ws$complete) && !anyNA(old_ws$complete) && 
                identical(old_ws$complete, complete.cases(original$data[, unique(c(base, cols)), drop = FALSE])) && 
                sum(old_ws$complete) == 95314L && identical(gs_hash(as.character(original$data$eid[old_ws$complete])), 
                birth_private_setting("PRIVATE_HASH_14")), "ORIGINAL_B0_IDENTITY_INDEX_AND_COMPLETE_MASK_CHANGED")
            gs_assert(identical(old_ws$pair_binding$target_hash, original$target_hash) && identical(old_ws$pair_binding$actual_index_sha256, 
                pins$bootstrap) && identical(gs_hash(as.numeric(old_ws$w)), old_ws$pair_binding$selection_weight_hash) && 
                identical(gs_hash(as.integer(old_index$indices)), old_ws$pair_binding$sampling_hash), "ORIGINAL_B0_WEIGHT_INDEX_PROVENANCE_CHANGED")
            header <- names(utils::read.csv(paths$clinical, nrows = 0L, check.names = FALSE))
            raw_cols <- c("eid", "p30750_i0", "p30756_i0")
            gs_assert(all(raw_cols %in% header), "RAW_BASELINE_HBA1C_AND_REPORTABILITY_COLUMNS_MISSING")
            classes <- rep("NULL", length(header))
            classes[match(raw_cols, header)] <- "character"
            raw <- utils::read.csv(paths$clinical, colClasses = classes, check.names = FALSE)
            gs_assert(!anyDuplicated(raw$eid), "DUPLICATE_RAW_CLINICAL_IDENTITY")
            j <- match(as.character(original$data$eid), raw$eid)
            gs_assert(!anyNA(j), "RAW_CLINICAL_TARGET_ALIGNMENT_FAILED")
            mask <- gs_raw_hba1c_mask(raw$p30750_i0[j], raw$p30756_i0[j])
            h <- mask$value
            available <- mask$available
            exclude <- mask$exclude
            retain <- mask$retain
            observed <- original$data$v5__clinical__hba1c
            gs_assert(identical(original$members$clinical__hba1c$column, "v5__clinical__hba1c") && identical(is.na(observed), 
                !available) && all(observed[available] == h[available]) && sum(available) == 113685L && sum(!available) == 
                9665L && sum(exclude) == 719L && sum(old_ws$complete & exclude) == 575L, "OBSERVED_RAW_HBA1C_EXCLUSION_PROOF_CHANGED")
            target_rows <- which(retain)
            z <- original
            z$data <- original$data[target_rows, , drop = FALSE]
            gs_assert(identical(names(original$risks), c("IHD", "T2D")), "PAIRED_ORIGINAL_ENDPOINTS_REQUIRED")
            z$risks <- lapply(original$risks, function(r) {
                gs_assert(identical(as.character(r$eid), as.character(original$data$eid)) && all(c("entry", "exit", 
                  "event", "death_event") %in% names(r)) && all(is.finite(r$exit) & r$exit > 0 & r$entry == 0 & 
                  r$event %in% 0:1 & r$death_event %in% 0:1 & r$event + r$death_event <= 1), "ORIGINAL_FULLFOLLOWUP_RISK_INVALID")
                r[target_rows, , drop = FALSE]
            })
            z$schema <- "GLYCAEMIA_SENSITIVITY_BASE11_B0_INPUTS_V1_20261106"
            z$target_hash <- gs_hash(as.character(z$data$eid))
            z$cluster_membership_hash <- gs_hash(paste(z$data$eid, z$data$cluster_id, sep = ":"))
            gs_assert(nrow(z$data) == 122631L && identical(z$target_hash, birth_private_setting("PRIVATE_HASH_15")) && 
                length(unique(z$data$cluster_id)) == 112038L, "FILTERED_TARGET_HASH_OR_COUNT_CHANGED")
            aux <- gs_vec(z$selection_auxiliary_columns)
            required_aux <- unique(unlist(regmatches(z$selection_functions$b3_missing, gregexpr("aux_[A-Za-z0-9_]+", 
                z$selection_functions$b3_missing))))
            gs_assert(length(required_aux) > 0L && all(required_aux %in% aux) && all(aux %in% names(z$data)) && 
                !any(aux %in% c(base, cols)), "FROZEN_SELECTION_AUXILIARIES_CHANGED")
            se <- new.env(parent = globalenv())
            for (name in names(z$selection_functions)) se[[name]] <- eval(parse(text = z$selection_functions[[name]]), 
                envir = se)
            progress("REFITTING_SELECTION_ON_FILTERED_TARGET", list(N_target = nrow(z$data)))
            ws <- se$b3_missing(z$data, cols, base)
            gs_assert(identical(ws$complete, complete.cases(z$data[, unique(c(base, cols)), drop = FALSE])) && 
                identical(ws$complete, old_ws$complete[target_rows]) && sum(ws$complete) == 94739L && length(ws$w) == 
                sum(ws$complete) && all(is.finite(ws$w) & ws$w > 0) && !is.null(ws$fit), "FRESH_TARGET_SELECTION_MODEL_OR_COMPLETE_MASK_INVALID")
            selection_fit <- ws$fit
            gs_assert(isTRUE(selection_fit$converged) && all(is.finite(stats::coef(selection_fit))) && stats::nobs(selection_fit) == 
                nrow(z$data), "SELECTION_MODEL_MUST_USE_FULL_FILTERED_TARGET")
            ws$design_weights <- rep(1, length(ws$w))
            ws$selection_weights <- as.numeric(ws$w)
            ws$selection_coefficients <- stats::coef(selection_fit)
            ws$selection_formula <- paste(deparse(stats::formula(selection_fit)), collapse = " ")
            ws$fit <- NULL
            d <- z$data[ws$complete, , drop = FALSE]
            risks <- lapply(z$risks, function(r) r[ws$complete, , drop = FALSE])
            sample_hash <- gs_hash(as.character(d$eid))
            gs_assert(identical(sample_hash, birth_private_setting("PRIVATE_HASH_16")) && all(d$v5__clinical__hba1c < 
                48) && sum(risks$IHD$event) == 6291L && sum(risks$T2D$event) == 3005L && sum(risks$IHD$death_event) == 
                4551L && sum(risks$T2D$death_event) == 4877L && sum(risks$IHD$event == 1 & risks$IHD$exit > 11) == 
                936L && sum(risks$T2D$event == 1 & risks$T2D$exit > 11) == 548L, "FILTERED_FIT_HASH_OR_FULLFOLLOWUP_ENDPOINT_COUNTS_CHANGED")
            gs_assert(!identical(as.numeric(ws$w), as.numeric(old_ws$w[retain[old_ws$complete]])), "OLD_IPW_SUBSETTING_INSTEAD_OF_NEW_TARGET_REFIT_REJECTED")
            clusters <- unique(as.character(z$data$cluster_id))
            index <- list(bootstrap_id = 0L, indices = seq_len(nrow(z$data)), cluster_multiplicity = rep.int(1L, 
                length(clusters)), target_hash = z$target_hash, cluster_membership_hash = z$cluster_membership_hash, 
                bootstrap_seed = 20261023L, RNG_kind = c("L'Ecuyer-CMRG", "Inversion", "Rejection"), bootstrap_sampling_performed = FALSE)
            z$sensitivity <- list(rule = "OBSERVED_RAW_BASELINE_HBA1C_GE48_REPORTABLE_EXCLUDED", original_target_hash = original$target_hash, 
                original_input_sha256 = pins$input, original_B0_index_sha256 = pins$bootstrap, original_B0_weights_sha256 = pins$weights, 
                raw_clinical_sha256 = pins$clinical, raw_value_field = "p30750_i0", raw_quality_field = "p30756_i0", 
                unknown_HbA1c_retained_in_target = TRUE, selection_refitted_on_filtered_target = TRUE, original_values_auxiliaries_and_fullfollowup_preserved = TRUE, 
                bootstrap_ids = 0L)
            files <- list(input = gs_new_rds(z, file.path(private, "LINE_INPUTS.private.rds")), index = gs_new_rds(index, 
                file.path(private, "bootstrap_0000.private.rds")), selection_fit = gs_new_rds(selection_fit, file.path(private, 
                "SELECTION_FIT.private.rds")))
            ws$pair_binding <- list(target_hash = z$target_hash, sample_hash = sample_hash, design_weight_hash = gs_hash(ws$design_weights), 
                selection_weight_hash = gs_hash(ws$selection_weights), sampling_hash = gs_hash(as.integer(index$indices)), 
                actual_index_sha256 = files$index$sha256)
            files$weights <- gs_new_rds(ws, file.path(private, "weights_0000.private.rds"))
            frozen <- readRDS(files$input$path)
            frozen_ws <- readRDS(files$weights$path)
            old_fit_ids <- as.character(original$data$eid[old_ws$complete])
            new_fit_ids <- as.character(d$eid)
            removed_rows <- old_ws$complete & !retain
            retained_fit_rows <- which(old_ws$complete & retain)
            columns_identical <- vapply(names(original$data), function(column) identical(frozen$data[[column]], 
                original$data[[column]][target_rows]), logical(1))
            risk_rows_identical <- vapply(names(original$risks), function(endpoint) identical(frozen$risks[[endpoint]], 
                original$risks[[endpoint]][target_rows, , drop = FALSE]) && identical(risks[[endpoint]], original$risks[[endpoint]][retained_fit_rows, 
                , drop = FALSE]), logical(1))
            definition_fields <- c("members", "base", "domains", "model_specifications")
            definition_hashes <- setNames(lapply(definition_fields, function(field) list(original = gs_hash(original[[field]]), 
                frozen = gs_hash(frozen[[field]]), identical = identical(original[[field]], frozen[[field]]))), 
                definition_fields)
            qa <- list(schema = "GLYCAEMIA_BASE11_B0_ACTUAL_SAMPLE_RISK_QA_V1", passed = TRUE, original_fit_N = length(old_fit_ids), 
                retained_fit_N = length(new_fit_ids), removed_fit_N = length(setdiff(old_fit_ids, new_fit_ids)), 
                added_fit_N = length(setdiff(new_fit_ids, old_fit_ids)), removed_non_high_N = sum(removed_rows & 
                  !exclude), filtered_target_N = nrow(frozen$data), retained_fit_exact_ordered_original_subset = identical(new_fit_ids, 
                  as.character(original$data$eid[retained_fit_rows])), all_retained_target_columns_identical = all(columns_identical), 
                retained_target_column_checks_N = length(columns_identical), changed_retained_target_column_N = sum(!columns_identical), 
                retained_entire_risk_rows_identical = as.list(risk_rows_identical), entry_exit_event_death_and_other_risk_fields_unchanged = all(risk_rows_identical), 
                scientific_definition_hashes = definition_hashes, roundtrip = list(input_data_hash_matches = identical(gs_hash(frozen$data), 
                  gs_hash(z$data)), input_target_risk_hash_matches = identical(gs_hash(frozen$risks), gs_hash(z$risks)), 
                  fit_sample_hash_matches = identical(gs_hash(as.character(frozen$data$eid[frozen_ws$complete])), 
                    sample_hash), fit_risk_hash_matches = identical(gs_hash(lapply(frozen$risks, function(r) r[frozen_ws$complete, 
                    , drop = FALSE])), gs_hash(risks)), weight_hash_matches = identical(gs_hash(as.numeric(frozen_ws$w)), 
                    gs_hash(as.numeric(ws$w))), complete_selection_payload_hash_matches = identical(gs_hash(frozen_ws), 
                    gs_hash(ws))), target_data_hash = gs_hash(frozen$data), target_risk_hash = gs_hash(frozen$risks), 
                fit_data_hash = gs_hash(d), fit_risk_hash = gs_hash(risks), fit_sample_hash = sample_hash, new_weight_hash = gs_hash(as.numeric(frozen_ws$w)), 
                participant_rows_in_reports = 0L)
            gs_assert(qa$original_fit_N == 95314L && qa$retained_fit_N == 94739L && qa$removed_fit_N == 575L && 
                qa$added_fit_N == 0L && qa$removed_non_high_N == 0L && qa$filtered_target_N == 122631L && isTRUE(qa$retained_fit_exact_ordered_original_subset) && 
                isTRUE(qa$all_retained_target_columns_identical) && isTRUE(qa$entry_exit_event_death_and_other_risk_fields_unchanged) && 
                all(vapply(definition_hashes, function(x) isTRUE(x$identical) && identical(x$original, x$frozen), 
                  logical(1))) && all(unlist(qa$roundtrip, use.names = FALSE)), "ACTUAL_SAMPLE_RISK_DIFF_OR_PRIVATE_FREEZE_ROUNDTRIP_FAILED")
            qa_file <- gs_new_json(qa, file.path(receipts, "SAMPLE_RISK_QA.json"))
            manifest <- list(schema = "GLYCAEMIA_SENSITIVITY_BASE11_B0_PREPARED_V1", complete = TRUE, utc = gs_utc(), 
                b = 0L, original_readonly_inputs = Map(function(path, sha256) list(path = path, sha256 = sha256), 
                  paths, pins), target_hash = z$target_hash, fit_sample_hash = sample_hash, pair_binding = ws$pair_binding, 
                data_hash = gs_hash(z$data), observed_fit_data_hash = gs_hash(d), risk_hash = gs_hash(risks), original_selected_definition_hash = gs_hash(original[c("members", 
                  "base", "domains", "model_specifications")]), filtered_selected_definition_hash = gs_hash(z[c("members", 
                  "base", "domains", "model_specifications")]), selection_source_hash = z$selection_source_hash, 
                selection_coefficient_hash = gs_hash(ws$selection_coefficients), weight_hash = gs_hash(as.numeric(ws$w)), 
                selection_formula = ws$selection_formula, N_original_target = 123350L, N_excluded = 719L, N_unknown_retained = 9665L, 
                N_target = nrow(z$data), N_original_fit = 95314L, N_fit = nrow(d), N_standardization = nrow(d), 
                N_target_families = length(clusters), endpoints = lapply(risks, gs_endpoint_counts), IPW = list(min = min(ws$w), 
                  median = stats::median(ws$w), max = max(ws$w), sum = sum(ws$w), ESS = sum(ws$w)^2/sum(ws$w^2), 
                  selection_fit_target_N = stats::nobs(selection_fit)), full_followup_unchanged = TRUE, report_horizon_years = 11L, 
                old_weight_reuse = FALSE, bootstrap_sampling_performed = FALSE, fresh_selection_refits = 1L, participant_rows_in_reports = 0L, 
                actual_sample_risk_qa = qa_file, files = files)
            gs_assert(identical(manifest$original_selected_definition_hash, manifest$filtered_selected_definition_hash), 
                "ROW_FILTER_ONLY_SELECTED_SCIENTIFIC_DEFINITION_CHANGED")
            gs_new_json(manifest, file.path(receipts, "PREPARATION_MANIFEST.json"))
            progress("PREPARATION_COMPLETE", list(N_target = nrow(z$data), N_fit = nrow(d)))
            list(z = z, d = d, w = as.numeric(ws$w), risks = risks, selection = ws, index = index, manifest = manifest, 
                files = files)
        }
    }), envir = context)
}

