# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

ff_p_assert <- function(ok, text) if (!isTRUE(ok)) stop(text, call. = FALSE)

ff_p_env <- function(name, default = "") {
    x <- Sys.getenv(name, "")
    if (nzchar(x)) 
        x
    else default
}

ff_p_sha <- function(path) {
    ff_p_assert(file.exists(path) && !dir.exists(path), paste("FF_REQUIRED_FILE_MISSING", path))
    digest::digest(file = path, algo = "sha256")
}

ff_p_cfg <- function(root) jsonlite::fromJSON(file.path(root, "config/runtime.json"), simplifyVector = FALSE)

ff_p_host <- function() ff_p_env("FF_HOST", if (Sys.info()[["sysname"]] == "Darwin") "macmini" else "winpc")

ff_p_rebase <- function(path) {
    prefix <- birth_private_setting("PRIVATE_INPUT")
    input <- ff_p_env("FF_INPUT_ROOT", birth_private_setting("PRIVATE_BIRTH"))
    if (startsWith(path, prefix)) 
        file.path(input, substring(path, nchar(prefix) + 1L))
    else path
}

ff_p_run <- function(cfg, create = FALSE) {
    host <- ff_p_host()
    run <- ff_p_env("FF_RUN_DIR", cfg$run_paths[[host]])
    ff_p_assert(grepl("/UKB_results/birth/border_fullfollowup_b100_20261027/run_[A-Za-z0-9_.-]+$", run) && !grepl("(^|/)[.][.](/|$)", 
        run) && (startsWith(run, birth_private_setting("PRIVATE_INPUT")) || startsWith(run, birth_private_setting("PRIVATE_INPUT"))), 
        "FF_ISOLATED_G_VERSIONED_RUN_REQUIRED")
    if (create && !dir.exists(run)) 
        ff_p_assert(dir.create(run, recursive = TRUE), "FF_NEW_RUN_CREATE_FAILED")
    normalizePath(run, mustWork = TRUE)
}

ff_p_R_hash <- function(root) {
    files <- sort(list.files(file.path(root, "f"), pattern = "[.]R$", full.names = TRUE))
    files <- c(files, sort(list.files(file.path(root, "src"), pattern = "[.]c$", full.names = TRUE)))
    digest::digest(paste(basename(files), vapply(files, ff_p_sha, ""), sep = ":", collapse = "\n"), algo = "sha256", 
        serialize = FALSE)
}

ff_p_auth <- function() ff_p_assert(identical(ff_p_env("FF_REAL_EXECUTION_AUTHORIZED"), "YES_FOR_THIS_RUN"), "FF_EXPLICIT_MANUAL_REAL_AUTHORIZATION_REQUIRED")

ff_p_exposure_alignment <- function(frozen, original) {
    ff_p_assert(is.numeric(frozen) && is.numeric(original) && length(frozen) == length(original) && length(frozen) > 
        0L && all(is.finite(frozen)) && all(is.finite(original)) && all(frozen >= 0 & frozen <= 1) && all(original >= 
        0 & original <= 1), "FF_ORIGINAL_CONTINUOUS_B_INVALID")
    difference <- max(abs(frozen - original))
    tolerance <- 8 * .Machine$double.eps
    ff_p_assert(difference <= tolerance, "FF_ORIGINAL_CONTINUOUS_B_CONFLICT")
    list(max_abs_difference = difference, nonidentical_double_values = sum(frozen != original), tolerance = tolerance, 
        tolerance_basis = "EIGHT_IEEE_DOUBLE_MACHINE_EPS_NO_VALUE_CHANGE", frozen_V5_B_preserved = TRUE)
}

ff_p_v5 <- function(root, cfg) {
    ff_p_auth()
    host <- ff_p_host()
    oldroot <- ff_p_env("FF_V5_ROOT", cfg$previous_code_paths[[host]])
    oldrun <- ff_p_env("FF_V5_RUN", cfg$previous_run_paths[[host]])
    adapter <- file.path(oldroot, "f/1.adapter.R")
    ff_p_assert(file.exists(adapter), paste("FF_EXISTING_V5_ADAPTER_MISSING", adapter))
    v5 <- new.env(parent = globalenv())
    sys.source(adapter, envir = v5)
    for (name in c("v5_bind", "v5_check", "v5_hash", "v5_sha", "v5_stage", "v5_utc")) assign(name, v5[[name]], 
        envir = environment(ff_p_v5))
    before <- Sys.getenv("V5_REAL_EXECUTION_AUTHORIZED", unset = NA_character_)
    on.exit(if (is.na(before)) Sys.unsetenv("V5_REAL_EXECUTION_AUTHORIZED") else Sys.setenv(V5_REAL_EXECUTION_AUTHORIZED = before), 
        add = TRUE)
    Sys.setenv(V5_REAL_EXECUTION_AUTHORIZED = "YES_FOR_THIS_RUN")
    ctx <- v5$v5_load(oldroot, oldrun)
    list(ctx = ctx, api = v5, root = oldroot, run = oldrun)
}

ff_p_endpoint <- function(parent_ids, original, risk, baseline, death, loss, admin) {
    j <- match(parent_ids, as.character(risk$eid))
    eligible <- !is.na(j)
    n <- length(parent_ids)
    r <- data.frame(eid = parent_ids, entry = rep(0, n), exit = rep(NA_real_, n), event = rep(NA_integer_, n), 
        death_event = rep(NA_integer_, n))
    r$entry[eligible] <- risk$.entry_time[j[eligible]]
    r$exit[eligible] <- risk$.exit_time[j[eligible]]
    r$event[eligible] <- risk$.event[j[eligible]]
    death_years <- as.numeric(death - baseline)/365.25
    loss_years <- as.numeric(loss - baseline)/365.25
    admin_years <- as.numeric(as.Date(admin) - baseline)/365.25
    r$death_event[eligible] <- as.integer(r$event[eligible] == 0L & !is.na(death_years[eligible]) & death_years[eligible] > 
        0 & death_years[eligible] == r$exit[eligible])
    reason <- rep("NOT_IN_ORIGINAL_RISKSET", n)
    reason[eligible & r$event == 1L] <- "TARGET_EVENT"
    reason[eligible & r$death_event == 1L] <- "DEATH_BEFORE_TARGET_EVENT"
    nd <- eligible & r$event == 0L & r$death_event == 0L
    nd[is.na(nd)] <- FALSE
    losscap <- !is.na(loss_years) & loss_years == r$exit
    admincap <- admin_years == r$exit
    reason[nd & losscap & !admincap] <- "RECORDED_NONDEATH_LOSS"
    reason[nd & !losscap & admincap] <- "ADMINISTRATIVE_CUTOFF"
    reason[nd & losscap & admincap] <- "LOSS_AND_ADMIN_SAME_DATE"
    ff_p_assert(all(reason[eligible] != "NOT_IN_ORIGINAL_RISKSET"), "FF_ORIGINAL_CENSOR_REASON_NOT_RECONSTRUCTED")
    list(risk = r, eligibility = eligible, censor_reason = reason)
}

ff_prepare <- function(root) {
    ff_p_auth()
    ff_p_assert(identical(ff_p_host(), "winpc") && Sys.info()[["sysname"]] != "Darwin", "FF_PARTICIPANT_PREPARATION_WINPC_ONLY")
    root <- normalizePath(root, mustWork = TRUE)
    cfg <- ff_p_cfg(root)
    run <- ff_p_run(cfg, TRUE)
    ff_p_assert(identical(as.integer(unlist(cfg$bootstrap_ids)), 0:100) && as.integer(cfg$max_bootstrap_id) == 
        100L, "FF_FIXED_B0_PLUS_1_TO_100_REQUIRED")
    private <- file.path(run, "private")
    if (!dir.exists(private)) 
        ff_p_assert(dir.create(private), "FF_PRIVATE_CREATE_FAILED")
    ff_p_assert(!file.exists(file.path(private, "FF_INPUT_LOCK.json")), "FF_PRESERVE_EXISTING_PREPARATION")
    attempt <- file.path(private, paste0("preparation_attempt_", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), 
        "_", Sys.getpid()))
    ff_p_assert(!dir.exists(attempt) && dir.create(attempt), "FF_UNIQUE_PREPARATION_ATTEMPT_REQUIRED")
    old <- ff_p_v5(root, cfg)
    ctx <- old$ctx
    v5 <- old$api
    scope_path <- file.path(root, cfg$active_approved_path)
    ff_p_assert(identical(ff_p_sha(scope_path), cfg$active_approved_csv_sha256), "FF_121_REGISTRY_SHA_CHANGED")
    scope <- read.csv(scope_path, colClasses = "character", check.names = FALSE)
    ff_p_assert(nrow(scope) == 121L && !anyDuplicated(scope$indicator_id) && !any(scope$category == "CMR") && all(scope$indicator_id %in% 
        ctx$scope$indicator_id), "FF_EXACT_EXISTING_121_NON_CMR_REQUIRED")
    main <- ff_p_rebase(cfg$previous_main$path)
    ff_p_assert(identical(ff_p_sha(main), cfg$previous_main$sha256), "FF_OLD_MAIN_SHA_CHANGED")
    incident_path <- ff_p_env("FF_INCIDENT_PREPARED", ff_p_rebase(cfg$incident$prepared))
    incident_sha <- ff_p_sha(incident_path)
    inc <- readRDS(incident_path)
    while (is.list(inc) && !is.null(inc$value)) inc <- inc$value
    ff_p_assert(all(c("data", "events", "baseline") %in% names(inc)), "FF_ORIGINAL_INCIDENT_SCHEMA_REQUIRED")
    code_root <- ff_p_env("FF_INCIDENT_CODE_ROOT", cfg$incident$code_root)
    cox_file <- file.path(code_root, "vendor/incident/cox.R")
    coverage_file <- file.path(code_root, "vendor/incident/provenance.R")
    cox_sha <- ff_p_sha(cox_file)
    coverage_sha <- ff_p_sha(coverage_file)
    ff_p_assert(grepl("attr\\(out, \"cutoff\"\\) <- \"2022-05-31\"", paste(readLines(coverage_file, warn = FALSE), 
        collapse = "\n")), "FF_FROZEN_UNIFORM_COVERAGE_CONTRACT_CHANGED")
    engine <- new.env(parent = globalenv())
    sys.source(cox_file, envir = engine)
    ids <- as.character(ctx$d_full$eid)
    match_source <- match(ids, as.character(inc$data$eid))
    ff_p_assert(!anyNA(match_source) && !anyDuplicated(as.character(inc$data$eid)), "FF_ORIGINAL_INCIDENT_NOT_COVERING_FROZEN_PARENT")
    sd <- inc$data[match_source, , drop = FALSE]
    for (name in c("age", "sex", "ethnicity", "center", "assessment_year", "S")) ff_p_assert(identical(as.character(ctx$d_full[[name]]), 
        as.character(sd[[name]])), paste("FF_ORIGINAL_COVARIATE_CONFLICT", name))
    exposure_alignment <- ff_p_exposure_alignment(as.numeric(ctx$d_full$B), as.numeric(sd$border_relative))
    admin <- cfg$incident$admin
    baseline <- as.Date(sd$date_attend)
    death <- as.Date(sd$date_death_uncapped)
    loss <- as.Date(sd$date_lost_registry)
    potential <- as.numeric(as.Date(admin) - baseline)
    ff_p_assert(all(is.finite(potential) & potential > 0), "FF_POTENTIAL_ADMIN_COVERAGE_MISSING")
    tau <- floor(min(potential/365.25))
    ff_p_assert(tau > 0, "FF_NO_POSITIVE_COMMON_ADMIN_REPORT_TIME")
    eps <- alignment <- coverage <- dates <- list()
    for (e in as.character(unlist(cfg$endpoints))) {
        ff_p_assert(e %in% names(inc$events) && e %in% names(inc$baseline), "FF_ORIGINAL_ENDPOINT_MISSING")
        raw <- engine$ip_riskset(inc$data, inc$events[[e]], inc$baseline[[e]], "primary", admin)$data
        z <- ff_p_endpoint(ids, inc, raw, baseline, death, loss, admin)
        oldep <- ctx$endpoints[[e]]
        oldrisk <- oldep$risk
        same <- z$eligibility & oldep$eligibility
        ok <- all(!oldep$eligibility | z$eligibility)
        rh <- digest::digest(z$risk, algo = "sha256")
        eps[[e]] <- list(status = if (ok) "AVAILABLE_FROZEN" else "BLOCKED_ALIGNMENT", reason = if (ok) "" else "V5_PARTICIPANTS_OUTSIDE_ORIGINAL_COMPLETED_INCIDENT_RISKSET_NO_SILENT_SAMPLE_CHANGE", 
            eligibility = oldep$eligibility, incident_eligibility = z$eligibility, risk = z$risk, endpoint_version = cfg$incident$definition, 
            risk_hash = rh, admin = admin, censor_reason = z$censor_reason, source = incident_path, source_sha256 = incident_sha)
        alignment[[e]] <- data.frame(endpoint_id = e, source_file = incident_path, source_sha256 = incident_sha, 
            frozen_incident_version = cfg$incident$definition, previous_V5_version = oldep$endpoint_version, N_original_incident_riskset = nrow(raw), 
            N_parent = length(ids), N_parent_ID_matched = sum(!is.na(match_source)), exposure_max_abs_difference = exposure_alignment$max_abs_difference, 
            exposure_numeric_tolerance = exposure_alignment$tolerance, frozen_V5_exposure_preserved = TRUE, N_V5_parent_eligible = sum(oldep$eligibility), 
            N_original_eligible_same_parent = sum(z$eligibility), N_risk_qualification_conflict = sum(oldep$eligibility & 
                !z$eligibility), N_same_ID_time_event_compared = sum(same), event_differences_vs_V5 = sum(z$risk$event[same] != 
                oldrisk$event[same]), time_differences_vs_V5 = sum(z$risk$exit[same] != oldrisk$exit[same]), original_events_same_ID = sum(z$risk$event[same]), 
            V5_events_same_ID = sum(oldrisk$event[same]), original_competing_deaths_same_ID = sum(z$risk$death_event[same]), 
            original_nondeath_censor_same_ID = sum(z$censor_reason[same] %in% c("RECORDED_NONDEATH_LOSS", "ADMINISTRATIVE_CUTOFF", 
                "LOSS_AND_ADMIN_SAME_DATE")), event_difference_vs_authoritative_incident = 0L, time_difference_vs_authoritative_incident = 0L, 
            dates_rule = "ORIGINAL_EVENT_ON_CENSOR_DATE_WINS;DEATH_ONLY_ADDED_WHERE_ORIGINAL_EVENT_ZERO", time_unit = "ELAPSED_DAYS_DIVIDED_BY_365.25_YEARS", 
            time_zero = "FROZEN_ORIGINAL_BASELINE_DATE", same_A_population_allowed = ok, Y_version_change_disclosed = !identical(oldep$endpoint_version, 
                cfg$incident$definition), mediator_missingness = "PER_INDICATOR_FROZEN_COMPLETE_RECORD_REQUIREMENTS_NOT_COMMON_121_INTERSECTION")
        coverage[[e]] <- data.frame(endpoint_id = e, coverage_end = admin, min_potential_admin_days = min(potential), 
            max_potential_admin_days = max(potential), raw_common_potential_years = min(potential/365.25), tau_report = tau, 
            source_contract = coverage_file, source_contract_sha256 = coverage_sha, coverage = "FROZEN_CONSERVATIVE_INSTITUTIONAL_UNIFORM_HOSPITAL_CUTOFF_NOT_INDIVIDUAL_LINKAGE_CERTIFICATION", 
            limitation = "ORIGINAL_UKB_EXTRACTION_RELEASE_ID_NOT_FOUND;EXACT_LOCAL_SNAPSHOT_HASHES_USED", uses_actual_event_death_loss_minimum = FALSE)
        dates[[e]] <- as.Date(inc$events[[e]][match_source])
    }
    cohort_alignment <- list()
    for (id in scope$indicator_id) {
        ca <- ctx$features[[id]]
        base <- unique(c(as.character(unlist(ctx$cfg$base_covariates)), ca$technical_covariates))
        cc <- stats::complete.cases(ctx$d_full[, unique(c(base, ca$column)), drop = FALSE])
        for (e in as.character(unlist(cfg$endpoints))) {
            target <- ctx$eligibility[[id]] & ctx$endpoints[[e]]$eligibility
            keep <- target & cc
            compare <- keep & eps[[e]]$incident_eligibility
            oldrisk <- ctx$endpoints[[e]]$risk
            newrisk <- eps[[e]]$risk
            cohort_alignment[[paste(id, e)]] <- data.frame(indicator_id = id, endpoint_id = e, N_frozen_parent = length(ids), 
                N_frozen_indicator_endpoint_target = sum(target), N_original_V5_complete_records = sum(keep), N_fullfollowup_same_ID_complete_records = sum(compare), 
                N_indicator_or_applicability_excluded = sum(!target), N_invalid_mediator_or_covariate_records = sum(target & 
                  !cc), N_original_incident_ineligible = sum(keep & !eps[[e]]$incident_eligibility), source_incident_events_in_complete_records = sum(newrisk$event[compare]), 
                source_competing_deaths_in_complete_records = sum(newrisk$death_event[compare]), event_differences_vs_V5_same_complete_ID = sum(newrisk$event[compare] != 
                  oldrisk$event[compare]), time_differences_vs_V5_same_complete_ID = sum(newrisk$exit[compare] != 
                  oldrisk$exit[compare]), event_time_differences_vs_original_incident = 0L, participants_ID_checked = TRUE, 
                A_target_or_complete_sample_changed = any(keep & !compare), exclusion_rule = "UNCHANGED_V5_INDICATOR_APPLICABILITY_AND_VALID_RECORDS_NO_121_INTERSECTION", 
                new_Y_definition = cfg$incident$definition, previous_Y_definition = ctx$endpoints[[e]]$endpoint_version)
        }
    }
    index_path <- file.path(attempt, "A_SOURCE_INDEX.json")
    oldwork <- ff_p_env("FF_V5_WORK", cfg$previous_work_paths$winpc)
    database <- file.path(oldwork, "coordinator/tasks.sqlite3")
    args <- c(file.path(root, "tools/export_a_source_index.py"), "--root", root, "--old-run", old$run, "--database", 
        database, "--output", index_path)
    audit <- file.path(private, "a_source_audit/A_SOURCE_INDEX.json")
    audit_seal <- sub("[.]json$", ".COMPLETE.json", audit)
    if (file.exists(audit_seal)) {
        accepted <- jsonlite::fromJSON(audit_seal, simplifyVector = FALSE)
        ff_p_assert(isTRUE(accepted$complete) && accepted$metadata_rows == 121L * 5L * 101L && identical(accepted$json_sha256, 
            ff_p_sha(audit)) && identical(accepted$old_input_lock_sha256, ctx$input_lock_sha256) && identical(accepted$approved_scope_sha256, 
            cfg$active_approved_csv_sha256) && identical(accepted$source_script_sha256, ff_p_sha(file.path(root, 
            "tools/export_a_source_index.py"))), "FF_PREVERIFIED_A_INDEX_BINDING_CHANGED")
        ff_p_assert(file.copy(audit, index_path, overwrite = FALSE), "FF_A_AUDIT_CHECKPOINT_COPY_FAILED")
        ff_p_assert(identical(ff_p_sha(index_path), accepted$json_sha256), "FF_A_AUDIT_CHECKPOINT_COPY_CHANGED")
    }
    else {
        meta <- system2(ff_p_env("FF_PYTHON", "python3"), vapply(args, shQuote, ""), stdout = TRUE, stderr = TRUE)
        ff_p_assert(is.null(attr(meta, "status")) || attr(meta, "status") == 0L, paste("FF_READONLY_A_INDEX_EXPORT_FAILED", 
            paste(tail(meta, 8L), collapse = "\n")))
    }
    ff_p_assert(file.exists(index_path), "FF_A_SOURCE_INDEX_NOT_CREATED")
    R_hash <- ff_p_R_hash(root)
    config_sha <- ff_p_sha(file.path(root, "config/runtime.json"))
    execution <- ff_p_env("FF_EXECUTION_SOURCE_HASH")
    ff_p_assert(grepl("^[a-f0-9]{64}$", execution), "FF_EXECUTION_SOURCE_HASH_REQUIRED")
    oldlp <- file.path(old$run, "private/V5_INPUT_LOCK.json")
    oldlock <- jsonlite::fromJSON(oldlp)
    ff_p_assert(identical(ff_p_sha(oldlp), ctx$input_lock_sha256) && identical(ff_p_sha(incident_path), incident_sha) && 
        identical(ff_p_sha(cox_file), cox_sha) && identical(ff_p_sha(coverage_file), coverage_sha), "FF_AUTHORITATIVE_SOURCE_CHANGED_DURING_PREPARATION")
    a_sha <- ff_p_sha(index_path)
    authorization <- file.path(root, "config/NATIVE_WEIGHTED_CIF_AUTHORIZATION_20261027.json")
    ff_p_assert(file.exists(authorization) && isTRUE(cfg$model$native_weighted_CIF_authorized), "FF_NATIVE_WEIGHTED_CIF_AUTHORIZATION_REQUIRED")
    authorization_sha <- ff_p_sha(authorization)
    method <- digest::digest(list(R_source_hash = R_hash, config_sha256 = config_sha, incident = incident_sha, 
        native_authorization = authorization_sha, A = "FROZEN_V5_PER_NUMBER_EMPIRICAL_SHIFT_OR_EXACT_CATEGORY_ENUMERATION", 
        weights = "ORIGINAL_SELECTION_IPW_NO_OLD_FIVE_YEAR_IPCW", tau = tau), algo = "sha256")
    frozen <- list(schema = "FF_THIN_FULLTIME_INPUTS_V1", endpoints = eps, tau_report = tau, source_incident_dates = list(baseline = baseline, 
        death = death, loss = loss, events = dates), alignment = do.call(rbind, alignment), coverage = do.call(rbind, 
        coverage), mediator_alignment = do.call(rbind, cohort_alignment), old_v5_input_lock_sha256 = ff_p_sha(oldlp), 
        old_v5_prepared_sha256 = oldlock$prepared_inputs_sha256, incident_source = incident_path, incident_sha256 = incident_sha, 
        original_cox_code_sha256 = cox_sha, A_index_sha256 = a_sha, N_parent = length(ids), parent_frame_hash = ctx$parent_frame_hash, 
        source_hash = R_hash, ff_method_hash = method, config_sha256 = config_sha, exposure_alignment = exposure_alignment)
    cache <- file.path(attempt, "FF_INPUTS.private.rds")
    saveRDS(frozen, cache, compress = FALSE)
    ff_p_assert(file.info(cache)$size < 500 * 1024^2, "FF_THIN_OUTCOME_CACHE_EXCEEDS_LOCAL_BULK_POLICY")
    v5$v5_csv(frozen$alignment, file.path(attempt, "INCIDENT_OUTCOME_ALIGNMENT.csv"))
    v5$v5_csv(frozen$mediator_alignment, file.path(attempt, "MEDIATOR_COHORT_ALIGNMENT.csv"))
    v5$v5_csv(frozen$coverage, file.path(attempt, "TAU_ADMIN_COVERAGE.csv"))
    lock <- list(complete = TRUE, models_run = 0L, run_id = cfg$run_id, prepared_inputs_sha256 = ff_p_sha(cache), 
        prepared_bytes = file.info(cache)$size, old_input_lock_sha256 = frozen$old_v5_input_lock_sha256, old_prepared_inputs_sha256 = frozen$old_v5_prepared_sha256, 
        incident_source_path = incident_path, incident_source_sha256 = incident_sha, original_cox_code_path = cox_file, 
        original_cox_code_sha256 = cox_sha, config_sha256 = config_sha, source_hash = R_hash, R_source_hash = R_hash, 
        execution_source_hash = execution, ff_method_hash = method, A_source_index_sha256 = a_sha, tau_report = tau, 
        N_parent = length(ids), N_indicators = 121L, N_endpoint_cells = 605L, max_bootstrap_id = 100L, scope_sha256 = ff_p_sha(scope_path), 
        previous_main_sha256 = ff_p_sha(main), native_weighted_authorization_sha256 = authorization_sha, exposure_alignment = exposure_alignment, 
        selection_auxiliary = "UNCHANGED_FROZEN_V5_E1_FIVE_YEAR_AUXILIARY_PREDICTORS_ONLY_NOT_NEW_Y", new_Y_authority = cfg$incident$definition, 
        old_outputs_modified = FALSE)
    v5$v5_json(lock, file.path(attempt, "FF_INPUT_LOCK.json"))
    for (name in setdiff(list.files(attempt), "FF_INPUT_LOCK.json")) v5$v5_publish(file.path(attempt, name), file.path(private, 
        name))
    v5$v5_publish(file.path(attempt, "FF_INPUT_LOCK.json"), file.path(private, "FF_INPUT_LOCK.json"))
    invisible(lock)
}

ff_load <- function(root) {
    ff_p_auth()
    root <- normalizePath(root, mustWork = TRUE)
    cfg <- ff_p_cfg(root)
    run <- ff_p_run(cfg)
    lp <- file.path(run, "private/FF_INPUT_LOCK.json")
    lock <- jsonlite::fromJSON(lp, simplifyVector = FALSE)
    ff_p_assert(isTRUE(lock$complete) && as.integer(lock$models_run) == 0L && identical(lock$R_source_hash, ff_p_R_hash(root)) && 
        identical(lock$config_sha256, ff_p_sha(file.path(root, "config/runtime.json"))), "FF_FROZEN_PREPARATION_OR_SOURCE_CHANGED")
    cache <- file.path(run, "private/FF_INPUTS.private.rds")
    ff_p_assert(identical(ff_p_sha(cache), lock$prepared_inputs_sha256), "FF_THIN_INPUT_CACHE_CHANGED")
    frozen <- readRDS(cache)
    old <- ff_p_v5(root, cfg)
    ctx <- old$ctx
    ff_p_assert(identical(ff_p_sha(file.path(old$run, "private/V5_INPUT_LOCK.json")), lock$old_input_lock_sha256) && 
        identical(ctx$parent_frame_hash, frozen$parent_frame_hash), "FF_ORIGINAL_V5_INPUT_OR_PARENT_CHANGED")
    scope <- read.csv(file.path(root, cfg$active_approved_path), colClasses = "character", check.names = FALSE)
    ff_p_assert(nrow(scope) == 121L && identical(ff_p_sha(file.path(root, cfg$active_approved_path)), lock$scope_sha256), 
        "FF_SCOPE_CHANGED")
    index_path <- file.path(run, "private/A_SOURCE_INDEX.json")
    ff_p_assert(identical(ff_p_sha(index_path), lock$A_source_index_sha256), "FF_A_REUSE_INDEX_CHANGED")
    idx <- jsonlite::fromJSON(index_path, simplifyVector = FALSE)
    ff_p_assert(identical(idx$old_input_lock_sha256, lock$old_input_lock_sha256) && length(idx$rows) == 121L * 
        5L * 101L, "FF_A_REUSE_INDEX_NOT_BOUND_TO_FROZEN_V5")
    for (i in seq_along(idx$rows)) for (name in c("path", "archive", "draw_artifact")) if (!is.null(idx$rows[[i]][[name]])) 
        idx$rows[[i]][[name]] <- ff_p_rebase(idx$rows[[i]][[name]])
    keys <- vapply(idx$rows, function(r) paste(r$indicator_id, r$bootstrap_id, r$endpoint_id, sep = "|"), "")
    ff_p_assert(!anyDuplicated(keys), "FF_DUPLICATED_A_SOURCE_CELL")
    ctx$a_reuse_index <- stats::setNames(idx$rows, keys)
    ctx$old_endpoints <- ctx$endpoints
    ctx$endpoints <- frozen$endpoints
    ctx$scope <- scope
    ctx$ff_config <- cfg
    ctx$tau_report <- frozen$tau_report
    ctx$old_run_id <- ctx$run_id
    ctx$run_id <- cfg$run_id
    ctx$old_run <- old$run
    ctx$run <- run
    ctx$ff_run <- run
    ctx$root <- root
    ctx$ff_method_hash <- lock$ff_method_hash
    ctx$source_hash <- lock$source_hash
    ctx$input_lock_sha256 <- ff_p_sha(lp)
    ctx$ff_input_lock_sha256 <- ctx$input_lock_sha256
    ctx$old_input_lock_sha256 <- lock$old_input_lock_sha256
    ctx$old_prepared_inputs_sha256 <- lock$old_prepared_inputs_sha256
    ctx$allow_missing_a_recompute <- TRUE
    ctx$ff_alignment <- frozen$alignment
    ctx$ff_coverage <- frozen$coverage
    ctx$ff_mediator_alignment <- frozen$mediator_alignment
    ctx$ff_dates <- frozen$source_incident_dates
    ctx
}

ff_draw <- function(ctx, bid) {
    ff_p_assert(is.numeric(bid) && length(bid) == 1L && bid == as.integer(bid) && bid %in% 0:100, "FF_FIXED_DRAW_ZERO_TO_100_REQUIRED")
    d <- file.path(ctx$old_run, "private/draws", sprintf("bootstrap_%04d", as.integer(bid)))
    marker <- jsonlite::fromJSON(file.path(d, "COMPLETE.json"), simplifyVector = FALSE)
    path <- file.path(d, "DRAW.private.rds")
    ff_p_assert(isTRUE(marker$complete) && as.integer(marker$bootstrap_id) == as.integer(bid) && identical(ff_p_sha(path), 
        marker$private_sha256), "FF_ORIGINAL_ACTUAL_DRAW_SEAL_CHANGED")
    saved <- readRDS(path)
    ff_p_assert(identical(saved$parent_frame_hash, ctx$parent_frame_hash) && identical(saved$cluster_order_hash, 
        ctx$cluster_order_hash) && identical(marker$parent_cluster_membership_hash, ctx$parent_cluster_membership_hash), 
        "FF_ACTUAL_DRAW_PARENT_CHANGED")
    rows <- saved$sampled_parent_rows
    ff_p_assert(is.integer(rows) && length(rows) > 0L && all(rows >= 1L & rows <= nrow(ctx$d_full)) && identical(ff_a_hash(rows), 
        marker$draw_hash) && identical(ff_a_hash(saved$cluster_multiplicity), marker$multiplicity_hash), "FF_ORIGINAL_ACTUAL_DRAW_CONTENT_CHANGED")
    list(ix = rows, multiplicity = saved$cluster_multiplicity, receipt = c(marker, list(actual_draw_private_sha256 = marker$private_sha256, 
        actual_draw_path = path, draw_source = "REUSED_RETAINED_ACTUAL_FULL_PARENT_INDICES_NO_NEW_RNG")))
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        ff_p_assert <- function(ok, text) if (!isTRUE(ok)) stop(text, call. = FALSE)
        ff_p_env <- function(name, default = "") {
            x <- Sys.getenv(name, "")
            if (nzchar(x)) x else default
        }
        ff_p_sha <- function(path) {
            ff_p_assert(file.exists(path) && !dir.exists(path), paste("FF_REQUIRED_FILE_MISSING", path))
            digest::digest(file = path, algo = "sha256")
        }
        ff_p_cfg <- function(root) jsonlite::fromJSON(file.path(root, "config/runtime.json"), simplifyVector = FALSE)
        ff_p_host <- function() ff_p_env("FF_HOST", if (Sys.info()[["sysname"]] == "Darwin") "macmini" else "winpc")
        ff_p_rebase <- function(path) {
            prefix <- birth_private_setting("PRIVATE_INPUT")
            input <- ff_p_env("FF_INPUT_ROOT", birth_private_setting("PRIVATE_BIRTH"))
            if (startsWith(path, prefix)) file.path(input, substring(path, nchar(prefix) + 1L)) else path
        }
        ff_p_run <- function(cfg, create = FALSE) {
            host <- ff_p_host()
            run <- ff_p_env("FF_RUN_DIR", cfg$run_paths[[host]])
            ff_p_assert(grepl("/UKB_results/birth/border_fullfollowup_b100_20261027/run_[A-Za-z0-9_.-]+$", run) && 
                !grepl("(^|/)[.][.](/|$)", run) && (startsWith(run, birth_private_setting("PRIVATE_INPUT")) || 
                startsWith(run, birth_private_setting("PRIVATE_INPUT"))), "FF_ISOLATED_G_VERSIONED_RUN_REQUIRED")
            if (create && !dir.exists(run)) ff_p_assert(dir.create(run, recursive = TRUE), "FF_NEW_RUN_CREATE_FAILED")
            normalizePath(run, mustWork = TRUE)
        }
        ff_p_R_hash <- function(root) {
            files <- sort(list.files(file.path(root, "f"), pattern = "[.]R$", full.names = TRUE))
            files <- c(files, sort(list.files(file.path(root, "src"), pattern = "[.]c$", full.names = TRUE)))
            digest::digest(paste(basename(files), vapply(files, ff_p_sha, ""), sep = ":", collapse = "\n"), algo = "sha256", 
                serialize = FALSE)
        }
        ff_p_auth <- function() ff_p_assert(identical(ff_p_env("FF_REAL_EXECUTION_AUTHORIZED"), "YES_FOR_THIS_RUN"), 
            "FF_EXPLICIT_MANUAL_REAL_AUTHORIZATION_REQUIRED")
        ff_p_exposure_alignment <- function(frozen, original) {
            ff_p_assert(is.numeric(frozen) && is.numeric(original) && length(frozen) == length(original) && length(frozen) > 
                0L && all(is.finite(frozen)) && all(is.finite(original)) && all(frozen >= 0 & frozen <= 1) && all(original >= 
                0 & original <= 1), "FF_ORIGINAL_CONTINUOUS_B_INVALID")
            difference <- max(abs(frozen - original))
            tolerance <- 8 * .Machine$double.eps
            ff_p_assert(difference <= tolerance, "FF_ORIGINAL_CONTINUOUS_B_CONFLICT")
            list(max_abs_difference = difference, nonidentical_double_values = sum(frozen != original), tolerance = tolerance, 
                tolerance_basis = "EIGHT_IEEE_DOUBLE_MACHINE_EPS_NO_VALUE_CHANGE", frozen_V5_B_preserved = TRUE)
        }
        ff_p_v5 <- function(root, cfg) {
            ff_p_auth()
            host <- ff_p_host()
            oldroot <- ff_p_env("FF_V5_ROOT", cfg$previous_code_paths[[host]])
            oldrun <- ff_p_env("FF_V5_RUN", cfg$previous_run_paths[[host]])
            adapter <- file.path(oldroot, "f/1.adapter.R")
            ff_p_assert(file.exists(adapter), paste("FF_EXISTING_V5_ADAPTER_MISSING", adapter))
            v5 <- new.env(parent = globalenv())
            sys.source(adapter, envir = v5)
            for (name in c("v5_bind", "v5_check", "v5_hash", "v5_sha", "v5_stage", "v5_utc")) assign(name, v5[[name]], 
                envir = environment(ff_p_v5))
            before <- Sys.getenv("V5_REAL_EXECUTION_AUTHORIZED", unset = NA_character_)
            on.exit(if (is.na(before)) Sys.unsetenv("V5_REAL_EXECUTION_AUTHORIZED") else Sys.setenv(V5_REAL_EXECUTION_AUTHORIZED = before), 
                add = TRUE)
            Sys.setenv(V5_REAL_EXECUTION_AUTHORIZED = "YES_FOR_THIS_RUN")
            ctx <- v5$v5_load(oldroot, oldrun)
            list(ctx = ctx, api = v5, root = oldroot, run = oldrun)
        }
        ff_p_endpoint <- function(parent_ids, original, risk, baseline, death, loss, admin) {
            j <- match(parent_ids, as.character(risk$eid))
            eligible <- !is.na(j)
            n <- length(parent_ids)
            r <- data.frame(eid = parent_ids, entry = rep(0, n), exit = rep(NA_real_, n), event = rep(NA_integer_, 
                n), death_event = rep(NA_integer_, n))
            r$entry[eligible] <- risk$.entry_time[j[eligible]]
            r$exit[eligible] <- risk$.exit_time[j[eligible]]
            r$event[eligible] <- risk$.event[j[eligible]]
            death_years <- as.numeric(death - baseline)/365.25
            loss_years <- as.numeric(loss - baseline)/365.25
            admin_years <- as.numeric(as.Date(admin) - baseline)/365.25
            r$death_event[eligible] <- as.integer(r$event[eligible] == 0L & !is.na(death_years[eligible]) & death_years[eligible] > 
                0 & death_years[eligible] == r$exit[eligible])
            reason <- rep("NOT_IN_ORIGINAL_RISKSET", n)
            reason[eligible & r$event == 1L] <- "TARGET_EVENT"
            reason[eligible & r$death_event == 1L] <- "DEATH_BEFORE_TARGET_EVENT"
            nd <- eligible & r$event == 0L & r$death_event == 0L
            nd[is.na(nd)] <- FALSE
            losscap <- !is.na(loss_years) & loss_years == r$exit
            admincap <- admin_years == r$exit
            reason[nd & losscap & !admincap] <- "RECORDED_NONDEATH_LOSS"
            reason[nd & !losscap & admincap] <- "ADMINISTRATIVE_CUTOFF"
            reason[nd & losscap & admincap] <- "LOSS_AND_ADMIN_SAME_DATE"
            ff_p_assert(all(reason[eligible] != "NOT_IN_ORIGINAL_RISKSET"), "FF_ORIGINAL_CENSOR_REASON_NOT_RECONSTRUCTED")
            list(risk = r, eligibility = eligible, censor_reason = reason)
        }
        ff_prepare <- function(root) {
            ff_p_auth()
            ff_p_assert(identical(ff_p_host(), "winpc") && Sys.info()[["sysname"]] != "Darwin", "FF_PARTICIPANT_PREPARATION_WINPC_ONLY")
            root <- normalizePath(root, mustWork = TRUE)
            cfg <- ff_p_cfg(root)
            run <- ff_p_run(cfg, TRUE)
            ff_p_assert(identical(as.integer(unlist(cfg$bootstrap_ids)), 0:100) && as.integer(cfg$max_bootstrap_id) == 
                100L, "FF_FIXED_B0_PLUS_1_TO_100_REQUIRED")
            private <- file.path(run, "private")
            if (!dir.exists(private)) ff_p_assert(dir.create(private), "FF_PRIVATE_CREATE_FAILED")
            ff_p_assert(!file.exists(file.path(private, "FF_INPUT_LOCK.json")), "FF_PRESERVE_EXISTING_PREPARATION")
            attempt <- file.path(private, paste0("preparation_attempt_", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), 
                "_", Sys.getpid()))
            ff_p_assert(!dir.exists(attempt) && dir.create(attempt), "FF_UNIQUE_PREPARATION_ATTEMPT_REQUIRED")
            old <- ff_p_v5(root, cfg)
            ctx <- old$ctx
            v5 <- old$api
            scope_path <- file.path(root, cfg$active_approved_path)
            ff_p_assert(identical(ff_p_sha(scope_path), cfg$active_approved_csv_sha256), "FF_121_REGISTRY_SHA_CHANGED")
            scope <- read.csv(scope_path, colClasses = "character", check.names = FALSE)
            ff_p_assert(nrow(scope) == 121L && !anyDuplicated(scope$indicator_id) && !any(scope$category == "CMR") && 
                all(scope$indicator_id %in% ctx$scope$indicator_id), "FF_EXACT_EXISTING_121_NON_CMR_REQUIRED")
            main <- ff_p_rebase(cfg$previous_main$path)
            ff_p_assert(identical(ff_p_sha(main), cfg$previous_main$sha256), "FF_OLD_MAIN_SHA_CHANGED")
            incident_path <- ff_p_env("FF_INCIDENT_PREPARED", ff_p_rebase(cfg$incident$prepared))
            incident_sha <- ff_p_sha(incident_path)
            inc <- readRDS(incident_path)
            while (is.list(inc) && !is.null(inc$value)) inc <- inc$value
            ff_p_assert(all(c("data", "events", "baseline") %in% names(inc)), "FF_ORIGINAL_INCIDENT_SCHEMA_REQUIRED")
            code_root <- ff_p_env("FF_INCIDENT_CODE_ROOT", cfg$incident$code_root)
            cox_file <- file.path(code_root, "vendor/incident/cox.R")
            coverage_file <- file.path(code_root, "vendor/incident/provenance.R")
            cox_sha <- ff_p_sha(cox_file)
            coverage_sha <- ff_p_sha(coverage_file)
            ff_p_assert(grepl("attr\\(out, \"cutoff\"\\) <- \"2022-05-31\"", paste(readLines(coverage_file, warn = FALSE), 
                collapse = "\n")), "FF_FROZEN_UNIFORM_COVERAGE_CONTRACT_CHANGED")
            engine <- new.env(parent = globalenv())
            sys.source(cox_file, envir = engine)
            ids <- as.character(ctx$d_full$eid)
            match_source <- match(ids, as.character(inc$data$eid))
            ff_p_assert(!anyNA(match_source) && !anyDuplicated(as.character(inc$data$eid)), "FF_ORIGINAL_INCIDENT_NOT_COVERING_FROZEN_PARENT")
            sd <- inc$data[match_source, , drop = FALSE]
            for (name in c("age", "sex", "ethnicity", "center", "assessment_year", "S")) ff_p_assert(identical(as.character(ctx$d_full[[name]]), 
                as.character(sd[[name]])), paste("FF_ORIGINAL_COVARIATE_CONFLICT", name))
            exposure_alignment <- ff_p_exposure_alignment(as.numeric(ctx$d_full$B), as.numeric(sd$border_relative))
            admin <- cfg$incident$admin
            baseline <- as.Date(sd$date_attend)
            death <- as.Date(sd$date_death_uncapped)
            loss <- as.Date(sd$date_lost_registry)
            potential <- as.numeric(as.Date(admin) - baseline)
            ff_p_assert(all(is.finite(potential) & potential > 0), "FF_POTENTIAL_ADMIN_COVERAGE_MISSING")
            tau <- floor(min(potential/365.25))
            ff_p_assert(tau > 0, "FF_NO_POSITIVE_COMMON_ADMIN_REPORT_TIME")
            eps <- alignment <- coverage <- dates <- list()
            for (e in as.character(unlist(cfg$endpoints))) {
                ff_p_assert(e %in% names(inc$events) && e %in% names(inc$baseline), "FF_ORIGINAL_ENDPOINT_MISSING")
                raw <- engine$ip_riskset(inc$data, inc$events[[e]], inc$baseline[[e]], "primary", admin)$data
                z <- ff_p_endpoint(ids, inc, raw, baseline, death, loss, admin)
                oldep <- ctx$endpoints[[e]]
                oldrisk <- oldep$risk
                same <- z$eligibility & oldep$eligibility
                ok <- all(!oldep$eligibility | z$eligibility)
                rh <- digest::digest(z$risk, algo = "sha256")
                eps[[e]] <- list(status = if (ok) "AVAILABLE_FROZEN" else "BLOCKED_ALIGNMENT", reason = if (ok) "" else "V5_PARTICIPANTS_OUTSIDE_ORIGINAL_COMPLETED_INCIDENT_RISKSET_NO_SILENT_SAMPLE_CHANGE", 
                  eligibility = oldep$eligibility, incident_eligibility = z$eligibility, risk = z$risk, endpoint_version = cfg$incident$definition, 
                  risk_hash = rh, admin = admin, censor_reason = z$censor_reason, source = incident_path, source_sha256 = incident_sha)
                alignment[[e]] <- data.frame(endpoint_id = e, source_file = incident_path, source_sha256 = incident_sha, 
                  frozen_incident_version = cfg$incident$definition, previous_V5_version = oldep$endpoint_version, 
                  N_original_incident_riskset = nrow(raw), N_parent = length(ids), N_parent_ID_matched = sum(!is.na(match_source)), 
                  exposure_max_abs_difference = exposure_alignment$max_abs_difference, exposure_numeric_tolerance = exposure_alignment$tolerance, 
                  frozen_V5_exposure_preserved = TRUE, N_V5_parent_eligible = sum(oldep$eligibility), N_original_eligible_same_parent = sum(z$eligibility), 
                  N_risk_qualification_conflict = sum(oldep$eligibility & !z$eligibility), N_same_ID_time_event_compared = sum(same), 
                  event_differences_vs_V5 = sum(z$risk$event[same] != oldrisk$event[same]), time_differences_vs_V5 = sum(z$risk$exit[same] != 
                    oldrisk$exit[same]), original_events_same_ID = sum(z$risk$event[same]), V5_events_same_ID = sum(oldrisk$event[same]), 
                  original_competing_deaths_same_ID = sum(z$risk$death_event[same]), original_nondeath_censor_same_ID = sum(z$censor_reason[same] %in% 
                    c("RECORDED_NONDEATH_LOSS", "ADMINISTRATIVE_CUTOFF", "LOSS_AND_ADMIN_SAME_DATE")), event_difference_vs_authoritative_incident = 0L, 
                  time_difference_vs_authoritative_incident = 0L, dates_rule = "ORIGINAL_EVENT_ON_CENSOR_DATE_WINS;DEATH_ONLY_ADDED_WHERE_ORIGINAL_EVENT_ZERO", 
                  time_unit = "ELAPSED_DAYS_DIVIDED_BY_365.25_YEARS", time_zero = "FROZEN_ORIGINAL_BASELINE_DATE", 
                  same_A_population_allowed = ok, Y_version_change_disclosed = !identical(oldep$endpoint_version, 
                    cfg$incident$definition), mediator_missingness = "PER_INDICATOR_FROZEN_COMPLETE_RECORD_REQUIREMENTS_NOT_COMMON_121_INTERSECTION")
                coverage[[e]] <- data.frame(endpoint_id = e, coverage_end = admin, min_potential_admin_days = min(potential), 
                  max_potential_admin_days = max(potential), raw_common_potential_years = min(potential/365.25), 
                  tau_report = tau, source_contract = coverage_file, source_contract_sha256 = coverage_sha, coverage = "FROZEN_CONSERVATIVE_INSTITUTIONAL_UNIFORM_HOSPITAL_CUTOFF_NOT_INDIVIDUAL_LINKAGE_CERTIFICATION", 
                  limitation = "ORIGINAL_UKB_EXTRACTION_RELEASE_ID_NOT_FOUND;EXACT_LOCAL_SNAPSHOT_HASHES_USED", 
                  uses_actual_event_death_loss_minimum = FALSE)
                dates[[e]] <- as.Date(inc$events[[e]][match_source])
            }
            cohort_alignment <- list()
            for (id in scope$indicator_id) {
                ca <- ctx$features[[id]]
                base <- unique(c(as.character(unlist(ctx$cfg$base_covariates)), ca$technical_covariates))
                cc <- stats::complete.cases(ctx$d_full[, unique(c(base, ca$column)), drop = FALSE])
                for (e in as.character(unlist(cfg$endpoints))) {
                  target <- ctx$eligibility[[id]] & ctx$endpoints[[e]]$eligibility
                  keep <- target & cc
                  compare <- keep & eps[[e]]$incident_eligibility
                  oldrisk <- ctx$endpoints[[e]]$risk
                  newrisk <- eps[[e]]$risk
                  cohort_alignment[[paste(id, e)]] <- data.frame(indicator_id = id, endpoint_id = e, N_frozen_parent = length(ids), 
                    N_frozen_indicator_endpoint_target = sum(target), N_original_V5_complete_records = sum(keep), 
                    N_fullfollowup_same_ID_complete_records = sum(compare), N_indicator_or_applicability_excluded = sum(!target), 
                    N_invalid_mediator_or_covariate_records = sum(target & !cc), N_original_incident_ineligible = sum(keep & 
                      !eps[[e]]$incident_eligibility), source_incident_events_in_complete_records = sum(newrisk$event[compare]), 
                    source_competing_deaths_in_complete_records = sum(newrisk$death_event[compare]), event_differences_vs_V5_same_complete_ID = sum(newrisk$event[compare] != 
                      oldrisk$event[compare]), time_differences_vs_V5_same_complete_ID = sum(newrisk$exit[compare] != 
                      oldrisk$exit[compare]), event_time_differences_vs_original_incident = 0L, participants_ID_checked = TRUE, 
                    A_target_or_complete_sample_changed = any(keep & !compare), exclusion_rule = "UNCHANGED_V5_INDICATOR_APPLICABILITY_AND_VALID_RECORDS_NO_121_INTERSECTION", 
                    new_Y_definition = cfg$incident$definition, previous_Y_definition = ctx$endpoints[[e]]$endpoint_version)
                }
            }
            index_path <- file.path(attempt, "A_SOURCE_INDEX.json")
            oldwork <- ff_p_env("FF_V5_WORK", cfg$previous_work_paths$winpc)
            database <- file.path(oldwork, "coordinator/tasks.sqlite3")
            args <- c(file.path(root, "tools/export_a_source_index.py"), "--root", root, "--old-run", old$run, 
                "--database", database, "--output", index_path)
            audit <- file.path(private, "a_source_audit/A_SOURCE_INDEX.json")
            audit_seal <- sub("[.]json$", ".COMPLETE.json", audit)
            if (file.exists(audit_seal)) {
                accepted <- jsonlite::fromJSON(audit_seal, simplifyVector = FALSE)
                ff_p_assert(isTRUE(accepted$complete) && accepted$metadata_rows == 121L * 5L * 101L && identical(accepted$json_sha256, 
                  ff_p_sha(audit)) && identical(accepted$old_input_lock_sha256, ctx$input_lock_sha256) && identical(accepted$approved_scope_sha256, 
                  cfg$active_approved_csv_sha256) && identical(accepted$source_script_sha256, ff_p_sha(file.path(root, 
                  "tools/export_a_source_index.py"))), "FF_PREVERIFIED_A_INDEX_BINDING_CHANGED")
                ff_p_assert(file.copy(audit, index_path, overwrite = FALSE), "FF_A_AUDIT_CHECKPOINT_COPY_FAILED")
                ff_p_assert(identical(ff_p_sha(index_path), accepted$json_sha256), "FF_A_AUDIT_CHECKPOINT_COPY_CHANGED")
            } else {
                meta <- system2(ff_p_env("FF_PYTHON", "python3"), vapply(args, shQuote, ""), stdout = TRUE, stderr = TRUE)
                ff_p_assert(is.null(attr(meta, "status")) || attr(meta, "status") == 0L, paste("FF_READONLY_A_INDEX_EXPORT_FAILED", 
                  paste(tail(meta, 8L), collapse = "\n")))
            }
            ff_p_assert(file.exists(index_path), "FF_A_SOURCE_INDEX_NOT_CREATED")
            R_hash <- ff_p_R_hash(root)
            config_sha <- ff_p_sha(file.path(root, "config/runtime.json"))
            execution <- ff_p_env("FF_EXECUTION_SOURCE_HASH")
            ff_p_assert(grepl("^[a-f0-9]{64}$", execution), "FF_EXECUTION_SOURCE_HASH_REQUIRED")
            oldlp <- file.path(old$run, "private/V5_INPUT_LOCK.json")
            oldlock <- jsonlite::fromJSON(oldlp)
            ff_p_assert(identical(ff_p_sha(oldlp), ctx$input_lock_sha256) && identical(ff_p_sha(incident_path), 
                incident_sha) && identical(ff_p_sha(cox_file), cox_sha) && identical(ff_p_sha(coverage_file), coverage_sha), 
                "FF_AUTHORITATIVE_SOURCE_CHANGED_DURING_PREPARATION")
            a_sha <- ff_p_sha(index_path)
            authorization <- file.path(root, "config/NATIVE_WEIGHTED_CIF_AUTHORIZATION_20261027.json")
            ff_p_assert(file.exists(authorization) && isTRUE(cfg$model$native_weighted_CIF_authorized), "FF_NATIVE_WEIGHTED_CIF_AUTHORIZATION_REQUIRED")
            authorization_sha <- ff_p_sha(authorization)
            method <- digest::digest(list(R_source_hash = R_hash, config_sha256 = config_sha, incident = incident_sha, 
                native_authorization = authorization_sha, A = "FROZEN_V5_PER_NUMBER_EMPIRICAL_SHIFT_OR_EXACT_CATEGORY_ENUMERATION", 
                weights = "ORIGINAL_SELECTION_IPW_NO_OLD_FIVE_YEAR_IPCW", tau = tau), algo = "sha256")
            frozen <- list(schema = "FF_THIN_FULLTIME_INPUTS_V1", endpoints = eps, tau_report = tau, source_incident_dates = list(baseline = baseline, 
                death = death, loss = loss, events = dates), alignment = do.call(rbind, alignment), coverage = do.call(rbind, 
                coverage), mediator_alignment = do.call(rbind, cohort_alignment), old_v5_input_lock_sha256 = ff_p_sha(oldlp), 
                old_v5_prepared_sha256 = oldlock$prepared_inputs_sha256, incident_source = incident_path, incident_sha256 = incident_sha, 
                original_cox_code_sha256 = cox_sha, A_index_sha256 = a_sha, N_parent = length(ids), parent_frame_hash = ctx$parent_frame_hash, 
                source_hash = R_hash, ff_method_hash = method, config_sha256 = config_sha, exposure_alignment = exposure_alignment)
            cache <- file.path(attempt, "FF_INPUTS.private.rds")
            saveRDS(frozen, cache, compress = FALSE)
            ff_p_assert(file.info(cache)$size < 500 * 1024^2, "FF_THIN_OUTCOME_CACHE_EXCEEDS_LOCAL_BULK_POLICY")
            v5$v5_csv(frozen$alignment, file.path(attempt, "INCIDENT_OUTCOME_ALIGNMENT.csv"))
            v5$v5_csv(frozen$mediator_alignment, file.path(attempt, "MEDIATOR_COHORT_ALIGNMENT.csv"))
            v5$v5_csv(frozen$coverage, file.path(attempt, "TAU_ADMIN_COVERAGE.csv"))
            lock <- list(complete = TRUE, models_run = 0L, run_id = cfg$run_id, prepared_inputs_sha256 = ff_p_sha(cache), 
                prepared_bytes = file.info(cache)$size, old_input_lock_sha256 = frozen$old_v5_input_lock_sha256, 
                old_prepared_inputs_sha256 = frozen$old_v5_prepared_sha256, incident_source_path = incident_path, 
                incident_source_sha256 = incident_sha, original_cox_code_path = cox_file, original_cox_code_sha256 = cox_sha, 
                config_sha256 = config_sha, source_hash = R_hash, R_source_hash = R_hash, execution_source_hash = execution, 
                ff_method_hash = method, A_source_index_sha256 = a_sha, tau_report = tau, N_parent = length(ids), 
                N_indicators = 121L, N_endpoint_cells = 605L, max_bootstrap_id = 100L, scope_sha256 = ff_p_sha(scope_path), 
                previous_main_sha256 = ff_p_sha(main), native_weighted_authorization_sha256 = authorization_sha, 
                exposure_alignment = exposure_alignment, selection_auxiliary = "UNCHANGED_FROZEN_V5_E1_FIVE_YEAR_AUXILIARY_PREDICTORS_ONLY_NOT_NEW_Y", 
                new_Y_authority = cfg$incident$definition, old_outputs_modified = FALSE)
            v5$v5_json(lock, file.path(attempt, "FF_INPUT_LOCK.json"))
            for (name in setdiff(list.files(attempt), "FF_INPUT_LOCK.json")) v5$v5_publish(file.path(attempt, name), 
                file.path(private, name))
            v5$v5_publish(file.path(attempt, "FF_INPUT_LOCK.json"), file.path(private, "FF_INPUT_LOCK.json"))
            invisible(lock)
        }
        ff_load <- function(root) {
            ff_p_auth()
            root <- normalizePath(root, mustWork = TRUE)
            cfg <- ff_p_cfg(root)
            run <- ff_p_run(cfg)
            lp <- file.path(run, "private/FF_INPUT_LOCK.json")
            lock <- jsonlite::fromJSON(lp, simplifyVector = FALSE)
            ff_p_assert(isTRUE(lock$complete) && as.integer(lock$models_run) == 0L && identical(lock$R_source_hash, 
                ff_p_R_hash(root)) && identical(lock$config_sha256, ff_p_sha(file.path(root, "config/runtime.json"))), 
                "FF_FROZEN_PREPARATION_OR_SOURCE_CHANGED")
            cache <- file.path(run, "private/FF_INPUTS.private.rds")
            ff_p_assert(identical(ff_p_sha(cache), lock$prepared_inputs_sha256), "FF_THIN_INPUT_CACHE_CHANGED")
            frozen <- readRDS(cache)
            old <- ff_p_v5(root, cfg)
            ctx <- old$ctx
            ff_p_assert(identical(ff_p_sha(file.path(old$run, "private/V5_INPUT_LOCK.json")), lock$old_input_lock_sha256) && 
                identical(ctx$parent_frame_hash, frozen$parent_frame_hash), "FF_ORIGINAL_V5_INPUT_OR_PARENT_CHANGED")
            scope <- read.csv(file.path(root, cfg$active_approved_path), colClasses = "character", check.names = FALSE)
            ff_p_assert(nrow(scope) == 121L && identical(ff_p_sha(file.path(root, cfg$active_approved_path)), lock$scope_sha256), 
                "FF_SCOPE_CHANGED")
            index_path <- file.path(run, "private/A_SOURCE_INDEX.json")
            ff_p_assert(identical(ff_p_sha(index_path), lock$A_source_index_sha256), "FF_A_REUSE_INDEX_CHANGED")
            idx <- jsonlite::fromJSON(index_path, simplifyVector = FALSE)
            ff_p_assert(identical(idx$old_input_lock_sha256, lock$old_input_lock_sha256) && length(idx$rows) == 
                121L * 5L * 101L, "FF_A_REUSE_INDEX_NOT_BOUND_TO_FROZEN_V5")
            for (i in seq_along(idx$rows)) for (name in c("path", "archive", "draw_artifact")) if (!is.null(idx$rows[[i]][[name]])) idx$rows[[i]][[name]] <- ff_p_rebase(idx$rows[[i]][[name]])
            keys <- vapply(idx$rows, function(r) paste(r$indicator_id, r$bootstrap_id, r$endpoint_id, sep = "|"), 
                "")
            ff_p_assert(!anyDuplicated(keys), "FF_DUPLICATED_A_SOURCE_CELL")
            ctx$a_reuse_index <- stats::setNames(idx$rows, keys)
            ctx$old_endpoints <- ctx$endpoints
            ctx$endpoints <- frozen$endpoints
            ctx$scope <- scope
            ctx$ff_config <- cfg
            ctx$tau_report <- frozen$tau_report
            ctx$old_run_id <- ctx$run_id
            ctx$run_id <- cfg$run_id
            ctx$old_run <- old$run
            ctx$run <- run
            ctx$ff_run <- run
            ctx$root <- root
            ctx$ff_method_hash <- lock$ff_method_hash
            ctx$source_hash <- lock$source_hash
            ctx$input_lock_sha256 <- ff_p_sha(lp)
            ctx$ff_input_lock_sha256 <- ctx$input_lock_sha256
            ctx$old_input_lock_sha256 <- lock$old_input_lock_sha256
            ctx$old_prepared_inputs_sha256 <- lock$old_prepared_inputs_sha256
            ctx$allow_missing_a_recompute <- TRUE
            ctx$ff_alignment <- frozen$alignment
            ctx$ff_coverage <- frozen$coverage
            ctx$ff_mediator_alignment <- frozen$mediator_alignment
            ctx$ff_dates <- frozen$source_incident_dates
            ctx
        }
        ff_draw <- function(ctx, bid) {
            ff_p_assert(is.numeric(bid) && length(bid) == 1L && bid == as.integer(bid) && bid %in% 0:100, "FF_FIXED_DRAW_ZERO_TO_100_REQUIRED")
            d <- file.path(ctx$old_run, "private/draws", sprintf("bootstrap_%04d", as.integer(bid)))
            marker <- jsonlite::fromJSON(file.path(d, "COMPLETE.json"), simplifyVector = FALSE)
            path <- file.path(d, "DRAW.private.rds")
            ff_p_assert(isTRUE(marker$complete) && as.integer(marker$bootstrap_id) == as.integer(bid) && identical(ff_p_sha(path), 
                marker$private_sha256), "FF_ORIGINAL_ACTUAL_DRAW_SEAL_CHANGED")
            saved <- readRDS(path)
            ff_p_assert(identical(saved$parent_frame_hash, ctx$parent_frame_hash) && identical(saved$cluster_order_hash, 
                ctx$cluster_order_hash) && identical(marker$parent_cluster_membership_hash, ctx$parent_cluster_membership_hash), 
                "FF_ACTUAL_DRAW_PARENT_CHANGED")
            rows <- saved$sampled_parent_rows
            ff_p_assert(is.integer(rows) && length(rows) > 0L && all(rows >= 1L & rows <= nrow(ctx$d_full)) && 
                identical(ff_a_hash(rows), marker$draw_hash) && identical(ff_a_hash(saved$cluster_multiplicity), 
                marker$multiplicity_hash), "FF_ORIGINAL_ACTUAL_DRAW_CONTENT_CHANGED")
            list(ix = rows, multiplicity = saved$cluster_multiplicity, receipt = c(marker, list(actual_draw_private_sha256 = marker$private_sha256, 
                actual_draw_path = path, draw_source = "REUSED_RETAINED_ACTUAL_FULL_PARENT_INDICES_NO_NEW_RNG")))
        }
    }), envir = context)
}

