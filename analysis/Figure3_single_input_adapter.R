# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

v5_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

v5_env <- function(k, default = "") {
    x <- Sys.getenv(k, "")
    if (nzchar(x)) 
        x
    else default
}

v5_sha <- function(p) {
    v5_assert(file.exists(p) && !dir.exists(p), paste("REQUIRED_FILE_MISSING", p))
    digest::digest(file = p, algo = "sha256")
}

v5_hash <- function(x) digest::digest(paste0(nchar(as.character(x)), ":", as.character(x), "\n", collapse = ""), 
    algo = "sha256", serialize = FALSE)

v5_check <- function(p, h) {
    v5_assert(is.character(h) && length(h) == 1L && grepl("^[a-f0-9]{64}$", h), "EXACT_SOURCE_SHA256_REQUIRED")
    v5_assert(identical(v5_sha(p), h), paste("FROZEN_SOURCE_SHA256_MISMATCH", p))
}

v5_bind <- function(xs) {
    xs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, xs)
    if (!length(xs)) 
        return(data.frame())
    ns <- unique(unlist(lapply(xs, names), use.names = FALSE))
    do.call(rbind, lapply(xs, function(x) {
        for (k in setdiff(ns, names(x))) x[[k]] <- NA
        x[, ns, drop = FALSE]
    }))
}

v5_json <- function(x, p) {
    v5_assert(!file.exists(p) && !dir.exists(p), paste("REFUSE_OVERWRITE", p))
    jsonlite::write_json(x, p, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null", digits = 17)
}

v5_csv <- function(x, p) {
    v5_assert(!file.exists(p) && !dir.exists(p), paste("REFUSE_OVERWRITE", p))
    write.csv(x, p, row.names = FALSE, na = "NA")
}

v5_publish <- function(source, dest) {
    h <- v5_sha(source)
    if (file.exists(dest)) {
        v5_check(dest, h)
        return(invisible(h))
    }
    v5_assert(!dir.exists(dest), "PUBLICATION_DESTINATION_IS_DIRECTORY")
    linked <- suppressWarnings(file.link(source, dest))
    if (!linked) 
        v5_assert(file.copy(source, dest, overwrite = FALSE, copy.date = TRUE), "ADDITIVE_PUBLICATION_FAILED_STAGING_RETAINED")
    v5_check(dest, h)
    invisible(h)
}

v5_deps <- function(fit = FALSE) {
    for (p in c("jsonlite", "digest", "survival", if (fit) c("brglm2", "nnet"))) v5_assert(requireNamespace(p, 
        quietly = TRUE), paste("DEPENDENCY_MISSING", p))
    if (fit) 
        v5_assert(as.character(packageVersion("brglm2")) == "1.1.1", "PINNED_BRGLM2_1_1_1_REQUIRED")
}

v5_id <- function(id) {
    v5_assert(is.numeric(id) && length(id) == 1L && is.finite(id) && id == as.integer(id) && id %in% 0:100, "FIXED_V5_ID_ZERO_TO_100_REQUIRED")
    as.integer(id)
}

v5_utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")

v5_scope <- function(root) {
    p <- file.path(root, "input/APPROVED_165_ORIGINAL.csv")
    d <- read.csv(p, colClasses = "character", check.names = FALSE, na.strings = "", fileEncoding = "UTF-8")
    v5_assert(nrow(d) == 165L && !anyNA(d$indicator_id) && !anyDuplicated(d$indicator_id) && sum(d$type == "binary") == 
        13L && sum(d$type == "categorical") == 15L && sum(d$type == "continuous") == 137L && sum(d$category == 
        "CMR") == 44L, "EXACT_APPROVED_165_13_15_137_44_REQUIRED")
    cfg <- jsonlite::fromJSON(file.path(root, "config/runtime.json"), simplifyVector = FALSE)
    v5_assert(identical(v5_sha(p), cfg$approved_csv_sha256) && cfg$max_bootstrap_id == 100L && !isTRUE(cfg$stage2_approved), 
        "V5_SCOPE_HASH_OR_AUTHORIZATION_CHANGED")
    d
}

v5_sources <- function(root) {
    v5_deps()
    v9 <- normalizePath(v5_env("B3_V9_ROOT", dirname(dirname(root))), mustWork = TRUE)
    ca <- new.env(parent = globalenv())
    catroot <- file.path(v9, "extensions/categorical_cpu_gpu_v2")
    sys.source(file.path(catroot, "R/reference_adapter.R"), envir = ca)
    sys.source(file.path(catroot, "R/prepare_and_cache.R"), envir = ca)
    loaded <- ca$c40_sources(catroot)
    b18root <- v5_env("V5_B18_ROOT", file.path(dirname(v9), "border_gcomp_18x5_completion_20261024_v1"))
    sys.source(file.path(b18root, "f/1.adapter.R"), envir = ca)
    list(engine = loaded$engine, discrete = ca, source_sha256 = c(loaded$source_sha256, categorical_adapter = v5_sha(file.path(catroot, 
        "R/reference_adapter.R")), categorical_preparer = v5_sha(file.path(catroot, "R/prepare_and_cache.R")), 
        b18_adapter = v5_sha(file.path(b18root, "f/1.adapter.R"))), v9 = v9)
}

v5_run <- function(run, create = FALSE) {
    p <- gsub(birth_private_setting("PRIVATE_INPUT"), "/", run)
    v5_assert(nzchar(p) && !grepl("(^|/)[.][.](/|$)", p) && grepl("/border_165_typewise_b100_v5/run_[A-Za-z0-9_.-]+$", 
        p), "ISOLATED_V5_VERSIONED_RUN_REQUIRED")
    v5_assert(startsWith(p, birth_private_setting("PRIVATE_INPUT")) || grepl("^[Gg]:/UKB_results/birth/", p) || 
        grepl("PRIVATE_CONTEXT_REQUIRED", p), "G_RESULT_ROOT_REQUIRED_NO_MAC_LOCAL_CACHE")
    if (create && !dir.exists(run)) 
        v5_assert(dir.create(run, recursive = TRUE, showWarnings = FALSE), "RUN_CREATE_FAILED")
    normalizePath(run, mustWork = TRUE)
}

v5_rebase <- function(p) {
    r <- v5_env("V5_INPUT_ROOT")
    if (nzchar(r) && startsWith(p, birth_private_setting("PRIVATE_INPUT"))) 
        file.path(r, substring(p, nchar(birth_private_setting("PRIVATE_INPUT")) + 1L))
    else p
}

v5_read_source <- function(desc) {
    v5_assert(is.list(desc) && !is.null(desc$path) && !is.null(desc$sha256), "SOURCE_DESCRIPTOR_PATH_SHA_REQUIRED")
    p <- v5_rebase(desc$path)
    v5_check(p, desc$sha256)
    z <- if (identical(desc$format, "json")) 
        jsonlite::fromJSON(p, simplifyVector = FALSE)
    else if (identical(desc$format, "csv")) {
        v5_assert(requireNamespace("data.table", quietly = TRUE), "DATA_TABLE_REQUIRED")
        as.data.frame(data.table::fread(p, select = if (is.null(desc$columns)) 
            NULL
        else unlist(desc$columns), showProgress = FALSE))
    }
    else readRDS(p)
    if (isTRUE(desc$unwrap)) 
        while (is.list(z) && !is.null(z$value)) {
            if (!is.null(z$value_hash)) 
                v5_assert(identical(digest::digest(z$value, algo = "sha256"), z$value_hash), "WRAPPED_SOURCE_CONTENT_CHANGED")
            z <- z$value
        }
    if (!is.null(desc$data_path)) 
        for (k in unlist(desc$data_path, use.names = FALSE)) {
            v5_assert(!is.null(z[[k]]), paste("SOURCE_DATA_PATH_MISSING", k))
            z <- z[[k]]
        }
    v5_check(p, desc$sha256)
    z
}

v5_join <- function(data, ids, require_all = TRUE) {
    v5_assert(is.data.frame(data) && "eid" %in% names(data), "SOURCE_DATA_EID_REQUIRED")
    sid <- as.character(data$eid)
    v5_assert(!anyNA(sid) && !anyDuplicated(sid), "SOURCE_IDS_NOT_UNIQUE")
    j <- match(ids, sid)
    if (require_all) 
        v5_assert(!anyNA(j), "SOURCE_NOT_COVERING_FROZEN_PARENT")
    j
}

v5_vendor_features <- function(spec) {
    v5_assert(requireNamespace("data.table", quietly = TRUE), "DATA_TABLE_REQUIRED")
    v5_assert(is.list(spec) && all(c("root", "out", "core", "files") %in% names(spec)), "VERIFIED_VENDOR_FEATURE_SOURCE_SPEC_REQUIRED")
    pins <- spec$files
    for (p in pins) v5_check(v5_rebase(p$path), p$sha256)
    e <- new.env(parent = globalenv())
    e$root <- v5_rebase(spec$root)
    e$out <- v5_rebase(spec$out)
    e$feature_cache <- new.env(parent = emptyenv())
    e$num <- function(x) suppressWarnings(as.numeric(as.character(x)))
    e$assert <- v5_assert
    e$fread <- data.table::fread
    e$read_json <- jsonlite::read_json
    e$cp_load <- function(path, strict = TRUE) {
        pp <- Filter(function(p) identical(normalizePath(v5_rebase(p$path), mustWork = TRUE), normalizePath(path, 
            mustWork = TRUE)), pins)
        v5_assert(length(pp) == 1L, paste("UNPINNED_VENDOR_CACHE", path))
        v5_read_source(c(pp[[1]], list(unwrap = TRUE)))
    }
    e$core <- v5_read_source(spec$core)
    dictfile <- file.path(e$root, "vendor/reused/8.clinical_omics.R")
    allowed <- "v2_feature_dictionary"
    for (expr in parse(file = dictfile)) if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) && 
        as.character(expr[[2]]) %in% allowed) 
        eval(expr, envir = e)
    v5_assert(is.function(e$v2_feature_dictionary), "FROZEN_FEATURE_DICTIONARY_FUNCTION_REQUIRED")
    e$dict <- function() e$v2_feature_dictionary(e$root)
    sys.source(file.path(e$root, "vendor/features.R"), envir = e)
    v5_assert(is.function(e$load_feature), "FROZEN_FEATURE_LOADER_REQUIRED")
    e$dictionary <- e$dict()
    e$verify <- function() {
        for (p in pins) v5_check(v5_rebase(p$path), p$sha256)
        invisible(TRUE)
    }
    e
}

v5_pure_risk <- function(v9) {
    rp <- file.path(dirname(v9), "border_population_divergence_explanation_v1/f/1.prepare.R")
    v5_check(rp, birth_private_setting("PRIVATE_HASH_6"))
    pure <- new.env(parent = globalenv())
    pure$pd_assert <- v5_assert
    pure$pd_hash <- function(x) digest::digest(x, algo = "sha256")
    sys.source(rp, envir = pure)
    pure
}

v5_dated_bridge <- function(parent, legacy, b18, pure) {
    ids <- as.character(b18$d_full$eid)
    j <- v5_join(parent$data, ids)
    p <- parent$data[j, , drop = FALSE]
    lj <- v5_join(legacy$data, ids)
    v5_assert(all(c("date_attend", "date_death_uncapped", "date_lost_registry", "record_observation_eligible", 
        "P_common") %in% names(p)) && all(!is.na(p$P_common) & p$P_common) && all(!is.na(p$record_observation_eligible) & 
        p$record_observation_eligible), "FROZEN_COMMON_TARGET_OBSERVATION_PROOF_REQUIRED")
    v5_assert(is.character(parent$baseline_note) && grepl("baseline union|BASELINE.*UNION", parent$baseline_note, 
        ignore.case = TRUE), "FROZEN_ALL_E1_BASELINE_FREE_TARGET_EVIDENCE_REQUIRED")
    ans <- data.frame(eid = ids, record_observation_eligible = p$record_observation_eligible)
    audit <- list()
    for (e in c("250", "250.2", "411", "411.3", "411.4")) {
        if (e %in% c("250", "411.3")) {
            v5_assert(e %in% names(legacy$events) && e %in% names(legacy$baseline) && inherits(legacy$events[[e]], 
                "Date"), "LEGACY_DM_AG_ACTUAL_DATES_REQUIRED")
            event <- legacy$events[[e]][lj]
            st <- as.integer(legacy$baseline[[e]][lj])
            source <- "EXACT_LEGACY_E0_EVENT_DATE_VECTOR"
        }
        else {
            r <- parent$risk$E1[[e]]
            v5_assert(is.data.frame(r) && identical(as.character(r$eid), as.character(parent$data$eid)) && all(b18$endpoints[[e]]$eligibility), 
                "FROZEN_E1_FULL_ELIGIBLE_PARENT_RISK_REQUIRED")
            r <- r[j, , drop = FALSE]
            days <- r$exit * 365.25
            v5_assert(all(is.finite(days)) && max(abs(days - round(days))) < 1e-07, "FROZEN_E1_EXIT_NOT_EXACT_INTEGER_DAY_NO_DATE_INFERENCE")
            event <- as.Date(rep(NA_real_, length(ids)), origin = "1970-01-01")
            ev <- r$event == 1L
            event[ev] <- as.Date(p$date_attend[ev]) + round(days[ev])
            st <- rep(0L, length(ids))
            source <- "EXACT_INTEGER_DAY_RECONSTRUCTION_OF_FROZEN_E1_OBSERVED_EVENT_DATE"
        }
        rr <- pure$pd_data_risk(p, event, st, p$record_observation_eligible, b18$cfg$admin)
        ep <- b18$endpoints[[e]]
        ok <- identical(rr$eligible, ep$eligibility) && isTRUE(all.equal(rr$data, ep$risk, tolerance = 1e-12, check.attributes = FALSE))
        v5_assert(ok, paste("DATED_BRIDGE_DOES_NOT_RECONSTRUCT_FROZEN_BASELINE_ENDPOINT", e))
        key <- gsub("[.]", "_", e)
        ans[[paste0("event__", key)]] <- event
        ans[[paste0("status__", key)]] <- st
        audit[[e]] <- data.frame(endpoint_id = e, source = source, baseline_risk_reconstructed_identically = ok, 
            N_parent = length(ids), N_dated_events = sum(!is.na(event)), endpoint_version = ep$endpoint_version, 
            admin = b18$cfg$admin, stringsAsFactors = FALSE)
    }
    list(data = ans, audit = v5_bind(audit))
}

v5_recover_social <- function(spec, loaded, ids, old, baseline_date) {
    v5_assert(all(c("raw", "dictionary", "encoder") %in% names(spec)), "EXACT_SOCIAL_RECOVERY_PINS_REQUIRED")
    raw <- v5_read_source(spec$raw)
    defs <- v5_read_source(spec$dictionary)
    v5_check(v5_rebase(spec$encoder$path), spec$encoder$sha256)
    encoder <- loaded$discrete$c40p_encoder(v5_rebase(spec$encoder$path))
    j <- v5_join(raw, ids)
    oj <- v5_join(old, ids, FALSE)
    v5_assert("p53_i0" %in% names(raw) && all(!is.na(as.Date(raw$p53_i0[j])) & as.Date(raw$p53_i0[j]) == as.Date(baseline_date)), 
        "SOCIAL_RAW_ACTUAL_BASELINE_DATES_CONFLICT")
    out <- data.frame(eid = ids)
    audit <- list()
    targets <- loaded$discrete$c40p_targets(raw[j, , drop = FALSE], encoder)
    for (id in unlist(spec$indicators, use.names = FALSE)) {
        def <- Filter(function(x) identical(x$id, id), defs)
        v5_assert(length(def) == 1L && identical(def[[1]]$type, "continuous"), "EXACT_FROZEN_SOCIAL_CONTINUOUS_DEFINITION_REQUIRED")
        def <- def[[1]]
        v5_assert(def$field %in% names(raw), "SOCIAL_RAW_FIELD_MISSING")
        x <- encoder$encode_marker(raw[j, , drop = FALSE], def)$value
        v5_assert(is.numeric(x) && length(x) == length(ids), "FROZEN_SOCIAL_ENCODER_NOT_NUMERIC")
        overlap <- !is.na(oj)
        v5_assert(id %in% names(old) && isTRUE(all.equal(x[overlap], old[[id]][oj[overlap]], tolerance = 0, check.attributes = FALSE)), 
            paste("FROZEN_SOCIAL_RECOVERY_OVERLAP_CHANGED", id))
        out[[id]] <- x
        tg <- targets[[def$eligible]]
        v5_assert(is.logical(tg) && !anyNA(tg), "FROZEN_SOCIAL_APPLICABILITY_TARGET_REQUIRED")
        out[[paste0("eligible__", id)]] <- tg
        audit[[id]] <- data.frame(indicator_id = id, N_parent = length(ids), N_target = sum(tg), target_rule = def$eligible, 
            N_old_cache_overlap = sum(overlap), N_new_parent_rows_recovered = sum(!overlap & !is.na(x)), N_valid = sum(!is.na(x)), 
            overlap_values_exact = TRUE, source = "EXISTING_FROZEN_CONTINUOUS_ENCODER_FULL_PARENT_RAW", stringsAsFactors = FALSE)
    }
    list(data = out, audit = v5_bind(audit))
}

v5_stage <- function(ctx, stage, endpoint = NULL, detail = NULL) {
    m <- ctx$monitor
    if (!is.environment(m)) 
        return(invisible(NULL))
    now <- as.numeric(Sys.time())
    mono <- proc.time()[["elapsed"]]
    if (is.null(endpoint)) 
        endpoint <- if (is.null(m$endpoint_id)) 
            "SHARED"
        else m$endpoint_id
    m$stage <- stage
    m$endpoint_id <- endpoint
    m$stage_started_epoch <- now
    m$stage_started_monotonic <- mono
    m$last_progress_epoch <- now
    m$sequence <- if (is.null(m$sequence)) 
        1L
    else m$sequence + 1L
    event <- list(schema = "V5_WORKER_STAGE_V1", run_id = ctx$run_id, worker_id = v5_env("V5_WORKER_SLOT"), host = v5_env("V5_HOST"), 
        pid = Sys.getpid(), process_start = m$process_start, indicator_id = m$indicator_id, bootstrap_id = m$bootstrap_id, 
        attempt_id = m$attempt_id, endpoint_id = endpoint, stage = stage, stage_detail = detail, stage_started_utc = v5_utc(), 
        stage_started_epoch = now, stage_started_monotonic = mono, stage_elapsed_s = 0, package_started_epoch = m$package_started_epoch, 
        package_elapsed_s = if (is.null(m$package_started_monotonic)) NA_real_ else mono - m$package_started_monotonic, 
        last_progress_epoch = now, sequence = m$sequence, endpoint_status = m$endpoint_status)
    ev <- v5_env("V5_WORKER_EVENTS")
    if (nzchar(ev)) {
        dir.create(dirname(ev), recursive = TRUE, showWarnings = FALSE)
        cat(jsonlite::toJSON(event, auto_unbox = TRUE, na = "null", null = "null", digits = 17), "\n", file = ev, 
            append = TRUE, sep = "")
    }
    sp <- v5_env("V5_WORKER_STATE")
    if (nzchar(sp)) {
        dir.create(dirname(sp), recursive = TRUE, showWarnings = FALSE)
        tmp <- paste0(sp, ".", Sys.getpid(), ".pending")
        jsonlite::write_json(event, tmp, auto_unbox = TRUE, na = "null", null = "null", digits = 17)
        v5_assert(file.rename(tmp, sp), "WORKER_STATE_ATOMIC_PUBLICATION_FAILED")
    }
    invisible(event)
}

v5_frame_endpoint <- function(ctx) {
    fs <- sys.frames()
    for (i in rev(seq_along(fs))) if (exists("e", envir = fs[[i]], inherits = FALSE)) {
        e <- get("e", envir = fs[[i]], inherits = FALSE)
        if (is.character(e) && length(e) == 1L && e %in% ctx$cfg$outcomes) 
            return(e)
    }
    if (is.environment(ctx$monitor) && !is.null(ctx$monitor$endpoint_id)) 
        ctx$monitor$endpoint_id
    else "SHARED"
}

v5_install_hooks <- function(ctx) {
    wrap <- function(env, name, stage) {
        original <- env[[name]]
        v5_assert(is.function(original), paste("FROZEN_FUNCTION_REQUIRED", name))
        env[[name]] <- local({
            fn <- original
            st <- stage
            nm <- name
            function(...) {
                ep <- v5_frame_endpoint(ctx)
                v5_stage(ctx, st, ep)
                value <- fn(...)
                if (nm == "b3_missing" && is.environment(ctx$monitor)) {
                  aa <- list(...)
                  d <- aa[[1]]
                  pr <- match(as.character(d$eid), as.character(ctx$d_full$eid))
                  if (!anyNA(pr)) {
                    key <- v5_hash(pr)
                    ctx$monitor$selection_weight_hashes[[key]] <- list(selection_weight_hash = v5_hash(format(value$w, 
                      digits = 17, scientific = TRUE)), selection_sample_hash = v5_hash(pr[value$complete]), selection_target_hash = key)
                  }
                }
                if (nm == "c40_exact_integral" && ep %in% ctx$cfg$outcomes) 
                  ctx$monitor$endpoint_status[[ep]] <- "COMPUTED_NOT_SAVED"
                value
            }
        })
    }
    hooks <- c(b3_missing = "SELECTION_WEIGHT", b3_censor = "CENSOR_MODEL", b3_logit = "OUTCOME_FIT", b3_predict_scenarios = "STANDARDIZE")
    for (k in names(hooks)) wrap(ctx$engine, k, hooks[[k]])
    wrap(ctx$discrete, "c40_discrete_fit", "MEDIATOR_FIT")
    wrap(ctx$discrete, "c40_exact_integral", "STANDARDIZE")
    ctx
}

v5_landmark <- function(ids, base, visit, event_date, status, observed, admin, pure, endpoint_version) {
    required <- c("date", "age", "center", "date_death_uncapped", "date_lost_registry")
    v5_assert(all(required %in% names(visit)) && nrow(visit) == length(ids), "ACTUAL_IMAGING_VISIT_DATA_REQUIRED")
    d <- base
    origin <- as.Date(visit$date)
    v5_assert(inherits(event_date, "Date"), "ACTUAL_FROZEN_EVENT_DATES_REQUIRED")
    d$date_attend <- origin
    d$date_death_uncapped <- as.Date(visit$date_death_uncapped)
    d$date_lost_registry <- as.Date(visit$date_lost_registry)
    d$age <- as.numeric(visit$age)
    d$center <- factor(visit$center)
    d$assessment_year <- factor(format(origin, "%Y"))
    r <- pure$pd_data_risk(d, event_date, status, observed, admin)
    eligible <- r$eligible & is.finite(d$age) & !is.na(d$center) & !is.na(d$assessment_year)
    known <- r$data$exit >= 5 | (r$data$exit <= 5 & (r$data$event == 1 | r$data$death_event == 1))
    support <- any(eligible & is.finite(r$data$exit) & r$data$exit >= 5)
    list(data = d, endpoint = list(status = if (support) "AVAILABLE_FROZEN" else "HOLD", reason = if (support) "" else "NO_FIVE_YEAR_LANDMARK_FOLLOWUP_SUPPORT_ORIGINAL_ADMIN_UNCHANGED", 
        risk = r$data, eligibility = eligible, endpoint_version = paste0(endpoint_version, "_INSTANCE2_LANDMARK"), 
        risk_hash = v5_hash(c(endpoint_version, as.character(origin), as.character(event_date), as.character(admin)))), 
        audit = data.frame(endpoint_version = endpoint_version, time_origin = "ACTUAL_INSTANCE2_DATE", admin = admin, 
            N_parent = length(ids), N_eligible = sum(eligible), N_known_5y = sum(eligible & known, na.rm = TRUE), 
            N_followup_ge5 = sum(eligible & r$data$exit >= 5, na.rm = TRUE), N_prevalent_at_imaging = sum(!is.na(event_date) & 
                event_date <= origin, na.rm = TRUE), N_events5 = sum(eligible & r$data$event == 1 & r$data$exit <= 
                5, na.rm = TRUE), status = if (support) "AVAILABLE_FROZEN" else "HOLD"))
}

v5_prepare <- function(root, run = v5_env("V5_RUN_DIR")) {
    v5_assert(identical(v5_env("V5_REAL_EXECUTION_AUTHORIZED"), "YES_FOR_THIS_RUN"), "EXPLICIT_V5_EXECUTION_AUTHORIZATION_REQUIRED")
    v5_assert(Sys.info()[["sysname"]] != "Darwin" && v5_env("V5_HOST") == "winpc", "PARTICIPANT_PREPARATION_WINPC_ONLY")
    root <- normalizePath(root, mustWork = TRUE)
    run <- v5_run(run, TRUE)
    v5_deps()
    scope <- v5_scope(root)
    loaded <- v5_sources(root)
    sp <- file.path(root, "config/adapter_sources.json")
    v5_assert(file.exists(sp), "VERIFIED_ADAPTER_SOURCES_JSON_REQUIRED")
    src <- jsonlite::fromJSON(sp, simplifyVector = FALSE)
    final_private <- file.path(run, "private")
    if (!dir.exists(final_private)) 
        v5_assert(dir.create(final_private), "PRIVATE_ROOT_CREATE_FAILED")
    v5_assert(!file.exists(file.path(final_private, "V5_INPUT_LOCK.json")), "PRESERVE_EXISTING_COMPLETE_PREPARATION")
    private <- file.path(final_private, paste0("preparation_attempt_", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), 
        "_", Sys.getpid()))
    v5_assert(!dir.exists(private) && dir.create(private), "UNIQUE_PREPARATION_ATTEMPT_REQUIRED")
    social_path <- file.path(root, "config/social_recovery.json")
    social_sha <- NULL
    if (file.exists(social_path)) {
        social_sha <- v5_sha(social_path)
        src$social_recovery <- jsonlite::fromJSON(social_path, simplifyVector = FALSE)
        src$social_recovery$dictionary$format <- "json"
        src$social_recovery$indicators <- names(Filter(function(x) identical(x$source, "social"), src$features))
    }
    b18 <- v5_read_source(src$b18_cache)
    catctx <- v5_read_source(src$categorical_cache)
    ids <- as.character(b18$d_full$eid)
    v5_assert(length(ids) == 123350L && identical(ids, as.character(catctx$d_full$eid)), "EXACT_PARENT_AND_CATEGORICAL_BRIDGE_REQUIRED")
    for (k in c("parent_frame_hash", "cluster_order_hash", "parent_cluster_membership_hash")) v5_assert(identical(b18[[k]], 
        catctx[[k]]), paste("PARENT_IDENTITY_CONFLICT", k))
    cfg <- b18$cfg
    cfg$outcomes <- c("250", "250.2", "411", "411.3", "411.4")
    cfg$base_covariates <- as.character(unlist(cfg$base_covariates))
    v5_assert(identical(cfg$base_covariates, c("B", "age", "sex", "ethnicity", "center", "assessment_year", "S")) && 
        cfg$horizon == 5 && cfg$admin == "2022-05-31" && cfg$seed == 20261023, "FROZEN_STATISTICAL_CONTRACT_CHANGED")
    d <- b18$d_full
    features <- variables <- eligibility <- list()
    metadata <- list()
    sources <- list()
    if (!is.null(src$sources)) 
        for (k in names(src$sources)) sources[[k]] <- v5_read_source(src$sources[[k]])
    vendor <- if (!is.null(src$vendor_features)) 
        v5_vendor_features(src$vendor_features)
    else NULL
    if (!is.null(vendor)) {
        md <- as.data.frame(vendor$core$metadata)
        cd <- as.data.frame(vendor$core$data)
        mj <- v5_join(md, as.character(cd$eid))
        md <- md[mj, , drop = FALSE]
        for (k in intersect(c("date_attend", "date_death_uncapped", "date_lost_registry"), names(cd))) md[[k]] <- cd[[k]]
        sources$vendor_visit <- md
        sources$vendor_core_data <- cd
    }
    social_audit <- NULL
    if (!is.null(src$social_recovery)) {
        v5_assert(!is.null(sources$social) && !is.null(sources$vendor_core_data), "SOCIAL_CACHE_AND_PARENT_DATES_REQUIRED")
        cd <- sources$vendor_core_data
        bj <- v5_join(cd, ids)
        recovered <- v5_recover_social(src$social_recovery, loaded, ids, sources$social, cd$date_attend[bj])
        sources$social_recovered <- recovered$data
        social_audit <- recovered$audit
        for (id in unlist(src$social_recovery$indicators, use.names = FALSE)) {
            v5_assert(identical(src$features[[id]]$source, "social"), "SOCIAL_RECOVERY_ONLY_OLD_SOCIAL_FEATURES")
            src$features[[id]]$source <- "social_recovered"
            src$features[[id]]$target_source <- "social_recovered"
            src$features[[id]]$target_column <- paste0("eligible__", id)
        }
    }
    dated_audit <- NULL
    if (!is.null(src$dated_frozen_bridge)) {
        db <- src$dated_frozen_bridge
        dp <- v5_read_source(db$parent_cache)
        dl <- v5_read_source(db$legacy_core)
        bridge <- v5_dated_bridge(dp, dl, b18, v5_pure_risk(loaded$v9))
        sources$dated_bridge <- bridge$data
        dated_audit <- bridge$audit
        if (is.null(src$dated_endpoints)) {
            src$dated_endpoints <- setNames(lapply(cfg$outcomes, function(e) {
                key <- gsub("[.]", "_", e)
                list(source = "dated_bridge", event_column = paste0("event__", key), baseline_column = paste0("status__", 
                  key), observation_column = "record_observation_eligible")
            }), cfg$outcomes)
        }
    }
    profile <- if (!is.null(src$profile_cache)) 
        v5_read_source(src$profile_cache)
    else NULL
    for (i in seq_len(nrow(scope))) {
        row <- scope[i, , drop = FALSE]
        id <- row$indicator_id
        is_cmr <- row$category == "CMR"
        spec <- src$features[[id]]
        if (row$type != "continuous") {
            ca <- Filter(function(x) identical(x$indicator_id, id), catctx$cfg$indicators)
            v5_assert(length(ca) == 1L && identical(ca[[1]]$kind, row$type), "FROZEN_CATEGORICAL_ID_TYPE_MISMATCH")
            ca <- ca[[1]]
            d[[ca$column]] <- catctx$d_full[[ca$column]]
            variables[[id]] <- catctx$variables[[id]]
            eligibility[[id]] <- catctx$eligibility[[id]]
            features[[id]] <- ca
        }
        else {
            old <- Filter(function(x) identical(x$feature_id, id), b18$cfg$features)
            pm <- if (!is.null(profile)) 
                profile$feature_metadata[profile$feature_metadata$feature_id == id, , drop = FALSE]
            else data.frame()
            if (length(old) == 1L) {
                ca <- old[[1]]
                x <- d[[ca$column]]
                target <- b18$eligibility[[ca$target_key]]
                ca$inverse_transform <- if (ca$prepared_transform == "log") 
                  "exp"
                else "identity"
                status <- "AVAILABLE_FROZEN"
                reason <- ""
            }
            else if (!is.null(vendor) && id %in% vendor$dictionary$mediator_id) {
                vr <- as.list(vendor$dictionary[vendor$dictionary$mediator_id == id, ])
                v5_assert(length(vr$mediator_id) == 1L, "UNIQUE_EXACT_VENDOR_FEATURE_REQUIRED")
                vf <- vendor$load_feature(vendor$core, vr)
                j <- v5_join(vf$data, ids)
                x <- vf$data$marker[j]
                fitted_unit <- if (isTRUE(vr$indexed)) 
                  paste0(vr$unit, "/m2")
                else vr$unit
                if (identical(vr$transform, "log")) 
                  fitted_unit <- paste0("log(", fitted_unit, ")")
                ca <- list(feature_id = id, candidate_id = id, column = paste0("v5__", make.names(id)), platform = row$category, 
                  target_key = id, raw_unit = vr$unit, prepared_unit = fitted_unit, prepared_transform = vr$transform, 
                  preprocessing = if (isTRUE(vr$indexed)) "same_instance2_BSA_division" else "FROZEN_FEATURE_LOADER_QC", 
                  inverse_transform = if (vr$transform == "log") "exp" else "identity", original_extra_covariates = as.character(vf$extra), 
                  source_qc = vf$source_qc, source_column = vr$column)
                target <- if (startsWith(id, "nmr__")) 
                  b18$eligibility$nmr_measured_i0
                else rep(TRUE, length(ids))
                status <- "AVAILABLE_FROZEN"
                reason <- ""
                if (is.null(spec) || !isTRUE(spec$physical_support_verified)) {
                  status <- "HOLD"
                  reason <- "PREPARED_VALUES_AVAILABLE_PHYSICAL_SUPPORT_SPEC_NOT_YET_VERIFIED"
                }
                else for (k in c("physical_lower_bound", "physical_upper_bound")) if (!is.null(spec[[k]])) 
                  ca[[k]] <- as.numeric(spec[[k]])
                ca$technical_covariates <- character()
                if (!is.null(spec$m3_technical_covariates)) 
                  for (k in unlist(spec$m3_technical_covariates)) {
                    v5_assert(k %in% vf$extra && k %in% names(vf$data), "TECHNICAL_COVARIATE_NOT_IN_FROZEN_LOADER")
                    nm <- paste0("tech__", make.names(id), "__", k)
                    d[[nm]] <- vf$data[[k]][j]
                    ca$technical_covariates <- c(ca$technical_covariates, nm)
                  }
                if (!is_cmr) {
                  v5_assert(all(is.na(vf$data$measurement_date[j]) | as.Date(vf$data$measurement_date[j]) == as.Date(sources$vendor_core_data$date_attend[match(ids, 
                    as.character(sources$vendor_core_data$eid))])), "BASELINE_MEASUREMENT_DATE_CONFLICT")
                }
            }
            else if (!is.null(spec)) {
                v5_assert(all(c("source", "column", "prepared_unit", "prepared_transform", "inverse_transform") %in% 
                  names(spec)), paste("EXPLICIT_SCALE_SPEC_REQUIRED", id))
                sd <- sources[[spec$source]]
                v5_assert(is.data.frame(sd) && spec$column %in% names(sd), paste("EXPLICIT_SOURCE_COLUMN_MISSING", 
                  id))
                j <- v5_join(sd, ids, FALSE)
                x <- sd[[spec$column]][j]
                v5_assert(is.numeric(x), "PREPARED_CONTINUOUS_NUMERIC_REQUIRED_NO_REENCODING")
                ca <- c(list(feature_id = id, candidate_id = id, column = paste0("v5__", make.names(id)), platform = row$category, 
                  target_key = id), spec)
                ca$column <- paste0("v5__", make.names(id))
                target <- rep(TRUE, length(ids))
                if (!is.null(spec$target_source)) {
                  td <- sources[[spec$target_source]]
                  tj <- v5_join(td, ids)
                  target <- td[[spec$target_column]][tj]
                  v5_assert(is.logical(target) && !anyNA(target), "EXPLICIT_TARGET_LOGICAL_REQUIRED")
                }
                status <- "AVAILABLE_FROZEN"
                reason <- ""
            }
            else if (nrow(pm) == 1L) {
                ca <- list(feature_id = id, candidate_id = id, column = paste0("v5__", make.names(id)), platform = row$category, 
                  target_key = id, prepared_unit = pm$unit, prepared_transform = pm$transform, inverse_transform = if (pm$transform == 
                    "log") "exp" else "identity")
                j <- v5_join(profile$data, ids)
                x <- profile$data[[pm$column]][j]
                target <- if (startsWith(id, "nmr__")) 
                  profile$eligibility$nmr_measured_i0[j]
                else rep(TRUE, length(ids))
                status <- "HOLD"
                reason <- "PREPARED_VALUES_AVAILABLE_PHYSICAL_SUPPORT_SPEC_NOT_YET_VERIFIED"
            }
            else {
                ca <- list(feature_id = id, candidate_id = id, column = paste0("v5__", make.names(id)), platform = row$category, 
                  target_key = id, prepared_unit = NA_character_, prepared_transform = NA_character_, inverse_transform = NA_character_)
                x <- rep(NA_real_, length(ids))
                target <- rep(TRUE, length(ids))
                status <- "HOLD"
                reason <- "EXACT_PREPARED_SOURCE_COLUMN_UNIT_TRANSFORM_NOT_YET_VERIFIED"
            }
            if (!is.null(spec) && length(old) == 1L) 
                for (k in c("physical_lower_bound", "physical_upper_bound")) if (!is.null(spec[[k]])) 
                  v5_assert(identical(as.numeric(spec[[k]]), as.numeric(ca[[k]])), paste("ORIGINAL_FROZEN_BOUND_CHANGED", 
                    id, k))
            d[[ca$column]] <- x
            features[[id]] <- ca
            eligibility[[id]] <- target
            variables[[id]] <- list(status = status, reason = reason, measurement_instance = row$measurement_instance, 
                source_label_retained = row$source_label_retained)
        }
        ca <- features[[id]]
        metadata[[id]] <- data.frame(indicator_id = id, feature_id = id, type = row$type, domain = row$category, 
            time_origin = if (is_cmr) 
                "imaging"
            else "baseline", column = ca$column, raw_unit = if (is.null(ca$raw_unit)) 
                if (is.null(ca$prepared_unit)) 
                  "category"
                else ca$prepared_unit
            else ca$raw_unit, prepared_unit = if (is.null(ca$prepared_unit)) 
                "category"
            else ca$prepared_unit, preprocessing = if (is.null(ca$preprocessing)) 
                "INHERITED_PREPARED_VALUES"
            else ca$preprocessing, prepared_transform = if (is.null(ca$prepared_transform)) 
                "FROZEN_FACTOR"
            else ca$prepared_transform, inverse_transform = if (is.null(ca$inverse_transform)) 
                "NOT_APPLICABLE"
            else ca$inverse_transform, status = variables[[id]]$status, reason = if (is.null(variables[[id]]$reason)) 
                ""
            else variables[[id]]$reason, N_target = sum(eligibility[[id]]), N_observed = sum(eligibility[[id]] & 
                !is.na(d[[ca$column]])), stringsAsFactors = FALSE)
    }
    imaging <- list(status = "HOLD", reason = "ACTUAL_INSTANCE2_VISIT_AND_FROZEN_DATED_ENDPOINTS_NOT_YET_VERIFIED")
    imaging_audit <- list()
    if (!is.null(src$imaging)) {
        im <- src$imaging
        vd <- sources[[im$source]]
        j <- v5_join(vd, ids, FALSE)
        cols <- c(date = im$date_column, age = im$age_column, center = im$center_column, date_death_uncapped = im$death_column, 
            date_lost_registry = im$loss_column)
        v5_assert(all(cols %in% names(vd)), "EXPLICIT_IMAGING_VISIT_COLUMNS_MISSING")
        visit <- vd[j, cols, drop = FALSE]
        names(visit) <- names(cols)
        pure <- v5_pure_risk(loaded$v9)
        ieps <- list()
        imd <- NULL
        for (e in cfg$outcomes) {
            es <- src$dated_endpoints[[e]]
            v5_assert(!is.null(es), "FIVE_FROZEN_DATED_ENDPOINT_SPECS_REQUIRED")
            ed <- sources[[es$source]]
            ej <- v5_join(ed, ids)
            z <- v5_landmark(ids, d, visit, as.Date(ed[[es$event_column]][ej]), as.integer(ed[[es$baseline_column]][ej]), 
                as.logical(ed[[es$observation_column]][ej]), cfg$admin, pure, b18$endpoints[[e]]$endpoint_version)
            ieps[[e]] <- z$endpoint
            imaging_audit[[e]] <- z$audit
            imd <- z$data
        }
        for (e in c("250.2", "411")) {
            r <- ieps[[e]]$risk
            k <- gsub("[.]", "_", e)
            imd[[paste0("aux_t_", k)]] <- pmin(r$exit, 5)
            imd[[paste0("aux_event_", k)]] <- as.integer(r$event == 1 & r$exit <= 5)
            imd[[paste0("aux_death_", k)]] <- as.integer(r$death_event == 1 & r$exit <= 5)
        }
        imaging <- list(status = "AVAILABLE_FROZEN", data = imd, endpoints = ieps, reason = "", covariate_time = "ACTUAL_INSTANCE2_AGE_YEAR_CENTER", 
            admin = cfg$admin)
    }
    for (id in scope$indicator_id[scope$category == "CMR"]) if (imaging$status != "AVAILABLE_FROZEN") {
        variables[[id]]$status <- "HOLD"
        variables[[id]]$reason <- imaging$reason
        metadata[[id]]$status <- "HOLD"
        metadata[[id]]$reason <- imaging$reason
    }
    if (!is.null(vendor)) 
        vendor$verify()
    cfg$features <- unname(features[scope$indicator_id[scope$type == "continuous"]])
    cfg$candidates <- cfg$features
    cfg$indicators <- unname(features[scope$indicator_id[scope$type != "continuous"]])
    cfg$outcomes <- c("250", "250.2", "411", "411.3", "411.4")
    ctx <- list(cfg = cfg, d_full = d, features = features, variables = variables, eligibility = eligibility, endpoints = b18$endpoints, 
        imaging = imaging, scope = scope, run_id = jsonlite::fromJSON(file.path(root, "config/runtime.json"))$run_id, 
        parent_frame_hash = b18$parent_frame_hash, cluster_order_hash = b18$cluster_order_hash, parent_cluster_membership_hash = b18$parent_cluster_membership_hash, 
        scope_hash = v5_sha(file.path(root, "input/APPROVED_165_ORIGINAL.csv")), method_hash = v5_hash(loaded$source_sha256), 
        source_hash = v5_hash(vapply(sort(list.files(file.path(root, "f"), full.names = TRUE, pattern = "[.]R$")), 
            v5_sha, "")))
    cache <- file.path(private, "V5_INPUTS.private.rds")
    saveRDS(ctx, cache, compress = FALSE)
    v5_csv(v5_bind(metadata), file.path(private, "FEATURE_INPUT_PREFLIGHT.csv"))
    if (length(imaging_audit)) 
        v5_csv(v5_bind(imaging_audit), file.path(private, "CMR_LANDMARK_PREFLIGHT.csv"))
    if (!is.null(dated_audit)) 
        v5_csv(dated_audit, file.path(private, "DATED_ENDPOINT_RECONSTRUCTION_PREFLIGHT.csv"))
    if (!is.null(social_audit)) 
        v5_csv(social_audit, file.path(private, "SOCIAL_FULL_PARENT_RECOVERY_PREFLIGHT.csv"))
    lock <- list(complete = TRUE, models_run = 0L, run_id = ctx$run_id, prepared_inputs_sha256 = v5_sha(cache), 
        approved_csv_sha256 = ctx$scope_hash, adapter_sources_sha256 = v5_sha(sp), social_recovery_sha256 = social_sha, 
        parent_source_sha256 = as.list(loaded$source_sha256), adapter_source_hash = ctx$source_hash, method_hash = ctx$method_hash, 
        parent_frame_hash = ctx$parent_frame_hash, cluster_order_hash = ctx$cluster_order_hash, parent_cluster_membership_hash = ctx$parent_cluster_membership_hash, 
        cache_sha256 = b18$cfg$cache_sha256, cohort_sha256 = b18$cfg$cohort_sha256, N_parent = length(ids), N_indicators = 165L, 
        N_continuous = 137L, N_discrete = 28L, participant_data_location = "G_ONLY_NO_MAC_BULK_COPY", prepared_bytes = as.numeric(file.info(cache)$size), 
        old_outputs_modified = FALSE)
    v5_json(lock, file.path(private, "V5_INPUT_LOCK.json"))
    files <- setdiff(list.files(private, full.names = TRUE), file.path(private, "V5_INPUT_LOCK.json"))
    for (p in files) v5_publish(p, file.path(final_private, basename(p)))
    v5_publish(file.path(private, "V5_INPUT_LOCK.json"), file.path(final_private, "V5_INPUT_LOCK.json"))
    lock
}

v5_load <- function(root, run = v5_env("V5_RUN_DIR")) {
    v5_assert(identical(v5_env("V5_REAL_EXECUTION_AUTHORIZED"), "YES_FOR_THIS_RUN"), "EXPLICIT_V5_EXECUTION_AUTHORIZATION_REQUIRED")
    v5_deps(TRUE)
    root <- normalizePath(root, mustWork = TRUE)
    run <- v5_run(run)
    scope <- v5_scope(root)
    loaded <- v5_sources(root)
    lp <- file.path(run, "private/V5_INPUT_LOCK.json")
    lock <- jsonlite::fromJSON(lp, simplifyVector = FALSE)
    v5_assert(isTRUE(lock$complete) && lock$models_run == 0L, "PREPARED_V5_INPUT_LOCK_REQUIRED")
    cache <- file.path(run, "private/V5_INPUTS.private.rds")
    v5_check(cache, lock$prepared_inputs_sha256)
    v5_assert(identical(v5_hash(vapply(sort(list.files(file.path(root, "f"), full.names = TRUE, pattern = "[.]R$")), 
        v5_sha, "")), lock$adapter_source_hash), "ADAPTER_SOURCE_CHANGED_AFTER_PREPARE")
    v5_assert(identical(as.character(unlist(lock$parent_source_sha256)), as.character(loaded$source_sha256)), "FROZEN_ENGINE_CHANGED_AFTER_PREPARE")
    v5_check(file.path(root, "config/adapter_sources.json"), lock$adapter_sources_sha256)
    if (!is.null(lock$social_recovery_sha256)) 
        v5_check(file.path(root, "config/social_recovery.json"), lock$social_recovery_sha256)
    ctx <- readRDS(cache)
    v5_assert(identical(ctx$scope_hash, v5_sha(file.path(root, "input/APPROVED_165_ORIGINAL.csv"))), "PREPARED_APPROVED_SCOPE_CHANGED")
    ctx$engine <- loaded$engine
    ctx$discrete <- loaded$discrete
    ctx$run <- run
    ctx$root <- root
    ctx$input_lock_sha256 <- v5_sha(lp)
    ctx$task_cache <- new.env(parent = emptyenv())
    ctx$matrix_cache <- new.env(parent = emptyenv())
    ctx$monitor <- new.env(parent = emptyenv())
    ctx$monitor$process_start <- v5_utc()
    ctx$monitor$endpoint_status <- as.list(setNames(rep("NOT_ATTEMPTED", 5), ctx$cfg$outcomes))
    v5_install_hooks(ctx)
}

v5_task_dir <- function(ctx, p) {
    v5_assert(is.character(p) && length(p) == 1L && !grepl("(^|/)[.][.](/|$)", gsub(birth_private_setting("PRIVATE_INPUT"), 
        "/", p)), "UNSAFE_TASK_PATH")
    parent <- dirname(p)
    ancestor <- parent
    parts <- character()
    while (!dir.exists(ancestor) && !identical(dirname(ancestor), ancestor)) {
        parts <- c(basename(ancestor), parts)
        ancestor <- dirname(ancestor)
    }
    ancestor <- normalizePath(ancestor, mustWork = TRUE)
    run <- normalizePath(ctx$run, mustWork = TRUE)
    v5_assert(identical(ancestor, run) || startsWith(paste0(ancestor, "/"), paste0(run, "/")), "TASK_OUTSIDE_NEW_V5_RUN")
    if (length(parts)) 
        parent <- do.call(file.path, as.list(c(ancestor, parts)))
    else parent <- ancestor
    if (!dir.exists(parent)) 
        dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    rp <- normalizePath(parent, mustWork = TRUE)
    v5_assert(startsWith(paste0(rp, "/"), paste0(run, "/")), "TASK_OUTSIDE_NEW_V5_RUN")
    dest <- file.path(rp, basename(p))
    v5_assert(!file.exists(dest) && !dir.exists(dest) && dir.create(dest), "PRESERVE_EXISTING_ATTEMPT_NO_OVERWRITE")
    dest
}

v5_failed_status <- function(reason) {
    if (grepl("cannot allocate|cannot open|permission|no space|write error|unused argument|could not find|object .* not found|subscript|transport|socket|FILE_MISSING|SHA256|PATH|SOURCE_CHANGED", 
        reason, ignore.case = TRUE)) 
        "FAILED_TECH"
    else "FAILED_STAT"
}

v5_base_row <- function(ctx, id, bid, e) {
    s <- ctx$scope[ctx$scope$indicator_id == id, , drop = FALSE]
    q <- data.frame(indicator_id = id, feature_id = id, label_zh = s$label_zh, type = s$type, domain = s$category, 
        endpoint_id = e, endpoint = e, endpoint_label = switch(e, `250` = "DM", `250.2` = "T2D", `411` = "IHD", 
            `411.3` = "AG", `411.4` = "CAS"), bootstrap_id = bid, time_origin = if (s$category == "CMR") 
            "imaging"
        else "baseline", horizon_years = 5L, status = "HOLD", reason = "", model_hash = ctx$method_hash, method_hash = ctx$method_hash, 
        input_hash = ctx$input_lock_sha256, scope_hash = ctx$scope_hash, draw_hash = NA_character_, target_hash = NA_character_, 
        sample_hash = NA_character_, stringsAsFactors = FALSE)
    for (k in c("R00", "R10", "R11", "TE", "DE", "IE", "PM", "PM_signed", "PM_signed_percent", "N_target", "N_complete", 
        "N_known_5y", "events_5y", "death_5y", "censor_5y")) q[[k]] <- NA_real_
    q
}

v5_compute <- function(ctx, indicator_id, bootstrap_id, output_dir, attempt_id) {
    bid <- v5_id(bootstrap_id)
    v5_assert(indicator_id %in% ctx$scope$indicator_id, "INDICATOR_OUTSIDE_APPROVED_165")
    m <- ctx$monitor
    m$indicator_id <- indicator_id
    m$bootstrap_id <- bid
    m$attempt_id <- attempt_id
    m$selection_weight_hashes <- list()
    m$package_started_epoch <- as.numeric(Sys.time())
    m$package_started_monotonic <- proc.time()[["elapsed"]]
    m$endpoint_status <- as.list(setNames(rep("NOT_ATTEMPTED", 5), ctx$cfg$outcomes))
    v5_stage(ctx, "TASK_CLAIMED", "SHARED")
    v5_stage(ctx, "LOAD_CACHE", "SHARED")
    scope <- ctx$scope[ctx$scope$indicator_id == indicator_id, , drop = FALSE]
    ca <- ctx$features[[indicator_id]]
    var <- ctx$variables[[indicator_id]]
    rows <- parameters <- diagnostics <- weights <- list()
    receipt <- NULL
    if (scope$type != "continuous" && identical(var$status, "AVAILABLE_FROZEN")) {
        v5_stage(ctx, "DRAW_PREPARE", "SHARED")
        ans <- ctx$discrete$c40_compute(ctx, indicator_id, bid, output_dir, backend = "reference")
        raw <- read.csv(file.path(output_dir, "RESULTS.csv"), stringsAsFactors = FALSE)
        receipt <- jsonlite::fromJSON(file.path(output_dir, "DRAW.json"), simplifyVector = FALSE)
        for (e in ctx$cfg$outcomes) {
            q <- v5_base_row(ctx, indicator_id, bid, e)
            r <- raw[as.character(raw$endpoint_id) == e, , drop = FALSE]
            for (k in names(r)) q[[k]] <- r[[k]]
            q$status <- if (r$status == "VALID") 
                "VALID"
            else if (grepl("BLOCKED|HOLD", r$status)) 
                "HOLD"
            else v5_failed_status(r$reason)
            wh <- m$selection_weight_hashes[[q$target_hash]]
            if (!is.null(wh)) 
                for (k in names(wh)) q[[k]] <- wh[[k]]
            rows[[e]] <- q
        }
    }
    else {
        td <- v5_task_dir(ctx, output_dir)
        v5_stage(ctx, "DRAW_PREPARE", "SHARED")
        actual <- ctx$discrete$c40_draw(ctx, bid)
        receipt <- actual$receipt
        receipt$indicator_id <- indicator_id
        receipt$method_hash <- ctx$method_hash
        receipt$input_hash <- ctx$input_lock_sha256
        for (k in names(ctx$discrete$c40_persist_draw(ctx, actual))) receipt[[k]] <- ctx$discrete$c40_persist_draw(ctx, 
            actual)[[k]]
        ix <- actual$ix
        cache <- new.env(parent = emptyenv())
        is_cmr <- scope$category == "CMR"
        dd_full <- if (is_cmr && ctx$imaging$status == "AVAILABLE_FROZEN") 
            ctx$imaging$data
        else ctx$d_full
        eps <- if (is_cmr && ctx$imaging$status == "AVAILABLE_FROZEN") 
            ctx$imaging$endpoints
        else ctx$endpoints
        warnings <- list()
        for (e in ctx$cfg$outcomes) {
            m$endpoint_id <- e
            q <- v5_base_row(ctx, indicator_id, bid, e)
            q$draw_hash <- receipt$draw_hash
            ep <- eps[[e]]
            cellwarnings <- character()
            if (var$status != "AVAILABLE_FROZEN") {
                q$reason <- var$reason
            }
            else if (ep$status != "AVAILABLE_FROZEN") {
                q$reason <- if (is.null(ep$reason)) 
                  "FROZEN_ENDPOINT_UNAVAILABLE"
                else ep$reason
            }
            else {
                oldwarning <- getOption("b3.warning_sink")
                options(b3.warning_sink = function(message) cellwarnings <<- c(cellwarnings, message))
                got <- tryCatch(withCallingHandlers({
                  target <- ctx$eligibility[[indicator_id]] & ep$eligibility
                  pr <- ix[target[ix]]
                  v5_assert(length(pr) > 0L, "NO_ELIGIBLE_TARGET")
                  d <- dd_full[pr, , drop = FALSE]
                  r <- ep$risk[pr, , drop = FALSE]
                  base <- unique(c(ctx$cfg$base_covariates, ca$technical_covariates))
                  q$target_hash <- v5_hash(pr)
                  q$N_target <- length(pr)
                  key <- v5_hash(c(indicator_id, bid, q$target_hash, ctx$method_hash, ctx$input_lock_sha256, ca$column, 
                    base))
                  if (exists(key, envir = cache, inherits = FALSE)) {
                    nu <- get(key, envir = cache)
                  }
                  else {
                    v5_stage(ctx, "SELECTION_WEIGHT", e)
                    mi <- ctx$engine$b3_missing(d, ca$column, base)
                    v5_stage(ctx, "MEDIATOR_FIT", "SHARED", list(target_hash = q$target_hash))
                    a <- ctx$engine$fp_a_path(d, ca$column, base, mi, ctx, paste0("single_", ca$candidate_id), 
                      bid)
                    physical <- ctx$engine$fp_physical(d, ca$column, mi, a$alpha, ctx)
                    v5_assert(all(physical$status == "VALID"), "SCENARIO_OUTSIDE_PHYSICAL_DOMAIN_NO_CLIPPING")
                    nu <- list(mi = mi, a = a, physical = physical)
                    assign(key, nu, envir = cache)
                  }
                  cc <- which(nu$mi$complete)
                  q$sample_hash <- v5_hash(pr[cc])
                  q$standardization_sample_hash <- q$sample_hash
                  q$N_complete <- length(cc)
                  q$prepared_unit <- ca$prepared_unit
                  q$prepared_transform <- ca$prepared_transform
                  q$inverse_transform <- ca$inverse_transform
                  q$endpoint_version <- ep$endpoint_version
                  q$risk_hash <- ep$risk_hash
                  q$target_id <- if (is_cmr) 
                    "INSTANCE2_SURVIVOR_ENDPOINT_DISEASE_FREE"
                  else ca$target_key
                  q$shared_nuisance_key <- key
                  q$selection_weight_hash <- v5_hash(format(nu$mi$w, digits = 17, scientific = TRUE))
                  q$selection_sample_hash <- q$sample_hash
                  q$selection_target_hash <- q$target_hash
                  q$formula_hash <- v5_hash(c(base, ca$column, "FROZEN_SELECTION_T2D_IHD_AUX", "FROZEN_RIDGE_CENSOR_THETA1", 
                    "M3_AS_MEAN_WEIGHTED_INTERCEPT_ONLY", "SPLINE_5_35_65_95"))
                  fit <- ctx$engine$fp_m3_effect(d, r, ca$column, base, nu$mi, diagnostics = bid == 0L, reuse_alpha = nu$a$alpha)
                  res <- fit$result
                  v5_assert(nrow(res) == 1L && all(is.finite(unlist(res[1, c("R00", "R10", "R11", "TE", "DE", "IE")]))), 
                    "NONFINITE_M3_EFFECT")
                  res$numerical_status <- res$status
                  res$status <- NULL
                  for (k in names(res)) q[[k]] <- res[[k]]
                  q$status <- "VALID"
                  q$reason <- ""
                  parameters[[e]] <- fit$parameters
                  diagnostics[[e]] <- v5_bind(list(nu$physical, nu$a$rows, fit$support, fit$residual_diagnostics, 
                    fit$censor_PH))
                  weights[[e]] <- ctx$engine$fp_weight_summary(nu$mi$w, "SELECTION_IPW", e)
                  TRUE
                }, warning = function(w) {
                  cellwarnings <<- c(cellwarnings, conditionMessage(w))
                  invokeRestart("muffleWarning")
                }), error = function(err) {
                  q$status <<- v5_failed_status(conditionMessage(err))
                  q$reason <<- conditionMessage(err)
                  FALSE
                })
                options(b3.warning_sink = oldwarning)
            }
            ww <- unique(c(cellwarnings, if ("warnings" %in% names(q)) q$warnings else character()))
            q$warning <- paste(ww, collapse = " | ")
            warnings[[e]] <- data.frame(endpoint_id = rep(e, length(ww)), message = ww, stringsAsFactors = FALSE)
            rows[[e]] <- q
            m$endpoint_status[[e]] <- q$status
            v5_stage(ctx, "WRITE_CELL", e)
            v5_csv(q, file.path(td, paste0("CELL_", gsub("[.]", "_", e), ".csv")))
        }
        v5_stage(ctx, "WRITE_PACKAGE", "SHARED")
        v5_json(receipt, file.path(td, "DRAW.json"))
        saveRDS(parameters, file.path(td, "PARAMETERS.rds"), compress = FALSE)
        v5_csv(v5_bind(diagnostics), file.path(td, "DIAGNOSTICS.csv"))
        v5_csv(v5_bind(weights), file.path(td, "WEIGHTS.csv"))
        v5_csv(v5_bind(warnings), file.path(td, "WARNINGS.csv"))
    }
    for (e in names(rows)) {
        q <- rows[[e]]
        q$PM_signed <- q$PM
        q$PM_signed_percent <- 100 * q$PM
        q$censor_5y <- if ("early_censor" %in% names(q)) 
            q$early_censor
        else NA_real_
        if (!"target_id" %in% names(q)) 
            q$target_id <- if ("target_key" %in% names(q)) 
                q$target_key
            else NA_character_
        for (k in c("TE", "DE", "IE")) q[[paste0(k, "_per_1000")]] <- 1000 * q[[k]]
        for (k in c("cluster_order_hash", "parent_frame_hash", "parent_cluster_membership_hash", "multiplicity_hash")) q[[k]] <- receipt[[k]]
        q$source_hash <- ctx$source_hash
        q$attempt_id <- attempt_id
        q$host <- v5_env("V5_HOST")
        q$source_label_retained <- scope$source_label_retained
        rows[[e]] <- q
        m$endpoint_status[[e]] <- q$status
        v5_stage(ctx, "WRITE_CELL", e)
        v5_csv(q, file.path(output_dir, paste0("V5_CELL_", gsub("[.]", "_", e), ".csv")))
    }
    v5_stage(ctx, "WRITE_PACKAGE", "SHARED")
    res <- v5_bind(rows)
    v5_assert(nrow(res) == 5L && setequal(as.character(res$endpoint_id), ctx$cfg$outcomes), "FIVE_INDEPENDENT_ENDPOINT_CELLS_REQUIRED")
    out <- normalizePath(output_dir, mustWork = TRUE)
    v5_csv(res, file.path(out, "V5_RESULTS.csv"))
    v5_csv(data.frame(total_seconds = proc.time()[["elapsed"]] - m$package_started_monotonic, pid = Sys.getpid(), 
        host = v5_env("V5_HOST")), file.path(out, "V5_TIMING.csv"))
    files <- list.files(out, full.names = TRUE)
    files <- files[!file.info(files)$isdir]
    seal <- list(complete = TRUE, indicator_id = indicator_id, bootstrap_id = bid, attempt_id = attempt_id, method_hash = ctx$method_hash, 
        input_hash = ctx$input_lock_sha256, draw_hash = receipt$draw_hash, artifacts_sha256 = as.list(setNames(vapply(files, 
            v5_sha, ""), basename(files))), cell_status = as.list(setNames(res$status, res$endpoint_id)), total_seconds = proc.time()[["elapsed"]] - 
            m$package_started_monotonic)
    v5_json(seal, file.path(out, "V5_COMPLETE.json"))
    v5_stage(ctx, "TASK_FINISHED", "SHARED")
    list(ok = TRUE, indicator_id = indicator_id, bootstrap_id = bid, output_dir = out, rows = unname(lapply(seq_len(nrow(res)), 
        function(i) as.list(res[i, , drop = FALSE]))), complete_sha256 = v5_sha(file.path(out, "V5_COMPLETE.json")), 
        total_seconds = seal$total_seconds)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        v5_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
        v5_env <- function(k, default = "") {
            x <- Sys.getenv(k, "")
            if (nzchar(x)) x else default
        }
        v5_sha <- function(p) {
            v5_assert(file.exists(p) && !dir.exists(p), paste("REQUIRED_FILE_MISSING", p))
            digest::digest(file = p, algo = "sha256")
        }
        v5_hash <- function(x) digest::digest(paste0(nchar(as.character(x)), ":", as.character(x), "\n", collapse = ""), 
            algo = "sha256", serialize = FALSE)
        v5_check <- function(p, h) {
            v5_assert(is.character(h) && length(h) == 1L && grepl("^[a-f0-9]{64}$", h), "EXACT_SOURCE_SHA256_REQUIRED")
            v5_assert(identical(v5_sha(p), h), paste("FROZEN_SOURCE_SHA256_MISMATCH", p))
        }
        v5_bind <- function(xs) {
            xs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, xs)
            if (!length(xs)) return(data.frame())
            ns <- unique(unlist(lapply(xs, names), use.names = FALSE))
            do.call(rbind, lapply(xs, function(x) {
                for (k in setdiff(ns, names(x))) x[[k]] <- NA
                x[, ns, drop = FALSE]
            }))
        }
        v5_json <- function(x, p) {
            v5_assert(!file.exists(p) && !dir.exists(p), paste("REFUSE_OVERWRITE", p))
            jsonlite::write_json(x, p, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null", digits = 17)
        }
        v5_csv <- function(x, p) {
            v5_assert(!file.exists(p) && !dir.exists(p), paste("REFUSE_OVERWRITE", p))
            write.csv(x, p, row.names = FALSE, na = "NA")
        }
        v5_publish <- function(source, dest) {
            h <- v5_sha(source)
            if (file.exists(dest)) {
                v5_check(dest, h)
                return(invisible(h))
            }
            v5_assert(!dir.exists(dest), "PUBLICATION_DESTINATION_IS_DIRECTORY")
            linked <- suppressWarnings(file.link(source, dest))
            if (!linked) v5_assert(file.copy(source, dest, overwrite = FALSE, copy.date = TRUE), "ADDITIVE_PUBLICATION_FAILED_STAGING_RETAINED")
            v5_check(dest, h)
            invisible(h)
        }
        v5_deps <- function(fit = FALSE) {
            for (p in c("jsonlite", "digest", "survival", if (fit) c("brglm2", "nnet"))) v5_assert(requireNamespace(p, 
                quietly = TRUE), paste("DEPENDENCY_MISSING", p))
            if (fit) v5_assert(as.character(packageVersion("brglm2")) == "1.1.1", "PINNED_BRGLM2_1_1_1_REQUIRED")
        }
        v5_id <- function(id) {
            v5_assert(is.numeric(id) && length(id) == 1L && is.finite(id) && id == as.integer(id) && id %in% 0:100, 
                "FIXED_V5_ID_ZERO_TO_100_REQUIRED")
            as.integer(id)
        }
        v5_utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
        v5_scope <- function(root) {
            p <- file.path(root, "input/APPROVED_165_ORIGINAL.csv")
            d <- read.csv(p, colClasses = "character", check.names = FALSE, na.strings = "", fileEncoding = "UTF-8")
            v5_assert(nrow(d) == 165L && !anyNA(d$indicator_id) && !anyDuplicated(d$indicator_id) && sum(d$type == 
                "binary") == 13L && sum(d$type == "categorical") == 15L && sum(d$type == "continuous") == 137L && 
                sum(d$category == "CMR") == 44L, "EXACT_APPROVED_165_13_15_137_44_REQUIRED")
            cfg <- jsonlite::fromJSON(file.path(root, "config/runtime.json"), simplifyVector = FALSE)
            v5_assert(identical(v5_sha(p), cfg$approved_csv_sha256) && cfg$max_bootstrap_id == 100L && !isTRUE(cfg$stage2_approved), 
                "V5_SCOPE_HASH_OR_AUTHORIZATION_CHANGED")
            d
        }
        v5_sources <- function(root) {
            v5_deps()
            v9 <- normalizePath(v5_env("B3_V9_ROOT", dirname(dirname(root))), mustWork = TRUE)
            ca <- new.env(parent = globalenv())
            catroot <- file.path(v9, "extensions/categorical_cpu_gpu_v2")
            sys.source(file.path(catroot, "R/reference_adapter.R"), envir = ca)
            sys.source(file.path(catroot, "R/prepare_and_cache.R"), envir = ca)
            loaded <- ca$c40_sources(catroot)
            b18root <- v5_env("V5_B18_ROOT", file.path(dirname(v9), "border_gcomp_18x5_completion_20261024_v1"))
            sys.source(file.path(b18root, "f/1.adapter.R"), envir = ca)
            list(engine = loaded$engine, discrete = ca, source_sha256 = c(loaded$source_sha256, categorical_adapter = v5_sha(file.path(catroot, 
                "R/reference_adapter.R")), categorical_preparer = v5_sha(file.path(catroot, "R/prepare_and_cache.R")), 
                b18_adapter = v5_sha(file.path(b18root, "f/1.adapter.R"))), v9 = v9)
        }
        v5_run <- function(run, create = FALSE) {
            p <- gsub(birth_private_setting("PRIVATE_INPUT"), "/", run)
            v5_assert(nzchar(p) && !grepl("(^|/)[.][.](/|$)", p) && grepl("/border_165_typewise_b100_v5/run_[A-Za-z0-9_.-]+$", 
                p), "ISOLATED_V5_VERSIONED_RUN_REQUIRED")
            v5_assert(startsWith(p, birth_private_setting("PRIVATE_INPUT")) || grepl("^[Gg]:/UKB_results/birth/", 
                p) || grepl("PRIVATE_CONTEXT_REQUIRED", p), "G_RESULT_ROOT_REQUIRED_NO_MAC_LOCAL_CACHE")
            if (create && !dir.exists(run)) v5_assert(dir.create(run, recursive = TRUE, showWarnings = FALSE), 
                "RUN_CREATE_FAILED")
            normalizePath(run, mustWork = TRUE)
        }
        v5_rebase <- function(p) {
            r <- v5_env("V5_INPUT_ROOT")
            if (nzchar(r) && startsWith(p, birth_private_setting("PRIVATE_INPUT"))) file.path(r, substring(p, nchar(birth_private_setting("PRIVATE_INPUT")) + 
                1L)) else p
        }
        v5_read_source <- function(desc) {
            v5_assert(is.list(desc) && !is.null(desc$path) && !is.null(desc$sha256), "SOURCE_DESCRIPTOR_PATH_SHA_REQUIRED")
            p <- v5_rebase(desc$path)
            v5_check(p, desc$sha256)
            z <- if (identical(desc$format, "json")) jsonlite::fromJSON(p, simplifyVector = FALSE) else if (identical(desc$format, 
                "csv")) {
                v5_assert(requireNamespace("data.table", quietly = TRUE), "DATA_TABLE_REQUIRED")
                as.data.frame(data.table::fread(p, select = if (is.null(desc$columns)) NULL else unlist(desc$columns), 
                  showProgress = FALSE))
            } else readRDS(p)
            if (isTRUE(desc$unwrap)) while (is.list(z) && !is.null(z$value)) {
                if (!is.null(z$value_hash)) v5_assert(identical(digest::digest(z$value, algo = "sha256"), z$value_hash), 
                  "WRAPPED_SOURCE_CONTENT_CHANGED")
                z <- z$value
            }
            if (!is.null(desc$data_path)) for (k in unlist(desc$data_path, use.names = FALSE)) {
                v5_assert(!is.null(z[[k]]), paste("SOURCE_DATA_PATH_MISSING", k))
                z <- z[[k]]
            }
            v5_check(p, desc$sha256)
            z
        }
        v5_join <- function(data, ids, require_all = TRUE) {
            v5_assert(is.data.frame(data) && "eid" %in% names(data), "SOURCE_DATA_EID_REQUIRED")
            sid <- as.character(data$eid)
            v5_assert(!anyNA(sid) && !anyDuplicated(sid), "SOURCE_IDS_NOT_UNIQUE")
            j <- match(ids, sid)
            if (require_all) v5_assert(!anyNA(j), "SOURCE_NOT_COVERING_FROZEN_PARENT")
            j
        }
        v5_vendor_features <- function(spec) {
            v5_assert(requireNamespace("data.table", quietly = TRUE), "DATA_TABLE_REQUIRED")
            v5_assert(is.list(spec) && all(c("root", "out", "core", "files") %in% names(spec)), "VERIFIED_VENDOR_FEATURE_SOURCE_SPEC_REQUIRED")
            pins <- spec$files
            for (p in pins) v5_check(v5_rebase(p$path), p$sha256)
            e <- new.env(parent = globalenv())
            e$root <- v5_rebase(spec$root)
            e$out <- v5_rebase(spec$out)
            e$feature_cache <- new.env(parent = emptyenv())
            e$num <- function(x) suppressWarnings(as.numeric(as.character(x)))
            e$assert <- v5_assert
            e$fread <- data.table::fread
            e$read_json <- jsonlite::read_json
            e$cp_load <- function(path, strict = TRUE) {
                pp <- Filter(function(p) identical(normalizePath(v5_rebase(p$path), mustWork = TRUE), normalizePath(path, 
                  mustWork = TRUE)), pins)
                v5_assert(length(pp) == 1L, paste("UNPINNED_VENDOR_CACHE", path))
                v5_read_source(c(pp[[1]], list(unwrap = TRUE)))
            }
            e$core <- v5_read_source(spec$core)
            dictfile <- file.path(e$root, "vendor/reused/8.clinical_omics.R")
            allowed <- "v2_feature_dictionary"
            for (expr in parse(file = dictfile)) if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) && 
                as.character(expr[[2]]) %in% allowed) eval(expr, envir = e)
            v5_assert(is.function(e$v2_feature_dictionary), "FROZEN_FEATURE_DICTIONARY_FUNCTION_REQUIRED")
            e$dict <- function() e$v2_feature_dictionary(e$root)
            sys.source(file.path(e$root, "vendor/features.R"), envir = e)
            v5_assert(is.function(e$load_feature), "FROZEN_FEATURE_LOADER_REQUIRED")
            e$dictionary <- e$dict()
            e$verify <- function() {
                for (p in pins) v5_check(v5_rebase(p$path), p$sha256)
                invisible(TRUE)
            }
            e
        }
        v5_pure_risk <- function(v9) {
            rp <- file.path(dirname(v9), "border_population_divergence_explanation_v1/f/1.prepare.R")
            v5_check(rp, birth_private_setting("PRIVATE_HASH_13"))
            pure <- new.env(parent = globalenv())
            pure$pd_assert <- v5_assert
            pure$pd_hash <- function(x) digest::digest(x, algo = "sha256")
            sys.source(rp, envir = pure)
            pure
        }
        v5_dated_bridge <- function(parent, legacy, b18, pure) {
            ids <- as.character(b18$d_full$eid)
            j <- v5_join(parent$data, ids)
            p <- parent$data[j, , drop = FALSE]
            lj <- v5_join(legacy$data, ids)
            v5_assert(all(c("date_attend", "date_death_uncapped", "date_lost_registry", "record_observation_eligible", 
                "P_common") %in% names(p)) && all(!is.na(p$P_common) & p$P_common) && all(!is.na(p$record_observation_eligible) & 
                p$record_observation_eligible), "FROZEN_COMMON_TARGET_OBSERVATION_PROOF_REQUIRED")
            v5_assert(is.character(parent$baseline_note) && grepl("baseline union|BASELINE.*UNION", parent$baseline_note, 
                ignore.case = TRUE), "FROZEN_ALL_E1_BASELINE_FREE_TARGET_EVIDENCE_REQUIRED")
            ans <- data.frame(eid = ids, record_observation_eligible = p$record_observation_eligible)
            audit <- list()
            for (e in c("250", "250.2", "411", "411.3", "411.4")) {
                if (e %in% c("250", "411.3")) {
                  v5_assert(e %in% names(legacy$events) && e %in% names(legacy$baseline) && inherits(legacy$events[[e]], 
                    "Date"), "LEGACY_DM_AG_ACTUAL_DATES_REQUIRED")
                  event <- legacy$events[[e]][lj]
                  st <- as.integer(legacy$baseline[[e]][lj])
                  source <- "EXACT_LEGACY_E0_EVENT_DATE_VECTOR"
                } else {
                  r <- parent$risk$E1[[e]]
                  v5_assert(is.data.frame(r) && identical(as.character(r$eid), as.character(parent$data$eid)) && 
                    all(b18$endpoints[[e]]$eligibility), "FROZEN_E1_FULL_ELIGIBLE_PARENT_RISK_REQUIRED")
                  r <- r[j, , drop = FALSE]
                  days <- r$exit * 365.25
                  v5_assert(all(is.finite(days)) && max(abs(days - round(days))) < 1e-07, "FROZEN_E1_EXIT_NOT_EXACT_INTEGER_DAY_NO_DATE_INFERENCE")
                  event <- as.Date(rep(NA_real_, length(ids)), origin = "1970-01-01")
                  ev <- r$event == 1L
                  event[ev] <- as.Date(p$date_attend[ev]) + round(days[ev])
                  st <- rep(0L, length(ids))
                  source <- "EXACT_INTEGER_DAY_RECONSTRUCTION_OF_FROZEN_E1_OBSERVED_EVENT_DATE"
                }
                rr <- pure$pd_data_risk(p, event, st, p$record_observation_eligible, b18$cfg$admin)
                ep <- b18$endpoints[[e]]
                ok <- identical(rr$eligible, ep$eligibility) && isTRUE(all.equal(rr$data, ep$risk, tolerance = 1e-12, 
                  check.attributes = FALSE))
                v5_assert(ok, paste("DATED_BRIDGE_DOES_NOT_RECONSTRUCT_FROZEN_BASELINE_ENDPOINT", e))
                key <- gsub("[.]", "_", e)
                ans[[paste0("event__", key)]] <- event
                ans[[paste0("status__", key)]] <- st
                audit[[e]] <- data.frame(endpoint_id = e, source = source, baseline_risk_reconstructed_identically = ok, 
                  N_parent = length(ids), N_dated_events = sum(!is.na(event)), endpoint_version = ep$endpoint_version, 
                  admin = b18$cfg$admin, stringsAsFactors = FALSE)
            }
            list(data = ans, audit = v5_bind(audit))
        }
        v5_recover_social <- function(spec, loaded, ids, old, baseline_date) {
            v5_assert(all(c("raw", "dictionary", "encoder") %in% names(spec)), "EXACT_SOCIAL_RECOVERY_PINS_REQUIRED")
            raw <- v5_read_source(spec$raw)
            defs <- v5_read_source(spec$dictionary)
            v5_check(v5_rebase(spec$encoder$path), spec$encoder$sha256)
            encoder <- loaded$discrete$c40p_encoder(v5_rebase(spec$encoder$path))
            j <- v5_join(raw, ids)
            oj <- v5_join(old, ids, FALSE)
            v5_assert("p53_i0" %in% names(raw) && all(!is.na(as.Date(raw$p53_i0[j])) & as.Date(raw$p53_i0[j]) == 
                as.Date(baseline_date)), "SOCIAL_RAW_ACTUAL_BASELINE_DATES_CONFLICT")
            out <- data.frame(eid = ids)
            audit <- list()
            targets <- loaded$discrete$c40p_targets(raw[j, , drop = FALSE], encoder)
            for (id in unlist(spec$indicators, use.names = FALSE)) {
                def <- Filter(function(x) identical(x$id, id), defs)
                v5_assert(length(def) == 1L && identical(def[[1]]$type, "continuous"), "EXACT_FROZEN_SOCIAL_CONTINUOUS_DEFINITION_REQUIRED")
                def <- def[[1]]
                v5_assert(def$field %in% names(raw), "SOCIAL_RAW_FIELD_MISSING")
                x <- encoder$encode_marker(raw[j, , drop = FALSE], def)$value
                v5_assert(is.numeric(x) && length(x) == length(ids), "FROZEN_SOCIAL_ENCODER_NOT_NUMERIC")
                overlap <- !is.na(oj)
                v5_assert(id %in% names(old) && isTRUE(all.equal(x[overlap], old[[id]][oj[overlap]], tolerance = 0, 
                  check.attributes = FALSE)), paste("FROZEN_SOCIAL_RECOVERY_OVERLAP_CHANGED", id))
                out[[id]] <- x
                tg <- targets[[def$eligible]]
                v5_assert(is.logical(tg) && !anyNA(tg), "FROZEN_SOCIAL_APPLICABILITY_TARGET_REQUIRED")
                out[[paste0("eligible__", id)]] <- tg
                audit[[id]] <- data.frame(indicator_id = id, N_parent = length(ids), N_target = sum(tg), target_rule = def$eligible, 
                  N_old_cache_overlap = sum(overlap), N_new_parent_rows_recovered = sum(!overlap & !is.na(x)), 
                  N_valid = sum(!is.na(x)), overlap_values_exact = TRUE, source = "EXISTING_FROZEN_CONTINUOUS_ENCODER_FULL_PARENT_RAW", 
                  stringsAsFactors = FALSE)
            }
            list(data = out, audit = v5_bind(audit))
        }
        v5_stage <- function(ctx, stage, endpoint = NULL, detail = NULL) {
            m <- ctx$monitor
            if (!is.environment(m)) return(invisible(NULL))
            now <- as.numeric(Sys.time())
            mono <- proc.time()[["elapsed"]]
            if (is.null(endpoint)) endpoint <- if (is.null(m$endpoint_id)) "SHARED" else m$endpoint_id
            m$stage <- stage
            m$endpoint_id <- endpoint
            m$stage_started_epoch <- now
            m$stage_started_monotonic <- mono
            m$last_progress_epoch <- now
            m$sequence <- if (is.null(m$sequence)) 1L else m$sequence + 1L
            event <- list(schema = "V5_WORKER_STAGE_V1", run_id = ctx$run_id, worker_id = v5_env("V5_WORKER_SLOT"), 
                host = v5_env("V5_HOST"), pid = Sys.getpid(), process_start = m$process_start, indicator_id = m$indicator_id, 
                bootstrap_id = m$bootstrap_id, attempt_id = m$attempt_id, endpoint_id = endpoint, stage = stage, 
                stage_detail = detail, stage_started_utc = v5_utc(), stage_started_epoch = now, stage_started_monotonic = mono, 
                stage_elapsed_s = 0, package_started_epoch = m$package_started_epoch, package_elapsed_s = if (is.null(m$package_started_monotonic)) NA_real_ else mono - 
                  m$package_started_monotonic, last_progress_epoch = now, sequence = m$sequence, endpoint_status = m$endpoint_status)
            ev <- v5_env("V5_WORKER_EVENTS")
            if (nzchar(ev)) {
                dir.create(dirname(ev), recursive = TRUE, showWarnings = FALSE)
                cat(jsonlite::toJSON(event, auto_unbox = TRUE, na = "null", null = "null", digits = 17), "\n", 
                  file = ev, append = TRUE, sep = "")
            }
            sp <- v5_env("V5_WORKER_STATE")
            if (nzchar(sp)) {
                dir.create(dirname(sp), recursive = TRUE, showWarnings = FALSE)
                tmp <- paste0(sp, ".", Sys.getpid(), ".pending")
                jsonlite::write_json(event, tmp, auto_unbox = TRUE, na = "null", null = "null", digits = 17)
                v5_assert(file.rename(tmp, sp), "WORKER_STATE_ATOMIC_PUBLICATION_FAILED")
            }
            invisible(event)
        }
        v5_frame_endpoint <- function(ctx) {
            fs <- sys.frames()
            for (i in rev(seq_along(fs))) if (exists("e", envir = fs[[i]], inherits = FALSE)) {
                e <- get("e", envir = fs[[i]], inherits = FALSE)
                if (is.character(e) && length(e) == 1L && e %in% ctx$cfg$outcomes) return(e)
            }
            if (is.environment(ctx$monitor) && !is.null(ctx$monitor$endpoint_id)) ctx$monitor$endpoint_id else "SHARED"
        }
        v5_install_hooks <- function(ctx) {
            wrap <- function(env, name, stage) {
                original <- env[[name]]
                v5_assert(is.function(original), paste("FROZEN_FUNCTION_REQUIRED", name))
                env[[name]] <- local({
                  fn <- original
                  st <- stage
                  nm <- name
                  function(...) {
                    ep <- v5_frame_endpoint(ctx)
                    v5_stage(ctx, st, ep)
                    value <- fn(...)
                    if (nm == "b3_missing" && is.environment(ctx$monitor)) {
                      aa <- list(...)
                      d <- aa[[1]]
                      pr <- match(as.character(d$eid), as.character(ctx$d_full$eid))
                      if (!anyNA(pr)) {
                        key <- v5_hash(pr)
                        ctx$monitor$selection_weight_hashes[[key]] <- list(selection_weight_hash = v5_hash(format(value$w, 
                          digits = 17, scientific = TRUE)), selection_sample_hash = v5_hash(pr[value$complete]), 
                          selection_target_hash = key)
                      }
                    }
                    if (nm == "c40_exact_integral" && ep %in% ctx$cfg$outcomes) ctx$monitor$endpoint_status[[ep]] <- "COMPUTED_NOT_SAVED"
                    value
                  }
                })
            }
            hooks <- c(b3_missing = "SELECTION_WEIGHT", b3_censor = "CENSOR_MODEL", b3_logit = "OUTCOME_FIT", b3_predict_scenarios = "STANDARDIZE")
            for (k in names(hooks)) wrap(ctx$engine, k, hooks[[k]])
            wrap(ctx$discrete, "c40_discrete_fit", "MEDIATOR_FIT")
            wrap(ctx$discrete, "c40_exact_integral", "STANDARDIZE")
            ctx
        }
        v5_landmark <- function(ids, base, visit, event_date, status, observed, admin, pure, endpoint_version) {
            required <- c("date", "age", "center", "date_death_uncapped", "date_lost_registry")
            v5_assert(all(required %in% names(visit)) && nrow(visit) == length(ids), "ACTUAL_IMAGING_VISIT_DATA_REQUIRED")
            d <- base
            origin <- as.Date(visit$date)
            v5_assert(inherits(event_date, "Date"), "ACTUAL_FROZEN_EVENT_DATES_REQUIRED")
            d$date_attend <- origin
            d$date_death_uncapped <- as.Date(visit$date_death_uncapped)
            d$date_lost_registry <- as.Date(visit$date_lost_registry)
            d$age <- as.numeric(visit$age)
            d$center <- factor(visit$center)
            d$assessment_year <- factor(format(origin, "%Y"))
            r <- pure$pd_data_risk(d, event_date, status, observed, admin)
            eligible <- r$eligible & is.finite(d$age) & !is.na(d$center) & !is.na(d$assessment_year)
            known <- r$data$exit >= 5 | (r$data$exit <= 5 & (r$data$event == 1 | r$data$death_event == 1))
            support <- any(eligible & is.finite(r$data$exit) & r$data$exit >= 5)
            list(data = d, endpoint = list(status = if (support) "AVAILABLE_FROZEN" else "HOLD", reason = if (support) "" else "NO_FIVE_YEAR_LANDMARK_FOLLOWUP_SUPPORT_ORIGINAL_ADMIN_UNCHANGED", 
                risk = r$data, eligibility = eligible, endpoint_version = paste0(endpoint_version, "_INSTANCE2_LANDMARK"), 
                risk_hash = v5_hash(c(endpoint_version, as.character(origin), as.character(event_date), as.character(admin)))), 
                audit = data.frame(endpoint_version = endpoint_version, time_origin = "ACTUAL_INSTANCE2_DATE", 
                  admin = admin, N_parent = length(ids), N_eligible = sum(eligible), N_known_5y = sum(eligible & 
                    known, na.rm = TRUE), N_followup_ge5 = sum(eligible & r$data$exit >= 5, na.rm = TRUE), N_prevalent_at_imaging = sum(!is.na(event_date) & 
                    event_date <= origin, na.rm = TRUE), N_events5 = sum(eligible & r$data$event == 1 & r$data$exit <= 
                    5, na.rm = TRUE), status = if (support) "AVAILABLE_FROZEN" else "HOLD"))
        }
        v5_prepare <- function(root, run = v5_env("V5_RUN_DIR")) {
            v5_assert(identical(v5_env("V5_REAL_EXECUTION_AUTHORIZED"), "YES_FOR_THIS_RUN"), "EXPLICIT_V5_EXECUTION_AUTHORIZATION_REQUIRED")
            v5_assert(Sys.info()[["sysname"]] != "Darwin" && v5_env("V5_HOST") == "winpc", "PARTICIPANT_PREPARATION_WINPC_ONLY")
            root <- normalizePath(root, mustWork = TRUE)
            run <- v5_run(run, TRUE)
            v5_deps()
            scope <- v5_scope(root)
            loaded <- v5_sources(root)
            sp <- file.path(root, "config/adapter_sources.json")
            v5_assert(file.exists(sp), "VERIFIED_ADAPTER_SOURCES_JSON_REQUIRED")
            src <- jsonlite::fromJSON(sp, simplifyVector = FALSE)
            final_private <- file.path(run, "private")
            if (!dir.exists(final_private)) v5_assert(dir.create(final_private), "PRIVATE_ROOT_CREATE_FAILED")
            v5_assert(!file.exists(file.path(final_private, "V5_INPUT_LOCK.json")), "PRESERVE_EXISTING_COMPLETE_PREPARATION")
            private <- file.path(final_private, paste0("preparation_attempt_", format(Sys.time(), "%Y%m%dT%H%M%S", 
                tz = "UTC"), "_", Sys.getpid()))
            v5_assert(!dir.exists(private) && dir.create(private), "UNIQUE_PREPARATION_ATTEMPT_REQUIRED")
            social_path <- file.path(root, "config/social_recovery.json")
            social_sha <- NULL
            if (file.exists(social_path)) {
                social_sha <- v5_sha(social_path)
                src$social_recovery <- jsonlite::fromJSON(social_path, simplifyVector = FALSE)
                src$social_recovery$dictionary$format <- "json"
                src$social_recovery$indicators <- names(Filter(function(x) identical(x$source, "social"), src$features))
            }
            b18 <- v5_read_source(src$b18_cache)
            catctx <- v5_read_source(src$categorical_cache)
            ids <- as.character(b18$d_full$eid)
            v5_assert(length(ids) == 123350L && identical(ids, as.character(catctx$d_full$eid)), "EXACT_PARENT_AND_CATEGORICAL_BRIDGE_REQUIRED")
            for (k in c("parent_frame_hash", "cluster_order_hash", "parent_cluster_membership_hash")) v5_assert(identical(b18[[k]], 
                catctx[[k]]), paste("PARENT_IDENTITY_CONFLICT", k))
            cfg <- b18$cfg
            cfg$outcomes <- c("250", "250.2", "411", "411.3", "411.4")
            cfg$base_covariates <- as.character(unlist(cfg$base_covariates))
            v5_assert(identical(cfg$base_covariates, c("B", "age", "sex", "ethnicity", "center", "assessment_year", 
                "S")) && cfg$horizon == 5 && cfg$admin == "2022-05-31" && cfg$seed == 20261023, "FROZEN_STATISTICAL_CONTRACT_CHANGED")
            d <- b18$d_full
            features <- variables <- eligibility <- list()
            metadata <- list()
            sources <- list()
            if (!is.null(src$sources)) for (k in names(src$sources)) sources[[k]] <- v5_read_source(src$sources[[k]])
            vendor <- if (!is.null(src$vendor_features)) v5_vendor_features(src$vendor_features) else NULL
            if (!is.null(vendor)) {
                md <- as.data.frame(vendor$core$metadata)
                cd <- as.data.frame(vendor$core$data)
                mj <- v5_join(md, as.character(cd$eid))
                md <- md[mj, , drop = FALSE]
                for (k in intersect(c("date_attend", "date_death_uncapped", "date_lost_registry"), names(cd))) md[[k]] <- cd[[k]]
                sources$vendor_visit <- md
                sources$vendor_core_data <- cd
            }
            social_audit <- NULL
            if (!is.null(src$social_recovery)) {
                v5_assert(!is.null(sources$social) && !is.null(sources$vendor_core_data), "SOCIAL_CACHE_AND_PARENT_DATES_REQUIRED")
                cd <- sources$vendor_core_data
                bj <- v5_join(cd, ids)
                recovered <- v5_recover_social(src$social_recovery, loaded, ids, sources$social, cd$date_attend[bj])
                sources$social_recovered <- recovered$data
                social_audit <- recovered$audit
                for (id in unlist(src$social_recovery$indicators, use.names = FALSE)) {
                  v5_assert(identical(src$features[[id]]$source, "social"), "SOCIAL_RECOVERY_ONLY_OLD_SOCIAL_FEATURES")
                  src$features[[id]]$source <- "social_recovered"
                  src$features[[id]]$target_source <- "social_recovered"
                  src$features[[id]]$target_column <- paste0("eligible__", id)
                }
            }
            dated_audit <- NULL
            if (!is.null(src$dated_frozen_bridge)) {
                db <- src$dated_frozen_bridge
                dp <- v5_read_source(db$parent_cache)
                dl <- v5_read_source(db$legacy_core)
                bridge <- v5_dated_bridge(dp, dl, b18, v5_pure_risk(loaded$v9))
                sources$dated_bridge <- bridge$data
                dated_audit <- bridge$audit
                if (is.null(src$dated_endpoints)) {
                  src$dated_endpoints <- setNames(lapply(cfg$outcomes, function(e) {
                    key <- gsub("[.]", "_", e)
                    list(source = "dated_bridge", event_column = paste0("event__", key), baseline_column = paste0("status__", 
                      key), observation_column = "record_observation_eligible")
                  }), cfg$outcomes)
                }
            }
            profile <- if (!is.null(src$profile_cache)) v5_read_source(src$profile_cache) else NULL
            for (i in seq_len(nrow(scope))) {
                row <- scope[i, , drop = FALSE]
                id <- row$indicator_id
                is_cmr <- row$category == "CMR"
                spec <- src$features[[id]]
                if (row$type != "continuous") {
                  ca <- Filter(function(x) identical(x$indicator_id, id), catctx$cfg$indicators)
                  v5_assert(length(ca) == 1L && identical(ca[[1]]$kind, row$type), "FROZEN_CATEGORICAL_ID_TYPE_MISMATCH")
                  ca <- ca[[1]]
                  d[[ca$column]] <- catctx$d_full[[ca$column]]
                  variables[[id]] <- catctx$variables[[id]]
                  eligibility[[id]] <- catctx$eligibility[[id]]
                  features[[id]] <- ca
                } else {
                  old <- Filter(function(x) identical(x$feature_id, id), b18$cfg$features)
                  pm <- if (!is.null(profile)) profile$feature_metadata[profile$feature_metadata$feature_id == 
                    id, , drop = FALSE] else data.frame()
                  if (length(old) == 1L) {
                    ca <- old[[1]]
                    x <- d[[ca$column]]
                    target <- b18$eligibility[[ca$target_key]]
                    ca$inverse_transform <- if (ca$prepared_transform == "log") "exp" else "identity"
                    status <- "AVAILABLE_FROZEN"
                    reason <- ""
                  } else if (!is.null(vendor) && id %in% vendor$dictionary$mediator_id) {
                    vr <- as.list(vendor$dictionary[vendor$dictionary$mediator_id == id, ])
                    v5_assert(length(vr$mediator_id) == 1L, "UNIQUE_EXACT_VENDOR_FEATURE_REQUIRED")
                    vf <- vendor$load_feature(vendor$core, vr)
                    j <- v5_join(vf$data, ids)
                    x <- vf$data$marker[j]
                    fitted_unit <- if (isTRUE(vr$indexed)) paste0(vr$unit, "/m2") else vr$unit
                    if (identical(vr$transform, "log")) fitted_unit <- paste0("log(", fitted_unit, ")")
                    ca <- list(feature_id = id, candidate_id = id, column = paste0("v5__", make.names(id)), platform = row$category, 
                      target_key = id, raw_unit = vr$unit, prepared_unit = fitted_unit, prepared_transform = vr$transform, 
                      preprocessing = if (isTRUE(vr$indexed)) "same_instance2_BSA_division" else "FROZEN_FEATURE_LOADER_QC", 
                      inverse_transform = if (vr$transform == "log") "exp" else "identity", original_extra_covariates = as.character(vf$extra), 
                      source_qc = vf$source_qc, source_column = vr$column)
                    target <- if (startsWith(id, "nmr__")) b18$eligibility$nmr_measured_i0 else rep(TRUE, length(ids))
                    status <- "AVAILABLE_FROZEN"
                    reason <- ""
                    if (is.null(spec) || !isTRUE(spec$physical_support_verified)) {
                      status <- "HOLD"
                      reason <- "PREPARED_VALUES_AVAILABLE_PHYSICAL_SUPPORT_SPEC_NOT_YET_VERIFIED"
                    } else for (k in c("physical_lower_bound", "physical_upper_bound")) if (!is.null(spec[[k]])) ca[[k]] <- as.numeric(spec[[k]])
                    ca$technical_covariates <- character()
                    if (!is.null(spec$m3_technical_covariates)) for (k in unlist(spec$m3_technical_covariates)) {
                      v5_assert(k %in% vf$extra && k %in% names(vf$data), "TECHNICAL_COVARIATE_NOT_IN_FROZEN_LOADER")
                      nm <- paste0("tech__", make.names(id), "__", k)
                      d[[nm]] <- vf$data[[k]][j]
                      ca$technical_covariates <- c(ca$technical_covariates, nm)
                    }
                    if (!is_cmr) {
                      v5_assert(all(is.na(vf$data$measurement_date[j]) | as.Date(vf$data$measurement_date[j]) == 
                        as.Date(sources$vendor_core_data$date_attend[match(ids, as.character(sources$vendor_core_data$eid))])), 
                        "BASELINE_MEASUREMENT_DATE_CONFLICT")
                    }
                  } else if (!is.null(spec)) {
                    v5_assert(all(c("source", "column", "prepared_unit", "prepared_transform", "inverse_transform") %in% 
                      names(spec)), paste("EXPLICIT_SCALE_SPEC_REQUIRED", id))
                    sd <- sources[[spec$source]]
                    v5_assert(is.data.frame(sd) && spec$column %in% names(sd), paste("EXPLICIT_SOURCE_COLUMN_MISSING", 
                      id))
                    j <- v5_join(sd, ids, FALSE)
                    x <- sd[[spec$column]][j]
                    v5_assert(is.numeric(x), "PREPARED_CONTINUOUS_NUMERIC_REQUIRED_NO_REENCODING")
                    ca <- c(list(feature_id = id, candidate_id = id, column = paste0("v5__", make.names(id)), platform = row$category, 
                      target_key = id), spec)
                    ca$column <- paste0("v5__", make.names(id))
                    target <- rep(TRUE, length(ids))
                    if (!is.null(spec$target_source)) {
                      td <- sources[[spec$target_source]]
                      tj <- v5_join(td, ids)
                      target <- td[[spec$target_column]][tj]
                      v5_assert(is.logical(target) && !anyNA(target), "EXPLICIT_TARGET_LOGICAL_REQUIRED")
                    }
                    status <- "AVAILABLE_FROZEN"
                    reason <- ""
                  } else if (nrow(pm) == 1L) {
                    ca <- list(feature_id = id, candidate_id = id, column = paste0("v5__", make.names(id)), platform = row$category, 
                      target_key = id, prepared_unit = pm$unit, prepared_transform = pm$transform, inverse_transform = if (pm$transform == 
                        "log") "exp" else "identity")
                    j <- v5_join(profile$data, ids)
                    x <- profile$data[[pm$column]][j]
                    target <- if (startsWith(id, "nmr__")) profile$eligibility$nmr_measured_i0[j] else rep(TRUE, 
                      length(ids))
                    status <- "HOLD"
                    reason <- "PREPARED_VALUES_AVAILABLE_PHYSICAL_SUPPORT_SPEC_NOT_YET_VERIFIED"
                  } else {
                    ca <- list(feature_id = id, candidate_id = id, column = paste0("v5__", make.names(id)), platform = row$category, 
                      target_key = id, prepared_unit = NA_character_, prepared_transform = NA_character_, inverse_transform = NA_character_)
                    x <- rep(NA_real_, length(ids))
                    target <- rep(TRUE, length(ids))
                    status <- "HOLD"
                    reason <- "EXACT_PREPARED_SOURCE_COLUMN_UNIT_TRANSFORM_NOT_YET_VERIFIED"
                  }
                  if (!is.null(spec) && length(old) == 1L) for (k in c("physical_lower_bound", "physical_upper_bound")) if (!is.null(spec[[k]])) v5_assert(identical(as.numeric(spec[[k]]), 
                    as.numeric(ca[[k]])), paste("ORIGINAL_FROZEN_BOUND_CHANGED", id, k))
                  d[[ca$column]] <- x
                  features[[id]] <- ca
                  eligibility[[id]] <- target
                  variables[[id]] <- list(status = status, reason = reason, measurement_instance = row$measurement_instance, 
                    source_label_retained = row$source_label_retained)
                }
                ca <- features[[id]]
                metadata[[id]] <- data.frame(indicator_id = id, feature_id = id, type = row$type, domain = row$category, 
                  time_origin = if (is_cmr) "imaging" else "baseline", column = ca$column, raw_unit = if (is.null(ca$raw_unit)) if (is.null(ca$prepared_unit)) "category" else ca$prepared_unit else ca$raw_unit, 
                  prepared_unit = if (is.null(ca$prepared_unit)) "category" else ca$prepared_unit, preprocessing = if (is.null(ca$preprocessing)) "INHERITED_PREPARED_VALUES" else ca$preprocessing, 
                  prepared_transform = if (is.null(ca$prepared_transform)) "FROZEN_FACTOR" else ca$prepared_transform, 
                  inverse_transform = if (is.null(ca$inverse_transform)) "NOT_APPLICABLE" else ca$inverse_transform, 
                  status = variables[[id]]$status, reason = if (is.null(variables[[id]]$reason)) "" else variables[[id]]$reason, 
                  N_target = sum(eligibility[[id]]), N_observed = sum(eligibility[[id]] & !is.na(d[[ca$column]])), 
                  stringsAsFactors = FALSE)
            }
            imaging <- list(status = "HOLD", reason = "ACTUAL_INSTANCE2_VISIT_AND_FROZEN_DATED_ENDPOINTS_NOT_YET_VERIFIED")
            imaging_audit <- list()
            if (!is.null(src$imaging)) {
                im <- src$imaging
                vd <- sources[[im$source]]
                j <- v5_join(vd, ids, FALSE)
                cols <- c(date = im$date_column, age = im$age_column, center = im$center_column, date_death_uncapped = im$death_column, 
                  date_lost_registry = im$loss_column)
                v5_assert(all(cols %in% names(vd)), "EXPLICIT_IMAGING_VISIT_COLUMNS_MISSING")
                visit <- vd[j, cols, drop = FALSE]
                names(visit) <- names(cols)
                pure <- v5_pure_risk(loaded$v9)
                ieps <- list()
                imd <- NULL
                for (e in cfg$outcomes) {
                  es <- src$dated_endpoints[[e]]
                  v5_assert(!is.null(es), "FIVE_FROZEN_DATED_ENDPOINT_SPECS_REQUIRED")
                  ed <- sources[[es$source]]
                  ej <- v5_join(ed, ids)
                  z <- v5_landmark(ids, d, visit, as.Date(ed[[es$event_column]][ej]), as.integer(ed[[es$baseline_column]][ej]), 
                    as.logical(ed[[es$observation_column]][ej]), cfg$admin, pure, b18$endpoints[[e]]$endpoint_version)
                  ieps[[e]] <- z$endpoint
                  imaging_audit[[e]] <- z$audit
                  imd <- z$data
                }
                for (e in c("250.2", "411")) {
                  r <- ieps[[e]]$risk
                  k <- gsub("[.]", "_", e)
                  imd[[paste0("aux_t_", k)]] <- pmin(r$exit, 5)
                  imd[[paste0("aux_event_", k)]] <- as.integer(r$event == 1 & r$exit <= 5)
                  imd[[paste0("aux_death_", k)]] <- as.integer(r$death_event == 1 & r$exit <= 5)
                }
                imaging <- list(status = "AVAILABLE_FROZEN", data = imd, endpoints = ieps, reason = "", covariate_time = "ACTUAL_INSTANCE2_AGE_YEAR_CENTER", 
                  admin = cfg$admin)
            }
            for (id in scope$indicator_id[scope$category == "CMR"]) if (imaging$status != "AVAILABLE_FROZEN") {
                variables[[id]]$status <- "HOLD"
                variables[[id]]$reason <- imaging$reason
                metadata[[id]]$status <- "HOLD"
                metadata[[id]]$reason <- imaging$reason
            }
            if (!is.null(vendor)) vendor$verify()
            cfg$features <- unname(features[scope$indicator_id[scope$type == "continuous"]])
            cfg$candidates <- cfg$features
            cfg$indicators <- unname(features[scope$indicator_id[scope$type != "continuous"]])
            cfg$outcomes <- c("250", "250.2", "411", "411.3", "411.4")
            ctx <- list(cfg = cfg, d_full = d, features = features, variables = variables, eligibility = eligibility, 
                endpoints = b18$endpoints, imaging = imaging, scope = scope, run_id = jsonlite::fromJSON(file.path(root, 
                  "config/runtime.json"))$run_id, parent_frame_hash = b18$parent_frame_hash, cluster_order_hash = b18$cluster_order_hash, 
                parent_cluster_membership_hash = b18$parent_cluster_membership_hash, scope_hash = v5_sha(file.path(root, 
                  "input/APPROVED_165_ORIGINAL.csv")), method_hash = v5_hash(loaded$source_sha256), source_hash = v5_hash(vapply(sort(list.files(file.path(root, 
                  "f"), full.names = TRUE, pattern = "[.]R$")), v5_sha, "")))
            cache <- file.path(private, "V5_INPUTS.private.rds")
            saveRDS(ctx, cache, compress = FALSE)
            v5_csv(v5_bind(metadata), file.path(private, "FEATURE_INPUT_PREFLIGHT.csv"))
            if (length(imaging_audit)) v5_csv(v5_bind(imaging_audit), file.path(private, "CMR_LANDMARK_PREFLIGHT.csv"))
            if (!is.null(dated_audit)) v5_csv(dated_audit, file.path(private, "DATED_ENDPOINT_RECONSTRUCTION_PREFLIGHT.csv"))
            if (!is.null(social_audit)) v5_csv(social_audit, file.path(private, "SOCIAL_FULL_PARENT_RECOVERY_PREFLIGHT.csv"))
            lock <- list(complete = TRUE, models_run = 0L, run_id = ctx$run_id, prepared_inputs_sha256 = v5_sha(cache), 
                approved_csv_sha256 = ctx$scope_hash, adapter_sources_sha256 = v5_sha(sp), social_recovery_sha256 = social_sha, 
                parent_source_sha256 = as.list(loaded$source_sha256), adapter_source_hash = ctx$source_hash, method_hash = ctx$method_hash, 
                parent_frame_hash = ctx$parent_frame_hash, cluster_order_hash = ctx$cluster_order_hash, parent_cluster_membership_hash = ctx$parent_cluster_membership_hash, 
                cache_sha256 = b18$cfg$cache_sha256, cohort_sha256 = b18$cfg$cohort_sha256, N_parent = length(ids), 
                N_indicators = 165L, N_continuous = 137L, N_discrete = 28L, participant_data_location = "G_ONLY_NO_MAC_BULK_COPY", 
                prepared_bytes = as.numeric(file.info(cache)$size), old_outputs_modified = FALSE)
            v5_json(lock, file.path(private, "V5_INPUT_LOCK.json"))
            files <- setdiff(list.files(private, full.names = TRUE), file.path(private, "V5_INPUT_LOCK.json"))
            for (p in files) v5_publish(p, file.path(final_private, basename(p)))
            v5_publish(file.path(private, "V5_INPUT_LOCK.json"), file.path(final_private, "V5_INPUT_LOCK.json"))
            lock
        }
        v5_load <- function(root, run = v5_env("V5_RUN_DIR")) {
            v5_assert(identical(v5_env("V5_REAL_EXECUTION_AUTHORIZED"), "YES_FOR_THIS_RUN"), "EXPLICIT_V5_EXECUTION_AUTHORIZATION_REQUIRED")
            v5_deps(TRUE)
            root <- normalizePath(root, mustWork = TRUE)
            run <- v5_run(run)
            scope <- v5_scope(root)
            loaded <- v5_sources(root)
            lp <- file.path(run, "private/V5_INPUT_LOCK.json")
            lock <- jsonlite::fromJSON(lp, simplifyVector = FALSE)
            v5_assert(isTRUE(lock$complete) && lock$models_run == 0L, "PREPARED_V5_INPUT_LOCK_REQUIRED")
            cache <- file.path(run, "private/V5_INPUTS.private.rds")
            v5_check(cache, lock$prepared_inputs_sha256)
            v5_assert(identical(v5_hash(vapply(sort(list.files(file.path(root, "f"), full.names = TRUE, pattern = "[.]R$")), 
                v5_sha, "")), lock$adapter_source_hash), "ADAPTER_SOURCE_CHANGED_AFTER_PREPARE")
            v5_assert(identical(as.character(unlist(lock$parent_source_sha256)), as.character(loaded$source_sha256)), 
                "FROZEN_ENGINE_CHANGED_AFTER_PREPARE")
            v5_check(file.path(root, "config/adapter_sources.json"), lock$adapter_sources_sha256)
            if (!is.null(lock$social_recovery_sha256)) v5_check(file.path(root, "config/social_recovery.json"), 
                lock$social_recovery_sha256)
            ctx <- readRDS(cache)
            v5_assert(identical(ctx$scope_hash, v5_sha(file.path(root, "input/APPROVED_165_ORIGINAL.csv"))), "PREPARED_APPROVED_SCOPE_CHANGED")
            ctx$engine <- loaded$engine
            ctx$discrete <- loaded$discrete
            ctx$run <- run
            ctx$root <- root
            ctx$input_lock_sha256 <- v5_sha(lp)
            ctx$task_cache <- new.env(parent = emptyenv())
            ctx$matrix_cache <- new.env(parent = emptyenv())
            ctx$monitor <- new.env(parent = emptyenv())
            ctx$monitor$process_start <- v5_utc()
            ctx$monitor$endpoint_status <- as.list(setNames(rep("NOT_ATTEMPTED", 5), ctx$cfg$outcomes))
            v5_install_hooks(ctx)
        }
        v5_task_dir <- function(ctx, p) {
            v5_assert(is.character(p) && length(p) == 1L && !grepl("(^|/)[.][.](/|$)", gsub(birth_private_setting("PRIVATE_INPUT"), 
                "/", p)), "UNSAFE_TASK_PATH")
            parent <- dirname(p)
            ancestor <- parent
            parts <- character()
            while (!dir.exists(ancestor) && !identical(dirname(ancestor), ancestor)) {
                parts <- c(basename(ancestor), parts)
                ancestor <- dirname(ancestor)
            }
            ancestor <- normalizePath(ancestor, mustWork = TRUE)
            run <- normalizePath(ctx$run, mustWork = TRUE)
            v5_assert(identical(ancestor, run) || startsWith(paste0(ancestor, "/"), paste0(run, "/")), "TASK_OUTSIDE_NEW_V5_RUN")
            if (length(parts)) parent <- do.call(file.path, as.list(c(ancestor, parts))) else parent <- ancestor
            if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
            rp <- normalizePath(parent, mustWork = TRUE)
            v5_assert(startsWith(paste0(rp, "/"), paste0(run, "/")), "TASK_OUTSIDE_NEW_V5_RUN")
            dest <- file.path(rp, basename(p))
            v5_assert(!file.exists(dest) && !dir.exists(dest) && dir.create(dest), "PRESERVE_EXISTING_ATTEMPT_NO_OVERWRITE")
            dest
        }
        v5_failed_status <- function(reason) {
            if (grepl("cannot allocate|cannot open|permission|no space|write error|unused argument|could not find|object .* not found|subscript|transport|socket|FILE_MISSING|SHA256|PATH|SOURCE_CHANGED", 
                reason, ignore.case = TRUE)) "FAILED_TECH" else "FAILED_STAT"
        }
        v5_base_row <- function(ctx, id, bid, e) {
            s <- ctx$scope[ctx$scope$indicator_id == id, , drop = FALSE]
            q <- data.frame(indicator_id = id, feature_id = id, label_zh = s$label_zh, type = s$type, domain = s$category, 
                endpoint_id = e, endpoint = e, endpoint_label = switch(e, `250` = "DM", `250.2` = "T2D", `411` = "IHD", 
                  `411.3` = "AG", `411.4` = "CAS"), bootstrap_id = bid, time_origin = if (s$category == "CMR") "imaging" else "baseline", 
                horizon_years = 5L, status = "HOLD", reason = "", model_hash = ctx$method_hash, method_hash = ctx$method_hash, 
                input_hash = ctx$input_lock_sha256, scope_hash = ctx$scope_hash, draw_hash = NA_character_, target_hash = NA_character_, 
                sample_hash = NA_character_, stringsAsFactors = FALSE)
            for (k in c("R00", "R10", "R11", "TE", "DE", "IE", "PM", "PM_signed", "PM_signed_percent", "N_target", 
                "N_complete", "N_known_5y", "events_5y", "death_5y", "censor_5y")) q[[k]] <- NA_real_
            q
        }
        v5_compute <- function(ctx, indicator_id, bootstrap_id, output_dir, attempt_id) {
            bid <- v5_id(bootstrap_id)
            v5_assert(indicator_id %in% ctx$scope$indicator_id, "INDICATOR_OUTSIDE_APPROVED_165")
            m <- ctx$monitor
            m$indicator_id <- indicator_id
            m$bootstrap_id <- bid
            m$attempt_id <- attempt_id
            m$selection_weight_hashes <- list()
            m$package_started_epoch <- as.numeric(Sys.time())
            m$package_started_monotonic <- proc.time()[["elapsed"]]
            m$endpoint_status <- as.list(setNames(rep("NOT_ATTEMPTED", 5), ctx$cfg$outcomes))
            v5_stage(ctx, "TASK_CLAIMED", "SHARED")
            v5_stage(ctx, "LOAD_CACHE", "SHARED")
            scope <- ctx$scope[ctx$scope$indicator_id == indicator_id, , drop = FALSE]
            ca <- ctx$features[[indicator_id]]
            var <- ctx$variables[[indicator_id]]
            rows <- parameters <- diagnostics <- weights <- list()
            receipt <- NULL
            if (scope$type != "continuous" && identical(var$status, "AVAILABLE_FROZEN")) {
                v5_stage(ctx, "DRAW_PREPARE", "SHARED")
                ans <- ctx$discrete$c40_compute(ctx, indicator_id, bid, output_dir, backend = "reference")
                raw <- read.csv(file.path(output_dir, "RESULTS.csv"), stringsAsFactors = FALSE)
                receipt <- jsonlite::fromJSON(file.path(output_dir, "DRAW.json"), simplifyVector = FALSE)
                for (e in ctx$cfg$outcomes) {
                  q <- v5_base_row(ctx, indicator_id, bid, e)
                  r <- raw[as.character(raw$endpoint_id) == e, , drop = FALSE]
                  for (k in names(r)) q[[k]] <- r[[k]]
                  q$status <- if (r$status == "VALID") "VALID" else if (grepl("BLOCKED|HOLD", r$status)) "HOLD" else v5_failed_status(r$reason)
                  wh <- m$selection_weight_hashes[[q$target_hash]]
                  if (!is.null(wh)) for (k in names(wh)) q[[k]] <- wh[[k]]
                  rows[[e]] <- q
                }
            } else {
                td <- v5_task_dir(ctx, output_dir)
                v5_stage(ctx, "DRAW_PREPARE", "SHARED")
                actual <- ctx$discrete$c40_draw(ctx, bid)
                receipt <- actual$receipt
                receipt$indicator_id <- indicator_id
                receipt$method_hash <- ctx$method_hash
                receipt$input_hash <- ctx$input_lock_sha256
                for (k in names(ctx$discrete$c40_persist_draw(ctx, actual))) receipt[[k]] <- ctx$discrete$c40_persist_draw(ctx, 
                  actual)[[k]]
                ix <- actual$ix
                cache <- new.env(parent = emptyenv())
                is_cmr <- scope$category == "CMR"
                dd_full <- if (is_cmr && ctx$imaging$status == "AVAILABLE_FROZEN") ctx$imaging$data else ctx$d_full
                eps <- if (is_cmr && ctx$imaging$status == "AVAILABLE_FROZEN") ctx$imaging$endpoints else ctx$endpoints
                warnings <- list()
                for (e in ctx$cfg$outcomes) {
                  m$endpoint_id <- e
                  q <- v5_base_row(ctx, indicator_id, bid, e)
                  q$draw_hash <- receipt$draw_hash
                  ep <- eps[[e]]
                  cellwarnings <- character()
                  if (var$status != "AVAILABLE_FROZEN") {
                    q$reason <- var$reason
                  } else if (ep$status != "AVAILABLE_FROZEN") {
                    q$reason <- if (is.null(ep$reason)) "FROZEN_ENDPOINT_UNAVAILABLE" else ep$reason
                  } else {
                    oldwarning <- getOption("b3.warning_sink")
                    options(b3.warning_sink = function(message) cellwarnings <<- c(cellwarnings, message))
                    got <- tryCatch(withCallingHandlers({
                      target <- ctx$eligibility[[indicator_id]] & ep$eligibility
                      pr <- ix[target[ix]]
                      v5_assert(length(pr) > 0L, "NO_ELIGIBLE_TARGET")
                      d <- dd_full[pr, , drop = FALSE]
                      r <- ep$risk[pr, , drop = FALSE]
                      base <- unique(c(ctx$cfg$base_covariates, ca$technical_covariates))
                      q$target_hash <- v5_hash(pr)
                      q$N_target <- length(pr)
                      key <- v5_hash(c(indicator_id, bid, q$target_hash, ctx$method_hash, ctx$input_lock_sha256, 
                        ca$column, base))
                      if (exists(key, envir = cache, inherits = FALSE)) {
                        nu <- get(key, envir = cache)
                      } else {
                        v5_stage(ctx, "SELECTION_WEIGHT", e)
                        mi <- ctx$engine$b3_missing(d, ca$column, base)
                        v5_stage(ctx, "MEDIATOR_FIT", "SHARED", list(target_hash = q$target_hash))
                        a <- ctx$engine$fp_a_path(d, ca$column, base, mi, ctx, paste0("single_", ca$candidate_id), 
                          bid)
                        physical <- ctx$engine$fp_physical(d, ca$column, mi, a$alpha, ctx)
                        v5_assert(all(physical$status == "VALID"), "SCENARIO_OUTSIDE_PHYSICAL_DOMAIN_NO_CLIPPING")
                        nu <- list(mi = mi, a = a, physical = physical)
                        assign(key, nu, envir = cache)
                      }
                      cc <- which(nu$mi$complete)
                      q$sample_hash <- v5_hash(pr[cc])
                      q$standardization_sample_hash <- q$sample_hash
                      q$N_complete <- length(cc)
                      q$prepared_unit <- ca$prepared_unit
                      q$prepared_transform <- ca$prepared_transform
                      q$inverse_transform <- ca$inverse_transform
                      q$endpoint_version <- ep$endpoint_version
                      q$risk_hash <- ep$risk_hash
                      q$target_id <- if (is_cmr) "INSTANCE2_SURVIVOR_ENDPOINT_DISEASE_FREE" else ca$target_key
                      q$shared_nuisance_key <- key
                      q$selection_weight_hash <- v5_hash(format(nu$mi$w, digits = 17, scientific = TRUE))
                      q$selection_sample_hash <- q$sample_hash
                      q$selection_target_hash <- q$target_hash
                      q$formula_hash <- v5_hash(c(base, ca$column, "FROZEN_SELECTION_T2D_IHD_AUX", "FROZEN_RIDGE_CENSOR_THETA1", 
                        "M3_AS_MEAN_WEIGHTED_INTERCEPT_ONLY", "SPLINE_5_35_65_95"))
                      fit <- ctx$engine$fp_m3_effect(d, r, ca$column, base, nu$mi, diagnostics = bid == 0L, reuse_alpha = nu$a$alpha)
                      res <- fit$result
                      v5_assert(nrow(res) == 1L && all(is.finite(unlist(res[1, c("R00", "R10", "R11", "TE", "DE", 
                        "IE")]))), "NONFINITE_M3_EFFECT")
                      res$numerical_status <- res$status
                      res$status <- NULL
                      for (k in names(res)) q[[k]] <- res[[k]]
                      q$status <- "VALID"
                      q$reason <- ""
                      parameters[[e]] <- fit$parameters
                      diagnostics[[e]] <- v5_bind(list(nu$physical, nu$a$rows, fit$support, fit$residual_diagnostics, 
                        fit$censor_PH))
                      weights[[e]] <- ctx$engine$fp_weight_summary(nu$mi$w, "SELECTION_IPW", e)
                      TRUE
                    }, warning = function(w) {
                      cellwarnings <<- c(cellwarnings, conditionMessage(w))
                      invokeRestart("muffleWarning")
                    }), error = function(err) {
                      q$status <<- v5_failed_status(conditionMessage(err))
                      q$reason <<- conditionMessage(err)
                      FALSE
                    })
                    options(b3.warning_sink = oldwarning)
                  }
                  ww <- unique(c(cellwarnings, if ("warnings" %in% names(q)) q$warnings else character()))
                  q$warning <- paste(ww, collapse = " | ")
                  warnings[[e]] <- data.frame(endpoint_id = rep(e, length(ww)), message = ww, stringsAsFactors = FALSE)
                  rows[[e]] <- q
                  m$endpoint_status[[e]] <- q$status
                  v5_stage(ctx, "WRITE_CELL", e)
                  v5_csv(q, file.path(td, paste0("CELL_", gsub("[.]", "_", e), ".csv")))
                }
                v5_stage(ctx, "WRITE_PACKAGE", "SHARED")
                v5_json(receipt, file.path(td, "DRAW.json"))
                saveRDS(parameters, file.path(td, "PARAMETERS.rds"), compress = FALSE)
                v5_csv(v5_bind(diagnostics), file.path(td, "DIAGNOSTICS.csv"))
                v5_csv(v5_bind(weights), file.path(td, "WEIGHTS.csv"))
                v5_csv(v5_bind(warnings), file.path(td, "WARNINGS.csv"))
            }
            for (e in names(rows)) {
                q <- rows[[e]]
                q$PM_signed <- q$PM
                q$PM_signed_percent <- 100 * q$PM
                q$censor_5y <- if ("early_censor" %in% names(q)) q$early_censor else NA_real_
                if (!"target_id" %in% names(q)) q$target_id <- if ("target_key" %in% names(q)) q$target_key else NA_character_
                for (k in c("TE", "DE", "IE")) q[[paste0(k, "_per_1000")]] <- 1000 * q[[k]]
                for (k in c("cluster_order_hash", "parent_frame_hash", "parent_cluster_membership_hash", "multiplicity_hash")) q[[k]] <- receipt[[k]]
                q$source_hash <- ctx$source_hash
                q$attempt_id <- attempt_id
                q$host <- v5_env("V5_HOST")
                q$source_label_retained <- scope$source_label_retained
                rows[[e]] <- q
                m$endpoint_status[[e]] <- q$status
                v5_stage(ctx, "WRITE_CELL", e)
                v5_csv(q, file.path(output_dir, paste0("V5_CELL_", gsub("[.]", "_", e), ".csv")))
            }
            v5_stage(ctx, "WRITE_PACKAGE", "SHARED")
            res <- v5_bind(rows)
            v5_assert(nrow(res) == 5L && setequal(as.character(res$endpoint_id), ctx$cfg$outcomes), "FIVE_INDEPENDENT_ENDPOINT_CELLS_REQUIRED")
            out <- normalizePath(output_dir, mustWork = TRUE)
            v5_csv(res, file.path(out, "V5_RESULTS.csv"))
            v5_csv(data.frame(total_seconds = proc.time()[["elapsed"]] - m$package_started_monotonic, pid = Sys.getpid(), 
                host = v5_env("V5_HOST")), file.path(out, "V5_TIMING.csv"))
            files <- list.files(out, full.names = TRUE)
            files <- files[!file.info(files)$isdir]
            seal <- list(complete = TRUE, indicator_id = indicator_id, bootstrap_id = bid, attempt_id = attempt_id, 
                method_hash = ctx$method_hash, input_hash = ctx$input_lock_sha256, draw_hash = receipt$draw_hash, 
                artifacts_sha256 = as.list(setNames(vapply(files, v5_sha, ""), basename(files))), cell_status = as.list(setNames(res$status, 
                  res$endpoint_id)), total_seconds = proc.time()[["elapsed"]] - m$package_started_monotonic)
            v5_json(seal, file.path(out, "V5_COMPLETE.json"))
            v5_stage(ctx, "TASK_FINISHED", "SHARED")
            list(ok = TRUE, indicator_id = indicator_id, bootstrap_id = bid, output_dir = out, rows = unname(lapply(seq_len(nrow(res)), 
                function(i) as.list(res[i, , drop = FALSE]))), complete_sha256 = v5_sha(file.path(out, "V5_COMPLETE.json")), 
                total_seconds = seal$total_seconds)
        }
    }), envir = context)
}

