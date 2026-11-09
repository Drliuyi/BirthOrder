# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

bi_cont_value <- function(x, default = NULL) if (is.null(x) || !length(x)) default else x

bi_cont_flag <- function(x) isTRUE(x) || identical(as.character(x), "TRUE")

bi_cont_save <- function(ctx, key, path, table) {
    directory <- if (key == "glycaemia") 
        "biochemistry"
    else key
    bi_save(ctx, file.path("raw", directory, path), table)
}

bi_cont_analyses <- function(feature, sensitivities) {
    selected <- Filter(function(s) is.null(s$features) || feature$feature_id %in% unlist(s$features), sensitivities)
    ids <- c("primary", vapply(selected, function(s) as.character(s$id), character(1)))
    if (bi_cont_flag(feature$bsa_index)) 
        ids <- c(ids, "bsa_indexed")
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) 
        stop("Sensitivity IDs must be nonempty and unique")
    ids
}

bi_cont_workers <- function(ctx, cfg, data, shared_gb = 0) {
    requested <- as.integer(bi_cont_value(cfg$max_workers, bi_cont_value(ctx$cfg$resources$max_workers, 1L)))
    if (length(requested) != 1L || !is.finite(requested) || requested < 1L) 
        stop("max_workers must be a positive integer")
    effective <- min(requested, 4L)
    reason <- "Forked exposure fits share one read-only feature matrix; master writes checkpoints"
    estimate <- max(0.25, as.numeric(object.size(data)) * 10/1024^3)
    reserve <- as.numeric(bi_cont_value(ctx$cfg$resources$minimum_free_gb, 2))
    ram_limit <- as.numeric(bi_cont_value(ctx$cfg$resources$ram_limit_gb, 24))
    if (!is.finite(ram_limit) || ram_limit <= 0 || !is.finite(shared_gb) || shared_gb < 0) 
        stop("Invalid continuous RAM budget")
    cap <- floor((ram_limit - shared_gb - reserve)/estimate)
    if (cap < 1L) 
        stop("CONTINUOUS_RAM_LIMIT_EXCEEDED_BY_SHARED_INPUTS")
    effective <- min(effective, cap)
    if (.Platform$OS.type != "unix") {
        effective <- 1L
        reason <- "Native Windows cannot fork; use reviewed WSL for parallel execution"
    }
    else if (file.exists("/proc/meminfo")) {
        line <- readLines("/proc/meminfo", warn = FALSE)
        available <- as.numeric(gsub("[^0-9]", "", line[grepl("^MemAvailable:", line)]))/1024^2
        if (length(available) != 1L || !is.finite(available)) 
            stop("Cannot validate available RAM for fork workers")
        cap <- floor((available - reserve)/estimate)
        if (cap < 1L) 
            stop("INSUFFICIENT_FREE_MEMORY_FOR_CONTINUOUS_WORKERS")
        effective <- min(effective, cap)
    }
    list(requested = requested, effective = effective, estimated_worker_gb = estimate, shared_input_gb = shared_gb, 
        ram_limit_gb = ram_limit, reason = reason)
}

bi_cont_map <- function(x, fun, workers) {
    if (workers <= 1L) 
        return(lapply(x, fun))
    ans <- parallel::mclapply(x, fun, mc.cores = min(workers, length(x)), mc.preschedule = FALSE, mc.set.seed = FALSE, 
        mc.cleanup = TRUE)
    if (any(vapply(ans, inherits, logical(1), "try-error"))) 
        stop("Continuous worker error: ", paste(unlist(ans[vapply(ans, inherits, logical(1), "try-error")]), collapse = "; "))
    ans
}

bi_cont_mapping_columns <- function(x, source) {
    if (!is.list(x)) 
        return(character())
    here <- if (identical(x$source, source) && !is.null(x$column)) 
        as.character(x$column)
    else character()
    unique(c(here, unlist(lapply(x, bi_cont_mapping_columns, source = source), use.names = FALSE)))
}

bi_cont_source_dialect <- function(path, ctx = NULL) {
    answer <- list(dialect = "standard", id_column = "eid")
    entries <- if (is.null(ctx)) 
        list()
    else bi_cont_value(ctx$cfg$source_dialects, list())
    if (!is.list(entries)) 
        stop("source_dialects must be a list of input/dialect/id_column records")
    if (!length(entries)) 
        return(answer)
    paths <- vapply(entries, function(entry) {
        if (!is.list(entry) || length(entry$input) != 1L || is.na(entry$input) || !nzchar(entry$input)) 
            stop("source_dialects entry requires one resolved input path")
        if (length(entry$dialect) != 1L || is.na(entry$dialect) || !entry$dialect %in% c("standard", "backslash", 
            "literal")) 
            stop("Invalid explicit source dialect")
        if (length(entry$id_column) != 1L || is.na(entry$id_column) || !nzchar(entry$id_column)) 
            stop("source_dialects entry requires one ID column")
        normalizePath(entry$input, winslash = "/", mustWork = FALSE)
    }, character(1))
    if (anyDuplicated(paths)) 
        stop("Duplicate resolved inputs in source_dialects")
    selected <- which(paths == normalizePath(path, winslash = "/", mustWork = TRUE))
    if (length(selected)) 
        answer <- entries[[selected]][c("dialect", "id_column")]
    answer
}

bi_cont_read <- function(path, columns, ctx = NULL) {
    if (!file.exists(path)) 
        stop("Continuous source missing: ", path)
    columns <- unique(as.character(columns[!is.na(columns) & nzchar(columns)]))
    if (grepl("[.]rds$", path, ignore.case = TRUE)) {
        d <- bi_read_table(path)
        return(d[, intersect(columns, names(d)), drop = FALSE])
    }
    dialect <- bi_cont_source_dialect(path, ctx)
    if (grepl("[.]gz$", path, ignore.case = TRUE) || dialect$dialect != "standard") {
        root <- if (!is.null(ctx) && !is.null(ctx$root)) 
            ctx$root
        else Sys.getenv("BI_ROOT")
        helper <- file.path(root, "f", "streaming_tsv_select.py")
        if (length(root) != 1L || !nzchar(root) || !file.exists(helper)) 
            stop("Streaming gzip projection helper unavailable; set BI_ROOT or ctx$root")
        python <- Sys.getenv("BI_PYTHON", unset = "python3")
        if (!nzchar(Sys.which(python)) && !file.exists(python)) 
            stop("Python 3 required for streaming gzip projection")
        delimiter <- if (grepl("[.]csv([.]gz)?$", path, ignore.case = TRUE)) 
            "comma"
        else "tab"
        args <- c(python, helper, "--input", normalizePath(path, mustWork = TRUE), "--columns-json", jsonlite::toJSON(columns, 
            auto_unbox = FALSE), "--delimiter", delimiter, "--dialect", dialect$dialect, "--id-column", dialect$id_column, 
            "--allow-missing", "--progress-every", "50000")
        command <- paste(vapply(args, shQuote, character(1)), collapse = " ")
        header_error <- tempfile("projection_header_diagnostic_")
        body_error <- tempfile("projection_body_diagnostic_")
        completion <- tempfile("projection_completion_", fileext = ".json")
        on.exit(unlink(c(header_error, body_error, completion)), add = TRUE)
        diagnostic <- function(p) if (file.exists(p)) 
            paste(readLines(p, warn = FALSE), collapse = "\n")
        else ""
        header_text <- suppressWarnings(system(paste(command, shQuote("--header-only"), "2>", shQuote(header_error)), 
            intern = TRUE))
        header_status <- attr(header_text, "status")
        if (!is.null(header_status) && header_status != 0L) 
            stop("SOURCE_PROJECTION_FAILED: ", diagnostic(header_error), call. = FALSE)
        header <- names(utils::read.csv(text = paste(header_text, collapse = "\n"), nrows = 0L, check.names = FALSE))
        body_command <- paste(command, shQuote("--audit-json"), shQuote(completion), shQuote("--diagnostic-log"), 
            shQuote(body_error))
        consumer_warnings <- character()
        projected <- tryCatch(withCallingHandlers(as.data.frame(data.table::fread(cmd = body_command, sep = ",", 
            header = TRUE, col.names = header, colClasses = if (dialect$id_column %in% header) 
                setNames("character", dialect$id_column)
            else NULL, check.names = FALSE, data.table = FALSE, showProgress = FALSE)), warning = function(w) {
            consumer_warnings <<- c(consumer_warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }), error = identity)
        details <- diagnostic(body_error)
        if (!file.exists(completion) || grepl("[continuous_projection] ERROR", details, fixed = TRUE)) {
            if (!nzchar(details)) 
                details <- "Parser completion audit absent; source stream did not complete"
            stop("SOURCE_PROJECTION_FAILED: ", details, call. = FALSE)
        }
        audit <- tryCatch(jsonlite::read_json(completion, simplifyVector = TRUE), error = function(e) NULL)
        if (is.null(audit) || !identical(audit$status, "PASS") || !isTRUE(audit$full_scan) || isTRUE(audit$header_only) || 
            isTRUE(audit$validation_only) || !is.null(audit$max_records) || !identical(audit$dialect, dialect$dialect) || 
            !identical(audit$id_column, dialect$id_column) || length(audit$rows) != 1L || !is.numeric(audit$rows) || 
            !is.finite(audit$rows) || audit$rows < 0 || length(audit$unique_ids) != 1L || audit$rows != audit$unique_ids) {
            stop("SOURCE_PROJECTION_FAILED: invalid completion audit; action=verify_helper_version_and_rerun_read_validation", 
                call. = FALSE)
        }
        if (inherits(projected, "error") || length(consumer_warnings)) {
            stop("PROJECTED_CSV_CONSUMER_FAILED after successful source validation; error_class=", if (inherits(projected, 
                "error")) 
                class(projected)[1L]
            else "none", " warning_count=", length(consumer_warnings), "; no participant output accepted; action=verify_R_data_table_adapter", 
                call. = FALSE)
        }
        if (nrow(projected) != audit$rows || ncol(projected) != audit$selected_columns || !identical(names(projected), 
            header)) 
            stop("PROJECTED_CSV_SHAPE_MISMATCH", call. = FALSE)
        quote_probe <- data.table::fread(text = "x\n\"a\"\"\"\"b\"\n", sep = ",", header = TRUE, colClasses = "character", 
            data.table = FALSE)[[1L]]
        if (identical(quote_probe, "a\"\"\"\"b")) {
            for (nm in names(projected)[vapply(projected, is.character, logical(1))]) {
                projected[[nm]] <- gsub("\"\"", "\"", projected[[nm]], fixed = TRUE)
            }
        }
        else if (!identical(quote_probe, "a\"\"b")) 
            stop("PROJECTED_CSV_QUOTE_DECODING_UNSUPPORTED", call. = FALSE)
        if (dialect$id_column %in% names(projected)) {
            ids <- projected[[dialect$id_column]]
            if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) 
                stop("PROJECTED_CSV_INVALID_ID_KEYS", call. = FALSE)
        }
        return(projected)
    }
    header <- names(data.table::fread(path, nrows = 0L, check.names = FALSE, data.table = FALSE, showProgress = FALSE))
    selected <- intersect(columns, header)
    if (!length(selected)) 
        stop("No configured columns present in continuous source")
    as.data.frame(data.table::fread(path, select = selected, check.names = FALSE, data.table = FALSE, showProgress = FALSE))
}

bi_cont_numbers <- function(x) {
    if (is.null(x) || !length(x) || all(is.na(x)) || identical(x, "")) 
        return(numeric())
    z <- suppressWarnings(as.numeric(unlist(strsplit(paste(x, collapse = ";"), "[;,]"))))
    if (anyNA(z)) 
        stop("Non-numeric configured missing code or bound")
    z
}

bi_cont_dictionary <- function(x) {
    if (is.character(x) && length(x) == 1L) 
        x <- bi_read_table(x)
    if (is.data.frame(x)) {
        rows <- lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
    }
    else if (is.list(x)) 
        rows <- x
    else stop("dictionary must be a table path or feature records")
    if (!length(rows)) 
        stop("No feature dictionary records")
    required <- c("feature_id", "column", "label", "unit", "family", "kind", "include", "qc_verified", "transform")
    defaults <- list(primary_id = "", missing_codes = "", valid_min = NA_real_, valid_max = NA_real_, winsor_lower = NA_real_, 
        winsor_upper = NA_real_, qc_column = "", qc_allowed = "", qc_source = "features", numeric_parse = "require_numeric", 
        bsa_index = FALSE, visit = "")
    out <- lapply(rows, function(r) {
        absent <- setdiff(required, names(r))
        if (length(absent)) 
            stop("Dictionary missing fields: ", paste(absent, collapse = ","))
        for (nm in names(defaults)) if (is.null(r[[nm]]) || !length(r[[nm]])) 
            r[[nm]] <- defaults[[nm]]
        for (nm in c("primary_id", "missing_codes", "qc_column", "qc_allowed", "qc_source", "visit")) {
            if (length(r[[nm]]) != 1L && nm != "missing_codes") 
                stop("Optional dictionary field must be scalar: ", nm)
            if (all(is.na(r[[nm]]))) 
                r[[nm]] <- ""
        }
        for (nm in c("feature_id", "column", "label", "unit", "family", "kind", "transform")) {
            if (length(r[[nm]]) != 1L || is.na(r[[nm]]) || !nzchar(as.character(r[[nm]]))) 
                stop("Empty dictionary field: ", nm)
            r[[nm]] <- as.character(r[[nm]])
        }
        r$include <- bi_cont_flag(r$include)
        r$qc_verified <- bi_cont_flag(r$qc_verified)
        r$bsa_index <- bi_cont_flag(r$bsa_index)
        if (!r$transform %in% c("identity", "log", "log1p", "rank_normal")) 
            stop("Unsupported locked transformation: ", r$transform)
        if (!r$numeric_parse %in% c("require_numeric", "non_numeric_to_missing")) 
            stop("Unsupported explicit numeric parsing rule")
        if (!r$kind %in% c("primary", "derived")) 
            stop("kind must be primary or derived")
        if (!nzchar(as.character(r$primary_id))) 
            r$primary_id <- r$feature_id
        r$missing_codes <- paste(bi_cont_numbers(r$missing_codes), collapse = ";")
        for (nm in c("valid_min", "valid_max", "winsor_lower", "winsor_upper")) {
            original <- r[[nm]]
            r[[nm]] <- suppressWarnings(as.numeric(r[[nm]]))
            if (length(r[[nm]]) != 1L) 
                stop("Bound must be scalar: ", nm)
            if (is.na(r[[nm]]) && !is.na(original) && nzchar(as.character(original))) 
                stop("Invalid numeric dictionary bound: ", nm)
        }
        if (is.finite(r$valid_min) && is.finite(r$valid_max) && r$valid_min >= r$valid_max) 
            stop("Invalid feature QC range")
        w <- c(r$winsor_lower, r$winsor_upper)
        if (any(is.finite(w)) && (!all(is.finite(w)) || w[1L] < 0 || w[2L] > 1 || w[1L] >= w[2L])) 
            stop("Invalid winsor quantiles")
        r[unique(c(required, names(defaults)))]
    })
    ans <- bi_bind(lapply(out, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
    for (nm in c("feature_id", "column", "primary_id")) if (anyDuplicated(ans[[nm]][ans$include])) 
        stop("Duplicate included ", nm)
    ans
}

bi_cont_gate <- function(cfg) {
    if (is.null(cfg) || !bi_cont_flag(cfg$enabled)) 
        return("Module not enabled in locked configuration")
    if (!bi_cont_flag(cfg$verified)) 
        return("Source/dictionary linkage has not been verified")
    if (is.null(cfg$source_version) || !nzchar(cfg$source_version)) 
        return("Verified source version is required")
    if (!bi_cont_flag(cfg$metadata$verified) || is.null(cfg$metadata$visit) || !nzchar(cfg$metadata$visit)) 
        return("Measurement visit metadata has not been verified")
    if (!bi_cont_flag(cfg$qc$verified) || is.null(cfg$qc$description) || !nzchar(cfg$qc$description)) 
        return("Assay/image QC and missing-code semantics have not been verified")
    if (!identical(cfg$qc$input_imputation, "none")) 
        return("Unimputed input provenance has not been verified; preimputed values cannot be silently treated as observed measurements")
    if (is.null(cfg$dictionary)) 
        return("Feature dictionary is not connected")
    for (nm in c("age", "center", "year")) if (is.null(cfg$metadata[[nm]]$column)) 
        return(paste("Visit-specific mapping missing:", nm))
    NULL
}

bi_cont_lookup <- function(mapping, cohort, features, index, label) {
    metadata <- attr(features, "bi_metadata_tables")
    if (is.null(mapping$source) || !mapping$source %in% c("cohort", "features", names(metadata)) || is.null(mapping$column)) 
        stop("Invalid explicit mapping: ", label)
    src <- if (mapping$source == "cohort") 
        cohort
    else if (mapping$source == "features") 
        features
    else metadata[[mapping$source]]$data
    if (!mapping$column %in% names(src)) 
        stop("Configured column absent for ", label, ": ", mapping$column)
    if (mapping$source == "cohort") 
        src[[mapping$column]]
    else if (mapping$source == "features") 
        src[[mapping$column]][index]
    else src[[mapping$column]][metadata[[mapping$source]]$index]
}

bi_cont_rule <- function(rule, cohort, features, index) {
    z <- bi_cont_lookup(rule, cohort, features, index, "QC/filter")
    if (is.null(rule$allowed) || !length(rule$allowed)) 
        stop("QC/filter allowed values must be explicit")
    !is.na(z) & as.character(z) %in% as.character(unlist(rule$allowed))
}

bi_cont_technical <- function(mapping, cohort, features, index, label) {
    z <- bi_cont_lookup(mapping, cohort, features, index, label)
    if (identical(mapping$type, "factor")) 
        return(factor(z))
    if (identical(mapping$type, "date_month")) {
        if (inherits(z, "POSIXt")) 
            z <- as.Date(z, tz = "UTC")
        if (!inherits(z, "Date")) {
            if (!is.character(z)) 
                stop("Assay dates must be Date or ISO date strings: ", mapping$column)
            parsed <- suppressWarnings(as.Date(z, format = "%Y-%m-%d"))
            if (any(!is.na(z) & nzchar(z) & is.na(parsed))) 
                stop("Unparseable configured assay dates: ", mapping$column)
            z <- parsed
        }
        return(factor(format(z, "%Y-%m")))
    }
    if (!identical(mapping$type, "numeric") || !is.numeric(z)) 
        stop("Technical mapping must specify factor, date_month, or an actual numeric column: ", mapping$column)
    z[!is.finite(z)] <- NA_real_
    if (!is.null(mapping$missing_codes)) 
        z[z %in% bi_cont_numbers(mapping$missing_codes)] <- NA_real_
    if (!is.null(mapping$valid_min)) 
        z[z < as.numeric(mapping$valid_min)] <- NA_real_
    if (!is.null(mapping$valid_max)) 
        z[z > as.numeric(mapping$valid_max)] <- NA_real_
    z
}

bi_cont_load <- function(ctx, cfg, key) {
    cohort <- bi_cohort(ctx)
    input <- bi_cont_value(cfg$input, "cohort")
    dict <- bi_cont_dictionary(cfg$dictionary)
    id <- bi_cont_value(cfg$id_column, "eid")
    columns <- unique(c(if (id != "rownames") id, dict$column[dict$include], dict$qc_column[dict$include & dict$qc_source == 
        "features"], bi_cont_mapping_columns(cfg, "features"), unlist(cfg$select_columns)))
    features <- if (identical(input, "cohort")) 
        cohort
    else bi_cont_read(input, columns, ctx)
    if (is.matrix(features)) 
        features <- as.data.frame(features, stringsAsFactors = FALSE)
    if (!is.data.frame(features)) 
        stop("Continuous input must be an explicitly keyed wide table")
    ids <- if (identical(id, "rownames")) 
        rownames(features)
    else {
        if (!id %in% names(features)) 
            stop("Feature ID column missing: ", id)
        as.character(features[[id]])
    }
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) 
        stop("Missing or duplicate feature participant keys")
    if (!"eid" %in% names(cohort) || anyNA(cohort$eid) || anyDuplicated(cohort$eid)) 
        stop("Invalid cohort participant keys")
    index <- match(as.character(cohort$eid), ids)
    metadata <- list()
    tables <- bi_cont_value(cfg$metadata_tables, list())
    if (length(tables) && (is.null(names(tables)) || any(!nzchar(names(tables))) || anyDuplicated(names(tables)) || 
        any(names(tables) %in% c("cohort", "features")))) 
        stop("Metadata table source names must be explicit, unique and nonreserved")
    for (name in names(tables)) {
        spec <- tables[[name]]
        if (!bi_cont_flag(spec$verified) || is.null(spec$input) || is.null(spec$id_column)) 
            stop("Unverified metadata source: ", name)
        selected <- unique(c(spec$id_column, unlist(spec$select_columns), bi_cont_mapping_columns(cfg, name), dict$qc_column[dict$include & 
            dict$qc_source == name]))
        md <- bi_cont_read(spec$input, selected, ctx)
        if (!spec$id_column %in% names(md)) 
            stop("Metadata participant key absent: ", name)
        mid <- as.character(md[[spec$id_column]])
        if (anyNA(mid) || any(!nzchar(mid)) || anyDuplicated(mid)) 
            stop("Duplicate/missing metadata participant keys: ", name)
        metadata[[name]] <- list(data = md, index = match(as.character(cohort$eid), mid))
    }
    attr(features, "bi_metadata_tables") <- metadata
    model_fields <- c("eid", "border12", "border1n", "border_rank", "border_relative", "S", "age", "sex", "ethnicity", 
        "center", "assessment_year")
    if (length(setdiff(model_fields, names(cohort)))) 
        stop("Frozen continuous model fields missing")
    d <- cohort[, model_fields, drop = FALSE]
    d$age <- bi_cont_lookup(cfg$metadata$age, cohort, features, index, "measurement age")
    if (!is.numeric(d$age)) 
        stop("Verified measurement age must be numeric years")
    d$age[!is.finite(d$age)] <- NA_real_
    d$center <- factor(bi_cont_lookup(cfg$metadata$center, cohort, features, index, "measurement center"))
    d$assessment_year <- factor(bi_cont_lookup(cfg$metadata$year, cohort, features, index, "measurement year"))
    extra <- character()
    mappings <- c(bi_cont_value(cfg$metadata$batches, list()), bi_cont_value(cfg$metadata$extra_covariates, list()))
    for (i in seq_along(mappings)) {
        m <- mappings[[i]]
        nm <- paste0("bi_technical_", i)
        d[[nm]] <- bi_cont_technical(m, cohort, features, index, nm)
        extra <- c(extra, nm)
    }
    keep <- !is.na(index)
    for (r in bi_cont_value(cfg$qc$sample_rules, list())) keep <- keep & bi_cont_rule(r, cohort, features, index)
    dict$exclusion_reason <- ifelse(!dict$include, "EXCLUDED_PREDECLARED", ifelse(!dict$qc_verified, "QC_UNVERIFIED", 
        ifelse(!dict$column %in% names(features), "COLUMN_UNAVAILABLE", "")))
    bad_visit <- nzchar(dict$visit) & dict$visit != cfg$metadata$visit
    dict$exclusion_reason[bad_visit] <- "VISIT_MISMATCH"
    if (key %in% c("nmr", "olink")) 
        dict$family <- paste0(key, "_all")
    list(cohort = cohort, features = features, index = index, data = d, extra = extra, keep = keep, dictionary = dict, 
        flow = data.frame(module = key, cohort_N = nrow(cohort), source_N = length(ids), matched_N = sum(!is.na(index)), 
            unmatched_source_N = sum(!ids %in% as.character(cohort$eid)), sample_qc_pass_N = sum(keep), measurement_visit = cfg$metadata$visit, 
            source_version = cfg$source_version, stringsAsFactors = FALSE))
}

bi_cont_prepare <- function(v, feature, keep, feature_qc = NULL) {
    invalid_text <- 0L
    if (!is.numeric(v)) {
        if (!identical(feature$numeric_parse, "non_numeric_to_missing") || !is.character(v)) 
            stop("Configured continuous feature is not numeric: ", feature$column)
        parsed <- suppressWarnings(as.numeric(v))
        invalid_text <- sum(!is.na(v) & nzchar(v) & is.na(parsed))
        v <- parsed
    }
    n_input <- sum(!is.na(v))
    v[!is.finite(v)] <- NA_real_
    v[v %in% bi_cont_numbers(feature$missing_codes)] <- NA_real_
    if (is.finite(feature$valid_min)) 
        v[v < feature$valid_min] <- NA_real_
    if (is.finite(feature$valid_max)) 
        v[v > feature$valid_max] <- NA_real_
    if (!is.null(feature_qc)) 
        keep <- keep & feature_qc
    v[!keep] <- NA_real_
    if (all(is.finite(c(feature$winsor_lower, feature$winsor_upper))) && any(!is.na(v))) {
        q <- stats::quantile(v, probs = c(feature$winsor_lower, feature$winsor_upper), na.rm = TRUE, names = FALSE, 
            type = 7)
        v <- pmax(q[1L], pmin(q[2L], v))
    }
    raw <- v
    invalid <- rep(FALSE, length(v))
    if (feature$transform == "log") 
        invalid <- !is.na(v) & v <= 0
    if (feature$transform == "log1p") 
        invalid <- !is.na(v) & v <= -1
    v[invalid] <- NA_real_
    raw[invalid] <- NA_real_
    if (feature$transform == "log") 
        v <- log(v)
    if (feature$transform == "log1p") 
        v <- log1p(v)
    if (feature$transform == "rank_normal") {
        ok <- !is.na(v)
        v[ok] <- stats::qnorm((rank(v[ok], ties.method = "average") - 0.5)/sum(ok))
    }
    list(y = v, raw = raw, qc = data.frame(feature_id = feature$feature_id, input_nonmissing_N = n_input, retained_measurement_N = sum(!is.na(v)), 
        transformation_invalid_N = sum(invalid), invalid_numeric_text_N = invalid_text, transformation = feature$transform, 
        original_unit = feature$unit, stringsAsFactors = FALSE))
}

bi_cont_good <- function(r) nrow(r) == 1L && r$status %in% c("PASS", "OK") && is.finite(r$beta) && is.finite(r$se) && 
    is.finite(r$p)

bi_cont_reference <- function(ctx, d, exposure, extra, result) {
    if (!bi_cont_good(result)) {
        expected <- grepl("INSUFFICIENT_CONTINUOUS_N|CONSTANT_EXPOSURE", as.character(result$error))
        return(list(status = if (expected) "NOT_ESTIMABLE" else "FAIL", reason = paste(result$status, result$error), 
            delta = NA_real_))
    }
    form <- bi_formula(exposure, extra = extra)
    s <- bi_sample(d, form, exposure)
    reference <- stats::lm(s$formula, data = s$data, na.action = stats::na.fail)
    co <- summary(reference)$coefficients
    if (!exposure %in% rownames(co)) 
        return(list(status = "FAIL", reason = "Reference coefficient absent", delta = NA_real_))
    expected <- c(co[exposure, 1L], co[exposure, 2L], co[exposure, 4L])
    observed <- c(result$beta, result$se, result$p)
    delta <- max(abs(expected - observed)/pmax(1, abs(expected)))
    sd_y <- stats::sd(s$data$y)
    z <- s$data
    z$y <- z$y/sd_y
    scaled <- summary(stats::lm(s$formula, data = z, na.action = stats::na.fail))$coefficients[exposure, 1L]
    correct <- is.finite(delta) && delta < 1e-07 && nrow(s$data) == result$N && is.finite(sd_y) && sd_y > 0 && 
        abs(scaled - result$beta/sd_y) < 1e-07
    list(status = if (correct) "PASS" else "FAIL", reason = if (correct) "Gaussian/lm and standardized-beta equality" else "Independent numeric/sample mismatch", 
        delta = delta)
}

bi_cont_empty_result <- function(exposure, status, reason) {
    data.frame(exposure = exposure, N = NA_integer_, cases = NA_integer_, noncases = NA_integer_, beta = NA_real_, 
        se = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_, p_raw = NA_real_, status = status, 
        warning = "", error = reason, formula = "", rank = NA_integer_, columns = NA_integer_, sample_hash = "", 
        stringsAsFactors = FALSE)
}

bi_cont_fit_one <- function(ctx, d, exposure, extra, feature, raw, key, analysis, check) {
    result <- bi_fit(ctx, d, exposure, extra = extra, family = "gaussian")
    result$trait <- feature$feature_id
    result$feature_id <- feature$feature_id
    result$label <- feature$label
    result$primary_id <- feature$primary_id
    result$analysis <- analysis
    result$kind <- feature$kind
    result$original_unit <- feature$unit
    result$transform <- feature$transform
    result$effect_unit <- if (feature$transform == "identity") 
        feature$unit
    else paste0(feature$transform, "(", feature$unit, ")")
    result$family_id <- paste(feature$family, analysis, sep = ":")
    result$test_id <- paste(key, feature$feature_id, analysis, exposure, sep = ":")
    result$sd_outcome <- NA_real_
    result$standardized_beta <- result$standardized_CI_low <- result$standardized_CI_high <- NA_real_
    result$original_beta <- result$original_CI_low <- result$original_CI_high <- result$original_p <- NA_real_
    result$original_status <- as.character(result$status)
    if (bi_cont_good(result)) {
        sample <- bi_sample(d, bi_formula(exposure, extra = extra), exposure)$data
        result$sd_outcome <- stats::sd(sample$y)
        if (is.finite(result$sd_outcome) && result$sd_outcome > 0) {
            result$standardized_beta <- result$beta/result$sd_outcome
            result$standardized_CI_low <- result$CI_low/result$sd_outcome
            result$standardized_CI_high <- result$CI_high/result$sd_outcome
        }
        raw_result <- result
        if (feature$transform != "identity") {
            raw_d <- d
            raw_d$y <- raw
            raw_result <- bi_fit(ctx, raw_d, exposure, extra = extra, family = "gaussian")
            if (bi_cont_good(raw_result) && (!identical(as.character(raw_result$sample_hash), as.character(result$sample_hash)) || 
                raw_result$N != result$N)) 
                stop("Original/transformed fit sample mismatch")
        }
        result$original_status <- as.character(raw_result$status)
        if (bi_cont_good(raw_result)) {
            result$original_beta <- raw_result$beta
            result$original_CI_low <- raw_result$CI_low
            result$original_CI_high <- raw_result$CI_high
            result$original_p <- raw_result$p
        }
    }
    reference <- if (check) 
        bi_cont_reference(ctx, d, exposure, extra, result)
    else list(status = "NOT_REQUESTED", reason = "", delta = NA_real_)
    list(result = result, selfcheck = data.frame(feature_id = feature$feature_id, exposure = exposure, status = reference$status, 
        reason = reference$reason, relative_max_difference = reference$delta, stringsAsFactors = FALSE))
}

bi_cont_sensitivity <- function(ctx, spec, dat, cfg) {
    if (!bi_cont_flag(spec$verified)) 
        return(list(available = FALSE, reason = "Sensitivity field/definition not verified"))
    if (identical(spec$type, "outcome_absent")) {
        y <- bi_outcomes(ctx)
        if (is.null(spec$code) || !spec$code %in% names(y)) 
            return(list(available = FALSE, reason = "Outcome status unavailable"))
        if (anyDuplicated(y$eid)) 
            stop("Duplicate outcome participant keys")
        value <- y[[spec$code]][match(dat$cohort$eid, y$eid)]
        if (any(!is.na(value) & !value %in% c(0, 1))) 
            stop("Outcome must preserve verified 0/1/NA")
        keep <- !is.na(value) & value == 0
    }
    else if (identical(spec$type, "keep_values")) {
        keep <- bi_cont_rule(spec, dat$cohort, dat$features, dat$index)
    }
    else stop("Unsupported predeclared sensitivity type")
    list(available = TRUE, keep = keep, reason = "")
}

bi_cont_selection <- function(ctx, key, cfg, dat) {
    ans <- list()
    if (!is.null(cfg$metadata$baseline_date) && !is.null(cfg$metadata$measurement_date)) {
        baseline <- bi_cont_lookup(cfg$metadata$baseline_date, dat$cohort, dat$features, dat$index, "baseline date")
        measurement <- bi_cont_lookup(cfg$metadata$measurement_date, dat$cohort, dat$features, dat$index, "measurement date")
        if (!inherits(baseline, "Date") || !inherits(measurement, "Date")) 
            stop("Verified selection dates must use Date class")
        delay <- as.numeric(measurement - baseline)/365.25
        good <- dat$keep & is.finite(delay) & delay >= 0
        ans[[length(ans) + 1L]] <- data.frame(metric = "baseline_to_measurement_years", N = sum(good), mean = if (any(good)) 
            mean(delay[good])
        else NA_real_, sd = if (sum(good) > 1) 
            stats::sd(delay[good])
        else NA_real_, invalid_negative_N = sum(dat$keep & !is.na(delay) & delay < 0), status = "PASS", stringsAsFactors = FALSE)
    }
    else if (key == "imaging") 
        ans[[length(ans) + 1L]] <- data.frame(metric = "baseline_to_measurement_years", N = NA_integer_, mean = NA_real_, 
            sd = NA_real_, invalid_negative_N = NA_integer_, status = "UNAVAILABLE_DATES")
    for (nm in names(bi_cont_value(cfg$metadata$prior_disease, list()))) {
        m <- cfg$metadata$prior_disease[[nm]]
        v <- bi_cont_lookup(m, dat$cohort, dat$features, dat$index, paste("prior disease", nm))
        if (any(!is.na(v) & !v %in% c(0, 1))) 
            stop("Prior disease status must use verified 0/1/NA")
        ans[[length(ans) + 1L]] <- data.frame(metric = paste0("disease_present_at_measurement_", nm), N = sum(dat$keep & 
            !is.na(v)), cases = sum(dat$keep & !is.na(v) & v == 1), status = "PASS")
    }
    if (length(ans)) 
        bi_cont_save(ctx, key, "selection_qc.csv", bi_bind(ans))
}

bi_feature_stage <- function(ctx, key) {
    cfg <- bi_module_config(ctx, key)
    gate <- bi_cont_gate(cfg)
    if (!is.null(gate)) 
        return(bi_skip(ctx, gate))
    exposures <- c("border12", "border1n", "border_rank", "border_relative")
    loaded <- tryCatch(bi_cont_load(ctx, cfg, key), error = function(e) e)
    if (inherits(loaded, "error")) 
        return(list(status = "FAILED", reason = paste("Adapter validation:", conditionMessage(loaded))))
    dat <- loaded
    dict <- dat$dictionary
    bi_cont_save(ctx, key, "feature_dictionary.csv", dict)
    bi_cont_save(ctx, key, "cohort_flow.csv", dat$flow)
    bi_cont_selection(ctx, key, cfg, dat)
    chosen <- which(dict$exclusion_reason == "")
    if (!length(chosen)) 
        return(bi_skip(ctx, "No predeclared, QC-verified, available features"))
    chunk_size <- as.integer(bi_cont_value(cfg$chunk_size, 25L))
    if (length(chunk_size) != 1L || !is.finite(chunk_size) || chunk_size < 1L || chunk_size > 1000L) 
        stop("chunk_size must be 1..1000")
    limit <- as.numeric(bi_cont_value(cfg$max_elapsed_seconds, 86400))
    if (length(limit) != 1L || !is.finite(limit) || limit <= 0) 
        stop("max_elapsed_seconds must be positive and finite")
    shared_gb <- (as.numeric(object.size(dat$cohort)) + as.numeric(object.size(dat$features)) + as.numeric(object.size(dat$data)))/1024^3
    workers <- bi_cont_workers(ctx, cfg, dat$data, shared_gb)
    bi_log(ctx, "continuous_resources", requested_workers = workers$requested, effective_workers = workers$effective, 
        estimated_worker_gb = workers$estimated_worker_gb, shared_input_gb = workers$shared_input_gb, ram_limit_gb = workers$ram_limit_gb, 
        reason = workers$reason, chunk_size = chunk_size)
    started <- proc.time()[["elapsed"]]
    results <- qc <- checks <- sensitivity_audit <- list()
    planned <- dict[chosen, , drop = FALSE]
    first_n <- min(10L, nrow(planned))
    stopped <- ""
    sensitivities <- bi_cont_value(cfg$sensitivities, list())
    registry <- list()
    for (i in seq_len(nrow(planned))) {
        f <- as.list(planned[i, , drop = FALSE])
        for (a in bi_cont_analyses(f, sensitivities)) for (e in exposures) {
            registry[[length(registry) + 1L]] <- data.frame(feature_id = f$feature_id, analysis = a, exposure = e, 
                family_id = paste(f$family, a, sep = ":"), test_id = paste(key, f$feature_id, a, e, sep = ":"))
        }
    }
    bi_cont_save(ctx, key, "planned_tests.csv", bi_bind(registry))
    for (i in seq_len(nrow(planned))) {
        feature <- as.list(planned[i, , drop = FALSE])
        if (i > first_n && !nzchar(stopped) && i == first_n + 1L) {
            ctab <- bi_bind(checks)
            if (any(ctab$status == "FAIL") || !any(ctab$status == "PASS")) 
                stopped <- "NUMERIC_SELFCHECK_FAILED"
        }
        if (!nzchar(stopped) && proc.time()[["elapsed"]] - started > limit) 
            stop("MODULE_TIME_LIMIT: completed feature caches preserved; resume with the same configuration")
        analyses <- bi_cont_analyses(feature, sensitivities)
        if (nzchar(stopped)) {
            for (analysis in analyses) for (exposure in exposures) {
                r <- bi_cont_empty_result(exposure, stopped, "Planned row retained; fitting not attempted")
                r$trait <- r$feature_id <- feature$feature_id
                r$analysis <- analysis
                r$family_id <- paste(feature$family, analysis, sep = ":")
                r$test_id <- paste(key, feature$feature_id, analysis, exposure, sep = ":")
                results[[length(results) + 1L]] <- r
            }
            next
        }
        job <- bi_checkpoint(ctx, paste0("continuous_", key, "_", i), function() tryCatch({
            feature_qc <- NULL
            if (nzchar(feature$qc_column)) {
                allowed <- strsplit(as.character(feature$qc_allowed), "[;,]")[[1L]]
                if (!length(allowed) || any(!nzchar(allowed))) 
                  stop("Explicit feature QC pass codes required")
                v <- bi_cont_lookup(list(source = feature$qc_source, column = feature$qc_column), dat$cohort, dat$features, 
                  dat$index, "feature reportability/QC")
                feature_qc <- !is.na(v) & as.character(v) %in% allowed
            }
            prepared <- bi_cont_prepare(dat$features[[feature$column]][dat$index], feature, dat$keep, feature_qc)
            d <- dat$data
            d$y <- prepared$y
            extra <- dat$extra
            specific <- bi_cont_value(cfg$feature_metadata[[feature$feature_id]], list())
            specific_mappings <- c(bi_cont_value(specific$batches, list()), bi_cont_value(specific$extra_covariates, 
                list()))
            for (j in seq_along(specific_mappings)) {
                nm <- paste0("bi_assay_technical_", j)
                d[[nm]] <- bi_cont_technical(specific_mappings[[j]], dat$cohort, dat$features, dat$index, nm)
                extra <- c(extra, nm)
            }
            output <- selfcheck <- audits <- list()
            for (analysis in analyses) {
                keep <- rep(TRUE, nrow(d))
                reason <- ""
                current <- d
                raw <- prepared$raw
                f <- feature
                if (!analysis %in% c("primary", "bsa_indexed")) {
                  spec <- sensitivities[[match(analysis, vapply(sensitivities, function(s) as.character(s$id), 
                    character(1)))]]
                  rule <- bi_cont_sensitivity(ctx, spec, dat, cfg)
                  if (!rule$available) 
                    reason <- rule$reason
                  else keep <- rule$keep
                }
                if (analysis == "bsa_indexed") {
                  m <- cfg$metadata$bsa
                  if (is.null(m) || !bi_cont_flag(m$verified) || !identical(m$unit, "m2")) 
                    reason <- "Verified same-visit BSA in m2 unavailable"
                  else {
                    bsa <- bi_cont_lookup(m, dat$cohort, dat$features, dat$index, "same-visit BSA")
                    if (!is.numeric(bsa)) 
                      stop("BSA must be numeric m2")
                    bsa[!is.finite(bsa) | bsa <= 0] <- NA_real_
                    indexed <- prepared$raw/bsa
                    f$unit <- paste0(feature$unit, "/m2")
                    f$valid_min <- f$valid_max <- f$winsor_lower <- f$winsor_upper <- NA_real_
                    f$missing_codes <- ""
                    indexed <- bi_cont_prepare(indexed, f, dat$keep)
                    current$y <- indexed$y
                    raw <- indexed$raw
                  }
                }
                audits[[length(audits) + 1L]] <- data.frame(feature_id = feature$feature_id, analysis = analysis, 
                  input_N = sum(!is.na(d$y)), retained_N = if (nzchar(reason)) 
                    NA_integer_
                  else sum(keep & !is.na(current$y)), status = if (nzchar(reason)) 
                    "SKIPPED_UNAVAILABLE"
                  else "PASS", reason = reason)
                fits <- bi_cont_map(exposures, function(exposure) {
                  if (nzchar(reason)) {
                    r <- bi_cont_empty_result(exposure, "SKIPPED_UNAVAILABLE", reason)
                    r$trait <- r$feature_id <- feature$feature_id
                    r$analysis <- analysis
                    r$family_id <- paste(feature$family, analysis, sep = ":")
                    r$test_id <- paste(key, feature$feature_id, analysis, exposure, sep = ":")
                    list(result = r, selfcheck = data.frame())
                  }
                  else {
                    bi_cont_fit_one(ctx, current[keep, , drop = FALSE], exposure, extra, f, raw[keep], key, analysis, 
                      check = i <= first_n && analysis == "primary")
                  }
                }, if (nzchar(reason)) 
                  1L
                else workers$effective)
                for (fit in fits) {
                  output[[length(output) + 1L]] <- fit$result
                  if (i <= first_n && analysis == "primary" && nrow(fit$selfcheck)) 
                    selfcheck[[length(selfcheck) + 1L]] <- fit$selfcheck
                }
            }
            list(result = bi_bind(output), qc = prepared$qc, checks = bi_bind(selfcheck), audit = bi_bind(audits))
        }, error = function(e) {
            failed <- list()
            for (a in analyses) for (ex in exposures) {
                r <- bi_cont_empty_result(ex, "FEATURE_ERROR", conditionMessage(e))
                r$trait <- r$feature_id <- feature$feature_id
                r$analysis <- a
                r$family_id <- paste(feature$family, a, sep = ":")
                r$test_id <- paste(key, feature$feature_id, a, ex, sep = ":")
                failed[[length(failed) + 1L]] <- r
            }
            list(result = bi_bind(failed), qc = data.frame(feature_id = feature$feature_id, status = "FEATURE_ERROR", 
                error = conditionMessage(e)), checks = if (i <= first_n) data.frame(feature_id = feature$feature_id, 
                exposure = exposures, status = "FAIL", reason = conditionMessage(e), relative_max_difference = NA_real_) else data.frame(), 
                audit = data.frame(feature_id = feature$feature_id, analysis = analyses, status = "FEATURE_ERROR", 
                  reason = conditionMessage(e)))
        }))
        results[[length(results) + 1L]] <- job$result
        qc[[length(qc) + 1L]] <- job$qc
        if (nrow(job$checks)) 
            checks[[length(checks) + 1L]] <- job$checks
        sensitivity_audit[[length(sensitivity_audit) + 1L]] <- job$audit
        if (i%%chunk_size == 0L || i == first_n || i == nrow(planned)) 
            bi_log(ctx, "continuous_chunk_complete", module = key, completed_features = i, planned_features = nrow(planned))
    }
    table <- bi_bind(results)
    table$data_snapshot_id <- cfg$source_version
    table$measurement_visit <- cfg$metadata$visit
    table$exposure_unit <- unname(linear_units[table$exposure])
    dictionary_index <- match(table$feature_id, planned$feature_id)
    for (nm in c("label", "primary_id", "kind", "transform")) {
        if (!nm %in% names(table)) 
            table[[nm]] <- NA_character_
        absent <- is.na(table[[nm]])
        table[[nm]][absent] <- planned[[nm]][dictionary_index[absent]]
    }
    if (!"original_unit" %in% names(table)) 
        table$original_unit <- NA_character_
    absent <- is.na(table$original_unit)
    table$original_unit[absent] <- planned$unit[dictionary_index[absent]]
    check_table <- bi_bind(checks)
    numeric_failure <- !nrow(check_table) || any(check_table$status == "FAIL") || !any(check_table$status == "PASS")
    if (numeric_failure) {
        table$p <- NA_real_
        table$status[table$status %in% c("PASS", "OK")] <- "NUMERIC_SELFCHECK_REVIEW"
        table$inference <- "NUMERIC_SELFCHECK_FAILED_NO_PUBLISHED_INFERENCE"
    }
    table <- bi_adjust(ctx, table)
    secondary <- table
    secondary$family_id <- paste(secondary$family_id, secondary$exposure, "secondary", sep = ":")
    secondary <- bi_adjust(ctx, secondary)
    if (!"q" %in% names(table) || !"q" %in% names(secondary)) 
        stop("bi_adjust must return column q")
    table$q_branch_secondary <- secondary$q[match(table$test_id, secondary$test_id)]
    bi_cont_save(ctx, key, "results.csv", table)
    bi_cont_save(ctx, key, "feature_qc.csv", bi_bind(qc))
    bi_cont_save(ctx, key, "numeric_selfcheck.csv", bi_bind(checks))
    bi_cont_save(ctx, key, "sensitivity_flow.csv", bi_bind(sensitivity_audit))
    if (key == "glycaemia") {
        bi_cont_save(ctx, key, "glucose_hba1c_results.csv", table[grepl("^glycaemia_primary:", table$family_id), 
            , drop = FALSE])
        bi_cont_save(ctx, key, "biochemistry_results.csv", table[!grepl("^glycaemia_primary:", table$family_id), 
            , drop = FALSE])
    }
    failures <- sum(!table$status %in% c("PASS", "OK"))
    partial <- failures > 0 || numeric_failure || any(dict$include & nzchar(dict$exclusion_reason))
    list(status = if (partial) "PARTIAL" else "PASS", reason = paste("features", nrow(planned), "planned tests", 
        nrow(table), "non-passing", failures, if (numeric_failure) "numeric selfcheck unavailable/failed" else "numeric selfcheck passed"), 
        features = nrow(planned), planned_tests = nrow(table))
}

bi_stage_glycaemia <- function(ctx) bi_feature_stage(ctx, "glycaemia")

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        bi_cont_value <- function(x, default = NULL) if (is.null(x) || !length(x)) default else x
        bi_cont_flag <- function(x) isTRUE(x) || identical(as.character(x), "TRUE")
        bi_cont_save <- function(ctx, key, path, table) {
            directory <- if (key == "glycaemia") "biochemistry" else key
            bi_save(ctx, file.path("raw", directory, path), table)
        }
        bi_cont_analyses <- function(feature, sensitivities) {
            selected <- Filter(function(s) is.null(s$features) || feature$feature_id %in% unlist(s$features), sensitivities)
            ids <- c("primary", vapply(selected, function(s) as.character(s$id), character(1)))
            if (bi_cont_flag(feature$bsa_index)) ids <- c(ids, "bsa_indexed")
            if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Sensitivity IDs must be nonempty and unique")
            ids
        }
        bi_cont_workers <- function(ctx, cfg, data, shared_gb = 0) {
            requested <- as.integer(bi_cont_value(cfg$max_workers, bi_cont_value(ctx$cfg$resources$max_workers, 
                1L)))
            if (length(requested) != 1L || !is.finite(requested) || requested < 1L) stop("max_workers must be a positive integer")
            effective <- min(requested, 4L)
            reason <- "Forked exposure fits share one read-only feature matrix; master writes checkpoints"
            estimate <- max(0.25, as.numeric(object.size(data)) * 10/1024^3)
            reserve <- as.numeric(bi_cont_value(ctx$cfg$resources$minimum_free_gb, 2))
            ram_limit <- as.numeric(bi_cont_value(ctx$cfg$resources$ram_limit_gb, 24))
            if (!is.finite(ram_limit) || ram_limit <= 0 || !is.finite(shared_gb) || shared_gb < 0) stop("Invalid continuous RAM budget")
            cap <- floor((ram_limit - shared_gb - reserve)/estimate)
            if (cap < 1L) stop("CONTINUOUS_RAM_LIMIT_EXCEEDED_BY_SHARED_INPUTS")
            effective <- min(effective, cap)
            if (.Platform$OS.type != "unix") {
                effective <- 1L
                reason <- "Native Windows cannot fork; use reviewed WSL for parallel execution"
            } else if (file.exists("/proc/meminfo")) {
                line <- readLines("/proc/meminfo", warn = FALSE)
                available <- as.numeric(gsub("[^0-9]", "", line[grepl("^MemAvailable:", line)]))/1024^2
                if (length(available) != 1L || !is.finite(available)) stop("Cannot validate available RAM for fork workers")
                cap <- floor((available - reserve)/estimate)
                if (cap < 1L) stop("INSUFFICIENT_FREE_MEMORY_FOR_CONTINUOUS_WORKERS")
                effective <- min(effective, cap)
            }
            list(requested = requested, effective = effective, estimated_worker_gb = estimate, shared_input_gb = shared_gb, 
                ram_limit_gb = ram_limit, reason = reason)
        }
        bi_cont_map <- function(x, fun, workers) {
            if (workers <= 1L) return(lapply(x, fun))
            ans <- parallel::mclapply(x, fun, mc.cores = min(workers, length(x)), mc.preschedule = FALSE, mc.set.seed = FALSE, 
                mc.cleanup = TRUE)
            if (any(vapply(ans, inherits, logical(1), "try-error"))) stop("Continuous worker error: ", paste(unlist(ans[vapply(ans, 
                inherits, logical(1), "try-error")]), collapse = "; "))
            ans
        }
        bi_cont_mapping_columns <- function(x, source) {
            if (!is.list(x)) return(character())
            here <- if (identical(x$source, source) && !is.null(x$column)) as.character(x$column) else character()
            unique(c(here, unlist(lapply(x, bi_cont_mapping_columns, source = source), use.names = FALSE)))
        }
        bi_cont_source_dialect <- function(path, ctx = NULL) {
            answer <- list(dialect = "standard", id_column = "eid")
            entries <- if (is.null(ctx)) list() else bi_cont_value(ctx$cfg$source_dialects, list())
            if (!is.list(entries)) stop("source_dialects must be a list of input/dialect/id_column records")
            if (!length(entries)) return(answer)
            paths <- vapply(entries, function(entry) {
                if (!is.list(entry) || length(entry$input) != 1L || is.na(entry$input) || !nzchar(entry$input)) stop("source_dialects entry requires one resolved input path")
                if (length(entry$dialect) != 1L || is.na(entry$dialect) || !entry$dialect %in% c("standard", "backslash", 
                  "literal")) stop("Invalid explicit source dialect")
                if (length(entry$id_column) != 1L || is.na(entry$id_column) || !nzchar(entry$id_column)) stop("source_dialects entry requires one ID column")
                normalizePath(entry$input, winslash = "/", mustWork = FALSE)
            }, character(1))
            if (anyDuplicated(paths)) stop("Duplicate resolved inputs in source_dialects")
            selected <- which(paths == normalizePath(path, winslash = "/", mustWork = TRUE))
            if (length(selected)) answer <- entries[[selected]][c("dialect", "id_column")]
            answer
        }
        bi_cont_read <- function(path, columns, ctx = NULL) {
            if (!file.exists(path)) stop("Continuous source missing: ", path)
            columns <- unique(as.character(columns[!is.na(columns) & nzchar(columns)]))
            if (grepl("[.]rds$", path, ignore.case = TRUE)) {
                d <- bi_read_table(path)
                return(d[, intersect(columns, names(d)), drop = FALSE])
            }
            dialect <- bi_cont_source_dialect(path, ctx)
            if (grepl("[.]gz$", path, ignore.case = TRUE) || dialect$dialect != "standard") {
                root <- if (!is.null(ctx) && !is.null(ctx$root)) ctx$root else Sys.getenv("BI_ROOT")
                helper <- file.path(root, "f", "streaming_tsv_select.py")
                if (length(root) != 1L || !nzchar(root) || !file.exists(helper)) stop("Streaming gzip projection helper unavailable; set BI_ROOT or ctx$root")
                python <- Sys.getenv("BI_PYTHON", unset = "python3")
                if (!nzchar(Sys.which(python)) && !file.exists(python)) stop("Python 3 required for streaming gzip projection")
                delimiter <- if (grepl("[.]csv([.]gz)?$", path, ignore.case = TRUE)) "comma" else "tab"
                args <- c(python, helper, "--input", normalizePath(path, mustWork = TRUE), "--columns-json", jsonlite::toJSON(columns, 
                  auto_unbox = FALSE), "--delimiter", delimiter, "--dialect", dialect$dialect, "--id-column", dialect$id_column, 
                  "--allow-missing", "--progress-every", "50000")
                command <- paste(vapply(args, shQuote, character(1)), collapse = " ")
                header_error <- tempfile("projection_header_diagnostic_")
                body_error <- tempfile("projection_body_diagnostic_")
                completion <- tempfile("projection_completion_", fileext = ".json")
                on.exit(unlink(c(header_error, body_error, completion)), add = TRUE)
                diagnostic <- function(p) if (file.exists(p)) paste(readLines(p, warn = FALSE), collapse = "\n") else ""
                header_text <- suppressWarnings(system(paste(command, shQuote("--header-only"), "2>", shQuote(header_error)), 
                  intern = TRUE))
                header_status <- attr(header_text, "status")
                if (!is.null(header_status) && header_status != 0L) stop("SOURCE_PROJECTION_FAILED: ", diagnostic(header_error), 
                  call. = FALSE)
                header <- names(utils::read.csv(text = paste(header_text, collapse = "\n"), nrows = 0L, check.names = FALSE))
                body_command <- paste(command, shQuote("--audit-json"), shQuote(completion), shQuote("--diagnostic-log"), 
                  shQuote(body_error))
                consumer_warnings <- character()
                projected <- tryCatch(withCallingHandlers(as.data.frame(data.table::fread(cmd = body_command, sep = ",", 
                  header = TRUE, col.names = header, colClasses = if (dialect$id_column %in% header) setNames("character", 
                    dialect$id_column) else NULL, check.names = FALSE, data.table = FALSE, showProgress = FALSE)), 
                  warning = function(w) {
                    consumer_warnings <<- c(consumer_warnings, conditionMessage(w))
                    invokeRestart("muffleWarning")
                  }), error = identity)
                details <- diagnostic(body_error)
                if (!file.exists(completion) || grepl("[continuous_projection] ERROR", details, fixed = TRUE)) {
                  if (!nzchar(details)) details <- "Parser completion audit absent; source stream did not complete"
                  stop("SOURCE_PROJECTION_FAILED: ", details, call. = FALSE)
                }
                audit <- tryCatch(jsonlite::read_json(completion, simplifyVector = TRUE), error = function(e) NULL)
                if (is.null(audit) || !identical(audit$status, "PASS") || !isTRUE(audit$full_scan) || isTRUE(audit$header_only) || 
                  isTRUE(audit$validation_only) || !is.null(audit$max_records) || !identical(audit$dialect, dialect$dialect) || 
                  !identical(audit$id_column, dialect$id_column) || length(audit$rows) != 1L || !is.numeric(audit$rows) || 
                  !is.finite(audit$rows) || audit$rows < 0 || length(audit$unique_ids) != 1L || audit$rows != audit$unique_ids) {
                  stop("SOURCE_PROJECTION_FAILED: invalid completion audit; action=verify_helper_version_and_rerun_read_validation", 
                    call. = FALSE)
                }
                if (inherits(projected, "error") || length(consumer_warnings)) {
                  stop("PROJECTED_CSV_CONSUMER_FAILED after successful source validation; error_class=", if (inherits(projected, 
                    "error")) class(projected)[1L] else "none", " warning_count=", length(consumer_warnings), "; no participant output accepted; action=verify_R_data_table_adapter", 
                    call. = FALSE)
                }
                if (nrow(projected) != audit$rows || ncol(projected) != audit$selected_columns || !identical(names(projected), 
                  header)) stop("PROJECTED_CSV_SHAPE_MISMATCH", call. = FALSE)
                quote_probe <- data.table::fread(text = "x\n\"a\"\"\"\"b\"\n", sep = ",", header = TRUE, colClasses = "character", 
                  data.table = FALSE)[[1L]]
                if (identical(quote_probe, "a\"\"\"\"b")) {
                  for (nm in names(projected)[vapply(projected, is.character, logical(1))]) {
                    projected[[nm]] <- gsub("\"\"", "\"", projected[[nm]], fixed = TRUE)
                  }
                } else if (!identical(quote_probe, "a\"\"b")) stop("PROJECTED_CSV_QUOTE_DECODING_UNSUPPORTED", 
                  call. = FALSE)
                if (dialect$id_column %in% names(projected)) {
                  ids <- projected[[dialect$id_column]]
                  if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) stop("PROJECTED_CSV_INVALID_ID_KEYS", 
                    call. = FALSE)
                }
                return(projected)
            }
            header <- names(data.table::fread(path, nrows = 0L, check.names = FALSE, data.table = FALSE, showProgress = FALSE))
            selected <- intersect(columns, header)
            if (!length(selected)) stop("No configured columns present in continuous source")
            as.data.frame(data.table::fread(path, select = selected, check.names = FALSE, data.table = FALSE, showProgress = FALSE))
        }
        bi_cont_numbers <- function(x) {
            if (is.null(x) || !length(x) || all(is.na(x)) || identical(x, "")) return(numeric())
            z <- suppressWarnings(as.numeric(unlist(strsplit(paste(x, collapse = ";"), "[;,]"))))
            if (anyNA(z)) stop("Non-numeric configured missing code or bound")
            z
        }
        bi_cont_dictionary <- function(x) {
            if (is.character(x) && length(x) == 1L) x <- bi_read_table(x)
            if (is.data.frame(x)) {
                rows <- lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
            } else if (is.list(x)) rows <- x else stop("dictionary must be a table path or feature records")
            if (!length(rows)) stop("No feature dictionary records")
            required <- c("feature_id", "column", "label", "unit", "family", "kind", "include", "qc_verified", 
                "transform")
            defaults <- list(primary_id = "", missing_codes = "", valid_min = NA_real_, valid_max = NA_real_, winsor_lower = NA_real_, 
                winsor_upper = NA_real_, qc_column = "", qc_allowed = "", qc_source = "features", numeric_parse = "require_numeric", 
                bsa_index = FALSE, visit = "")
            out <- lapply(rows, function(r) {
                absent <- setdiff(required, names(r))
                if (length(absent)) stop("Dictionary missing fields: ", paste(absent, collapse = ","))
                for (nm in names(defaults)) if (is.null(r[[nm]]) || !length(r[[nm]])) r[[nm]] <- defaults[[nm]]
                for (nm in c("primary_id", "missing_codes", "qc_column", "qc_allowed", "qc_source", "visit")) {
                  if (length(r[[nm]]) != 1L && nm != "missing_codes") stop("Optional dictionary field must be scalar: ", 
                    nm)
                  if (all(is.na(r[[nm]]))) r[[nm]] <- ""
                }
                for (nm in c("feature_id", "column", "label", "unit", "family", "kind", "transform")) {
                  if (length(r[[nm]]) != 1L || is.na(r[[nm]]) || !nzchar(as.character(r[[nm]]))) stop("Empty dictionary field: ", 
                    nm)
                  r[[nm]] <- as.character(r[[nm]])
                }
                r$include <- bi_cont_flag(r$include)
                r$qc_verified <- bi_cont_flag(r$qc_verified)
                r$bsa_index <- bi_cont_flag(r$bsa_index)
                if (!r$transform %in% c("identity", "log", "log1p", "rank_normal")) stop("Unsupported locked transformation: ", 
                  r$transform)
                if (!r$numeric_parse %in% c("require_numeric", "non_numeric_to_missing")) stop("Unsupported explicit numeric parsing rule")
                if (!r$kind %in% c("primary", "derived")) stop("kind must be primary or derived")
                if (!nzchar(as.character(r$primary_id))) r$primary_id <- r$feature_id
                r$missing_codes <- paste(bi_cont_numbers(r$missing_codes), collapse = ";")
                for (nm in c("valid_min", "valid_max", "winsor_lower", "winsor_upper")) {
                  original <- r[[nm]]
                  r[[nm]] <- suppressWarnings(as.numeric(r[[nm]]))
                  if (length(r[[nm]]) != 1L) stop("Bound must be scalar: ", nm)
                  if (is.na(r[[nm]]) && !is.na(original) && nzchar(as.character(original))) stop("Invalid numeric dictionary bound: ", 
                    nm)
                }
                if (is.finite(r$valid_min) && is.finite(r$valid_max) && r$valid_min >= r$valid_max) stop("Invalid feature QC range")
                w <- c(r$winsor_lower, r$winsor_upper)
                if (any(is.finite(w)) && (!all(is.finite(w)) || w[1L] < 0 || w[2L] > 1 || w[1L] >= w[2L])) stop("Invalid winsor quantiles")
                r[unique(c(required, names(defaults)))]
            })
            ans <- bi_bind(lapply(out, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
            for (nm in c("feature_id", "column", "primary_id")) if (anyDuplicated(ans[[nm]][ans$include])) stop("Duplicate included ", 
                nm)
            ans
        }
        bi_cont_gate <- function(cfg) {
            if (is.null(cfg) || !bi_cont_flag(cfg$enabled)) return("Module not enabled in locked configuration")
            if (!bi_cont_flag(cfg$verified)) return("Source/dictionary linkage has not been verified")
            if (is.null(cfg$source_version) || !nzchar(cfg$source_version)) return("Verified source version is required")
            if (!bi_cont_flag(cfg$metadata$verified) || is.null(cfg$metadata$visit) || !nzchar(cfg$metadata$visit)) return("Measurement visit metadata has not been verified")
            if (!bi_cont_flag(cfg$qc$verified) || is.null(cfg$qc$description) || !nzchar(cfg$qc$description)) return("Assay/image QC and missing-code semantics have not been verified")
            if (!identical(cfg$qc$input_imputation, "none")) return("Unimputed input provenance has not been verified; preimputed values cannot be silently treated as observed measurements")
            if (is.null(cfg$dictionary)) return("Feature dictionary is not connected")
            for (nm in c("age", "center", "year")) if (is.null(cfg$metadata[[nm]]$column)) return(paste("Visit-specific mapping missing:", 
                nm))
            NULL
        }
        bi_cont_lookup <- function(mapping, cohort, features, index, label) {
            metadata <- attr(features, "bi_metadata_tables")
            if (is.null(mapping$source) || !mapping$source %in% c("cohort", "features", names(metadata)) || is.null(mapping$column)) stop("Invalid explicit mapping: ", 
                label)
            src <- if (mapping$source == "cohort") cohort else if (mapping$source == "features") features else metadata[[mapping$source]]$data
            if (!mapping$column %in% names(src)) stop("Configured column absent for ", label, ": ", mapping$column)
            if (mapping$source == "cohort") src[[mapping$column]] else if (mapping$source == "features") src[[mapping$column]][index] else src[[mapping$column]][metadata[[mapping$source]]$index]
        }
        bi_cont_rule <- function(rule, cohort, features, index) {
            z <- bi_cont_lookup(rule, cohort, features, index, "QC/filter")
            if (is.null(rule$allowed) || !length(rule$allowed)) stop("QC/filter allowed values must be explicit")
            !is.na(z) & as.character(z) %in% as.character(unlist(rule$allowed))
        }
        bi_cont_technical <- function(mapping, cohort, features, index, label) {
            z <- bi_cont_lookup(mapping, cohort, features, index, label)
            if (identical(mapping$type, "factor")) return(factor(z))
            if (identical(mapping$type, "date_month")) {
                if (inherits(z, "POSIXt")) z <- as.Date(z, tz = "UTC")
                if (!inherits(z, "Date")) {
                  if (!is.character(z)) stop("Assay dates must be Date or ISO date strings: ", mapping$column)
                  parsed <- suppressWarnings(as.Date(z, format = "%Y-%m-%d"))
                  if (any(!is.na(z) & nzchar(z) & is.na(parsed))) stop("Unparseable configured assay dates: ", 
                    mapping$column)
                  z <- parsed
                }
                return(factor(format(z, "%Y-%m")))
            }
            if (!identical(mapping$type, "numeric") || !is.numeric(z)) stop("Technical mapping must specify factor, date_month, or an actual numeric column: ", 
                mapping$column)
            z[!is.finite(z)] <- NA_real_
            if (!is.null(mapping$missing_codes)) z[z %in% bi_cont_numbers(mapping$missing_codes)] <- NA_real_
            if (!is.null(mapping$valid_min)) z[z < as.numeric(mapping$valid_min)] <- NA_real_
            if (!is.null(mapping$valid_max)) z[z > as.numeric(mapping$valid_max)] <- NA_real_
            z
        }
        bi_cont_load <- function(ctx, cfg, key) {
            cohort <- bi_cohort(ctx)
            input <- bi_cont_value(cfg$input, "cohort")
            dict <- bi_cont_dictionary(cfg$dictionary)
            id <- bi_cont_value(cfg$id_column, "eid")
            columns <- unique(c(if (id != "rownames") id, dict$column[dict$include], dict$qc_column[dict$include & 
                dict$qc_source == "features"], bi_cont_mapping_columns(cfg, "features"), unlist(cfg$select_columns)))
            features <- if (identical(input, "cohort")) cohort else bi_cont_read(input, columns, ctx)
            if (is.matrix(features)) features <- as.data.frame(features, stringsAsFactors = FALSE)
            if (!is.data.frame(features)) stop("Continuous input must be an explicitly keyed wide table")
            ids <- if (identical(id, "rownames")) rownames(features) else {
                if (!id %in% names(features)) stop("Feature ID column missing: ", id)
                as.character(features[[id]])
            }
            if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Missing or duplicate feature participant keys")
            if (!"eid" %in% names(cohort) || anyNA(cohort$eid) || anyDuplicated(cohort$eid)) stop("Invalid cohort participant keys")
            index <- match(as.character(cohort$eid), ids)
            metadata <- list()
            tables <- bi_cont_value(cfg$metadata_tables, list())
            if (length(tables) && (is.null(names(tables)) || any(!nzchar(names(tables))) || anyDuplicated(names(tables)) || 
                any(names(tables) %in% c("cohort", "features")))) stop("Metadata table source names must be explicit, unique and nonreserved")
            for (name in names(tables)) {
                spec <- tables[[name]]
                if (!bi_cont_flag(spec$verified) || is.null(spec$input) || is.null(spec$id_column)) stop("Unverified metadata source: ", 
                  name)
                selected <- unique(c(spec$id_column, unlist(spec$select_columns), bi_cont_mapping_columns(cfg, 
                  name), dict$qc_column[dict$include & dict$qc_source == name]))
                md <- bi_cont_read(spec$input, selected, ctx)
                if (!spec$id_column %in% names(md)) stop("Metadata participant key absent: ", name)
                mid <- as.character(md[[spec$id_column]])
                if (anyNA(mid) || any(!nzchar(mid)) || anyDuplicated(mid)) stop("Duplicate/missing metadata participant keys: ", 
                  name)
                metadata[[name]] <- list(data = md, index = match(as.character(cohort$eid), mid))
            }
            attr(features, "bi_metadata_tables") <- metadata
            model_fields <- c("eid", "border12", "border1n", "border_rank", "border_relative", "S", "age", "sex", 
                "ethnicity", "center", "assessment_year")
            if (length(setdiff(model_fields, names(cohort)))) stop("Frozen continuous model fields missing")
            d <- cohort[, model_fields, drop = FALSE]
            d$age <- bi_cont_lookup(cfg$metadata$age, cohort, features, index, "measurement age")
            if (!is.numeric(d$age)) stop("Verified measurement age must be numeric years")
            d$age[!is.finite(d$age)] <- NA_real_
            d$center <- factor(bi_cont_lookup(cfg$metadata$center, cohort, features, index, "measurement center"))
            d$assessment_year <- factor(bi_cont_lookup(cfg$metadata$year, cohort, features, index, "measurement year"))
            extra <- character()
            mappings <- c(bi_cont_value(cfg$metadata$batches, list()), bi_cont_value(cfg$metadata$extra_covariates, 
                list()))
            for (i in seq_along(mappings)) {
                m <- mappings[[i]]
                nm <- paste0("bi_technical_", i)
                d[[nm]] <- bi_cont_technical(m, cohort, features, index, nm)
                extra <- c(extra, nm)
            }
            keep <- !is.na(index)
            for (r in bi_cont_value(cfg$qc$sample_rules, list())) keep <- keep & bi_cont_rule(r, cohort, features, 
                index)
            dict$exclusion_reason <- ifelse(!dict$include, "EXCLUDED_PREDECLARED", ifelse(!dict$qc_verified, "QC_UNVERIFIED", 
                ifelse(!dict$column %in% names(features), "COLUMN_UNAVAILABLE", "")))
            bad_visit <- nzchar(dict$visit) & dict$visit != cfg$metadata$visit
            dict$exclusion_reason[bad_visit] <- "VISIT_MISMATCH"
            if (key %in% c("nmr", "olink")) dict$family <- paste0(key, "_all")
            list(cohort = cohort, features = features, index = index, data = d, extra = extra, keep = keep, dictionary = dict, 
                flow = data.frame(module = key, cohort_N = nrow(cohort), source_N = length(ids), matched_N = sum(!is.na(index)), 
                  unmatched_source_N = sum(!ids %in% as.character(cohort$eid)), sample_qc_pass_N = sum(keep), measurement_visit = cfg$metadata$visit, 
                  source_version = cfg$source_version, stringsAsFactors = FALSE))
        }
        bi_cont_prepare <- function(v, feature, keep, feature_qc = NULL) {
            invalid_text <- 0L
            if (!is.numeric(v)) {
                if (!identical(feature$numeric_parse, "non_numeric_to_missing") || !is.character(v)) stop("Configured continuous feature is not numeric: ", 
                  feature$column)
                parsed <- suppressWarnings(as.numeric(v))
                invalid_text <- sum(!is.na(v) & nzchar(v) & is.na(parsed))
                v <- parsed
            }
            n_input <- sum(!is.na(v))
            v[!is.finite(v)] <- NA_real_
            v[v %in% bi_cont_numbers(feature$missing_codes)] <- NA_real_
            if (is.finite(feature$valid_min)) v[v < feature$valid_min] <- NA_real_
            if (is.finite(feature$valid_max)) v[v > feature$valid_max] <- NA_real_
            if (!is.null(feature_qc)) keep <- keep & feature_qc
            v[!keep] <- NA_real_
            if (all(is.finite(c(feature$winsor_lower, feature$winsor_upper))) && any(!is.na(v))) {
                q <- stats::quantile(v, probs = c(feature$winsor_lower, feature$winsor_upper), na.rm = TRUE, names = FALSE, 
                  type = 7)
                v <- pmax(q[1L], pmin(q[2L], v))
            }
            raw <- v
            invalid <- rep(FALSE, length(v))
            if (feature$transform == "log") invalid <- !is.na(v) & v <= 0
            if (feature$transform == "log1p") invalid <- !is.na(v) & v <= -1
            v[invalid] <- NA_real_
            raw[invalid] <- NA_real_
            if (feature$transform == "log") v <- log(v)
            if (feature$transform == "log1p") v <- log1p(v)
            if (feature$transform == "rank_normal") {
                ok <- !is.na(v)
                v[ok] <- stats::qnorm((rank(v[ok], ties.method = "average") - 0.5)/sum(ok))
            }
            list(y = v, raw = raw, qc = data.frame(feature_id = feature$feature_id, input_nonmissing_N = n_input, 
                retained_measurement_N = sum(!is.na(v)), transformation_invalid_N = sum(invalid), invalid_numeric_text_N = invalid_text, 
                transformation = feature$transform, original_unit = feature$unit, stringsAsFactors = FALSE))
        }
        bi_cont_good <- function(r) nrow(r) == 1L && r$status %in% c("PASS", "OK") && is.finite(r$beta) && is.finite(r$se) && 
            is.finite(r$p)
        bi_cont_reference <- function(ctx, d, exposure, extra, result) {
            if (!bi_cont_good(result)) {
                expected <- grepl("INSUFFICIENT_CONTINUOUS_N|CONSTANT_EXPOSURE", as.character(result$error))
                return(list(status = if (expected) "NOT_ESTIMABLE" else "FAIL", reason = paste(result$status, result$error), 
                  delta = NA_real_))
            }
            form <- bi_formula(exposure, extra = extra)
            s <- bi_sample(d, form, exposure)
            reference <- stats::lm(s$formula, data = s$data, na.action = stats::na.fail)
            co <- summary(reference)$coefficients
            if (!exposure %in% rownames(co)) return(list(status = "FAIL", reason = "Reference coefficient absent", 
                delta = NA_real_))
            expected <- c(co[exposure, 1L], co[exposure, 2L], co[exposure, 4L])
            observed <- c(result$beta, result$se, result$p)
            delta <- max(abs(expected - observed)/pmax(1, abs(expected)))
            sd_y <- stats::sd(s$data$y)
            z <- s$data
            z$y <- z$y/sd_y
            scaled <- summary(stats::lm(s$formula, data = z, na.action = stats::na.fail))$coefficients[exposure, 
                1L]
            correct <- is.finite(delta) && delta < 1e-07 && nrow(s$data) == result$N && is.finite(sd_y) && sd_y > 
                0 && abs(scaled - result$beta/sd_y) < 1e-07
            list(status = if (correct) "PASS" else "FAIL", reason = if (correct) "Gaussian/lm and standardized-beta equality" else "Independent numeric/sample mismatch", 
                delta = delta)
        }
        bi_cont_empty_result <- function(exposure, status, reason) {
            data.frame(exposure = exposure, N = NA_integer_, cases = NA_integer_, noncases = NA_integer_, beta = NA_real_, 
                se = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_, p_raw = NA_real_, 
                status = status, warning = "", error = reason, formula = "", rank = NA_integer_, columns = NA_integer_, 
                sample_hash = "", stringsAsFactors = FALSE)
        }
        bi_cont_fit_one <- function(ctx, d, exposure, extra, feature, raw, key, analysis, check) {
            result <- bi_fit(ctx, d, exposure, extra = extra, family = "gaussian")
            result$trait <- feature$feature_id
            result$feature_id <- feature$feature_id
            result$label <- feature$label
            result$primary_id <- feature$primary_id
            result$analysis <- analysis
            result$kind <- feature$kind
            result$original_unit <- feature$unit
            result$transform <- feature$transform
            result$effect_unit <- if (feature$transform == "identity") feature$unit else paste0(feature$transform, 
                "(", feature$unit, ")")
            result$family_id <- paste(feature$family, analysis, sep = ":")
            result$test_id <- paste(key, feature$feature_id, analysis, exposure, sep = ":")
            result$sd_outcome <- NA_real_
            result$standardized_beta <- result$standardized_CI_low <- result$standardized_CI_high <- NA_real_
            result$original_beta <- result$original_CI_low <- result$original_CI_high <- result$original_p <- NA_real_
            result$original_status <- as.character(result$status)
            if (bi_cont_good(result)) {
                sample <- bi_sample(d, bi_formula(exposure, extra = extra), exposure)$data
                result$sd_outcome <- stats::sd(sample$y)
                if (is.finite(result$sd_outcome) && result$sd_outcome > 0) {
                  result$standardized_beta <- result$beta/result$sd_outcome
                  result$standardized_CI_low <- result$CI_low/result$sd_outcome
                  result$standardized_CI_high <- result$CI_high/result$sd_outcome
                }
                raw_result <- result
                if (feature$transform != "identity") {
                  raw_d <- d
                  raw_d$y <- raw
                  raw_result <- bi_fit(ctx, raw_d, exposure, extra = extra, family = "gaussian")
                  if (bi_cont_good(raw_result) && (!identical(as.character(raw_result$sample_hash), as.character(result$sample_hash)) || 
                    raw_result$N != result$N)) stop("Original/transformed fit sample mismatch")
                }
                result$original_status <- as.character(raw_result$status)
                if (bi_cont_good(raw_result)) {
                  result$original_beta <- raw_result$beta
                  result$original_CI_low <- raw_result$CI_low
                  result$original_CI_high <- raw_result$CI_high
                  result$original_p <- raw_result$p
                }
            }
            reference <- if (check) bi_cont_reference(ctx, d, exposure, extra, result) else list(status = "NOT_REQUESTED", 
                reason = "", delta = NA_real_)
            list(result = result, selfcheck = data.frame(feature_id = feature$feature_id, exposure = exposure, 
                status = reference$status, reason = reference$reason, relative_max_difference = reference$delta, 
                stringsAsFactors = FALSE))
        }
        bi_cont_sensitivity <- function(ctx, spec, dat, cfg) {
            if (!bi_cont_flag(spec$verified)) return(list(available = FALSE, reason = "Sensitivity field/definition not verified"))
            if (identical(spec$type, "outcome_absent")) {
                y <- bi_outcomes(ctx)
                if (is.null(spec$code) || !spec$code %in% names(y)) return(list(available = FALSE, reason = "Outcome status unavailable"))
                if (anyDuplicated(y$eid)) stop("Duplicate outcome participant keys")
                value <- y[[spec$code]][match(dat$cohort$eid, y$eid)]
                if (any(!is.na(value) & !value %in% c(0, 1))) stop("Outcome must preserve verified 0/1/NA")
                keep <- !is.na(value) & value == 0
            } else if (identical(spec$type, "keep_values")) {
                keep <- bi_cont_rule(spec, dat$cohort, dat$features, dat$index)
            } else stop("Unsupported predeclared sensitivity type")
            list(available = TRUE, keep = keep, reason = "")
        }
        bi_cont_selection <- function(ctx, key, cfg, dat) {
            ans <- list()
            if (!is.null(cfg$metadata$baseline_date) && !is.null(cfg$metadata$measurement_date)) {
                baseline <- bi_cont_lookup(cfg$metadata$baseline_date, dat$cohort, dat$features, dat$index, "baseline date")
                measurement <- bi_cont_lookup(cfg$metadata$measurement_date, dat$cohort, dat$features, dat$index, 
                  "measurement date")
                if (!inherits(baseline, "Date") || !inherits(measurement, "Date")) stop("Verified selection dates must use Date class")
                delay <- as.numeric(measurement - baseline)/365.25
                good <- dat$keep & is.finite(delay) & delay >= 0
                ans[[length(ans) + 1L]] <- data.frame(metric = "baseline_to_measurement_years", N = sum(good), 
                  mean = if (any(good)) mean(delay[good]) else NA_real_, sd = if (sum(good) > 1) stats::sd(delay[good]) else NA_real_, 
                  invalid_negative_N = sum(dat$keep & !is.na(delay) & delay < 0), status = "PASS", stringsAsFactors = FALSE)
            } else if (key == "imaging") ans[[length(ans) + 1L]] <- data.frame(metric = "baseline_to_measurement_years", 
                N = NA_integer_, mean = NA_real_, sd = NA_real_, invalid_negative_N = NA_integer_, status = "UNAVAILABLE_DATES")
            for (nm in names(bi_cont_value(cfg$metadata$prior_disease, list()))) {
                m <- cfg$metadata$prior_disease[[nm]]
                v <- bi_cont_lookup(m, dat$cohort, dat$features, dat$index, paste("prior disease", nm))
                if (any(!is.na(v) & !v %in% c(0, 1))) stop("Prior disease status must use verified 0/1/NA")
                ans[[length(ans) + 1L]] <- data.frame(metric = paste0("disease_present_at_measurement_", nm), N = sum(dat$keep & 
                  !is.na(v)), cases = sum(dat$keep & !is.na(v) & v == 1), status = "PASS")
            }
            if (length(ans)) bi_cont_save(ctx, key, "selection_qc.csv", bi_bind(ans))
        }
        bi_feature_stage <- function(ctx, key) {
            cfg <- bi_module_config(ctx, key)
            gate <- bi_cont_gate(cfg)
            if (!is.null(gate)) return(bi_skip(ctx, gate))
            exposures <- c("border12", "border1n", "border_rank", "border_relative")
            loaded <- tryCatch(bi_cont_load(ctx, cfg, key), error = function(e) e)
            if (inherits(loaded, "error")) return(list(status = "FAILED", reason = paste("Adapter validation:", 
                conditionMessage(loaded))))
            dat <- loaded
            dict <- dat$dictionary
            bi_cont_save(ctx, key, "feature_dictionary.csv", dict)
            bi_cont_save(ctx, key, "cohort_flow.csv", dat$flow)
            bi_cont_selection(ctx, key, cfg, dat)
            chosen <- which(dict$exclusion_reason == "")
            if (!length(chosen)) return(bi_skip(ctx, "No predeclared, QC-verified, available features"))
            chunk_size <- as.integer(bi_cont_value(cfg$chunk_size, 25L))
            if (length(chunk_size) != 1L || !is.finite(chunk_size) || chunk_size < 1L || chunk_size > 1000L) stop("chunk_size must be 1..1000")
            limit <- as.numeric(bi_cont_value(cfg$max_elapsed_seconds, 86400))
            if (length(limit) != 1L || !is.finite(limit) || limit <= 0) stop("max_elapsed_seconds must be positive and finite")
            shared_gb <- (as.numeric(object.size(dat$cohort)) + as.numeric(object.size(dat$features)) + as.numeric(object.size(dat$data)))/1024^3
            workers <- bi_cont_workers(ctx, cfg, dat$data, shared_gb)
            bi_log(ctx, "continuous_resources", requested_workers = workers$requested, effective_workers = workers$effective, 
                estimated_worker_gb = workers$estimated_worker_gb, shared_input_gb = workers$shared_input_gb, ram_limit_gb = workers$ram_limit_gb, 
                reason = workers$reason, chunk_size = chunk_size)
            started <- proc.time()[["elapsed"]]
            results <- qc <- checks <- sensitivity_audit <- list()
            planned <- dict[chosen, , drop = FALSE]
            first_n <- min(10L, nrow(planned))
            stopped <- ""
            sensitivities <- bi_cont_value(cfg$sensitivities, list())
            registry <- list()
            for (i in seq_len(nrow(planned))) {
                f <- as.list(planned[i, , drop = FALSE])
                for (a in bi_cont_analyses(f, sensitivities)) for (e in exposures) {
                  registry[[length(registry) + 1L]] <- data.frame(feature_id = f$feature_id, analysis = a, exposure = e, 
                    family_id = paste(f$family, a, sep = ":"), test_id = paste(key, f$feature_id, a, e, sep = ":"))
                }
            }
            bi_cont_save(ctx, key, "planned_tests.csv", bi_bind(registry))
            for (i in seq_len(nrow(planned))) {
                feature <- as.list(planned[i, , drop = FALSE])
                if (i > first_n && !nzchar(stopped) && i == first_n + 1L) {
                  ctab <- bi_bind(checks)
                  if (any(ctab$status == "FAIL") || !any(ctab$status == "PASS")) stopped <- "NUMERIC_SELFCHECK_FAILED"
                }
                if (!nzchar(stopped) && proc.time()[["elapsed"]] - started > limit) stop("MODULE_TIME_LIMIT: completed feature caches preserved; resume with the same configuration")
                analyses <- bi_cont_analyses(feature, sensitivities)
                if (nzchar(stopped)) {
                  for (analysis in analyses) for (exposure in exposures) {
                    r <- bi_cont_empty_result(exposure, stopped, "Planned row retained; fitting not attempted")
                    r$trait <- r$feature_id <- feature$feature_id
                    r$analysis <- analysis
                    r$family_id <- paste(feature$family, analysis, sep = ":")
                    r$test_id <- paste(key, feature$feature_id, analysis, exposure, sep = ":")
                    results[[length(results) + 1L]] <- r
                  }
                  next
                }
                job <- bi_checkpoint(ctx, paste0("continuous_", key, "_", i), function() tryCatch({
                  feature_qc <- NULL
                  if (nzchar(feature$qc_column)) {
                    allowed <- strsplit(as.character(feature$qc_allowed), "[;,]")[[1L]]
                    if (!length(allowed) || any(!nzchar(allowed))) stop("Explicit feature QC pass codes required")
                    v <- bi_cont_lookup(list(source = feature$qc_source, column = feature$qc_column), dat$cohort, 
                      dat$features, dat$index, "feature reportability/QC")
                    feature_qc <- !is.na(v) & as.character(v) %in% allowed
                  }
                  prepared <- bi_cont_prepare(dat$features[[feature$column]][dat$index], feature, dat$keep, feature_qc)
                  d <- dat$data
                  d$y <- prepared$y
                  extra <- dat$extra
                  specific <- bi_cont_value(cfg$feature_metadata[[feature$feature_id]], list())
                  specific_mappings <- c(bi_cont_value(specific$batches, list()), bi_cont_value(specific$extra_covariates, 
                    list()))
                  for (j in seq_along(specific_mappings)) {
                    nm <- paste0("bi_assay_technical_", j)
                    d[[nm]] <- bi_cont_technical(specific_mappings[[j]], dat$cohort, dat$features, dat$index, nm)
                    extra <- c(extra, nm)
                  }
                  output <- selfcheck <- audits <- list()
                  for (analysis in analyses) {
                    keep <- rep(TRUE, nrow(d))
                    reason <- ""
                    current <- d
                    raw <- prepared$raw
                    f <- feature
                    if (!analysis %in% c("primary", "bsa_indexed")) {
                      spec <- sensitivities[[match(analysis, vapply(sensitivities, function(s) as.character(s$id), 
                        character(1)))]]
                      rule <- bi_cont_sensitivity(ctx, spec, dat, cfg)
                      if (!rule$available) reason <- rule$reason else keep <- rule$keep
                    }
                    if (analysis == "bsa_indexed") {
                      m <- cfg$metadata$bsa
                      if (is.null(m) || !bi_cont_flag(m$verified) || !identical(m$unit, "m2")) reason <- "Verified same-visit BSA in m2 unavailable" else {
                        bsa <- bi_cont_lookup(m, dat$cohort, dat$features, dat$index, "same-visit BSA")
                        if (!is.numeric(bsa)) stop("BSA must be numeric m2")
                        bsa[!is.finite(bsa) | bsa <= 0] <- NA_real_
                        indexed <- prepared$raw/bsa
                        f$unit <- paste0(feature$unit, "/m2")
                        f$valid_min <- f$valid_max <- f$winsor_lower <- f$winsor_upper <- NA_real_
                        f$missing_codes <- ""
                        indexed <- bi_cont_prepare(indexed, f, dat$keep)
                        current$y <- indexed$y
                        raw <- indexed$raw
                      }
                    }
                    audits[[length(audits) + 1L]] <- data.frame(feature_id = feature$feature_id, analysis = analysis, 
                      input_N = sum(!is.na(d$y)), retained_N = if (nzchar(reason)) NA_integer_ else sum(keep & 
                        !is.na(current$y)), status = if (nzchar(reason)) "SKIPPED_UNAVAILABLE" else "PASS", reason = reason)
                    fits <- bi_cont_map(exposures, function(exposure) {
                      if (nzchar(reason)) {
                        r <- bi_cont_empty_result(exposure, "SKIPPED_UNAVAILABLE", reason)
                        r$trait <- r$feature_id <- feature$feature_id
                        r$analysis <- analysis
                        r$family_id <- paste(feature$family, analysis, sep = ":")
                        r$test_id <- paste(key, feature$feature_id, analysis, exposure, sep = ":")
                        list(result = r, selfcheck = data.frame())
                      } else {
                        bi_cont_fit_one(ctx, current[keep, , drop = FALSE], exposure, extra, f, raw[keep], key, 
                          analysis, check = i <= first_n && analysis == "primary")
                      }
                    }, if (nzchar(reason)) 1L else workers$effective)
                    for (fit in fits) {
                      output[[length(output) + 1L]] <- fit$result
                      if (i <= first_n && analysis == "primary" && nrow(fit$selfcheck)) selfcheck[[length(selfcheck) + 
                        1L]] <- fit$selfcheck
                    }
                  }
                  list(result = bi_bind(output), qc = prepared$qc, checks = bi_bind(selfcheck), audit = bi_bind(audits))
                }, error = function(e) {
                  failed <- list()
                  for (a in analyses) for (ex in exposures) {
                    r <- bi_cont_empty_result(ex, "FEATURE_ERROR", conditionMessage(e))
                    r$trait <- r$feature_id <- feature$feature_id
                    r$analysis <- a
                    r$family_id <- paste(feature$family, a, sep = ":")
                    r$test_id <- paste(key, feature$feature_id, a, ex, sep = ":")
                    failed[[length(failed) + 1L]] <- r
                  }
                  list(result = bi_bind(failed), qc = data.frame(feature_id = feature$feature_id, status = "FEATURE_ERROR", 
                    error = conditionMessage(e)), checks = if (i <= first_n) data.frame(feature_id = feature$feature_id, 
                    exposure = exposures, status = "FAIL", reason = conditionMessage(e), relative_max_difference = NA_real_) else data.frame(), 
                    audit = data.frame(feature_id = feature$feature_id, analysis = analyses, status = "FEATURE_ERROR", 
                      reason = conditionMessage(e)))
                }))
                results[[length(results) + 1L]] <- job$result
                qc[[length(qc) + 1L]] <- job$qc
                if (nrow(job$checks)) checks[[length(checks) + 1L]] <- job$checks
                sensitivity_audit[[length(sensitivity_audit) + 1L]] <- job$audit
                if (i%%chunk_size == 0L || i == first_n || i == nrow(planned)) bi_log(ctx, "continuous_chunk_complete", 
                  module = key, completed_features = i, planned_features = nrow(planned))
            }
            table <- bi_bind(results)
            table$data_snapshot_id <- cfg$source_version
            table$measurement_visit <- cfg$metadata$visit
            table$exposure_unit <- unname(linear_units[table$exposure])
            dictionary_index <- match(table$feature_id, planned$feature_id)
            for (nm in c("label", "primary_id", "kind", "transform")) {
                if (!nm %in% names(table)) table[[nm]] <- NA_character_
                absent <- is.na(table[[nm]])
                table[[nm]][absent] <- planned[[nm]][dictionary_index[absent]]
            }
            if (!"original_unit" %in% names(table)) table$original_unit <- NA_character_
            absent <- is.na(table$original_unit)
            table$original_unit[absent] <- planned$unit[dictionary_index[absent]]
            check_table <- bi_bind(checks)
            numeric_failure <- !nrow(check_table) || any(check_table$status == "FAIL") || !any(check_table$status == 
                "PASS")
            if (numeric_failure) {
                table$p <- NA_real_
                table$status[table$status %in% c("PASS", "OK")] <- "NUMERIC_SELFCHECK_REVIEW"
                table$inference <- "NUMERIC_SELFCHECK_FAILED_NO_PUBLISHED_INFERENCE"
            }
            table <- bi_adjust(ctx, table)
            secondary <- table
            secondary$family_id <- paste(secondary$family_id, secondary$exposure, "secondary", sep = ":")
            secondary <- bi_adjust(ctx, secondary)
            if (!"q" %in% names(table) || !"q" %in% names(secondary)) stop("bi_adjust must return column q")
            table$q_branch_secondary <- secondary$q[match(table$test_id, secondary$test_id)]
            bi_cont_save(ctx, key, "results.csv", table)
            bi_cont_save(ctx, key, "feature_qc.csv", bi_bind(qc))
            bi_cont_save(ctx, key, "numeric_selfcheck.csv", bi_bind(checks))
            bi_cont_save(ctx, key, "sensitivity_flow.csv", bi_bind(sensitivity_audit))
            if (key == "glycaemia") {
                bi_cont_save(ctx, key, "glucose_hba1c_results.csv", table[grepl("^glycaemia_primary:", table$family_id), 
                  , drop = FALSE])
                bi_cont_save(ctx, key, "biochemistry_results.csv", table[!grepl("^glycaemia_primary:", table$family_id), 
                  , drop = FALSE])
            }
            failures <- sum(!table$status %in% c("PASS", "OK"))
            partial <- failures > 0 || numeric_failure || any(dict$include & nzchar(dict$exclusion_reason))
            list(status = if (partial) "PARTIAL" else "PASS", reason = paste("features", nrow(planned), "planned tests", 
                nrow(table), "non-passing", failures, if (numeric_failure) "numeric selfcheck unavailable/failed" else "numeric selfcheck passed"), 
                features = nrow(planned), planned_tests = nrow(table))
        }
        bi_stage_glycaemia <- function(ctx) bi_feature_stage(ctx, "glycaemia")
    }), envir = context)
}

