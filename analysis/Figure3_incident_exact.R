# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

bi_isa_codes <- function() c("411", "411.3", "411.4", "351", "250", "250.2")

bi_isa_exposures <- function() c("border12", "border1n", "border_rank", "border_relative")

bi_isa_ids <- function(x, label) {
    bi_assert(is.atomic(x) && !is.list(x), paste("INCIDENT_AUDIT_INVALID_IDS", label))
    if (is.numeric(x)) 
        bi_assert(all(is.finite(x) & x == floor(x) & abs(x) <= 2^53), paste("INCIDENT_AUDIT_INEXACT_IDS", label))
    z <- as.character(x)
    bi_assert(!anyNA(z) && all(nzchar(z)) && all(z == trimws(z)) && !anyDuplicated(z), paste("INCIDENT_AUDIT_DUPLICATE_OR_MISSING_IDS", 
        label))
    z
}

bi_isa_dates <- function(x) {
    if (inherits(x, "Date")) 
        return(as.Date(x))
    bi_assert(!is.numeric(x) && !is.logical(x), "INCIDENT_AUDIT_UNTYPED_DATE")
    s <- as.character(x)
    missing <- is.na(s) | s == ""
    valid <- missing | grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", s)
    bi_assert(all(valid), "INCIDENT_AUDIT_NON_ISO_DATE")
    ans <- as.Date(s, format = "%Y-%m-%d")
    bi_assert(all(missing | (!is.na(ans) & format(ans, "%Y-%m-%d") == s)), "INCIDENT_AUDIT_INVALID_CALENDAR_DATE")
    ans
}

bi_isa_normalize <- function(x) {
    s <- toupper(trimws(as.character(x)))
    s[is.na(s) | s %in% c("", "NA", "NULL", "NAN")] <- NA_character_
    s <- gsub(".", "", s, fixed = TRUE)
    bi_assert(all(is.na(s) | grepl("^[A-Z0-9]+$", s)), "INCIDENT_AUDIT_MALFORMED_DIAGNOSIS_CODE")
    s
}

bi_isa_canonical_hash <- function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    bi_assert(!anyNA(x), "INCIDENT_AUDIT_MAP_MISSING_VALUE")
    x <- x[do.call(order, c(lapply(x, as.character), list(method = "radix"))), , drop = FALSE]
    digest::digest(paste(apply(x, 1L, paste, collapse = "\t"), collapse = "\n"), algo = "sha256", serialize = FALSE)
}

bi_isa_map <- function(direct, rollup) {
    bi_assert(identical(names(direct), c("vocabulary_id", "code", "phecode")) && identical(names(rollup), c("code", 
        "phecode_unrolled")), "INCIDENT_AUDIT_MAP_SCHEMA_MISMATCH")
    bi_assert(!anyNA(direct) && !anyNA(rollup), "INCIDENT_AUDIT_MAP_MISSING_VALUE")
    z <- merge(direct, rollup, by.x = "phecode", by.y = "code", sort = FALSE)
    z <- z[z$vocabulary_id == "ICD10CM" & z$phecode_unrolled %in% bi_isa_codes(), , drop = FALSE]
    keys <- gsub(".", "", toupper(trimws(as.character(z$code))), fixed = TRUE)
    bi_assert(all(grepl("^[A-Z0-9]+(-[A-Z0-9]+)?$", keys)), "INCIDENT_AUDIT_UNSUPPORTED_DICTIONARY_TOKEN")
    z <- unique(data.frame(trait = as.character(z$phecode_unrolled), source_code = as.character(z$code), normalized_code = keys, 
        stringsAsFactors = FALSE))
    bi_assert(setequal(unique(z$trait), bi_isa_codes()), "INCIDENT_AUDIT_SIX_MAPS_REQUIRED")
    z
}

bi_isa_pair <- function(code, dates) {
    bi_assert(is.data.frame(code) && is.data.frame(dates) && ncol(code) >= 2L && identical(dim(code), dim(dates)) && 
        identical(names(code)[1], "eid") && identical(names(dates)[1], "eid"), "INCIDENT_AUDIT_PAIR_SCHEMA_MISMATCH")
    ci <- bi_isa_ids(code$eid, "code")
    di <- bi_isa_ids(dates$eid, "date")
    bi_assert(identical(ci, di), "INCIDENT_AUDIT_CODE_DATE_ID_ORDER_MISMATCH")
    bi_assert(!anyDuplicated(names(code)) && !anyDuplicated(names(dates)) && all(grepl("^i[0-9]+$", names(code)[-1])) && 
        all(grepl("^d[0-9]+$", names(dates)[-1])) && identical(sub("^i", "", names(code)[-1]), sub("^d", "", names(dates)[-1])), 
        "INCIDENT_AUDIT_ARRAY_POSITION_MISMATCH")
    invisible(ci)
}

bi_isa_field_summary <- function(cohort, fields) {
    bi_assert("date_attend" %in% names(cohort), "INCIDENT_AUDIT_ENTRY_FIELD_MISSING")
    at <- bi_isa_dates(cohort$date_attend)
    bi_bind(lapply(fields, function(f) {
        present <- f %in% names(cohort)
        tt <- if (present) 
            bi_isa_dates(cohort[[f]])
        else rep(as.Date(NA), nrow(cohort))
        data.frame(field = f, present = present, cohort_n = nrow(cohort), nonmissing_count = if (present) 
            sum(!is.na(tt))
        else NA_integer_, on_or_before_entry_count = if (present) 
            sum(tt <= at, na.rm = TRUE)
        else NA_integer_, after_entry_unfiltered_count = if (present) 
            sum(tt > at, na.rm = TRUE)
        else NA_integer_, limitation = "Candidate source only; not PheCode-specific incidence, risk-set events or follow-up coverage", 
            stringsAsFactors = FALSE)
    }))
}

bi_isa_extract <- function(code, dates, mapping) {
    bi_isa_pair(code, dates)
    out <- list()
    for (k in bi_isa_codes()) {
        allowed <- unique(mapping$normalized_code[mapping$trait == k])
        bi_assert(length(allowed) > 0L, "INCIDENT_AUDIT_EMPTY_EXACT_MAP")
        first <- rep(NA_real_, nrow(code))
        observed <- undated <- logical(nrow(code))
        pairs <- 0L
        for (j in seq.int(2L, ncol(code))) {
            cc <- bi_isa_normalize(code[[j]])
            hit <- !is.na(cc) & cc %in% allowed
            if (!any(hit)) 
                next
            tt <- rep(NA_real_, nrow(code))
            tt[hit] <- as.numeric(bi_isa_dates(dates[[j]][hit]))
            observed <- observed | hit
            undated <- undated | (hit & is.na(tt))
            pairs <- pairs + sum(hit)
            first <- pmin(first, tt, na.rm = TRUE)
            first[!is.finite(first)] <- NA_real_
        }
        out[[k]] <- list(first = as.Date(first, origin = "1970-01-01"), observed = observed, undated = undated, 
            pairs = pairs)
    }
    out
}

bi_isa_candidate_adapter <- function(cohort, code, dates, mapping) {
    ids <- bi_isa_pair(code, dates)
    aid <- bi_isa_ids(cohort$eid, "adapter cohort")
    at <- bi_isa_dates(cohort$date_attend)
    ix <- match(aid, ids)
    events <- bi_isa_extract(code, dates, mapping)
    out <- list()
    for (k in bi_isa_codes()) {
        z <- events[[k]]
        known <- !is.na(ix) & !z$undated[ix]
        known[is.na(known)] <- FALSE
        event <- z$first[ix]
        event[!known] <- as.Date(NA)
        seen <- ifelse(known, as.integer(z$observed[ix]), NA_integer_)
        prior <- ifelse(known & !is.na(at), as.integer(!is.na(event) & event <= at), NA_integer_)
        slug <- gsub(".", "_", k, fixed = TRUE)
        p <- paste0("phe_", slug)
        dd <- data.frame(eid = aid, stringsAsFactors = FALSE)
        dd[[paste0(p, "_first_hospital_date")]] <- event
        dd[[paste0(p, "_hospital_event_observed")]] <- seen
        dd[[paste0(p, "_hospital_prior_only")]] <- prior
        dd[[paste0(p, "_approved_all_source_baseline")]] <- rep(NA_integer_, length(aid))
        out[[k]] <- dd
    }
    out
}

bi_isa_config_template <- function() {
    outcomes <- setNames(lapply(bi_isa_codes(), function(k) {
        p <- paste0("phe_", gsub(".", "_", k, fixed = TRUE))
        list(enabled = TRUE, definition_verified = FALSE, event_dates_verified = FALSE, baseline_exclusion_verified = FALSE, 
            followup_coverage_verified = FALSE, definition_id = paste0("CANDIDATE_EXACT_HOSPITAL_PHECODE_", k), 
            source_id = "ICD10_PHEWAS_EXACT_DIRECT_PLUS_ROLLUP_AUDIT_ONLY", review_id = "", definition_manifest = "", 
            event_semantics = "first_ever_recorded", baseline_semantics = "approved_all_available_history", event_date_field = paste0(p, 
                "_first_hospital_date"), event_observed_field = paste0(p, "_hospital_event_observed"), baseline_disease_field = paste0(p, 
                "_approved_all_source_baseline"), input_path = "", id_field = "eid", coverage_start_field = "", 
            coverage_end_field = "", admin_censor_field = "", warning = "Template only: hospital_prior_only must NOT populate approved_all_source_baseline; preserve unverified gates")
    }), bi_isa_codes())
    list(enabled = TRUE, protocol_verified = FALSE, protocol_id = "", time_scale = "attained_age", age_time_scale_verified = FALSE, 
        censoring_verified = FALSE, date_precision_verified = FALSE, date_review_id = "", date_fields = list(birth = "birth_date", 
            entry = "date_attend", death = "", loss = "date_lost"), birth_date_rule = "", tie_rule = "event_on_censor_counts", 
        sensitivity_s8 = TRUE, outcomes = outcomes, warning = "Raw death source and region-specific completeness must be reviewed; no current FOD maximum date is a cutoff")
}

bi_isa_scan <- function(cohort, code, dates, flags, mapping) {
    ids <- bi_isa_pair(code, dates)
    aid <- bi_isa_ids(cohort$eid, "frozen all")
    bi_assert("id" %in% names(flags), "INCIDENT_AUDIT_PHEWAS_ID_FIELD_REQUIRED")
    pid <- bi_isa_ids(flags$id, "PheWAS")
    ix <- match(ids, aid)
    px <- match(ids, pid)
    bi_assert(all(bi_isa_codes() %in% names(flags)), "INCIDENT_AUDIT_PHEWAS_COLUMNS_MISSING")
    at <- bi_isa_dates(cohort$date_attend)[ix]
    rows <- list()
    events <- bi_isa_extract(code, dates, mapping)
    for (k in bi_isa_codes()) {
        mm <- mapping[mapping$trait == k, , drop = FALSE]
        allowed <- unique(mm$normalized_code)
        z <- events[[k]]
        first <- as.numeric(z$first)
        observed <- z$observed
        undated <- z$undated
        pairs <- z$pairs
        old <- flags[[k]][px]
        known <- !is.na(px)
        bi_assert(all(is.na(old) | old %in% c(0, 1)), "INCIDENT_AUDIT_INVALID_PHEWAS_FLAGS")
        rows[[k]] <- data.frame(trait = k, source = "hospital_ICD10_exact_PheWAS_direct_plus_rollup", source_rows = nrow(code), 
            source_positions = ncol(code) - 1L, cohort_n = nrow(cohort), source_cohort_overlap_count = sum(!is.na(ix)), 
            cohort_without_source_count = sum(!aid %in% ids), map_rows = nrow(mm), unique_exact_codes = length(allowed), 
            positive_people_count = sum(observed), positive_code_date_pair_count = pairs, dated_positive_people_count = sum(!is.na(first)), 
            undated_positive_people_count = sum(undated), on_or_before_entry_count = sum(first <= as.numeric(at), 
                na.rm = TRUE), after_entry_unfiltered_count = sum(first > as.numeric(at), na.rm = TRUE), phewas_compared_count = sum(known), 
            phewas_positive_count = sum(old == 1, na.rm = TRUE), phewas_flag_disagreement_count = sum(observed[known] != 
                (old[known] == 1)), status = "AUDIT_ONLY_NOT_APPROVED", limitation = "Post-baseline raw first-code count only; no baseline all-source exclusion, censoring, risk set or regression", 
            stringsAsFactors = FALSE)
    }
    bi_bind(rows)
}

bi_isa_blockers <- function(m, code) {
    o <- m$outcomes[[code]] %||% list()
    checks <- c(PROTOCOL_NOT_APPROVED = !isTRUE(m$protocol_verified), EXACT_OUTCOME_DEFINITION_NOT_APPROVED = !isTRUE(o$definition_verified), 
        EVENT_DATE_ADAPTER_NOT_APPROVED = !isTRUE(o$event_dates_verified), BASELINE_ALL_SOURCE_EXCLUSION_NOT_APPROVED = !isTRUE(o$baseline_exclusion_verified), 
        SOURCE_REGION_COMPLETENESS_NOT_APPROVED = !isTRUE(o$followup_coverage_verified), CENSORING_NOT_APPROVED = !isTRUE(m$censoring_verified), 
        DATE_PRECISION_NOT_APPROVED = !isTRUE(m$date_precision_verified), ATTAINED_AGE_NOT_APPROVED = !isTRUE(m$age_time_scale_verified), 
        DEFINITION_MANIFEST_MISSING = !bi_temporal_text(o$definition_manifest), EVENT_STATUS_BASELINE_BINDINGS_MISSING = !all(vapply(list(o$event_date_field, 
            o$event_observed_field, o$baseline_disease_field), bi_temporal_text, logical(1))))
    paste(names(checks)[checks], collapse = ";")
}

bi_isa_nofit_matrix <- function(m) {
    z <- expand.grid(exposure = bi_isa_exposures(), trait = bi_isa_codes(), variant = c("primary", "S_le_8"), stringsAsFactors = FALSE)
    z$test_id <- paste(z$exposure, z$trait, z$variant, sep = "__")
    z$status <- "AUDIT_ONLY_NOT_FITTED"
    z$blockers <- vapply(z$trait, function(k) bi_isa_blockers(m, k), character(1))
    z$HR <- z$CI_low <- z$CI_high <- z$p <- z$q <- NA_real_
    z$release_allowed <- FALSE
    z$interpretation <- "Not estimated: neither a negative association nor Cox failure; this audit cannot approve a fit"
    bi_assert(nrow(z) == 48L && !anyDuplicated(z$test_id), "INCIDENT_AUDIT_48_SLOT_MATRIX_FAILURE")
    z
}

bi_stage_incident_source_audit <- function(ctx) {
    bi_assert(ctx$stage %in% c("incident_audit", "incident_source_audit") || isTRUE(ctx$synthetic), "INCIDENT_SOURCE_AUDIT_MANUAL_STAGE_REQUIRED")
    bi_assert(isTRUE(ctx$real_authorized) || isTRUE(ctx$synthetic), "REAL_DATA_ACCESS_NOT_AUTHORIZED")
    m <- bi_module_config(ctx, "incident")
    ac <- bi_module_config(ctx, "incident_source_audit")
    invpath <- file.path(ctx$root, "config", "incident_source_inventory.csv")
    bi_assert(file.exists(invpath), paste("INCIDENT_SOURCE_INVENTORY_MISSING", invpath))
    inventory <- read.csv(invpath, colClasses = "character", check.names = FALSE)
    bi_assert(identical(inventory$trait, bi_isa_codes()), "INCIDENT_SOURCE_INVENTORY_SCOPE_MISMATCH")
    allpath <- ctx$cfg$inputs$all_rds
    phepath <- ctx$cfg$inputs$phewas_rds
    rdir <- dirname(phepath)
    paths <- c(all_rds = allpath, phewas_rds = phepath, icd10_code = ac$code_path %||% file.path(rdir, "icd10.code0.rds"), 
        icd10_date = ac$date_path %||% file.path(rdir, "icd10.date0.rds"))
    cat("INCIDENT_SOURCE_AUDIT_ROOT=", ctx$root, "\nALL_RDS=", allpath, "\nPHE_RDATA=", rdir, "\nOUT=", ctx$out, 
        "\nAUDIT_ONLY_NO_MODELS\n", sep = "")
    bi_assert(all(file.exists(paths)), paste("INCIDENT_AUDIT_MISSING_INPUT", paste(paths[!file.exists(paths)], 
        collapse = ",")))
    dest <- normalizePath(ctx$out, mustWork = FALSE)
    for (p in paths) {
        protected <- normalizePath(dirname(p), mustWork = TRUE)
        bi_assert(!identical(dest, protected) && !startsWith(dest, paste0(protected, "/")), "INCIDENT_AUDIT_OUTPUT_OVERLAPS_INPUT")
    }
    before <- lapply(paths, function(p) list(bytes = unname(file.info(p)$size), mtime = unname(as.numeric(file.info(p)$mtime)), 
        sha256 = bi_hash(p, TRUE)))
    expected <- ctx$cfg$input_sha256
    for (k in intersect(names(expected), names(paths))) bi_assert(identical(before[[k]]$sha256, expected[[k]]), 
        paste("INCIDENT_AUDIT_INPUT_SHA256_MISMATCH", k))
    bi_resource_guard(ctx)
    cohort <- as.data.frame(readRDS(allpath))
    bi_isa_ids(cohort$eid, "all")
    fields <- unique(c("birth_date", "date_attend", "date_death", "date_lost", unlist(strsplit(inventory$candidate_fod_fields, 
        ";", fixed = TRUE))))
    field_audit <- bi_isa_field_summary(cohort, fields[nzchar(fields)])
    cohort <- cohort[, c("eid", "date_attend"), drop = FALSE]
    bi_assert(requireNamespace("PheWAS", quietly = TRUE), "INCIDENT_AUDIT_PHEWAS_PACKAGE_REQUIRED")
    maps <- new.env(parent = emptyenv())
    utils::data(list = c("phecode_map", "phecode_rollup_map"), package = "PheWAS", envir = maps)
    direct <- as.data.frame(maps$phecode_map)
    rollup <- as.data.frame(maps$phecode_rollup_map)
    mh <- bi_isa_canonical_hash(direct)
    rh <- bi_isa_canonical_hash(rollup)
    bi_assert(identical(mh, birth_private_setting("PRIVATE_HASH_1")) && identical(rh, birth_private_setting("PRIVATE_HASH_2")), 
        "INCIDENT_AUDIT_PHEWAS_MAP_CHANGED_REVIEW_REQUIRED")
    mapping <- bi_isa_map(direct, rollup)
    pair_audit <- bi_isa_scan(cohort, as.data.frame(readRDS(paths["icd10_code"])), as.data.frame(readRDS(paths["icd10_date"])), 
        as.data.frame(readRDS(phepath)), mapping)
    pair_audit$map_sha256 <- mh
    pair_audit$rollup_sha256 <- rh
    for (k in names(paths)) {
        now <- file.info(paths[[k]])
        bi_assert(identical(unname(now$size), before[[k]]$bytes) && identical(unname(as.numeric(now$mtime)), before[[k]]$mtime), 
            "INCIDENT_AUDIT_SOURCE_CHANGED_DURING_READ")
    }
    source_audit <- data.frame(source = names(paths), path = unname(paths), bytes = vapply(before, `[[`, numeric(1), 
        "bytes"), sha256 = vapply(before, `[[`, character(1), "sha256"), stringsAsFactors = FALSE)
    bi_save(ctx, "raw/incident_source_audit/source_manifest.csv", source_audit)
    bi_save(ctx, "raw/incident_source_audit/source_inventory.csv", inventory)
    bi_save(ctx, "raw/incident_source_audit/field_availability.csv", field_audit)
    bi_save(ctx, "raw/incident_source_audit/exact_phecode_date_counts.csv", pair_audit)
    bi_save(ctx, "raw/incident_source_audit/nofit_matrix.csv", bi_isa_nofit_matrix(m))
    bi_save(ctx, "raw/incident_source_audit/unapproved_config_template.json", jsonlite::toJSON(bi_isa_config_template(), 
        auto_unbox = TRUE, pretty = TRUE, null = "null"))
    list(status = "AUDIT_ONLY_NOT_FITTED", reason = "48 nofit slots retained; source agreement/counts cannot approve outcome, baseline exclusion or censoring; no risk set or model created")
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        bi_isa_codes <- function() c("411", "411.3", "411.4", "351", "250", "250.2")
        bi_isa_exposures <- function() c("border12", "border1n", "border_rank", "border_relative")
        bi_isa_ids <- function(x, label) {
            bi_assert(is.atomic(x) && !is.list(x), paste("INCIDENT_AUDIT_INVALID_IDS", label))
            if (is.numeric(x)) bi_assert(all(is.finite(x) & x == floor(x) & abs(x) <= 2^53), paste("INCIDENT_AUDIT_INEXACT_IDS", 
                label))
            z <- as.character(x)
            bi_assert(!anyNA(z) && all(nzchar(z)) && all(z == trimws(z)) && !anyDuplicated(z), paste("INCIDENT_AUDIT_DUPLICATE_OR_MISSING_IDS", 
                label))
            z
        }
        bi_isa_dates <- function(x) {
            if (inherits(x, "Date")) return(as.Date(x))
            bi_assert(!is.numeric(x) && !is.logical(x), "INCIDENT_AUDIT_UNTYPED_DATE")
            s <- as.character(x)
            missing <- is.na(s) | s == ""
            valid <- missing | grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", s)
            bi_assert(all(valid), "INCIDENT_AUDIT_NON_ISO_DATE")
            ans <- as.Date(s, format = "%Y-%m-%d")
            bi_assert(all(missing | (!is.na(ans) & format(ans, "%Y-%m-%d") == s)), "INCIDENT_AUDIT_INVALID_CALENDAR_DATE")
            ans
        }
        bi_isa_normalize <- function(x) {
            s <- toupper(trimws(as.character(x)))
            s[is.na(s) | s %in% c("", "NA", "NULL", "NAN")] <- NA_character_
            s <- gsub(".", "", s, fixed = TRUE)
            bi_assert(all(is.na(s) | grepl("^[A-Z0-9]+$", s)), "INCIDENT_AUDIT_MALFORMED_DIAGNOSIS_CODE")
            s
        }
        bi_isa_canonical_hash <- function(x) {
            x <- as.data.frame(x, stringsAsFactors = FALSE)
            bi_assert(!anyNA(x), "INCIDENT_AUDIT_MAP_MISSING_VALUE")
            x <- x[do.call(order, c(lapply(x, as.character), list(method = "radix"))), , drop = FALSE]
            digest::digest(paste(apply(x, 1L, paste, collapse = "\t"), collapse = "\n"), algo = "sha256", serialize = FALSE)
        }
        bi_isa_map <- function(direct, rollup) {
            bi_assert(identical(names(direct), c("vocabulary_id", "code", "phecode")) && identical(names(rollup), 
                c("code", "phecode_unrolled")), "INCIDENT_AUDIT_MAP_SCHEMA_MISMATCH")
            bi_assert(!anyNA(direct) && !anyNA(rollup), "INCIDENT_AUDIT_MAP_MISSING_VALUE")
            z <- merge(direct, rollup, by.x = "phecode", by.y = "code", sort = FALSE)
            z <- z[z$vocabulary_id == "ICD10CM" & z$phecode_unrolled %in% bi_isa_codes(), , drop = FALSE]
            keys <- gsub(".", "", toupper(trimws(as.character(z$code))), fixed = TRUE)
            bi_assert(all(grepl("^[A-Z0-9]+(-[A-Z0-9]+)?$", keys)), "INCIDENT_AUDIT_UNSUPPORTED_DICTIONARY_TOKEN")
            z <- unique(data.frame(trait = as.character(z$phecode_unrolled), source_code = as.character(z$code), 
                normalized_code = keys, stringsAsFactors = FALSE))
            bi_assert(setequal(unique(z$trait), bi_isa_codes()), "INCIDENT_AUDIT_SIX_MAPS_REQUIRED")
            z
        }
        bi_isa_pair <- function(code, dates) {
            bi_assert(is.data.frame(code) && is.data.frame(dates) && ncol(code) >= 2L && identical(dim(code), dim(dates)) && 
                identical(names(code)[1], "eid") && identical(names(dates)[1], "eid"), "INCIDENT_AUDIT_PAIR_SCHEMA_MISMATCH")
            ci <- bi_isa_ids(code$eid, "code")
            di <- bi_isa_ids(dates$eid, "date")
            bi_assert(identical(ci, di), "INCIDENT_AUDIT_CODE_DATE_ID_ORDER_MISMATCH")
            bi_assert(!anyDuplicated(names(code)) && !anyDuplicated(names(dates)) && all(grepl("^i[0-9]+$", names(code)[-1])) && 
                all(grepl("^d[0-9]+$", names(dates)[-1])) && identical(sub("^i", "", names(code)[-1]), sub("^d", 
                "", names(dates)[-1])), "INCIDENT_AUDIT_ARRAY_POSITION_MISMATCH")
            invisible(ci)
        }
        bi_isa_field_summary <- function(cohort, fields) {
            bi_assert("date_attend" %in% names(cohort), "INCIDENT_AUDIT_ENTRY_FIELD_MISSING")
            at <- bi_isa_dates(cohort$date_attend)
            bi_bind(lapply(fields, function(f) {
                present <- f %in% names(cohort)
                tt <- if (present) bi_isa_dates(cohort[[f]]) else rep(as.Date(NA), nrow(cohort))
                data.frame(field = f, present = present, cohort_n = nrow(cohort), nonmissing_count = if (present) sum(!is.na(tt)) else NA_integer_, 
                  on_or_before_entry_count = if (present) sum(tt <= at, na.rm = TRUE) else NA_integer_, after_entry_unfiltered_count = if (present) sum(tt > 
                    at, na.rm = TRUE) else NA_integer_, limitation = "Candidate source only; not PheCode-specific incidence, risk-set events or follow-up coverage", 
                  stringsAsFactors = FALSE)
            }))
        }
        bi_isa_extract <- function(code, dates, mapping) {
            bi_isa_pair(code, dates)
            out <- list()
            for (k in bi_isa_codes()) {
                allowed <- unique(mapping$normalized_code[mapping$trait == k])
                bi_assert(length(allowed) > 0L, "INCIDENT_AUDIT_EMPTY_EXACT_MAP")
                first <- rep(NA_real_, nrow(code))
                observed <- undated <- logical(nrow(code))
                pairs <- 0L
                for (j in seq.int(2L, ncol(code))) {
                  cc <- bi_isa_normalize(code[[j]])
                  hit <- !is.na(cc) & cc %in% allowed
                  if (!any(hit)) next
                  tt <- rep(NA_real_, nrow(code))
                  tt[hit] <- as.numeric(bi_isa_dates(dates[[j]][hit]))
                  observed <- observed | hit
                  undated <- undated | (hit & is.na(tt))
                  pairs <- pairs + sum(hit)
                  first <- pmin(first, tt, na.rm = TRUE)
                  first[!is.finite(first)] <- NA_real_
                }
                out[[k]] <- list(first = as.Date(first, origin = "1970-01-01"), observed = observed, undated = undated, 
                  pairs = pairs)
            }
            out
        }
        bi_isa_candidate_adapter <- function(cohort, code, dates, mapping) {
            ids <- bi_isa_pair(code, dates)
            aid <- bi_isa_ids(cohort$eid, "adapter cohort")
            at <- bi_isa_dates(cohort$date_attend)
            ix <- match(aid, ids)
            events <- bi_isa_extract(code, dates, mapping)
            out <- list()
            for (k in bi_isa_codes()) {
                z <- events[[k]]
                known <- !is.na(ix) & !z$undated[ix]
                known[is.na(known)] <- FALSE
                event <- z$first[ix]
                event[!known] <- as.Date(NA)
                seen <- ifelse(known, as.integer(z$observed[ix]), NA_integer_)
                prior <- ifelse(known & !is.na(at), as.integer(!is.na(event) & event <= at), NA_integer_)
                slug <- gsub(".", "_", k, fixed = TRUE)
                p <- paste0("phe_", slug)
                dd <- data.frame(eid = aid, stringsAsFactors = FALSE)
                dd[[paste0(p, "_first_hospital_date")]] <- event
                dd[[paste0(p, "_hospital_event_observed")]] <- seen
                dd[[paste0(p, "_hospital_prior_only")]] <- prior
                dd[[paste0(p, "_approved_all_source_baseline")]] <- rep(NA_integer_, length(aid))
                out[[k]] <- dd
            }
            out
        }
        bi_isa_config_template <- function() {
            outcomes <- setNames(lapply(bi_isa_codes(), function(k) {
                p <- paste0("phe_", gsub(".", "_", k, fixed = TRUE))
                list(enabled = TRUE, definition_verified = FALSE, event_dates_verified = FALSE, baseline_exclusion_verified = FALSE, 
                  followup_coverage_verified = FALSE, definition_id = paste0("CANDIDATE_EXACT_HOSPITAL_PHECODE_", 
                    k), source_id = "ICD10_PHEWAS_EXACT_DIRECT_PLUS_ROLLUP_AUDIT_ONLY", review_id = "", definition_manifest = "", 
                  event_semantics = "first_ever_recorded", baseline_semantics = "approved_all_available_history", 
                  event_date_field = paste0(p, "_first_hospital_date"), event_observed_field = paste0(p, "_hospital_event_observed"), 
                  baseline_disease_field = paste0(p, "_approved_all_source_baseline"), input_path = "", id_field = "eid", 
                  coverage_start_field = "", coverage_end_field = "", admin_censor_field = "", warning = "Template only: hospital_prior_only must NOT populate approved_all_source_baseline; preserve unverified gates")
            }), bi_isa_codes())
            list(enabled = TRUE, protocol_verified = FALSE, protocol_id = "", time_scale = "attained_age", age_time_scale_verified = FALSE, 
                censoring_verified = FALSE, date_precision_verified = FALSE, date_review_id = "", date_fields = list(birth = "birth_date", 
                  entry = "date_attend", death = "", loss = "date_lost"), birth_date_rule = "", tie_rule = "event_on_censor_counts", 
                sensitivity_s8 = TRUE, outcomes = outcomes, warning = "Raw death source and region-specific completeness must be reviewed; no current FOD maximum date is a cutoff")
        }
        bi_isa_scan <- function(cohort, code, dates, flags, mapping) {
            ids <- bi_isa_pair(code, dates)
            aid <- bi_isa_ids(cohort$eid, "frozen all")
            bi_assert("id" %in% names(flags), "INCIDENT_AUDIT_PHEWAS_ID_FIELD_REQUIRED")
            pid <- bi_isa_ids(flags$id, "PheWAS")
            ix <- match(ids, aid)
            px <- match(ids, pid)
            bi_assert(all(bi_isa_codes() %in% names(flags)), "INCIDENT_AUDIT_PHEWAS_COLUMNS_MISSING")
            at <- bi_isa_dates(cohort$date_attend)[ix]
            rows <- list()
            events <- bi_isa_extract(code, dates, mapping)
            for (k in bi_isa_codes()) {
                mm <- mapping[mapping$trait == k, , drop = FALSE]
                allowed <- unique(mm$normalized_code)
                z <- events[[k]]
                first <- as.numeric(z$first)
                observed <- z$observed
                undated <- z$undated
                pairs <- z$pairs
                old <- flags[[k]][px]
                known <- !is.na(px)
                bi_assert(all(is.na(old) | old %in% c(0, 1)), "INCIDENT_AUDIT_INVALID_PHEWAS_FLAGS")
                rows[[k]] <- data.frame(trait = k, source = "hospital_ICD10_exact_PheWAS_direct_plus_rollup", source_rows = nrow(code), 
                  source_positions = ncol(code) - 1L, cohort_n = nrow(cohort), source_cohort_overlap_count = sum(!is.na(ix)), 
                  cohort_without_source_count = sum(!aid %in% ids), map_rows = nrow(mm), unique_exact_codes = length(allowed), 
                  positive_people_count = sum(observed), positive_code_date_pair_count = pairs, dated_positive_people_count = sum(!is.na(first)), 
                  undated_positive_people_count = sum(undated), on_or_before_entry_count = sum(first <= as.numeric(at), 
                    na.rm = TRUE), after_entry_unfiltered_count = sum(first > as.numeric(at), na.rm = TRUE), phewas_compared_count = sum(known), 
                  phewas_positive_count = sum(old == 1, na.rm = TRUE), phewas_flag_disagreement_count = sum(observed[known] != 
                    (old[known] == 1)), status = "AUDIT_ONLY_NOT_APPROVED", limitation = "Post-baseline raw first-code count only; no baseline all-source exclusion, censoring, risk set or regression", 
                  stringsAsFactors = FALSE)
            }
            bi_bind(rows)
        }
        bi_isa_blockers <- function(m, code) {
            o <- m$outcomes[[code]] %||% list()
            checks <- c(PROTOCOL_NOT_APPROVED = !isTRUE(m$protocol_verified), EXACT_OUTCOME_DEFINITION_NOT_APPROVED = !isTRUE(o$definition_verified), 
                EVENT_DATE_ADAPTER_NOT_APPROVED = !isTRUE(o$event_dates_verified), BASELINE_ALL_SOURCE_EXCLUSION_NOT_APPROVED = !isTRUE(o$baseline_exclusion_verified), 
                SOURCE_REGION_COMPLETENESS_NOT_APPROVED = !isTRUE(o$followup_coverage_verified), CENSORING_NOT_APPROVED = !isTRUE(m$censoring_verified), 
                DATE_PRECISION_NOT_APPROVED = !isTRUE(m$date_precision_verified), ATTAINED_AGE_NOT_APPROVED = !isTRUE(m$age_time_scale_verified), 
                DEFINITION_MANIFEST_MISSING = !bi_temporal_text(o$definition_manifest), EVENT_STATUS_BASELINE_BINDINGS_MISSING = !all(vapply(list(o$event_date_field, 
                  o$event_observed_field, o$baseline_disease_field), bi_temporal_text, logical(1))))
            paste(names(checks)[checks], collapse = ";")
        }
        bi_isa_nofit_matrix <- function(m) {
            z <- expand.grid(exposure = bi_isa_exposures(), trait = bi_isa_codes(), variant = c("primary", "S_le_8"), 
                stringsAsFactors = FALSE)
            z$test_id <- paste(z$exposure, z$trait, z$variant, sep = "__")
            z$status <- "AUDIT_ONLY_NOT_FITTED"
            z$blockers <- vapply(z$trait, function(k) bi_isa_blockers(m, k), character(1))
            z$HR <- z$CI_low <- z$CI_high <- z$p <- z$q <- NA_real_
            z$release_allowed <- FALSE
            z$interpretation <- "Not estimated: neither a negative association nor Cox failure; this audit cannot approve a fit"
            bi_assert(nrow(z) == 48L && !anyDuplicated(z$test_id), "INCIDENT_AUDIT_48_SLOT_MATRIX_FAILURE")
            z
        }
        bi_stage_incident_source_audit <- function(ctx) {
            bi_assert(ctx$stage %in% c("incident_audit", "incident_source_audit") || isTRUE(ctx$synthetic), "INCIDENT_SOURCE_AUDIT_MANUAL_STAGE_REQUIRED")
            bi_assert(isTRUE(ctx$real_authorized) || isTRUE(ctx$synthetic), "REAL_DATA_ACCESS_NOT_AUTHORIZED")
            m <- bi_module_config(ctx, "incident")
            ac <- bi_module_config(ctx, "incident_source_audit")
            invpath <- file.path(ctx$root, "config", "incident_source_inventory.csv")
            bi_assert(file.exists(invpath), paste("INCIDENT_SOURCE_INVENTORY_MISSING", invpath))
            inventory <- read.csv(invpath, colClasses = "character", check.names = FALSE)
            bi_assert(identical(inventory$trait, bi_isa_codes()), "INCIDENT_SOURCE_INVENTORY_SCOPE_MISMATCH")
            allpath <- ctx$cfg$inputs$all_rds
            phepath <- ctx$cfg$inputs$phewas_rds
            rdir <- dirname(phepath)
            paths <- c(all_rds = allpath, phewas_rds = phepath, icd10_code = ac$code_path %||% file.path(rdir, 
                "icd10.code0.rds"), icd10_date = ac$date_path %||% file.path(rdir, "icd10.date0.rds"))
            cat("INCIDENT_SOURCE_AUDIT_ROOT=", ctx$root, "\nALL_RDS=", allpath, "\nPHE_RDATA=", rdir, "\nOUT=", 
                ctx$out, "\nAUDIT_ONLY_NO_MODELS\n", sep = "")
            bi_assert(all(file.exists(paths)), paste("INCIDENT_AUDIT_MISSING_INPUT", paste(paths[!file.exists(paths)], 
                collapse = ",")))
            dest <- normalizePath(ctx$out, mustWork = FALSE)
            for (p in paths) {
                protected <- normalizePath(dirname(p), mustWork = TRUE)
                bi_assert(!identical(dest, protected) && !startsWith(dest, paste0(protected, "/")), "INCIDENT_AUDIT_OUTPUT_OVERLAPS_INPUT")
            }
            before <- lapply(paths, function(p) list(bytes = unname(file.info(p)$size), mtime = unname(as.numeric(file.info(p)$mtime)), 
                sha256 = bi_hash(p, TRUE)))
            expected <- ctx$cfg$input_sha256
            for (k in intersect(names(expected), names(paths))) bi_assert(identical(before[[k]]$sha256, expected[[k]]), 
                paste("INCIDENT_AUDIT_INPUT_SHA256_MISMATCH", k))
            bi_resource_guard(ctx)
            cohort <- as.data.frame(readRDS(allpath))
            bi_isa_ids(cohort$eid, "all")
            fields <- unique(c("birth_date", "date_attend", "date_death", "date_lost", unlist(strsplit(inventory$candidate_fod_fields, 
                ";", fixed = TRUE))))
            field_audit <- bi_isa_field_summary(cohort, fields[nzchar(fields)])
            cohort <- cohort[, c("eid", "date_attend"), drop = FALSE]
            bi_assert(requireNamespace("PheWAS", quietly = TRUE), "INCIDENT_AUDIT_PHEWAS_PACKAGE_REQUIRED")
            maps <- new.env(parent = emptyenv())
            utils::data(list = c("phecode_map", "phecode_rollup_map"), package = "PheWAS", envir = maps)
            direct <- as.data.frame(maps$phecode_map)
            rollup <- as.data.frame(maps$phecode_rollup_map)
            mh <- bi_isa_canonical_hash(direct)
            rh <- bi_isa_canonical_hash(rollup)
            bi_assert(identical(mh, birth_private_setting("PRIVATE_HASH_3")) && identical(rh, birth_private_setting("PRIVATE_HASH_4")), 
                "INCIDENT_AUDIT_PHEWAS_MAP_CHANGED_REVIEW_REQUIRED")
            mapping <- bi_isa_map(direct, rollup)
            pair_audit <- bi_isa_scan(cohort, as.data.frame(readRDS(paths["icd10_code"])), as.data.frame(readRDS(paths["icd10_date"])), 
                as.data.frame(readRDS(phepath)), mapping)
            pair_audit$map_sha256 <- mh
            pair_audit$rollup_sha256 <- rh
            for (k in names(paths)) {
                now <- file.info(paths[[k]])
                bi_assert(identical(unname(now$size), before[[k]]$bytes) && identical(unname(as.numeric(now$mtime)), 
                  before[[k]]$mtime), "INCIDENT_AUDIT_SOURCE_CHANGED_DURING_READ")
            }
            source_audit <- data.frame(source = names(paths), path = unname(paths), bytes = vapply(before, `[[`, 
                numeric(1), "bytes"), sha256 = vapply(before, `[[`, character(1), "sha256"), stringsAsFactors = FALSE)
            bi_save(ctx, "raw/incident_source_audit/source_manifest.csv", source_audit)
            bi_save(ctx, "raw/incident_source_audit/source_inventory.csv", inventory)
            bi_save(ctx, "raw/incident_source_audit/field_availability.csv", field_audit)
            bi_save(ctx, "raw/incident_source_audit/exact_phecode_date_counts.csv", pair_audit)
            bi_save(ctx, "raw/incident_source_audit/nofit_matrix.csv", bi_isa_nofit_matrix(m))
            bi_save(ctx, "raw/incident_source_audit/unapproved_config_template.json", jsonlite::toJSON(bi_isa_config_template(), 
                auto_unbox = TRUE, pretty = TRUE, null = "null"))
            list(status = "AUDIT_ONLY_NOT_FITTED", reason = "48 nofit slots retained; source agreement/counts cannot approve outcome, baseline exclusion or censoring; no risk set or model created")
        }
    }), envir = context)
}

