# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

ev_inc_date <- function(x, label = "date") {
    if (inherits(x, "Date")) 
        return(x)
    if (inherits(x, "POSIXt")) 
        return(as.Date(x))
    z <- as.character(x)
    missing <- is.na(z) | toupper(trimws(z)) %in% c("", "NA", "NAN", "NULL")
    valid <- !missing & grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", z)
    if (any(!missing & !valid)) 
        stop(label, ": non-ISO/unknown date precision")
    ans <- rep(as.Date(NA), length(z))
    ans[valid] <- as.Date(z[valid], format = "%Y-%m-%d")
    if (any(valid & is.na(ans))) 
        stop(label, ": invalid calendar date")
    ans
}

ev_inc_code <- function(x) {
    z <- toupper(gsub("[.[:space:]]", "", as.character(x)))
    z[is.na(x) | z %in% c("", "NA", "NAN", "NULL")] <- NA_character_
    z
}

ev_inc_evidence <- function(x) {
    is.list(x) && isTRUE(x$verified) && length(x$evidence) == 1L && !is.na(x$evidence) && nzchar(trimws(as.character(x$evidence)))
}

ev_inc_contract_example <- function(indir) {
    list(hospital = list(verified = FALSE, evidence = "", position_verified = FALSE, no_record_means_no_hospital_event = FALSE, 
        code_path = file.path(indir, "Rdata", "icd10.code0.rds"), date_path = file.path(indir, "Rdata", "icd10.date0.rds")), 
        coverage = list(verified = FALSE, evidence = "", region_column = "", dates = c(England = NA_character_, 
            Scotland = NA_character_, Wales = NA_character_)), censor = list(verified = FALSE, evidence = "", death_column = "date_death", 
            lost_column = "date_lost", missing_death_means_alive = FALSE, missing_lost_means_not_lost = FALSE, 
            same_day_event = "unverified", age_precision = "unverified"), srd = list(verified = FALSE, evidence = "", 
            baseline_positions = character(), no_record_means_no_report = FALSE, code_path = file.path(indir, "Rdata", 
                "srd.code0.rds")), opcs = list(verified = FALSE, evidence = "", available = TRUE, ihd_codes = character(), 
            position_verified = FALSE, no_record_means_no_procedure = FALSE, code_path = file.path(indir, "Rdata", 
                "opcs4.code0.rds"), date_path = file.path(indir, "Rdata", "opcs4.date0.rds")), other_prior_ihd = list(verified = FALSE, 
            evidence = "", reviewed_sources = character(), status_columns = character(), date_columns = character(), 
            absent_sources = character(), date_missing_means_no_prior_record = FALSE), t2d = list(verified = FALSE, 
            evidence = "", algorithm = "", baseline_status_column = "", baseline_free_values = "no_diabetes", type1_values = "type1", 
            type2_values = "type2", unknown_values = c("uncertain", "unknown"), baseline_only_verified = FALSE, 
            inputs_and_rules_verified = FALSE))
}

ev_inc_specs <- function() {
    data.frame(endpoint = c("T2D", "IHD", "MI", "angina", "atherosclerosis"), phecode_context = c("250.2", "411", 
        "411.2", "411.3", "411.4"), endpoint_name = c("First hospital T2D record among baseline-unidentified diabetes", 
        "First hospital broad IHD record among baseline-unidentified IHD", "First hospital acute MI record among baseline-unidentified IHD", 
        "First hospital angina record among baseline-unidentified IHD", "First hospital coronary atherosclerosis record among baseline-unidentified IHD"), 
        pattern = c("^E11[0-9A-Z]*$", "^I2[0-5][0-9A-Z]*$", "^I2[12][0-9A-Z]*$", "^I20[0-9A-Z]*$", "^I251[0-9A-Z]*$"), 
        definition = c("Hospital ICD-10 E11 (all subcodes); not all-source incident capture", "Hospital ICD-10 I20-I25 (all subcodes); not hard MI", 
            "Hospital ICD-10 I21/I22 (all subcodes); not PheCode 411.2 or legacy date_mi", "Hospital ICD-10 I20 (all subcodes)", 
            "Hospital ICD-10 I25.1 (all subcodes)"), test_family = c("longitudinal_core", "longitudinal_core", 
            rep("longitudinal_secondary", 3)), planned_tests = c(2L, 2L, 3L, 3L, 3L), definition_reference = c("https://doi.org/10.1371/journal.pone.0162388", 
            rep("https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=41270", 4)), stringsAsFactors = FALSE)
}

ev_inc_read_pair <- function(code_path, date_path, eids, positions = NULL) {
    for (p in c(code_path, date_path)) if (length(p) != 1L || !file.exists(p)) 
        stop("Missing validated source: ", paste(p, collapse = ","))
    codes <- as.data.frame(readRDS(code_path), stringsAsFactors = FALSE)
    dates <- as.data.frame(readRDS(date_path), stringsAsFactors = FALSE)
    if (!"eid" %in% names(codes) || !"eid" %in% names(dates)) 
        stop("Source lacks eid")
    if (anyNA(codes$eid) || anyNA(dates$eid) || anyDuplicated(codes$eid) || anyDuplicated(dates$eid)) 
        stop("Source eid missing/duplicate")
    if (!identical(as.character(codes$eid), as.character(dates$eid))) 
        stop("Code/date source eid orders differ; no unchecked positional join")
    cn <- setdiff(names(codes), "eid")
    dn <- setdiff(names(dates), "eid")
    if (!length(cn) || !all(grepl("^i[0-9]+$", cn)) || !all(grepl("^d[0-9]+$", dn)) || !identical(sub("^i", "", 
        cn), sub("^d", "", dn))) 
        stop("Code/date array positions not identical iN/dN")
    if (!is.null(positions)) {
        if (anyDuplicated(positions) || any(!positions %in% cn)) 
            stop("Invalid source positions")
        dn <- paste0("d", sub("^i", "", positions))
        cn <- positions
    }
    idx <- match(as.character(eids), as.character(codes$eid))
    list(codes = codes[idx, cn, drop = FALSE], dates = dates[idx, dn, drop = FALSE], linked = !is.na(idx), positions = cn)
}

ev_inc_scan <- function(pair, pattern, entry) {
    n <- length(entry)
    first <- rep(as.Date(NA), n)
    prior <- rep(FALSE, n)
    undated <- rep(FALSE, n)
    for (j in seq_len(ncol(pair$codes))) {
        code <- ev_inc_code(pair$codes[[j]])
        positive <- !is.na(code) & grepl(pattern, code, perl = TRUE)
        if (!any(positive)) 
            next
        date <- rep(as.Date(NA), n)
        date[positive] <- ev_inc_date(pair$dates[[j]][positive], paste0("source date position ", j))
        undated <- undated | (positive & is.na(date))
        prior <- prior | (positive & !is.na(date) & !is.na(entry) & date <= entry)
        newer <- positive & !is.na(date) & (is.na(first) | date < first)
        first[newer] <- date[newer]
    }
    list(first = first, prior = prior, undated = undated, known = pair$linked & !undated)
}

ev_inc_srd <- function(contract, eids) {
    p <- contract$code_path
    if (length(p) != 1L || !file.exists(p)) 
        stop("Missing baseline SRD code source")
    x <- as.data.frame(readRDS(p), stringsAsFactors = FALSE)
    if (!"eid" %in% names(x) || anyNA(x$eid) || anyDuplicated(x$eid)) 
        stop("Baseline SRD source eid missing/duplicate")
    cols <- contract$baseline_positions
    if (!length(cols) || anyDuplicated(cols) || any(!cols %in% names(x))) 
        stop("Baseline SRD instance-0 positions absent/unverified")
    idx <- match(as.character(eids), as.character(x$eid))
    present <- !is.na(idx)
    ihd <- dm <- unknown <- rep(FALSE, length(eids))
    for (col in cols) {
        code <- as.character(x[[col]][idx])
        vals <- suppressWarnings(as.numeric(code))
        positive <- !is.na(vals)
        ihd <- ihd | (positive & vals %in% c(1074, 1075))
        dm <- dm | (positive & vals %in% c(1220, 1222, 1223))
        unknown <- unknown | (positive & vals %in% c(-1, -3))
        bad <- !is.na(code) & nzchar(trimws(code)) & !toupper(trimws(code)) %in% c("NA", "NAN", "NULL") & is.na(vals)
        if (any(bad)) 
            stop("Unparsed baseline SRD code: source needs audit")
    }
    list(ihd = ihd, dm = dm, known = present & !unknown, unknown = unknown)
}

ev_inc_column <- function(ctx, name) {
    if (length(name) != 1L || is.na(name) || !nzchar(name)) 
        stop("Unconfigured column")
    if (name %in% names(ctx$d)) 
        return(ctx$d[[name]])
    if (!is.null(ctx$raw_all) && name %in% names(ctx$raw_all)) {
        raw <- ctx$raw_all
        if (anyNA(raw$eid) || anyDuplicated(raw$eid)) 
            stop("raw_all eid not unique")
        return(raw[[name]][match(as.character(ctx$d$eid), as.character(raw$eid))])
    }
    stop("Missing audited column: ", name)
}

ev_inc_censor <- function(ctx, cfg, entry) {
    region <- as.character(ev_inc_column(ctx, cfg$coverage$region_column))
    cutoffs <- cfg$coverage$dates
    if (is.null(names(cutoffs)) || anyDuplicated(names(cutoffs))) 
        stop("Coverage dates must be region-named")
    coverage <- ev_inc_date(unname(cutoffs[match(region, names(cutoffs))]), "coverage")
    admin <- ev_inc_date(ctx$config$admin_follow_end, "administrative cutoff")
    if (length(admin) != 1L || is.na(admin) || admin > as.Date("2023-04-01")) 
        stop("Administrative cutoff absent or later than locked 2023-04-01")
    death <- ev_inc_date(ev_inc_column(ctx, cfg$censor$death_column), "death")
    lost <- ev_inc_date(ev_inc_column(ctx, cfg$censor$lost_column), "lost")
    stopifnot(length(death) == length(entry), length(lost) == length(entry))
    end <- pmin(rep(admin, length(entry)), coverage, death, lost, na.rm = TRUE)
    end[is.na(coverage)] <- as.Date(NA)
    list(end = end, death = death, lost = lost, coverage = coverage, known = !is.na(coverage), version = paste("admin<=2023-04-01;audited uniform conservative hospital horizon;", 
        "min(event,death,lost,coverage,admin);same-day", cfg$censor$same_day_event), missing_region_n = sum(is.na(coverage)))
}

ev_inc_add_years <- function(x, years = 2L) {
    ans <- rep(as.Date(NA), length(x))
    valid <- !is.na(x)
    yr <- as.integer(format(x[valid], "%Y")) + as.integer(years)
    md <- format(x[valid], "%m-%d")
    leap <- (yr%%4L == 0L & yr%%100L != 0L) | yr%%400L == 0L
    md[md == "02-29" & !leap] <- "02-28"
    ans[valid] <- as.Date(paste0(yr, "-", md), format = "%Y-%m-%d")
    ans
}

ev_inc_risk <- function(d, entry, event_date, censor, prior, known, landmark = FALSE, same_day = "event") {
    stopifnot(nrow(d) == length(entry), length(entry) == length(event_date), length(entry) == length(censor$end), 
        length(entry) == length(prior), length(entry) == length(known))
    if (anyNA(prior) || anyNA(known)) 
        stop("Prior/known masks cannot silently contain NA")
    if (!"S" %in% names(d)) 
        stop("Missing locked full-sibship size")
    first_last <- !is.na(d$X) & d$X %in% c(0, 1)
    supported <- sort(unique(d$S[first_last & !is.na(d$S) & d$S >= 2 & d$S == floor(d$S)]))
    supported <- supported[vapply(supported, function(s) all(c(0, 1) %in% d$X[first_last & !is.na(d$S) & d$S == 
        s]), logical(1))]
    mask <- first_last
    flow <- list()
    add <- function(step, keep) {
        before <- sum(mask)
        mask <<- mask & keep
        flow[[length(flow) + 1L]] <<- data.frame(step = step, n_before = before, n_excluded = before - sum(mask), 
            n_remaining = sum(mask), stringsAsFactors = FALSE)
    }
    add("outcome_independent_supported_exact_sibship_size", !is.na(d$S) & d$S %in% supported)
    add("baseline_source_observable", known)
    add("baseline_without_identified_disease", !prior)
    add("baseline_date_present", !is.na(entry))
    add("valid_age_at_entry", !is.na(d$age_entry) & is.finite(d$age_entry) & d$age_entry > 0)
    add("regional_followup_coverage_known", censor$known)
    add("positive_followup", !is.na(censor$end) & !is.na(entry) & censor$end > entry)
    start <- entry
    if (landmark) {
        start <- ev_inc_add_years(entry, 2L)
        alive <- is.na(censor$death) | censor$death > start
        observable <- !is.na(censor$end) & censor$end > start
        event_free <- is.na(event_date) | event_date > start
        add("alive_at_2calendar_year_landmark", alive & !is.na(start))
        add("observable_beyond_landmark", observable)
        add("event_free_through_landmark_day", event_free)
    }
    event <- !is.na(event_date) & !is.na(start) & !is.na(censor$end) & event_date > start & event_date <= censor$end
    if (same_day == "death") 
        event <- event & (is.na(censor$death) | event_date < censor$death)
    if (!same_day %in% c("event", "death")) 
        stop("Unverified event/death same-day rule")
    exit <- censor$end
    exit[event] <- event_date[event]
    out <- d[mask, , drop = FALSE]
    out$event <- as.integer(event[mask])
    out$followup_years <- as.numeric(exit[mask] - start[mask])/365.25
    out$age_entry <- d$age_entry[mask] + as.numeric(start[mask] - entry[mask])/365.25
    out$age_exit <- d$age_entry[mask] + as.numeric(exit[mask] - entry[mask])/365.25
    if (anyDuplicated(out$eid) || anyNA(out$eid) || any(out$followup_years <= 0) || any(out$age_exit <= out$age_entry)) 
        stop("Invalid Cox interval/eid after risk-set construction")
    list(data = out, flow = do.call(rbind, flow), n_eligible = sum(mask), events_after_death = sum(mask & !is.na(event_date) & 
        !is.na(censor$death) & event_date > censor$death & event_date > entry, na.rm = TRUE), events_after_death_all = sum(!is.na(event_date) & 
        !is.na(censor$death) & event_date > censor$death & event_date > entry, na.rm = TRUE), event_on_death = sum(mask & 
        !is.na(event_date) & !is.na(censor$death) & event_date == censor$death & event_date > entry, na.rm = TRUE))
}

ev_inc_terms <- function(ctx, model) {
    if (model == "time_on_study") {
        if (is.null(ctx$primary_terms) || !length(ctx$primary_terms)) 
            stop("Missing frozen primary terms")
        return(paste(unique(c("X", ctx$primary_terms)), collapse = " + "))
    }
    if (length(ctx$birth_term) != 1L || !nzchar(ctx$birth_term)) 
        stop("Missing frozen birth spline")
    paste(c("X", "factor(S)", "factor(sex)", "factor(ethnicity)", "factor(center)", ctx$birth_term), collapse = " + ")
}

ev_inc_fit <- function(ctx, dat, spec, model, censor_version) {
    fit_file <- file.path(ctx$out, "cache", paste0("incident_fit_", spec$endpoint, "_", model, ".rds"))
    if (file.exists(fit_file)) 
        stop("Protected existing fit without a matching task checkpoint: ", fit_file)
    response <- if (model == "time_on_study") 
        "survival::Surv(followup_years, event)"
    else "survival::Surv(age_entry, age_exit, event)"
    formula <- stats::as.formula(paste(response, "~", ev_inc_terms(ctx, model)))
    environment(formula) <- globalenv()
    mm <- stats::model.matrix(formula, dat)
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    if (!ncol(mm) || qr(mm)$rank < ncol(mm)) 
        stop("Cox design rank deficient; no silent term removal")
    cond <- tryCatch(kappa(mm, exact = FALSE), error = function(e) NA_real_)
    if (length(unique(dat$X)) != 2L || sum(dat$event) < 20L || sum(dat$event == 0L) < 20L) 
        stop("Prespecified estimability rule fails: both birth groups, >=20 events and non-events required")
    warnings <- character()
    fit <- withCallingHandlers(survival::coxph(formula, data = dat, ties = "efron", x = TRUE, y = TRUE, singular.ok = FALSE, 
        na.action = stats::na.fail), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    tab <- summary(fit)$coefficients
    if (!"X" %in% rownames(tab) || any(!is.finite(stats::coef(fit))) || any(!is.finite(sqrt(diag(stats::vcov(fit))))) || 
        any(grepl("infinite|converg|iteration|singular", warnings, ignore.case = TRUE))) 
        stop("Cox unstable/non-convergent/infinite estimate: ", paste(warnings, collapse = "; "))
    if (length(fit$n) != 1L || fit$n != nrow(dat) || fit$n != length(unique(dat$eid)) || fit$nevent != sum(dat$event) || 
        nrow(fit$x) != nrow(dat)) 
        stop("Cox actual N/events/eid invariant failed")
    beta <- tab["X", "coef"]
    se <- tab["X", "se(coef)"]
    p <- tab["X", "Pr(>|z|)"]
    if (!all(is.finite(c(beta, se, p))) || se <= 0) 
        stop("X coefficient/SE/P invalid")
    ph <- tryCatch(survival::cox.zph(fit, terms = TRUE, global = TRUE), error = identity)
    ph_status <- "diagnostic_failed"
    ph_x <- ph_global <- NA_real_
    ph_table <- NULL
    if (!inherits(ph, "error")) {
        ph_table <- data.frame(endpoint = spec$endpoint, model = model, term = rownames(ph$table), ph$table, row.names = NULL, 
            check.names = FALSE)
        if ("X" %in% rownames(ph$table)) 
            ph_x <- ph$table["X", "p"]
        if ("GLOBAL" %in% rownames(ph$table)) 
            ph_global <- ph$table["GLOBAL", "p"]
        ph_status <- if (is.finite(ph_x) && is.finite(ph_global)) 
            "diagnostics_available_review_plot"
        else "diagnostic_incomplete"
        ph_file <- file.path(ctx$out, "figures", paste0("Schoenfeld_", spec$endpoint, "_", model, ".pdf"))
        dir.create(dirname(ph_file), recursive = TRUE, showWarnings = FALSE)
        if (file.exists(ph_file)) 
            stop("Protected existing PH figure: ", ph_file)
        tmp <- paste0(ph_file, ".tmp.", Sys.getpid())
        grDevices::pdf(tmp, width = 6.4, height = 5)
        tryCatch({
            if ("X" %in% colnames(ph$y)) 
                graphics::plot(ph, var = match("X", colnames(ph$y)), main = paste(spec$endpoint, model, "X PH diagnostic"))
            else graphics::plot.new()
        }, finally = grDevices::dev.off())
        if (!file.rename(tmp, ph_file)) 
            stop("Cannot install PH figure: ", ph_file)
    }
    else warnings <- c(warnings, paste("cox.zph:", conditionMessage(ph)))
    stats <- data.frame(endpoint = spec$endpoint, endpoint_name = spec$endpoint_name, outcome = spec$phecode_context, 
        metric = "HR", exposure = "X", scope = "hospital_first_record", model = model, primary_flag = model == 
            "primary_age", time_scale = if (model == "time_on_study") 
            "time_on_study"
        else "attained_age", landmark_years = if (model == "landmark_2y") 
            2L
        else 0L, status = "ok", technical_status = if (ph_status == "diagnostic_failed") 
            "needs_review"
        else "fit_complete", reason = "", formula = paste(deparse(formula, width.cutoff = 500), collapse = " "), 
        N = as.integer(fit$n), model_rows = nrow(fit$x), unique_eid = length(unique(dat$eid)), sample_hash = digest::digest(sort(as.character(dat$eid)), 
            algo = "sha256"), events = as.integer(fit$nevent), N_first = sum(dat$X == 0), N_last = sum(dat$X == 
            1), events_first = sum(dat$X == 0 & dat$event == 1), events_last = sum(dat$X == 1 & dat$event == 1), 
        person_years_first = sum(dat$followup_years[dat$X == 0]), person_years_last = sum(dat$followup_years[dat$X == 
            1]), followup_median = stats::median(dat$followup_years), followup_q25 = unname(stats::quantile(dat$followup_years, 
            0.25)), followup_q75 = unname(stats::quantile(dat$followup_years, 0.75)), beta = beta, SE = se, HR = exp(beta), 
        CI_low = exp(beta - 1.96 * se), CI_high = exp(beta + 1.96 * se), P = p, p_adjusted = NA_real_, test_family = spec$test_family, 
        planned_tests = spec$planned_tests, correction = if (spec$test_family == "longitudinal_core") 
            "Holm n=2"
        else "BH n=3", ph_X_P = ph_x, ph_global_P = ph_global, ph_status = ph_status, condition_number = cond, 
        censor_version = censor_version, warnings = paste(unique(warnings), collapse = "; "), inference_scope = "Post-discovery internal robustness; cause-specific hospital-record HR, not cumulative risk", 
        stringsAsFactors = FALSE)
    ev_save(fit, fit_file)
    list(table = stats, ph = ph_table)
}

ev_inc_blank <- function(spec, model, status, reason) {
    data.frame(endpoint = spec$endpoint, endpoint_name = spec$endpoint_name, model = model, outcome = spec$phecode_context, 
        metric = "HR", exposure = "X", scope = "hospital_first_record", primary_flag = model == "primary_age", 
        time_scale = if (model == "time_on_study") 
            "time_on_study"
        else "attained_age", landmark_years = if (model == "landmark_2y") 
            2L
        else 0L, status = status, technical_status = if (status == "blocked") 
            "data_insufficient"
        else "needs_repair", reason = reason, N = NA_integer_, model_rows = NA_integer_, unique_eid = NA_integer_, 
        events = NA_integer_, N_first = NA_integer_, N_last = NA_integer_, events_first = NA_integer_, events_last = NA_integer_, 
        person_years_first = NA_real_, person_years_last = NA_real_, followup_median = NA_real_, followup_q25 = NA_real_, 
        followup_q75 = NA_real_, beta = NA_real_, SE = NA_real_, HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, 
        P = NA_real_, p_adjusted = NA_real_, test_family = spec$test_family, planned_tests = spec$planned_tests, 
        correction = if (spec$test_family == "longitudinal_core") 
            "Holm n=2"
        else "BH n=3", ph_X_P = NA_real_, ph_global_P = NA_real_, ph_status = "not_evaluated", condition_number = NA_real_, 
        censor_version = "unverified", warnings = "", inference_scope = "Not evaluated; not evidence of no association", 
        stringsAsFactors = FALSE)
}

ev_inc_gate <- function(cfg, endpoint) {
    reasons <- character()
    need <- function(ok, msg) if (!isTRUE(ok)) 
        reasons <<- c(reasons, msg)
    h <- cfg$hospital
    c <- cfg$coverage
    z <- cfg$censor
    need(ev_inc_evidence(h) && isTRUE(h$position_verified) && isTRUE(h$no_record_means_no_hospital_event), "Audit hospital source/version, iN/dN alignment and no-record semantics")
    need(ev_inc_evidence(c) && length(c$region_column) == 1L && nzchar(c$region_column) && !is.null(names(c$dates)), 
        "Supply verified local hospital completeness cutoff per source/region, with region column")
    need(ev_inc_evidence(z) && isTRUE(z$missing_death_means_alive) && isTRUE(z$missing_lost_means_not_lost) && 
        z$same_day_event %in% c("event", "death") && length(z$age_precision) == 1L && !z$age_precision %in% c("", 
        "unverified"), "Verify uncapped death/lost provenance, absent-record semantics, same-day rule and age precision")
    if (endpoint == "T2D") {
        t <- cfg$t2d
        need(ev_inc_evidence(t) && isTRUE(t$baseline_only_verified) && isTRUE(t$inputs_and_rules_verified) && length(t$algorithm) == 
            1L && nzchar(t$algorithm) && length(t$baseline_status_column) == 1L && nzchar(t$baseline_status_column) && 
            length(t$baseline_free_values) > 0L && length(t$type1_values) > 0L && length(t$type2_values) > 0L && 
            length(t$unknown_values) > 0L, "Verify Eastwood appendix or explicit local baseline DM algorithm/inputs; distinguish type1/type2/uncertain; do not use dm.yes")
    }
    else {
        s <- cfg$srd
        o <- cfg$opcs
        p <- cfg$other_prior_ihd
        need(ev_inc_evidence(s) && length(s$baseline_positions) > 0L && isTRUE(s$no_record_means_no_report), "Audit instance-0 SRD positions and missing semantics for prior 1074/1075")
        need(ev_inc_evidence(o) && identical(o$available, FALSE) || ev_inc_evidence(o) && isTRUE(o$available) && 
            isTRUE(o$position_verified) && length(o$ihd_codes) > 0L, "Review OPCS prior coronary disease/revascularisation codes, or document genuinely unavailable source")
        need(!isTRUE(o$available) || isTRUE(o$no_record_means_no_procedure), "Verify absent OPCS array entries mean no recorded procedure, not unobservable data")
        need(ev_inc_evidence(p) && length(p$reviewed_sources) > 0L, "Review ICD9/GP/UKB algorithm/self-report prior-IHD evidence; specify verified flags/dates or unavailable sources")
        need(!length(p$date_columns) || isTRUE(p$date_missing_means_no_prior_record), "Verify missing derived prior-IHD dates mean no recorded disease; otherwise use 0/1/NA source status flags")
    }
    unique(reasons)
}

ev_inc_endpoint <- function(ctx, cfg, spec, hospital, censor, srd = NULL, opcs = NULL) {
    d <- as.data.frame(ctx$d)
    entry <- ev_inc_date(ev_inc_column(ctx, "date_attend"), "baseline measurement")
    event <- ev_inc_scan(hospital, spec$pattern, entry)
    prior <- event$prior
    known <- event$known
    baseline_counts <- list(hospital_endpoint = event$prior, hospital_undated = event$undated)
    if (spec$endpoint == "T2D") {
        dm <- ev_inc_scan(hospital, "^E1[0-4][0-9A-Z]*$", entry)
        t <- cfg$t2d
        status <- as.character(ev_inc_column(ctx, t$baseline_status_column))
        free_values <- as.character(t$baseline_free_values)
        excluded_values <- unique(c(as.character(t$type1_values), as.character(t$type2_values), as.character(t$unknown_values)))
        if (!length(free_values) || anyNA(free_values) || any(free_values %in% excluded_values)) 
            stop("Baseline DM algorithm status sets overlap/empty")
        free <- !is.na(status) & status %in% free_values
        diagnosed <- !is.na(status) & status %in% c(t$type1_values, t$type2_values)
        prior <- prior | dm$prior | diagnosed
        known <- known & dm$known & (free | diagnosed)
        baseline_counts$hospital_any_dm <- dm$prior
        baseline_counts$algorithm_not_verified_diabetes_free <- !free
        baseline_counts$algorithm_type1 <- !is.na(status) & status %in% t$type1_values
        baseline_counts$algorithm_type2 <- !is.na(status) & status %in% t$type2_values
        baseline_counts$algorithm_uncertain <- is.na(status) | status %in% t$unknown_values
    }
    else {
        ihd <- ev_inc_scan(hospital, "^I2[0-5][0-9A-Z]*$", entry)
        prior <- prior | ihd$prior | srd$ihd
        known <- known & ihd$known & srd$known
        baseline_counts$hospital_any_ihd <- ihd$prior
        baseline_counts$baseline_srd_ihd <- srd$ihd
        baseline_counts$baseline_srd_unknown <- !srd$known
        if (!is.null(opcs)) {
            codes <- ev_inc_code(cfg$opcs$ihd_codes)
            if (anyNA(codes) || any(!grepl("^[A-Z][0-9A-Z]{1,6}$", codes))) 
                stop("Prior OPCS mapping requires exact codes, not unaudited regular expressions")
            regex <- paste0("^(?:", paste(codes, collapse = "|"), ")$")
            surgery <- ev_inc_scan(opcs, regex, entry)
            prior <- prior | surgery$prior
            known <- known & surgery$known
            baseline_counts$prior_coronary_procedure <- surgery$prior
            baseline_counts$undated_coronary_procedure <- surgery$undated
        }
        for (col in cfg$other_prior_ihd$status_columns) {
            v <- ev_inc_column(ctx, col)
            if (any(!is.na(v) & !v %in% c(0, 1, FALSE, TRUE))) 
                stop("Prior-IHD flag not 0/1/NA: ", col)
            prior <- prior | (!is.na(v) & v == 1)
            known <- known & !is.na(v)
            baseline_counts[[paste0("prior_flag_", col)]] <- !is.na(v) & v == 1
        }
        for (col in cfg$other_prior_ihd$date_columns) {
            date <- ev_inc_date(ev_inc_column(ctx, col), col)
            positive <- !is.na(date) & !is.na(entry) & date <= entry
            prior <- prior | positive
            baseline_counts[[paste0("prior_date_", col)]] <- positive
        }
    }
    risk <- ev_inc_risk(d, entry, event$first, censor, prior, known, same_day = cfg$censor$same_day_event)
    if (risk$events_after_death > 0L) 
        stop("Endpoint dates after verified death: ", risk$events_after_death, "; source reconciliation required")
    union_formula <- stats::as.formula(paste("~", ev_inc_terms(ctx, "primary_age"), "+", ev_inc_terms(ctx, "time_on_study")))
    base <- risk$data
    mf <- stats::model.frame(union_formula, base, na.action = stats::na.pass)
    complete <- stats::complete.cases(mf)
    for (v in names(mf)) if (is.numeric(mf[[v]])) {
        z <- mf[[v]]
        finite <- if (is.matrix(z)) 
            rowSums(!is.finite(z)) == 0L
        else is.finite(z)
        complete <- complete & finite
    }
    if (!is.logical(complete) || !is.null(dim(complete)) || length(complete) != nrow(base)) 
        stop("Cox complete-case mask must be one logical value per participant")
    extra <- data.frame(step = "main_time_covariate_union_complete", n_before = nrow(base), n_excluded = sum(!complete), 
        n_remaining = sum(complete), stringsAsFactors = FALSE)
    base <- droplevels(base[complete, , drop = FALSE])
    flow <- rbind(risk$flow, extra)
    flow$endpoint <- spec$endpoint
    flow$model <- "primary_age_and_time_on_study"
    baseline <- data.frame(step = paste0("baseline_source_count_nonmutuallyexclusive:", names(baseline_counts)), 
        n_before = sum(!is.na(d$X)), n_excluded = vapply(baseline_counts, function(z) sum(z & !is.na(d$X), na.rm = TRUE), 
            integer(1)), n_remaining = NA_integer_, endpoint = spec$endpoint, model = "baseline_audit", stringsAsFactors = FALSE)
    rows <- list()
    ph <- list()
    for (model in c("primary_age", "time_on_study", "landmark_2y")) {
        model_data <- base
        if (model == "landmark_2y") {
            lrisk <- ev_inc_risk(d, entry, event$first, censor, prior, known, landmark = TRUE, same_day = cfg$censor$same_day_event)
            model_data <- lrisk$data[as.character(lrisk$data$eid) %in% as.character(base$eid), , drop = FALSE]
            model_data <- droplevels(model_data)
            lflow <- lrisk$flow
            lflow$endpoint <- spec$endpoint
            lflow$model <- model
            flow <- rbind(flow, lflow, data.frame(step = "baseline_union_complete_landmark_intersection", n_before = nrow(lrisk$data), 
                n_excluded = nrow(lrisk$data) - nrow(model_data), n_remaining = nrow(model_data), endpoint = spec$endpoint, 
                model = model))
        }
        result <- tryCatch(ev_run_task(ctx, paste0("incident_", spec$endpoint, "_", model), function() ev_inc_fit(ctx, 
            model_data, spec, model, censor$version)), error = function(e) {
            if (inherits(e, "ev_budget_end")) 
                stop(e)
            e
        })
        if (inherits(result, "error")) {
            rows[[model]] <- ev_inc_blank(spec, model, if (inherits(result, "ev_budget_end")) 
                "deferred"
            else "failed", conditionMessage(result))
        }
        else if (is.null(result) || !is.list(result) || is.null(result$table)) {
            rows[[model]] <- ev_inc_blank(spec, model, "deferred", "Task budget/checkpoint returned no model result")
        }
        else {
            rows[[model]] <- result$table
            ph[[model]] <- result$ph
        }
    }
    list(table = ev_bind(rows), flow = rbind(flow, baseline), ph = ev_bind(ph), date_audit = data.frame(endpoint = spec$endpoint, 
        events_after_death = risk$events_after_death, events_after_death_all = risk$events_after_death_all, event_on_death = risk$event_on_death, 
        missing_region_n = censor$missing_region_n))
}

ev_incident <- function(ctx, audit_only = FALSE) {
    existing_table <- file.path(ctx$out, "tables", "incident_core.csv")
    if (audit_only && file.exists(existing_table)) {
        previous <- utils::read.csv(existing_table, stringsAsFactors = FALSE, check.names = FALSE)
        if (("HR" %in% names(previous) && any(is.finite(suppressWarnings(as.numeric(previous$HR))))) || ("status" %in% 
            names(previous) && any(previous$status == "ok", na.rm = TRUE))) 
            stop("Audit-only cannot replace completed incident results; use a new output directory")
    }
    specs <- ev_inc_specs()
    cfg <- ctx$config$incident
    if (is.null(cfg)) 
        cfg <- ev_inc_contract_example(ctx$indir)
    manifests <- tables <- flows <- ph <- date_audit <- list()
    shared <- NULL
    shared_error <- NULL
    for (i in seq_len(nrow(specs))) {
        spec <- specs[i, , drop = FALSE]
        reasons <- ev_inc_gate(cfg, spec$endpoint)
        models <- c("primary_age", "time_on_study", "landmark_2y")
        gate_status <- if (length(reasons)) 
            "blocked"
        else "audited_contract_ready"
        if (!length(reasons) && audit_only) {
            gate_status <- "audited_contract_not_runtime_checked_not_run"
            tables[[spec$endpoint]] <- ev_bind(lapply(models, function(model) ev_inc_blank(spec, model, "not_run", 
                "Audit-only invocation; no source arrays read and no models fitted")))
            flows[[spec$endpoint]] <- data.frame(endpoint = spec$endpoint, model = "all", step = "audit_only_no_risk_set_constructed", 
                n_before = sum(!is.na(ctx$d$X)), n_excluded = NA_integer_, n_remaining = NA_integer_, reason = "No model requested")
        }
        if (!length(reasons) && !audit_only) {
            if (is.null(shared) && is.null(shared_error)) {
                shared <- tryCatch({
                  entry <- ev_inc_date(ev_inc_column(ctx, "date_attend"), "baseline measurement")
                  list(hospital = ev_inc_read_pair(cfg$hospital$code_path, cfg$hospital$date_path, ctx$d$eid), 
                    censor = ev_inc_censor(ctx, cfg, entry))
                }, error = function(e) {
                  if (inherits(e, "ev_budget_end")) 
                    stop(e)
                  e
                })
                if (inherits(shared, "error")) {
                  shared_error <- conditionMessage(shared)
                  shared <- NULL
                }
            }
            if (!is.null(shared_error)) 
                reasons <- shared_error
            if (!length(reasons)) {
                result <- tryCatch({
                  srd <- opcs <- NULL
                  if (spec$endpoint != "T2D") {
                    if (is.null(shared$srd)) 
                      shared$srd <- ev_inc_srd(cfg$srd, ctx$d$eid)
                    srd <- shared$srd
                    if (isTRUE(cfg$opcs$available)) {
                      if (is.null(shared$opcs)) 
                        shared$opcs <- ev_inc_read_pair(cfg$opcs$code_path, cfg$opcs$date_path, ctx$d$eid)
                      opcs <- shared$opcs
                    }
                  }
                  ev_inc_endpoint(ctx, cfg, spec, shared$hospital, shared$censor, srd, opcs)
                }, error = function(e) {
                  if (inherits(e, "ev_budget_end")) 
                    stop(e)
                  e
                })
                if (inherits(result, "error")) 
                  reasons <- conditionMessage(result)
                else {
                  tables[[spec$endpoint]] <- result$table
                  flows[[spec$endpoint]] <- result$flow
                  ph[[spec$endpoint]] <- result$ph
                  date_audit[[spec$endpoint]] <- result$date_audit
                }
            }
        }
        if (length(reasons)) {
            gate_status <- "blocked"
            tables[[spec$endpoint]] <- ev_bind(lapply(models, function(model) ev_inc_blank(spec, model, "blocked", 
                paste(reasons, collapse = "; "))))
            flows[[spec$endpoint]] <- data.frame(endpoint = spec$endpoint, model = "all", step = "blocked_before_risk_set", 
                n_before = sum(!is.na(ctx$d$X)), n_excluded = NA_integer_, n_remaining = NA_integer_, reason = paste(reasons, 
                  collapse = "; "))
        }
        manifests[[spec$endpoint]] <- data.frame(endpoint = spec$endpoint, endpoint_name = spec$endpoint_name, 
            definition = spec$definition, normalized_code_regex = spec$pattern, source = "Hospital ICD-10 summary diagnosis/date arrays", 
            source_code_path = if (length(cfg$hospital$code_path)) 
                cfg$hospital$code_path
            else "unconfigured", source_date_path = if (length(cfg$hospital$date_path)) 
                cfg$hospital$date_path
            else "unconfigured", source_evidence = paste(cfg$hospital$evidence, collapse = "; "), srd_baseline_positions = paste(cfg$srd$baseline_positions, 
                collapse = "; "), opcs_prior_ihd_exact_codes = paste(cfg$opcs$ihd_codes, collapse = "; "), other_prior_ihd_reviewed_sources = paste(cfg$other_prior_ihd$reviewed_sources, 
                collapse = "; "), other_prior_ihd_absent_sources = paste(cfg$other_prior_ihd$absent_sources, collapse = "; "), 
            other_prior_ihd_status_columns = paste(cfg$other_prior_ihd$status_columns, collapse = "; "), other_prior_ihd_date_columns = paste(cfg$other_prior_ihd$date_columns, 
                collapse = "; "), definition_reference = spec$definition_reference, baseline_exclusion = if (spec$endpoint == 
                "T2D") 
                "Verified baseline DM algorithm: exclude type1/type2/uncertain/unknown; prior E10-E14 and undated DM excluded"
            else "Any identified baseline IHD: hospital I20-I25 + baseline SRD 1074/1075 + verified procedures/other available sources; undated unknown excluded", 
            baseline_algorithm = if (spec$endpoint == "T2D") 
                paste(cfg$t2d$algorithm, collapse = "; ")
            else "Audited multi-source baseline IHD exclusion; not a clinical adjudication algorithm", baseline_evidence = if (spec$endpoint == 
                "T2D") 
                paste(cfg$t2d$evidence, collapse = "; ")
            else paste(c(cfg$srd$evidence, cfg$opcs$evidence, cfg$other_prior_ihd$evidence), collapse = "; "), 
            time_origin = "Verified baseline exposure measurement date; future-response backfill prohibited", age_method = paste(cfg$censor$age_precision, 
                collapse = "; "), end_rule = "First event/death/lost/verified source-region coverage/administrative cutoff, whichever earliest", 
            administrative_cap = "2023-04-01; earlier verified local source coverage takes precedence", coverage_evidence = paste(cfg$coverage$evidence, 
                collapse = "; "), censor_evidence = paste(cfg$censor$evidence, collapse = "; "), same_day_rule = paste(cfg$censor$same_day_event, 
                collapse = "; "), planned_tests = spec$planned_tests, test_family = spec$test_family, status = gate_status, 
            reason = paste(reasons, collapse = "; "), interpretation = "Post-discovery internal robustness; first hospital record, not complete incidence or causal inference", 
            stringsAsFactors = FALSE)
    }
    tab <- ev_bind(tables)
    for (model in c("primary_age", "time_on_study", "landmark_2y")) for (family in c("longitudinal_core", "longitudinal_secondary")) {
        ix <- which(tab$model == model & tab$test_family == family & is.finite(tab$P))
        if (length(ix)) 
            tab$p_adjusted[ix] <- stats::p.adjust(tab$P[ix], method = if (family == "longitudinal_core") 
                "holm"
            else "BH", n = if (family == "longitudinal_core") 
                2L
            else 3L)
    }
    manifest <- ev_bind(manifests)
    flow <- ev_bind(flows)
    ev_write(manifest, file.path(ctx$out, "tables", "outcome_manifest.csv"))
    ev_write(tab, file.path(ctx$out, "tables", "incident_core.csv"))
    ev_write(flow, file.path(ctx$out, "tables", "incident_flow.csv"))
    ph_out <- ev_bind(ph)
    if (!nrow(ph_out)) 
        ph_out <- data.frame(endpoint = character(), model = character(), term = character(), chisq = numeric(), 
            df = numeric(), p = numeric())
    dates_out <- ev_bind(date_audit)
    if (!nrow(dates_out)) 
        dates_out <- data.frame(endpoint = character(), events_after_death = integer(), events_after_death_all = integer(), 
            event_on_death = integer(), missing_region_n = integer())
    ev_write(ph_out, file.path(ctx$out, "tables", "incident_ph.csv"))
    ev_write(dates_out, file.path(ctx$out, "audit", "incident_date_audit.csv"))
    ev_record_status(ctx, tab)
    invisible(list(table = tab, manifest = manifest, flow = flow))
}

ev_incident_audit <- function(ctx) ev_incident(ctx, audit_only = TRUE)

ev_incident_tests <- function() {
    d <- data.frame(eid = as.character(1:6), X = c(0, 1, 0, 1, 0, 1), age_entry = rep(50, 6), S = rep(2L, 6))
    entry <- as.Date(rep("2020-01-01", 6))
    event <- as.Date(c("2019-01-01", "2021-01-01", "2022-01-01", "2023-01-01", NA, "2025-01-01"))
    censor <- list(end = as.Date(rep("2024-01-01", 6)), death = rep(as.Date(NA), 6), known = rep(TRUE, 6))
    r <- ev_inc_risk(d, entry, event, censor, prior = c(TRUE, rep(FALSE, 5)), known = rep(TRUE, 6))
    stopifnot(!"1" %in% r$data$eid, nrow(r$data) == 5L, sum(r$data$event) == 3L, r$data$event[r$data$eid == "6"] == 
        0L, length(unique(r$data$eid)) == nrow(r$data))
    l <- ev_inc_risk(d, entry, event, censor, prior = c(TRUE, rep(FALSE, 5)), known = rep(TRUE, 6), landmark = TRUE)
    stopifnot(identical(l$data$eid, c("4", "5", "6")), nrow(l$data) == 3L, sum(l$data$event) == 1L, isTRUE(all.equal(l$data$followup_years[l$data$eid == 
        "5"], 730/365.25)), isTRUE(all.equal(l$data$age_entry[1], 50 + 731/365.25)), ev_inc_add_years(as.Date("2020-02-29")) == 
        as.Date("2022-02-28"))
    unknown <- ev_inc_risk(d, entry, event, censor, prior = c(TRUE, rep(FALSE, 5)), known = c(TRUE, TRUE, TRUE, 
        TRUE, FALSE, TRUE))
    stopifnot(!"5" %in% unknown$data$eid)
    toy <- data.frame(eid = as.character(seq_len(120)), X = rep(0:1, 60), age_entry = rep(40, 120), age_exit = 40 + 
        (seq_len(120)%%17 + 1)/10, event = as.integer(seq_len(120)%%3 == 0))
    fit <- survival::coxph(survival::Surv(age_entry, age_exit, event) ~ X, data = toy, x = TRUE, y = TRUE, na.action = stats::na.fail)
    stopifnot(fit$n == nrow(toy), fit$n == length(unique(toy$eid)), fit$nevent == sum(toy$event), nrow(fit$x) == 
        nrow(toy))
    invisible(TRUE)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        ev_inc_date <- function(x, label = "date") {
            if (inherits(x, "Date")) return(x)
            if (inherits(x, "POSIXt")) return(as.Date(x))
            z <- as.character(x)
            missing <- is.na(z) | toupper(trimws(z)) %in% c("", "NA", "NAN", "NULL")
            valid <- !missing & grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", z)
            if (any(!missing & !valid)) stop(label, ": non-ISO/unknown date precision")
            ans <- rep(as.Date(NA), length(z))
            ans[valid] <- as.Date(z[valid], format = "%Y-%m-%d")
            if (any(valid & is.na(ans))) stop(label, ": invalid calendar date")
            ans
        }
        ev_inc_code <- function(x) {
            z <- toupper(gsub("[.[:space:]]", "", as.character(x)))
            z[is.na(x) | z %in% c("", "NA", "NAN", "NULL")] <- NA_character_
            z
        }
        ev_inc_evidence <- function(x) {
            is.list(x) && isTRUE(x$verified) && length(x$evidence) == 1L && !is.na(x$evidence) && nzchar(trimws(as.character(x$evidence)))
        }
        ev_inc_contract_example <- function(indir) {
            list(hospital = list(verified = FALSE, evidence = "", position_verified = FALSE, no_record_means_no_hospital_event = FALSE, 
                code_path = file.path(indir, "Rdata", "icd10.code0.rds"), date_path = file.path(indir, "Rdata", 
                  "icd10.date0.rds")), coverage = list(verified = FALSE, evidence = "", region_column = "", dates = c(England = NA_character_, 
                Scotland = NA_character_, Wales = NA_character_)), censor = list(verified = FALSE, evidence = "", 
                death_column = "date_death", lost_column = "date_lost", missing_death_means_alive = FALSE, missing_lost_means_not_lost = FALSE, 
                same_day_event = "unverified", age_precision = "unverified"), srd = list(verified = FALSE, evidence = "", 
                baseline_positions = character(), no_record_means_no_report = FALSE, code_path = file.path(indir, 
                  "Rdata", "srd.code0.rds")), opcs = list(verified = FALSE, evidence = "", available = TRUE, ihd_codes = character(), 
                position_verified = FALSE, no_record_means_no_procedure = FALSE, code_path = file.path(indir, "Rdata", 
                  "opcs4.code0.rds"), date_path = file.path(indir, "Rdata", "opcs4.date0.rds")), other_prior_ihd = list(verified = FALSE, 
                evidence = "", reviewed_sources = character(), status_columns = character(), date_columns = character(), 
                absent_sources = character(), date_missing_means_no_prior_record = FALSE), t2d = list(verified = FALSE, 
                evidence = "", algorithm = "", baseline_status_column = "", baseline_free_values = "no_diabetes", 
                type1_values = "type1", type2_values = "type2", unknown_values = c("uncertain", "unknown"), baseline_only_verified = FALSE, 
                inputs_and_rules_verified = FALSE))
        }
        ev_inc_specs <- function() {
            data.frame(endpoint = c("T2D", "IHD", "MI", "angina", "atherosclerosis"), phecode_context = c("250.2", 
                "411", "411.2", "411.3", "411.4"), endpoint_name = c("First hospital T2D record among baseline-unidentified diabetes", 
                "First hospital broad IHD record among baseline-unidentified IHD", "First hospital acute MI record among baseline-unidentified IHD", 
                "First hospital angina record among baseline-unidentified IHD", "First hospital coronary atherosclerosis record among baseline-unidentified IHD"), 
                pattern = c("^E11[0-9A-Z]*$", "^I2[0-5][0-9A-Z]*$", "^I2[12][0-9A-Z]*$", "^I20[0-9A-Z]*$", "^I251[0-9A-Z]*$"), 
                definition = c("Hospital ICD-10 E11 (all subcodes); not all-source incident capture", "Hospital ICD-10 I20-I25 (all subcodes); not hard MI", 
                  "Hospital ICD-10 I21/I22 (all subcodes); not PheCode 411.2 or legacy date_mi", "Hospital ICD-10 I20 (all subcodes)", 
                  "Hospital ICD-10 I25.1 (all subcodes)"), test_family = c("longitudinal_core", "longitudinal_core", 
                  rep("longitudinal_secondary", 3)), planned_tests = c(2L, 2L, 3L, 3L, 3L), definition_reference = c("https://doi.org/10.1371/journal.pone.0162388", 
                  rep("https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=41270", 4)), stringsAsFactors = FALSE)
        }
        ev_inc_read_pair <- function(code_path, date_path, eids, positions = NULL) {
            for (p in c(code_path, date_path)) if (length(p) != 1L || !file.exists(p)) stop("Missing validated source: ", 
                paste(p, collapse = ","))
            codes <- as.data.frame(readRDS(code_path), stringsAsFactors = FALSE)
            dates <- as.data.frame(readRDS(date_path), stringsAsFactors = FALSE)
            if (!"eid" %in% names(codes) || !"eid" %in% names(dates)) stop("Source lacks eid")
            if (anyNA(codes$eid) || anyNA(dates$eid) || anyDuplicated(codes$eid) || anyDuplicated(dates$eid)) stop("Source eid missing/duplicate")
            if (!identical(as.character(codes$eid), as.character(dates$eid))) stop("Code/date source eid orders differ; no unchecked positional join")
            cn <- setdiff(names(codes), "eid")
            dn <- setdiff(names(dates), "eid")
            if (!length(cn) || !all(grepl("^i[0-9]+$", cn)) || !all(grepl("^d[0-9]+$", dn)) || !identical(sub("^i", 
                "", cn), sub("^d", "", dn))) stop("Code/date array positions not identical iN/dN")
            if (!is.null(positions)) {
                if (anyDuplicated(positions) || any(!positions %in% cn)) stop("Invalid source positions")
                dn <- paste0("d", sub("^i", "", positions))
                cn <- positions
            }
            idx <- match(as.character(eids), as.character(codes$eid))
            list(codes = codes[idx, cn, drop = FALSE], dates = dates[idx, dn, drop = FALSE], linked = !is.na(idx), 
                positions = cn)
        }
        ev_inc_scan <- function(pair, pattern, entry) {
            n <- length(entry)
            first <- rep(as.Date(NA), n)
            prior <- rep(FALSE, n)
            undated <- rep(FALSE, n)
            for (j in seq_len(ncol(pair$codes))) {
                code <- ev_inc_code(pair$codes[[j]])
                positive <- !is.na(code) & grepl(pattern, code, perl = TRUE)
                if (!any(positive)) next
                date <- rep(as.Date(NA), n)
                date[positive] <- ev_inc_date(pair$dates[[j]][positive], paste0("source date position ", j))
                undated <- undated | (positive & is.na(date))
                prior <- prior | (positive & !is.na(date) & !is.na(entry) & date <= entry)
                newer <- positive & !is.na(date) & (is.na(first) | date < first)
                first[newer] <- date[newer]
            }
            list(first = first, prior = prior, undated = undated, known = pair$linked & !undated)
        }
        ev_inc_srd <- function(contract, eids) {
            p <- contract$code_path
            if (length(p) != 1L || !file.exists(p)) stop("Missing baseline SRD code source")
            x <- as.data.frame(readRDS(p), stringsAsFactors = FALSE)
            if (!"eid" %in% names(x) || anyNA(x$eid) || anyDuplicated(x$eid)) stop("Baseline SRD source eid missing/duplicate")
            cols <- contract$baseline_positions
            if (!length(cols) || anyDuplicated(cols) || any(!cols %in% names(x))) stop("Baseline SRD instance-0 positions absent/unverified")
            idx <- match(as.character(eids), as.character(x$eid))
            present <- !is.na(idx)
            ihd <- dm <- unknown <- rep(FALSE, length(eids))
            for (col in cols) {
                code <- as.character(x[[col]][idx])
                vals <- suppressWarnings(as.numeric(code))
                positive <- !is.na(vals)
                ihd <- ihd | (positive & vals %in% c(1074, 1075))
                dm <- dm | (positive & vals %in% c(1220, 1222, 1223))
                unknown <- unknown | (positive & vals %in% c(-1, -3))
                bad <- !is.na(code) & nzchar(trimws(code)) & !toupper(trimws(code)) %in% c("NA", "NAN", "NULL") & 
                  is.na(vals)
                if (any(bad)) stop("Unparsed baseline SRD code: source needs audit")
            }
            list(ihd = ihd, dm = dm, known = present & !unknown, unknown = unknown)
        }
        ev_inc_column <- function(ctx, name) {
            if (length(name) != 1L || is.na(name) || !nzchar(name)) stop("Unconfigured column")
            if (name %in% names(ctx$d)) return(ctx$d[[name]])
            if (!is.null(ctx$raw_all) && name %in% names(ctx$raw_all)) {
                raw <- ctx$raw_all
                if (anyNA(raw$eid) || anyDuplicated(raw$eid)) stop("raw_all eid not unique")
                return(raw[[name]][match(as.character(ctx$d$eid), as.character(raw$eid))])
            }
            stop("Missing audited column: ", name)
        }
        ev_inc_censor <- function(ctx, cfg, entry) {
            region <- as.character(ev_inc_column(ctx, cfg$coverage$region_column))
            cutoffs <- cfg$coverage$dates
            if (is.null(names(cutoffs)) || anyDuplicated(names(cutoffs))) stop("Coverage dates must be region-named")
            coverage <- ev_inc_date(unname(cutoffs[match(region, names(cutoffs))]), "coverage")
            admin <- ev_inc_date(ctx$config$admin_follow_end, "administrative cutoff")
            if (length(admin) != 1L || is.na(admin) || admin > as.Date("2023-04-01")) stop("Administrative cutoff absent or later than locked 2023-04-01")
            death <- ev_inc_date(ev_inc_column(ctx, cfg$censor$death_column), "death")
            lost <- ev_inc_date(ev_inc_column(ctx, cfg$censor$lost_column), "lost")
            stopifnot(length(death) == length(entry), length(lost) == length(entry))
            end <- pmin(rep(admin, length(entry)), coverage, death, lost, na.rm = TRUE)
            end[is.na(coverage)] <- as.Date(NA)
            list(end = end, death = death, lost = lost, coverage = coverage, known = !is.na(coverage), version = paste("admin<=2023-04-01;audited uniform conservative hospital horizon;", 
                "min(event,death,lost,coverage,admin);same-day", cfg$censor$same_day_event), missing_region_n = sum(is.na(coverage)))
        }
        ev_inc_add_years <- function(x, years = 2L) {
            ans <- rep(as.Date(NA), length(x))
            valid <- !is.na(x)
            yr <- as.integer(format(x[valid], "%Y")) + as.integer(years)
            md <- format(x[valid], "%m-%d")
            leap <- (yr%%4L == 0L & yr%%100L != 0L) | yr%%400L == 0L
            md[md == "02-29" & !leap] <- "02-28"
            ans[valid] <- as.Date(paste0(yr, "-", md), format = "%Y-%m-%d")
            ans
        }
        ev_inc_risk <- function(d, entry, event_date, censor, prior, known, landmark = FALSE, same_day = "event") {
            stopifnot(nrow(d) == length(entry), length(entry) == length(event_date), length(entry) == length(censor$end), 
                length(entry) == length(prior), length(entry) == length(known))
            if (anyNA(prior) || anyNA(known)) stop("Prior/known masks cannot silently contain NA")
            if (!"S" %in% names(d)) stop("Missing locked full-sibship size")
            first_last <- !is.na(d$X) & d$X %in% c(0, 1)
            supported <- sort(unique(d$S[first_last & !is.na(d$S) & d$S >= 2 & d$S == floor(d$S)]))
            supported <- supported[vapply(supported, function(s) all(c(0, 1) %in% d$X[first_last & !is.na(d$S) & 
                d$S == s]), logical(1))]
            mask <- first_last
            flow <- list()
            add <- function(step, keep) {
                before <- sum(mask)
                mask <<- mask & keep
                flow[[length(flow) + 1L]] <<- data.frame(step = step, n_before = before, n_excluded = before - 
                  sum(mask), n_remaining = sum(mask), stringsAsFactors = FALSE)
            }
            add("outcome_independent_supported_exact_sibship_size", !is.na(d$S) & d$S %in% supported)
            add("baseline_source_observable", known)
            add("baseline_without_identified_disease", !prior)
            add("baseline_date_present", !is.na(entry))
            add("valid_age_at_entry", !is.na(d$age_entry) & is.finite(d$age_entry) & d$age_entry > 0)
            add("regional_followup_coverage_known", censor$known)
            add("positive_followup", !is.na(censor$end) & !is.na(entry) & censor$end > entry)
            start <- entry
            if (landmark) {
                start <- ev_inc_add_years(entry, 2L)
                alive <- is.na(censor$death) | censor$death > start
                observable <- !is.na(censor$end) & censor$end > start
                event_free <- is.na(event_date) | event_date > start
                add("alive_at_2calendar_year_landmark", alive & !is.na(start))
                add("observable_beyond_landmark", observable)
                add("event_free_through_landmark_day", event_free)
            }
            event <- !is.na(event_date) & !is.na(start) & !is.na(censor$end) & event_date > start & event_date <= 
                censor$end
            if (same_day == "death") event <- event & (is.na(censor$death) | event_date < censor$death)
            if (!same_day %in% c("event", "death")) stop("Unverified event/death same-day rule")
            exit <- censor$end
            exit[event] <- event_date[event]
            out <- d[mask, , drop = FALSE]
            out$event <- as.integer(event[mask])
            out$followup_years <- as.numeric(exit[mask] - start[mask])/365.25
            out$age_entry <- d$age_entry[mask] + as.numeric(start[mask] - entry[mask])/365.25
            out$age_exit <- d$age_entry[mask] + as.numeric(exit[mask] - entry[mask])/365.25
            if (anyDuplicated(out$eid) || anyNA(out$eid) || any(out$followup_years <= 0) || any(out$age_exit <= 
                out$age_entry)) stop("Invalid Cox interval/eid after risk-set construction")
            list(data = out, flow = do.call(rbind, flow), n_eligible = sum(mask), events_after_death = sum(mask & 
                !is.na(event_date) & !is.na(censor$death) & event_date > censor$death & event_date > entry, na.rm = TRUE), 
                events_after_death_all = sum(!is.na(event_date) & !is.na(censor$death) & event_date > censor$death & 
                  event_date > entry, na.rm = TRUE), event_on_death = sum(mask & !is.na(event_date) & !is.na(censor$death) & 
                  event_date == censor$death & event_date > entry, na.rm = TRUE))
        }
        ev_inc_terms <- function(ctx, model) {
            if (model == "time_on_study") {
                if (is.null(ctx$primary_terms) || !length(ctx$primary_terms)) stop("Missing frozen primary terms")
                return(paste(unique(c("X", ctx$primary_terms)), collapse = " + "))
            }
            if (length(ctx$birth_term) != 1L || !nzchar(ctx$birth_term)) stop("Missing frozen birth spline")
            paste(c("X", "factor(S)", "factor(sex)", "factor(ethnicity)", "factor(center)", ctx$birth_term), collapse = " + ")
        }
        ev_inc_fit <- function(ctx, dat, spec, model, censor_version) {
            fit_file <- file.path(ctx$out, "cache", paste0("incident_fit_", spec$endpoint, "_", model, ".rds"))
            if (file.exists(fit_file)) stop("Protected existing fit without a matching task checkpoint: ", fit_file)
            response <- if (model == "time_on_study") "survival::Surv(followup_years, event)" else "survival::Surv(age_entry, age_exit, event)"
            formula <- stats::as.formula(paste(response, "~", ev_inc_terms(ctx, model)))
            environment(formula) <- globalenv()
            mm <- stats::model.matrix(formula, dat)
            mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
            if (!ncol(mm) || qr(mm)$rank < ncol(mm)) stop("Cox design rank deficient; no silent term removal")
            cond <- tryCatch(kappa(mm, exact = FALSE), error = function(e) NA_real_)
            if (length(unique(dat$X)) != 2L || sum(dat$event) < 20L || sum(dat$event == 0L) < 20L) stop("Prespecified estimability rule fails: both birth groups, >=20 events and non-events required")
            warnings <- character()
            fit <- withCallingHandlers(survival::coxph(formula, data = dat, ties = "efron", x = TRUE, y = TRUE, 
                singular.ok = FALSE, na.action = stats::na.fail), warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
            tab <- summary(fit)$coefficients
            if (!"X" %in% rownames(tab) || any(!is.finite(stats::coef(fit))) || any(!is.finite(sqrt(diag(stats::vcov(fit))))) || 
                any(grepl("infinite|converg|iteration|singular", warnings, ignore.case = TRUE))) stop("Cox unstable/non-convergent/infinite estimate: ", 
                paste(warnings, collapse = "; "))
            if (length(fit$n) != 1L || fit$n != nrow(dat) || fit$n != length(unique(dat$eid)) || fit$nevent != 
                sum(dat$event) || nrow(fit$x) != nrow(dat)) stop("Cox actual N/events/eid invariant failed")
            beta <- tab["X", "coef"]
            se <- tab["X", "se(coef)"]
            p <- tab["X", "Pr(>|z|)"]
            if (!all(is.finite(c(beta, se, p))) || se <= 0) stop("X coefficient/SE/P invalid")
            ph <- tryCatch(survival::cox.zph(fit, terms = TRUE, global = TRUE), error = identity)
            ph_status <- "diagnostic_failed"
            ph_x <- ph_global <- NA_real_
            ph_table <- NULL
            if (!inherits(ph, "error")) {
                ph_table <- data.frame(endpoint = spec$endpoint, model = model, term = rownames(ph$table), ph$table, 
                  row.names = NULL, check.names = FALSE)
                if ("X" %in% rownames(ph$table)) ph_x <- ph$table["X", "p"]
                if ("GLOBAL" %in% rownames(ph$table)) ph_global <- ph$table["GLOBAL", "p"]
                ph_status <- if (is.finite(ph_x) && is.finite(ph_global)) "diagnostics_available_review_plot" else "diagnostic_incomplete"
                ph_file <- file.path(ctx$out, "figures", paste0("Schoenfeld_", spec$endpoint, "_", model, ".pdf"))
                dir.create(dirname(ph_file), recursive = TRUE, showWarnings = FALSE)
                if (file.exists(ph_file)) stop("Protected existing PH figure: ", ph_file)
                tmp <- paste0(ph_file, ".tmp.", Sys.getpid())
                grDevices::pdf(tmp, width = 6.4, height = 5)
                tryCatch({
                  if ("X" %in% colnames(ph$y)) graphics::plot(ph, var = match("X", colnames(ph$y)), main = paste(spec$endpoint, 
                    model, "X PH diagnostic")) else graphics::plot.new()
                }, finally = grDevices::dev.off())
                if (!file.rename(tmp, ph_file)) stop("Cannot install PH figure: ", ph_file)
            } else warnings <- c(warnings, paste("cox.zph:", conditionMessage(ph)))
            stats <- data.frame(endpoint = spec$endpoint, endpoint_name = spec$endpoint_name, outcome = spec$phecode_context, 
                metric = "HR", exposure = "X", scope = "hospital_first_record", model = model, primary_flag = model == 
                  "primary_age", time_scale = if (model == "time_on_study") "time_on_study" else "attained_age", 
                landmark_years = if (model == "landmark_2y") 2L else 0L, status = "ok", technical_status = if (ph_status == 
                  "diagnostic_failed") "needs_review" else "fit_complete", reason = "", formula = paste(deparse(formula, 
                  width.cutoff = 500), collapse = " "), N = as.integer(fit$n), model_rows = nrow(fit$x), unique_eid = length(unique(dat$eid)), 
                sample_hash = digest::digest(sort(as.character(dat$eid)), algo = "sha256"), events = as.integer(fit$nevent), 
                N_first = sum(dat$X == 0), N_last = sum(dat$X == 1), events_first = sum(dat$X == 0 & dat$event == 
                  1), events_last = sum(dat$X == 1 & dat$event == 1), person_years_first = sum(dat$followup_years[dat$X == 
                  0]), person_years_last = sum(dat$followup_years[dat$X == 1]), followup_median = stats::median(dat$followup_years), 
                followup_q25 = unname(stats::quantile(dat$followup_years, 0.25)), followup_q75 = unname(stats::quantile(dat$followup_years, 
                  0.75)), beta = beta, SE = se, HR = exp(beta), CI_low = exp(beta - 1.96 * se), CI_high = exp(beta + 
                  1.96 * se), P = p, p_adjusted = NA_real_, test_family = spec$test_family, planned_tests = spec$planned_tests, 
                correction = if (spec$test_family == "longitudinal_core") "Holm n=2" else "BH n=3", ph_X_P = ph_x, 
                ph_global_P = ph_global, ph_status = ph_status, condition_number = cond, censor_version = censor_version, 
                warnings = paste(unique(warnings), collapse = "; "), inference_scope = "Post-discovery internal robustness; cause-specific hospital-record HR, not cumulative risk", 
                stringsAsFactors = FALSE)
            ev_save(fit, fit_file)
            list(table = stats, ph = ph_table)
        }
        ev_inc_blank <- function(spec, model, status, reason) {
            data.frame(endpoint = spec$endpoint, endpoint_name = spec$endpoint_name, model = model, outcome = spec$phecode_context, 
                metric = "HR", exposure = "X", scope = "hospital_first_record", primary_flag = model == "primary_age", 
                time_scale = if (model == "time_on_study") "time_on_study" else "attained_age", landmark_years = if (model == 
                  "landmark_2y") 2L else 0L, status = status, technical_status = if (status == "blocked") "data_insufficient" else "needs_repair", 
                reason = reason, N = NA_integer_, model_rows = NA_integer_, unique_eid = NA_integer_, events = NA_integer_, 
                N_first = NA_integer_, N_last = NA_integer_, events_first = NA_integer_, events_last = NA_integer_, 
                person_years_first = NA_real_, person_years_last = NA_real_, followup_median = NA_real_, followup_q25 = NA_real_, 
                followup_q75 = NA_real_, beta = NA_real_, SE = NA_real_, HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, 
                P = NA_real_, p_adjusted = NA_real_, test_family = spec$test_family, planned_tests = spec$planned_tests, 
                correction = if (spec$test_family == "longitudinal_core") "Holm n=2" else "BH n=3", ph_X_P = NA_real_, 
                ph_global_P = NA_real_, ph_status = "not_evaluated", condition_number = NA_real_, censor_version = "unverified", 
                warnings = "", inference_scope = "Not evaluated; not evidence of no association", stringsAsFactors = FALSE)
        }
        ev_inc_gate <- function(cfg, endpoint) {
            reasons <- character()
            need <- function(ok, msg) if (!isTRUE(ok)) reasons <<- c(reasons, msg)
            h <- cfg$hospital
            c <- cfg$coverage
            z <- cfg$censor
            need(ev_inc_evidence(h) && isTRUE(h$position_verified) && isTRUE(h$no_record_means_no_hospital_event), 
                "Audit hospital source/version, iN/dN alignment and no-record semantics")
            need(ev_inc_evidence(c) && length(c$region_column) == 1L && nzchar(c$region_column) && !is.null(names(c$dates)), 
                "Supply verified local hospital completeness cutoff per source/region, with region column")
            need(ev_inc_evidence(z) && isTRUE(z$missing_death_means_alive) && isTRUE(z$missing_lost_means_not_lost) && 
                z$same_day_event %in% c("event", "death") && length(z$age_precision) == 1L && !z$age_precision %in% 
                c("", "unverified"), "Verify uncapped death/lost provenance, absent-record semantics, same-day rule and age precision")
            if (endpoint == "T2D") {
                t <- cfg$t2d
                need(ev_inc_evidence(t) && isTRUE(t$baseline_only_verified) && isTRUE(t$inputs_and_rules_verified) && 
                  length(t$algorithm) == 1L && nzchar(t$algorithm) && length(t$baseline_status_column) == 1L && 
                  nzchar(t$baseline_status_column) && length(t$baseline_free_values) > 0L && length(t$type1_values) > 
                  0L && length(t$type2_values) > 0L && length(t$unknown_values) > 0L, "Verify Eastwood appendix or explicit local baseline DM algorithm/inputs; distinguish type1/type2/uncertain; do not use dm.yes")
            } else {
                s <- cfg$srd
                o <- cfg$opcs
                p <- cfg$other_prior_ihd
                need(ev_inc_evidence(s) && length(s$baseline_positions) > 0L && isTRUE(s$no_record_means_no_report), 
                  "Audit instance-0 SRD positions and missing semantics for prior 1074/1075")
                need(ev_inc_evidence(o) && identical(o$available, FALSE) || ev_inc_evidence(o) && isTRUE(o$available) && 
                  isTRUE(o$position_verified) && length(o$ihd_codes) > 0L, "Review OPCS prior coronary disease/revascularisation codes, or document genuinely unavailable source")
                need(!isTRUE(o$available) || isTRUE(o$no_record_means_no_procedure), "Verify absent OPCS array entries mean no recorded procedure, not unobservable data")
                need(ev_inc_evidence(p) && length(p$reviewed_sources) > 0L, "Review ICD9/GP/UKB algorithm/self-report prior-IHD evidence; specify verified flags/dates or unavailable sources")
                need(!length(p$date_columns) || isTRUE(p$date_missing_means_no_prior_record), "Verify missing derived prior-IHD dates mean no recorded disease; otherwise use 0/1/NA source status flags")
            }
            unique(reasons)
        }
        ev_inc_endpoint <- function(ctx, cfg, spec, hospital, censor, srd = NULL, opcs = NULL) {
            d <- as.data.frame(ctx$d)
            entry <- ev_inc_date(ev_inc_column(ctx, "date_attend"), "baseline measurement")
            event <- ev_inc_scan(hospital, spec$pattern, entry)
            prior <- event$prior
            known <- event$known
            baseline_counts <- list(hospital_endpoint = event$prior, hospital_undated = event$undated)
            if (spec$endpoint == "T2D") {
                dm <- ev_inc_scan(hospital, "^E1[0-4][0-9A-Z]*$", entry)
                t <- cfg$t2d
                status <- as.character(ev_inc_column(ctx, t$baseline_status_column))
                free_values <- as.character(t$baseline_free_values)
                excluded_values <- unique(c(as.character(t$type1_values), as.character(t$type2_values), as.character(t$unknown_values)))
                if (!length(free_values) || anyNA(free_values) || any(free_values %in% excluded_values)) stop("Baseline DM algorithm status sets overlap/empty")
                free <- !is.na(status) & status %in% free_values
                diagnosed <- !is.na(status) & status %in% c(t$type1_values, t$type2_values)
                prior <- prior | dm$prior | diagnosed
                known <- known & dm$known & (free | diagnosed)
                baseline_counts$hospital_any_dm <- dm$prior
                baseline_counts$algorithm_not_verified_diabetes_free <- !free
                baseline_counts$algorithm_type1 <- !is.na(status) & status %in% t$type1_values
                baseline_counts$algorithm_type2 <- !is.na(status) & status %in% t$type2_values
                baseline_counts$algorithm_uncertain <- is.na(status) | status %in% t$unknown_values
            } else {
                ihd <- ev_inc_scan(hospital, "^I2[0-5][0-9A-Z]*$", entry)
                prior <- prior | ihd$prior | srd$ihd
                known <- known & ihd$known & srd$known
                baseline_counts$hospital_any_ihd <- ihd$prior
                baseline_counts$baseline_srd_ihd <- srd$ihd
                baseline_counts$baseline_srd_unknown <- !srd$known
                if (!is.null(opcs)) {
                  codes <- ev_inc_code(cfg$opcs$ihd_codes)
                  if (anyNA(codes) || any(!grepl("^[A-Z][0-9A-Z]{1,6}$", codes))) stop("Prior OPCS mapping requires exact codes, not unaudited regular expressions")
                  regex <- paste0("^(?:", paste(codes, collapse = "|"), ")$")
                  surgery <- ev_inc_scan(opcs, regex, entry)
                  prior <- prior | surgery$prior
                  known <- known & surgery$known
                  baseline_counts$prior_coronary_procedure <- surgery$prior
                  baseline_counts$undated_coronary_procedure <- surgery$undated
                }
                for (col in cfg$other_prior_ihd$status_columns) {
                  v <- ev_inc_column(ctx, col)
                  if (any(!is.na(v) & !v %in% c(0, 1, FALSE, TRUE))) stop("Prior-IHD flag not 0/1/NA: ", col)
                  prior <- prior | (!is.na(v) & v == 1)
                  known <- known & !is.na(v)
                  baseline_counts[[paste0("prior_flag_", col)]] <- !is.na(v) & v == 1
                }
                for (col in cfg$other_prior_ihd$date_columns) {
                  date <- ev_inc_date(ev_inc_column(ctx, col), col)
                  positive <- !is.na(date) & !is.na(entry) & date <= entry
                  prior <- prior | positive
                  baseline_counts[[paste0("prior_date_", col)]] <- positive
                }
            }
            risk <- ev_inc_risk(d, entry, event$first, censor, prior, known, same_day = cfg$censor$same_day_event)
            if (risk$events_after_death > 0L) stop("Endpoint dates after verified death: ", risk$events_after_death, 
                "; source reconciliation required")
            union_formula <- stats::as.formula(paste("~", ev_inc_terms(ctx, "primary_age"), "+", ev_inc_terms(ctx, 
                "time_on_study")))
            base <- risk$data
            mf <- stats::model.frame(union_formula, base, na.action = stats::na.pass)
            complete <- stats::complete.cases(mf)
            for (v in names(mf)) if (is.numeric(mf[[v]])) {
                z <- mf[[v]]
                finite <- if (is.matrix(z)) rowSums(!is.finite(z)) == 0L else is.finite(z)
                complete <- complete & finite
            }
            if (!is.logical(complete) || !is.null(dim(complete)) || length(complete) != nrow(base)) stop("Cox complete-case mask must be one logical value per participant")
            extra <- data.frame(step = "main_time_covariate_union_complete", n_before = nrow(base), n_excluded = sum(!complete), 
                n_remaining = sum(complete), stringsAsFactors = FALSE)
            base <- droplevels(base[complete, , drop = FALSE])
            flow <- rbind(risk$flow, extra)
            flow$endpoint <- spec$endpoint
            flow$model <- "primary_age_and_time_on_study"
            baseline <- data.frame(step = paste0("baseline_source_count_nonmutuallyexclusive:", names(baseline_counts)), 
                n_before = sum(!is.na(d$X)), n_excluded = vapply(baseline_counts, function(z) sum(z & !is.na(d$X), 
                  na.rm = TRUE), integer(1)), n_remaining = NA_integer_, endpoint = spec$endpoint, model = "baseline_audit", 
                stringsAsFactors = FALSE)
            rows <- list()
            ph <- list()
            for (model in c("primary_age", "time_on_study", "landmark_2y")) {
                model_data <- base
                if (model == "landmark_2y") {
                  lrisk <- ev_inc_risk(d, entry, event$first, censor, prior, known, landmark = TRUE, same_day = cfg$censor$same_day_event)
                  model_data <- lrisk$data[as.character(lrisk$data$eid) %in% as.character(base$eid), , drop = FALSE]
                  model_data <- droplevels(model_data)
                  lflow <- lrisk$flow
                  lflow$endpoint <- spec$endpoint
                  lflow$model <- model
                  flow <- rbind(flow, lflow, data.frame(step = "baseline_union_complete_landmark_intersection", 
                    n_before = nrow(lrisk$data), n_excluded = nrow(lrisk$data) - nrow(model_data), n_remaining = nrow(model_data), 
                    endpoint = spec$endpoint, model = model))
                }
                result <- tryCatch(ev_run_task(ctx, paste0("incident_", spec$endpoint, "_", model), function() ev_inc_fit(ctx, 
                  model_data, spec, model, censor$version)), error = function(e) {
                  if (inherits(e, "ev_budget_end")) stop(e)
                  e
                })
                if (inherits(result, "error")) {
                  rows[[model]] <- ev_inc_blank(spec, model, if (inherits(result, "ev_budget_end")) "deferred" else "failed", 
                    conditionMessage(result))
                } else if (is.null(result) || !is.list(result) || is.null(result$table)) {
                  rows[[model]] <- ev_inc_blank(spec, model, "deferred", "Task budget/checkpoint returned no model result")
                } else {
                  rows[[model]] <- result$table
                  ph[[model]] <- result$ph
                }
            }
            list(table = ev_bind(rows), flow = rbind(flow, baseline), ph = ev_bind(ph), date_audit = data.frame(endpoint = spec$endpoint, 
                events_after_death = risk$events_after_death, events_after_death_all = risk$events_after_death_all, 
                event_on_death = risk$event_on_death, missing_region_n = censor$missing_region_n))
        }
        ev_incident <- function(ctx, audit_only = FALSE) {
            existing_table <- file.path(ctx$out, "tables", "incident_core.csv")
            if (audit_only && file.exists(existing_table)) {
                previous <- utils::read.csv(existing_table, stringsAsFactors = FALSE, check.names = FALSE)
                if (("HR" %in% names(previous) && any(is.finite(suppressWarnings(as.numeric(previous$HR))))) || 
                  ("status" %in% names(previous) && any(previous$status == "ok", na.rm = TRUE))) stop("Audit-only cannot replace completed incident results; use a new output directory")
            }
            specs <- ev_inc_specs()
            cfg <- ctx$config$incident
            if (is.null(cfg)) cfg <- ev_inc_contract_example(ctx$indir)
            manifests <- tables <- flows <- ph <- date_audit <- list()
            shared <- NULL
            shared_error <- NULL
            for (i in seq_len(nrow(specs))) {
                spec <- specs[i, , drop = FALSE]
                reasons <- ev_inc_gate(cfg, spec$endpoint)
                models <- c("primary_age", "time_on_study", "landmark_2y")
                gate_status <- if (length(reasons)) "blocked" else "audited_contract_ready"
                if (!length(reasons) && audit_only) {
                  gate_status <- "audited_contract_not_runtime_checked_not_run"
                  tables[[spec$endpoint]] <- ev_bind(lapply(models, function(model) ev_inc_blank(spec, model, "not_run", 
                    "Audit-only invocation; no source arrays read and no models fitted")))
                  flows[[spec$endpoint]] <- data.frame(endpoint = spec$endpoint, model = "all", step = "audit_only_no_risk_set_constructed", 
                    n_before = sum(!is.na(ctx$d$X)), n_excluded = NA_integer_, n_remaining = NA_integer_, reason = "No model requested")
                }
                if (!length(reasons) && !audit_only) {
                  if (is.null(shared) && is.null(shared_error)) {
                    shared <- tryCatch({
                      entry <- ev_inc_date(ev_inc_column(ctx, "date_attend"), "baseline measurement")
                      list(hospital = ev_inc_read_pair(cfg$hospital$code_path, cfg$hospital$date_path, ctx$d$eid), 
                        censor = ev_inc_censor(ctx, cfg, entry))
                    }, error = function(e) {
                      if (inherits(e, "ev_budget_end")) stop(e)
                      e
                    })
                    if (inherits(shared, "error")) {
                      shared_error <- conditionMessage(shared)
                      shared <- NULL
                    }
                  }
                  if (!is.null(shared_error)) reasons <- shared_error
                  if (!length(reasons)) {
                    result <- tryCatch({
                      srd <- opcs <- NULL
                      if (spec$endpoint != "T2D") {
                        if (is.null(shared$srd)) shared$srd <- ev_inc_srd(cfg$srd, ctx$d$eid)
                        srd <- shared$srd
                        if (isTRUE(cfg$opcs$available)) {
                          if (is.null(shared$opcs)) shared$opcs <- ev_inc_read_pair(cfg$opcs$code_path, cfg$opcs$date_path, 
                            ctx$d$eid)
                          opcs <- shared$opcs
                        }
                      }
                      ev_inc_endpoint(ctx, cfg, spec, shared$hospital, shared$censor, srd, opcs)
                    }, error = function(e) {
                      if (inherits(e, "ev_budget_end")) stop(e)
                      e
                    })
                    if (inherits(result, "error")) reasons <- conditionMessage(result) else {
                      tables[[spec$endpoint]] <- result$table
                      flows[[spec$endpoint]] <- result$flow
                      ph[[spec$endpoint]] <- result$ph
                      date_audit[[spec$endpoint]] <- result$date_audit
                    }
                  }
                }
                if (length(reasons)) {
                  gate_status <- "blocked"
                  tables[[spec$endpoint]] <- ev_bind(lapply(models, function(model) ev_inc_blank(spec, model, "blocked", 
                    paste(reasons, collapse = "; "))))
                  flows[[spec$endpoint]] <- data.frame(endpoint = spec$endpoint, model = "all", step = "blocked_before_risk_set", 
                    n_before = sum(!is.na(ctx$d$X)), n_excluded = NA_integer_, n_remaining = NA_integer_, reason = paste(reasons, 
                      collapse = "; "))
                }
                manifests[[spec$endpoint]] <- data.frame(endpoint = spec$endpoint, endpoint_name = spec$endpoint_name, 
                  definition = spec$definition, normalized_code_regex = spec$pattern, source = "Hospital ICD-10 summary diagnosis/date arrays", 
                  source_code_path = if (length(cfg$hospital$code_path)) cfg$hospital$code_path else "unconfigured", 
                  source_date_path = if (length(cfg$hospital$date_path)) cfg$hospital$date_path else "unconfigured", 
                  source_evidence = paste(cfg$hospital$evidence, collapse = "; "), srd_baseline_positions = paste(cfg$srd$baseline_positions, 
                    collapse = "; "), opcs_prior_ihd_exact_codes = paste(cfg$opcs$ihd_codes, collapse = "; "), 
                  other_prior_ihd_reviewed_sources = paste(cfg$other_prior_ihd$reviewed_sources, collapse = "; "), 
                  other_prior_ihd_absent_sources = paste(cfg$other_prior_ihd$absent_sources, collapse = "; "), 
                  other_prior_ihd_status_columns = paste(cfg$other_prior_ihd$status_columns, collapse = "; "), 
                  other_prior_ihd_date_columns = paste(cfg$other_prior_ihd$date_columns, collapse = "; "), definition_reference = spec$definition_reference, 
                  baseline_exclusion = if (spec$endpoint == "T2D") "Verified baseline DM algorithm: exclude type1/type2/uncertain/unknown; prior E10-E14 and undated DM excluded" else "Any identified baseline IHD: hospital I20-I25 + baseline SRD 1074/1075 + verified procedures/other available sources; undated unknown excluded", 
                  baseline_algorithm = if (spec$endpoint == "T2D") paste(cfg$t2d$algorithm, collapse = "; ") else "Audited multi-source baseline IHD exclusion; not a clinical adjudication algorithm", 
                  baseline_evidence = if (spec$endpoint == "T2D") paste(cfg$t2d$evidence, collapse = "; ") else paste(c(cfg$srd$evidence, 
                    cfg$opcs$evidence, cfg$other_prior_ihd$evidence), collapse = "; "), time_origin = "Verified baseline exposure measurement date; future-response backfill prohibited", 
                  age_method = paste(cfg$censor$age_precision, collapse = "; "), end_rule = "First event/death/lost/verified source-region coverage/administrative cutoff, whichever earliest", 
                  administrative_cap = "2023-04-01; earlier verified local source coverage takes precedence", coverage_evidence = paste(cfg$coverage$evidence, 
                    collapse = "; "), censor_evidence = paste(cfg$censor$evidence, collapse = "; "), same_day_rule = paste(cfg$censor$same_day_event, 
                    collapse = "; "), planned_tests = spec$planned_tests, test_family = spec$test_family, status = gate_status, 
                  reason = paste(reasons, collapse = "; "), interpretation = "Post-discovery internal robustness; first hospital record, not complete incidence or causal inference", 
                  stringsAsFactors = FALSE)
            }
            tab <- ev_bind(tables)
            for (model in c("primary_age", "time_on_study", "landmark_2y")) for (family in c("longitudinal_core", 
                "longitudinal_secondary")) {
                ix <- which(tab$model == model & tab$test_family == family & is.finite(tab$P))
                if (length(ix)) tab$p_adjusted[ix] <- stats::p.adjust(tab$P[ix], method = if (family == "longitudinal_core") "holm" else "BH", 
                  n = if (family == "longitudinal_core") 2L else 3L)
            }
            manifest <- ev_bind(manifests)
            flow <- ev_bind(flows)
            ev_write(manifest, file.path(ctx$out, "tables", "outcome_manifest.csv"))
            ev_write(tab, file.path(ctx$out, "tables", "incident_core.csv"))
            ev_write(flow, file.path(ctx$out, "tables", "incident_flow.csv"))
            ph_out <- ev_bind(ph)
            if (!nrow(ph_out)) ph_out <- data.frame(endpoint = character(), model = character(), term = character(), 
                chisq = numeric(), df = numeric(), p = numeric())
            dates_out <- ev_bind(date_audit)
            if (!nrow(dates_out)) dates_out <- data.frame(endpoint = character(), events_after_death = integer(), 
                events_after_death_all = integer(), event_on_death = integer(), missing_region_n = integer())
            ev_write(ph_out, file.path(ctx$out, "tables", "incident_ph.csv"))
            ev_write(dates_out, file.path(ctx$out, "audit", "incident_date_audit.csv"))
            ev_record_status(ctx, tab)
            invisible(list(table = tab, manifest = manifest, flow = flow))
        }
        ev_incident_audit <- function(ctx) ev_incident(ctx, audit_only = TRUE)
        ev_incident_tests <- function() {
            d <- data.frame(eid = as.character(1:6), X = c(0, 1, 0, 1, 0, 1), age_entry = rep(50, 6), S = rep(2L, 
                6))
            entry <- as.Date(rep("2020-01-01", 6))
            event <- as.Date(c("2019-01-01", "2021-01-01", "2022-01-01", "2023-01-01", NA, "2025-01-01"))
            censor <- list(end = as.Date(rep("2024-01-01", 6)), death = rep(as.Date(NA), 6), known = rep(TRUE, 
                6))
            r <- ev_inc_risk(d, entry, event, censor, prior = c(TRUE, rep(FALSE, 5)), known = rep(TRUE, 6))
            stopifnot(!"1" %in% r$data$eid, nrow(r$data) == 5L, sum(r$data$event) == 3L, r$data$event[r$data$eid == 
                "6"] == 0L, length(unique(r$data$eid)) == nrow(r$data))
            l <- ev_inc_risk(d, entry, event, censor, prior = c(TRUE, rep(FALSE, 5)), known = rep(TRUE, 6), landmark = TRUE)
            stopifnot(identical(l$data$eid, c("4", "5", "6")), nrow(l$data) == 3L, sum(l$data$event) == 1L, isTRUE(all.equal(l$data$followup_years[l$data$eid == 
                "5"], 730/365.25)), isTRUE(all.equal(l$data$age_entry[1], 50 + 731/365.25)), ev_inc_add_years(as.Date("2020-02-29")) == 
                as.Date("2022-02-28"))
            unknown <- ev_inc_risk(d, entry, event, censor, prior = c(TRUE, rep(FALSE, 5)), known = c(TRUE, TRUE, 
                TRUE, TRUE, FALSE, TRUE))
            stopifnot(!"5" %in% unknown$data$eid)
            toy <- data.frame(eid = as.character(seq_len(120)), X = rep(0:1, 60), age_entry = rep(40, 120), age_exit = 40 + 
                (seq_len(120)%%17 + 1)/10, event = as.integer(seq_len(120)%%3 == 0))
            fit <- survival::coxph(survival::Surv(age_entry, age_exit, event) ~ X, data = toy, x = TRUE, y = TRUE, 
                na.action = stats::na.fail)
            stopifnot(fit$n == nrow(toy), fit$n == length(unique(toy$eid)), fit$nevent == sum(toy$event), nrow(fit$x) == 
                nrow(toy))
            invisible(TRUE)
        }
    }), envir = context)
}

