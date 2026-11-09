# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

ip_source_ids <- function(x, label = "source") {
    x <- as.data.frame(x)
    key <- intersect(c("eid", "id"), names(x))
    bi_assert(length(key) == 1L, paste("SOURCE_ID_COLUMN_INVALID", label))
    names(x)[names(x) == key] <- "eid"
    x$eid <- bi_isa_ids(x$eid, label)
    bi_assert(!anyDuplicated(names(x)), paste("SOURCE_DUPLICATE_COLUMNS", label))
    x
}

ip_source_num <- function(x) if (is.numeric(x) || is.logical(x)) as.numeric(x) else bi_num(x)

ip_source_read <- function(path, label) {
    x <- if (grepl("[.]rds$", path, ignore.case = TRUE)) 
        readRDS(path)
    else data.table::fread(path, data.table = FALSE, check.names = FALSE, showProgress = FALSE)
    ip_source_ids(x, label)
}

ip_source_manifest <- function(paths, authenticated = NULL) {
    if (!is.null(authenticated)) {
        authenticated <- as.data.frame(authenticated)
        if ("key" %in% names(authenticated) && !"role" %in% names(authenticated)) 
            names(authenticated)[names(authenticated) == "key"] <- "role"
        bi_assert(all(c("role", "path", "sha256", "bytes", "mtime") %in% names(authenticated)) && !anyDuplicated(authenticated$role), 
            "AUTHENTICATED_SOURCE_MANIFEST_INVALID")
    }
    bi_bind(lapply(names(paths), function(role) {
        p <- normalizePath(paths[[role]], mustWork = TRUE)
        b <- file.info(p)
        if (is.null(authenticated)) 
            sha <- digest::digest(file = p, algo = "sha256", serialize = FALSE)
        else {
            i <- match(role, authenticated$role)
            bi_assert(!is.na(i) && identical(p, normalizePath(authenticated$path[i], mustWork = TRUE)) && b$size == 
                authenticated$bytes[i] && abs(as.numeric(b$mtime) - authenticated$mtime[i]) < 1e-05 && grepl("^[0-9a-f]{64}$", 
                authenticated$sha256[i]), paste("AUTHENTICATED_SOURCE_CHANGED", role))
            sha <- as.character(authenticated$sha256[i])
        }
        a <- file.info(p)
        bi_assert(identical(b$size, a$size) && identical(b$mtime, a$mtime), paste("SOURCE_CHANGED_DURING_HASH", 
            role))
        data.frame(role = role, path = p, bytes = unname(a$size), mtime = as.numeric(a$mtime), sha256 = sha)
    }))
}

ip_source_pair <- function(code, dates, eids, label) {
    bi_isa_pair(code, dates)
    ix <- match(eids, as.character(code$eid))
    list(codes = code[ix, -1L, drop = FALSE], dates = dates[ix, -1L, drop = FALSE], linked = !is.na(ix), label = label)
}

ip_source_scan <- function(pair, entry, birth_year, allowed = NULL, pattern = NULL) {
    bi_assert(xor(is.null(allowed), is.null(pattern)), "EXACT_OR_ANCHORED_SCAN_REQUIRED")
    n <- length(entry)
    first <- rep(as.Date(NA), n)
    prior <- undated <- observed <- rep(FALSE, n)
    audit_as_of <- as.Date("2026-10-16")
    pseudo <- as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"))
    yr <- ip_source_num(birth_year)
    birth_floor <- rep(as.Date(NA), n)
    good <- is.finite(yr) & yr == floor(yr) & yr >= 1900 & yr <= 2100
    birth_floor[good] <- as.Date(paste0(as.integer(yr[good]), "-01-01"))
    bad_n <- post_n <- 0L
    for (j in seq_len(ncol(pair$codes))) {
        cc <- bi_isa_normalize(pair$codes[[j]])
        hit <- !is.na(cc) & if (is.null(pattern)) 
            cc %in% allowed
        else grepl(pattern, cc, perl = TRUE)
        if (!any(hit)) 
            next
        dd <- rep(as.Date(NA), n)
        dd[hit] <- bi_isa_dates(pair$dates[[j]][hit])
        bad <- hit & !is.na(dd) & (dd %in% pseudo | dd > audit_as_of | (!is.na(birth_floor) & dd < birth_floor))
        bad[is.na(bad)] <- FALSE
        bad_n <- bad_n + sum(bad)
        dd[bad] <- as.Date(NA)
        observed <- observed | hit
        undated <- undated | (hit & is.na(dd))
        prior <- prior | (hit & !is.na(dd) & !is.na(entry) & dd <= entry)
        post_n <- post_n + sum(hit & !is.na(dd) & !is.na(entry) & dd > entry)
        take <- hit & !is.na(dd) & (is.na(first) | dd < first)
        first[take] <- dd[take]
    }
    known <- pair$linked & !undated & !is.na(entry)
    flag <- as.integer(prior)
    flag[!known & !prior] <- NA_integer_
    list(first = first, prior = flag, observed = observed, undated = undated, invalid_date_count = bad_n, after_entry_pair_count = post_n, 
        known = known, audit_as_of = audit_as_of)
}

ip_history_union <- function(...) {
    z <- list(...)
    bi_assert(length(z) > 0L && length(unique(lengths(z))) == 1L, "HISTORY_LENGTH_MISMATCH")
    bi_assert(all(vapply(z, function(x) all(is.na(x) | x %in% 0:1), logical(1))), "INVALID_HISTORY_FLAG")
    positive <- Reduce(`|`, lapply(z, function(x) !is.na(x) & x == 1L))
    unknown <- Reduce(`|`, lapply(z, is.na))
    out <- as.integer(positive)
    out[unknown & !positive] <- NA_integer_
    out
}

ip_source_exposure_audit <- function(d) {
    ex <- c("border12", "border1n", "border_rank", "border_relative")
    covars <- c("age", "sex", "ethnicity", "center", "assessment_year")
    req <- c("eid", "cnt_older", "cnt_brother", "cnt_sister", "adopted_num", "eligible_border", "sibling_n", "sibship_size", 
        "birth_rank", "S", ex, covars)
    bi_assert(all(req %in% names(d)), "EXPOSURE_AUDIT_FIELDS_MISSING")
    bi_isa_ids(d$eid, "exposure audit")
    same <- function(x, y) identical(is.na(x), is.na(y)) && all(x == y, na.rm = TRUE)
    older <- ip_source_num(d$cnt_older)
    brother <- ip_source_num(d$cnt_brother)
    sister <- ip_source_num(d$cnt_sister)
    for (z in list(older, brother, sister)) bi_assert(all(is.na(z) | (is.finite(z) & z >= 0 & z == floor(z))), 
        "UNCLEAN_FROZEN_SIBLING_COUNT")
    adopted <- ip_source_num(d$adopted_num)
    bi_assert(all(is.na(adopted) | adopted %in% 0:1), "INVALID_FROZEN_ADOPTION_CODE")
    siblings <- brother + sister
    S <- siblings + 1
    rank <- older + 1
    eligible <- !is.na(older) & !is.na(brother) & !is.na(sister) & older <= siblings & (is.na(adopted) | adopted != 
        1)
    eligible[is.na(eligible)] <- FALSE
    bi_assert(identical(as.logical(d$eligible_border), eligible) && same(ip_source_num(d$sibling_n), siblings) && 
        same(ip_source_num(d$sibship_size), S) && same(ip_source_num(d$S), S) && same(ip_source_num(d$birth_rank), 
        rank), "FROZEN_EXPOSURE_DERIVATION_MISMATCH")
    expected <- setNames(rep(list(rep(NA_real_, nrow(d))), length(ex)), ex)
    two <- eligible & !is.na(S) & S == 2
    expected$border12[two] <- rank[two] - 1
    multiple <- eligible & !is.na(S) & S >= 2
    first <- multiple & !is.na(rank) & rank == 1
    last <- multiple & !is.na(rank) & rank == S
    expected$border1n[first] <- 0
    expected$border1n[last] <- 1
    expected$border_rank[multiple] <- rank[multiple]
    expected$border_relative[multiple] <- (rank[multiple] - 1)/(S[multiple] - 1)
    masks <- lapply(ex, function(k) {
        bi_assert(same(ip_source_num(d[[k]]), expected[[k]]), paste("FROZEN_EXPOSURE_VALUE_MISMATCH", k))
        !is.na(expected[[k]]) & is.finite(expected[[k]])
    })
    names(masks) <- ex
    selected <- Reduce(`|`, masks)
    bi_assert(identical(selected, multiple), "FROZEN_EXPOSURE_ELIGIBILITY_MISMATCH")
    bi_assert(any(selected), "NO_ELIGIBLE_FROZEN_EXPOSURE_ROWS")
    cc_masks <- list()
    rows <- list()
    for (k in ex) {
        variables <- c(covars, if (k != "border12") "S")
        cc <- stats::complete.cases(d[, variables, drop = FALSE])
        for (v in variables) if (is.numeric(d[[v]])) 
            cc <- cc & is.finite(d[[v]])
        cc_masks[[k]] <- masks[[k]] & cc
        rows[[k]] <- data.frame(exposure = k, input_count = nrow(d), raw_eligible_count = sum(masks[[k]]), covariate_complete_count = sum(cc_masks[[k]]), 
            covariate_incomplete_count = sum(masks[[k]] & !cc), eligibility_excluded_count = sum(!masks[[k]]), 
            covariates = paste(variables, collapse = ";"), stage = "Before endpoint baseline exclusions and censoring; not final model N", 
            stringsAsFactors = FALSE)
    }
    any_cc <- Reduce(`|`, cc_masks)
    rows$any_branch <- data.frame(exposure = "any_branch", input_count = nrow(d), raw_eligible_count = sum(selected), 
        covariate_complete_count = sum(any_cc), covariate_incomplete_count = sum(selected & !any_cc), eligibility_excluded_count = sum(!selected), 
        covariates = "Union of the four branch-specific covariate-complete sets", stage = "Before endpoint baseline exclusions and censoring; not final model N", 
        stringsAsFactors = FALSE)
    list(selected = selected, counts = bi_bind(rows))
}

ip_prepare <- function(ctx) {
    bi_assert(isTRUE(ctx$synthetic) || identical(Sys.info()[["sysname"]], "Linux"), "PARTICIPANT_ADAPTER_WINPC_WSL_ONLY")
    roles <- c("all", "phewas", "baseline", "operations", "icd10_code", "icd10_date", "icd9_code", "icd9_date", 
        "opcs4_code", "opcs4_date", "firstocc", "gp", "srd")
    bi_assert(all(roles %in% names(ctx$paths)) && all(file.exists(unlist(ctx$paths[roles]))), "REQUIRED_SOURCE_PATH_MISSING")
    bi_assert(isTRUE(ctx$synthetic) || !is.null(ctx$input_manifest), "AUTHENTICATED_INPUT_MANIFEST_REQUIRED")
    paths <- ctx$paths[roles]
    before <- ip_source_manifest(paths, ctx$input_manifest)
    ip_csv(ctx, "audit/source_adapter_manifest.csv", before)
    ip_log(ctx, "Authenticate full source rosters and exact six-PheCode mapping")
    all <- ip_source_read(paths$all, "frozen all")
    raw <- ip_source_read(paths$baseline, "full baseline extract")
    ops <- ip_source_read(paths$operations, "baseline operations extract")
    hc <- ip_source_read(paths$icd10_code, "ICD10 codes")
    hd <- ip_source_read(paths$icd10_date, "ICD10 dates")
    bi_isa_pair(hc, hd)
    bi_assert(setequal(raw$eid, hc$eid), "FULL_RAW_HOSPITAL_ROSTER_MISMATCH")
    bi_assert(setequal(ops$eid, raw$eid), "OPERATION_RAW_ROSTER_MISMATCH")
    bi_assert(all(hc$eid %in% all$eid), "HOSPITAL_PARTICIPANT_ABSENT_FROM_FROZEN_ALL")
    bi_assert(all(c("date_attend", "birth_year") %in% names(all)), "BASELINE_DATE_OR_BIRTH_YEAR_MISSING")
    all$date_attend <- bi_isa_dates(all$date_attend)
    all$birth_year <- ip_source_num(all$birth_year)
    bi_assert("p31_i0" %in% names(raw), "RAW_BASELINE_SEX_MISSING")
    ax <- match(raw$eid, all$eid)
    sx <- ip_source_num(all$sex[ax])
    sr <- ip_source_num(raw$p31_i0)
    bi_assert(!any(!is.na(sx) & !is.na(sr) & sx != sr), "RAW_FROZEN_SEX_ALIGNMENT_MISMATCH")
    maps <- new.env(parent = emptyenv())
    utils::data(list = c("phecode_map", "phecode_map_icd10", "phecode_rollup_map"), package = "PheWAS", envir = maps)
    cm_direct <- as.data.frame(maps$phecode_map)
    bi_assert(identical(bi_isa_canonical_hash(cm_direct), birth_private_setting("PRIVATE_HASH_1")), "ORIGINAL_ICD9CM_HISTORY_MAP_CHANGED")
    direct <- as.data.frame(maps$phecode_map_icd10)
    rollup <- as.data.frame(maps$phecode_rollup_map)
    mh <- bi_isa_canonical_hash(direct)
    rh <- bi_isa_canonical_hash(rollup)
    bi_assert(identical(mh, birth_private_setting("PRIVATE_HASH_2")) && identical(rh, birth_private_setting("PRIVATE_HASH_3")), 
        "PINNED_PHECODE_MAPPING_CHANGED")
    mapping <- sc_uk_map(direct, rollup)
    exact <- bi_isa_extract(hc, hd, mapping)
    flags <- ip_source_read(paths$phewas, "reconciliation-only PheWAS")
    bi_assert(all(bi_isa_codes() %in% names(flags)), "SIX_RECONCILIATION_FLAGS_MISSING")
    fi <- match(hc$eid, flags$eid)
    audit <- list()
    for (k in bi_isa_codes()) {
        old <- ip_source_num(flags[[k]][fi])
        bi_assert(all(is.na(old) | old %in% 0:1), "INVALID_RECONCILIATION_FLAG")
        comparable <- !is.na(old)
        delta <- sum(exact[[k]]$observed[comparable] != (old[comparable] == 1L))
        bi_assert(sum(comparable) > 0L, paste("NO_OLD_PHENOTYPE_COMPARISON", k))
        audit[[k]] <- data.frame(trait = k, compared_count = sum(comparable), disagreement_count = delta, hospital_positive_count = sum(exact[[k]]$observed), 
            hospital_undated_count = sum(exact[[k]]$undated), flag_missing_count = sum(!comparable), flag_missing_used_for_exclusion = FALSE, 
            exact_ICD10_code_count = sum(mapping$trait == k), direct_sha256 = mh, rollup_sha256 = rh, comparison = "New ICD10 map versus old US ICD10CM; differences reported, not forced to zero")
    }
    ip_csv(ctx, "audit/exact_flag_reconciliation.csv", bi_bind(audit))
    ip_csv(ctx, "audit/exact_icd10_mapping.csv", mapping)
    d <- linear_covars(make_birth_order(all))
    eligibility <- ip_source_exposure_audit(d)
    ip_csv(ctx, "audit/exposure_eligibility_counts.csv", eligibility$counts)
    selected <- eligibility$selected
    d <- d[selected, , drop = FALSE]
    d$eid <- as.character(d$eid)
    ids <- d$eid
    entry <- bi_isa_dates(d$date_attend)
    d$date_attend <- entry
    ci <- match(ids, hc$eid)
    ri <- match(ids, raw$eid)
    hospital_known <- !is.na(ci)
    raw_known <- !is.na(ri)
    ip_csv(ctx, "audit/source_roster_counts.csv", data.frame(source = c("frozen_all", "full_raw", "hospital", "operations", 
        "exposure_selected", "exposure_without_hospital", "exposure_without_raw"), count = c(nrow(all), nrow(raw), 
        nrow(hc), nrow(ops), nrow(d), sum(!hospital_known), sum(!raw_known))))
    rm(all, flags)
    registry <- ev_prov_registry_dates(raw, ids[raw_known])
    ix <- match(ids, registry$eid)
    d$date_death_uncapped <- registry$date_death_uncapped[ix]
    d$date_lost_registry <- registry$date_lost_registry[ix]
    d$registry_source_row_present <- raw_known
    d$hospital_coverage_group <- ifelse(hospital_known & raw_known, "UKB_institutional_uniform", NA_character_)
    p10 <- ip_source_pair(hc, hd, ids, "ICD10")
    events <- baseline_exact <- list()
    source_audit <- list()
    for (k in bi_isa_codes()) {
        z <- ip_source_scan(p10, entry, d$birth_year, allowed = mapping$normalized_code[mapping$trait == k])
        events[[k]] <- z$first
        baseline_exact[[k]] <- z$prior
        events[[k]][!z$known] <- as.Date(NA)
        source_audit[[paste0("ICD10_", k)]] <- data.frame(source = paste0("exact_ICD10_", k), prior_count = sum(z$prior == 
            1L, na.rm = TRUE), unknown_count = sum(is.na(z$prior)), invalid_date_count = z$invalid_date_count, 
            after_entry_pair_count = z$after_entry_pair_count)
    }
    ihd10 <- ip_source_scan(p10, entry, d$birth_year, pattern = "^I2[0-5][0-9A-Z]*$")
    dm10 <- ip_source_scan(p10, entry, d$birth_year, pattern = "^E1[0-4][0-9A-Z]*$")
    rm(hc, hd, exact, p10)
    ip_log(ctx, "Construct reviewed coronary, diabetes and exact-351 baseline history")
    c9 <- ip_source_read(paths$icd9_code, "ICD9 codes")
    d9 <- ip_source_read(paths$icd9_date, "ICD9 dates")
    p9 <- ip_source_pair(c9, d9, ids, "ICD9")
    rolled <- merge(cm_direct, rollup, by.x = "phecode", by.y = "code", sort = FALSE)
    m9 <- unique(rolled[rolled$vocabulary_id == "ICD9CM" & rolled$phecode_unrolled %in% c("250", "351"), c("code", 
        "phecode_unrolled")])
    m9$normalized_code <- bi_isa_normalize(m9$code)
    bi_assert(sum(m9$phecode_unrolled == "250") == 66L && sum(m9$phecode_unrolled == "351") == 18L, "EXACT_ICD9_HISTORY_MAPPING_CHANGED")
    ip_csv(ctx, "audit/exact_icd9_baseline_mapping.csv", m9)
    ihd9 <- ip_source_scan(p9, entry, d$birth_year, pattern = "^41[0-4][0-9A-Z]*$")
    dm9 <- ip_source_scan(p9, entry, d$birth_year, allowed = m9$normalized_code[m9$phecode_unrolled == "250"])
    neuro9 <- ip_source_scan(p9, entry, d$birth_year, allowed = m9$normalized_code[m9$phecode_unrolled == "351"])
    rm(c9, d9, p9)
    oc <- ip_source_read(paths$opcs4_code, "OPCS codes")
    od <- ip_source_read(paths$opcs4_date, "OPCS dates")
    surgery <- ip_source_scan(ip_source_pair(oc, od, ids, "OPCS4"), entry, d$birth_year, allowed = ev_prov_opcs_coronary_codes())
    rm(oc, od)
    srd <- ev_inc_srd(list(code_path = paths$srd, baseline_positions = paste0("i", 1:37)), ids)
    srd_flag <- as.integer(srd$ihd)
    srd_flag[!srd$known & !srd$ihd] <- NA_integer_
    firstocc <- ev_prov_firstocc_prior(ip_source_read(paths$firstocc, "first occurrence"), ids, entry)
    gp <- ev_prov_gp_positive_prior(ip_source_read(paths$gp, "partial GP"), ids, entry)
    ts <- ev_prov_touchscreen_prior(ops, ids)
    ev_ids <- ip_source_ids
    ev_num <- ip_source_num
    opfun <- ev_operation_prior
    environment(opfun) <- environment()
    op <- opfun(ops, ids)
    dmcontract <- list(verified = TRUE, baseline_only_verified = TRUE, arrays_verified = TRUE, medication_dictionary_verified = TRUE, 
        evidence = "Pinned conservative baseline adaptation; authenticated full baseline extract and runtime gates", 
        algorithm = "Eastwood2016_Fig2_S1_conservative_baseline_v1")
    dm <- ev_diabetes_baseline(raw, ids, dmcontract)
    bi_assert(isTRUE(dm$verified) && identical(as.character(dm$eid), ids), "DIABETES_BASELINE_GATE_FAILED")
    ip_csv(ctx, "audit/diabetes_source_gates.csv", dm$gates)
    ip_csv(ctx, "audit/diabetes_baseline_counts.csv", dm$audit)
    dpcontract <- list(verified = TRUE, baseline_comparison_verified = TRUE, evidence = "Pinned positive-only GP SRD and first-occurrence definitions; uncertainty excluded", 
        source_columns = ev_diabetes_prior_specs()$column)
    dp <- ev_diabetes_prior_sources(d, ids, entry, birth_year = d$birth_year, contract = dpcontract)
    bi_assert(isTRUE(dp$verified), "DIABETES_DATED_PRIOR_GATE_FAILED")
    ip_csv(ctx, "audit/diabetes_dated_prior_counts.csv", dp$audit)
    dmflag <- ifelse(dm$status == "no_diabetes", 0L, ifelse(dm$status %in% c("type1", "type2"), 1L, NA_integer_))
    dpflag <- as.integer(dp$prior)
    dpflag[dp$timing_uncertain & !dp$prior] <- NA_integer_
    coronary <- ip_history_union(ihd10$prior, ihd9$prior, srd_flag, surgery$prior, firstocc$prior_ihd_ukb_firstocc, 
        gp$prior_ihd_gp_available_positive, ts$prior_ihd_touchscreen, op$prior_ihd_operation_report)
    diabetes <- ip_history_union(baseline_exact[["250"]], dm10$prior, dm9$prior, dmflag, dpflag)
    baseline <- list()
    for (k in c("411", "411.3", "411.4")) baseline[[k]] <- ip_history_union(baseline_exact[[k]], coronary)
    for (k in c("250", "250.2")) baseline[[k]] <- ip_history_union(baseline_exact[[k]], diabetes)
    baseline[["351"]] <- ip_history_union(baseline_exact[["351"]], neuro9$prior)
    baseline <- baseline[bi_isa_codes()]
    events <- events[bi_isa_codes()]
    for (k in names(baseline)) baseline[[k]][!hospital_known | !raw_known] <- NA_integer_
    extra <- list(broad_ICD10_IHD = ihd10, broad_ICD10_DM = dm10, ICD9_IHD = ihd9, exact_ICD9_DM = dm9, exact_ICD9_351 = neuro9, 
        OPCS_coronary = surgery)
    for (k in names(extra)) {
        z <- extra[[k]]
        source_audit[[k]] <- data.frame(source = k, prior_count = sum(z$prior == 1L, na.rm = TRUE), unknown_count = sum(is.na(z$prior)), 
            invalid_date_count = z$invalid_date_count, after_entry_pair_count = z$after_entry_pair_count)
    }
    supplementary <- list(baseline_SRD_IHD = srd_flag, firstocc_IHD = firstocc$prior_ihd_ukb_firstocc, GP_available_positive_IHD = gp$prior_ihd_gp_available_positive, 
        touchscreen_IHD = ts$prior_ihd_touchscreen, baseline_operation_IHD = op$prior_ihd_operation_report, conservative_DM_algorithm = dmflag, 
        available_dated_DM_history = dpflag)
    for (k in names(supplementary)) source_audit[[k]] <- data.frame(source = k, prior_count = sum(supplementary[[k]] == 
        1L, na.rm = TRUE), unknown_count = sum(is.na(supplementary[[k]])), invalid_date_count = NA_integer_, after_entry_pair_count = NA_integer_)
    ip_csv(ctx, "audit/baseline_component_counts.csv", bi_bind(source_audit))
    ip_csv(ctx, "audit/endpoint_adapter_counts.csv", bi_bind(lapply(names(baseline), function(k) data.frame(trait = k, 
        exposure_population_count = nrow(d), baseline_positive_count = sum(baseline[[k]] == 1L, na.rm = TRUE), 
        baseline_unknown_count = sum(is.na(baseline[[k]])), baseline_unidentified_count = sum(baseline[[k]] == 
            0L, na.rm = TRUE), dated_hospital_event_count = sum(!is.na(events[[k]])), baseline_scope = if (k == 
            "351") 
            "Exact mapped ICD10 and ICD9 history only; no unverified SRD GP or broad neuro substitution"
        else if (k %in% c("250", "250.2")) 
            "Conservative baseline DM algorithm plus mapped and broad hospital DM and available dated history"
        else "Exact target plus broad IHD ICD10 ICD9 SRD OPCS touchscreen operations firstocc and available GP", 
        estimand = "First recorded hospital event among participants without identified baseline disease; not complete clinical incidence"))))
    after <- ip_source_manifest(paths, before)
    bi_assert(identical(before, after), "SOURCE_CHANGED_DURING_ADAPTER")
    ip_csv(ctx, "audit/source_adapter_after_manifest.csv", after)
    bi_assert(all(vapply(events, function(x) inherits(x, "Date") && length(x) == nrow(d), logical(1))) && all(vapply(baseline, 
        function(x) length(x) == nrow(d) && all(is.na(x) | x %in% 0:1), logical(1))), "ADAPTER_OUTPUT_SCHEMA_INVALID")
    ip_log(ctx, paste("Source adapter complete; eligible exposure rows", nrow(d), "six exact hospital endpoints"))
    list(data = d, events = events, baseline = baseline)
}

sc_uk_map <- function(direct, rollup) {
    stopifnot(identical(names(direct), c("vocabulary_id", "code", "phecode")), all(direct$vocabulary_id == "ICD10"), 
        identical(bi_isa_canonical_hash(direct), birth_private_setting("PRIVATE_HASH_4")))
    z <- merge(direct, rollup, by.x = "phecode", by.y = "code", sort = FALSE)
    z <- z[z$phecode_unrolled %in% bi_isa_codes(), ]
    ans <- unique(data.frame(trait = as.character(z$phecode_unrolled), source_code = as.character(z$code), normalized_code = bi_isa_normalize(z$code)))
    stopifnot(setequal(unique(ans$trait), bi_isa_codes()))
    expected <- c(`250` = 37L, `250.2` = 23L, `351` = 23L, `411` = 43L, `411.3` = 4L, `411.4` = 7L)
    stopifnot(all(table(ans$trait)[names(expected)] == expected))
    ans
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        ip_source_ids <- function(x, label = "source") {
            x <- as.data.frame(x)
            key <- intersect(c("eid", "id"), names(x))
            bi_assert(length(key) == 1L, paste("SOURCE_ID_COLUMN_INVALID", label))
            names(x)[names(x) == key] <- "eid"
            x$eid <- bi_isa_ids(x$eid, label)
            bi_assert(!anyDuplicated(names(x)), paste("SOURCE_DUPLICATE_COLUMNS", label))
            x
        }
        ip_source_num <- function(x) if (is.numeric(x) || is.logical(x)) as.numeric(x) else bi_num(x)
        ip_source_read <- function(path, label) {
            x <- if (grepl("[.]rds$", path, ignore.case = TRUE)) readRDS(path) else data.table::fread(path, data.table = FALSE, 
                check.names = FALSE, showProgress = FALSE)
            ip_source_ids(x, label)
        }
        ip_source_manifest <- function(paths, authenticated = NULL) {
            if (!is.null(authenticated)) {
                authenticated <- as.data.frame(authenticated)
                if ("key" %in% names(authenticated) && !"role" %in% names(authenticated)) names(authenticated)[names(authenticated) == 
                  "key"] <- "role"
                bi_assert(all(c("role", "path", "sha256", "bytes", "mtime") %in% names(authenticated)) && !anyDuplicated(authenticated$role), 
                  "AUTHENTICATED_SOURCE_MANIFEST_INVALID")
            }
            bi_bind(lapply(names(paths), function(role) {
                p <- normalizePath(paths[[role]], mustWork = TRUE)
                b <- file.info(p)
                if (is.null(authenticated)) sha <- digest::digest(file = p, algo = "sha256", serialize = FALSE) else {
                  i <- match(role, authenticated$role)
                  bi_assert(!is.na(i) && identical(p, normalizePath(authenticated$path[i], mustWork = TRUE)) && 
                    b$size == authenticated$bytes[i] && abs(as.numeric(b$mtime) - authenticated$mtime[i]) < 1e-05 && 
                    grepl("^[0-9a-f]{64}$", authenticated$sha256[i]), paste("AUTHENTICATED_SOURCE_CHANGED", role))
                  sha <- as.character(authenticated$sha256[i])
                }
                a <- file.info(p)
                bi_assert(identical(b$size, a$size) && identical(b$mtime, a$mtime), paste("SOURCE_CHANGED_DURING_HASH", 
                  role))
                data.frame(role = role, path = p, bytes = unname(a$size), mtime = as.numeric(a$mtime), sha256 = sha)
            }))
        }
        ip_source_pair <- function(code, dates, eids, label) {
            bi_isa_pair(code, dates)
            ix <- match(eids, as.character(code$eid))
            list(codes = code[ix, -1L, drop = FALSE], dates = dates[ix, -1L, drop = FALSE], linked = !is.na(ix), 
                label = label)
        }
        ip_source_scan <- function(pair, entry, birth_year, allowed = NULL, pattern = NULL) {
            bi_assert(xor(is.null(allowed), is.null(pattern)), "EXACT_OR_ANCHORED_SCAN_REQUIRED")
            n <- length(entry)
            first <- rep(as.Date(NA), n)
            prior <- undated <- observed <- rep(FALSE, n)
            audit_as_of <- as.Date("2026-10-16")
            pseudo <- as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"))
            yr <- ip_source_num(birth_year)
            birth_floor <- rep(as.Date(NA), n)
            good <- is.finite(yr) & yr == floor(yr) & yr >= 1900 & yr <= 2100
            birth_floor[good] <- as.Date(paste0(as.integer(yr[good]), "-01-01"))
            bad_n <- post_n <- 0L
            for (j in seq_len(ncol(pair$codes))) {
                cc <- bi_isa_normalize(pair$codes[[j]])
                hit <- !is.na(cc) & if (is.null(pattern)) cc %in% allowed else grepl(pattern, cc, perl = TRUE)
                if (!any(hit)) next
                dd <- rep(as.Date(NA), n)
                dd[hit] <- bi_isa_dates(pair$dates[[j]][hit])
                bad <- hit & !is.na(dd) & (dd %in% pseudo | dd > audit_as_of | (!is.na(birth_floor) & dd < birth_floor))
                bad[is.na(bad)] <- FALSE
                bad_n <- bad_n + sum(bad)
                dd[bad] <- as.Date(NA)
                observed <- observed | hit
                undated <- undated | (hit & is.na(dd))
                prior <- prior | (hit & !is.na(dd) & !is.na(entry) & dd <= entry)
                post_n <- post_n + sum(hit & !is.na(dd) & !is.na(entry) & dd > entry)
                take <- hit & !is.na(dd) & (is.na(first) | dd < first)
                first[take] <- dd[take]
            }
            known <- pair$linked & !undated & !is.na(entry)
            flag <- as.integer(prior)
            flag[!known & !prior] <- NA_integer_
            list(first = first, prior = flag, observed = observed, undated = undated, invalid_date_count = bad_n, 
                after_entry_pair_count = post_n, known = known, audit_as_of = audit_as_of)
        }
        ip_history_union <- function(...) {
            z <- list(...)
            bi_assert(length(z) > 0L && length(unique(lengths(z))) == 1L, "HISTORY_LENGTH_MISMATCH")
            bi_assert(all(vapply(z, function(x) all(is.na(x) | x %in% 0:1), logical(1))), "INVALID_HISTORY_FLAG")
            positive <- Reduce(`|`, lapply(z, function(x) !is.na(x) & x == 1L))
            unknown <- Reduce(`|`, lapply(z, is.na))
            out <- as.integer(positive)
            out[unknown & !positive] <- NA_integer_
            out
        }
        ip_source_exposure_audit <- function(d) {
            ex <- c("border12", "border1n", "border_rank", "border_relative")
            covars <- c("age", "sex", "ethnicity", "center", "assessment_year")
            req <- c("eid", "cnt_older", "cnt_brother", "cnt_sister", "adopted_num", "eligible_border", "sibling_n", 
                "sibship_size", "birth_rank", "S", ex, covars)
            bi_assert(all(req %in% names(d)), "EXPOSURE_AUDIT_FIELDS_MISSING")
            bi_isa_ids(d$eid, "exposure audit")
            same <- function(x, y) identical(is.na(x), is.na(y)) && all(x == y, na.rm = TRUE)
            older <- ip_source_num(d$cnt_older)
            brother <- ip_source_num(d$cnt_brother)
            sister <- ip_source_num(d$cnt_sister)
            for (z in list(older, brother, sister)) bi_assert(all(is.na(z) | (is.finite(z) & z >= 0 & z == floor(z))), 
                "UNCLEAN_FROZEN_SIBLING_COUNT")
            adopted <- ip_source_num(d$adopted_num)
            bi_assert(all(is.na(adopted) | adopted %in% 0:1), "INVALID_FROZEN_ADOPTION_CODE")
            siblings <- brother + sister
            S <- siblings + 1
            rank <- older + 1
            eligible <- !is.na(older) & !is.na(brother) & !is.na(sister) & older <= siblings & (is.na(adopted) | 
                adopted != 1)
            eligible[is.na(eligible)] <- FALSE
            bi_assert(identical(as.logical(d$eligible_border), eligible) && same(ip_source_num(d$sibling_n), siblings) && 
                same(ip_source_num(d$sibship_size), S) && same(ip_source_num(d$S), S) && same(ip_source_num(d$birth_rank), 
                rank), "FROZEN_EXPOSURE_DERIVATION_MISMATCH")
            expected <- setNames(rep(list(rep(NA_real_, nrow(d))), length(ex)), ex)
            two <- eligible & !is.na(S) & S == 2
            expected$border12[two] <- rank[two] - 1
            multiple <- eligible & !is.na(S) & S >= 2
            first <- multiple & !is.na(rank) & rank == 1
            last <- multiple & !is.na(rank) & rank == S
            expected$border1n[first] <- 0
            expected$border1n[last] <- 1
            expected$border_rank[multiple] <- rank[multiple]
            expected$border_relative[multiple] <- (rank[multiple] - 1)/(S[multiple] - 1)
            masks <- lapply(ex, function(k) {
                bi_assert(same(ip_source_num(d[[k]]), expected[[k]]), paste("FROZEN_EXPOSURE_VALUE_MISMATCH", k))
                !is.na(expected[[k]]) & is.finite(expected[[k]])
            })
            names(masks) <- ex
            selected <- Reduce(`|`, masks)
            bi_assert(identical(selected, multiple), "FROZEN_EXPOSURE_ELIGIBILITY_MISMATCH")
            bi_assert(any(selected), "NO_ELIGIBLE_FROZEN_EXPOSURE_ROWS")
            cc_masks <- list()
            rows <- list()
            for (k in ex) {
                variables <- c(covars, if (k != "border12") "S")
                cc <- stats::complete.cases(d[, variables, drop = FALSE])
                for (v in variables) if (is.numeric(d[[v]])) cc <- cc & is.finite(d[[v]])
                cc_masks[[k]] <- masks[[k]] & cc
                rows[[k]] <- data.frame(exposure = k, input_count = nrow(d), raw_eligible_count = sum(masks[[k]]), 
                  covariate_complete_count = sum(cc_masks[[k]]), covariate_incomplete_count = sum(masks[[k]] & 
                    !cc), eligibility_excluded_count = sum(!masks[[k]]), covariates = paste(variables, collapse = ";"), 
                  stage = "Before endpoint baseline exclusions and censoring; not final model N", stringsAsFactors = FALSE)
            }
            any_cc <- Reduce(`|`, cc_masks)
            rows$any_branch <- data.frame(exposure = "any_branch", input_count = nrow(d), raw_eligible_count = sum(selected), 
                covariate_complete_count = sum(any_cc), covariate_incomplete_count = sum(selected & !any_cc), eligibility_excluded_count = sum(!selected), 
                covariates = "Union of the four branch-specific covariate-complete sets", stage = "Before endpoint baseline exclusions and censoring; not final model N", 
                stringsAsFactors = FALSE)
            list(selected = selected, counts = bi_bind(rows))
        }
        ip_prepare <- function(ctx) {
            bi_assert(isTRUE(ctx$synthetic) || identical(Sys.info()[["sysname"]], "Linux"), "PARTICIPANT_ADAPTER_WINPC_WSL_ONLY")
            roles <- c("all", "phewas", "baseline", "operations", "icd10_code", "icd10_date", "icd9_code", "icd9_date", 
                "opcs4_code", "opcs4_date", "firstocc", "gp", "srd")
            bi_assert(all(roles %in% names(ctx$paths)) && all(file.exists(unlist(ctx$paths[roles]))), "REQUIRED_SOURCE_PATH_MISSING")
            bi_assert(isTRUE(ctx$synthetic) || !is.null(ctx$input_manifest), "AUTHENTICATED_INPUT_MANIFEST_REQUIRED")
            paths <- ctx$paths[roles]
            before <- ip_source_manifest(paths, ctx$input_manifest)
            ip_csv(ctx, "audit/source_adapter_manifest.csv", before)
            ip_log(ctx, "Authenticate full source rosters and exact six-PheCode mapping")
            all <- ip_source_read(paths$all, "frozen all")
            raw <- ip_source_read(paths$baseline, "full baseline extract")
            ops <- ip_source_read(paths$operations, "baseline operations extract")
            hc <- ip_source_read(paths$icd10_code, "ICD10 codes")
            hd <- ip_source_read(paths$icd10_date, "ICD10 dates")
            bi_isa_pair(hc, hd)
            bi_assert(setequal(raw$eid, hc$eid), "FULL_RAW_HOSPITAL_ROSTER_MISMATCH")
            bi_assert(setequal(ops$eid, raw$eid), "OPERATION_RAW_ROSTER_MISMATCH")
            bi_assert(all(hc$eid %in% all$eid), "HOSPITAL_PARTICIPANT_ABSENT_FROM_FROZEN_ALL")
            bi_assert(all(c("date_attend", "birth_year") %in% names(all)), "BASELINE_DATE_OR_BIRTH_YEAR_MISSING")
            all$date_attend <- bi_isa_dates(all$date_attend)
            all$birth_year <- ip_source_num(all$birth_year)
            bi_assert("p31_i0" %in% names(raw), "RAW_BASELINE_SEX_MISSING")
            ax <- match(raw$eid, all$eid)
            sx <- ip_source_num(all$sex[ax])
            sr <- ip_source_num(raw$p31_i0)
            bi_assert(!any(!is.na(sx) & !is.na(sr) & sx != sr), "RAW_FROZEN_SEX_ALIGNMENT_MISMATCH")
            maps <- new.env(parent = emptyenv())
            utils::data(list = c("phecode_map", "phecode_map_icd10", "phecode_rollup_map"), package = "PheWAS", 
                envir = maps)
            cm_direct <- as.data.frame(maps$phecode_map)
            bi_assert(identical(bi_isa_canonical_hash(cm_direct), birth_private_setting("PRIVATE_HASH_5")), "ORIGINAL_ICD9CM_HISTORY_MAP_CHANGED")
            direct <- as.data.frame(maps$phecode_map_icd10)
            rollup <- as.data.frame(maps$phecode_rollup_map)
            mh <- bi_isa_canonical_hash(direct)
            rh <- bi_isa_canonical_hash(rollup)
            bi_assert(identical(mh, birth_private_setting("PRIVATE_HASH_6")) && identical(rh, birth_private_setting("PRIVATE_HASH_7")), 
                "PINNED_PHECODE_MAPPING_CHANGED")
            mapping <- sc_uk_map(direct, rollup)
            exact <- bi_isa_extract(hc, hd, mapping)
            flags <- ip_source_read(paths$phewas, "reconciliation-only PheWAS")
            bi_assert(all(bi_isa_codes() %in% names(flags)), "SIX_RECONCILIATION_FLAGS_MISSING")
            fi <- match(hc$eid, flags$eid)
            audit <- list()
            for (k in bi_isa_codes()) {
                old <- ip_source_num(flags[[k]][fi])
                bi_assert(all(is.na(old) | old %in% 0:1), "INVALID_RECONCILIATION_FLAG")
                comparable <- !is.na(old)
                delta <- sum(exact[[k]]$observed[comparable] != (old[comparable] == 1L))
                bi_assert(sum(comparable) > 0L, paste("NO_OLD_PHENOTYPE_COMPARISON", k))
                audit[[k]] <- data.frame(trait = k, compared_count = sum(comparable), disagreement_count = delta, 
                  hospital_positive_count = sum(exact[[k]]$observed), hospital_undated_count = sum(exact[[k]]$undated), 
                  flag_missing_count = sum(!comparable), flag_missing_used_for_exclusion = FALSE, exact_ICD10_code_count = sum(mapping$trait == 
                    k), direct_sha256 = mh, rollup_sha256 = rh, comparison = "New ICD10 map versus old US ICD10CM; differences reported, not forced to zero")
            }
            ip_csv(ctx, "audit/exact_flag_reconciliation.csv", bi_bind(audit))
            ip_csv(ctx, "audit/exact_icd10_mapping.csv", mapping)
            d <- linear_covars(make_birth_order(all))
            eligibility <- ip_source_exposure_audit(d)
            ip_csv(ctx, "audit/exposure_eligibility_counts.csv", eligibility$counts)
            selected <- eligibility$selected
            d <- d[selected, , drop = FALSE]
            d$eid <- as.character(d$eid)
            ids <- d$eid
            entry <- bi_isa_dates(d$date_attend)
            d$date_attend <- entry
            ci <- match(ids, hc$eid)
            ri <- match(ids, raw$eid)
            hospital_known <- !is.na(ci)
            raw_known <- !is.na(ri)
            ip_csv(ctx, "audit/source_roster_counts.csv", data.frame(source = c("frozen_all", "full_raw", "hospital", 
                "operations", "exposure_selected", "exposure_without_hospital", "exposure_without_raw"), count = c(nrow(all), 
                nrow(raw), nrow(hc), nrow(ops), nrow(d), sum(!hospital_known), sum(!raw_known))))
            rm(all, flags)
            registry <- ev_prov_registry_dates(raw, ids[raw_known])
            ix <- match(ids, registry$eid)
            d$date_death_uncapped <- registry$date_death_uncapped[ix]
            d$date_lost_registry <- registry$date_lost_registry[ix]
            d$registry_source_row_present <- raw_known
            d$hospital_coverage_group <- ifelse(hospital_known & raw_known, "UKB_institutional_uniform", NA_character_)
            p10 <- ip_source_pair(hc, hd, ids, "ICD10")
            events <- baseline_exact <- list()
            source_audit <- list()
            for (k in bi_isa_codes()) {
                z <- ip_source_scan(p10, entry, d$birth_year, allowed = mapping$normalized_code[mapping$trait == 
                  k])
                events[[k]] <- z$first
                baseline_exact[[k]] <- z$prior
                events[[k]][!z$known] <- as.Date(NA)
                source_audit[[paste0("ICD10_", k)]] <- data.frame(source = paste0("exact_ICD10_", k), prior_count = sum(z$prior == 
                  1L, na.rm = TRUE), unknown_count = sum(is.na(z$prior)), invalid_date_count = z$invalid_date_count, 
                  after_entry_pair_count = z$after_entry_pair_count)
            }
            ihd10 <- ip_source_scan(p10, entry, d$birth_year, pattern = "^I2[0-5][0-9A-Z]*$")
            dm10 <- ip_source_scan(p10, entry, d$birth_year, pattern = "^E1[0-4][0-9A-Z]*$")
            rm(hc, hd, exact, p10)
            ip_log(ctx, "Construct reviewed coronary, diabetes and exact-351 baseline history")
            c9 <- ip_source_read(paths$icd9_code, "ICD9 codes")
            d9 <- ip_source_read(paths$icd9_date, "ICD9 dates")
            p9 <- ip_source_pair(c9, d9, ids, "ICD9")
            rolled <- merge(cm_direct, rollup, by.x = "phecode", by.y = "code", sort = FALSE)
            m9 <- unique(rolled[rolled$vocabulary_id == "ICD9CM" & rolled$phecode_unrolled %in% c("250", "351"), 
                c("code", "phecode_unrolled")])
            m9$normalized_code <- bi_isa_normalize(m9$code)
            bi_assert(sum(m9$phecode_unrolled == "250") == 66L && sum(m9$phecode_unrolled == "351") == 18L, "EXACT_ICD9_HISTORY_MAPPING_CHANGED")
            ip_csv(ctx, "audit/exact_icd9_baseline_mapping.csv", m9)
            ihd9 <- ip_source_scan(p9, entry, d$birth_year, pattern = "^41[0-4][0-9A-Z]*$")
            dm9 <- ip_source_scan(p9, entry, d$birth_year, allowed = m9$normalized_code[m9$phecode_unrolled == 
                "250"])
            neuro9 <- ip_source_scan(p9, entry, d$birth_year, allowed = m9$normalized_code[m9$phecode_unrolled == 
                "351"])
            rm(c9, d9, p9)
            oc <- ip_source_read(paths$opcs4_code, "OPCS codes")
            od <- ip_source_read(paths$opcs4_date, "OPCS dates")
            surgery <- ip_source_scan(ip_source_pair(oc, od, ids, "OPCS4"), entry, d$birth_year, allowed = ev_prov_opcs_coronary_codes())
            rm(oc, od)
            srd <- ev_inc_srd(list(code_path = paths$srd, baseline_positions = paste0("i", 1:37)), ids)
            srd_flag <- as.integer(srd$ihd)
            srd_flag[!srd$known & !srd$ihd] <- NA_integer_
            firstocc <- ev_prov_firstocc_prior(ip_source_read(paths$firstocc, "first occurrence"), ids, entry)
            gp <- ev_prov_gp_positive_prior(ip_source_read(paths$gp, "partial GP"), ids, entry)
            ts <- ev_prov_touchscreen_prior(ops, ids)
            ev_ids <- ip_source_ids
            ev_num <- ip_source_num
            opfun <- ev_operation_prior
            environment(opfun) <- environment()
            op <- opfun(ops, ids)
            dmcontract <- list(verified = TRUE, baseline_only_verified = TRUE, arrays_verified = TRUE, medication_dictionary_verified = TRUE, 
                evidence = "Pinned conservative baseline adaptation; authenticated full baseline extract and runtime gates", 
                algorithm = "Eastwood2016_Fig2_S1_conservative_baseline_v1")
            dm <- ev_diabetes_baseline(raw, ids, dmcontract)
            bi_assert(isTRUE(dm$verified) && identical(as.character(dm$eid), ids), "DIABETES_BASELINE_GATE_FAILED")
            ip_csv(ctx, "audit/diabetes_source_gates.csv", dm$gates)
            ip_csv(ctx, "audit/diabetes_baseline_counts.csv", dm$audit)
            dpcontract <- list(verified = TRUE, baseline_comparison_verified = TRUE, evidence = "Pinned positive-only GP SRD and first-occurrence definitions; uncertainty excluded", 
                source_columns = ev_diabetes_prior_specs()$column)
            dp <- ev_diabetes_prior_sources(d, ids, entry, birth_year = d$birth_year, contract = dpcontract)
            bi_assert(isTRUE(dp$verified), "DIABETES_DATED_PRIOR_GATE_FAILED")
            ip_csv(ctx, "audit/diabetes_dated_prior_counts.csv", dp$audit)
            dmflag <- ifelse(dm$status == "no_diabetes", 0L, ifelse(dm$status %in% c("type1", "type2"), 1L, NA_integer_))
            dpflag <- as.integer(dp$prior)
            dpflag[dp$timing_uncertain & !dp$prior] <- NA_integer_
            coronary <- ip_history_union(ihd10$prior, ihd9$prior, srd_flag, surgery$prior, firstocc$prior_ihd_ukb_firstocc, 
                gp$prior_ihd_gp_available_positive, ts$prior_ihd_touchscreen, op$prior_ihd_operation_report)
            diabetes <- ip_history_union(baseline_exact[["250"]], dm10$prior, dm9$prior, dmflag, dpflag)
            baseline <- list()
            for (k in c("411", "411.3", "411.4")) baseline[[k]] <- ip_history_union(baseline_exact[[k]], coronary)
            for (k in c("250", "250.2")) baseline[[k]] <- ip_history_union(baseline_exact[[k]], diabetes)
            baseline[["351"]] <- ip_history_union(baseline_exact[["351"]], neuro9$prior)
            baseline <- baseline[bi_isa_codes()]
            events <- events[bi_isa_codes()]
            for (k in names(baseline)) baseline[[k]][!hospital_known | !raw_known] <- NA_integer_
            extra <- list(broad_ICD10_IHD = ihd10, broad_ICD10_DM = dm10, ICD9_IHD = ihd9, exact_ICD9_DM = dm9, 
                exact_ICD9_351 = neuro9, OPCS_coronary = surgery)
            for (k in names(extra)) {
                z <- extra[[k]]
                source_audit[[k]] <- data.frame(source = k, prior_count = sum(z$prior == 1L, na.rm = TRUE), unknown_count = sum(is.na(z$prior)), 
                  invalid_date_count = z$invalid_date_count, after_entry_pair_count = z$after_entry_pair_count)
            }
            supplementary <- list(baseline_SRD_IHD = srd_flag, firstocc_IHD = firstocc$prior_ihd_ukb_firstocc, 
                GP_available_positive_IHD = gp$prior_ihd_gp_available_positive, touchscreen_IHD = ts$prior_ihd_touchscreen, 
                baseline_operation_IHD = op$prior_ihd_operation_report, conservative_DM_algorithm = dmflag, available_dated_DM_history = dpflag)
            for (k in names(supplementary)) source_audit[[k]] <- data.frame(source = k, prior_count = sum(supplementary[[k]] == 
                1L, na.rm = TRUE), unknown_count = sum(is.na(supplementary[[k]])), invalid_date_count = NA_integer_, 
                after_entry_pair_count = NA_integer_)
            ip_csv(ctx, "audit/baseline_component_counts.csv", bi_bind(source_audit))
            ip_csv(ctx, "audit/endpoint_adapter_counts.csv", bi_bind(lapply(names(baseline), function(k) data.frame(trait = k, 
                exposure_population_count = nrow(d), baseline_positive_count = sum(baseline[[k]] == 1L, na.rm = TRUE), 
                baseline_unknown_count = sum(is.na(baseline[[k]])), baseline_unidentified_count = sum(baseline[[k]] == 
                  0L, na.rm = TRUE), dated_hospital_event_count = sum(!is.na(events[[k]])), baseline_scope = if (k == 
                  "351") "Exact mapped ICD10 and ICD9 history only; no unverified SRD GP or broad neuro substitution" else if (k %in% 
                  c("250", "250.2")) "Conservative baseline DM algorithm plus mapped and broad hospital DM and available dated history" else "Exact target plus broad IHD ICD10 ICD9 SRD OPCS touchscreen operations firstocc and available GP", 
                estimand = "First recorded hospital event among participants without identified baseline disease; not complete clinical incidence"))))
            after <- ip_source_manifest(paths, before)
            bi_assert(identical(before, after), "SOURCE_CHANGED_DURING_ADAPTER")
            ip_csv(ctx, "audit/source_adapter_after_manifest.csv", after)
            bi_assert(all(vapply(events, function(x) inherits(x, "Date") && length(x) == nrow(d), logical(1))) && 
                all(vapply(baseline, function(x) length(x) == nrow(d) && all(is.na(x) | x %in% 0:1), logical(1))), 
                "ADAPTER_OUTPUT_SCHEMA_INVALID")
            ip_log(ctx, paste("Source adapter complete; eligible exposure rows", nrow(d), "six exact hospital endpoints"))
            list(data = d, events = events, baseline = baseline)
        }
        sc_uk_map <- function(direct, rollup) {
            stopifnot(identical(names(direct), c("vocabulary_id", "code", "phecode")), all(direct$vocabulary_id == 
                "ICD10"), identical(bi_isa_canonical_hash(direct), birth_private_setting("PRIVATE_HASH_8")))
            z <- merge(direct, rollup, by.x = "phecode", by.y = "code", sort = FALSE)
            z <- z[z$phecode_unrolled %in% bi_isa_codes(), ]
            ans <- unique(data.frame(trait = as.character(z$phecode_unrolled), source_code = as.character(z$code), 
                normalized_code = bi_isa_normalize(z$code)))
            stopifnot(setequal(unique(ans$trait), bi_isa_codes()))
            expected <- c(`250` = 37L, `250.2` = 23L, `351` = 23L, `411` = 43L, `411.3` = 4L, `411.4` = 7L)
            stopifnot(all(table(ans$trait)[names(expected)] == expected))
            ans
        }
    }), envir = context)
}

