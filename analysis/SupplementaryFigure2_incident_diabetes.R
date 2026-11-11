# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

ev_dm_contract_example <- function() list(verified = FALSE, baseline_only_verified = FALSE, arrays_verified = FALSE, 
    medication_dictionary_verified = FALSE, evidence = "", algorithm = "Eastwood2016_Fig2_S1_conservative_baseline_v1")

ev_dm_medication_dictionary <- function() {
    code <- c(1140883066, 1140884600, 1140874686, 1140868902, 1140874646, 1140874658, 1140874664, 1140874674, 1140874706, 
        1140874718, 1140874744, 1141152590, 1141157284, 1141153254, 1141168660, 1141171646, 1141173882, 1141177600, 
        1141189090, 1140868908, 1140874650, 1140874652, 1140874660, 1140874666, 1140874678, 1140874680, 1140874690, 
        1140874712, 1140874716, 1140874724, 1140874726, 1140874728, 1140874732, 1140874736, 1140874740, 1140874746, 
        1140874754, 1140874756, 1140882964, 1141156984, 1141153262, 1141168668, 1141169504, 1141171508, 1141171652, 
        1141173786, 1141177606, 1141189094)
    meaning <- c("insulin product", "metformin", "glucophage 500mg tablet", "acarbose", "glipizide", "gliquidone", 
        "tolazamide", "tolbutamide", "chlorpropamide", "glibenclamide", "gliclazide", "glimepiride", "glipizide product", 
        "troglitazone", "repaglinide", "pioglitazone", "nateglinide", "rosiglitazone", "rosiglitazone 1mg / metformin 500mg tablet", 
        "glucobay 50mg tablet", "glibenese 5mg tablet", "minodiab 2.5mg tablet", "glurenorm 30mg tablet", "tolanase 100mg tablet", 
        "glyconon 500mg tablet", "rastinon 500mg tablet", "orabet 500mg tablet", "diabinese 100mg tablet", "glymese 250mg tablet", 
        "daonil 5mg tablet", "semi-daonil 2.5mg tablet", "euglucon 2.5mg tablet", "malix 2.5mg tablet", "diabetamide 2.5mg tablet", 
        "calabren 2.5mg tablet", "diamicron 80mg tablet", "guar gum", "guarem 5g/sachet granules", "oral hypoglycaemic", 
        "amaryl 1mg tablet", "romozin 200mg tablet", "novonorm 0.5mg tablet", "diaglyk 80mg tablet", "vivazide 80mg tablet", 
        "actos 15mg tablet", "starlix 60mg tablet", "avandia 4mg tablet", "avandamet 1mg / 500mg tablet")
    data.frame(code = code, meaning = meaning, category = c("insulin", rep("metformin", 2L), rep("nonmetformin_oral", 
        16L), rep("diabetes_medication_or_uncertain", 29L)), stringsAsFactors = FALSE)
}

ev_dm_num <- function(x, label) {
    z <- trimws(as.character(x))
    miss <- is.na(x) | toupper(z) %in% c("", "NA", "NAN", "NULL")
    ok <- miss | grepl("^-?[0-9]+(?:\\.[0-9]+)?$", z, perl = TRUE)
    if (any(!ok)) 
        stop("Unparsed baseline numeric field: ", label)
    ans <- suppressWarnings(as.numeric(z))
    ans[miss] <- NA_real_
    ans
}

ev_dm_columns <- function(raw, field, array = FALSE) {
    pat <- if (array) 
        paste0("^p", field, "_i0_a[0-9]+$")
    else paste0("^p", field, "_i0(?:_a[0-9]+)?$")
    out <- grep(pat, names(raw), value = TRUE, perl = TRUE)
    if (array && length(out)) 
        out <- out[order(as.integer(sub(".*_a", "", out)))]
    out
}

ev_dm_array <- function(raw, cols, count, label) {
    n <- nrow(raw)
    vals <- lapply(cols, function(x) ev_dm_num(raw[[x]], x))
    used <- integer(n)
    unknown <- rep(FALSE, n)
    for (x in vals) {
        used <- used + as.integer(!is.na(x))
        unknown <- unknown | (!is.na(x) & (x < 0 | x != floor(x)))
    }
    complete <- !is.na(count) & count >= 0 & count == floor(count) & count == used & !unknown
    list(values = vals, complete = complete, used = used, has = function(codes) Reduce(`|`, lapply(vals, function(x) !is.na(x) & 
        x %in% codes), init = rep(FALSE, n)))
}

ev_dm_touchscreen_insulin <- function(raw, cols, allowed) {
    n <- nrow(raw)
    insulin <- observed <- unknown <- none <- other <- rep(FALSE, n)
    for (col in cols) {
        zz <- trimws(as.character(raw[[col]]))
        zz[is.na(zz)] <- ""
        tokens <- strsplit(zz, "[|,;[:space:]]+", perl = TRUE)
        for (i in which(nzchar(zz) & !toupper(zz) %in% c("NA", "NAN", "NULL"))) {
            t <- tokens[[i]]
            t <- t[nzchar(t)]
            if (any(!grepl("^-?[0-9]+$", t))) 
                stop("Unparsed insulin touchscreen choices: ", col)
            t <- as.numeric(t)
            observed[i] <- TRUE
            unknown[i] <- unknown[i] | any(!t %in% c(allowed, -7))
            insulin[i] <- insulin[i] | 3 %in% t
            none[i] <- none[i] | -7 %in% t
            other[i] <- other[i] | any(t > 0)
        }
    }
    list(insulin = insulin, known = observed & !unknown & !(none & other))
}

ev_diabetes_baseline <- function(raw, eids = NULL, contract = ev_dm_contract_example()) {
    raw <- as.data.frame(raw, stringsAsFactors = FALSE)
    if (!"eid" %in% names(raw) || anyNA(raw$eid) || anyDuplicated(raw$eid) || anyDuplicated(names(raw))) 
        stop("Baseline raw input requires unique nonmissing eid and column names")
    if (is.null(eids)) 
        eids <- raw$eid
    if (anyNA(eids) || anyDuplicated(eids)) 
        stop("Requested baseline eid not unique/nonmissing")
    ix <- match(as.character(eids), as.character(raw$eid))
    linked <- !is.na(ix)
    raw <- raw[ix, , drop = FALSE]
    n <- length(eids)
    scalar_ids <- c(31, 135, 137, 2443, 4041, 2976, 2986, 21000)
    scalar <- paste0("p", scalar_ids, "_i0")
    arrays <- lapply(c(20002, 20003, 20009), function(id) ev_dm_columns(raw, id, TRUE))
    names(arrays) <- c("nurse", "medication", "nurse_age")
    ts <- list(female = ev_dm_columns(raw, 6153), male = ev_dm_columns(raw, 6177))
    missing <- c(setdiff(scalar, names(raw)), names(arrays)[!lengths(arrays)], names(ts)[!lengths(ts)])
    positions_match <- length(arrays$nurse) > 0L && identical(sub("^p20002", "", arrays$nurse), sub("^p20009", 
        "", arrays$nurse_age))
    evidence <- if (is.null(contract$evidence)) 
        ""
    else as.character(contract$evidence)
    gates <- data.frame(gate = c("source_contract", "baseline_only", "all_source_arrays", "medication_dictionary", 
        "required_fields_present", "nurse_code_age_positions"), passed = c(isTRUE(contract$verified) && length(evidence) == 
        1L && !is.na(evidence) && nzchar(trimws(evidence)), isTRUE(contract$baseline_only_verified), isTRUE(contract$arrays_verified), 
        isTRUE(contract$medication_dictionary_verified), !length(missing), positions_match), evidence = c(rep(if (length(evidence) == 
        1L && !is.na(evidence)) evidence else "", 4L), paste(missing, collapse = ";"), "20002_i0_aN and 20009_i0_aN must match exactly"), 
        stringsAsFactors = FALSE)
    if (!all(gates$passed)) 
        return(list(status = rep("unknown", n), classification = rep("not_adjudicated", n), eid = eids, verified = FALSE, 
            gates = gates, audit = data.frame(item = "blocked_gate", n = sum(!gates$passed)), flags = NULL, scope = "BLOCKED: required baseline fields or source contract unaudited"))
    get <- function(id) ev_dm_num(raw[[paste0("p", id, "_i0")]], paste0("p", id))
    sex <- get(31)
    doctor <- get(2443)
    gest <- get(4041)
    early <- get(2986)
    for (z in list(doctor, gest, early)) if (any(!is.na(z) & !z %in% c(0, 1, -1, -3, -2))) 
        stop("Unexpected baseline diabetes Yes/No code; audit coding before use")
    ni <- ev_dm_array(raw, arrays$nurse, get(135), "nurse")
    rx <- ev_dm_array(raw, arrays$medication, get(137), "medication")
    female <- !is.na(sex) & sex == 0
    male <- !is.na(sex) & sex == 1
    sf <- ev_dm_touchscreen_insulin(raw, ts$female, 1:5)
    sm <- ev_dm_touchscreen_insulin(raw, ts$male, 1:3)
    insulin_ts <- (female & sf$insulin) | (male & sm$insulin)
    ts_known <- (female & sf$known) | (male & sm$known)
    meds <- ev_dm_medication_dictionary()
    insulin_ni <- rx$has(meds$code[meds$category == "insulin"])
    metformin <- rx$has(meds$code[meds$category == "metformin"]) | rx$has(1141189090)
    oral <- rx$has(meds$code[meds$category == "nonmetformin_oral"])
    rx_uncertain <- rx$has(c(99999, meds$code[meds$category == "diabetes_medication_or_uncertain"]))
    ni_uncertain <- ni$has(99999)
    ni_t1_raw <- ni$has(1222)
    ni_t2_raw <- ni$has(1223)
    both_types <- ni_t1_raw & ni_t2_raw
    ni_t1 <- ni_t1_raw & !both_types
    ni_any <- ni$has(c(1220, 1222, 1223))
    ni_gdm <- ni$has(1221)
    complications <- ni$has(c(1276, 1468, 1607))
    gdm_ts <- female & !is.na(doctor) & doctor == 1 & !is.na(gest) & gest == 1
    gdm_any <- female & (gdm_ts | ni_gdm)
    ages <- lapply(arrays$nurse_age, function(x) ev_dm_num(raw[[x]], x))
    age_ni <- age_gdm <- rep(NA_real_, n)
    age_conflict <- rep(FALSE, n)
    for (j in seq_along(ni$values)) {
        cc <- ni$values[[j]]
        aa <- ages[[j]]
        valid <- !is.na(aa) & aa >= 0 & is.finite(aa)
        dm <- !is.na(cc) & cc %in% c(1220, 1222, 1223) & valid
        age_conflict <- age_conflict | (dm & !is.na(age_ni) & abs(aa - age_ni) > 1e-08)
        fill <- dm & is.na(age_ni)
        age_ni[fill] <- aa[fill]
        gg <- !is.na(cc) & cc == 1221 & valid & is.na(age_gdm)
        age_gdm[gg] <- aa[gg]
    }
    age_ts <- get(2976)
    age_ts[!is.na(age_ts) & age_ts < 0] <- NA_real_
    age <- ifelse(!is.na(age_ni), age_ni, age_ts)
    ethnic <- get(21000)
    european <- !is.na(ethnic) & ethnic %in% c(1, 1001, 1002, 1003)
    sa_ac <- !is.na(ethnic) & ethnic %in% c(3001, 3002, 3003, 4001)
    young <- !is.na(age) & ((european & age <= 36) | (sa_ac & age <= 30))
    age_route_known <- is.na(age) | european | sa_ac
    base_known <- linked & ni$complete & rx$complete & ts_known & !rx_uncertain & !ni_uncertain
    present <- ni_any | gdm_any | insulin_ts | insulin_ni | metformin | oral
    classification <- rep("unknown", n)
    terminal <- rep("unresolved", n)
    assign <- function(mask, value, rule) {
        mask[is.na(mask)] <- FALSE
        classification[mask] <<- value
        terminal[mask] <<- rule
    }
    assign(base_known & !present, "diabetes_unlikely", "1.1:no_NI_DM_GDM_or_diabetes_medication")
    no_meds <- !insulin_ts & !insulin_ni & !metformin & !oral
    gdm_only <- base_known & present & gdm_any & !ni_any & no_meds & (gdm_ts | (ni_gdm & !is.na(age_gdm) & age_gdm < 
        50))
    assign(gdm_only, "possible_gestational", "1.2:gestational_only")
    active <- base_known & present & !gdm_only & !age_conflict
    type1_arm <- active & !oral & age_route_known & young & (insulin_ts | insulin_ni | (!is.na(early) & early == 
        1) | ni_t1)
    assign(type1_arm, "possible_type1", "2.2:insufficient_probable_type1_evidence")
    assign(type1_arm & (ni_t1 | (!is.na(early) & early == 1 & (insulin_ts | insulin_ni))), "probable_type1", "2.1_or_2.2")
    early_known <- !is.na(early) & early %in% c(0, 1)
    type2_arm <- active & (oral | (age_route_known & !young) | (age_route_known & young & !insulin_ts & !insulin_ni & 
        !ni_t1 & early_known & early == 0))
    only_metformin <- metformin & !insulin_ni & !oral
    assign(type2_arm & only_metformin & !ni_any, "diabetes_unlikely", "3.2:metformin_without_NI_DM")
    rest <- type2_arm & !(only_metformin & !ni_any)
    assign(rest & oral, "probable_type2", "3.3:nonmetformin_oral")
    assign(rest & !oral & !insulin_ni, "probable_type2", "3.4:no_NI_insulin")
    assign(rest & !oral & insulin_ni & ni_t1, "probable_type1", "3.5:NI_type1")
    assign(rest & !oral & insulin_ni & !ni_t1, "possible_type2", "3.5:NI_insulin_without_NI_type1")
    status <- rep("unknown", n)
    status[classification == "probable_type1"] <- "type1"
    status[classification == "probable_type2"] <- "type2"
    positive_or_uncertain <- ni_any | ni_gdm | gdm_ts | (!is.na(doctor) & doctor == 1) | insulin_ts | insulin_ni | 
        metformin | oral | rx_uncertain | complications | (!is.na(early) & early == 1)
    status[positive_or_uncertain & !status %in% c("type1", "type2")] <- "any_diabetes_or_uncertain"
    status[both_types | age_conflict] <- "any_diabetes_or_uncertain"
    free <- base_known & classification == "diabetes_unlikely" & !is.na(doctor) & doctor == 0 & !positive_or_uncertain
    status[free] <- "no_diabetes"
    flags <- data.frame(eid = eids, baseline_status = status, classification = classification, terminal_rule = terminal, 
        raw_linked = linked, nurse_complete = ni$complete, medication_complete = rx$complete, touchscreen_insulin_known = ts_known, 
        nurse_any_diabetes = ni_any, nurse_type1_type2_conflict = both_types, nurse_diagnostic_age_conflict = age_conflict, 
        doctor_yes = !is.na(doctor) & doctor == 1, gestational_report = ni_gdm | gdm_ts, insulin_TS = insulin_ts, 
        insulin_NI = insulin_ni, metformin_NI = metformin, nonmetformin_oral_NI = oral, medication_uncertain = rx_uncertain, 
        diabetic_complication = complications, diagnostic_age = age, age_routing_supported = age_route_known, stringsAsFactors = FALSE)
    tab <- function(x, prefix) {
        tt <- table(x, useNA = "always")
        data.frame(item = paste0(prefix, names(tt)), n = as.integer(tt), stringsAsFactors = FALSE)
    }
    audit <- rbind(tab(status, "status:"), tab(classification, "classification:"), data.frame(item = c("unlinked", 
        "nurse_incomplete", "medication_incomplete", "TS_insulin_unknown", "medication_uncertain", "both_types_NI", 
        "NI_age_conflict", "unsupported_ethnicity_age_routing", "doctor_TS_only_not_free"), n = c(sum(!linked), 
        sum(!ni$complete), sum(!rx$complete), sum(!ts_known), sum(rx_uncertain), sum(both_types), sum(age_conflict), 
        sum(active & !age_route_known), sum(classification == "diabetes_unlikely" & !is.na(doctor) & doctor == 
            1))))
    list(status = status, classification = classification, eid = eids, flags = flags, audit = audit, gates = gates, 
        verified = TRUE, algorithm = contract$algorithm, scope = paste("Conservative Eastwood-based baseline self-report adaptation; probable", 
            "type1/type2 separated, possible/GDM/TS-only/metformin-only/complications", "excluded from the diabetes-free risk set; not clinical absence or full Eastwood", 
            "reproduction. Hospital E11 outcomes remain first-record endpoints, not all-source incidence."))
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        ev_dm_contract_example <- function() list(verified = FALSE, baseline_only_verified = FALSE, arrays_verified = FALSE, 
            medication_dictionary_verified = FALSE, evidence = "", algorithm = "Eastwood2016_Fig2_S1_conservative_baseline_v1")
        ev_dm_medication_dictionary <- function() {
            code <- c(1140883066, 1140884600, 1140874686, 1140868902, 1140874646, 1140874658, 1140874664, 1140874674, 
                1140874706, 1140874718, 1140874744, 1141152590, 1141157284, 1141153254, 1141168660, 1141171646, 
                1141173882, 1141177600, 1141189090, 1140868908, 1140874650, 1140874652, 1140874660, 1140874666, 
                1140874678, 1140874680, 1140874690, 1140874712, 1140874716, 1140874724, 1140874726, 1140874728, 
                1140874732, 1140874736, 1140874740, 1140874746, 1140874754, 1140874756, 1140882964, 1141156984, 
                1141153262, 1141168668, 1141169504, 1141171508, 1141171652, 1141173786, 1141177606, 1141189094)
            meaning <- c("insulin product", "metformin", "glucophage 500mg tablet", "acarbose", "glipizide", "gliquidone", 
                "tolazamide", "tolbutamide", "chlorpropamide", "glibenclamide", "gliclazide", "glimepiride", "glipizide product", 
                "troglitazone", "repaglinide", "pioglitazone", "nateglinide", "rosiglitazone", "rosiglitazone 1mg / metformin 500mg tablet", 
                "glucobay 50mg tablet", "glibenese 5mg tablet", "minodiab 2.5mg tablet", "glurenorm 30mg tablet", 
                "tolanase 100mg tablet", "glyconon 500mg tablet", "rastinon 500mg tablet", "orabet 500mg tablet", 
                "diabinese 100mg tablet", "glymese 250mg tablet", "daonil 5mg tablet", "semi-daonil 2.5mg tablet", 
                "euglucon 2.5mg tablet", "malix 2.5mg tablet", "diabetamide 2.5mg tablet", "calabren 2.5mg tablet", 
                "diamicron 80mg tablet", "guar gum", "guarem 5g/sachet granules", "oral hypoglycaemic", "amaryl 1mg tablet", 
                "romozin 200mg tablet", "novonorm 0.5mg tablet", "diaglyk 80mg tablet", "vivazide 80mg tablet", 
                "actos 15mg tablet", "starlix 60mg tablet", "avandia 4mg tablet", "avandamet 1mg / 500mg tablet")
            data.frame(code = code, meaning = meaning, category = c("insulin", rep("metformin", 2L), rep("nonmetformin_oral", 
                16L), rep("diabetes_medication_or_uncertain", 29L)), stringsAsFactors = FALSE)
        }
        ev_dm_num <- function(x, label) {
            z <- trimws(as.character(x))
            miss <- is.na(x) | toupper(z) %in% c("", "NA", "NAN", "NULL")
            ok <- miss | grepl("^-?[0-9]+(?:\\.[0-9]+)?$", z, perl = TRUE)
            if (any(!ok)) stop("Unparsed baseline numeric field: ", label)
            ans <- suppressWarnings(as.numeric(z))
            ans[miss] <- NA_real_
            ans
        }
        ev_dm_columns <- function(raw, field, array = FALSE) {
            pat <- if (array) paste0("^p", field, "_i0_a[0-9]+$") else paste0("^p", field, "_i0(?:_a[0-9]+)?$")
            out <- grep(pat, names(raw), value = TRUE, perl = TRUE)
            if (array && length(out)) out <- out[order(as.integer(sub(".*_a", "", out)))]
            out
        }
        ev_dm_array <- function(raw, cols, count, label) {
            n <- nrow(raw)
            vals <- lapply(cols, function(x) ev_dm_num(raw[[x]], x))
            used <- integer(n)
            unknown <- rep(FALSE, n)
            for (x in vals) {
                used <- used + as.integer(!is.na(x))
                unknown <- unknown | (!is.na(x) & (x < 0 | x != floor(x)))
            }
            complete <- !is.na(count) & count >= 0 & count == floor(count) & count == used & !unknown
            list(values = vals, complete = complete, used = used, has = function(codes) Reduce(`|`, lapply(vals, 
                function(x) !is.na(x) & x %in% codes), init = rep(FALSE, n)))
        }
        ev_dm_touchscreen_insulin <- function(raw, cols, allowed) {
            n <- nrow(raw)
            insulin <- observed <- unknown <- none <- other <- rep(FALSE, n)
            for (col in cols) {
                zz <- trimws(as.character(raw[[col]]))
                zz[is.na(zz)] <- ""
                tokens <- strsplit(zz, "[|,;[:space:]]+", perl = TRUE)
                for (i in which(nzchar(zz) & !toupper(zz) %in% c("NA", "NAN", "NULL"))) {
                  t <- tokens[[i]]
                  t <- t[nzchar(t)]
                  if (any(!grepl("^-?[0-9]+$", t))) stop("Unparsed insulin touchscreen choices: ", col)
                  t <- as.numeric(t)
                  observed[i] <- TRUE
                  unknown[i] <- unknown[i] | any(!t %in% c(allowed, -7))
                  insulin[i] <- insulin[i] | 3 %in% t
                  none[i] <- none[i] | -7 %in% t
                  other[i] <- other[i] | any(t > 0)
                }
            }
            list(insulin = insulin, known = observed & !unknown & !(none & other))
        }
        ev_diabetes_baseline <- function(raw, eids = NULL, contract = ev_dm_contract_example()) {
            raw <- as.data.frame(raw, stringsAsFactors = FALSE)
            if (!"eid" %in% names(raw) || anyNA(raw$eid) || anyDuplicated(raw$eid) || anyDuplicated(names(raw))) stop("Baseline raw input requires unique nonmissing eid and column names")
            if (is.null(eids)) eids <- raw$eid
            if (anyNA(eids) || anyDuplicated(eids)) stop("Requested baseline eid not unique/nonmissing")
            ix <- match(as.character(eids), as.character(raw$eid))
            linked <- !is.na(ix)
            raw <- raw[ix, , drop = FALSE]
            n <- length(eids)
            scalar_ids <- c(31, 135, 137, 2443, 4041, 2976, 2986, 21000)
            scalar <- paste0("p", scalar_ids, "_i0")
            arrays <- lapply(c(20002, 20003, 20009), function(id) ev_dm_columns(raw, id, TRUE))
            names(arrays) <- c("nurse", "medication", "nurse_age")
            ts <- list(female = ev_dm_columns(raw, 6153), male = ev_dm_columns(raw, 6177))
            missing <- c(setdiff(scalar, names(raw)), names(arrays)[!lengths(arrays)], names(ts)[!lengths(ts)])
            positions_match <- length(arrays$nurse) > 0L && identical(sub("^p20002", "", arrays$nurse), sub("^p20009", 
                "", arrays$nurse_age))
            evidence <- if (is.null(contract$evidence)) "" else as.character(contract$evidence)
            gates <- data.frame(gate = c("source_contract", "baseline_only", "all_source_arrays", "medication_dictionary", 
                "required_fields_present", "nurse_code_age_positions"), passed = c(isTRUE(contract$verified) && 
                length(evidence) == 1L && !is.na(evidence) && nzchar(trimws(evidence)), isTRUE(contract$baseline_only_verified), 
                isTRUE(contract$arrays_verified), isTRUE(contract$medication_dictionary_verified), !length(missing), 
                positions_match), evidence = c(rep(if (length(evidence) == 1L && !is.na(evidence)) evidence else "", 
                4L), paste(missing, collapse = ";"), "20002_i0_aN and 20009_i0_aN must match exactly"), stringsAsFactors = FALSE)
            if (!all(gates$passed)) return(list(status = rep("unknown", n), classification = rep("not_adjudicated", 
                n), eid = eids, verified = FALSE, gates = gates, audit = data.frame(item = "blocked_gate", n = sum(!gates$passed)), 
                flags = NULL, scope = "BLOCKED: required baseline fields or source contract unaudited"))
            get <- function(id) ev_dm_num(raw[[paste0("p", id, "_i0")]], paste0("p", id))
            sex <- get(31)
            doctor <- get(2443)
            gest <- get(4041)
            early <- get(2986)
            for (z in list(doctor, gest, early)) if (any(!is.na(z) & !z %in% c(0, 1, -1, -3, -2))) stop("Unexpected baseline diabetes Yes/No code; audit coding before use")
            ni <- ev_dm_array(raw, arrays$nurse, get(135), "nurse")
            rx <- ev_dm_array(raw, arrays$medication, get(137), "medication")
            female <- !is.na(sex) & sex == 0
            male <- !is.na(sex) & sex == 1
            sf <- ev_dm_touchscreen_insulin(raw, ts$female, 1:5)
            sm <- ev_dm_touchscreen_insulin(raw, ts$male, 1:3)
            insulin_ts <- (female & sf$insulin) | (male & sm$insulin)
            ts_known <- (female & sf$known) | (male & sm$known)
            meds <- ev_dm_medication_dictionary()
            insulin_ni <- rx$has(meds$code[meds$category == "insulin"])
            metformin <- rx$has(meds$code[meds$category == "metformin"]) | rx$has(1141189090)
            oral <- rx$has(meds$code[meds$category == "nonmetformin_oral"])
            rx_uncertain <- rx$has(c(99999, meds$code[meds$category == "diabetes_medication_or_uncertain"]))
            ni_uncertain <- ni$has(99999)
            ni_t1_raw <- ni$has(1222)
            ni_t2_raw <- ni$has(1223)
            both_types <- ni_t1_raw & ni_t2_raw
            ni_t1 <- ni_t1_raw & !both_types
            ni_any <- ni$has(c(1220, 1222, 1223))
            ni_gdm <- ni$has(1221)
            complications <- ni$has(c(1276, 1468, 1607))
            gdm_ts <- female & !is.na(doctor) & doctor == 1 & !is.na(gest) & gest == 1
            gdm_any <- female & (gdm_ts | ni_gdm)
            ages <- lapply(arrays$nurse_age, function(x) ev_dm_num(raw[[x]], x))
            age_ni <- age_gdm <- rep(NA_real_, n)
            age_conflict <- rep(FALSE, n)
            for (j in seq_along(ni$values)) {
                cc <- ni$values[[j]]
                aa <- ages[[j]]
                valid <- !is.na(aa) & aa >= 0 & is.finite(aa)
                dm <- !is.na(cc) & cc %in% c(1220, 1222, 1223) & valid
                age_conflict <- age_conflict | (dm & !is.na(age_ni) & abs(aa - age_ni) > 1e-08)
                fill <- dm & is.na(age_ni)
                age_ni[fill] <- aa[fill]
                gg <- !is.na(cc) & cc == 1221 & valid & is.na(age_gdm)
                age_gdm[gg] <- aa[gg]
            }
            age_ts <- get(2976)
            age_ts[!is.na(age_ts) & age_ts < 0] <- NA_real_
            age <- ifelse(!is.na(age_ni), age_ni, age_ts)
            ethnic <- get(21000)
            european <- !is.na(ethnic) & ethnic %in% c(1, 1001, 1002, 1003)
            sa_ac <- !is.na(ethnic) & ethnic %in% c(3001, 3002, 3003, 4001)
            young <- !is.na(age) & ((european & age <= 36) | (sa_ac & age <= 30))
            age_route_known <- is.na(age) | european | sa_ac
            base_known <- linked & ni$complete & rx$complete & ts_known & !rx_uncertain & !ni_uncertain
            present <- ni_any | gdm_any | insulin_ts | insulin_ni | metformin | oral
            classification <- rep("unknown", n)
            terminal <- rep("unresolved", n)
            assign <- function(mask, value, rule) {
                mask[is.na(mask)] <- FALSE
                classification[mask] <<- value
                terminal[mask] <<- rule
            }
            assign(base_known & !present, "diabetes_unlikely", "1.1:no_NI_DM_GDM_or_diabetes_medication")
            no_meds <- !insulin_ts & !insulin_ni & !metformin & !oral
            gdm_only <- base_known & present & gdm_any & !ni_any & no_meds & (gdm_ts | (ni_gdm & !is.na(age_gdm) & 
                age_gdm < 50))
            assign(gdm_only, "possible_gestational", "1.2:gestational_only")
            active <- base_known & present & !gdm_only & !age_conflict
            type1_arm <- active & !oral & age_route_known & young & (insulin_ts | insulin_ni | (!is.na(early) & 
                early == 1) | ni_t1)
            assign(type1_arm, "possible_type1", "2.2:insufficient_probable_type1_evidence")
            assign(type1_arm & (ni_t1 | (!is.na(early) & early == 1 & (insulin_ts | insulin_ni))), "probable_type1", 
                "2.1_or_2.2")
            early_known <- !is.na(early) & early %in% c(0, 1)
            type2_arm <- active & (oral | (age_route_known & !young) | (age_route_known & young & !insulin_ts & 
                !insulin_ni & !ni_t1 & early_known & early == 0))
            only_metformin <- metformin & !insulin_ni & !oral
            assign(type2_arm & only_metformin & !ni_any, "diabetes_unlikely", "3.2:metformin_without_NI_DM")
            rest <- type2_arm & !(only_metformin & !ni_any)
            assign(rest & oral, "probable_type2", "3.3:nonmetformin_oral")
            assign(rest & !oral & !insulin_ni, "probable_type2", "3.4:no_NI_insulin")
            assign(rest & !oral & insulin_ni & ni_t1, "probable_type1", "3.5:NI_type1")
            assign(rest & !oral & insulin_ni & !ni_t1, "possible_type2", "3.5:NI_insulin_without_NI_type1")
            status <- rep("unknown", n)
            status[classification == "probable_type1"] <- "type1"
            status[classification == "probable_type2"] <- "type2"
            positive_or_uncertain <- ni_any | ni_gdm | gdm_ts | (!is.na(doctor) & doctor == 1) | insulin_ts | insulin_ni | 
                metformin | oral | rx_uncertain | complications | (!is.na(early) & early == 1)
            status[positive_or_uncertain & !status %in% c("type1", "type2")] <- "any_diabetes_or_uncertain"
            status[both_types | age_conflict] <- "any_diabetes_or_uncertain"
            free <- base_known & classification == "diabetes_unlikely" & !is.na(doctor) & doctor == 0 & !positive_or_uncertain
            status[free] <- "no_diabetes"
            flags <- data.frame(eid = eids, baseline_status = status, classification = classification, terminal_rule = terminal, 
                raw_linked = linked, nurse_complete = ni$complete, medication_complete = rx$complete, touchscreen_insulin_known = ts_known, 
                nurse_any_diabetes = ni_any, nurse_type1_type2_conflict = both_types, nurse_diagnostic_age_conflict = age_conflict, 
                doctor_yes = !is.na(doctor) & doctor == 1, gestational_report = ni_gdm | gdm_ts, insulin_TS = insulin_ts, 
                insulin_NI = insulin_ni, metformin_NI = metformin, nonmetformin_oral_NI = oral, medication_uncertain = rx_uncertain, 
                diabetic_complication = complications, diagnostic_age = age, age_routing_supported = age_route_known, 
                stringsAsFactors = FALSE)
            tab <- function(x, prefix) {
                tt <- table(x, useNA = "always")
                data.frame(item = paste0(prefix, names(tt)), n = as.integer(tt), stringsAsFactors = FALSE)
            }
            audit <- rbind(tab(status, "status:"), tab(classification, "classification:"), data.frame(item = c("unlinked", 
                "nurse_incomplete", "medication_incomplete", "TS_insulin_unknown", "medication_uncertain", "both_types_NI", 
                "NI_age_conflict", "unsupported_ethnicity_age_routing", "doctor_TS_only_not_free"), n = c(sum(!linked), 
                sum(!ni$complete), sum(!rx$complete), sum(!ts_known), sum(rx_uncertain), sum(both_types), sum(age_conflict), 
                sum(active & !age_route_known), sum(classification == "diabetes_unlikely" & !is.na(doctor) & doctor == 
                  1))))
            list(status = status, classification = classification, eid = eids, flags = flags, audit = audit, gates = gates, 
                verified = TRUE, algorithm = contract$algorithm, scope = paste("Conservative Eastwood-based baseline self-report adaptation; probable", 
                  "type1/type2 separated, possible/GDM/TS-only/metformin-only/complications", "excluded from the diabetes-free risk set; not clinical absence or full Eastwood", 
                  "reproduction. Hospital E11 outcomes remain first-record endpoints, not all-source incidence."))
        }
    }), envir = context)
}

