# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

nmrD_moments <- function(x, w) {
    nmrA_assert(is.numeric(x) && length(x) == length(w) && all(is.finite(x)) && all(is.finite(w) & w > 0), "NMR_DIAGNOSTIC_INVALID_VALUES_OR_WEIGHTS")
    w <- w/sum(w)
    mean <- sum(w * x)
    sd <- sqrt(sum(w * (x - mean)^2))
    o <- order(x)
    values <- x[o]
    cumulative <- cumsum(w[o])
    cumulative[length(cumulative)] <- 1
    quantile <- function(p) values[which(cumulative >= p)[1L]]
    c(N = length(x), min = min(x), max = max(x), zero_n = sum(x == 0), negative_n = sum(x < 0), mean = mean, SD = sd, 
        q50 = quantile(0.5), q99 = quantile(0.99), q999 = quantile(0.999))
}

nmrD_weighted_rank <- function(x, w) {
    o <- order(x)
    group <- cumsum(c(TRUE, diff(x[o]) != 0))
    mass <- as.numeric(rowsum(w[o], group, reorder = FALSE))/sum(w)
    midpoint <- cumsum(mass) - mass/2
    answer <- numeric(length(x))
    answer[o] <- midpoint[group]
    answer
}

nmrD_spearman <- function(x, y, w) {
    a <- nmrD_weighted_rank(x, w)
    b <- nmrD_weighted_rank(y, w)
    w <- w/sum(w)
    a <- a - sum(w * a)
    b <- b - sum(w * b)
    denominator <- sqrt(sum(w * a^2) * sum(w * b^2))
    if (denominator == 0) 
        return(NA_real_)
    sum(w * a * b)/denominator
}

nmrD_default_pairs <- function() list(c("nmr__nmr_23400", "nmr__nmr_23415"), c("nmr__nmr_23400", "nmr__nmr_23406"), 
    c("nmr__nmr_23406", "nmr__nmr_23418"), c("nmr__nmr_23406", "nmr__nmr_23568"), c("nmr__nmr_23568", "nmr__nmr_23569"), 
    c("nmr__nmr_23427", "nmr__nmr_23430"), c("nmr__nmr_23427", "nmr__nmr_23523"), c("nmr__nmr_23524", "nmr__nmr_23526"), 
    c("nmr__nmr_23527", "nmr__nmr_23611"), c("nmr__nmr_23506", "nmr__nmr_23596"), c("nmr__nmr_23535", "nmr__nmr_23617"), 
    c("nmr__nmr_23568", "nmr__nmr_23640"), c("nmr__nmr_23404", "nmr__nmr_23532"), c("nmr__nmr_23400", "nmr__nmr_23427"), 
    c("nmr__nmr_23451", "nmr__nmr_23459"), c("bmi", "nmr__nmr_23404"), c("bmi", "nmr__nmr_23435"), c("waist", "nmr__nmr_23400"))

nmrD_table <- function(values, w, kind, block = NA_character_) {
    if (!length(values)) 
        return(data.frame())
    do.call(rbind, lapply(names(values), function(key) data.frame(block = block, name = key, kind = kind, as.list(nmrD_moments(values[[key]], 
        w)), check.names = FALSE, stringsAsFactors = FALSE)))
}

nmr_observed_summary <- function(fit, df, w, pairs = nmrD_default_pairs()) {
    nmrA_assert(is.data.frame(df) && length(w) == nrow(df) && all(is.finite(w) & w > 0) && is.finite(sum(w)), "NMR_DIAGNOSTIC_INVALID_DATA_OR_WEIGHTS")
    raw <- nmrA_raw_values(df, fit$members[fit$nmr_ids])
    raw_statistics <- nmrD_table(raw, w, "ORIGINAL_REPORTED_UNITS")
    physical_failures <- list()
    component_tables <- list()
    normalized_tables <- list()
    closure_tables <- list()
    auxiliary_tables <- list()
    roundtrips <- list()
    patterns <- list()
    for (step in fit$extension) {
        if (!identical(step$type, "composition")) 
            next
        block <- step$block
        encoded <- tryCatch(nmrA_roundtrip(block, raw, nrow(df)), error = identity)
        if (inherits(encoded, "error")) {
            physical_failures[[block$id]] <- conditionMessage(encoded)
            next
        }
        enc <- encoded$encoded
        component_tables[[block$id]] <- nmrD_table(as.list(as.data.frame(enc$parts)), w, "DISJOINT_PHYSICAL_OR_REPORTED_BASIS_PART", 
            block$id)
        positive <- enc$total > 0
        if (all(positive)) 
            normalized_tables[[block$id]] <- nmrD_table(as.list(as.data.frame(enc$parts/enc$total)), w, "EXACTLY_CLOSED_PHYSICAL_PROPORTION", 
                block$id)
        if (identical(block$component_unit, "reported_percent")) 
            closure_tables[[block$id]] <- nmrD_table(list(reported_part_sum_minus_100 = enc$total - 100), w, if (is.null(block$fixed_total)) 
                "RETAINED_REPORTED_CLOSURE_DEVIATION"
            else "SOURCE_DEFINED_FIXED_CLOSURE_NUMERICAL_AUDIT", block$id)
        auxiliary_tables[[block$id]] <- nmrD_table(enc$aux, w, "RETAINED_RECORDED_RELATION_DISCREPANCY_NOT_NEW_MOLECULAR_PART", 
            block$id)
        roundtrips[[block$id]] <- encoded[-1L]
        patterns[[block$id]] <- list(unweighted_counts = as.list(table(enc$mask)), weighted_mass = as.list(tapply(w/sum(w), 
            enc$mask, sum)))
    }
    for (id in fit$nmr_ids) {
        field <- sub("^nmr__nmr_", "", id)
        bounds <- nmrA_member_bounds(field, fit$members, fit$raw_support)
        value <- raw[[paste0("F", field)]]
        if (any(value < bounds[1L] | value > bounds[2L])) 
            physical_failures[[paste0("RAW_", field)]] <- "REPORTED_RAW_VALUE_OUTSIDE_SOURCE_PHYSICAL_SUPPORT"
    }
    pair_rows <- lapply(pairs, function(pair) {
        if (!all(pair %in% names(fit$members))) 
            return(NULL)
        columns <- vapply(fit$members[pair], `[[`, "", "column")
        if (!all(columns %in% names(df))) 
            return(NULL)
        data.frame(left = pair[1L], right = pair[2L], rho = nmrD_spearman(df[[columns[1L]]], df[[columns[2L]]], 
            w), stringsAsFactors = FALSE)
    })
    bind <- function(x) if (length(Filter(function(z) nrow(z) > 0, x))) 
        do.call(rbind, x)
    else data.frame()
    list(schema = "NMR_RECONSTRUCTION_ENGINEERING_SUMMARY_V1", N = nrow(df), sum_weight = sum(w), weight_ESS = sum(w)^2/sum(w^2), 
        physical_pass = !length(physical_failures), physical_failures = physical_failures, raw = raw_statistics, 
        parts = bind(component_tables), physical_proportions = bind(normalized_tables), reported_closure = bind(closure_tables), 
        reported_discrepancy = bind(auxiliary_tables), zero_patterns = patterns, roundtrip = roundtrips, selected_spearman = if (length(Filter(Negate(is.null), 
            pair_rows))) do.call(rbind, pair_rows) else data.frame(), definitions = list(weighted_SD = "sqrt(sum(normalized_w*(x-weighted_mean)^2))", 
            quantiles = "WEIGHTED_EMPIRICAL_INVERSE_CDF_NO_INTERPOLATION_SAME_FOR_OBSERVED_AND_GENERATED", spearman = "WEIGHTED_CORRELATION_OF_WEIGHTED_MID_DISTRIBUTION_RANKS", 
            weight_ESS_is_row_weight_diagnostic_not_independent_participant_count = TRUE, exact_recorded_source_ratios_claimed = FALSE, 
            positive_parts_closed_separately_from_retained_reported_deviations = TRUE))
}

nmr_compare_summaries <- function(observed, generated, ratio_limits = c(0.5, 2), near_zero = 1e-08, rho_delta_warning = 0.35) {
    nmrA_assert(length(ratio_limits) == 2L && ratio_limits[1L] > 0 && ratio_limits[1L] < 1 && ratio_limits[2L] > 
        1 && near_zero > 0, "NMR_DIAGNOSTIC_THRESHOLDS_INVALID")
    columns <- c("mean", "SD", "q99", "q999")
    o <- observed$raw
    g <- generated$raw
    ix <- match(o$name, g$name)
    nmrA_assert(!anyNA(ix), "NMR_DIAGNOSTIC_RAW_INFORMATION_MISMATCH")
    g <- g[ix, , drop = FALSE]
    rows <- list()
    for (metric in columns) for (j in seq_len(nrow(o))) {
        eligible <- abs(o[[metric]][j]) > near_zero
        ratio <- if (eligible) 
            g[[metric]][j]/o[[metric]][j]
        else NA_real_
        rows[[length(rows) + 1L]] <- data.frame(field = o$name[j], metric = metric, observed = o[[metric]][j], 
            generated = g[[metric]][j], delta = g[[metric]][j] - o[[metric]][j], ratio = ratio, ratio_eligible = eligible, 
            catastrophic_ratio_flag = eligible && (!is.finite(ratio) || ratio < ratio_limits[1L] || ratio > ratio_limits[2L]), 
            stringsAsFactors = FALSE)
    }
    comparison <- do.call(rbind, rows)
    rho <- if (nrow(observed$selected_spearman) && nrow(generated$selected_spearman)) 
        merge(observed$selected_spearman, generated$selected_spearman, by = c("left", "right"), suffixes = c("_observed", 
            "_generated"))
    else data.frame()
    if (nrow(rho)) {
        rho$delta <- rho$rho_generated - rho$rho_observed
        rho$structure_warning <- is.finite(rho$delta) & abs(rho$delta) > rho_delta_warning
    }
    flags <- comparison[comparison$catastrophic_ratio_flag, , drop = FALSE]
    list(schema = "NMR_OBSERVED_B_ENGINEERING_COMPARISON_V1", physical_pass = isTRUE(generated$physical_pass), 
        comparison = comparison, selected_spearman = rho, catastrophic_flags = flags, engineering_pass = isTRUE(generated$physical_pass) && 
            !nrow(flags), thresholds = list(ratio_limits = ratio_limits, absolute_near_zero = near_zero, rho_absolute_delta_warning = rho_delta_warning, 
            rho_warning_is_distribution_gate = FALSE), interpretation = paste("Broad predefined engineering distortion screen, not equality of distributions,", 
            "not a statistical validation or precise 99.9th-percentile guarantee,", "and never selected using IE, endpoint significance, or outcomes."))
}

nmr_generated_summary <- function(fit, C, w, engine = fit$engine, Ksmall = 3L, seedfixed = 20261103L, uniforms = NULL, 
    include_exposure_arms = TRUE, pairs = nmrD_default_pairs(), ratio_limits = c(0.5, 2), near_zero = 1e-08) {
    nmrA_assert(is.data.frame(C) && "B" %in% names(C) && length(w) == nrow(C) && length(Ksmall) == 1L && Ksmall >= 
        1L && Ksmall == as.integer(Ksmall) && all(is.finite(w) & w > 0) && is.finite(sum(w)), "NMR_SIMULATION_DIAGNOSTIC_INPUT_INVALID")
    before <- nmrA_rng()
    on.exit({
        do.call(RNGkind, as.list(before$kind))
        if (before$present) assign(".Random.seed", before$seed, envir = .GlobalEnv) else if (exists(".Random.seed", 
            .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    if (is.null(uniforms)) 
        set.seed(seedfixed)
    arms <- c("observed_B", if (isTRUE(include_exposure_arms)) c("B0", "B1"))
    raw_columns <- vapply(fit$members[fit$nmr_ids], `[[`, "", "column")
    pair_ids <- intersect(unique(unlist(pairs)), names(fit$members))
    keep <- unique(c("B", raw_columns, vapply(fit$members[pair_ids], `[[`, "", "column")))
    pooled <- setNames(lapply(arms, function(x) vector("list", Ksmall)), arms)
    failure <- NULL
    for (k in seq_len(Ksmall)) {
        u <- if (is.null(uniforms)) 
            matrix(stats::runif(nrow(C) * fit$uniform_dim), nrow(C), byrow = TRUE)
        else if (is.function(uniforms)) 
            uniforms(k)
        else uniforms[[k]]
        nmrA_assert(is.matrix(u) && nrow(u) == nrow(C) && ncol(u) == fit$uniform_dim && all(is.finite(u) & u >= 
            0 & u < 1), "NMR_DIAGNOSTIC_PROVIDED_UNIFORM_INVALID")
        for (arm in arms) {
            b <- switch(arm, observed_B = C$B, B0 = 0, B1 = 1)
            result <- tryCatch(draw_nmr_composition_A(fit, C, b, u, engine = engine), error = identity)
            if (inherits(result, "error")) {
                failure <- list(replicate = k, arm = arm, message = conditionMessage(result), skipped_or_retried_draws = 0L)
                break
            }
            pooled[[arm]][[k]] <- result[intersect(keep, names(result))]
        }
        if (!is.null(failure)) 
            break
    }
    if (!is.null(failure)) 
        return(list(schema = "NMR_PAIRED_ENGINEERING_SIMULATION_V1", diagnostic_complete = FALSE, physical_pass = FALSE, 
            engineering_pass = FALSE, failure = failure, official_MC_integrations = 0L, bootstrap_draws = 0L, fitted_models = 0L, 
            retries = 0L))
    summaries <- lapply(pooled, function(chunks) nmr_observed_summary(fit, do.call(rbind, chunks), rep(w/Ksmall, 
        Ksmall), pairs))
    available_observed <- all(raw_columns %in% names(C))
    comparison <- if (available_observed) 
        nmr_compare_summaries(nmr_observed_summary(fit, C, w, pairs), summaries$observed_B, ratio_limits = ratio_limits, 
            near_zero = near_zero)
    else NULL
    list(schema = "NMR_PAIRED_ENGINEERING_SIMULATION_V1", diagnostic_complete = TRUE, Ksmall = as.integer(Ksmall), 
        seedfixed = if (is.null(uniforms)) seedfixed else NULL, original_participants = nrow(C), simulated_rows_per_arm = nrow(C) * 
            Ksmall, provided_uniforms = !is.null(uniforms), paired_same_uniforms_across_arms = TRUE, fixed_SELECTED_BASE11_uniform_prefix = 11L, 
        uniform_dim = fit$uniform_dim, physical_pass = all(vapply(summaries, `[[`, logical(1), "physical_pass")), 
        engineering_pass = if (is.null(comparison)) NA else comparison$engineering_pass, summaries = summaries, 
        observed_B_comparison = comparison, exposure_arm_distribution_differences_not_gated = TRUE, official_MC_integrations = 0L, 
        bootstrap_draws = 0L, fitted_models = 0L, interpretation = paste("Additional deterministic engineering simulation only; not a numbered bootstrap", 
            "or an official K50 estimate. Ksmall tail summaries are descriptive and Monte Carlo imprecise.", "A broad ratio flag requires review; distributional identity is not required."))
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        nmrD_moments <- function(x, w) {
            nmrA_assert(is.numeric(x) && length(x) == length(w) && all(is.finite(x)) && all(is.finite(w) & w > 
                0), "NMR_DIAGNOSTIC_INVALID_VALUES_OR_WEIGHTS")
            w <- w/sum(w)
            mean <- sum(w * x)
            sd <- sqrt(sum(w * (x - mean)^2))
            o <- order(x)
            values <- x[o]
            cumulative <- cumsum(w[o])
            cumulative[length(cumulative)] <- 1
            quantile <- function(p) values[which(cumulative >= p)[1L]]
            c(N = length(x), min = min(x), max = max(x), zero_n = sum(x == 0), negative_n = sum(x < 0), mean = mean, 
                SD = sd, q50 = quantile(0.5), q99 = quantile(0.99), q999 = quantile(0.999))
        }
        nmrD_weighted_rank <- function(x, w) {
            o <- order(x)
            group <- cumsum(c(TRUE, diff(x[o]) != 0))
            mass <- as.numeric(rowsum(w[o], group, reorder = FALSE))/sum(w)
            midpoint <- cumsum(mass) - mass/2
            answer <- numeric(length(x))
            answer[o] <- midpoint[group]
            answer
        }
        nmrD_spearman <- function(x, y, w) {
            a <- nmrD_weighted_rank(x, w)
            b <- nmrD_weighted_rank(y, w)
            w <- w/sum(w)
            a <- a - sum(w * a)
            b <- b - sum(w * b)
            denominator <- sqrt(sum(w * a^2) * sum(w * b^2))
            if (denominator == 0) return(NA_real_)
            sum(w * a * b)/denominator
        }
        nmrD_default_pairs <- function() list(c("nmr__nmr_23400", "nmr__nmr_23415"), c("nmr__nmr_23400", "nmr__nmr_23406"), 
            c("nmr__nmr_23406", "nmr__nmr_23418"), c("nmr__nmr_23406", "nmr__nmr_23568"), c("nmr__nmr_23568", "nmr__nmr_23569"), 
            c("nmr__nmr_23427", "nmr__nmr_23430"), c("nmr__nmr_23427", "nmr__nmr_23523"), c("nmr__nmr_23524", "nmr__nmr_23526"), 
            c("nmr__nmr_23527", "nmr__nmr_23611"), c("nmr__nmr_23506", "nmr__nmr_23596"), c("nmr__nmr_23535", "nmr__nmr_23617"), 
            c("nmr__nmr_23568", "nmr__nmr_23640"), c("nmr__nmr_23404", "nmr__nmr_23532"), c("nmr__nmr_23400", "nmr__nmr_23427"), 
            c("nmr__nmr_23451", "nmr__nmr_23459"), c("bmi", "nmr__nmr_23404"), c("bmi", "nmr__nmr_23435"), c("waist", 
                "nmr__nmr_23400"))
        nmrD_table <- function(values, w, kind, block = NA_character_) {
            if (!length(values)) return(data.frame())
            do.call(rbind, lapply(names(values), function(key) data.frame(block = block, name = key, kind = kind, 
                as.list(nmrD_moments(values[[key]], w)), check.names = FALSE, stringsAsFactors = FALSE)))
        }
        nmr_observed_summary <- function(fit, df, w, pairs = nmrD_default_pairs()) {
            nmrA_assert(is.data.frame(df) && length(w) == nrow(df) && all(is.finite(w) & w > 0) && is.finite(sum(w)), 
                "NMR_DIAGNOSTIC_INVALID_DATA_OR_WEIGHTS")
            raw <- nmrA_raw_values(df, fit$members[fit$nmr_ids])
            raw_statistics <- nmrD_table(raw, w, "ORIGINAL_REPORTED_UNITS")
            physical_failures <- list()
            component_tables <- list()
            normalized_tables <- list()
            closure_tables <- list()
            auxiliary_tables <- list()
            roundtrips <- list()
            patterns <- list()
            for (step in fit$extension) {
                if (!identical(step$type, "composition")) next
                block <- step$block
                encoded <- tryCatch(nmrA_roundtrip(block, raw, nrow(df)), error = identity)
                if (inherits(encoded, "error")) {
                  physical_failures[[block$id]] <- conditionMessage(encoded)
                  next
                }
                enc <- encoded$encoded
                component_tables[[block$id]] <- nmrD_table(as.list(as.data.frame(enc$parts)), w, "DISJOINT_PHYSICAL_OR_REPORTED_BASIS_PART", 
                  block$id)
                positive <- enc$total > 0
                if (all(positive)) normalized_tables[[block$id]] <- nmrD_table(as.list(as.data.frame(enc$parts/enc$total)), 
                  w, "EXACTLY_CLOSED_PHYSICAL_PROPORTION", block$id)
                if (identical(block$component_unit, "reported_percent")) closure_tables[[block$id]] <- nmrD_table(list(reported_part_sum_minus_100 = enc$total - 
                  100), w, if (is.null(block$fixed_total)) "RETAINED_REPORTED_CLOSURE_DEVIATION" else "SOURCE_DEFINED_FIXED_CLOSURE_NUMERICAL_AUDIT", 
                  block$id)
                auxiliary_tables[[block$id]] <- nmrD_table(enc$aux, w, "RETAINED_RECORDED_RELATION_DISCREPANCY_NOT_NEW_MOLECULAR_PART", 
                  block$id)
                roundtrips[[block$id]] <- encoded[-1L]
                patterns[[block$id]] <- list(unweighted_counts = as.list(table(enc$mask)), weighted_mass = as.list(tapply(w/sum(w), 
                  enc$mask, sum)))
            }
            for (id in fit$nmr_ids) {
                field <- sub("^nmr__nmr_", "", id)
                bounds <- nmrA_member_bounds(field, fit$members, fit$raw_support)
                value <- raw[[paste0("F", field)]]
                if (any(value < bounds[1L] | value > bounds[2L])) physical_failures[[paste0("RAW_", field)]] <- "REPORTED_RAW_VALUE_OUTSIDE_SOURCE_PHYSICAL_SUPPORT"
            }
            pair_rows <- lapply(pairs, function(pair) {
                if (!all(pair %in% names(fit$members))) return(NULL)
                columns <- vapply(fit$members[pair], `[[`, "", "column")
                if (!all(columns %in% names(df))) return(NULL)
                data.frame(left = pair[1L], right = pair[2L], rho = nmrD_spearman(df[[columns[1L]]], df[[columns[2L]]], 
                  w), stringsAsFactors = FALSE)
            })
            bind <- function(x) if (length(Filter(function(z) nrow(z) > 0, x))) do.call(rbind, x) else data.frame()
            list(schema = "NMR_RECONSTRUCTION_ENGINEERING_SUMMARY_V1", N = nrow(df), sum_weight = sum(w), weight_ESS = sum(w)^2/sum(w^2), 
                physical_pass = !length(physical_failures), physical_failures = physical_failures, raw = raw_statistics, 
                parts = bind(component_tables), physical_proportions = bind(normalized_tables), reported_closure = bind(closure_tables), 
                reported_discrepancy = bind(auxiliary_tables), zero_patterns = patterns, roundtrip = roundtrips, 
                selected_spearman = if (length(Filter(Negate(is.null), pair_rows))) do.call(rbind, pair_rows) else data.frame(), 
                definitions = list(weighted_SD = "sqrt(sum(normalized_w*(x-weighted_mean)^2))", quantiles = "WEIGHTED_EMPIRICAL_INVERSE_CDF_NO_INTERPOLATION_SAME_FOR_OBSERVED_AND_GENERATED", 
                  spearman = "WEIGHTED_CORRELATION_OF_WEIGHTED_MID_DISTRIBUTION_RANKS", weight_ESS_is_row_weight_diagnostic_not_independent_participant_count = TRUE, 
                  exact_recorded_source_ratios_claimed = FALSE, positive_parts_closed_separately_from_retained_reported_deviations = TRUE))
        }
        nmr_compare_summaries <- function(observed, generated, ratio_limits = c(0.5, 2), near_zero = 1e-08, rho_delta_warning = 0.35) {
            nmrA_assert(length(ratio_limits) == 2L && ratio_limits[1L] > 0 && ratio_limits[1L] < 1 && ratio_limits[2L] > 
                1 && near_zero > 0, "NMR_DIAGNOSTIC_THRESHOLDS_INVALID")
            columns <- c("mean", "SD", "q99", "q999")
            o <- observed$raw
            g <- generated$raw
            ix <- match(o$name, g$name)
            nmrA_assert(!anyNA(ix), "NMR_DIAGNOSTIC_RAW_INFORMATION_MISMATCH")
            g <- g[ix, , drop = FALSE]
            rows <- list()
            for (metric in columns) for (j in seq_len(nrow(o))) {
                eligible <- abs(o[[metric]][j]) > near_zero
                ratio <- if (eligible) g[[metric]][j]/o[[metric]][j] else NA_real_
                rows[[length(rows) + 1L]] <- data.frame(field = o$name[j], metric = metric, observed = o[[metric]][j], 
                  generated = g[[metric]][j], delta = g[[metric]][j] - o[[metric]][j], ratio = ratio, ratio_eligible = eligible, 
                  catastrophic_ratio_flag = eligible && (!is.finite(ratio) || ratio < ratio_limits[1L] || ratio > 
                    ratio_limits[2L]), stringsAsFactors = FALSE)
            }
            comparison <- do.call(rbind, rows)
            rho <- if (nrow(observed$selected_spearman) && nrow(generated$selected_spearman)) merge(observed$selected_spearman, 
                generated$selected_spearman, by = c("left", "right"), suffixes = c("_observed", "_generated")) else data.frame()
            if (nrow(rho)) {
                rho$delta <- rho$rho_generated - rho$rho_observed
                rho$structure_warning <- is.finite(rho$delta) & abs(rho$delta) > rho_delta_warning
            }
            flags <- comparison[comparison$catastrophic_ratio_flag, , drop = FALSE]
            list(schema = "NMR_OBSERVED_B_ENGINEERING_COMPARISON_V1", physical_pass = isTRUE(generated$physical_pass), 
                comparison = comparison, selected_spearman = rho, catastrophic_flags = flags, engineering_pass = isTRUE(generated$physical_pass) && 
                  !nrow(flags), thresholds = list(ratio_limits = ratio_limits, absolute_near_zero = near_zero, 
                  rho_absolute_delta_warning = rho_delta_warning, rho_warning_is_distribution_gate = FALSE), interpretation = paste("Broad predefined engineering distortion screen, not equality of distributions,", 
                  "not a statistical validation or precise 99.9th-percentile guarantee,", "and never selected using IE, endpoint significance, or outcomes."))
        }
        nmr_generated_summary <- function(fit, C, w, engine = fit$engine, Ksmall = 3L, seedfixed = 20261103L, uniforms = NULL, 
            include_exposure_arms = TRUE, pairs = nmrD_default_pairs(), ratio_limits = c(0.5, 2), near_zero = 1e-08) {
            nmrA_assert(is.data.frame(C) && "B" %in% names(C) && length(w) == nrow(C) && length(Ksmall) == 1L && 
                Ksmall >= 1L && Ksmall == as.integer(Ksmall) && all(is.finite(w) & w > 0) && is.finite(sum(w)), 
                "NMR_SIMULATION_DIAGNOSTIC_INPUT_INVALID")
            before <- nmrA_rng()
            on.exit({
                do.call(RNGkind, as.list(before$kind))
                if (before$present) assign(".Random.seed", before$seed, envir = .GlobalEnv) else if (exists(".Random.seed", 
                  .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
            }, add = TRUE)
            if (is.null(uniforms)) set.seed(seedfixed)
            arms <- c("observed_B", if (isTRUE(include_exposure_arms)) c("B0", "B1"))
            raw_columns <- vapply(fit$members[fit$nmr_ids], `[[`, "", "column")
            pair_ids <- intersect(unique(unlist(pairs)), names(fit$members))
            keep <- unique(c("B", raw_columns, vapply(fit$members[pair_ids], `[[`, "", "column")))
            pooled <- setNames(lapply(arms, function(x) vector("list", Ksmall)), arms)
            failure <- NULL
            for (k in seq_len(Ksmall)) {
                u <- if (is.null(uniforms)) matrix(stats::runif(nrow(C) * fit$uniform_dim), nrow(C), byrow = TRUE) else if (is.function(uniforms)) uniforms(k) else uniforms[[k]]
                nmrA_assert(is.matrix(u) && nrow(u) == nrow(C) && ncol(u) == fit$uniform_dim && all(is.finite(u) & 
                  u >= 0 & u < 1), "NMR_DIAGNOSTIC_PROVIDED_UNIFORM_INVALID")
                for (arm in arms) {
                  b <- switch(arm, observed_B = C$B, B0 = 0, B1 = 1)
                  result <- tryCatch(draw_nmr_composition_A(fit, C, b, u, engine = engine), error = identity)
                  if (inherits(result, "error")) {
                    failure <- list(replicate = k, arm = arm, message = conditionMessage(result), skipped_or_retried_draws = 0L)
                    break
                  }
                  pooled[[arm]][[k]] <- result[intersect(keep, names(result))]
                }
                if (!is.null(failure)) break
            }
            if (!is.null(failure)) return(list(schema = "NMR_PAIRED_ENGINEERING_SIMULATION_V1", diagnostic_complete = FALSE, 
                physical_pass = FALSE, engineering_pass = FALSE, failure = failure, official_MC_integrations = 0L, 
                bootstrap_draws = 0L, fitted_models = 0L, retries = 0L))
            summaries <- lapply(pooled, function(chunks) nmr_observed_summary(fit, do.call(rbind, chunks), rep(w/Ksmall, 
                Ksmall), pairs))
            available_observed <- all(raw_columns %in% names(C))
            comparison <- if (available_observed) nmr_compare_summaries(nmr_observed_summary(fit, C, w, pairs), 
                summaries$observed_B, ratio_limits = ratio_limits, near_zero = near_zero) else NULL
            list(schema = "NMR_PAIRED_ENGINEERING_SIMULATION_V1", diagnostic_complete = TRUE, Ksmall = as.integer(Ksmall), 
                seedfixed = if (is.null(uniforms)) seedfixed else NULL, original_participants = nrow(C), simulated_rows_per_arm = nrow(C) * 
                  Ksmall, provided_uniforms = !is.null(uniforms), paired_same_uniforms_across_arms = TRUE, fixed_SELECTED_BASE11_uniform_prefix = 11L, 
                uniform_dim = fit$uniform_dim, physical_pass = all(vapply(summaries, `[[`, logical(1), "physical_pass")), 
                engineering_pass = if (is.null(comparison)) NA else comparison$engineering_pass, summaries = summaries, 
                observed_B_comparison = comparison, exposure_arm_distribution_differences_not_gated = TRUE, official_MC_integrations = 0L, 
                bootstrap_draws = 0L, fitted_models = 0L, interpretation = paste("Additional deterministic engineering simulation only; not a numbered bootstrap", 
                  "or an official K50 estimate. Ksmall tail summaries are descriptive and Monte Carlo imprecise.", 
                  "A broad ratio flag requires review; distributional identity is not required."))
        }
    }), envir = context)
}

