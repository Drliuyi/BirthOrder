# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

body_tqeD_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

body_tqeD_summary <- function(A11, d, w, engine, methods, reference) {
    ids <- c("bmi", "weight", "fat_free_mass", "fat_mass", "hip", "waist")
    cols <- vapply(A11$members[ids], `[[`, "", "column")
    values <- setNames(as.list(d[cols]), ids)
    W <- values$weight
    F <- values$fat_mass
    L <- values$fat_free_mass
    T <- F + L
    q <- F/T
    e <- W - T
    all_finite <- all(vapply(values, function(x) all(is.finite(x)), logical(1)))
    positive <- all(W > 0 & F > 0 & L > 0 & T > 0 & is.finite(q) & q > 0 & q < 1)
    grain <- reference$measurement_grain
    arithmetic <- 32 * .Machine$double.eps * pmax(1, abs(W), abs(F), abs(L), abs(T))
    frequencies <- NULL
    state_pass <- TRUE
    precision_pass <- TRUE
    if (!is.null(grain)) {
        grid <- round(e/grain)
        grid[grid == 0] <- 0
        remainder <- abs(e - grid * grain)
        precision_pass <- all(remainder <= arithmetic)
        known <- grid %in% reference$e_grid_states
        state_pass <- all(known)
        states <- sort(unique(c(reference$e_grid_states, grid)))
        frequencies <- do.call(rbind, lapply(states, function(s) {
            take <- grid == s
            data.frame(integer_grid_state = s, e_value = s * grain, n = sum(take), weighted_mass = sum(w[take])/sum(w), 
                observed_state = s %in% reference$e_grid_states)
        }))
        max_precision_residue <- max(remainder)
        component_excess <- pmax(F - W, L - W, 0)
        total_excess <- pmax(T - W, 0)
        containment_pass <- all(component_excess <= grain + arithmetic & total_excess <= grain + arithmetic)
        containment <- list(measurement_grain = grain, maximum_component_excess = max(component_excess), maximum_total_excess = max(total_excess), 
            maximum_component_excess_grains = max(component_excess)/grain, maximum_total_excess_grains = max(total_excess)/grain, 
            allowed_excess_grains = 1, arithmetic_roundoff_bound_max = max(arithmetic), passed = containment_pass, 
            criterion = "RAW_COMPONENTS_AND_TOTAL_NOT_ABOVE_WEIGHT_BY_MORE_THAN_ONE_OBSERVED_MEASUREMENT_GRAIN_PLUS_MACHINE_ROUNDOFF")
    }
    else {
        max_precision_residue <- NULL
        containment_pass <- TRUE
        containment <- list(measurement_grain = NULL, passed = TRUE, criterion = "NO_DECIMAL_GRAIN_IN_OBSERVED_RECORDS; NO_INVENTED_CONTAINMENT_TOLERANCE", 
            empirical_e_law = TRUE)
    }
    pairs <- utils::combn(ids, 2, simplify = FALSE)
    correlations <- do.call(rbind, lapply(pairs, function(pair) data.frame(left = pair[1L], right = pair[2L], rho = methods$nmrD_spearman(values[[pair[1L]]], 
        values[[pair[2L]]], w), stringsAsFactors = FALSE)))
    list(N = nrow(d), sum_weight = sum(w), weight_ESS = sum(w)^2/sum(w^2), physical_pass = all_finite && positive && 
        precision_pass && state_pass && containment_pass, raw = methods$nmrD_table(values, w, "BODY_ORIGINAL_REPORTED_UNITS"), 
        coordinates = methods$nmrD_table(list(T = T, q = q, e = e, F_over_W = F/W, L_over_W = L/W), w, "BODY_DERIVED_DIAGNOSTIC_NOT_NEW_MEDIATORS"), 
        selected_spearman = correlations, support = list(all_finite = all_finite, W_F_L_T_strict_positive = positive, 
            min_T = min(T), q_range = range(q), q_strictly_between_zero_one = all(q > 0 & q < 1), e_range = range(e), 
            e_measurement_grain = grain, e_precision_pass = precision_pass, e_states_within_observed_alphabet = state_pass, 
            e_max_precision_residue = max_precision_residue, containment = containment), e_class_frequency = frequencies)
}

body_tqe_diagnostics <- function(A11, d, w, engine, methods, seed, width = 11L, K = 3L) {
    body_tqeD_assert(is.data.frame(d) && "B" %in% names(d) && length(w) == nrow(d) && all(is.finite(w) & w > 0) && 
        is.finite(sum(w)) && K == as.integer(K) && K >= 1L && width %in% c(11L, 15L, 49L) && length(A11$members) == 
        11L && identical(A11$body_tqe$revision, "SELECTED44_BODY_TQE_V1"), "BODY_TQE_DIAGNOSTIC_EXACT_A11_DATA_WEIGHT_WIDTH_REQUIRED")
    needed <- c("nmrD_table", "nmrD_spearman", "nmr_compare_summaries", "nmrA_rng")
    body_tqeD_assert(all(vapply(needed, function(name) is.function(methods[[name]]), logical(1))), "BODY_TQE_DIAGNOSTIC_EXISTING_DESCRIPTIVE_HELPERS_REQUIRED")
    before <- methods$nmrA_rng()
    on.exit({
        do.call(RNGkind, as.list(before$kind))
        if (before$present) assign(".Random.seed", before$seed, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, 
            inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    encoded <- engine$body_tqe_encode(d, A11$members)
    precision <- encoded$precision$e_decimal_digits
    grain <- if (is.na(precision)) 
        NULL
    else 10^(-precision)
    grid <- if (is.null(grain)) 
        NULL
    else sort(unique(round(encoded$z$e/grain)))
    if (!is.null(grid)) 
        grid[grid == 0] <- 0
    reference <- list(measurement_grain = grain, e_grid_states = grid, measured_precision = encoded$precision, 
        source = "EXACT_FROZEN_LINE_B0000_COMPLETE_RECORDS", decimal_state_rounding_is_binary_arithmetic_canonicalization_not_data_edit = TRUE)
    observed <- body_tqeD_summary(A11, d, w, engine, methods, reference)
    arms <- c("observed_B", "B0", "B1")
    ids <- c("bmi", "weight", "fat_free_mass", "fat_mass", "hip", "waist")
    cols <- vapply(A11$members[ids], `[[`, "", "column")
    pooled <- setNames(lapply(arms, function(arm) vector("list", K)), arms)
    set.seed(as.integer(seed))
    failure <- NULL
    for (k in seq_len(K)) {
        U <- matrix(stats::runif(nrow(d) * width), nrow(d), byrow = TRUE)
        for (arm in arms) {
            b <- switch(arm, observed_B = d$B, B0 = 0, B1 = 1)
            draw <- tryCatch(engine$c39_draw_a(A11, d, b, U[, seq_len(11L), drop = FALSE]), error = identity)
            if (inherits(draw, "error")) {
                failure <- list(replicate = k, arm = arm, message = conditionMessage(draw), retries = 0L)
                break
            }
            pooled[[arm]][[k]] <- draw[cols]
        }
        if (!is.null(failure)) 
            break
    }
    if (!is.null(failure)) 
        return(list(schema = "BODY_TQE_PAIRED_ENGINEERING_DIAGNOSTIC_V1", diagnostic_complete = FALSE, physical_pass = FALSE, 
            catastrophic_screen_pass = FALSE, passed = FALSE, observed = observed, reference = reference, failure = failure, 
            official_MC_integrations = 0L, bootstrap_draws = 0L, fitted_models = 0L, retries = 0L))
    summaries <- lapply(pooled, function(chunks) body_tqeD_summary(A11, do.call(rbind, chunks), rep(w/K, K), engine, 
        methods, reference))
    comparison <- methods$nmr_compare_summaries(observed, summaries$observed_B)
    physical <- isTRUE(observed$physical_pass) && all(vapply(summaries, `[[`, logical(1), "physical_pass"))
    list(schema = "BODY_TQE_PAIRED_ENGINEERING_DIAGNOSTIC_V1", diagnostic_complete = TRUE, actual_rows = nrow(d), 
        Ksmall = as.integer(K), seed = as.integer(seed), uniform_dim = as.integer(width), uniform_prefix_slots = seq_len(11L), 
        paired_same_U_for_all_three_B_arms = TRUE, raw_members_unchanged = ids, derived_coordinates_are_diagnostics_not_members = TRUE, 
        reference = reference, observed = observed, summaries = summaries, observed_B_comparison = comparison, 
        physical_pass = physical, catastrophic_screen_pass = isTRUE(comparison$engineering_pass), passed = physical && 
            isTRUE(comparison$engineering_pass), exposure_arm_distribution_differences_not_gated = TRUE, ratio_screen = list(limits = c(0.5, 
            2), role = "EXISTING_BROAD_CATASTROPHIC_DISTORTION_REVIEW_ONLY_NOT_EQUIVALENCE", rho_delta_warning = 0.35, 
            rho_warning_is_gate = FALSE), official_MC_integrations = 0L, bootstrap_draws = 0L, fitted_models = 0L, 
        retries = 0L, rng_restored_on_exit = TRUE, interpretation = paste("Actual original BODY distributions, weighted mean/SD/q99/q999 and all 15 pairwise", 
            "weighted Spearman relationships are descriptive. Three draws give Monte Carlo-imprecise tails;", "the pre-existing broad ratio screen is not distributional equivalence or biological validation.", 
            "Physical containment uses the measured decimal grain, not an invented exact W=F+L rule."))
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        body_tqeD_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
        body_tqeD_summary <- function(A11, d, w, engine, methods, reference) {
            ids <- c("bmi", "weight", "fat_free_mass", "fat_mass", "hip", "waist")
            cols <- vapply(A11$members[ids], `[[`, "", "column")
            values <- setNames(as.list(d[cols]), ids)
            W <- values$weight
            F <- values$fat_mass
            L <- values$fat_free_mass
            T <- F + L
            q <- F/T
            e <- W - T
            all_finite <- all(vapply(values, function(x) all(is.finite(x)), logical(1)))
            positive <- all(W > 0 & F > 0 & L > 0 & T > 0 & is.finite(q) & q > 0 & q < 1)
            grain <- reference$measurement_grain
            arithmetic <- 32 * .Machine$double.eps * pmax(1, abs(W), abs(F), abs(L), abs(T))
            frequencies <- NULL
            state_pass <- TRUE
            precision_pass <- TRUE
            if (!is.null(grain)) {
                grid <- round(e/grain)
                grid[grid == 0] <- 0
                remainder <- abs(e - grid * grain)
                precision_pass <- all(remainder <= arithmetic)
                known <- grid %in% reference$e_grid_states
                state_pass <- all(known)
                states <- sort(unique(c(reference$e_grid_states, grid)))
                frequencies <- do.call(rbind, lapply(states, function(s) {
                  take <- grid == s
                  data.frame(integer_grid_state = s, e_value = s * grain, n = sum(take), weighted_mass = sum(w[take])/sum(w), 
                    observed_state = s %in% reference$e_grid_states)
                }))
                max_precision_residue <- max(remainder)
                component_excess <- pmax(F - W, L - W, 0)
                total_excess <- pmax(T - W, 0)
                containment_pass <- all(component_excess <= grain + arithmetic & total_excess <= grain + arithmetic)
                containment <- list(measurement_grain = grain, maximum_component_excess = max(component_excess), 
                  maximum_total_excess = max(total_excess), maximum_component_excess_grains = max(component_excess)/grain, 
                  maximum_total_excess_grains = max(total_excess)/grain, allowed_excess_grains = 1, arithmetic_roundoff_bound_max = max(arithmetic), 
                  passed = containment_pass, criterion = "RAW_COMPONENTS_AND_TOTAL_NOT_ABOVE_WEIGHT_BY_MORE_THAN_ONE_OBSERVED_MEASUREMENT_GRAIN_PLUS_MACHINE_ROUNDOFF")
            } else {
                max_precision_residue <- NULL
                containment_pass <- TRUE
                containment <- list(measurement_grain = NULL, passed = TRUE, criterion = "NO_DECIMAL_GRAIN_IN_OBSERVED_RECORDS; NO_INVENTED_CONTAINMENT_TOLERANCE", 
                  empirical_e_law = TRUE)
            }
            pairs <- utils::combn(ids, 2, simplify = FALSE)
            correlations <- do.call(rbind, lapply(pairs, function(pair) data.frame(left = pair[1L], right = pair[2L], 
                rho = methods$nmrD_spearman(values[[pair[1L]]], values[[pair[2L]]], w), stringsAsFactors = FALSE)))
            list(N = nrow(d), sum_weight = sum(w), weight_ESS = sum(w)^2/sum(w^2), physical_pass = all_finite && 
                positive && precision_pass && state_pass && containment_pass, raw = methods$nmrD_table(values, 
                w, "BODY_ORIGINAL_REPORTED_UNITS"), coordinates = methods$nmrD_table(list(T = T, q = q, e = e, 
                F_over_W = F/W, L_over_W = L/W), w, "BODY_DERIVED_DIAGNOSTIC_NOT_NEW_MEDIATORS"), selected_spearman = correlations, 
                support = list(all_finite = all_finite, W_F_L_T_strict_positive = positive, min_T = min(T), q_range = range(q), 
                  q_strictly_between_zero_one = all(q > 0 & q < 1), e_range = range(e), e_measurement_grain = grain, 
                  e_precision_pass = precision_pass, e_states_within_observed_alphabet = state_pass, e_max_precision_residue = max_precision_residue, 
                  containment = containment), e_class_frequency = frequencies)
        }
        body_tqe_diagnostics <- function(A11, d, w, engine, methods, seed, width = 11L, K = 3L) {
            body_tqeD_assert(is.data.frame(d) && "B" %in% names(d) && length(w) == nrow(d) && all(is.finite(w) & 
                w > 0) && is.finite(sum(w)) && K == as.integer(K) && K >= 1L && width %in% c(11L, 15L, 49L) && 
                length(A11$members) == 11L && identical(A11$body_tqe$revision, "SELECTED44_BODY_TQE_V1"), "BODY_TQE_DIAGNOSTIC_EXACT_A11_DATA_WEIGHT_WIDTH_REQUIRED")
            needed <- c("nmrD_table", "nmrD_spearman", "nmr_compare_summaries", "nmrA_rng")
            body_tqeD_assert(all(vapply(needed, function(name) is.function(methods[[name]]), logical(1))), "BODY_TQE_DIAGNOSTIC_EXISTING_DESCRIPTIVE_HELPERS_REQUIRED")
            before <- methods$nmrA_rng()
            on.exit({
                do.call(RNGkind, as.list(before$kind))
                if (before$present) assign(".Random.seed", before$seed, .GlobalEnv) else if (exists(".Random.seed", 
                  .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
            }, add = TRUE)
            encoded <- engine$body_tqe_encode(d, A11$members)
            precision <- encoded$precision$e_decimal_digits
            grain <- if (is.na(precision)) NULL else 10^(-precision)
            grid <- if (is.null(grain)) NULL else sort(unique(round(encoded$z$e/grain)))
            if (!is.null(grid)) grid[grid == 0] <- 0
            reference <- list(measurement_grain = grain, e_grid_states = grid, measured_precision = encoded$precision, 
                source = "EXACT_FROZEN_LINE_B0000_COMPLETE_RECORDS", decimal_state_rounding_is_binary_arithmetic_canonicalization_not_data_edit = TRUE)
            observed <- body_tqeD_summary(A11, d, w, engine, methods, reference)
            arms <- c("observed_B", "B0", "B1")
            ids <- c("bmi", "weight", "fat_free_mass", "fat_mass", "hip", "waist")
            cols <- vapply(A11$members[ids], `[[`, "", "column")
            pooled <- setNames(lapply(arms, function(arm) vector("list", K)), arms)
            set.seed(as.integer(seed))
            failure <- NULL
            for (k in seq_len(K)) {
                U <- matrix(stats::runif(nrow(d) * width), nrow(d), byrow = TRUE)
                for (arm in arms) {
                  b <- switch(arm, observed_B = d$B, B0 = 0, B1 = 1)
                  draw <- tryCatch(engine$c39_draw_a(A11, d, b, U[, seq_len(11L), drop = FALSE]), error = identity)
                  if (inherits(draw, "error")) {
                    failure <- list(replicate = k, arm = arm, message = conditionMessage(draw), retries = 0L)
                    break
                  }
                  pooled[[arm]][[k]] <- draw[cols]
                }
                if (!is.null(failure)) break
            }
            if (!is.null(failure)) return(list(schema = "BODY_TQE_PAIRED_ENGINEERING_DIAGNOSTIC_V1", diagnostic_complete = FALSE, 
                physical_pass = FALSE, catastrophic_screen_pass = FALSE, passed = FALSE, observed = observed, reference = reference, 
                failure = failure, official_MC_integrations = 0L, bootstrap_draws = 0L, fitted_models = 0L, retries = 0L))
            summaries <- lapply(pooled, function(chunks) body_tqeD_summary(A11, do.call(rbind, chunks), rep(w/K, 
                K), engine, methods, reference))
            comparison <- methods$nmr_compare_summaries(observed, summaries$observed_B)
            physical <- isTRUE(observed$physical_pass) && all(vapply(summaries, `[[`, logical(1), "physical_pass"))
            list(schema = "BODY_TQE_PAIRED_ENGINEERING_DIAGNOSTIC_V1", diagnostic_complete = TRUE, actual_rows = nrow(d), 
                Ksmall = as.integer(K), seed = as.integer(seed), uniform_dim = as.integer(width), uniform_prefix_slots = seq_len(11L), 
                paired_same_U_for_all_three_B_arms = TRUE, raw_members_unchanged = ids, derived_coordinates_are_diagnostics_not_members = TRUE, 
                reference = reference, observed = observed, summaries = summaries, observed_B_comparison = comparison, 
                physical_pass = physical, catastrophic_screen_pass = isTRUE(comparison$engineering_pass), passed = physical && 
                  isTRUE(comparison$engineering_pass), exposure_arm_distribution_differences_not_gated = TRUE, 
                ratio_screen = list(limits = c(0.5, 2), role = "EXISTING_BROAD_CATASTROPHIC_DISTORTION_REVIEW_ONLY_NOT_EQUIVALENCE", 
                  rho_delta_warning = 0.35, rho_warning_is_gate = FALSE), official_MC_integrations = 0L, bootstrap_draws = 0L, 
                fitted_models = 0L, retries = 0L, rng_restored_on_exit = TRUE, interpretation = paste("Actual original BODY distributions, weighted mean/SD/q99/q999 and all 15 pairwise", 
                  "weighted Spearman relationships are descriptive. Three draws give Monte Carlo-imprecise tails;", 
                  "the pre-existing broad ratio screen is not distributional equivalence or biological validation.", 
                  "Physical containment uses the measured decimal grain, not an invented exact W=F+L rule."))
        }
    }), envir = context)
}

