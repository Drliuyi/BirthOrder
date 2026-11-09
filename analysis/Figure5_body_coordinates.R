# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

body_tqe_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

body_tqe_revision <- function() "SELECTED44_BODY_TQE_V1"

body_tqe_ids <- function() c("smoking", "bmi", "weight", "fat_free_mass", "fat_mass", "hip", "waist", "clinical__hdl", 
    "clinical__tg", "clinical__glucose", "clinical__hba1c")

body_tqe_mass_ids <- function() c("weight", "fat_free_mass", "fat_mass")

body_tqe_names <- function() c(logT = ".body_tqe_logT", logitq = ".body_tqe_logitq", e = ".body_tqe_e", affineT = ".body_tqe_affineT", 
    affineF = ".body_tqe_affineF", affinee = ".body_tqe_affinee")

body_tqe_hash <- function(x) {
    body_tqe_assert(requireNamespace("digest", quietly = TRUE), "BODY_TQE_DIGEST_REQUIRED")
    digest::digest(x, algo = "sha256", serialize = TRUE, serializeVersion = 2L)
}

body_tqe_binding <- function(d, members, base, w) {
    cols <- vapply(members[seq_len(11L)], `[[`, "", "column")
    body_tqe_hash(list(revision = body_tqe_revision(), members = members[seq_len(11L)], base = base, data = d[unique(c(base, 
        cols))], eid = if ("eid" %in% names(d)) as.character(d$eid), weights = as.numeric(w)))
}

body_tqe_precision <- function(x) {
    tolerance <- 32 * .Machine$double.eps * pmax(1, abs(x))
    for (k in 0:10) if (all(abs(x - round(x, k)) <= tolerance)) 
        return(k)
    NA_integer_
}

body_tqe_columns <- function(members) vapply(members[body_tqe_mass_ids()], `[[`, "", "column")

body_tqe_encode <- function(d, members) {
    cols <- body_tqe_columns(members)
    W <- d[[cols[["weight"]]]]
    L <- d[[cols[["fat_free_mass"]]]]
    F <- d[[cols[["fat_mass"]]]]
    body_tqe_assert(all(is.finite(W) & W > 0 & is.finite(L) & L > 0 & is.finite(F) & F > 0), "BODY_TQE_STRICTLY_POSITIVE_RAW_W_F_L_REQUIRED_NO_REPAIR")
    T <- F + L
    q <- F/T
    e <- W - T
    body_tqe_assert(all(is.finite(T) & T > 0 & is.finite(q) & q > 0 & q < 1), "BODY_TQE_INVALID_TOTAL_OR_FRACTION_NO_REPAIR")
    z <- list(logT = log(T), logitq = log(F) - log(L), e = e)
    F2 <- exp(z$logT) * stats::plogis(z$logitq)
    L2 <- exp(z$logT) * stats::plogis(-z$logitq)
    W2 <- F2 + L2 + e
    error <- c(W = max(abs(W2 - W)/pmax(1, abs(W))), F = max(abs(F2 - F)/pmax(1, abs(F))), L = max(abs(L2 - L)/pmax(1, 
        abs(L))))
    body_tqe_assert(all(is.finite(error) & error <= 1e-10), "BODY_TQE_RAW_ROUNDTRIP_LOSES_RECORDED_VARIATION")
    digits <- vapply(list(W, L, F), body_tqe_precision, integer(1))
    precision <- if (anyNA(digits)) 
        NA_integer_
    else max(digits)
    states <- if (is.na(precision)) 
        NULL
    else round(e * 10^precision)/10^precision
    if (!is.null(states)) 
        states[states == 0] <- 0
    state_error <- if (is.null(states)) 
        NA_real_
    else max(abs(states - e)/pmax(1, abs(W)))
    if (!is.null(states)) 
        body_tqe_assert(state_error <= 1e-10, "BODY_TQE_DECIMAL_STATE_COLLAPSE_LOSES_RECORDED_VARIATION")
    list(z = z, total = T, fraction = q, columns = cols, e_states = states, precision = list(raw_decimal_digits = digits, 
        e_decimal_digits = precision, e_binary_unique_n = length(unique(e)), e_state_n = if (is.null(states)) NULL else length(unique(states)), 
        max_scaled_state_residue = state_error, arithmetic_tolerance = "32*double.eps*pmax(1,abs(raw)); state roundtrip <=1e-10"), 
        roundtrip = list(max_scaled_error = error, tolerance = 1e-10))
}

body_tqe_scale <- function(x, w) {
    center <- sum(w * x)/sum(w)
    scale <- sqrt(sum(w * (x - center)^2)/sum(w))
    list(center = center, scale = if (scale == 0) 1 else scale)
}

body_tqe_standardize <- function(x, ref) (x - ref$center)/ref$scale

body_tqe_reference_copy <- function(reference) {
    body_tqe_assert(is.list(reference) && inherits(reference$formula, "formula"), "BODY_TQE_UPSTREAM_FORMULA_REFERENCE_REQUIRED")
    formula <- reference$formula
    environment(formula) <- baseenv()
    reference$formula <- formula
    reference
}

body_tqe_affine_prepare <- function(d, previous, columns, w = NULL, reference = NULL) {
    if (!all(columns %in% previous)) 
        return(list(data = d, previous = previous, reference = NULL))
    keys <- body_tqe_names()[c("affineT", "affineF", "affinee")]
    W <- d[[columns[["weight"]]]]
    L <- d[[columns[["fat_free_mass"]]]]
    F <- d[[columns[["fat_mass"]]]]
    values <- list(F + L, F, W - (F + L))
    if (is.null(reference)) 
        reference <- list(columns = columns, keys = keys, scale = lapply(values, body_tqe_scale, w = w), basis = "INVERTIBLE_AFFINE_T_F_e_SAME_RAW_W_L_F_PREDICTOR_SPAN")
    body_tqe_assert(identical(reference$columns, columns) && identical(reference$keys, keys), "BODY_TQE_AFFINE_PREDICTOR_REFERENCE_DRIFT")
    for (j in seq_along(keys)) d[[keys[[j]]]] <- body_tqe_standardize(values[[j]], reference$scale[[j]])
    previous[match(columns, previous)] <- unname(keys)
    list(data = d, previous = previous, reference = reference)
}

body_tqe_raw_fit <- function(d, member, base, previous, w, columns, continuous_rule, count_rule) {
    prepared <- body_tqe_affine_prepare(d, previous, columns, w)
    design <- c39_a_design(prepared$data, base, prepared$previous)
    terms <- attr(stats::terms(design$reference$formula), "term.labels")
    callback <- getOption("c39.a_checkpoint_callback")
    options(c39.a_checkpoint_callback = NULL)
    on.exit(options(c39.a_checkpoint_callback = callback), add = TRUE)
    a <- .body_tqe_original_fit_a(prepared$data, setNames(list(member), member$indicator_id), terms, w, continuous_rule, 
        count_rule)
    f <- a$models[[1L]]
    f$design_reference <- design$reference
    f$formula <- design$reference$formula
    f$predecessor_columns <- previous
    f$body_tqe_affine <- prepared$reference
    f
}

body_tqe_coordinate_fit <- function(d, y, base, previous, w, id, member) {
    design <- c39_a_design(d, base, previous)
    f <- c39_linear_fit(design$X, y, w, paste0("body_tqe__", id), "EMPIRICAL_RESIDUAL_UNTRUNCATED")
    f$design_reference <- design$reference
    f$columns <- colnames(design$X)
    f$formula <- design$reference$formula
    f$predecessor_columns <- previous
    f$member <- member
    f$body_tqe_role <- id
    f$body_tqe_scale <- body_tqe_scale(y, w)
    f
}

body_tqe_e_fit <- function(d, encoded, base, previous, w, member) {
    design <- c39_a_design(d, base, previous)
    states <- encoded$e_states
    if (is.null(states)) {
        f <- body_tqe_coordinate_fit(d, encoded$z$e, base, previous, w, "e", member)
        f$e_law <- "WEIGHTED_EMPIRICAL_RESIDUAL_CONDITIONAL_ON_W_POSITIVE"
        return(f)
    }
    values <- sort(unique(states))
    body_tqe_assert(length(values) <= 256L, "BODY_TQE_MEASURED_DISCRETE_e_OVER_256_STATES_REQUIRES_REVIEW")
    if (length(values) == 1L) 
        f <- list(kind = "body_tqe_category", values = values, coefficients = NULL, warnings = character())
    else {
        X <- design$X
        y <- factor(match(states, values), levels = seq_along(values))
        body_tqe_assert(requireNamespace("nnet", quietly = TRUE) && qr(X)$rank == ncol(X), "BODY_TQE_e_CATEGORY_DESIGN_RANK_DEFICIENT")
        warnings <- character()
        fit <- withCallingHandlers(nnet::multinom(y ~ X - 1, weights = w, trace = FALSE, maxit = 1000L, MaxNWts = 100000L, 
            model = FALSE), warning = function(e) {
            warnings <<- c(warnings, conditionMessage(e))
            invokeRestart("muffleWarning")
        })
        body_tqe_assert(identical(fit$convergence, 0L) && all(is.finite(stats::coef(fit))), "BODY_TQE_e_CATEGORY_NOT_CONVERGED")
        f <- list(kind = "body_tqe_category", values = values, coefficients = stats::coef(fit), warnings = unique(warnings), 
            convergence = fit$convergence)
    }
    f$design_reference <- design$reference
    f$columns <- colnames(design$X)
    f$formula <- design$reference$formula
    f$predecessor_columns <- previous
    f$member <- member
    f$body_tqe_role <- "e"
    f$e_law <- "MEASURED_DECIMAL_STATES_WEIGHTED_CONDITIONAL_MULTINOMIAL"
    f$body_tqe_scale <- body_tqe_scale(encoded$z$e, w)
    f
}

body_tqe_fit_a <- function(d, members, base, w, continuous_rule = "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL", 
    count_rule = "FIXED_SUPPORT_MULTINOMIAL_0_7", prefix_fit = getOption("body_tqe.prefix_A"), upstream_fit = getOption("body_tqe.upstream_A")) {
    ids <- vapply(members, `[[`, "", "indicator_id")
    body_tqe_assert(length(ids) >= 11L && identical(unname(ids[seq_len(11L)]), body_tqe_ids()) && identical(names(members), 
        unname(ids)), "BODY_TQE_FROZEN_11_PREFIX_REQUIRED")
    body_tqe_assert(length(w) == nrow(d) && all(is.finite(w) & w > 0) && is.finite(sum(w)), "BODY_TQE_EXACT_SELECTION_WEIGHTS_INVALID")
    body_tqe_assert(!any(body_tqe_names() %in% names(d)), "BODY_TQE_RESERVED_COLUMNS_EXIST")
    .body_tqe_original_support(d, members)
    binding <- body_tqe_binding(d, members, base, w)
    encoded <- body_tqe_encode(d, members)
    columns <- encoded$columns
    internal <- body_tqe_names()
    original_previous <- character()
    models <- list()
    warnings <- character()
    start <- 1L
    upstream_provenance <- NULL
    if (!is.null(prefix_fit)) {
        body_tqe_assert(identical(prefix_fit$body_tqe$revision, body_tqe_revision()) && identical(prefix_fit$body_tqe$prefix_binding, 
            binding) && length(prefix_fit$members) == 11L && identical(prefix_fit$members, members[seq_len(11L)]), 
            "BODY_TQE_PREFIX_SAMPLE_WEIGHTS_MEMBERS_BINDING_MISMATCH")
        models <- prefix_fit$models
        warnings <- prefix_fit$warnings
        original_previous <- vapply(members[seq_len(11L)], `[[`, "", "column")
        start <- 12L
        upstream_provenance <- prefix_fit$body_tqe$upstream_provenance
    }
    else if (!is.null(upstream_fit)) {
        upstream_ids <- body_tqe_ids()[1:2]
        body_tqe_assert(identical(upstream_fit$prefix_binding, binding) && identical(upstream_fit$members, members[1:2]) && 
            identical(names(upstream_fit$models), upstream_ids) && is.list(upstream_fit$provenance) && length(upstream_fit$provenance) > 
            0L, "BODY_TQE_OLD_UPSTREAM_INPUT_WEIGHTS_MEMBERS_PROVENANCE_MISMATCH")
        before <- character()
        for (id in upstream_ids) {
            f <- upstream_fit$models[[id]]
            expected_design <- c39_a_design(d, base, before)
            expected <- expected_design$reference
            old_design <- c39_a_design(d, NULL, NULL, f$design_reference)
            body_tqe_assert(identical(f$member, members[[id]]) && identical(body_tqe_reference_copy(f$design_reference), 
                body_tqe_reference_copy(expected)) && identical(old_design$X, expected_design$X) && identical(f$predecessor_columns, 
                before) && is.null(f$body_tqe_role), paste("BODY_TQE_OLD_UPSTREAM_CONDITIONAL_DESIGN_CHANGED", 
                id))
            before <- c(before, members[[id]]$column)
        }
        body_tqe_assert(identical(upstream_fit$models$bmi$log_a_revision$revision, "CORE39_BMI_LOG_FAT_LINEAR_V1"), 
            "BODY_TQE_OLD_UPSTREAM_APPROVED_LOG_BMI_REQUIRED")
        models <- upstream_fit$models
        warnings <- unique(unlist(lapply(models, `[[`, "warnings")))
        original_previous <- before
        start <- 3L
        upstream_provenance <- upstream_fit$provenance
    }
    if (start <= length(members)) 
        for (j in seq.int(start, length(members))) {
            m <- members[[j]]
            id <- m$indicator_id
            if (id == "weight") {
                f <- body_tqe_coordinate_fit(d, encoded$z$logT, base, original_previous, w, "logT", m)
                d[[internal[["logT"]]]] <- body_tqe_standardize(encoded$z$logT, f$body_tqe_scale)
            }
            else if (id == "fat_free_mass") {
                upstream <- vapply(members[1:2], `[[`, "", "column")
                previous <- c(upstream, internal[["logT"]])
                f <- body_tqe_coordinate_fit(d, encoded$z$logitq, base, previous, w, "logitq", m)
                d[[internal[["logitq"]]]] <- body_tqe_standardize(encoded$z$logitq, f$body_tqe_scale)
            }
            else if (id == "fat_mass") {
                upstream <- vapply(members[1:2], `[[`, "", "column")
                f <- body_tqe_e_fit(d, encoded, base, c(upstream, internal[c("logT", "logitq")]), w, m)
            }
            else f <- body_tqe_raw_fit(d, m, base, original_previous, w, columns, continuous_rule, count_rule)
            models[[id]] <- f
            warnings <- c(warnings, f$warnings)
            original_previous <- c(original_previous, m$column)
            checkpoint <- getOption("c39.a_checkpoint_callback")
            if (is.function(checkpoint)) 
                checkpoint(id, f, j)
        }
    list(models = models, members = members, warnings = unique(warnings), conditional_order = names(models), continuous_rule = continuous_rule, 
        model_revision = body_tqe_revision(), body_tqe = list(revision = body_tqe_revision(), prefix_binding = binding, 
            precision = encoded$precision, roundtrip = encoded$roundtrip, uniform_slot_mapping = c(weight = "logT", 
                fat_free_mass = "logitq", fat_mass = "e"), raw_columns = columns, old_upstream_ids = body_tqe_ids()[1:2], 
            old_upstream_exact_reused = !is.null(upstream_provenance), upstream_provenance = upstream_provenance, 
            downstream_predictor_law = "UNCHANGED_RAW_W_L_F_LINEAR_SPAN_AFFINE_BASIS_ONLY", e_law = models$fat_mass$e_law, 
            no_height = TRUE, no_clipping_no_retry_no_reseeding = TRUE), residual_assumption = "Joint constrained T,q,e reconstruction; all other original response laws retained")
}

body_tqe_category_draw <- function(f, X, u, lower, upper, id, audit) {
    k <- length(f$values)
    n <- nrow(X)
    if (k == 1L) 
        p <- matrix(1, n, 1L)
    else if (k == 2L) {
        pp <- stats::plogis(as.numeric(X %*% f$coefficients))
        p <- cbind(1 - pp, pp)
    }
    else {
        lp <- cbind(0, X %*% t(f$coefficients))
        lp <- lp - apply(lp, 1L, max)
        p <- exp(lp)
        p <- p/rowSums(p)
    }
    legal <- outer(lower, f$values, `<`) & outer(upper, f$values, `>=`)
    eligible_n <- rowSums(legal)
    retained_mass <- rowSums(p * legal)
    if (is.function(audit)) 
        audit(id, eligible_n, retained_mass, retained_mass^2/rowSums((p * legal)^2))
    body_tqe_assert(all(is.finite(p) & p >= 0 & p <= 1) && all(retained_mass > 0), "BODY_TQE_e_CATEGORY_CONDITIONAL_SUPPORT_EMPTY_NO_RETRY")
    p <- p * legal/retained_mass
    cumulative <- p
    if (k > 1L) 
        for (j in 2:k) cumulative[, j] <- cumulative[, j - 1L] + p[, j]
    cumulative[, k] <- 1
    ix <- 1L + rowSums(u >= cumulative[, -k, drop = FALSE])
    value <- f$values[ix]
    body_tqe_assert(all(value > lower & value <= upper), "BODY_TQE_e_CATEGORY_DRAW_INVALID_NO_CLIPPING")
    value
}

body_tqe_draw_a <- function(a, C, b, u, audit = NULL) {
    if (is.null(a$body_tqe)) 
        return(.body_tqe_original_draw_a(a, C, b, u, audit))
    body_tqe_assert(identical(a$body_tqe$revision, body_tqe_revision()) && is.matrix(u) && nrow(u) == nrow(C) && 
        ncol(u) == length(a$models) && all(is.finite(u) & u >= 0 & u < 1) && length(b) %in% c(1L, nrow(C)) && all(is.finite(b) & 
        b >= 0 & b <= 1), "BODY_TQE_INVALID_JOINT_DRAW_INPUT")
    d <- C
    d$B <- b
    internal <- body_tqe_names()
    columns <- a$body_tqe$raw_columns
    T <- q <- NULL
    for (j in seq_along(a$models)) {
        f <- a$models[[j]]
        role <- f$body_tqe_role
        id <- f$member$indicator_id
        if (is.null(role)) {
            prepared <- body_tqe_affine_prepare(d, f$predecessor_columns, columns, reference = f$body_tqe_affine)
            single <- list(models = setNames(list(f), id), members = setNames(list(f$member), id))
            d <- .body_tqe_original_draw_a(single, prepared$data, b, u[, j, drop = FALSE], audit)
        }
        else {
            X <- c39_a_design(d, NULL, NULL, f$design_reference)$X
            if (role != "e") {
                value <- c39_residual_draw(f, X, u[, j], -Inf, Inf, FALSE, paste0("body_tqe__", role), audit)
                d[[internal[[role]]]] <- body_tqe_standardize(value, f$body_tqe_scale)
                if (role == "logT") {
                  T <- exp(value)
                  body_tqe_assert(all(is.finite(T) & T > 0), "BODY_TQE_TOTAL_EXP_INVALID_NO_CLIPPING")
                }
                else {
                  q <- stats::plogis(value)
                  F <- T * q
                  L <- T * stats::plogis(-value)
                  body_tqe_assert(all(is.finite(F) & F > 0 & is.finite(L) & L > 0), "BODY_TQE_FRACTION_INVERSE_NONPOSITIVE_NO_CLIPPING")
                  d[[columns[["fat_mass"]]]] <- F
                  d[[columns[["fat_free_mass"]]]] <- L
                  T <- F + L
                }
            }
            else {
                Wmember <- a$members$weight
                lower <- pmax(-T, Wmember$lower - T)
                upper <- Wmember$upper - T
                if (identical(f$kind, "body_tqe_category")) 
                  e <- body_tqe_category_draw(f, X, u[, j], lower, upper, "body_tqe__e", audit)
                else {
                  f$rule <- "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL"
                  e <- c39_residual_draw(f, X, u[, j], lower, upper, TRUE, "body_tqe__e", audit)
                }
                d[[columns[["weight"]]]] <- T + e
                body_tqe_assert(all(is.finite(d[[columns[["weight"]]]]) & d[[columns[["weight"]]]] > 0), "BODY_TQE_RECONSTRUCTED_WEIGHT_NONPOSITIVE_NO_CLIPPING")
            }
        }
    }
    body_tqe_support(d, a$members)
    added <- setdiff(unname(internal), names(C))
    d[intersect(added, names(d))] <- NULL
    d
}

body_tqe_support <- function(d, members) {
    .body_tqe_original_support(d, members)
    if (all(body_tqe_mass_ids() %in% names(members))) {
        columns <- body_tqe_columns(members)
        body_tqe_assert(all(vapply(d[columns], function(x) all(is.finite(x) & x > 0), logical(1))), "BODY_TQE_RAW_POSITIVE_SUPPORT_FAILURE_NO_REPAIR")
        T <- d[[columns[["fat_mass"]]]] + d[[columns[["fat_free_mass"]]]]
        e <- d[[columns[["weight"]]]] - T
        body_tqe_assert(all(is.finite(T) & T > 0 & is.finite(e) & T + e > 0), "BODY_TQE_RECONSTRUCTION_SUPPORT_FAILURE_NO_REPAIR")
    }
    invisible(TRUE)
}

body_tqe_install <- function(engine) {
    body_tqe_assert(is.environment(engine) && !identical(engine, .GlobalEnv) && is.function(engine$c39_fit_a) && 
        is.function(engine$c39_draw_a) && is.function(engine$c39_support) && is.function(engine$c39_linear_fit) && 
        exists(".c39_loga_installed", engine, inherits = FALSE), "BODY_TQE_INSTALL_MIXED_A_ISOLATED_ENGINE_FIRST")
    if (exists(".body_tqe_installed", engine, inherits = FALSE)) {
        body_tqe_assert(identical(engine$.body_tqe_installed, body_tqe_revision()) && identical(engine$c39_fit_a, 
            engine$body_tqe_fit_a) && identical(engine$c39_draw_a, engine$body_tqe_draw_a), "BODY_TQE_INSTALLATION_DRIFT")
        return(invisible(engine))
    }
    engine$.body_tqe_original_fit_a <- engine$c39_fit_a
    engine$.body_tqe_original_draw_a <- engine$c39_draw_a
    engine$.body_tqe_original_support <- engine$c39_support
    source_environment <- environment(body_tqe_install)
    functions <- ls(source_environment, pattern = "^body_tqe_")
    for (name in functions) if (is.function(get(name, source_environment, inherits = FALSE))) {
        f <- get(name, source_environment, inherits = FALSE)
        environment(f) <- engine
        assign(name, f, envir = engine)
    }
    engine$c39_fit_a <- engine$body_tqe_fit_a
    engine$c39_draw_a <- engine$body_tqe_draw_a
    engine$c39_support <- engine$body_tqe_support
    engine$fd_fit_a <- engine$c39_fit_a
    engine$fd_draw_a <- engine$c39_draw_a
    engine$fd_support <- engine$c39_support
    engine$.body_tqe_installed <- body_tqe_revision()
    invisible(engine)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        body_tqe_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
        body_tqe_revision <- function() "SELECTED44_BODY_TQE_V1"
        body_tqe_ids <- function() c("smoking", "bmi", "weight", "fat_free_mass", "fat_mass", "hip", "waist", "clinical__hdl", 
            "clinical__tg", "clinical__glucose", "clinical__hba1c")
        body_tqe_mass_ids <- function() c("weight", "fat_free_mass", "fat_mass")
        body_tqe_names <- function() c(logT = ".body_tqe_logT", logitq = ".body_tqe_logitq", e = ".body_tqe_e", 
            affineT = ".body_tqe_affineT", affineF = ".body_tqe_affineF", affinee = ".body_tqe_affinee")
        body_tqe_hash <- function(x) {
            body_tqe_assert(requireNamespace("digest", quietly = TRUE), "BODY_TQE_DIGEST_REQUIRED")
            digest::digest(x, algo = "sha256", serialize = TRUE, serializeVersion = 2L)
        }
        body_tqe_binding <- function(d, members, base, w) {
            cols <- vapply(members[seq_len(11L)], `[[`, "", "column")
            body_tqe_hash(list(revision = body_tqe_revision(), members = members[seq_len(11L)], base = base, data = d[unique(c(base, 
                cols))], eid = if ("eid" %in% names(d)) as.character(d$eid), weights = as.numeric(w)))
        }
        body_tqe_precision <- function(x) {
            tolerance <- 32 * .Machine$double.eps * pmax(1, abs(x))
            for (k in 0:10) if (all(abs(x - round(x, k)) <= tolerance)) return(k)
            NA_integer_
        }
        body_tqe_columns <- function(members) vapply(members[body_tqe_mass_ids()], `[[`, "", "column")
        body_tqe_encode <- function(d, members) {
            cols <- body_tqe_columns(members)
            W <- d[[cols[["weight"]]]]
            L <- d[[cols[["fat_free_mass"]]]]
            F <- d[[cols[["fat_mass"]]]]
            body_tqe_assert(all(is.finite(W) & W > 0 & is.finite(L) & L > 0 & is.finite(F) & F > 0), "BODY_TQE_STRICTLY_POSITIVE_RAW_W_F_L_REQUIRED_NO_REPAIR")
            T <- F + L
            q <- F/T
            e <- W - T
            body_tqe_assert(all(is.finite(T) & T > 0 & is.finite(q) & q > 0 & q < 1), "BODY_TQE_INVALID_TOTAL_OR_FRACTION_NO_REPAIR")
            z <- list(logT = log(T), logitq = log(F) - log(L), e = e)
            F2 <- exp(z$logT) * stats::plogis(z$logitq)
            L2 <- exp(z$logT) * stats::plogis(-z$logitq)
            W2 <- F2 + L2 + e
            error <- c(W = max(abs(W2 - W)/pmax(1, abs(W))), F = max(abs(F2 - F)/pmax(1, abs(F))), L = max(abs(L2 - 
                L)/pmax(1, abs(L))))
            body_tqe_assert(all(is.finite(error) & error <= 1e-10), "BODY_TQE_RAW_ROUNDTRIP_LOSES_RECORDED_VARIATION")
            digits <- vapply(list(W, L, F), body_tqe_precision, integer(1))
            precision <- if (anyNA(digits)) NA_integer_ else max(digits)
            states <- if (is.na(precision)) NULL else round(e * 10^precision)/10^precision
            if (!is.null(states)) states[states == 0] <- 0
            state_error <- if (is.null(states)) NA_real_ else max(abs(states - e)/pmax(1, abs(W)))
            if (!is.null(states)) body_tqe_assert(state_error <= 1e-10, "BODY_TQE_DECIMAL_STATE_COLLAPSE_LOSES_RECORDED_VARIATION")
            list(z = z, total = T, fraction = q, columns = cols, e_states = states, precision = list(raw_decimal_digits = digits, 
                e_decimal_digits = precision, e_binary_unique_n = length(unique(e)), e_state_n = if (is.null(states)) NULL else length(unique(states)), 
                max_scaled_state_residue = state_error, arithmetic_tolerance = "32*double.eps*pmax(1,abs(raw)); state roundtrip <=1e-10"), 
                roundtrip = list(max_scaled_error = error, tolerance = 1e-10))
        }
        body_tqe_scale <- function(x, w) {
            center <- sum(w * x)/sum(w)
            scale <- sqrt(sum(w * (x - center)^2)/sum(w))
            list(center = center, scale = if (scale == 0) 1 else scale)
        }
        body_tqe_standardize <- function(x, ref) (x - ref$center)/ref$scale
        body_tqe_reference_copy <- function(reference) {
            body_tqe_assert(is.list(reference) && inherits(reference$formula, "formula"), "BODY_TQE_UPSTREAM_FORMULA_REFERENCE_REQUIRED")
            formula <- reference$formula
            environment(formula) <- baseenv()
            reference$formula <- formula
            reference
        }
        body_tqe_affine_prepare <- function(d, previous, columns, w = NULL, reference = NULL) {
            if (!all(columns %in% previous)) return(list(data = d, previous = previous, reference = NULL))
            keys <- body_tqe_names()[c("affineT", "affineF", "affinee")]
            W <- d[[columns[["weight"]]]]
            L <- d[[columns[["fat_free_mass"]]]]
            F <- d[[columns[["fat_mass"]]]]
            values <- list(F + L, F, W - (F + L))
            if (is.null(reference)) reference <- list(columns = columns, keys = keys, scale = lapply(values, body_tqe_scale, 
                w = w), basis = "INVERTIBLE_AFFINE_T_F_e_SAME_RAW_W_L_F_PREDICTOR_SPAN")
            body_tqe_assert(identical(reference$columns, columns) && identical(reference$keys, keys), "BODY_TQE_AFFINE_PREDICTOR_REFERENCE_DRIFT")
            for (j in seq_along(keys)) d[[keys[[j]]]] <- body_tqe_standardize(values[[j]], reference$scale[[j]])
            previous[match(columns, previous)] <- unname(keys)
            list(data = d, previous = previous, reference = reference)
        }
        body_tqe_raw_fit <- function(d, member, base, previous, w, columns, continuous_rule, count_rule) {
            prepared <- body_tqe_affine_prepare(d, previous, columns, w)
            design <- c39_a_design(prepared$data, base, prepared$previous)
            terms <- attr(stats::terms(design$reference$formula), "term.labels")
            callback <- getOption("c39.a_checkpoint_callback")
            options(c39.a_checkpoint_callback = NULL)
            on.exit(options(c39.a_checkpoint_callback = callback), add = TRUE)
            a <- .body_tqe_original_fit_a(prepared$data, setNames(list(member), member$indicator_id), terms, w, 
                continuous_rule, count_rule)
            f <- a$models[[1L]]
            f$design_reference <- design$reference
            f$formula <- design$reference$formula
            f$predecessor_columns <- previous
            f$body_tqe_affine <- prepared$reference
            f
        }
        body_tqe_coordinate_fit <- function(d, y, base, previous, w, id, member) {
            design <- c39_a_design(d, base, previous)
            f <- c39_linear_fit(design$X, y, w, paste0("body_tqe__", id), "EMPIRICAL_RESIDUAL_UNTRUNCATED")
            f$design_reference <- design$reference
            f$columns <- colnames(design$X)
            f$formula <- design$reference$formula
            f$predecessor_columns <- previous
            f$member <- member
            f$body_tqe_role <- id
            f$body_tqe_scale <- body_tqe_scale(y, w)
            f
        }
        body_tqe_e_fit <- function(d, encoded, base, previous, w, member) {
            design <- c39_a_design(d, base, previous)
            states <- encoded$e_states
            if (is.null(states)) {
                f <- body_tqe_coordinate_fit(d, encoded$z$e, base, previous, w, "e", member)
                f$e_law <- "WEIGHTED_EMPIRICAL_RESIDUAL_CONDITIONAL_ON_W_POSITIVE"
                return(f)
            }
            values <- sort(unique(states))
            body_tqe_assert(length(values) <= 256L, "BODY_TQE_MEASURED_DISCRETE_e_OVER_256_STATES_REQUIRES_REVIEW")
            if (length(values) == 1L) f <- list(kind = "body_tqe_category", values = values, coefficients = NULL, 
                warnings = character()) else {
                X <- design$X
                y <- factor(match(states, values), levels = seq_along(values))
                body_tqe_assert(requireNamespace("nnet", quietly = TRUE) && qr(X)$rank == ncol(X), "BODY_TQE_e_CATEGORY_DESIGN_RANK_DEFICIENT")
                warnings <- character()
                fit <- withCallingHandlers(nnet::multinom(y ~ X - 1, weights = w, trace = FALSE, maxit = 1000L, 
                  MaxNWts = 100000L, model = FALSE), warning = function(e) {
                  warnings <<- c(warnings, conditionMessage(e))
                  invokeRestart("muffleWarning")
                })
                body_tqe_assert(identical(fit$convergence, 0L) && all(is.finite(stats::coef(fit))), "BODY_TQE_e_CATEGORY_NOT_CONVERGED")
                f <- list(kind = "body_tqe_category", values = values, coefficients = stats::coef(fit), warnings = unique(warnings), 
                  convergence = fit$convergence)
            }
            f$design_reference <- design$reference
            f$columns <- colnames(design$X)
            f$formula <- design$reference$formula
            f$predecessor_columns <- previous
            f$member <- member
            f$body_tqe_role <- "e"
            f$e_law <- "MEASURED_DECIMAL_STATES_WEIGHTED_CONDITIONAL_MULTINOMIAL"
            f$body_tqe_scale <- body_tqe_scale(encoded$z$e, w)
            f
        }
        body_tqe_fit_a <- function(d, members, base, w, continuous_rule = "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL", 
            count_rule = "FIXED_SUPPORT_MULTINOMIAL_0_7", prefix_fit = getOption("body_tqe.prefix_A"), upstream_fit = getOption("body_tqe.upstream_A")) {
            ids <- vapply(members, `[[`, "", "indicator_id")
            body_tqe_assert(length(ids) >= 11L && identical(unname(ids[seq_len(11L)]), body_tqe_ids()) && identical(names(members), 
                unname(ids)), "BODY_TQE_FROZEN_11_PREFIX_REQUIRED")
            body_tqe_assert(length(w) == nrow(d) && all(is.finite(w) & w > 0) && is.finite(sum(w)), "BODY_TQE_EXACT_SELECTION_WEIGHTS_INVALID")
            body_tqe_assert(!any(body_tqe_names() %in% names(d)), "BODY_TQE_RESERVED_COLUMNS_EXIST")
            .body_tqe_original_support(d, members)
            binding <- body_tqe_binding(d, members, base, w)
            encoded <- body_tqe_encode(d, members)
            columns <- encoded$columns
            internal <- body_tqe_names()
            original_previous <- character()
            models <- list()
            warnings <- character()
            start <- 1L
            upstream_provenance <- NULL
            if (!is.null(prefix_fit)) {
                body_tqe_assert(identical(prefix_fit$body_tqe$revision, body_tqe_revision()) && identical(prefix_fit$body_tqe$prefix_binding, 
                  binding) && length(prefix_fit$members) == 11L && identical(prefix_fit$members, members[seq_len(11L)]), 
                  "BODY_TQE_PREFIX_SAMPLE_WEIGHTS_MEMBERS_BINDING_MISMATCH")
                models <- prefix_fit$models
                warnings <- prefix_fit$warnings
                original_previous <- vapply(members[seq_len(11L)], `[[`, "", "column")
                start <- 12L
                upstream_provenance <- prefix_fit$body_tqe$upstream_provenance
            } else if (!is.null(upstream_fit)) {
                upstream_ids <- body_tqe_ids()[1:2]
                body_tqe_assert(identical(upstream_fit$prefix_binding, binding) && identical(upstream_fit$members, 
                  members[1:2]) && identical(names(upstream_fit$models), upstream_ids) && is.list(upstream_fit$provenance) && 
                  length(upstream_fit$provenance) > 0L, "BODY_TQE_OLD_UPSTREAM_INPUT_WEIGHTS_MEMBERS_PROVENANCE_MISMATCH")
                before <- character()
                for (id in upstream_ids) {
                  f <- upstream_fit$models[[id]]
                  expected_design <- c39_a_design(d, base, before)
                  expected <- expected_design$reference
                  old_design <- c39_a_design(d, NULL, NULL, f$design_reference)
                  body_tqe_assert(identical(f$member, members[[id]]) && identical(body_tqe_reference_copy(f$design_reference), 
                    body_tqe_reference_copy(expected)) && identical(old_design$X, expected_design$X) && identical(f$predecessor_columns, 
                    before) && is.null(f$body_tqe_role), paste("BODY_TQE_OLD_UPSTREAM_CONDITIONAL_DESIGN_CHANGED", 
                    id))
                  before <- c(before, members[[id]]$column)
                }
                body_tqe_assert(identical(upstream_fit$models$bmi$log_a_revision$revision, "CORE39_BMI_LOG_FAT_LINEAR_V1"), 
                  "BODY_TQE_OLD_UPSTREAM_APPROVED_LOG_BMI_REQUIRED")
                models <- upstream_fit$models
                warnings <- unique(unlist(lapply(models, `[[`, "warnings")))
                original_previous <- before
                start <- 3L
                upstream_provenance <- upstream_fit$provenance
            }
            if (start <= length(members)) for (j in seq.int(start, length(members))) {
                m <- members[[j]]
                id <- m$indicator_id
                if (id == "weight") {
                  f <- body_tqe_coordinate_fit(d, encoded$z$logT, base, original_previous, w, "logT", m)
                  d[[internal[["logT"]]]] <- body_tqe_standardize(encoded$z$logT, f$body_tqe_scale)
                } else if (id == "fat_free_mass") {
                  upstream <- vapply(members[1:2], `[[`, "", "column")
                  previous <- c(upstream, internal[["logT"]])
                  f <- body_tqe_coordinate_fit(d, encoded$z$logitq, base, previous, w, "logitq", m)
                  d[[internal[["logitq"]]]] <- body_tqe_standardize(encoded$z$logitq, f$body_tqe_scale)
                } else if (id == "fat_mass") {
                  upstream <- vapply(members[1:2], `[[`, "", "column")
                  f <- body_tqe_e_fit(d, encoded, base, c(upstream, internal[c("logT", "logitq")]), w, m)
                } else f <- body_tqe_raw_fit(d, m, base, original_previous, w, columns, continuous_rule, count_rule)
                models[[id]] <- f
                warnings <- c(warnings, f$warnings)
                original_previous <- c(original_previous, m$column)
                checkpoint <- getOption("c39.a_checkpoint_callback")
                if (is.function(checkpoint)) checkpoint(id, f, j)
            }
            list(models = models, members = members, warnings = unique(warnings), conditional_order = names(models), 
                continuous_rule = continuous_rule, model_revision = body_tqe_revision(), body_tqe = list(revision = body_tqe_revision(), 
                  prefix_binding = binding, precision = encoded$precision, roundtrip = encoded$roundtrip, uniform_slot_mapping = c(weight = "logT", 
                    fat_free_mass = "logitq", fat_mass = "e"), raw_columns = columns, old_upstream_ids = body_tqe_ids()[1:2], 
                  old_upstream_exact_reused = !is.null(upstream_provenance), upstream_provenance = upstream_provenance, 
                  downstream_predictor_law = "UNCHANGED_RAW_W_L_F_LINEAR_SPAN_AFFINE_BASIS_ONLY", e_law = models$fat_mass$e_law, 
                  no_height = TRUE, no_clipping_no_retry_no_reseeding = TRUE), residual_assumption = "Joint constrained T,q,e reconstruction; all other original response laws retained")
        }
        body_tqe_category_draw <- function(f, X, u, lower, upper, id, audit) {
            k <- length(f$values)
            n <- nrow(X)
            if (k == 1L) p <- matrix(1, n, 1L) else if (k == 2L) {
                pp <- stats::plogis(as.numeric(X %*% f$coefficients))
                p <- cbind(1 - pp, pp)
            } else {
                lp <- cbind(0, X %*% t(f$coefficients))
                lp <- lp - apply(lp, 1L, max)
                p <- exp(lp)
                p <- p/rowSums(p)
            }
            legal <- outer(lower, f$values, `<`) & outer(upper, f$values, `>=`)
            eligible_n <- rowSums(legal)
            retained_mass <- rowSums(p * legal)
            if (is.function(audit)) audit(id, eligible_n, retained_mass, retained_mass^2/rowSums((p * legal)^2))
            body_tqe_assert(all(is.finite(p) & p >= 0 & p <= 1) && all(retained_mass > 0), "BODY_TQE_e_CATEGORY_CONDITIONAL_SUPPORT_EMPTY_NO_RETRY")
            p <- p * legal/retained_mass
            cumulative <- p
            if (k > 1L) for (j in 2:k) cumulative[, j] <- cumulative[, j - 1L] + p[, j]
            cumulative[, k] <- 1
            ix <- 1L + rowSums(u >= cumulative[, -k, drop = FALSE])
            value <- f$values[ix]
            body_tqe_assert(all(value > lower & value <= upper), "BODY_TQE_e_CATEGORY_DRAW_INVALID_NO_CLIPPING")
            value
        }
        body_tqe_draw_a <- function(a, C, b, u, audit = NULL) {
            if (is.null(a$body_tqe)) return(.body_tqe_original_draw_a(a, C, b, u, audit))
            body_tqe_assert(identical(a$body_tqe$revision, body_tqe_revision()) && is.matrix(u) && nrow(u) == nrow(C) && 
                ncol(u) == length(a$models) && all(is.finite(u) & u >= 0 & u < 1) && length(b) %in% c(1L, nrow(C)) && 
                all(is.finite(b) & b >= 0 & b <= 1), "BODY_TQE_INVALID_JOINT_DRAW_INPUT")
            d <- C
            d$B <- b
            internal <- body_tqe_names()
            columns <- a$body_tqe$raw_columns
            T <- q <- NULL
            for (j in seq_along(a$models)) {
                f <- a$models[[j]]
                role <- f$body_tqe_role
                id <- f$member$indicator_id
                if (is.null(role)) {
                  prepared <- body_tqe_affine_prepare(d, f$predecessor_columns, columns, reference = f$body_tqe_affine)
                  single <- list(models = setNames(list(f), id), members = setNames(list(f$member), id))
                  d <- .body_tqe_original_draw_a(single, prepared$data, b, u[, j, drop = FALSE], audit)
                } else {
                  X <- c39_a_design(d, NULL, NULL, f$design_reference)$X
                  if (role != "e") {
                    value <- c39_residual_draw(f, X, u[, j], -Inf, Inf, FALSE, paste0("body_tqe__", role), audit)
                    d[[internal[[role]]]] <- body_tqe_standardize(value, f$body_tqe_scale)
                    if (role == "logT") {
                      T <- exp(value)
                      body_tqe_assert(all(is.finite(T) & T > 0), "BODY_TQE_TOTAL_EXP_INVALID_NO_CLIPPING")
                    } else {
                      q <- stats::plogis(value)
                      F <- T * q
                      L <- T * stats::plogis(-value)
                      body_tqe_assert(all(is.finite(F) & F > 0 & is.finite(L) & L > 0), "BODY_TQE_FRACTION_INVERSE_NONPOSITIVE_NO_CLIPPING")
                      d[[columns[["fat_mass"]]]] <- F
                      d[[columns[["fat_free_mass"]]]] <- L
                      T <- F + L
                    }
                  } else {
                    Wmember <- a$members$weight
                    lower <- pmax(-T, Wmember$lower - T)
                    upper <- Wmember$upper - T
                    if (identical(f$kind, "body_tqe_category")) e <- body_tqe_category_draw(f, X, u[, j], lower, 
                      upper, "body_tqe__e", audit) else {
                      f$rule <- "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL"
                      e <- c39_residual_draw(f, X, u[, j], lower, upper, TRUE, "body_tqe__e", audit)
                    }
                    d[[columns[["weight"]]]] <- T + e
                    body_tqe_assert(all(is.finite(d[[columns[["weight"]]]]) & d[[columns[["weight"]]]] > 0), "BODY_TQE_RECONSTRUCTED_WEIGHT_NONPOSITIVE_NO_CLIPPING")
                  }
                }
            }
            body_tqe_support(d, a$members)
            added <- setdiff(unname(internal), names(C))
            d[intersect(added, names(d))] <- NULL
            d
        }
        body_tqe_support <- function(d, members) {
            .body_tqe_original_support(d, members)
            if (all(body_tqe_mass_ids() %in% names(members))) {
                columns <- body_tqe_columns(members)
                body_tqe_assert(all(vapply(d[columns], function(x) all(is.finite(x) & x > 0), logical(1))), "BODY_TQE_RAW_POSITIVE_SUPPORT_FAILURE_NO_REPAIR")
                T <- d[[columns[["fat_mass"]]]] + d[[columns[["fat_free_mass"]]]]
                e <- d[[columns[["weight"]]]] - T
                body_tqe_assert(all(is.finite(T) & T > 0 & is.finite(e) & T + e > 0), "BODY_TQE_RECONSTRUCTION_SUPPORT_FAILURE_NO_REPAIR")
            }
            invisible(TRUE)
        }
        body_tqe_install <- function(engine) {
            body_tqe_assert(is.environment(engine) && !identical(engine, .GlobalEnv) && is.function(engine$c39_fit_a) && 
                is.function(engine$c39_draw_a) && is.function(engine$c39_support) && is.function(engine$c39_linear_fit) && 
                exists(".c39_loga_installed", engine, inherits = FALSE), "BODY_TQE_INSTALL_MIXED_A_ISOLATED_ENGINE_FIRST")
            if (exists(".body_tqe_installed", engine, inherits = FALSE)) {
                body_tqe_assert(identical(engine$.body_tqe_installed, body_tqe_revision()) && identical(engine$c39_fit_a, 
                  engine$body_tqe_fit_a) && identical(engine$c39_draw_a, engine$body_tqe_draw_a), "BODY_TQE_INSTALLATION_DRIFT")
                return(invisible(engine))
            }
            engine$.body_tqe_original_fit_a <- engine$c39_fit_a
            engine$.body_tqe_original_draw_a <- engine$c39_draw_a
            engine$.body_tqe_original_support <- engine$c39_support
            source_environment <- environment(body_tqe_install)
            functions <- ls(source_environment, pattern = "^body_tqe_")
            for (name in functions) if (is.function(get(name, source_environment, inherits = FALSE))) {
                f <- get(name, source_environment, inherits = FALSE)
                environment(f) <- engine
                assign(name, f, envir = engine)
            }
            engine$c39_fit_a <- engine$body_tqe_fit_a
            engine$c39_draw_a <- engine$body_tqe_draw_a
            engine$c39_support <- engine$body_tqe_support
            engine$fd_fit_a <- engine$c39_fit_a
            engine$fd_draw_a <- engine$c39_draw_a
            engine$fd_support <- engine$c39_support
            engine$.body_tqe_installed <- body_tqe_revision()
            invisible(engine)
        }
    }), envir = context)
}

