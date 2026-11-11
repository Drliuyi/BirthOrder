# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

nmrA_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

nmrA_revision <- function() "SELECTED44_BASE11_RAW29_COMPOSITION_RECONSTRUCTION_V1"

nmrA_rng <- function() list(kind = RNGkind(), present = exists(".Random.seed", .GlobalEnv, inherits = FALSE), seed = if (exists(".Random.seed", 
    .GlobalEnv, inherits = FALSE)) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL)

nmrA_eval <- function(expression, values, n) {
    walk <- function(x) {
        if (is.numeric(x)) 
            return(x)
        if (is.symbol(x)) {
            key <- as.character(x)
            nmrA_assert(key %in% names(values), paste("NMR_EXPRESSION_UNKNOWN_SYMBOL", key))
            return(values[[key]])
        }
        nmrA_assert(is.call(x), "NMR_EXPRESSION_LANGUAGE_INVALID")
        op <- as.character(x[[1L]])
        nmrA_assert(op %in% c("+", "-", "*", "/", "log", "exp", "("), paste("NMR_EXPRESSION_OPERATION_NOT_ALLOWED", 
            op))
        args <- lapply(as.list(x)[-1L], walk)
        nmrA_assert(length(args) %in% if (op %in% c("+", "-")) 
            1:2
        else if (op %in% c("*", "/")) 
            2L
        else 1L, "NMR_EXPRESSION_ARITY_INVALID")
        do.call(switch(op, `+` = `+`, `-` = `-`, `*` = `*`, `/` = `/`, log = log, exp = exp, `(` = identity), args)
    }
    ast <- parse(text = expression, keep.source = FALSE)
    nmrA_assert(length(ast) == 1L, "NMR_EXPRESSION_SINGLE_EXPRESSION_REQUIRED")
    answer <- walk(ast[[1L]])
    nmrA_assert(is.numeric(answer) && length(answer) %in% c(1L, n) && all(is.finite(answer)), paste("NMR_EXPRESSION_NONFINITE", 
        expression))
    rep_len(answer, n)
}

nmrA_basis <- function(D) {
    nmrA_assert(length(D) == 1L && D >= 1L, "NMR_ILR_DIMENSION_INVALID")
    if (D == 1L) 
        return(matrix(numeric(), 1L, 0L))
    H <- stats::contr.helmert(D)
    sweep(H, 2L, sqrt(colSums(H^2)), `/`)
}

nmrA_ilr <- function(parts) {
    nmrA_assert(is.matrix(parts) && is.numeric(parts) && ncol(parts) >= 1L && all(is.finite(parts) & parts > 0), 
        "NMR_ILR_STRICTLY_POSITIVE_PARTS_REQUIRED_NO_EPSILON")
    log(parts) %*% nmrA_basis(ncol(parts))
}

nmrA_inverse_ilr <- function(z, total, D = ncol(z) + 1L) {
    nmrA_assert(is.matrix(z) && ncol(z) == D - 1L && all(is.finite(z)) && length(total) %in% c(1L, nrow(z)) && 
        all(is.finite(total) & total > 0), "NMR_INVERSE_ILR_INVALID_INPUT")
    lp <- z %*% t(nmrA_basis(D))
    lp <- lp - apply(lp, 1L, max)
    p <- exp(lp)
    answer <- p/rowSums(p) * rep_len(total, nrow(z))
    nmrA_assert(all(is.finite(answer) & answer > 0), "NMR_INVERSE_ILR_UNDERFLOW_OR_OVERFLOW_NO_CLIPPING")
    answer
}

nmrA_raw_values <- function(data, members) {
    fields <- sub("^nmr__nmr_", "F", names(members))
    cols <- vapply(members, `[[`, "", "column")
    nmrA_assert(all(cols %in% names(data)) && !anyDuplicated(fields) && all(vapply(data[cols], is.numeric, logical(1))), 
        "NMR_RAW_NUMERIC_COLUMNS_REQUIRED")
    setNames(as.list(data[cols]), fields)
}

nmrA_block_fields <- function(block) if (identical(block$type, "scalar")) as.character(block$field) else names(block$raw)

nmrA_block_encode <- function(block, raw, n) {
    nmrA_assert(identical(block$type, "composition"), "NMR_COMPOSITION_BLOCK_REQUIRED")
    part_values <- list()
    values <- raw
    for (key in names(block$parts)) {
        part_values[[key]] <- nmrA_eval(block$parts[[key]], values, n)
        values[[key]] <- part_values[[key]]
    }
    parts <- do.call(cbind, part_values)
    colnames(parts) <- names(block$parts)
    nmrA_assert(ncol(parts) >= 2L && !anyDuplicated(colnames(parts)) && all(is.finite(parts) & parts >= 0), paste("NMR_NEGATIVE_BIOCHEMICAL_PART_NO_REPAIR", 
        block$id))
    total <- rowSums(parts)
    if (!is.null(block$fixed_total)) 
        nmrA_assert(all(abs(total - block$fixed_total) <= 512 * .Machine$double.eps * pmax(1, abs(total))), paste("NMR_SOURCE_FIXED_TOTAL_IDENTITY_FAILED", 
            block$id))
    masks <- apply(parts > 0, 1L, function(x) paste0(as.integer(x), collapse = ""))
    aux <- list()
    values <- c(raw, as.list(as.data.frame(parts)), list(TOTAL = total))
    for (key in names(block$aux)) {
        specification <- block$aux[[key]]
        expression <- if (is.character(specification)) 
            specification
        else specification$expression
        aux[[key]] <- nmrA_eval(expression, values, n)
        values[[key]] <- aux[[key]]
    }
    list(parts = parts, total = total, mask = masks, aux = aux)
}

nmrA_block_decode <- function(block, parts, aux, prior = list()) {
    n <- nrow(parts)
    values <- c(as.list(as.data.frame(parts)), list(TOTAL = rowSums(parts)), aux, prior)
    setNames(lapply(block$raw, nmrA_eval, values = values, n = n), paste0("F", names(block$raw)))
}

nmrA_roundtrip <- function(block, raw, n) {
    encoded <- nmrA_block_encode(block, raw, n)
    decoded <- nmrA_block_decode(block, encoded$parts, encoded$aux, raw)
    discrepancies <- vapply(names(decoded), function(key) max(abs(decoded[[key]] - raw[[key]])/pmax(1, abs(raw[[key]]))), 
        numeric(1))
    tolerance <- 1e-10
    nmrA_assert(all(discrepancies <= tolerance), paste("NMR_RAW_ROUNDTRIP_LOSES_RECORDED_VARIATION", block$id, 
        paste(names(discrepancies)[discrepancies > tolerance], collapse = ",")))
    list(encoded = encoded, max_scaled_error = max(discrepancies), arithmetic_scaled_tolerance = tolerance, exact_recorded_derived_identities_claimed = FALSE)
}

nmrA_status_fit <- function(X, status, weights, id) {
    levels <- sort(unique(status))
    y <- factor(status, levels = levels)
    if (length(levels) == 1L) 
        return(list(kind = "constant", levels = levels))
    nmrA_assert(requireNamespace("nnet", quietly = TRUE) && qr(X)$rank == ncol(X), paste("NMR_ZERO_MASK_DESIGN_RANK_DEFICIENT", 
        id))
    warnings <- character()
    fit <- withCallingHandlers(nnet::multinom(y ~ X - 1, weights = weights, trace = FALSE, maxit = 1000L, MaxNWts = 100000L, 
        model = FALSE), warning = function(e) {
        warnings <<- c(warnings, conditionMessage(e))
        invokeRestart("muffleWarning")
    })
    coef <- stats::coef(fit)
    nmrA_assert(identical(fit$convergence, 0L) && all(is.finite(coef)), paste("NMR_ZERO_MASK_MODEL_NOT_CONVERGED", 
        id))
    list(kind = "multinomial", levels = levels, coefficients = coef, warnings = unique(warnings))
}

nmrA_status_draw <- function(fit, X, u) {
    if (identical(fit$kind, "constant")) 
        return(rep(fit$levels, nrow(X)))
    if (length(fit$levels) == 2L) {
        pp <- stats::plogis(as.numeric(X %*% fit$coefficients))
        p <- cbind(1 - pp, pp)
    }
    else {
        lp <- cbind(0, X %*% t(fit$coefficients))
        lp <- lp - apply(lp, 1L, max)
        p <- exp(lp)
        p <- p/rowSums(p)
    }
    nmrA_assert(all(is.finite(p) & p >= 0 & p <= 1), "NMR_ZERO_MASK_PROBABILITY_INVALID")
    cumulative <- p
    if (ncol(p) > 1L) 
        for (j in 2:ncol(p)) cumulative[, j] <- cumulative[, j - 1L] + p[, j]
    fit$levels[1L + rowSums(u >= cumulative[, -ncol(p), drop = FALSE])]
}

nmrA_conditional_fit <- function(data, y, base, previous, weights, engine, id, lower = -Inf, upper = Inf, rule = "EMPIRICAL_RESIDUAL_UNTRUNCATED") {
    design <- engine$c39_a_design(data, base, previous)
    fit <- engine$c39_linear_fit(design$X, y, weights, id, rule)
    fit$design_reference <- design$reference
    fit$id <- id
    fit$lower <- lower
    fit$upper <- upper
    fit$predictor_center <- sum(weights * y)/sum(weights)
    fit$predictor_scale <- sqrt(sum(weights * (y - fit$predictor_center)^2)/sum(weights))
    if (fit$predictor_scale == 0) 
        fit$predictor_scale <- 1
    fit
}

nmrA_predictor <- function(fit, value) (value - fit$predictor_center)/fit$predictor_scale

nmrA_predictor_inverse <- function(fit, value) value * fit$predictor_scale + fit$predictor_center

nmrA_conditional_draw <- function(fit, data, u, engine, audit = NULL, lower = fit$lower, upper = fit$upper) {
    X <- engine$c39_a_design(data, NULL, NULL, fit$design_reference)$X
    f <- fit
    if (any(is.finite(lower)) || any(is.finite(upper))) 
        f$rule <- "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL"
    engine$c39_residual_draw(f, X, u, lower, upper, FALSE, fit$id, audit)
}

nmrA_member_bounds <- function(field, members, raw_support = NULL) {
    m <- members[[paste0("nmr__nmr_", field)]]
    support <- raw_support[[field]]
    c(max(m$lower, if (is.null(support$lower)) -Inf else support$lower), min(m$upper, if (is.null(support$upper)) Inf else support$upper))
}

nmrA_affine_bounds <- function(block, values0, values1, members, raw_support, variable, allowed_aux = character()) {
    n <- length(values0$TOTAL)
    lower <- rep(-Inf, n)
    upper <- rep(Inf, n)
    for (field in names(block$raw)) {
        expression <- block$raw[[field]]
        symbols <- all.names(parse(text = expression)[[1L]], functions = FALSE, unique = TRUE)
        unresolved <- setdiff(intersect(symbols, names(block$aux)), allowed_aux)
        if (length(unresolved)) 
            next
        a <- nmrA_eval(expression, values0, n)
        b <- nmrA_eval(expression, values1, n) - a
        bounds <- nmrA_member_bounds(field, members, raw_support)
        positive <- b > 0
        negative <- b < 0
        fixed <- b == 0
        lower[positive] <- pmax(lower[positive], (bounds[1L] - a[positive])/b[positive])
        upper[positive] <- pmin(upper[positive], (bounds[2L] - a[positive])/b[positive])
        lower[negative] <- pmax(lower[negative], (bounds[2L] - a[negative])/b[negative])
        upper[negative] <- pmin(upper[negative], (bounds[1L] - a[negative])/b[negative])
        nmrA_assert(all(!fixed | (a >= bounds[1L] & a <= bounds[2L])), paste("NMR_RECONSTRUCTION_FIXED_SUPPORT_FAILURE", 
            block$id, field, variable))
    }
    nmrA_assert(all(lower <= upper), paste("NMR_DYNAMIC_PHYSICAL_SUPPORT_EMPTY", block$id, variable))
    list(lower = lower, upper = upper)
}

nmrA_total_bounds <- function(block, proportions, raw, members, raw_support) {
    n <- nrow(proportions)
    if (!identical(block$component_unit, "reported_percent")) 
        return(list(lower = rep(-Inf, n), upper = rep(Inf, n)))
    p0 <- proportions * 0
    values0 <- c(as.list(as.data.frame(p0)), list(TOTAL = rep(0, n)), raw)
    values1 <- c(as.list(as.data.frame(proportions)), list(TOTAL = rep(1, n)), raw)
    bounds <- nmrA_affine_bounds(block, values0, values1, members, raw_support, "TOTAL")
    nmrA_assert(all(bounds$upper > 0), "NMR_POSITIVE_TOTAL_PHYSICAL_SUPPORT_EMPTY")
    list(lower = ifelse(bounds$lower > 0, log(bounds$lower), -Inf), upper = log(bounds$upper))
}

nmrA_aux_dynamic_bounds <- function(block, key, parts, auxiliary, raw, members, raw_support) {
    n <- nrow(parts)
    common <- c(as.list(as.data.frame(parts)), list(TOTAL = rowSums(parts)), raw, auxiliary)
    v0 <- common
    v1 <- common
    v0[[key]] <- rep(0, n)
    v1[[key]] <- rep(1, n)
    nmrA_affine_bounds(block, v0, v1, members, raw_support, key, c(names(auxiliary), key))
}

nmrA_aux_bounds <- function(specification, y) {
    if (is.character(specification)) 
        return(c(-Inf, Inf))
    nmrA_assert(!identical(specification$support, "observed_range"), "NMR_OBSERVED_EXTREMA_ARE_NOT_PHYSICAL_MEASUREMENT_SUPPORT")
    c(if (is.null(specification$lower)) -Inf else as.numeric(specification$lower), if (is.null(specification$upper)) Inf else as.numeric(specification$upper))
}

fit_nmr_composition_A <- function(data, members, base, weights, engine, blockspec, core_fit = NULL, core_callback = NULL) {
    rng <- nmrA_rng()
    nmrA_assert(is.data.frame(data) && is.list(members) && is.environment(engine) && length(weights) == nrow(data) && 
        all(is.finite(weights) & weights > 0) && is.finite(sum(weights)), "NMR_A_INVALID_INPUT")
    nmrA_assert(!any(grepl("^\\.nmrA_", names(data))), "NMR_A_RESERVED_COORDINATE_COLUMN_EXISTS")
    nmrA_assert(identical(blockspec$schema, "NMR_COMPOSITION_BLOCKSPEC_V1") && length(blockspec$blocks) > 0L, "NMR_REVIEWED_BLOCKSPEC_REQUIRED")
    core_ids <- names(members)[seq_len(11L)]
    nmr_ids <- names(members)[-(seq_len(11L))]
    expected <- if (isTRUE(blockspec$synthetic)) 
        length(nmr_ids)
    else 29L
    nmrA_assert(length(core_ids) == 11L && length(nmr_ids) == expected && all(grepl("^nmr__nmr_[0-9]+$", nmr_ids)), 
        "NMR_A_SELECTED_BASE11_PLUS_FROZEN29_REQUIRED")
    covered <- unlist(lapply(blockspec$blocks, nmrA_block_fields), use.names = FALSE)
    nmrA_assert(!anyDuplicated(covered) && setequal(covered, sub("^nmr__nmr_", "", nmr_ids)), "NMR_BLOCKS_MUST_PARTITION_ALL_RAW_INFORMATION_NO_SELECTION")
    core_members <- members[core_ids]
    extension_members <- members[nmr_ids]
    if (is.null(core_fit)) 
        core_fit <- engine$c39_fit_a(data, core_members, base, weights)
    nmrA_assert(identical(names(core_fit$models), core_ids) && identical(core_fit$members, core_members), "NMR_SELECTED_BASE11_PREFIX_FIT_BINDING_MISMATCH")
    if (is.function(core_callback)) 
        core_callback(core_fit)
    original <- data
    raw <- nmrA_raw_values(data, extension_members)
    previous <- vapply(core_members, `[[`, "", "column")
    steps <- list()
    coordinates <- character()
    manifest <- list()
    offset <- 11L
    for (block in blockspec$blocks) {
        nmrA_assert(is.character(block$id) && length(block$id) == 1L && grepl("^[A-Za-z][A-Za-z0-9_]*$", block$id), 
            "NMR_STABLE_BLOCK_ID_REQUIRED")
        prefix <- paste0(".nmrA_", block$id, "_")
        if (identical(block$type, "scalar")) {
            field <- paste0("F", block$field)
            key <- paste0(prefix, "linear")
            member <- extension_members[[paste0("nmr__nmr_", block$field)]]
            bounds <- nmrA_member_bounds(block$field, members, blockspec$raw_support)
            f <- nmrA_conditional_fit(data, raw[[field]], base, previous, weights, engine, key, bounds[1L], bounds[2L], 
                "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL")
            offset <- offset + 1L
            data[[key]] <- nmrA_predictor(f, raw[[field]])
            steps[[block$id]] <- list(type = "scalar", block = block, fit = f, key = key, uniform = offset)
            previous <- c(previous, key)
            coordinates <- c(coordinates, key)
            manifest[[block$id]] <- list(type = "scalar", raw_fields = block$field, family = "INHERITED_IDENTITY_LINEAR_WEIGHTED_EMPIRICAL_RESIDUAL", 
                uniform_columns = setNames(offset, key))
            next
        }
        nmrA_assert(identical(block$type, "composition"), "NMR_BLOCK_TYPE_UNREVIEWED")
        rt <- nmrA_roundtrip(block, raw, nrow(data))
        enc <- rt$encoded
        D <- ncol(enc$parts)
        mask_key <- paste0(prefix, "mask")
        log_key <- paste0(prefix, "logtotal")
        ilr_keys <- paste0(prefix, "ilr", seq_len(D - 1L))
        design <- engine$c39_a_design(data, base, previous)
        status <- nmrA_status_fit(design$X, enc$mask, weights, block$id)
        status$design_reference <- design$reference
        fixed_total <- !is.null(block$fixed_total)
        slot_names <- c(mask_key, if (!fixed_total) log_key, ilr_keys)
        slots <- setNames(offset + seq_along(slot_names), slot_names)
        offset <- max(slots)
        data[[mask_key]] <- factor(enc$mask, levels = status$levels)
        data[[log_key]] <- 0
        for (key in ilr_keys) data[[key]] <- 0
        pattern_fits <- list()
        for (pattern in status$levels) {
            ii <- enc$mask == pattern
            active <- which(strsplit(pattern, "", fixed = TRUE)[[1L]] == "1")
            pf <- list(active = active, observed_n = sum(ii), total = NULL, ilr = list())
            if (length(active)) {
                nmrA_assert(!identical(block$total_support, "observed_range"), "NMR_OBSERVED_EXTREMA_ARE_NOT_PHYSICAL_TOTAL_SUPPORT")
                y <- if (fixed_total) 
                  rep(log(block$fixed_total), sum(ii))
                else log(enc$total[ii])
                bounds <- c(-Inf, Inf)
                z <- nmrA_ilr(enc$parts[ii, active, drop = FALSE])
                for (j in seq_len(ncol(z))) {
                  key <- ilr_keys[j]
                  pred <- c(previous, ilr_keys[seq_len(j - 1L)])
                  pf$ilr[[key]] <- nmrA_conditional_fit(data[ii, , drop = FALSE], z[, j], base, pred, weights[ii], 
                    engine, key)
                  data[[key]][ii] <- nmrA_predictor(pf$ilr[[key]], z[, j])
                }
                if (!fixed_total) {
                  pf$total <- nmrA_conditional_fit(data[ii, , drop = FALSE], y, base, c(previous, names(pf$ilr)), 
                    weights[ii], engine, log_key)
                  data[[log_key]][ii] <- nmrA_predictor(pf$total, y)
                }
            }
            pattern_fits[[pattern]] <- pf
        }
        current <- c(if (!fixed_total) log_key, ilr_keys)
        if (length(status$levels) > 1L) 
            current <- c(mask_key, current)
        auxiliary <- list()
        for (key in names(enc$aux)) {
            column <- paste0(prefix, key)
            y <- enc$aux[[key]]
            bounds <- nmrA_aux_bounds(block$aux[[key]], y)
            f <- nmrA_conditional_fit(data, y, base, c(previous, current), weights, engine, column, bounds[1L], 
                bounds[2L], if (all(is.infinite(bounds))) 
                  "EMPIRICAL_RESIDUAL_UNTRUNCATED"
                else "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL")
            offset <- offset + 1L
            slots[column] <- offset
            data[[column]] <- nmrA_predictor(f, y)
            auxiliary[[key]] <- list(fit = f, column = column, uniform = offset, support_policy = if (is.character(block$aux[[key]])) "unbounded" else block$aux[[key]]$support)
            current <- c(current, column)
        }
        steps[[block$id]] <- list(type = "composition", block = block, status = status, patterns = pattern_fits, 
            mask_key = mask_key, log_key = log_key, ilr_keys = ilr_keys, auxiliary = auxiliary, slots = slots)
        previous <- c(previous, current)
        coordinates <- c(coordinates, names(slots))
        manifest[[block$id]] <- list(type = "real_biochemical_composition", raw_fields = names(block$raw), parts = names(block$parts), 
            observed_zero_patterns = status$levels, observed_pattern_n = as.list(table(enc$mask)), zero_policy = "EXPLICIT_CONDITIONAL_MASK_NO_EPSILON", 
            independent_positive_composition_dimensions = D - 1L, variable_total_conditioned_on_ilr = !fixed_total, 
            fixed_total = block$fixed_total, total_support = if (is.null(block$total_support)) "positive_unbounded" else block$total_support, 
            measurement_auxiliary = block$aux, observed_raw_roundtrip = rt[-1L], uniform_columns = as.list(slots), 
            fixed_mask_uniform_reserved = length(status$levels) == 1L)
    }
    nmrA_assert(identical(original, data[names(original)]), "NMR_A_MODIFIED_OBSERVED_RAW_OR_BASE_COLUMNS")
    list(revision = nmrA_revision(), base11 = core_fit, extension = steps, members = members, core_ids = core_ids, 
        nmr_ids = nmr_ids, base = base, uniform_dim = as.integer(offset), uniform_columns = c(setNames(seq_len(11L), 
            core_ids), setNames(seq.int(12L, offset), coordinates)), manifest = list(revision = nmrA_revision(), 
            raw_information_count = length(nmr_ids), common_SELECTED_BASE11_uniform_prefix = 11L, uniform_dim = as.integer(offset), 
            positive_interior_extension_coordinates = as.integer(offset - 11L - sum(vapply(steps, function(s) identical(s$type, 
                "composition"), logical(1)))), reserved_zero_mask_uniform_slots = sum(vapply(steps, function(s) identical(s$type, 
                "composition"), logical(1))), coordinate_basis = "ORTHONORMAL_HELMERT_ILR_OF_POSITIVE_SUBCOMPOSITION", 
            conditional_predecessors = "SELECTED_BASE11_RAW_THEN_INDEPENDENT_SELECTED29_EXTENSION_COORDINATES", 
            original_raw_Y_columns_retained = TRUE, no_PCA_no_ridge_no_correlation_selection = TRUE, predictor_conditioning = "DECLARED_REVERSIBLE_WEIGHTED_CENTER_AND_SCALE_EXTENSION_ONLY", 
            optimizer_RNG_consumed = !identical(rng, nmrA_rng()), optimizer_RNG = list(before = rng, after = nmrA_rng()), 
            measurement_relation_errors_not_silently_discarded = TRUE, blocks = manifest), raw_support = blockspec$raw_support, 
        engine = engine)
}

draw_nmr_composition_A <- function(a, C, b, uniforms, audit = NULL, engine = a$engine, required_ids = NULL) {
    rng <- nmrA_rng()
    on.exit(nmrA_assert(identical(rng, nmrA_rng()), "NMR_A_DRAW_CHANGED_RNG"), add = TRUE)
    nmrA_assert(identical(a$revision, nmrA_revision()) && is.matrix(uniforms) && nrow(uniforms) == nrow(C) && ncol(uniforms) == 
        a$uniform_dim && all(is.finite(uniforms) & uniforms >= 0 & uniforms < 1), "NMR_A_DRAW_UNIFORM_CONTRACT_INVALID")
    nmrA_assert(!any(grepl("^\\.nmrA_", names(C))), "NMR_A_RESERVED_COORDINATE_COLUMN_EXISTS")
    if (!is.null(required_ids)) {
        nmrA_assert(is.character(required_ids) && length(required_ids) > 0L && !anyNA(required_ids) && !anyDuplicated(required_ids) && 
            all(required_ids %in% names(a$members)), "NMR_A_REQUESTED_DOMAIN_IDS_INVALID")
        if (all(required_ids %in% a$core_ids)) {
            last <- max(match(required_ids, a$core_ids))
            ids <- a$core_ids[seq_len(last)]
            prefix <- a$base11
            prefix$models <- prefix$models[ids]
            prefix$members <- prefix$members[ids]
            prefix$conditional_order <- ids
            return(engine$c39_draw_a(prefix, C, b, uniforms[, seq_len(last), drop = FALSE], audit))
        }
    }
    d <- engine$c39_draw_a(a$base11, C, b, uniforms[, seq_len(11L), drop = FALSE], audit)
    raw <- list()
    for (step in a$extension) {
        block <- step$block
        if (identical(step$type, "scalar")) {
            value <- nmrA_conditional_draw(step$fit, d, uniforms[, step$uniform], engine, audit)
            d[[step$key]] <- nmrA_predictor(step$fit, value)
            raw[[paste0("F", block$field)]] <- value
        }
        else {
            X <- engine$c39_a_design(d, NULL, NULL, step$status$design_reference)$X
            mask <- nmrA_status_draw(step$status, X, uniforms[, step$slots[[step$mask_key]]])
            d[[step$mask_key]] <- factor(mask, levels = step$status$levels)
            d[[step$log_key]] <- 0
            for (key in step$ilr_keys) d[[key]] <- 0
            parts <- matrix(0, nrow(d), length(block$parts), dimnames = list(NULL, names(block$parts)))
            for (pattern in unique(mask)) {
                ii <- mask == pattern
                pf <- step$patterns[[pattern]]
                nmrA_assert(is.list(pf), "NMR_UNOBSERVED_ZERO_MASK_DRAWN")
                if (!length(pf$active)) 
                  next
                z <- matrix(numeric(), sum(ii), 0L)
                for (key in names(pf$ilr)) {
                  value <- nmrA_conditional_draw(pf$ilr[[key]], d[ii, , drop = FALSE], uniforms[ii, step$slots[[key]]], 
                    engine, audit)
                  d[[key]][ii] <- nmrA_predictor(pf$ilr[[key]], value)
                  z <- cbind(z, value)
                }
                proportions <- matrix(0, sum(ii), ncol(parts), dimnames = list(NULL, colnames(parts)))
                proportions[, pf$active] <- nmrA_inverse_ilr(z, 1, length(pf$active))
                bounds <- if (is.null(block$fixed_total)) 
                  nmrA_total_bounds(block, proportions, lapply(raw, `[`, ii), a$members, a$raw_support)
                else list(lower = rep(-Inf, sum(ii)), upper = rep(Inf, sum(ii)))
                logtotal <- if (!is.null(block$fixed_total)) 
                  rep(log(block$fixed_total), sum(ii))
                else nmrA_conditional_draw(pf$total, d[ii, , drop = FALSE], uniforms[ii, step$slots[[step$log_key]]], 
                  engine, audit, bounds$lower, bounds$upper)
                if (!is.null(pf$total)) 
                  d[[step$log_key]][ii] <- nmrA_predictor(pf$total, logtotal)
                parts[ii, ] <- proportions * exp(logtotal)
            }
            aux <- list()
            for (key in names(step$auxiliary)) {
                af <- step$auxiliary[[key]]
                bounds <- nmrA_aux_dynamic_bounds(block, key, parts, aux, raw, a$members, a$raw_support)
                value <- nmrA_conditional_draw(af$fit, d, uniforms[, af$uniform], engine, audit, pmax(af$fit$lower, 
                  bounds$lower), pmin(af$fit$upper, bounds$upper))
                aux[[key]] <- value
                d[[af$column]] <- nmrA_predictor(af$fit, value)
            }
            decoded <- nmrA_block_decode(block, parts, aux, raw)
            raw[names(decoded)] <- decoded
        }
        for (field in nmrA_block_fields(block)) {
            member <- a$members[[paste0("nmr__nmr_", field)]]
            d[[member$column]] <- raw[[paste0("F", field)]]
        }
    }
    engine$c39_support(d, a$members)
    for (id in a$nmr_ids) {
        field <- sub("^nmr__nmr_", "", id)
        bounds <- nmrA_member_bounds(field, a$members, a$raw_support)
        value <- d[[a$members[[id]]$column]]
        nmrA_assert(all(is.finite(value) & value >= bounds[1L] & value <= bounds[2L]), paste("NMR_RECONSTRUCTED_REPORTED_RAW_PHYSICAL_SUPPORT_FAILURE", 
            field))
    }
    nmrA_assert(all(vapply(d[vapply(a$members[a$nmr_ids], `[[`, "", "column")], function(x) all(is.finite(x)), 
        logical(1))), "NMR_RAW_RECONSTRUCTION_NONFINITE")
    d[, setdiff(names(d), grep("^\\.nmrA_", names(d), value = TRUE)), drop = FALSE]
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        nmrA_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
        nmrA_revision <- function() "SELECTED44_BASE11_RAW29_COMPOSITION_RECONSTRUCTION_V1"
        nmrA_rng <- function() list(kind = RNGkind(), present = exists(".Random.seed", .GlobalEnv, inherits = FALSE), 
            seed = if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL)
        nmrA_eval <- function(expression, values, n) {
            walk <- function(x) {
                if (is.numeric(x)) return(x)
                if (is.symbol(x)) {
                  key <- as.character(x)
                  nmrA_assert(key %in% names(values), paste("NMR_EXPRESSION_UNKNOWN_SYMBOL", key))
                  return(values[[key]])
                }
                nmrA_assert(is.call(x), "NMR_EXPRESSION_LANGUAGE_INVALID")
                op <- as.character(x[[1L]])
                nmrA_assert(op %in% c("+", "-", "*", "/", "log", "exp", "("), paste("NMR_EXPRESSION_OPERATION_NOT_ALLOWED", 
                  op))
                args <- lapply(as.list(x)[-1L], walk)
                nmrA_assert(length(args) %in% if (op %in% c("+", "-")) 1:2 else if (op %in% c("*", "/")) 2L else 1L, 
                  "NMR_EXPRESSION_ARITY_INVALID")
                do.call(switch(op, `+` = `+`, `-` = `-`, `*` = `*`, `/` = `/`, log = log, exp = exp, `(` = identity), 
                  args)
            }
            ast <- parse(text = expression, keep.source = FALSE)
            nmrA_assert(length(ast) == 1L, "NMR_EXPRESSION_SINGLE_EXPRESSION_REQUIRED")
            answer <- walk(ast[[1L]])
            nmrA_assert(is.numeric(answer) && length(answer) %in% c(1L, n) && all(is.finite(answer)), paste("NMR_EXPRESSION_NONFINITE", 
                expression))
            rep_len(answer, n)
        }
        nmrA_basis <- function(D) {
            nmrA_assert(length(D) == 1L && D >= 1L, "NMR_ILR_DIMENSION_INVALID")
            if (D == 1L) return(matrix(numeric(), 1L, 0L))
            H <- stats::contr.helmert(D)
            sweep(H, 2L, sqrt(colSums(H^2)), `/`)
        }
        nmrA_ilr <- function(parts) {
            nmrA_assert(is.matrix(parts) && is.numeric(parts) && ncol(parts) >= 1L && all(is.finite(parts) & parts > 
                0), "NMR_ILR_STRICTLY_POSITIVE_PARTS_REQUIRED_NO_EPSILON")
            log(parts) %*% nmrA_basis(ncol(parts))
        }
        nmrA_inverse_ilr <- function(z, total, D = ncol(z) + 1L) {
            nmrA_assert(is.matrix(z) && ncol(z) == D - 1L && all(is.finite(z)) && length(total) %in% c(1L, nrow(z)) && 
                all(is.finite(total) & total > 0), "NMR_INVERSE_ILR_INVALID_INPUT")
            lp <- z %*% t(nmrA_basis(D))
            lp <- lp - apply(lp, 1L, max)
            p <- exp(lp)
            answer <- p/rowSums(p) * rep_len(total, nrow(z))
            nmrA_assert(all(is.finite(answer) & answer > 0), "NMR_INVERSE_ILR_UNDERFLOW_OR_OVERFLOW_NO_CLIPPING")
            answer
        }
        nmrA_raw_values <- function(data, members) {
            fields <- sub("^nmr__nmr_", "F", names(members))
            cols <- vapply(members, `[[`, "", "column")
            nmrA_assert(all(cols %in% names(data)) && !anyDuplicated(fields) && all(vapply(data[cols], is.numeric, 
                logical(1))), "NMR_RAW_NUMERIC_COLUMNS_REQUIRED")
            setNames(as.list(data[cols]), fields)
        }
        nmrA_block_fields <- function(block) if (identical(block$type, "scalar")) as.character(block$field) else names(block$raw)
        nmrA_block_encode <- function(block, raw, n) {
            nmrA_assert(identical(block$type, "composition"), "NMR_COMPOSITION_BLOCK_REQUIRED")
            part_values <- list()
            values <- raw
            for (key in names(block$parts)) {
                part_values[[key]] <- nmrA_eval(block$parts[[key]], values, n)
                values[[key]] <- part_values[[key]]
            }
            parts <- do.call(cbind, part_values)
            colnames(parts) <- names(block$parts)
            nmrA_assert(ncol(parts) >= 2L && !anyDuplicated(colnames(parts)) && all(is.finite(parts) & parts >= 
                0), paste("NMR_NEGATIVE_BIOCHEMICAL_PART_NO_REPAIR", block$id))
            total <- rowSums(parts)
            if (!is.null(block$fixed_total)) nmrA_assert(all(abs(total - block$fixed_total) <= 512 * .Machine$double.eps * 
                pmax(1, abs(total))), paste("NMR_SOURCE_FIXED_TOTAL_IDENTITY_FAILED", block$id))
            masks <- apply(parts > 0, 1L, function(x) paste0(as.integer(x), collapse = ""))
            aux <- list()
            values <- c(raw, as.list(as.data.frame(parts)), list(TOTAL = total))
            for (key in names(block$aux)) {
                specification <- block$aux[[key]]
                expression <- if (is.character(specification)) specification else specification$expression
                aux[[key]] <- nmrA_eval(expression, values, n)
                values[[key]] <- aux[[key]]
            }
            list(parts = parts, total = total, mask = masks, aux = aux)
        }
        nmrA_block_decode <- function(block, parts, aux, prior = list()) {
            n <- nrow(parts)
            values <- c(as.list(as.data.frame(parts)), list(TOTAL = rowSums(parts)), aux, prior)
            setNames(lapply(block$raw, nmrA_eval, values = values, n = n), paste0("F", names(block$raw)))
        }
        nmrA_roundtrip <- function(block, raw, n) {
            encoded <- nmrA_block_encode(block, raw, n)
            decoded <- nmrA_block_decode(block, encoded$parts, encoded$aux, raw)
            discrepancies <- vapply(names(decoded), function(key) max(abs(decoded[[key]] - raw[[key]])/pmax(1, 
                abs(raw[[key]]))), numeric(1))
            tolerance <- 1e-10
            nmrA_assert(all(discrepancies <= tolerance), paste("NMR_RAW_ROUNDTRIP_LOSES_RECORDED_VARIATION", block$id, 
                paste(names(discrepancies)[discrepancies > tolerance], collapse = ",")))
            list(encoded = encoded, max_scaled_error = max(discrepancies), arithmetic_scaled_tolerance = tolerance, 
                exact_recorded_derived_identities_claimed = FALSE)
        }
        nmrA_status_fit <- function(X, status, weights, id) {
            levels <- sort(unique(status))
            y <- factor(status, levels = levels)
            if (length(levels) == 1L) return(list(kind = "constant", levels = levels))
            nmrA_assert(requireNamespace("nnet", quietly = TRUE) && qr(X)$rank == ncol(X), paste("NMR_ZERO_MASK_DESIGN_RANK_DEFICIENT", 
                id))
            warnings <- character()
            fit <- withCallingHandlers(nnet::multinom(y ~ X - 1, weights = weights, trace = FALSE, maxit = 1000L, 
                MaxNWts = 100000L, model = FALSE), warning = function(e) {
                warnings <<- c(warnings, conditionMessage(e))
                invokeRestart("muffleWarning")
            })
            coef <- stats::coef(fit)
            nmrA_assert(identical(fit$convergence, 0L) && all(is.finite(coef)), paste("NMR_ZERO_MASK_MODEL_NOT_CONVERGED", 
                id))
            list(kind = "multinomial", levels = levels, coefficients = coef, warnings = unique(warnings))
        }
        nmrA_status_draw <- function(fit, X, u) {
            if (identical(fit$kind, "constant")) return(rep(fit$levels, nrow(X)))
            if (length(fit$levels) == 2L) {
                pp <- stats::plogis(as.numeric(X %*% fit$coefficients))
                p <- cbind(1 - pp, pp)
            } else {
                lp <- cbind(0, X %*% t(fit$coefficients))
                lp <- lp - apply(lp, 1L, max)
                p <- exp(lp)
                p <- p/rowSums(p)
            }
            nmrA_assert(all(is.finite(p) & p >= 0 & p <= 1), "NMR_ZERO_MASK_PROBABILITY_INVALID")
            cumulative <- p
            if (ncol(p) > 1L) for (j in 2:ncol(p)) cumulative[, j] <- cumulative[, j - 1L] + p[, j]
            fit$levels[1L + rowSums(u >= cumulative[, -ncol(p), drop = FALSE])]
        }
        nmrA_conditional_fit <- function(data, y, base, previous, weights, engine, id, lower = -Inf, upper = Inf, 
            rule = "EMPIRICAL_RESIDUAL_UNTRUNCATED") {
            design <- engine$c39_a_design(data, base, previous)
            fit <- engine$c39_linear_fit(design$X, y, weights, id, rule)
            fit$design_reference <- design$reference
            fit$id <- id
            fit$lower <- lower
            fit$upper <- upper
            fit$predictor_center <- sum(weights * y)/sum(weights)
            fit$predictor_scale <- sqrt(sum(weights * (y - fit$predictor_center)^2)/sum(weights))
            if (fit$predictor_scale == 0) fit$predictor_scale <- 1
            fit
        }
        nmrA_predictor <- function(fit, value) (value - fit$predictor_center)/fit$predictor_scale
        nmrA_predictor_inverse <- function(fit, value) value * fit$predictor_scale + fit$predictor_center
        nmrA_conditional_draw <- function(fit, data, u, engine, audit = NULL, lower = fit$lower, upper = fit$upper) {
            X <- engine$c39_a_design(data, NULL, NULL, fit$design_reference)$X
            f <- fit
            if (any(is.finite(lower)) || any(is.finite(upper))) f$rule <- "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL"
            engine$c39_residual_draw(f, X, u, lower, upper, FALSE, fit$id, audit)
        }
        nmrA_member_bounds <- function(field, members, raw_support = NULL) {
            m <- members[[paste0("nmr__nmr_", field)]]
            support <- raw_support[[field]]
            c(max(m$lower, if (is.null(support$lower)) -Inf else support$lower), min(m$upper, if (is.null(support$upper)) Inf else support$upper))
        }
        nmrA_affine_bounds <- function(block, values0, values1, members, raw_support, variable, allowed_aux = character()) {
            n <- length(values0$TOTAL)
            lower <- rep(-Inf, n)
            upper <- rep(Inf, n)
            for (field in names(block$raw)) {
                expression <- block$raw[[field]]
                symbols <- all.names(parse(text = expression)[[1L]], functions = FALSE, unique = TRUE)
                unresolved <- setdiff(intersect(symbols, names(block$aux)), allowed_aux)
                if (length(unresolved)) next
                a <- nmrA_eval(expression, values0, n)
                b <- nmrA_eval(expression, values1, n) - a
                bounds <- nmrA_member_bounds(field, members, raw_support)
                positive <- b > 0
                negative <- b < 0
                fixed <- b == 0
                lower[positive] <- pmax(lower[positive], (bounds[1L] - a[positive])/b[positive])
                upper[positive] <- pmin(upper[positive], (bounds[2L] - a[positive])/b[positive])
                lower[negative] <- pmax(lower[negative], (bounds[2L] - a[negative])/b[negative])
                upper[negative] <- pmin(upper[negative], (bounds[1L] - a[negative])/b[negative])
                nmrA_assert(all(!fixed | (a >= bounds[1L] & a <= bounds[2L])), paste("NMR_RECONSTRUCTION_FIXED_SUPPORT_FAILURE", 
                  block$id, field, variable))
            }
            nmrA_assert(all(lower <= upper), paste("NMR_DYNAMIC_PHYSICAL_SUPPORT_EMPTY", block$id, variable))
            list(lower = lower, upper = upper)
        }
        nmrA_total_bounds <- function(block, proportions, raw, members, raw_support) {
            n <- nrow(proportions)
            if (!identical(block$component_unit, "reported_percent")) return(list(lower = rep(-Inf, n), upper = rep(Inf, 
                n)))
            p0 <- proportions * 0
            values0 <- c(as.list(as.data.frame(p0)), list(TOTAL = rep(0, n)), raw)
            values1 <- c(as.list(as.data.frame(proportions)), list(TOTAL = rep(1, n)), raw)
            bounds <- nmrA_affine_bounds(block, values0, values1, members, raw_support, "TOTAL")
            nmrA_assert(all(bounds$upper > 0), "NMR_POSITIVE_TOTAL_PHYSICAL_SUPPORT_EMPTY")
            list(lower = ifelse(bounds$lower > 0, log(bounds$lower), -Inf), upper = log(bounds$upper))
        }
        nmrA_aux_dynamic_bounds <- function(block, key, parts, auxiliary, raw, members, raw_support) {
            n <- nrow(parts)
            common <- c(as.list(as.data.frame(parts)), list(TOTAL = rowSums(parts)), raw, auxiliary)
            v0 <- common
            v1 <- common
            v0[[key]] <- rep(0, n)
            v1[[key]] <- rep(1, n)
            nmrA_affine_bounds(block, v0, v1, members, raw_support, key, c(names(auxiliary), key))
        }
        nmrA_aux_bounds <- function(specification, y) {
            if (is.character(specification)) return(c(-Inf, Inf))
            nmrA_assert(!identical(specification$support, "observed_range"), "NMR_OBSERVED_EXTREMA_ARE_NOT_PHYSICAL_MEASUREMENT_SUPPORT")
            c(if (is.null(specification$lower)) -Inf else as.numeric(specification$lower), if (is.null(specification$upper)) Inf else as.numeric(specification$upper))
        }
        fit_nmr_composition_A <- function(data, members, base, weights, engine, blockspec, core_fit = NULL, core_callback = NULL) {
            rng <- nmrA_rng()
            nmrA_assert(is.data.frame(data) && is.list(members) && is.environment(engine) && length(weights) == 
                nrow(data) && all(is.finite(weights) & weights > 0) && is.finite(sum(weights)), "NMR_A_INVALID_INPUT")
            nmrA_assert(!any(grepl("^\\.nmrA_", names(data))), "NMR_A_RESERVED_COORDINATE_COLUMN_EXISTS")
            nmrA_assert(identical(blockspec$schema, "NMR_COMPOSITION_BLOCKSPEC_V1") && length(blockspec$blocks) > 
                0L, "NMR_REVIEWED_BLOCKSPEC_REQUIRED")
            core_ids <- names(members)[seq_len(11L)]
            nmr_ids <- names(members)[-(seq_len(11L))]
            expected <- if (isTRUE(blockspec$synthetic)) length(nmr_ids) else 29L
            nmrA_assert(length(core_ids) == 11L && length(nmr_ids) == expected && all(grepl("^nmr__nmr_[0-9]+$", 
                nmr_ids)), "NMR_A_SELECTED_BASE11_PLUS_FROZEN29_REQUIRED")
            covered <- unlist(lapply(blockspec$blocks, nmrA_block_fields), use.names = FALSE)
            nmrA_assert(!anyDuplicated(covered) && setequal(covered, sub("^nmr__nmr_", "", nmr_ids)), "NMR_BLOCKS_MUST_PARTITION_ALL_RAW_INFORMATION_NO_SELECTION")
            core_members <- members[core_ids]
            extension_members <- members[nmr_ids]
            if (is.null(core_fit)) core_fit <- engine$c39_fit_a(data, core_members, base, weights)
            nmrA_assert(identical(names(core_fit$models), core_ids) && identical(core_fit$members, core_members), 
                "NMR_SELECTED_BASE11_PREFIX_FIT_BINDING_MISMATCH")
            if (is.function(core_callback)) core_callback(core_fit)
            original <- data
            raw <- nmrA_raw_values(data, extension_members)
            previous <- vapply(core_members, `[[`, "", "column")
            steps <- list()
            coordinates <- character()
            manifest <- list()
            offset <- 11L
            for (block in blockspec$blocks) {
                nmrA_assert(is.character(block$id) && length(block$id) == 1L && grepl("^[A-Za-z][A-Za-z0-9_]*$", 
                  block$id), "NMR_STABLE_BLOCK_ID_REQUIRED")
                prefix <- paste0(".nmrA_", block$id, "_")
                if (identical(block$type, "scalar")) {
                  field <- paste0("F", block$field)
                  key <- paste0(prefix, "linear")
                  member <- extension_members[[paste0("nmr__nmr_", block$field)]]
                  bounds <- nmrA_member_bounds(block$field, members, blockspec$raw_support)
                  f <- nmrA_conditional_fit(data, raw[[field]], base, previous, weights, engine, key, bounds[1L], 
                    bounds[2L], "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL")
                  offset <- offset + 1L
                  data[[key]] <- nmrA_predictor(f, raw[[field]])
                  steps[[block$id]] <- list(type = "scalar", block = block, fit = f, key = key, uniform = offset)
                  previous <- c(previous, key)
                  coordinates <- c(coordinates, key)
                  manifest[[block$id]] <- list(type = "scalar", raw_fields = block$field, family = "INHERITED_IDENTITY_LINEAR_WEIGHTED_EMPIRICAL_RESIDUAL", 
                    uniform_columns = setNames(offset, key))
                  next
                }
                nmrA_assert(identical(block$type, "composition"), "NMR_BLOCK_TYPE_UNREVIEWED")
                rt <- nmrA_roundtrip(block, raw, nrow(data))
                enc <- rt$encoded
                D <- ncol(enc$parts)
                mask_key <- paste0(prefix, "mask")
                log_key <- paste0(prefix, "logtotal")
                ilr_keys <- paste0(prefix, "ilr", seq_len(D - 1L))
                design <- engine$c39_a_design(data, base, previous)
                status <- nmrA_status_fit(design$X, enc$mask, weights, block$id)
                status$design_reference <- design$reference
                fixed_total <- !is.null(block$fixed_total)
                slot_names <- c(mask_key, if (!fixed_total) log_key, ilr_keys)
                slots <- setNames(offset + seq_along(slot_names), slot_names)
                offset <- max(slots)
                data[[mask_key]] <- factor(enc$mask, levels = status$levels)
                data[[log_key]] <- 0
                for (key in ilr_keys) data[[key]] <- 0
                pattern_fits <- list()
                for (pattern in status$levels) {
                  ii <- enc$mask == pattern
                  active <- which(strsplit(pattern, "", fixed = TRUE)[[1L]] == "1")
                  pf <- list(active = active, observed_n = sum(ii), total = NULL, ilr = list())
                  if (length(active)) {
                    nmrA_assert(!identical(block$total_support, "observed_range"), "NMR_OBSERVED_EXTREMA_ARE_NOT_PHYSICAL_TOTAL_SUPPORT")
                    y <- if (fixed_total) rep(log(block$fixed_total), sum(ii)) else log(enc$total[ii])
                    bounds <- c(-Inf, Inf)
                    z <- nmrA_ilr(enc$parts[ii, active, drop = FALSE])
                    for (j in seq_len(ncol(z))) {
                      key <- ilr_keys[j]
                      pred <- c(previous, ilr_keys[seq_len(j - 1L)])
                      pf$ilr[[key]] <- nmrA_conditional_fit(data[ii, , drop = FALSE], z[, j], base, pred, weights[ii], 
                        engine, key)
                      data[[key]][ii] <- nmrA_predictor(pf$ilr[[key]], z[, j])
                    }
                    if (!fixed_total) {
                      pf$total <- nmrA_conditional_fit(data[ii, , drop = FALSE], y, base, c(previous, names(pf$ilr)), 
                        weights[ii], engine, log_key)
                      data[[log_key]][ii] <- nmrA_predictor(pf$total, y)
                    }
                  }
                  pattern_fits[[pattern]] <- pf
                }
                current <- c(if (!fixed_total) log_key, ilr_keys)
                if (length(status$levels) > 1L) current <- c(mask_key, current)
                auxiliary <- list()
                for (key in names(enc$aux)) {
                  column <- paste0(prefix, key)
                  y <- enc$aux[[key]]
                  bounds <- nmrA_aux_bounds(block$aux[[key]], y)
                  f <- nmrA_conditional_fit(data, y, base, c(previous, current), weights, engine, column, bounds[1L], 
                    bounds[2L], if (all(is.infinite(bounds))) "EMPIRICAL_RESIDUAL_UNTRUNCATED" else "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL")
                  offset <- offset + 1L
                  slots[column] <- offset
                  data[[column]] <- nmrA_predictor(f, y)
                  auxiliary[[key]] <- list(fit = f, column = column, uniform = offset, support_policy = if (is.character(block$aux[[key]])) "unbounded" else block$aux[[key]]$support)
                  current <- c(current, column)
                }
                steps[[block$id]] <- list(type = "composition", block = block, status = status, patterns = pattern_fits, 
                  mask_key = mask_key, log_key = log_key, ilr_keys = ilr_keys, auxiliary = auxiliary, slots = slots)
                previous <- c(previous, current)
                coordinates <- c(coordinates, names(slots))
                manifest[[block$id]] <- list(type = "real_biochemical_composition", raw_fields = names(block$raw), 
                  parts = names(block$parts), observed_zero_patterns = status$levels, observed_pattern_n = as.list(table(enc$mask)), 
                  zero_policy = "EXPLICIT_CONDITIONAL_MASK_NO_EPSILON", independent_positive_composition_dimensions = D - 
                    1L, variable_total_conditioned_on_ilr = !fixed_total, fixed_total = block$fixed_total, total_support = if (is.null(block$total_support)) "positive_unbounded" else block$total_support, 
                  measurement_auxiliary = block$aux, observed_raw_roundtrip = rt[-1L], uniform_columns = as.list(slots), 
                  fixed_mask_uniform_reserved = length(status$levels) == 1L)
            }
            nmrA_assert(identical(original, data[names(original)]), "NMR_A_MODIFIED_OBSERVED_RAW_OR_BASE_COLUMNS")
            list(revision = nmrA_revision(), base11 = core_fit, extension = steps, members = members, core_ids = core_ids, 
                nmr_ids = nmr_ids, base = base, uniform_dim = as.integer(offset), uniform_columns = c(setNames(seq_len(11L), 
                  core_ids), setNames(seq.int(12L, offset), coordinates)), manifest = list(revision = nmrA_revision(), 
                  raw_information_count = length(nmr_ids), common_SELECTED_BASE11_uniform_prefix = 11L, uniform_dim = as.integer(offset), 
                  positive_interior_extension_coordinates = as.integer(offset - 11L - sum(vapply(steps, function(s) identical(s$type, 
                    "composition"), logical(1)))), reserved_zero_mask_uniform_slots = sum(vapply(steps, function(s) identical(s$type, 
                    "composition"), logical(1))), coordinate_basis = "ORTHONORMAL_HELMERT_ILR_OF_POSITIVE_SUBCOMPOSITION", 
                  conditional_predecessors = "SELECTED_BASE11_RAW_THEN_INDEPENDENT_SELECTED29_EXTENSION_COORDINATES", 
                  original_raw_Y_columns_retained = TRUE, no_PCA_no_ridge_no_correlation_selection = TRUE, predictor_conditioning = "DECLARED_REVERSIBLE_WEIGHTED_CENTER_AND_SCALE_EXTENSION_ONLY", 
                  optimizer_RNG_consumed = !identical(rng, nmrA_rng()), optimizer_RNG = list(before = rng, after = nmrA_rng()), 
                  measurement_relation_errors_not_silently_discarded = TRUE, blocks = manifest), raw_support = blockspec$raw_support, 
                engine = engine)
        }
        draw_nmr_composition_A <- function(a, C, b, uniforms, audit = NULL, engine = a$engine, required_ids = NULL) {
            rng <- nmrA_rng()
            on.exit(nmrA_assert(identical(rng, nmrA_rng()), "NMR_A_DRAW_CHANGED_RNG"), add = TRUE)
            nmrA_assert(identical(a$revision, nmrA_revision()) && is.matrix(uniforms) && nrow(uniforms) == nrow(C) && 
                ncol(uniforms) == a$uniform_dim && all(is.finite(uniforms) & uniforms >= 0 & uniforms < 1), "NMR_A_DRAW_UNIFORM_CONTRACT_INVALID")
            nmrA_assert(!any(grepl("^\\.nmrA_", names(C))), "NMR_A_RESERVED_COORDINATE_COLUMN_EXISTS")
            if (!is.null(required_ids)) {
                nmrA_assert(is.character(required_ids) && length(required_ids) > 0L && !anyNA(required_ids) && 
                  !anyDuplicated(required_ids) && all(required_ids %in% names(a$members)), "NMR_A_REQUESTED_DOMAIN_IDS_INVALID")
                if (all(required_ids %in% a$core_ids)) {
                  last <- max(match(required_ids, a$core_ids))
                  ids <- a$core_ids[seq_len(last)]
                  prefix <- a$base11
                  prefix$models <- prefix$models[ids]
                  prefix$members <- prefix$members[ids]
                  prefix$conditional_order <- ids
                  return(engine$c39_draw_a(prefix, C, b, uniforms[, seq_len(last), drop = FALSE], audit))
                }
            }
            d <- engine$c39_draw_a(a$base11, C, b, uniforms[, seq_len(11L), drop = FALSE], audit)
            raw <- list()
            for (step in a$extension) {
                block <- step$block
                if (identical(step$type, "scalar")) {
                  value <- nmrA_conditional_draw(step$fit, d, uniforms[, step$uniform], engine, audit)
                  d[[step$key]] <- nmrA_predictor(step$fit, value)
                  raw[[paste0("F", block$field)]] <- value
                } else {
                  X <- engine$c39_a_design(d, NULL, NULL, step$status$design_reference)$X
                  mask <- nmrA_status_draw(step$status, X, uniforms[, step$slots[[step$mask_key]]])
                  d[[step$mask_key]] <- factor(mask, levels = step$status$levels)
                  d[[step$log_key]] <- 0
                  for (key in step$ilr_keys) d[[key]] <- 0
                  parts <- matrix(0, nrow(d), length(block$parts), dimnames = list(NULL, names(block$parts)))
                  for (pattern in unique(mask)) {
                    ii <- mask == pattern
                    pf <- step$patterns[[pattern]]
                    nmrA_assert(is.list(pf), "NMR_UNOBSERVED_ZERO_MASK_DRAWN")
                    if (!length(pf$active)) next
                    z <- matrix(numeric(), sum(ii), 0L)
                    for (key in names(pf$ilr)) {
                      value <- nmrA_conditional_draw(pf$ilr[[key]], d[ii, , drop = FALSE], uniforms[ii, step$slots[[key]]], 
                        engine, audit)
                      d[[key]][ii] <- nmrA_predictor(pf$ilr[[key]], value)
                      z <- cbind(z, value)
                    }
                    proportions <- matrix(0, sum(ii), ncol(parts), dimnames = list(NULL, colnames(parts)))
                    proportions[, pf$active] <- nmrA_inverse_ilr(z, 1, length(pf$active))
                    bounds <- if (is.null(block$fixed_total)) nmrA_total_bounds(block, proportions, lapply(raw, 
                      `[`, ii), a$members, a$raw_support) else list(lower = rep(-Inf, sum(ii)), upper = rep(Inf, 
                      sum(ii)))
                    logtotal <- if (!is.null(block$fixed_total)) rep(log(block$fixed_total), sum(ii)) else nmrA_conditional_draw(pf$total, 
                      d[ii, , drop = FALSE], uniforms[ii, step$slots[[step$log_key]]], engine, audit, bounds$lower, 
                      bounds$upper)
                    if (!is.null(pf$total)) d[[step$log_key]][ii] <- nmrA_predictor(pf$total, logtotal)
                    parts[ii, ] <- proportions * exp(logtotal)
                  }
                  aux <- list()
                  for (key in names(step$auxiliary)) {
                    af <- step$auxiliary[[key]]
                    bounds <- nmrA_aux_dynamic_bounds(block, key, parts, aux, raw, a$members, a$raw_support)
                    value <- nmrA_conditional_draw(af$fit, d, uniforms[, af$uniform], engine, audit, pmax(af$fit$lower, 
                      bounds$lower), pmin(af$fit$upper, bounds$upper))
                    aux[[key]] <- value
                    d[[af$column]] <- nmrA_predictor(af$fit, value)
                  }
                  decoded <- nmrA_block_decode(block, parts, aux, raw)
                  raw[names(decoded)] <- decoded
                }
                for (field in nmrA_block_fields(block)) {
                  member <- a$members[[paste0("nmr__nmr_", field)]]
                  d[[member$column]] <- raw[[paste0("F", field)]]
                }
            }
            engine$c39_support(d, a$members)
            for (id in a$nmr_ids) {
                field <- sub("^nmr__nmr_", "", id)
                bounds <- nmrA_member_bounds(field, a$members, a$raw_support)
                value <- d[[a$members[[id]]$column]]
                nmrA_assert(all(is.finite(value) & value >= bounds[1L] & value <= bounds[2L]), paste("NMR_RECONSTRUCTED_REPORTED_RAW_PHYSICAL_SUPPORT_FAILURE", 
                  field))
            }
            nmrA_assert(all(vapply(d[vapply(a$members[a$nmr_ids], `[[`, "", "column")], function(x) all(is.finite(x)), 
                logical(1))), "NMR_RAW_RECONSTRUCTION_NONFINITE")
            d[, setdiff(names(d), grep("^\\.nmrA_", names(d), value = TRUE)), drop = FALSE]
        }
    }), envir = context)
}

