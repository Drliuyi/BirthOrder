# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

body_pmm_revision <- function() "SELECTED44_GROUPED_BODY_PREDICTIVE_MEAN_CELLS_V2"

body_pmm_ids <- function() c("weight", "fat_free_mass", "fat_mass", "hip", "waist")

body_pmm_spec <- function() list(method = "WEIGHTED_GROUPED_CONDITIONAL_HOTDECK_MULTIVARIATE_PREDICTIVE_MEAN_CELLS", 
    predictors = "ORIGINAL_B_C_SMOKING_BMI_ONLY", score_responses = c("logT", "logitq", "e", "hip", "waist"), fitting = "EXACT_SELECTION_WEIGHTED_LEAST_SQUARES", 
    score_metric = "WEIGHTED_SD_STANDARDIZED_PREDICTIVE_MEANS", tree = "ONCE_FITTED_BALANCED_KD_PARTITION_MAX_LOCAL_SCORE_RANGE_MEDIAN_DISTINCT_SPLITS", 
    maximum_leaf_rows = 128L, minimum_child_rows = 32L, minimum_unique_eid_per_leaf = 2L, minimum_unique_BODY_tuple_per_leaf = 2L, 
    donor_probability = "SELECTION_WEIGHT_NORMALIZED_WITHIN_LEAF", tuple = c("T", "q", "e", "hip", "waist"), uniform_slot = 3L, 
    unused_fixed_slots = 4:7, dimensions_unchanged = c(BASE = 11L, NMR = 49L, PROTEIN = 15L), assumptions = "Conditional exchangeability within coarsened predictive-mean cells; prediction-model/cell coarsening and finite donor-support limitations remain.", 
    not_guaranteed_unbiased = TRUE, not_exact_knn = TRUE, no_outcome_or_IE_tuning = TRUE, references = c("https://www.samhsa.gov/data/sites/default/files/NSDUH-PMNImputationMethods2015.pdf", 
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC3117228/", "https://www150.statcan.gc.ca/n1/pub/12-001-x/2005002/article/9041-eng.pdf", 
        "https://link.springer.com/article/10.1186/1471-2288-14-75"))

body_pmm_tree_accelerate <- function(tree) {
    tree$flat <- list(dim = vapply(tree$nodes, function(z) if (is.null(z$dim)) 0L else as.integer(z$dim), integer(1)), 
        cut = vapply(tree$nodes, function(z) if (is.null(z$cut)) NA_real_ else z$cut, numeric(1)), left = vapply(tree$nodes, 
            function(z) if (is.null(z$left)) 0L else as.integer(z$left), integer(1)), right = vapply(tree$nodes, 
            function(z) if (is.null(z$right)) 0L else as.integer(z$right), integer(1)), leaf = vapply(tree$nodes, 
            function(z) if (is.null(z$leaf)) 0L else as.integer(z$leaf), integer(1)))
    largest <- max(vapply(tree$leaves, `[[`, integer(1), "N"))
    cells <- largest * length(tree$leaves)
    tree$padded <- NULL
    if (largest <= 256L && cells <= 2000000L) {
        indices <- matrix(0L, length(tree$leaves), largest)
        cdf <- matrix(1, length(tree$leaves), largest)
        for (j in seq_along(tree$leaves)) {
            pool <- tree$leaves[[j]]
            cols <- seq_len(pool$N)
            indices[j, cols] <- pool$rows
            cdf[j, cols] <- pool$cdf[-1L]
        }
        tree$padded <- list(indices = indices, cdf = cdf)
    }
    tree$performance <- list(flat_depth_vectorized_query = TRUE, padded_donor_cdf = !is.null(tree$padded), largest_leaf = largest, 
        padding_limit_leaf = 256L, padding_limit_cells = 2000000L, original_nodes_donors_probabilities_unchanged = TRUE)
    tree
}

body_pmm_tree_fit <- function(scores, tuple, eid, w, spec) {
    nodes <- list()
    leaves <- list()
    n <- nrow(scores)
    make_leaf <- function(ix) {
        ix <- sort(ix)
        prob <- w[ix]/sum(w[ix])
        cdf <- c(0, cumsum(prob))
        cdf[length(cdf)] <- 1
        leaf <- length(leaves) + 1L
        unique_eid <- length(unique(eid[ix]))
        unique_tuple <- sum(!duplicated(as.data.frame(tuple[ix, , drop = FALSE])))
        body_tqe_assert(unique_eid >= spec$minimum_unique_eid_per_leaf && unique_tuple >= spec$minimum_unique_BODY_tuple_per_leaf, 
            "BODY_PMM_SINGLE_UNIQUE_DONOR_CELL_NOT_ALLOWED")
        leaves[[leaf]] <<- list(rows = ix, cdf = cdf, N = length(ix), unique_eid = unique_eid, unique_tuple = unique_tuple, 
            weight_sum = sum(w[ix]), ESS = sum(w[ix])^2/sum(w[ix]^2), score_min = apply(scores[ix, , drop = FALSE], 
                2, min), score_max = apply(scores[ix, , drop = FALSE], 2, max))
        list(leaf = leaf)
    }
    build <- function(ix) {
        node <- length(nodes) + 1L
        nodes[[node]] <<- list()
        if (length(ix) <= spec$maximum_leaf_rows) {
            nodes[[node]] <<- make_leaf(ix)
            return(node)
        }
        span <- apply(scores[ix, , drop = FALSE], 2, function(x) diff(range(x)))
        split <- NULL
        for (dim in order(span, decreasing = TRUE)) {
            if (span[dim] == 0) 
                next
            order_ix <- ix[order(scores[ix, dim], ix)]
            x <- scores[order_ix, dim]
            positions <- which(diff(x) > 0)
            positions <- positions[positions >= spec$minimum_child_rows & length(ix) - positions >= spec$minimum_child_rows]
            if (!length(positions)) 
                next
            k <- positions[which.min(abs(positions - length(ix)/2))]
            left <- order_ix[seq_len(k)]
            right <- order_ix[seq.int(k + 1L, length(ix))]
            if (length(unique(eid[left])) < 2L || length(unique(eid[right])) < 2L) 
                next
            if (sum(!duplicated(as.data.frame(tuple[left, , drop = FALSE]))) < 2L || sum(!duplicated(as.data.frame(tuple[right, 
                , drop = FALSE]))) < 2L) 
                next
            split <- list(dim = dim, cut = x[k], left = left, right = right)
            break
        }
        if (is.null(split)) {
            nodes[[node]] <<- make_leaf(ix)
            return(node)
        }
        left <- build(split$left)
        right <- build(split$right)
        nodes[[node]] <<- list(dim = split$dim, cut = split$cut, left = left, right = right)
        node
    }
    root <- build(seq_len(n))
    body_tqe_assert(length(leaves) > 1L, "BODY_PMM_NO_CONDITIONAL_PARTITION_UNCONDITIONAL_POOL_FORBIDDEN")
    stats <- do.call(rbind, lapply(leaves, function(z) c(N = z$N, unique_eid = z$unique_eid, unique_tuple = z$unique_tuple, 
        weight_sum = z$weight_sum, ESS = z$ESS)))
    body_pmm_tree_accelerate(list(nodes = nodes, leaves = leaves, root = root, leaf_statistics = stats, leaf_statistics_quantiles = apply(stats, 
        2, stats::quantile, probs = c(0, 0.01, 0.1, 0.5, 0.9, 0.99, 1)), leaf_count = length(leaves), node_count = length(nodes), 
        query_complexity = "DEPTH_VECTORIZED_FLAT_TREE_AND_BOUNDED_PADDED_CDF; NO_ALL_DONOR_SCAN_OR_TREE_REBUILD_AT_DRAW"))
}

body_pmm_tree_query <- function(tree, scores) {
    if (!is.null(tree$flat)) {
        answer <- integer(nrow(scores))
        node <- rep.int(tree$root, nrow(scores))
        active <- seq_len(nrow(scores))
        depth <- 0L
        while (length(active)) {
            depth <- depth + 1L
            body_tqe_assert(depth <= length(tree$nodes), "BODY_PMM_FLAT_TREE_CYCLE")
            at <- node[active]
            leaf <- tree$flat$leaf[at]
            done <- leaf > 0L
            answer[active[done]] <- leaf[done]
            active <- active[!done]
            at <- at[!done]
            if (length(active)) {
                left <- scores[cbind(active, tree$flat$dim[at])] <= tree$flat$cut[at]
                next_node <- tree$flat$right[at]
                next_node[left] <- tree$flat$left[at[left]]
                node[active] <- next_node
            }
        }
        body_tqe_assert(all(answer >= 1L), "BODY_PMM_FLAT_TREE_QUERY_FAILED")
        return(answer)
    }
    answer <- integer(nrow(scores))
    visit <- function(node, rows) {
        if (!length(rows)) 
            return(invisible(NULL))
        part <- tree$nodes[[node]]
        if (!is.null(part$leaf)) {
            answer[rows] <<- part$leaf
            return(invisible(NULL))
        }
        left <- scores[rows, part$dim] <= part$cut
        visit(part$left, rows[left])
        visit(part$right, rows[!left])
        invisible(NULL)
    }
    visit(tree$root, seq_len(nrow(scores)))
    body_tqe_assert(all(answer >= 1L), "BODY_PMM_TREE_QUERY_FAILED")
    answer
}

body_pmm_donor_query <- function(tree, leaf_ids, u) {
    if (!is.null(tree$padded)) {
        rank <- 1L + rowSums(u >= tree$padded$cdf[leaf_ids, , drop = FALSE])
        donor <- tree$padded$indices[cbind(leaf_ids, rank)]
        body_tqe_assert(all(donor >= 1L), "BODY_PMM_PADDED_DONOR_CDF_FAILURE_NO_RETRY")
        return(donor)
    }
    donor <- integer(length(leaf_ids))
    for (leaf in unique(leaf_ids)) {
        rows <- which(leaf_ids == leaf)
        pool <- tree$leaves[[leaf]]
        j <- findInterval(u[rows], pool$cdf, rightmost.closed = TRUE)
        body_tqe_assert(all(j >= 1L & j <= pool$N), "BODY_PMM_DONOR_CDF_FAILURE_NO_RETRY")
        donor[rows] <- pool$rows[j]
    }
    donor
}

body_pmm_group_fit <- function(d, members, base, w) {
    spec <- body_pmm_spec()
    enc <- body_tqe_encode(d, members)
    cols <- vapply(members[body_pmm_ids()], `[[`, "", "column")
    upstream <- vapply(members[c("smoking", "bmi")], `[[`, "", "column")
    design <- c39_a_design(d, base, upstream)
    X <- design$X
    Y <- cbind(logT = enc$z$logT, logitq = enc$z$logitq, e = enc$z$e, hip = d[[cols[["hip"]]]], waist = d[[cols[["waist"]]]])
    body_tqe_assert(qr(X)$rank == ncol(X), "BODY_PMM_MATCHING_DESIGN_RANK_DEFICIENT")
    fit <- stats::lm.wfit(X, Y, w)
    body_tqe_assert(fit$rank == ncol(X) && all(is.finite(fit$coefficients)), "BODY_PMM_MATCHING_WEIGHTED_FIT_FAILED")
    means <- X %*% fit$coefficients
    center <- colSums(means * w)/sum(w)
    scale <- sqrt(colSums(sweep(means, 2, center, "-")^2 * w)/sum(w))
    scale[scale == 0] <- 1
    scores <- sweep(sweep(means, 2, center, "-"), 2, scale, "/")
    tuple <- cbind(T = enc$total, q = enc$fraction, e = enc$z$e, hip = d[[cols[["hip"]]]], waist = d[[cols[["waist"]]]])
    eid <- if ("eid" %in% names(d)) 
        as.character(d$eid)
    else as.character(seq_len(nrow(d)))
    tree <- body_pmm_tree_fit(scores, tuple, eid, w, spec)
    list(revision = body_pmm_revision(), spec = spec, columns = cols, tuple = tuple, coefficients = fit$coefficients, 
        design_reference = design$reference, score_center = center, score_scale = scale, score_min = apply(scores, 
            2, min), score_max = apply(scores, 2, max), tree = tree, donor_N = nrow(d), unique_donor_eid = length(unique(eid)), 
        weight_ESS = sum(w)^2/sum(w^2), observed_precision = enc$precision, roundtrip = enc$roundtrip)
}

body_pmm_group_draw <- function(group, d, u, audit = NULL) {
    X <- c39_a_design(d, NULL, NULL, group$design_reference)$X
    scores <- sweep(sweep(X %*% group$coefficients, 2, group$score_center, "-"), 2, group$score_scale, "/")
    leaf_ids <- body_pmm_tree_query(group$tree, scores)
    donor <- body_pmm_donor_query(group$tree, leaf_ids, u)
    t <- group$tuple[donor, , drop = FALSE]
    F <- t[, "T"] * t[, "q"]
    L <- t[, "T"] * (1 - t[, "q"])
    W <- (F + L) + t[, "e"]
    body_tqe_assert(all(is.finite(F) & F > 0 & is.finite(L) & L > 0 & is.finite(W) & W > 0), "BODY_PMM_DONOR_RECONSTRUCTION_ILLEGAL_NO_REPAIR")
    values <- list(weight = W, fat_free_mass = L, fat_mass = F, hip = t[, "hip"], waist = t[, "waist"])
    for (id in names(values)) d[[group$columns[[id]]]] <- values[[id]]
    if (is.function(audit)) 
        audit("body_group_predictive_mean_cell", vapply(group$tree$leaves[leaf_ids], `[[`, integer(1), "N"), rep(1, 
            nrow(d)), vapply(group$tree$leaves[leaf_ids], `[[`, numeric(1), "ESS"))
    observer <- getOption("body_pmm.donor_audit_callback")
    if (is.function(observer)) {
        outside <- sweep(scores, 2, group$score_min, "<") | sweep(scores, 2, group$score_max, ">")
        observer(list(leaf_ids = leaf_ids, donor_rows = donor, scores = scores, outside_per_dimension = outside, 
            outside_global_score_range = rowSums(outside) > 0))
    }
    d
}

body_pmm_check_reused <- function(old, d, members, base, w, binding) {
    body_tqe_assert(is.list(old$A) && identical(old$binding, binding) && identical(old$A$body_tqe$prefix_binding, 
        binding) && identical(old$A$members, members[seq_len(11L)]) && is.list(old$provenance) && length(old$provenance) > 
        0L, "BODY_PMM_OLD_A11_EXACT_INPUT_WEIGHT_MEMBER_PROVENANCE_MISMATCH")
    ids <- names(members)[seq_len(11L)]
    cols <- body_tqe_columns(members)
    for (id in ids[c(1L, 2L, 8:11)]) {
        f <- old$A$models[[id]]
        before <- vapply(members[seq_len(match(id, ids) - 1L)], `[[`, "", "column")
        prepared <- body_tqe_affine_prepare(d, before, cols, reference = f$body_tqe_affine, w = w)
        expected <- c39_a_design(prepared$data, base, prepared$previous)
        actual <- c39_a_design(prepared$data, NULL, NULL, f$design_reference)
        body_tqe_assert(identical(f$member, members[[id]]) && identical(unname(f$predecessor_columns), unname(before)) && 
            identical(body_tqe_reference_copy(f$design_reference), body_tqe_reference_copy(expected$reference)) && 
            identical(actual$X, expected$X), "BODY_PMM_REUSED_UPSTREAM_OR_CLINICAL_DESIGN_CHANGED")
    }
    body_tqe_assert(identical(old$A$models$bmi$log_a_revision$revision, "CORE39_BMI_LOG_FAT_LINEAR_V1"), "BODY_PMM_REUSED_APPROVED_LOG_BMI_REQUIRED")
    invisible(TRUE)
}

body_pmm_fit_a <- function(d, members, base, w, continuous_rule = "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL", 
    count_rule = "FIXED_SUPPORT_MULTINOMIAL_0_7", prefix_fit = getOption("body_pmm.prefix_A"), old_fit = getOption("body_pmm.old_A11")) {
    ids <- vapply(members, `[[`, "", "indicator_id")
    body_tqe_assert(length(ids) >= 11L && identical(unname(ids[1:11]), body_tqe_ids()) && identical(names(members), 
        unname(ids)), "BODY_PMM_FROZEN_SELECTED11_PREFIX_REQUIRED")
    body_tqe_assert(length(w) == nrow(d) && all(is.finite(w) & w > 0) && is.finite(sum(w)), "BODY_PMM_INVALID_EXACT_SELECTION_WEIGHTS")
    body_tqe_support(d, members)
    binding <- body_tqe_binding(d, members, base, w)
    models <- list()
    previous <- character()
    warnings <- character()
    provenance <- NULL
    start <- 1L
    if (!is.null(prefix_fit)) {
        body_tqe_assert(identical(prefix_fit$body_pmm$revision, body_pmm_revision()) && identical(prefix_fit$body_pmm$prefix_binding, 
            binding) && identical(prefix_fit$members, members[1:11]), "BODY_PMM_NEW_PREFIX_BINDING_MISMATCH")
        models <- prefix_fit$models
        group <- prefix_fit$body_pmm$group
        start <- 12L
        previous <- vapply(members[1:11], `[[`, "", "column")
        warnings <- prefix_fit$warnings
        provenance <- prefix_fit$body_pmm$reused_old_provenance
    }
    else {
        if (!is.null(old_fit)) {
            body_pmm_check_reused(old_fit, d, members, base, w, binding)
            provenance <- old_fit$provenance
        }
        group <- NULL
    }
    if (start <= length(members)) 
        for (j in seq.int(start, length(members))) {
            member <- members[[j]]
            id <- member$indicator_id
            if (id %in% body_pmm_ids()) {
                if (id == "weight") 
                  group <- body_pmm_group_fit(d, members, base, w)
                f <- list(kind = "body_group_donor_slot", member = member, body_pmm_group_slot = TRUE, original_uniform_slot = j, 
                  uniform_used = id == "weight", e_law = "WHOLE_OBSERVED_TUPLE_CONDITIONAL_DONOR")
            }
            else if (!is.null(old_fit) && j %in% c(1L, 2L, 8:11)) 
                f <- old_fit$A$models[[id]]
            else f <- body_tqe_raw_fit(d, member, base, previous, w, body_tqe_columns(members), continuous_rule, 
                count_rule)
            models[[id]] <- f
            warnings <- c(warnings, f$warnings)
            previous <- c(previous, member$column)
            checkpoint <- getOption("c39.a_checkpoint_callback")
            if (is.function(checkpoint)) 
                checkpoint(id, f, j)
        }
    list(models = models, members = members, warnings = unique(warnings), conditional_order = names(models), continuous_rule = continuous_rule, 
        model_revision = body_pmm_revision(), body_pmm = list(revision = body_pmm_revision(), prefix_binding = binding, 
            group = group, reused_old_ids = if (is.null(provenance)) character() else body_tqe_ids()[c(1L, 2L, 
                8:11)], reused_old_provenance = provenance, unused_fixed_uniform_slots = 4:7, group_uniform_slot = 3L, 
            raw_members_unchanged = TRUE, no_height = TRUE), body_tqe = list(revision = body_tqe_revision(), generation_revision = body_pmm_revision(), 
            representation_only = TRUE, prefix_binding = binding, precision = group$observed_precision, roundtrip = group$roundtrip, 
            raw_columns = body_tqe_columns(members)), residual_assumption = "BODY5 common weighted conditional-cell donor; original upstream and clinical residual laws retained")
}

body_pmm_draw_a <- function(a, C, b, u, audit = NULL) {
    if (is.null(a$body_pmm)) 
        return(.body_pmm_previous_draw_a(a, C, b, u, audit))
    body_tqe_assert(identical(a$body_pmm$revision, body_pmm_revision()) && is.matrix(u) && nrow(u) == nrow(C) && 
        ncol(u) == length(a$models) && all(is.finite(u) & u >= 0 & u < 1) && length(b) %in% c(1L, nrow(C)) && all(is.finite(b) & 
        b >= 0 & b <= 1), "BODY_PMM_INVALID_RANDOM_INPUT")
    d <- C
    d$B <- b
    masscols <- a$body_pmm$group$columns[body_tqe_mass_ids()]
    for (j in seq_along(a$models)) {
        f <- a$models[[j]]
        id <- f$member$indicator_id
        if (isTRUE(f$body_pmm_group_slot)) {
            if (id == "weight") 
                d <- body_pmm_group_draw(a$body_pmm$group, d, u[, j], audit)
        }
        else {
            prepared <- body_tqe_affine_prepare(d, f$predecessor_columns, masscols, reference = f$body_tqe_affine)
            single <- list(models = setNames(list(f), id), members = setNames(list(f$member), id))
            d <- .body_tqe_original_draw_a(single, prepared$data, b, u[, j, drop = FALSE], audit)
        }
    }
    body_tqe_support(d, a$members)
    internal <- unname(body_tqe_names())
    d[intersect(setdiff(internal, names(C)), names(d))] <- NULL
    d
}

body_pmm_extend_saved_A <- function(prefix, old, d, base, w, provenance) {
    binding <- body_tqe_binding(d, old$members, base, w)
    body_tqe_assert(identical(prefix$body_pmm$revision, body_pmm_revision()) && identical(prefix$body_pmm$prefix_binding, 
        binding) && identical(old$body_tqe$prefix_binding, binding) && identical(prefix$members, old$members[1:11]) && 
        is.list(provenance) && length(provenance) > 0L, "BODY_PMM_SAVED_EXTENSION_INPUT_MEMBER_WEIGHT_PROVENANCE_MISMATCH")
    ids <- names(old$members)
    if (length(ids) > 11L) 
        for (j in seq.int(12L, length(ids))) {
            f <- old$models[[j]]
            previous <- vapply(old$members[seq_len(j - 1L)], `[[`, "", "column")
            prepared <- body_tqe_affine_prepare(d, previous, body_tqe_columns(old$members), w = w, reference = f$body_tqe_affine)
            expected <- c39_a_design(prepared$data, base, prepared$previous)
            actual <- c39_a_design(prepared$data, NULL, NULL, f$design_reference)
            body_tqe_assert(identical(f$member, old$members[[j]]) && identical(unname(f$predecessor_columns), unname(previous)) && 
                identical(body_tqe_reference_copy(f$design_reference), body_tqe_reference_copy(expected$reference)) && 
                identical(actual$X, expected$X), "BODY_PMM_SAVED_EXTENSION_ACTUAL_CONDITIONAL_DESIGN_CHANGED")
        }
    out <- old
    out$models[1:11] <- prefix$models
    out$body_pmm <- prefix$body_pmm
    out$body_tqe <- prefix$body_tqe
    out$model_revision <- body_pmm_revision()
    out$body_pmm$reused_extension_ids <- ids[-seq_len(11L)]
    out$body_pmm$reused_extension_provenance <- provenance
    body_tqe_assert(identical(out$models[-seq_len(11L)], old$models[-seq_len(11L)]), "BODY_PMM_SAVED_EXTENSION_OBJECT_CHANGED")
    out
}

body_pmm_diagnostics <- function(A11, d, w, engine, methods, seed, width = 11L, K = 3L) {
    previous <- getOption("body_pmm.donor_audit_callback")
    query <- list()
    options(body_pmm.donor_audit_callback = function(info) {
        query[[length(query) + 1L]] <<- list(call = length(query) + 1L, rows = length(info$leaf_ids), outside_global_score_range_n = sum(info$outside_global_score_range), 
            outside_global_score_range_weighted_mass = sum(w[info$outside_global_score_range])/sum(w), outside_n_by_score = colSums(info$outside_per_dimension), 
            distinct_selected_donor_rows = length(unique(info$donor_rows)), matched_leaf_count = length(unique(info$leaf_ids)))
        if (is.function(previous)) 
            previous(info)
    })
    on.exit(options(body_pmm.donor_audit_callback = previous), add = TRUE)
    out <- methods$body_tqe_diagnostics(A11, d, w, engine, methods, seed, width, K)
    group <- A11$body_pmm$group
    out$generation_revision <- body_pmm_revision()
    out$matching_spec <- group$spec
    out$donor_support <- list(N = group$donor_N, unique_eid = group$unique_donor_eid, selection_weight_ESS = group$weight_ESS, 
        leaf_count = group$tree$leaf_count, leaf_statistics_quantiles = group$tree$leaf_statistics_quantiles)
    out$query_support <- list(calls = query, call_order = "FOR_EACH_REPLICATE_OBSERVED_B_THEN_B0_THEN_B1", score_support = "MARGINAL_MIN_MAX_OF_TRAINING_PREDICTIVE_MEANS_NOT_A_CONVEX_HULL_TEST", 
        outside_handling = "DETERMINISTIC_BOUNDARY_CELL_ROUTING; NO_UNCONDITIONAL_FALLBACK_OR_SCORE_CLIPPING")
    out$no_new_ratio_or_correlation_threshold <- TRUE
    out
}

body_pmm_install <- function(engine) {
    body_tqe_assert(is.environment(engine) && exists(".body_tqe_installed", engine, inherits = FALSE), "BODY_PMM_INSTALL_TQE_ISOLATED_ENGINE_FIRST")
    if (exists(".body_pmm_installed", engine, inherits = FALSE)) {
        body_tqe_assert(identical(engine$.body_pmm_installed, body_pmm_revision()) && identical(engine$c39_fit_a, 
            engine$body_pmm_fit_a) && identical(engine$c39_draw_a, engine$body_pmm_draw_a), "BODY_PMM_INSTALLATION_DRIFT")
        return(invisible(engine))
    }
    engine$.body_pmm_previous_draw_a <- engine$c39_draw_a
    source_environment <- environment(body_pmm_install)
    for (name in ls(source_environment, pattern = "^body_pmm_")) {
        f <- get(name, source_environment, inherits = FALSE)
        if (is.function(f)) {
            environment(f) <- engine
            assign(name, f, engine)
        }
    }
    engine$c39_fit_a <- engine$body_pmm_fit_a
    engine$c39_draw_a <- engine$body_pmm_draw_a
    engine$fd_fit_a <- engine$c39_fit_a
    engine$fd_draw_a <- engine$c39_draw_a
    engine$.body_pmm_installed <- body_pmm_revision()
    invisible(engine)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        body_pmm_revision <- function() "SELECTED44_GROUPED_BODY_PREDICTIVE_MEAN_CELLS_V2"
        body_pmm_ids <- function() c("weight", "fat_free_mass", "fat_mass", "hip", "waist")
        body_pmm_spec <- function() list(method = "WEIGHTED_GROUPED_CONDITIONAL_HOTDECK_MULTIVARIATE_PREDICTIVE_MEAN_CELLS", 
            predictors = "ORIGINAL_B_C_SMOKING_BMI_ONLY", score_responses = c("logT", "logitq", "e", "hip", "waist"), 
            fitting = "EXACT_SELECTION_WEIGHTED_LEAST_SQUARES", score_metric = "WEIGHTED_SD_STANDARDIZED_PREDICTIVE_MEANS", 
            tree = "ONCE_FITTED_BALANCED_KD_PARTITION_MAX_LOCAL_SCORE_RANGE_MEDIAN_DISTINCT_SPLITS", maximum_leaf_rows = 128L, 
            minimum_child_rows = 32L, minimum_unique_eid_per_leaf = 2L, minimum_unique_BODY_tuple_per_leaf = 2L, 
            donor_probability = "SELECTION_WEIGHT_NORMALIZED_WITHIN_LEAF", tuple = c("T", "q", "e", "hip", "waist"), 
            uniform_slot = 3L, unused_fixed_slots = 4:7, dimensions_unchanged = c(BASE = 11L, NMR = 49L, PROTEIN = 15L), 
            assumptions = "Conditional exchangeability within coarsened predictive-mean cells; prediction-model/cell coarsening and finite donor-support limitations remain.", 
            not_guaranteed_unbiased = TRUE, not_exact_knn = TRUE, no_outcome_or_IE_tuning = TRUE, references = c("https://www.samhsa.gov/data/sites/default/files/NSDUH-PMNImputationMethods2015.pdf", 
                "https://pmc.ncbi.nlm.nih.gov/articles/PMC3117228/", "https://www150.statcan.gc.ca/n1/pub/12-001-x/2005002/article/9041-eng.pdf", 
                "https://link.springer.com/article/10.1186/1471-2288-14-75"))
        body_pmm_tree_accelerate <- function(tree) {
            tree$flat <- list(dim = vapply(tree$nodes, function(z) if (is.null(z$dim)) 0L else as.integer(z$dim), 
                integer(1)), cut = vapply(tree$nodes, function(z) if (is.null(z$cut)) NA_real_ else z$cut, numeric(1)), 
                left = vapply(tree$nodes, function(z) if (is.null(z$left)) 0L else as.integer(z$left), integer(1)), 
                right = vapply(tree$nodes, function(z) if (is.null(z$right)) 0L else as.integer(z$right), integer(1)), 
                leaf = vapply(tree$nodes, function(z) if (is.null(z$leaf)) 0L else as.integer(z$leaf), integer(1)))
            largest <- max(vapply(tree$leaves, `[[`, integer(1), "N"))
            cells <- largest * length(tree$leaves)
            tree$padded <- NULL
            if (largest <= 256L && cells <= 2000000L) {
                indices <- matrix(0L, length(tree$leaves), largest)
                cdf <- matrix(1, length(tree$leaves), largest)
                for (j in seq_along(tree$leaves)) {
                  pool <- tree$leaves[[j]]
                  cols <- seq_len(pool$N)
                  indices[j, cols] <- pool$rows
                  cdf[j, cols] <- pool$cdf[-1L]
                }
                tree$padded <- list(indices = indices, cdf = cdf)
            }
            tree$performance <- list(flat_depth_vectorized_query = TRUE, padded_donor_cdf = !is.null(tree$padded), 
                largest_leaf = largest, padding_limit_leaf = 256L, padding_limit_cells = 2000000L, original_nodes_donors_probabilities_unchanged = TRUE)
            tree
        }
        body_pmm_tree_fit <- function(scores, tuple, eid, w, spec) {
            nodes <- list()
            leaves <- list()
            n <- nrow(scores)
            make_leaf <- function(ix) {
                ix <- sort(ix)
                prob <- w[ix]/sum(w[ix])
                cdf <- c(0, cumsum(prob))
                cdf[length(cdf)] <- 1
                leaf <- length(leaves) + 1L
                unique_eid <- length(unique(eid[ix]))
                unique_tuple <- sum(!duplicated(as.data.frame(tuple[ix, , drop = FALSE])))
                body_tqe_assert(unique_eid >= spec$minimum_unique_eid_per_leaf && unique_tuple >= spec$minimum_unique_BODY_tuple_per_leaf, 
                  "BODY_PMM_SINGLE_UNIQUE_DONOR_CELL_NOT_ALLOWED")
                leaves[[leaf]] <<- list(rows = ix, cdf = cdf, N = length(ix), unique_eid = unique_eid, unique_tuple = unique_tuple, 
                  weight_sum = sum(w[ix]), ESS = sum(w[ix])^2/sum(w[ix]^2), score_min = apply(scores[ix, , drop = FALSE], 
                    2, min), score_max = apply(scores[ix, , drop = FALSE], 2, max))
                list(leaf = leaf)
            }
            build <- function(ix) {
                node <- length(nodes) + 1L
                nodes[[node]] <<- list()
                if (length(ix) <= spec$maximum_leaf_rows) {
                  nodes[[node]] <<- make_leaf(ix)
                  return(node)
                }
                span <- apply(scores[ix, , drop = FALSE], 2, function(x) diff(range(x)))
                split <- NULL
                for (dim in order(span, decreasing = TRUE)) {
                  if (span[dim] == 0) next
                  order_ix <- ix[order(scores[ix, dim], ix)]
                  x <- scores[order_ix, dim]
                  positions <- which(diff(x) > 0)
                  positions <- positions[positions >= spec$minimum_child_rows & length(ix) - positions >= spec$minimum_child_rows]
                  if (!length(positions)) next
                  k <- positions[which.min(abs(positions - length(ix)/2))]
                  left <- order_ix[seq_len(k)]
                  right <- order_ix[seq.int(k + 1L, length(ix))]
                  if (length(unique(eid[left])) < 2L || length(unique(eid[right])) < 2L) next
                  if (sum(!duplicated(as.data.frame(tuple[left, , drop = FALSE]))) < 2L || sum(!duplicated(as.data.frame(tuple[right, 
                    , drop = FALSE]))) < 2L) next
                  split <- list(dim = dim, cut = x[k], left = left, right = right)
                  break
                }
                if (is.null(split)) {
                  nodes[[node]] <<- make_leaf(ix)
                  return(node)
                }
                left <- build(split$left)
                right <- build(split$right)
                nodes[[node]] <<- list(dim = split$dim, cut = split$cut, left = left, right = right)
                node
            }
            root <- build(seq_len(n))
            body_tqe_assert(length(leaves) > 1L, "BODY_PMM_NO_CONDITIONAL_PARTITION_UNCONDITIONAL_POOL_FORBIDDEN")
            stats <- do.call(rbind, lapply(leaves, function(z) c(N = z$N, unique_eid = z$unique_eid, unique_tuple = z$unique_tuple, 
                weight_sum = z$weight_sum, ESS = z$ESS)))
            body_pmm_tree_accelerate(list(nodes = nodes, leaves = leaves, root = root, leaf_statistics = stats, 
                leaf_statistics_quantiles = apply(stats, 2, stats::quantile, probs = c(0, 0.01, 0.1, 0.5, 0.9, 
                  0.99, 1)), leaf_count = length(leaves), node_count = length(nodes), query_complexity = "DEPTH_VECTORIZED_FLAT_TREE_AND_BOUNDED_PADDED_CDF; NO_ALL_DONOR_SCAN_OR_TREE_REBUILD_AT_DRAW"))
        }
        body_pmm_tree_query <- function(tree, scores) {
            if (!is.null(tree$flat)) {
                answer <- integer(nrow(scores))
                node <- rep.int(tree$root, nrow(scores))
                active <- seq_len(nrow(scores))
                depth <- 0L
                while (length(active)) {
                  depth <- depth + 1L
                  body_tqe_assert(depth <= length(tree$nodes), "BODY_PMM_FLAT_TREE_CYCLE")
                  at <- node[active]
                  leaf <- tree$flat$leaf[at]
                  done <- leaf > 0L
                  answer[active[done]] <- leaf[done]
                  active <- active[!done]
                  at <- at[!done]
                  if (length(active)) {
                    left <- scores[cbind(active, tree$flat$dim[at])] <= tree$flat$cut[at]
                    next_node <- tree$flat$right[at]
                    next_node[left] <- tree$flat$left[at[left]]
                    node[active] <- next_node
                  }
                }
                body_tqe_assert(all(answer >= 1L), "BODY_PMM_FLAT_TREE_QUERY_FAILED")
                return(answer)
            }
            answer <- integer(nrow(scores))
            visit <- function(node, rows) {
                if (!length(rows)) return(invisible(NULL))
                part <- tree$nodes[[node]]
                if (!is.null(part$leaf)) {
                  answer[rows] <<- part$leaf
                  return(invisible(NULL))
                }
                left <- scores[rows, part$dim] <= part$cut
                visit(part$left, rows[left])
                visit(part$right, rows[!left])
                invisible(NULL)
            }
            visit(tree$root, seq_len(nrow(scores)))
            body_tqe_assert(all(answer >= 1L), "BODY_PMM_TREE_QUERY_FAILED")
            answer
        }
        body_pmm_donor_query <- function(tree, leaf_ids, u) {
            if (!is.null(tree$padded)) {
                rank <- 1L + rowSums(u >= tree$padded$cdf[leaf_ids, , drop = FALSE])
                donor <- tree$padded$indices[cbind(leaf_ids, rank)]
                body_tqe_assert(all(donor >= 1L), "BODY_PMM_PADDED_DONOR_CDF_FAILURE_NO_RETRY")
                return(donor)
            }
            donor <- integer(length(leaf_ids))
            for (leaf in unique(leaf_ids)) {
                rows <- which(leaf_ids == leaf)
                pool <- tree$leaves[[leaf]]
                j <- findInterval(u[rows], pool$cdf, rightmost.closed = TRUE)
                body_tqe_assert(all(j >= 1L & j <= pool$N), "BODY_PMM_DONOR_CDF_FAILURE_NO_RETRY")
                donor[rows] <- pool$rows[j]
            }
            donor
        }
        body_pmm_group_fit <- function(d, members, base, w) {
            spec <- body_pmm_spec()
            enc <- body_tqe_encode(d, members)
            cols <- vapply(members[body_pmm_ids()], `[[`, "", "column")
            upstream <- vapply(members[c("smoking", "bmi")], `[[`, "", "column")
            design <- c39_a_design(d, base, upstream)
            X <- design$X
            Y <- cbind(logT = enc$z$logT, logitq = enc$z$logitq, e = enc$z$e, hip = d[[cols[["hip"]]]], waist = d[[cols[["waist"]]]])
            body_tqe_assert(qr(X)$rank == ncol(X), "BODY_PMM_MATCHING_DESIGN_RANK_DEFICIENT")
            fit <- stats::lm.wfit(X, Y, w)
            body_tqe_assert(fit$rank == ncol(X) && all(is.finite(fit$coefficients)), "BODY_PMM_MATCHING_WEIGHTED_FIT_FAILED")
            means <- X %*% fit$coefficients
            center <- colSums(means * w)/sum(w)
            scale <- sqrt(colSums(sweep(means, 2, center, "-")^2 * w)/sum(w))
            scale[scale == 0] <- 1
            scores <- sweep(sweep(means, 2, center, "-"), 2, scale, "/")
            tuple <- cbind(T = enc$total, q = enc$fraction, e = enc$z$e, hip = d[[cols[["hip"]]]], waist = d[[cols[["waist"]]]])
            eid <- if ("eid" %in% names(d)) as.character(d$eid) else as.character(seq_len(nrow(d)))
            tree <- body_pmm_tree_fit(scores, tuple, eid, w, spec)
            list(revision = body_pmm_revision(), spec = spec, columns = cols, tuple = tuple, coefficients = fit$coefficients, 
                design_reference = design$reference, score_center = center, score_scale = scale, score_min = apply(scores, 
                  2, min), score_max = apply(scores, 2, max), tree = tree, donor_N = nrow(d), unique_donor_eid = length(unique(eid)), 
                weight_ESS = sum(w)^2/sum(w^2), observed_precision = enc$precision, roundtrip = enc$roundtrip)
        }
        body_pmm_group_draw <- function(group, d, u, audit = NULL) {
            X <- c39_a_design(d, NULL, NULL, group$design_reference)$X
            scores <- sweep(sweep(X %*% group$coefficients, 2, group$score_center, "-"), 2, group$score_scale, 
                "/")
            leaf_ids <- body_pmm_tree_query(group$tree, scores)
            donor <- body_pmm_donor_query(group$tree, leaf_ids, u)
            t <- group$tuple[donor, , drop = FALSE]
            F <- t[, "T"] * t[, "q"]
            L <- t[, "T"] * (1 - t[, "q"])
            W <- (F + L) + t[, "e"]
            body_tqe_assert(all(is.finite(F) & F > 0 & is.finite(L) & L > 0 & is.finite(W) & W > 0), "BODY_PMM_DONOR_RECONSTRUCTION_ILLEGAL_NO_REPAIR")
            values <- list(weight = W, fat_free_mass = L, fat_mass = F, hip = t[, "hip"], waist = t[, "waist"])
            for (id in names(values)) d[[group$columns[[id]]]] <- values[[id]]
            if (is.function(audit)) audit("body_group_predictive_mean_cell", vapply(group$tree$leaves[leaf_ids], 
                `[[`, integer(1), "N"), rep(1, nrow(d)), vapply(group$tree$leaves[leaf_ids], `[[`, numeric(1), 
                "ESS"))
            observer <- getOption("body_pmm.donor_audit_callback")
            if (is.function(observer)) {
                outside <- sweep(scores, 2, group$score_min, "<") | sweep(scores, 2, group$score_max, ">")
                observer(list(leaf_ids = leaf_ids, donor_rows = donor, scores = scores, outside_per_dimension = outside, 
                  outside_global_score_range = rowSums(outside) > 0))
            }
            d
        }
        body_pmm_check_reused <- function(old, d, members, base, w, binding) {
            body_tqe_assert(is.list(old$A) && identical(old$binding, binding) && identical(old$A$body_tqe$prefix_binding, 
                binding) && identical(old$A$members, members[seq_len(11L)]) && is.list(old$provenance) && length(old$provenance) > 
                0L, "BODY_PMM_OLD_A11_EXACT_INPUT_WEIGHT_MEMBER_PROVENANCE_MISMATCH")
            ids <- names(members)[seq_len(11L)]
            cols <- body_tqe_columns(members)
            for (id in ids[c(1L, 2L, 8:11)]) {
                f <- old$A$models[[id]]
                before <- vapply(members[seq_len(match(id, ids) - 1L)], `[[`, "", "column")
                prepared <- body_tqe_affine_prepare(d, before, cols, reference = f$body_tqe_affine, w = w)
                expected <- c39_a_design(prepared$data, base, prepared$previous)
                actual <- c39_a_design(prepared$data, NULL, NULL, f$design_reference)
                body_tqe_assert(identical(f$member, members[[id]]) && identical(unname(f$predecessor_columns), 
                  unname(before)) && identical(body_tqe_reference_copy(f$design_reference), body_tqe_reference_copy(expected$reference)) && 
                  identical(actual$X, expected$X), "BODY_PMM_REUSED_UPSTREAM_OR_CLINICAL_DESIGN_CHANGED")
            }
            body_tqe_assert(identical(old$A$models$bmi$log_a_revision$revision, "CORE39_BMI_LOG_FAT_LINEAR_V1"), 
                "BODY_PMM_REUSED_APPROVED_LOG_BMI_REQUIRED")
            invisible(TRUE)
        }
        body_pmm_fit_a <- function(d, members, base, w, continuous_rule = "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL", 
            count_rule = "FIXED_SUPPORT_MULTINOMIAL_0_7", prefix_fit = getOption("body_pmm.prefix_A"), old_fit = getOption("body_pmm.old_A11")) {
            ids <- vapply(members, `[[`, "", "indicator_id")
            body_tqe_assert(length(ids) >= 11L && identical(unname(ids[1:11]), body_tqe_ids()) && identical(names(members), 
                unname(ids)), "BODY_PMM_FROZEN_SELECTED11_PREFIX_REQUIRED")
            body_tqe_assert(length(w) == nrow(d) && all(is.finite(w) & w > 0) && is.finite(sum(w)), "BODY_PMM_INVALID_EXACT_SELECTION_WEIGHTS")
            body_tqe_support(d, members)
            binding <- body_tqe_binding(d, members, base, w)
            models <- list()
            previous <- character()
            warnings <- character()
            provenance <- NULL
            start <- 1L
            if (!is.null(prefix_fit)) {
                body_tqe_assert(identical(prefix_fit$body_pmm$revision, body_pmm_revision()) && identical(prefix_fit$body_pmm$prefix_binding, 
                  binding) && identical(prefix_fit$members, members[1:11]), "BODY_PMM_NEW_PREFIX_BINDING_MISMATCH")
                models <- prefix_fit$models
                group <- prefix_fit$body_pmm$group
                start <- 12L
                previous <- vapply(members[1:11], `[[`, "", "column")
                warnings <- prefix_fit$warnings
                provenance <- prefix_fit$body_pmm$reused_old_provenance
            } else {
                if (!is.null(old_fit)) {
                  body_pmm_check_reused(old_fit, d, members, base, w, binding)
                  provenance <- old_fit$provenance
                }
                group <- NULL
            }
            if (start <= length(members)) for (j in seq.int(start, length(members))) {
                member <- members[[j]]
                id <- member$indicator_id
                if (id %in% body_pmm_ids()) {
                  if (id == "weight") group <- body_pmm_group_fit(d, members, base, w)
                  f <- list(kind = "body_group_donor_slot", member = member, body_pmm_group_slot = TRUE, original_uniform_slot = j, 
                    uniform_used = id == "weight", e_law = "WHOLE_OBSERVED_TUPLE_CONDITIONAL_DONOR")
                } else if (!is.null(old_fit) && j %in% c(1L, 2L, 8:11)) f <- old_fit$A$models[[id]] else f <- body_tqe_raw_fit(d, 
                  member, base, previous, w, body_tqe_columns(members), continuous_rule, count_rule)
                models[[id]] <- f
                warnings <- c(warnings, f$warnings)
                previous <- c(previous, member$column)
                checkpoint <- getOption("c39.a_checkpoint_callback")
                if (is.function(checkpoint)) checkpoint(id, f, j)
            }
            list(models = models, members = members, warnings = unique(warnings), conditional_order = names(models), 
                continuous_rule = continuous_rule, model_revision = body_pmm_revision(), body_pmm = list(revision = body_pmm_revision(), 
                  prefix_binding = binding, group = group, reused_old_ids = if (is.null(provenance)) character() else body_tqe_ids()[c(1L, 
                    2L, 8:11)], reused_old_provenance = provenance, unused_fixed_uniform_slots = 4:7, group_uniform_slot = 3L, 
                  raw_members_unchanged = TRUE, no_height = TRUE), body_tqe = list(revision = body_tqe_revision(), 
                  generation_revision = body_pmm_revision(), representation_only = TRUE, prefix_binding = binding, 
                  precision = group$observed_precision, roundtrip = group$roundtrip, raw_columns = body_tqe_columns(members)), 
                residual_assumption = "BODY5 common weighted conditional-cell donor; original upstream and clinical residual laws retained")
        }
        body_pmm_draw_a <- function(a, C, b, u, audit = NULL) {
            if (is.null(a$body_pmm)) return(.body_pmm_previous_draw_a(a, C, b, u, audit))
            body_tqe_assert(identical(a$body_pmm$revision, body_pmm_revision()) && is.matrix(u) && nrow(u) == nrow(C) && 
                ncol(u) == length(a$models) && all(is.finite(u) & u >= 0 & u < 1) && length(b) %in% c(1L, nrow(C)) && 
                all(is.finite(b) & b >= 0 & b <= 1), "BODY_PMM_INVALID_RANDOM_INPUT")
            d <- C
            d$B <- b
            masscols <- a$body_pmm$group$columns[body_tqe_mass_ids()]
            for (j in seq_along(a$models)) {
                f <- a$models[[j]]
                id <- f$member$indicator_id
                if (isTRUE(f$body_pmm_group_slot)) {
                  if (id == "weight") d <- body_pmm_group_draw(a$body_pmm$group, d, u[, j], audit)
                } else {
                  prepared <- body_tqe_affine_prepare(d, f$predecessor_columns, masscols, reference = f$body_tqe_affine)
                  single <- list(models = setNames(list(f), id), members = setNames(list(f$member), id))
                  d <- .body_tqe_original_draw_a(single, prepared$data, b, u[, j, drop = FALSE], audit)
                }
            }
            body_tqe_support(d, a$members)
            internal <- unname(body_tqe_names())
            d[intersect(setdiff(internal, names(C)), names(d))] <- NULL
            d
        }
        body_pmm_extend_saved_A <- function(prefix, old, d, base, w, provenance) {
            binding <- body_tqe_binding(d, old$members, base, w)
            body_tqe_assert(identical(prefix$body_pmm$revision, body_pmm_revision()) && identical(prefix$body_pmm$prefix_binding, 
                binding) && identical(old$body_tqe$prefix_binding, binding) && identical(prefix$members, old$members[1:11]) && 
                is.list(provenance) && length(provenance) > 0L, "BODY_PMM_SAVED_EXTENSION_INPUT_MEMBER_WEIGHT_PROVENANCE_MISMATCH")
            ids <- names(old$members)
            if (length(ids) > 11L) for (j in seq.int(12L, length(ids))) {
                f <- old$models[[j]]
                previous <- vapply(old$members[seq_len(j - 1L)], `[[`, "", "column")
                prepared <- body_tqe_affine_prepare(d, previous, body_tqe_columns(old$members), w = w, reference = f$body_tqe_affine)
                expected <- c39_a_design(prepared$data, base, prepared$previous)
                actual <- c39_a_design(prepared$data, NULL, NULL, f$design_reference)
                body_tqe_assert(identical(f$member, old$members[[j]]) && identical(unname(f$predecessor_columns), 
                  unname(previous)) && identical(body_tqe_reference_copy(f$design_reference), body_tqe_reference_copy(expected$reference)) && 
                  identical(actual$X, expected$X), "BODY_PMM_SAVED_EXTENSION_ACTUAL_CONDITIONAL_DESIGN_CHANGED")
            }
            out <- old
            out$models[1:11] <- prefix$models
            out$body_pmm <- prefix$body_pmm
            out$body_tqe <- prefix$body_tqe
            out$model_revision <- body_pmm_revision()
            out$body_pmm$reused_extension_ids <- ids[-seq_len(11L)]
            out$body_pmm$reused_extension_provenance <- provenance
            body_tqe_assert(identical(out$models[-seq_len(11L)], old$models[-seq_len(11L)]), "BODY_PMM_SAVED_EXTENSION_OBJECT_CHANGED")
            out
        }
        body_pmm_diagnostics <- function(A11, d, w, engine, methods, seed, width = 11L, K = 3L) {
            previous <- getOption("body_pmm.donor_audit_callback")
            query <- list()
            options(body_pmm.donor_audit_callback = function(info) {
                query[[length(query) + 1L]] <<- list(call = length(query) + 1L, rows = length(info$leaf_ids), outside_global_score_range_n = sum(info$outside_global_score_range), 
                  outside_global_score_range_weighted_mass = sum(w[info$outside_global_score_range])/sum(w), outside_n_by_score = colSums(info$outside_per_dimension), 
                  distinct_selected_donor_rows = length(unique(info$donor_rows)), matched_leaf_count = length(unique(info$leaf_ids)))
                if (is.function(previous)) previous(info)
            })
            on.exit(options(body_pmm.donor_audit_callback = previous), add = TRUE)
            out <- methods$body_tqe_diagnostics(A11, d, w, engine, methods, seed, width, K)
            group <- A11$body_pmm$group
            out$generation_revision <- body_pmm_revision()
            out$matching_spec <- group$spec
            out$donor_support <- list(N = group$donor_N, unique_eid = group$unique_donor_eid, selection_weight_ESS = group$weight_ESS, 
                leaf_count = group$tree$leaf_count, leaf_statistics_quantiles = group$tree$leaf_statistics_quantiles)
            out$query_support <- list(calls = query, call_order = "FOR_EACH_REPLICATE_OBSERVED_B_THEN_B0_THEN_B1", 
                score_support = "MARGINAL_MIN_MAX_OF_TRAINING_PREDICTIVE_MEANS_NOT_A_CONVEX_HULL_TEST", outside_handling = "DETERMINISTIC_BOUNDARY_CELL_ROUTING; NO_UNCONDITIONAL_FALLBACK_OR_SCORE_CLIPPING")
            out$no_new_ratio_or_correlation_threshold <- TRUE
            out
        }
        body_pmm_install <- function(engine) {
            body_tqe_assert(is.environment(engine) && exists(".body_tqe_installed", engine, inherits = FALSE), 
                "BODY_PMM_INSTALL_TQE_ISOLATED_ENGINE_FIRST")
            if (exists(".body_pmm_installed", engine, inherits = FALSE)) {
                body_tqe_assert(identical(engine$.body_pmm_installed, body_pmm_revision()) && identical(engine$c39_fit_a, 
                  engine$body_pmm_fit_a) && identical(engine$c39_draw_a, engine$body_pmm_draw_a), "BODY_PMM_INSTALLATION_DRIFT")
                return(invisible(engine))
            }
            engine$.body_pmm_previous_draw_a <- engine$c39_draw_a
            source_environment <- environment(body_pmm_install)
            for (name in ls(source_environment, pattern = "^body_pmm_")) {
                f <- get(name, source_environment, inherits = FALSE)
                if (is.function(f)) {
                  environment(f) <- engine
                  assign(name, f, engine)
                }
            }
            engine$c39_fit_a <- engine$body_pmm_fit_a
            engine$c39_draw_a <- engine$body_pmm_draw_a
            engine$fd_fit_a <- engine$c39_fit_a
            engine$fd_draw_a <- engine$c39_draw_a
            engine$.body_pmm_installed <- body_pmm_revision()
            invisible(engine)
        }
    }), envir = context)
}

