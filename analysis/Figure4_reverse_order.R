# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

omics_product_stream <- function(C, weights, domains, draw_joint, predict_risk, validate_support, uniform_dim, 
    n_mc = 50L, block_size = 512L, seed = 1L, allocation_order = names(domains), include_reverse = TRUE, endpoints = c("IHD", 
        "T2D"), binding, checkpoint = NULL, save_checkpoint = NULL, progress = NULL) {
    spec <- new121_spec(C, weights, domains, allocation_order, endpoints)
    if (!all(vapply(list(draw_joint, predict_risk, validate_support), is.function, logical(1))) || (!is.null(save_checkpoint) && 
        !is.function(save_checkpoint)) || (!is.null(progress) && !is.function(progress))) 
        new121_abort("invalid callback")
    for (key in c("uniform_dim", "n_mc", "block_size", "seed")) {
        value <- get(key)
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 1 || value > .Machine$integer.max || 
            value != floor(value)) 
            new121_abort(paste(key, "must be a positive integer"))
    }
    if (!is.logical(include_reverse) || length(include_reverse) != 1L || is.na(include_reverse)) 
        new121_abort("include_reverse must be scalar logical")
    required <- c("panel_id", "bootstrap_id", "sample_hash", "population_hash", "model_spec_hash", "input_hash", 
        "rng_spec_hash")
    if (!is.list(binding) || !all(required %in% names(binding)) || anyDuplicated(names(binding))) 
        new121_abort(paste("scientific binding requires", paste(required, collapse = ", ")))
    if (!is.numeric(binding$bootstrap_id) || length(binding$bootstrap_id) != 1L || !is.finite(binding$bootstrap_id) || 
        binding$bootstrap_id < 0 || binding$bootstrap_id > 100 || binding$bootstrap_id != floor(binding$bootstrap_id)) 
        new121_abort("fixed bootstrap ID 0..100 required")
    for (key in setdiff(required, "bootstrap_id")) if (!is.character(binding[[key]]) || length(binding[[key]]) != 
        1L || is.na(binding[[key]]) || !nzchar(binding[[key]])) 
        new121_abort(paste("invalid scientific binding", key))
    binding <- binding[sort(required)]
    orders <- list(forward = allocation_order)
    if (include_reverse) 
        orders$reverse <- rev(allocation_order)
    zero <- lapply(orders, function(order) matrix(0, length(new121_scenarios(order)), length(endpoints), dimnames = list(new121_scenarios(order), 
        endpoints)))
    identity <- list(schema = "OMICS_PRODUCT_MC_REVERSE_ISOLATION_V1", binding = binding, domains = domains, orders = orders, 
        endpoints = endpoints, n_target = nrow(C), C_columns = names(C), uniform_dim = as.integer(uniform_dim), 
        n_mc = as.integer(n_mc), seed = as.integer(seed), weight_sum = sum(weights), stream_contract = "JOINT_PLUS_ONE_INDEPENDENT_STREAM_PER_DOMAIN_PAIRED_B0_B1", 
        reverse_failure_policy = "DROP_ALL_REVERSE_ESTIMATES_KEEP_COMPLETE_FORWARD_NO_JOINT_REPLAY")
    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) 
        old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
        do.call(RNGkind, as.list(old_kind))
        if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, 
            inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    stream_names <- c("JOINT", paste0("DOMAIN:", names(domains)))
    states <- list(get(".Random.seed", envir = .GlobalEnv))
    for (j in seq_len(length(domains))) states[[j + 1L]] <- parallel::nextRNGStream(states[[j]])
    names(states) <- stream_names
    initial_states <- states
    completed <- 0L
    mc_risks <- list()
    prediction_calls <- 0L
    resumed <- FALSE
    reverse_failure <- NULL
    if (!is.null(checkpoint)) {
        old_identity <- checkpoint$identity
        if (!is.list(old_identity)) 
            new121_abort("missing checkpoint identity")
        old_sum <- old_identity$weight_sum
        old_identity$weight_sum <- NULL
        compare_identity <- identity
        compare_identity$weight_sum <- NULL
        if (!identical(old_identity, compare_identity) || !isTRUE(all.equal(old_sum, identity$weight_sum, tolerance = 1e-10))) 
            new121_abort("checkpoint scientific/sample/stream identity mismatch")
        completed <- checkpoint$completed_mc
        mc_risks <- checkpoint$mc_risks
        states <- checkpoint$stream_states
        if (!is.numeric(completed) || length(completed) != 1L || !is.finite(completed) || completed < 0 || completed > 
            n_mc || completed != floor(completed) || length(mc_risks) != completed || !identical(names(states), 
            stream_names) || any(!vapply(states, function(s) is.integer(s) && length(s) == 7L && !anyNA(s), logical(1)))) 
            new121_abort("corrupt checkpoint: completed-pass accounting/streams")
        flag <- checkpoint$reverse_failure_recorded
        if (!is.logical(flag) || length(flag) != 1L || is.na(flag) || !identical(flag, !is.null(checkpoint$reverse_failure))) 
            new121_abort("checkpoint explicit reverse-failure flag missing or inconsistent")
        reverse_failure <- checkpoint$reverse_failure
        if (!is.null(reverse_failure) && (!include_reverse || !is.list(reverse_failure) || !identical(reverse_failure$allocation_order, 
            "reverse") || !is.character(reverse_failure$message) || length(reverse_failure$message) != 1L || !nzchar(reverse_failure$message) || 
            !is.numeric(reverse_failure$mc) || length(reverse_failure$mc) != 1L || reverse_failure$mc < 1L || reverse_failure$mc > 
            completed || reverse_failure$mc != floor(reverse_failure$mc))) 
            new121_abort("corrupt checkpoint reverse-failure diagnostic")
        active <- if (is.null(reverse_failure)) 
            names(orders)
        else "forward"
        for (pass in mc_risks) {
            if (!is.list(pass) || !identical(names(pass), active)) 
                new121_abort("checkpoint failed reverse estimates or order placeholders forbidden")
            for (ord in active) {
                x <- pass[[ord]]
                if (!is.matrix(x) || !identical(dimnames(x), dimnames(zero[[ord]])) || any(!is.finite(x)) || any(x < 
                  -1e-12 | x > 1 + 1e-12)) 
                  new121_abort("corrupt checkpoint risk means")
            }
        }
        prediction_calls <- checkpoint$prediction_calls
        if (length(prediction_calls) != 1L || !is.finite(prediction_calls) || prediction_calls < 0) 
            new121_abort("corrupt prediction count")
        resumed <- completed > 0L
    }
    next_uniforms <- function(stream, n) {
        assign(".Random.seed", states[[stream]], envir = .GlobalEnv)
        z <- matrix(runif(n * uniform_dim), nrow = n, byrow = TRUE)
        states[[stream]] <<- get(".Random.seed", envir = .GlobalEnv)
        z
    }
    deterministic <- function(callback) {
        before <- get(".Random.seed", envir = .GlobalEnv)
        value <- callback()
        if (!identical(before, get(".Random.seed", envir = .GlobalEnv))) 
            new121_abort("simulation/prediction/support callbacks must use supplied uniforms only")
        value
    }
    started <- proc.time()[["elapsed"]]
    if (completed < n_mc) 
        for (k in seq.int(completed + 1L, n_mc)) {
            sums <- zero
            for (lo in seq.int(1L, nrow(C), by = block_size)) {
                ii <- seq.int(lo, min(nrow(C), lo + block_size - 1L))
                cc <- C[ii, , drop = FALSE]
                ww <- weights[ii]
                context <- paste0("MC=", k, ";rows=", lo, "-", max(ii))
                check <- function(M, label) deterministic(function() new121_check_M(M, cc, spec$columns, validate_support, 
                  paste(context, label, sep = ";")))
                draw <- function(b, u, label) check(deterministic(function() draw_joint(b = b, C = cc, uniforms = u)), 
                  label)
                predict <- function(b, M, label) {
                  prediction_calls <<- prediction_calls + 1L
                  colSums(deterministic(function() new121_predict(b, M, cc, predict_risk, endpoints, paste(context, 
                    label, sep = ";"))) * ww)
                }
                u <- next_uniforms("JOINT", length(ii))
                Q0 <- draw(0, u, "jointQ0")
                Q1 <- draw(1, u, "jointQ1")
                overall <- rbind(R00 = predict(0, Q0, "R00"), R10 = predict(1, Q0, "R10"), R11 = predict(1, Q1, 
                  "R11"))
                P0 <- Q0
                P1 <- Q1
                for (g in names(domains)) {
                  ug <- next_uniforms(paste0("DOMAIN:", g), length(ii))
                  D0 <- draw(0, ug, paste0("domain=", g, ";Q0"))
                  D1 <- draw(1, ug, paste0("domain=", g, ";Q1"))
                  cols <- domains[[g]]
                  P0[, cols] <- D0[, cols, drop = FALSE]
                  P1[, cols] <- D1[, cols, drop = FALSE]
                }
                P0 <- check(P0, "all0_conditional_product")
                P1 <- check(P1, "all1_conditional_product")
                endpoints_product <- rbind(P0 = predict(1, P0, "P0"), P1 = predict(1, P1, "P1"))
                for (ord in names(orders)) {
                  if (ord == "reverse" && !is.null(reverse_failure)) 
                    next
                  attempt <- tryCatch({
                    sums[[ord]][c("R00", "R10", "R11"), ] <- sums[[ord]][c("R00", "R10", "R11"), ] + overall
                    sums[[ord]][c("P0", "P1"), ] <- sums[[ord]][c("P0", "P1"), ] + endpoints_product
                    M <- P0
                    for (j in seq_along(orders[[ord]])) {
                      g <- orders[[ord]][j]
                      cols <- domains[[g]]
                      M[, cols] <- P1[, cols, drop = FALSE]
                      scenario <- paste0("S:", g)
                      value <- if (j == length(orders[[ord]])) 
                        endpoints_product["P1", ]
                      else predict(1, check(M, paste(ord, scenario)), paste(ord, scenario))
                      sums[[ord]][scenario, ] <- sums[[ord]][scenario, ] + value
                    }
                    NULL
                  }, error = function(e) {
                    if (ord != "reverse") 
                      stop(e)
                    list(allocation_order = "reverse", mc = as.integer(k), rows_first = as.integer(lo), rows_last = as.integer(max(ii)), 
                      context = context, message = conditionMessage(e), condition_class = class(e), failed_stage = "REVERSE_CUMULATIVE_REPLACEMENT", 
                      disposition = "ALL_REVERSE_ESTIMATES_DROPPED_FORWARD_CONTINUES")
                  })
                  if (!is.null(attempt)) {
                    reverse_failure <- attempt
                    mc_risks <- lapply(mc_risks, function(pass) pass["forward"])
                  }
                }
                if (!is.null(progress)) 
                  progress(list(stage = "integrating", bootstrap_id = binding$bootstrap_id, mc = k, n_mc = n_mc, 
                    rows_complete = max(ii), n_target = nrow(C), prediction_calls = prediction_calls, reverse_failure_recorded = !is.null(reverse_failure), 
                    scenarios_per_block = 5L + sum(vapply(orders[if (is.null(reverse_failure)) names(orders) else "forward"], 
                      length, integer(1)) - 1L)))
            }
            active <- if (is.null(reverse_failure)) 
                names(orders)
            else "forward"
            mc_risks[[k]] <- lapply(sums[active], function(x) x/sum(weights))
            completed <- k
            state <- list(identity = identity, completed_mc = completed, mc_risks = mc_risks, stream_states = states, 
                prediction_calls = prediction_calls, reverse_failure_recorded = !is.null(reverse_failure), reverse_failure = reverse_failure)
            if (!is.null(save_checkpoint)) 
                save_checkpoint(state)
        }
    active <- if (is.null(reverse_failure)) 
        names(orders)
    else "forward"
    returned_orders <- orders[active]
    risks <- lapply(active, function(ord) Reduce(`+`, lapply(mc_risks, `[[`, ord))/n_mc)
    names(risks) <- active
    result <- new121_finish(risks, returned_orders, endpoints, mc_risks, list(schema = "NEW121_PRODUCT_RESULT_V1", 
        method = "conditional independent domain marginals; cumulative replacement", binding = binding, n_target = nrow(C), 
        weight_sum = sum(weights), n_mc = completed, block_size = block_size, uniform_dim = uniform_dim, seed = seed, 
        initial_stream_states = initial_states, final_stream_states = states, stream_contract = identity$stream_contract, 
        same_C_within_all_domain_draws = TRUE, independent_domain_streams = TRUE, paired_uniforms_at_exposure_0_1 = TRUE, 
        R_scenarios_full_joint = TRUE, shared_product_endpoints_both_orders = TRUE, resumed = resumed, completed_mc_on_entry = if (resumed) checkpoint$completed_mc else 0L, 
        prediction_calls = prediction_calls, model_fits = 0L, requested_orders = names(orders), returned_orders = active, 
        reverse_failure_recorded = !is.null(reverse_failure), reverse_failure = reverse_failure, reverse_failure_is_not_forward_failure = TRUE, 
        reverse_partial_estimates_or_placeholders_returned = FALSE, elapsed_seconds = proc.time()[["elapsed"]] - 
            started, units = "probability and explicit per_1000 columns", cross_device_probability_reference = 1e-06, 
        bitwise_cross_device_agreement_required = FALSE))
    result$reverse_failure <- reverse_failure
    result
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        omics_product_stream <- function(C, weights, domains, draw_joint, predict_risk, validate_support, uniform_dim, 
            n_mc = 50L, block_size = 512L, seed = 1L, allocation_order = names(domains), include_reverse = TRUE, 
            endpoints = c("IHD", "T2D"), binding, checkpoint = NULL, save_checkpoint = NULL, progress = NULL) {
            spec <- new121_spec(C, weights, domains, allocation_order, endpoints)
            if (!all(vapply(list(draw_joint, predict_risk, validate_support), is.function, logical(1))) || (!is.null(save_checkpoint) && 
                !is.function(save_checkpoint)) || (!is.null(progress) && !is.function(progress))) new121_abort("invalid callback")
            for (key in c("uniform_dim", "n_mc", "block_size", "seed")) {
                value <- get(key)
                if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 1 || value > .Machine$integer.max || 
                  value != floor(value)) new121_abort(paste(key, "must be a positive integer"))
            }
            if (!is.logical(include_reverse) || length(include_reverse) != 1L || is.na(include_reverse)) new121_abort("include_reverse must be scalar logical")
            required <- c("panel_id", "bootstrap_id", "sample_hash", "population_hash", "model_spec_hash", "input_hash", 
                "rng_spec_hash")
            if (!is.list(binding) || !all(required %in% names(binding)) || anyDuplicated(names(binding))) new121_abort(paste("scientific binding requires", 
                paste(required, collapse = ", ")))
            if (!is.numeric(binding$bootstrap_id) || length(binding$bootstrap_id) != 1L || !is.finite(binding$bootstrap_id) || 
                binding$bootstrap_id < 0 || binding$bootstrap_id > 100 || binding$bootstrap_id != floor(binding$bootstrap_id)) new121_abort("fixed bootstrap ID 0..100 required")
            for (key in setdiff(required, "bootstrap_id")) if (!is.character(binding[[key]]) || length(binding[[key]]) != 
                1L || is.na(binding[[key]]) || !nzchar(binding[[key]])) new121_abort(paste("invalid scientific binding", 
                key))
            binding <- binding[sort(required)]
            orders <- list(forward = allocation_order)
            if (include_reverse) orders$reverse <- rev(allocation_order)
            zero <- lapply(orders, function(order) matrix(0, length(new121_scenarios(order)), length(endpoints), 
                dimnames = list(new121_scenarios(order), endpoints)))
            identity <- list(schema = "OMICS_PRODUCT_MC_REVERSE_ISOLATION_V1", binding = binding, domains = domains, 
                orders = orders, endpoints = endpoints, n_target = nrow(C), C_columns = names(C), uniform_dim = as.integer(uniform_dim), 
                n_mc = as.integer(n_mc), seed = as.integer(seed), weight_sum = sum(weights), stream_contract = "JOINT_PLUS_ONE_INDEPENDENT_STREAM_PER_DOMAIN_PAIRED_B0_B1", 
                reverse_failure_policy = "DROP_ALL_REVERSE_ESTIMATES_KEEP_COMPLETE_FORWARD_NO_JOINT_REPLAY")
            old_kind <- RNGkind()
            had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
            if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
            on.exit({
                do.call(RNGkind, as.list(old_kind))
                if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", 
                  envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
            }, add = TRUE)
            RNGkind("L'Ecuyer-CMRG")
            set.seed(seed)
            stream_names <- c("JOINT", paste0("DOMAIN:", names(domains)))
            states <- list(get(".Random.seed", envir = .GlobalEnv))
            for (j in seq_len(length(domains))) states[[j + 1L]] <- parallel::nextRNGStream(states[[j]])
            names(states) <- stream_names
            initial_states <- states
            completed <- 0L
            mc_risks <- list()
            prediction_calls <- 0L
            resumed <- FALSE
            reverse_failure <- NULL
            if (!is.null(checkpoint)) {
                old_identity <- checkpoint$identity
                if (!is.list(old_identity)) new121_abort("missing checkpoint identity")
                old_sum <- old_identity$weight_sum
                old_identity$weight_sum <- NULL
                compare_identity <- identity
                compare_identity$weight_sum <- NULL
                if (!identical(old_identity, compare_identity) || !isTRUE(all.equal(old_sum, identity$weight_sum, 
                  tolerance = 1e-10))) new121_abort("checkpoint scientific/sample/stream identity mismatch")
                completed <- checkpoint$completed_mc
                mc_risks <- checkpoint$mc_risks
                states <- checkpoint$stream_states
                if (!is.numeric(completed) || length(completed) != 1L || !is.finite(completed) || completed < 0 || 
                  completed > n_mc || completed != floor(completed) || length(mc_risks) != completed || !identical(names(states), 
                  stream_names) || any(!vapply(states, function(s) is.integer(s) && length(s) == 7L && !anyNA(s), 
                  logical(1)))) new121_abort("corrupt checkpoint: completed-pass accounting/streams")
                flag <- checkpoint$reverse_failure_recorded
                if (!is.logical(flag) || length(flag) != 1L || is.na(flag) || !identical(flag, !is.null(checkpoint$reverse_failure))) new121_abort("checkpoint explicit reverse-failure flag missing or inconsistent")
                reverse_failure <- checkpoint$reverse_failure
                if (!is.null(reverse_failure) && (!include_reverse || !is.list(reverse_failure) || !identical(reverse_failure$allocation_order, 
                  "reverse") || !is.character(reverse_failure$message) || length(reverse_failure$message) != 1L || 
                  !nzchar(reverse_failure$message) || !is.numeric(reverse_failure$mc) || length(reverse_failure$mc) != 
                  1L || reverse_failure$mc < 1L || reverse_failure$mc > completed || reverse_failure$mc != floor(reverse_failure$mc))) new121_abort("corrupt checkpoint reverse-failure diagnostic")
                active <- if (is.null(reverse_failure)) names(orders) else "forward"
                for (pass in mc_risks) {
                  if (!is.list(pass) || !identical(names(pass), active)) new121_abort("checkpoint failed reverse estimates or order placeholders forbidden")
                  for (ord in active) {
                    x <- pass[[ord]]
                    if (!is.matrix(x) || !identical(dimnames(x), dimnames(zero[[ord]])) || any(!is.finite(x)) || 
                      any(x < -1e-12 | x > 1 + 1e-12)) new121_abort("corrupt checkpoint risk means")
                  }
                }
                prediction_calls <- checkpoint$prediction_calls
                if (length(prediction_calls) != 1L || !is.finite(prediction_calls) || prediction_calls < 0) new121_abort("corrupt prediction count")
                resumed <- completed > 0L
            }
            next_uniforms <- function(stream, n) {
                assign(".Random.seed", states[[stream]], envir = .GlobalEnv)
                z <- matrix(runif(n * uniform_dim), nrow = n, byrow = TRUE)
                states[[stream]] <<- get(".Random.seed", envir = .GlobalEnv)
                z
            }
            deterministic <- function(callback) {
                before <- get(".Random.seed", envir = .GlobalEnv)
                value <- callback()
                if (!identical(before, get(".Random.seed", envir = .GlobalEnv))) new121_abort("simulation/prediction/support callbacks must use supplied uniforms only")
                value
            }
            started <- proc.time()[["elapsed"]]
            if (completed < n_mc) for (k in seq.int(completed + 1L, n_mc)) {
                sums <- zero
                for (lo in seq.int(1L, nrow(C), by = block_size)) {
                  ii <- seq.int(lo, min(nrow(C), lo + block_size - 1L))
                  cc <- C[ii, , drop = FALSE]
                  ww <- weights[ii]
                  context <- paste0("MC=", k, ";rows=", lo, "-", max(ii))
                  check <- function(M, label) deterministic(function() new121_check_M(M, cc, spec$columns, validate_support, 
                    paste(context, label, sep = ";")))
                  draw <- function(b, u, label) check(deterministic(function() draw_joint(b = b, C = cc, uniforms = u)), 
                    label)
                  predict <- function(b, M, label) {
                    prediction_calls <<- prediction_calls + 1L
                    colSums(deterministic(function() new121_predict(b, M, cc, predict_risk, endpoints, paste(context, 
                      label, sep = ";"))) * ww)
                  }
                  u <- next_uniforms("JOINT", length(ii))
                  Q0 <- draw(0, u, "jointQ0")
                  Q1 <- draw(1, u, "jointQ1")
                  overall <- rbind(R00 = predict(0, Q0, "R00"), R10 = predict(1, Q0, "R10"), R11 = predict(1, Q1, 
                    "R11"))
                  P0 <- Q0
                  P1 <- Q1
                  for (g in names(domains)) {
                    ug <- next_uniforms(paste0("DOMAIN:", g), length(ii))
                    D0 <- draw(0, ug, paste0("domain=", g, ";Q0"))
                    D1 <- draw(1, ug, paste0("domain=", g, ";Q1"))
                    cols <- domains[[g]]
                    P0[, cols] <- D0[, cols, drop = FALSE]
                    P1[, cols] <- D1[, cols, drop = FALSE]
                  }
                  P0 <- check(P0, "all0_conditional_product")
                  P1 <- check(P1, "all1_conditional_product")
                  endpoints_product <- rbind(P0 = predict(1, P0, "P0"), P1 = predict(1, P1, "P1"))
                  for (ord in names(orders)) {
                    if (ord == "reverse" && !is.null(reverse_failure)) next
                    attempt <- tryCatch({
                      sums[[ord]][c("R00", "R10", "R11"), ] <- sums[[ord]][c("R00", "R10", "R11"), ] + overall
                      sums[[ord]][c("P0", "P1"), ] <- sums[[ord]][c("P0", "P1"), ] + endpoints_product
                      M <- P0
                      for (j in seq_along(orders[[ord]])) {
                        g <- orders[[ord]][j]
                        cols <- domains[[g]]
                        M[, cols] <- P1[, cols, drop = FALSE]
                        scenario <- paste0("S:", g)
                        value <- if (j == length(orders[[ord]])) endpoints_product["P1", ] else predict(1, check(M, 
                          paste(ord, scenario)), paste(ord, scenario))
                        sums[[ord]][scenario, ] <- sums[[ord]][scenario, ] + value
                      }
                      NULL
                    }, error = function(e) {
                      if (ord != "reverse") stop(e)
                      list(allocation_order = "reverse", mc = as.integer(k), rows_first = as.integer(lo), rows_last = as.integer(max(ii)), 
                        context = context, message = conditionMessage(e), condition_class = class(e), failed_stage = "REVERSE_CUMULATIVE_REPLACEMENT", 
                        disposition = "ALL_REVERSE_ESTIMATES_DROPPED_FORWARD_CONTINUES")
                    })
                    if (!is.null(attempt)) {
                      reverse_failure <- attempt
                      mc_risks <- lapply(mc_risks, function(pass) pass["forward"])
                    }
                  }
                  if (!is.null(progress)) progress(list(stage = "integrating", bootstrap_id = binding$bootstrap_id, 
                    mc = k, n_mc = n_mc, rows_complete = max(ii), n_target = nrow(C), prediction_calls = prediction_calls, 
                    reverse_failure_recorded = !is.null(reverse_failure), scenarios_per_block = 5L + sum(vapply(orders[if (is.null(reverse_failure)) names(orders) else "forward"], 
                      length, integer(1)) - 1L)))
                }
                active <- if (is.null(reverse_failure)) names(orders) else "forward"
                mc_risks[[k]] <- lapply(sums[active], function(x) x/sum(weights))
                completed <- k
                state <- list(identity = identity, completed_mc = completed, mc_risks = mc_risks, stream_states = states, 
                  prediction_calls = prediction_calls, reverse_failure_recorded = !is.null(reverse_failure), reverse_failure = reverse_failure)
                if (!is.null(save_checkpoint)) save_checkpoint(state)
            }
            active <- if (is.null(reverse_failure)) names(orders) else "forward"
            returned_orders <- orders[active]
            risks <- lapply(active, function(ord) Reduce(`+`, lapply(mc_risks, `[[`, ord))/n_mc)
            names(risks) <- active
            result <- new121_finish(risks, returned_orders, endpoints, mc_risks, list(schema = "NEW121_PRODUCT_RESULT_V1", 
                method = "conditional independent domain marginals; cumulative replacement", binding = binding, 
                n_target = nrow(C), weight_sum = sum(weights), n_mc = completed, block_size = block_size, uniform_dim = uniform_dim, 
                seed = seed, initial_stream_states = initial_states, final_stream_states = states, stream_contract = identity$stream_contract, 
                same_C_within_all_domain_draws = TRUE, independent_domain_streams = TRUE, paired_uniforms_at_exposure_0_1 = TRUE, 
                R_scenarios_full_joint = TRUE, shared_product_endpoints_both_orders = TRUE, resumed = resumed, 
                completed_mc_on_entry = if (resumed) checkpoint$completed_mc else 0L, prediction_calls = prediction_calls, 
                model_fits = 0L, requested_orders = names(orders), returned_orders = active, reverse_failure_recorded = !is.null(reverse_failure), 
                reverse_failure = reverse_failure, reverse_failure_is_not_forward_failure = TRUE, reverse_partial_estimates_or_placeholders_returned = FALSE, 
                elapsed_seconds = proc.time()[["elapsed"]] - started, units = "probability and explicit per_1000 columns", 
                cross_device_probability_reference = 1e-06, bitwise_cross_device_agreement_required = FALSE))
            result$reverse_failure <- reverse_failure
            result
        }
    }), envir = context)
}

