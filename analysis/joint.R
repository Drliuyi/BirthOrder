# NEW121 product-of-conditional-domain-marginals kernel. No data loading/fitting.
# R00/R10/R11 retain the complete joint law; product draws never shuffle rows.
joint_abort <- function(message) stop(message, call. = FALSE)

joint_spec <- function(C, weights, domains, allocation_order, endpoints) {
  if (!is.data.frame(C) || nrow(C) < 1L) joint_abort("nonempty C data.frame required")
  if (!is.numeric(weights) || length(weights) != nrow(C) ||
      any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0)
    joint_abort("finite nonnegative weights with positive total required")
  if (!is.list(domains) || !length(domains) || is.null(names(domains)) ||
      anyNA(names(domains)) || any(!nzchar(names(domains))) || anyDuplicated(names(domains)) ||
      any(!vapply(domains, function(x) is.character(x) && length(x) > 0L &&
                    !anyNA(x) && all(nzchar(x)), logical(1))))
    joint_abort("uniquely named, nonempty domain column lists required")
  columns <- unlist(domains, use.names = FALSE)
  if (anyDuplicated(columns) || anyDuplicated(allocation_order) ||
      !setequal(allocation_order, names(domains))) joint_abort("invalid domains/allocation order")
  if (!is.character(endpoints) || length(endpoints) != 2L || anyNA(endpoints) ||
      any(!nzchar(endpoints)) || anyDuplicated(endpoints)) joint_abort("two distinct endpoints required")
  list(columns = columns, domains = domains, order = allocation_order, endpoints = endpoints)
}

joint_check_M <- function(M, C, columns, validate_support, context) {
  if (is.matrix(M)) M <- as.data.frame(M, stringsAsFactors = FALSE)
  # A named list of column vectors is also accepted; arbitrary gate/derived
  # coordinate names are legal as long as they partition the registered domains.
  if (is.list(M) && !is.data.frame(M) && !is.null(names(M)) &&
      all(vapply(M, length, integer(1)) == nrow(C)))
    M <- as.data.frame(M, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.data.frame(M) || nrow(M) != nrow(C) || anyDuplicated(names(M)) ||
      !setequal(names(M), columns) || anyNA(M) ||
      any(vapply(M, function(x) is.numeric(x) && any(!is.finite(x)), logical(1))))
    joint_abort(paste("invalid complete mediator vector:", context))
  M <- M[, columns, drop = FALSE]
  verdict <- validate_support(M = M, C = C, context = context)
  if (!isTRUE(verdict)) joint_abort(paste(context, paste(verdict, collapse = "; "),
    "support failure: no clipping, selective redraw, or row exchange permitted"))
  M
}

joint_predict <- function(b, M, C, predict_risk, endpoints, context) {
  risk <- predict_risk(b = b, M = M, C = C)
  if (is.data.frame(risk)) risk <- as.matrix(risk)
  if (!is.matrix(risk) || !is.numeric(risk) || nrow(risk) != nrow(C) ||
      is.null(colnames(risk)) || anyDuplicated(colnames(risk)) ||
      !setequal(colnames(risk), endpoints) || any(!is.finite(risk)) ||
      any(risk < 0 | risk > 1)) joint_abort(paste("invalid probability:", context))
  risk[, endpoints, drop = FALSE]
}

joint_scenarios <- function(order) c("R00", "R10", "R11", "P0", "P1", paste0("S:", order))

joint_finish <- function(risks, orders, endpoints, mc_risks, audit) {
  scenario_rows <- component_rows <- paired_rows <- list()
  for (ord in names(orders)) {
    rr <- risks[[ord]]
    for (e in endpoints) {
      scenario_rows[[length(scenario_rows) + 1L]] <- data.frame(
        allocation_order = ord, endpoint = e, scenario = rownames(rr),
        risk = rr[, e], risk_per_1000 = 1000 * rr[, e], row.names = NULL)
      chain <- c(rr["P0", e], rr[paste0("S:", orders[[ord]]), e])
      iie <- diff(chain)
      dep <- (rr["R11", e] - rr["P1", e]) - (rr["R10", e] - rr["P0", e])
      effect <- c(rr["R11", e] - rr["R00", e], rr["R10", e] - rr["R00", e],
                  rr["R11", e] - rr["R10", e], iie, dep)
      if (abs(sum(iie) + dep - effect[3L]) > 1e-12 ||
          abs(effect[1L] - effect[2L] - effect[3L]) > 1e-12)
        joint_abort("same-vector decomposition failed to close")
      component_rows[[length(component_rows) + 1L]] <- data.frame(
        allocation_order = ord, endpoint = e,
        component_type = c("TE", "DE", "IE_joint", rep("IIE", length(iie)), "DEP"),
        domain_id = c(NA_character_, NA_character_, NA_character_, orders[[ord]], NA_character_),
        effect = unname(effect), effect_per_1000 = 1000 * unname(effect), row.names = NULL)
    }
  }
  components <- do.call(rbind, component_rows)
  for (ord in names(orders)) {
    a <- components[components$allocation_order == ord & components$endpoint == endpoints[1L], ]
    b <- components[components$allocation_order == ord & components$endpoint == endpoints[2L], ]
    paired_rows[[length(paired_rows) + 1L]] <- data.frame(
      allocation_order = ord, component_type = a$component_type, domain_id = a$domain_id,
      contrast = paste(endpoints, collapse = " - "), Delta = a$effect - b$effect,
      Delta_per_1000 = 1000 * (a$effect - b$effect), row.names = NULL)
  }
  # MC standard errors describe integration error conditional on these fits,
  # not sampling uncertainty. Paired MC effects retain common-random covariance.
  mc_effects <- lapply(seq_along(mc_risks), function(k) {
    z <- joint_finish_no_mc(mc_risks[[k]], orders, endpoints)
    z$mc <- k
    z
  })
  mc_effects <- do.call(rbind, mc_effects)
  keys <- unique(mc_effects[, c("allocation_order", "endpoint", "component_type", "domain_id")])
  mc_se <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    same_domain <- if (is.na(key$domain_id)) is.na(mc_effects$domain_id) else
      !is.na(mc_effects$domain_id) & mc_effects$domain_id == key$domain_id
    v <- mc_effects$effect[mc_effects$allocation_order == key$allocation_order &
      mc_effects$endpoint == key$endpoint & mc_effects$component_type == key$component_type & same_domain]
    key$mc_se_probability <- if (length(v) > 1L) sd(v) / sqrt(length(v)) else NA_real_
    key$mc_se_per_1000 <- 1000 * key$mc_se_probability
    key
  })
  structure(list(scenarios = do.call(rbind, scenario_rows), components = components,
    paired = do.call(rbind, paired_rows), risks = risks, orders = orders,
    mc_risks = mc_risks, mc_effects = mc_effects, mc_se = do.call(rbind, mc_se), audit = audit),
    class = "joint_product_result")
}

joint_finish_no_mc <- function(risks, orders, endpoints) {
  rows <- list()
  for (ord in names(orders)) for (e in endpoints) {
    rr <- risks[[ord]]
    values <- c(rr["R11", e] - rr["R00", e], rr["R10", e] - rr["R00", e],
      rr["R11", e] - rr["R10", e], diff(c(rr["P0", e], rr[paste0("S:", orders[[ord]]), e])),
      (rr["R11", e] - rr["P1", e]) - (rr["R10", e] - rr["P0", e]))
    rows[[length(rows) + 1L]] <- data.frame(allocation_order = ord, endpoint = e,
      component_type = c("TE", "DE", "IE_joint", rep("IIE", length(orders[[ord]])), "DEP"),
      domain_id = c(NA_character_, NA_character_, NA_character_, orders[[ord]], NA_character_),
      effect = unname(values), stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# draw_joint(b,C,uniforms) must deterministically draw a complete joint vector,
# with ONE scalar exposure b=0/1 for its entire conditional chain. Bootstrap ID
# is metadata and must NEVER be passed as this b. Each domain's own independent
# L'Ecuyer stream drives a full pair of draws; only its columns are extracted.
# predict_risk is the cached, already loaded stable DOUBLE exponential CIF
# callback. This kernel neither sources/compiles C nor changes its integration.
# binding hashes identify scientific input/encoding/sample/model/RNG semantics,
# NOT RDS bytes, formula environments, cache paths, or coefficient last bits.
# save_checkpoint(state) persists aggregate-only state atomically/versionedly;
# invoke it only on the worker's private output tree. It is called once/MC pass.
joint_product_stream <- function(C, weights, domains, draw_joint, predict_risk, validate_support,
    uniform_dim, n_mc = 50L, block_size = 512L, seed = 1L,
    allocation_order = names(domains), include_reverse = TRUE, endpoints = c("IHD", "T2D"),
    binding, checkpoint = NULL, save_checkpoint = NULL, progress = NULL) {
  spec <- joint_spec(C, weights, domains, allocation_order, endpoints)
  if (!all(vapply(list(draw_joint, predict_risk, validate_support), is.function, logical(1))) ||
      (!is.null(save_checkpoint) && !is.function(save_checkpoint)) ||
      (!is.null(progress) && !is.function(progress))) joint_abort("invalid callback")
  for (key in c("uniform_dim", "n_mc", "block_size", "seed")) {
    value <- get(key)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value > .Machine$integer.max || value != floor(value))
      joint_abort(paste(key, "must be a positive integer"))
  }
  if (!is.logical(include_reverse) || length(include_reverse) != 1L || is.na(include_reverse))
    joint_abort("include_reverse must be scalar logical")
  required <- c("panel_id", "bootstrap_id", "sample_hash", "population_hash", "model_spec_hash", "input_hash", "rng_spec_hash")
  if (!is.list(binding) || !all(required %in% names(binding)) || anyDuplicated(names(binding)))
    joint_abort(paste("scientific binding requires", paste(required, collapse = ", ")))
  if (!is.numeric(binding$bootstrap_id) || length(binding$bootstrap_id) != 1L ||
      !is.finite(binding$bootstrap_id) || binding$bootstrap_id < 0 || binding$bootstrap_id > 100 ||
      binding$bootstrap_id != floor(binding$bootstrap_id)) joint_abort("fixed bootstrap ID 0..100 required")
  for (key in setdiff(required, "bootstrap_id"))
    if (!is.character(binding[[key]]) || length(binding[[key]]) != 1L || is.na(binding[[key]]) ||
        !nzchar(binding[[key]])) joint_abort(paste("invalid scientific binding", key))
  binding <- binding[sort(required)]
  orders <- list(forward = allocation_order)
  if (include_reverse) orders$reverse <- rev(allocation_order)
  zero <- lapply(orders, function(order) matrix(0, length(joint_scenarios(order)), length(endpoints),
    dimnames = list(joint_scenarios(order), endpoints)))
  identity <- list(schema = "NEW121_PRODUCT_MC_CHECKPOINT_V1", binding = binding,
    domains = domains, orders = orders, endpoints = endpoints, n_target = nrow(C),
    C_columns = names(C), uniform_dim = as.integer(uniform_dim), n_mc = as.integer(n_mc), seed = as.integer(seed),
    # Not an exact-float comparison or a cross-device reproducibility gate.
    weight_sum = sum(weights), stream_contract = "JOINT_PLUS_ONE_INDEPENDENT_STREAM_PER_DOMAIN_PAIRED_B0_B1")
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
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
  if (!is.null(checkpoint)) {
    old_identity <- checkpoint$identity
    if (!is.list(old_identity)) joint_abort("missing checkpoint identity")
    old_sum <- old_identity$weight_sum
    old_identity$weight_sum <- NULL
    compare_identity <- identity
    compare_identity$weight_sum <- NULL
    if (!identical(old_identity, compare_identity) ||
        !isTRUE(all.equal(old_sum, identity$weight_sum, tolerance = 1e-10)))
      joint_abort("checkpoint scientific/sample/stream identity mismatch")
    completed <- checkpoint$completed_mc
    mc_risks <- checkpoint$mc_risks
    states <- checkpoint$stream_states
    if (!is.numeric(completed) || length(completed) != 1L || !is.finite(completed) ||
        completed < 0 || completed > n_mc || completed != floor(completed) ||
        length(mc_risks) != completed || !identical(names(states), stream_names) ||
        any(!vapply(states, function(s) is.integer(s) && length(s) == 7L && !anyNA(s), logical(1))))
      joint_abort("corrupt checkpoint: completed-pass accounting/streams")
    for (pass in mc_risks) for (ord in names(orders)) {
      x <- pass[[ord]]
      if (!is.matrix(x) || !identical(dimnames(x), dimnames(zero[[ord]])) ||
          any(!is.finite(x)) || any(x < -1e-12 | x > 1+1e-12)) joint_abort("corrupt checkpoint risk means")
    }
    prediction_calls <- checkpoint$prediction_calls
    if (length(prediction_calls) != 1L || !is.finite(prediction_calls) || prediction_calls < 0)
      joint_abort("corrupt prediction count")
    resumed <- completed > 0L
  }
  next_uniforms <- function(stream, n) {
    assign(".Random.seed", states[[stream]], envir = .GlobalEnv)
    # Row-major consumption keeps the participant/variate mapping independent
    # of block boundaries, including resumption on a different worker.
    z <- matrix(runif(n * uniform_dim), nrow = n, byrow = TRUE)
    states[[stream]] <<- get(".Random.seed", envir = .GlobalEnv)
    z
  }
  deterministic <- function(callback) {
    before <- get(".Random.seed", envir = .GlobalEnv)
    value <- callback()
    if (!identical(before, get(".Random.seed", envir = .GlobalEnv)))
      joint_abort("simulation/prediction/support callbacks must use supplied uniforms only")
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
      check <- function(M, label) deterministic(function()
        joint_check_M(M, cc, spec$columns, validate_support, paste(context, label, sep = ";")))
      draw <- function(b, u, label) check(deterministic(function()
        draw_joint(b = b, C = cc, uniforms = u)), label)
      predict <- function(b, M, label) {
        prediction_calls <<- prediction_calls + 1L
        colSums(deterministic(function() joint_predict(b, M, cc, predict_risk,
          endpoints, paste(context, label, sep = ";"))) * ww)
      }
      u <- next_uniforms("JOINT", length(ii))
      Q0 <- draw(0, u, "jointQ0")
      Q1 <- draw(1, u, "jointQ1")
      overall <- rbind(R00 = predict(0, Q0, "R00"), R10 = predict(1, Q0, "R10"),
        R11 = predict(1, Q1, "R11"))
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
        sums[[ord]][c("R00", "R10", "R11"), ] <- sums[[ord]][c("R00", "R10", "R11"), ] + overall
        sums[[ord]][c("P0", "P1"), ] <- sums[[ord]][c("P0", "P1"), ] + endpoints_product
        M <- P0
        for (j in seq_along(orders[[ord]])) {
          g <- orders[[ord]][j]
          cols <- domains[[g]]
          M[, cols] <- P1[, cols, drop = FALSE]
          scenario <- paste0("S:", g)
          value <- if (j == length(orders[[ord]])) endpoints_product["P1", ] else
            predict(1, check(M, paste(ord, scenario)), paste(ord, scenario))
          sums[[ord]][scenario, ] <- sums[[ord]][scenario, ] + value
        }
      }
      if (!is.null(progress)) progress(list(stage = "integrating", bootstrap_id = binding$bootstrap_id,
        mc = k, n_mc = n_mc, rows_complete = max(ii), n_target = nrow(C),
        prediction_calls = prediction_calls, scenarios_per_block = 5L +
          sum(vapply(orders, length, integer(1)) - 1L)))
    }
    # Commit whole MC passes only. A partial pass after interruption is replayed
    # with its original saved streams; no partially accumulated values are reused.
    mc_risks[[k]] <- lapply(sums, function(x) x / sum(weights))
    completed <- k
    state <- list(identity = identity, completed_mc = completed, mc_risks = mc_risks,
      stream_states = states, prediction_calls = prediction_calls)
    if (!is.null(save_checkpoint)) save_checkpoint(state)
  }
  risks <- lapply(names(orders), function(ord) Reduce(`+`, lapply(mc_risks, `[[`, ord)) / n_mc)
  names(risks) <- names(orders)
  joint_finish(risks, orders, endpoints, mc_risks, list(
    schema = "NEW121_PRODUCT_RESULT_V1", method = "conditional independent domain marginals; cumulative replacement",
    binding = binding, n_target = nrow(C), weight_sum = sum(weights), n_mc = completed,
    block_size = block_size, uniform_dim = uniform_dim, seed = seed,
    initial_stream_states = initial_states, final_stream_states = states,
    stream_contract = identity$stream_contract, same_C_within_all_domain_draws = TRUE,
    independent_domain_streams = TRUE, paired_uniforms_at_exposure_0_1 = TRUE,
    R_scenarios_full_joint = TRUE, shared_product_endpoints_both_orders = TRUE,
    resumed = resumed, completed_mc_on_entry = if (resumed) checkpoint$completed_mc else 0L,
    prediction_calls = prediction_calls, model_fits = 0L,
    elapsed_seconds = proc.time()[["elapsed"]] - started, units = "probability and explicit per_1000 columns",
    cross_device_probability_reference = 1e-6, bitwise_cross_device_agreement_required = FALSE))
}

