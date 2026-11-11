# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

s44pr_a_fail <- function(code, detail = "", rebuild_eligible = FALSE) {
    stop(structure(list(message = paste(code, detail), call = NULL, code = code, rebuild_eligible = rebuild_eligible), 
        class = c("s44pr_a_cache_error", "error", "condition")))
}

s44pr_a_assert <- function(ok, code) if (!isTRUE(ok)) s44pr_a_fail(code)

s44pr_a_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)

s44pr_a_vec <- function(x) as.character(unlist(x, use.names = FALSE))

s44pr_a_sha <- function(path) digest::digest(file = path, algo = "sha256")

s44pr_A_source_pins <- function() c(`R/body_group_pmm_A_v2.R` = birth_private_setting("PRIVATE_HASH_1"), `R/body_tqe_A.R` = birth_private_setting("PRIVATE_HASH_2"), 
    `vendor/core39/R/core39_models.R` = birth_private_setting("PRIVATE_HASH_3"), `vendor/core39/repairs/bmi_log_fat_linear_v1/R/mixed_a_revision.R` = birth_private_setting("PRIVATE_HASH_4"), 
    `vendor/border_fine_domain_behavior_joint_20261031_v1/f/00.core.R` = birth_private_setting("PRIVATE_HASH_5"), 
    `vendor/border_fine_domain_behavior_joint_20261031_v1/f/20.joint.R` = birth_private_setting("PRIVATE_HASH_6"), 
    `R/selected_nmr_composition_A.R` = birth_private_setting("PRIVATE_HASH_7"), `R/selected_nmr_blockspec.R` = birth_private_setting("PRIVATE_HASH_8"), 
    `vendor/NMR_BLOCKS_RAW54.source.json` = birth_private_setting("PRIVATE_HASH_9"), `R/nmr23591_logit_A_v2.R` = birth_private_setting("PRIVATE_HASH_10"))

s44pr_A_context <- function(z, ir, ws, ticket, input_sha256, index_sha256, weights_sha256, expected_identity = NULL) {
    s44pr_a_assert(requireNamespace("digest", quietly = TRUE), "DIGEST_REQUIRED_NO_INSTALL")
    i <- ticket$identity
    draw <- ticket$draw
    s44pr_a_assert(i$model_id %in% c("S44_PROT_REF11", "S44_PROT_PLUS15") && i$line == "PROTEIN" && as.integer(i$b) %in% 
        0:100 && as.integer(i$k) == 50L, "ORIGINAL_TWO_PROTEIN_B0_100_SCOPE_REQUIRED")
    if (!is.null(expected_identity)) {
        for (key in c("model_id", "b", "line", "input_hash", "population_hash", "sampling_hash", "draw_token", 
            "seed", "uniform_dim", "rng_spec_hash", "model_spec")) s44pr_a_assert(identical(expected_identity[[key]], 
            i[[key]]), paste0("A_IDENTITY_CHANGED:", key))
    }
    s44pr_a_assert(identical(input_sha256, i$input_hash) && identical(index_sha256, draw$index_sha256) && identical(weights_sha256, 
        draw$weights_sha256), "ONCE_VERIFIED_ORIGINAL_INPUT_INDEX_WEIGHT_SHA_REQUIRED")
    s44pr_a_assert(identical(z$schema, "SELECTED44_LINE_FROZEN_INPUTS_V1") && z$line == "PROTEIN" && identical(s44pr_a_hash(as.character(z$data$eid)), 
        z$target_hash) && identical(z$target_hash, i$population_hash), "ORIGINAL_TARGET_AND_INPUT_BINDING_CHANGED")
    target <- if (is.null(ir$targethash)) 
        ir$target_hash
    else ir$targethash
    s44pr_a_assert(identical(target, z$target_hash) && as.integer(ir$bootstrap_id) == as.integer(i$b) && identical(ir$cluster_membership_hash, 
        z$cluster_membership_hash), "ORIGINAL_BOOTSTRAP_FAMILIES_CHANGED")
    ix <- ir$indices
    s44pr_a_assert(is.numeric(ix) && !anyNA(ix) && all(ix == floor(ix) & ix >= 1 & ix <= nrow(z$data)) && is.logical(ws$complete) && 
        length(ws$complete) == length(ix) && !anyNA(ws$complete), "ORIGINAL_INDEX_OR_COMPLETENESS_MASK_INVALID")
    d <- z$data[ix[ws$complete], , drop = FALSE]
    w <- ws$w
    pb <- ws$pair_binding
    s44pr_a_assert(length(w) == nrow(d) && all(is.finite(w) & w > 0) && "cluster_id" %in% names(d) && !anyNA(d$cluster_id), 
        "ORIGINAL_A_ROWS_WEIGHTS_FAMILIES_INVALID")
    s44pr_a_assert(identical(s44pr_a_hash(as.integer(ix)), i$sampling_hash) && identical(draw$sampling_hash, i$sampling_hash) && 
        as.integer(draw$b) == as.integer(i$b) && as.integer(draw$seed) == as.integer(i$seed), "ORIGINAL_DRAW_AND_RNG_SEED_CHANGED")
    s44pr_a_assert(is.list(pb) && identical(pb$target_hash, z$target_hash) && identical(pb$sample_hash, s44pr_a_hash(as.character(d$eid))) && 
        identical(pb$design_weight_hash, s44pr_a_hash(ws$design_weights)) && identical(pb$selection_weight_hash, 
        s44pr_a_hash(ws$selection_weights)) && identical(pb$sampling_hash, i$sampling_hash) && identical(pb$actual_index_sha256, 
        index_sha256) && identical(as.numeric(w), as.numeric(ws$design_weights * ws$selection_weights)), "ORIGINAL_COMPLETE_MASK_AND_IPW_COMPONENT_BINDING_CHANGED")
    s44pr_a_assert(identical(draw$sample_hash, pb$sample_hash) && identical(draw$selection_weight_hash, pb$selection_weight_hash) && 
        as.integer(draw$N_fit) == nrow(d) && as.integer(draw$N_resampled) == length(ix), "ORIGINAL_DRAW_POPULATION_COUNTS_OR_SELECTION_WEIGHTS_CHANGED")
    s44pr_a_assert(setequal(names(z$risks), c("IHD", "T2D")), "ORIGINAL_TWO_ENDPOINT_RISK_MAPPING_REQUIRED")
    rows <- ix[ws$complete]
    risks <- lapply(z$risks, function(r) {
        s44pr_a_assert(nrow(r) == nrow(z$data) && all(c("event", "death_event", "exit") %in% names(r)), "ORIGINAL_FULLFOLLOWUP_RISK_ROW_COUNT_CHANGED")
        value <- r[rows, , drop = FALSE]
        s44pr_a_assert(!anyNA(value[, c("event", "death_event", "exit")]) && all(value$event %in% 0:1) && all(value$death_event %in% 
            0:1) && all(is.finite(value$exit) & value$exit >= 0), "ORIGINAL_RISK_FLAGS_OR_FOLLOWUP_INVALID")
        value
    })
    ids <- s44pr_a_vec(if (is.null(i$model_spec$mediator_ids)) 
        i$model_spec$members
    else i$model_spec$mediator_ids)
    expected_ids <- if (i$model_id == "S44_PROT_REF11") 
        s44pr_a_vec(z$core_ids)
    else names(z$members)
    s44pr_a_assert(identical(ids, expected_ids) && length(ids) == if (i$model_id == "S44_PROT_REF11") 
        11L
    else 15L, "ORIGINAL_MODEL_MEMBERS_AND_ORDER_CHANGED")
    members <- z$members[ids]
    cols <- vapply(members, `[[`, "", "column")
    domains <- lapply(z$domains, function(x) intersect(s44pr_a_vec(x), cols))
    domains <- Filter(length, domains)
    order <- s44pr_a_vec(i$model_spec$allocation_order)
    s44pr_a_assert(setequal(names(domains), order), "ORIGINAL_DOMAIN_ALLOCATION_CHANGED")
    domains <- domains[order]
    base <- s44pr_a_vec(z$base)
    s44pr_a_assert(identical(base, s44pr_a_vec(i$model_spec$base)) && as.integer(i$uniform_dim) == 15L, "ORIGINAL_BASE_OR_PROTEIN15_COMMON_RNG_WIDTH_CHANGED")
    list(z = z, ir = ir, ws = ws, ticket = ticket, identity = i, draw = draw, d = d, w = w, risks = risks, pair_binding = pb, 
        model_ids = ids, members = members, domains = domains, base = base, family_hash = s44pr_a_hash(as.character(d$cluster_id)), 
        risk_row_mapping_hash = s44pr_a_hash(list(rows = as.integer(rows), eid = as.character(d$eid), risks = lapply(risks, 
            function(r) r[, c("event", "death_event", "exit"), drop = FALSE]))), rng = list(kind = "L'Ecuyer-CMRG", 
            normal = "Inversion", sample = "Rejection", seed = as.integer(i$seed), uniform_dim = 15L, independent_domain_streams = TRUE, 
            rng_spec_hash = i$rng_spec_hash))
}

s44pr_A_original_key <- function(context, source_root, methods) {
    i <- context$identity
    s44pr_a_assert(is.function(methods$s44_nmr_blockspec), "FROZEN_BLOCKSPEC_FUNCTION_REQUIRED_FOR_ORIGINAL_A_KEY")
    spec <- methods$s44_nmr_blockspec(jsonlite::fromJSON(file.path(source_root, "vendor/NMR_BLOCKS_RAW54.source.json"), 
        simplifyVector = FALSE))
    composition <- s44pr_a_hash(list(spec = spec, A = s44pr_a_sha(file.path(source_root, "R/selected_nmr_composition_A.R")), 
        block_adapter = s44pr_a_sha(file.path(source_root, "R/selected_nmr_blockspec.R"))))
    s44pr_a_hash(list(input = i$input_hash, population = i$population_hash, sample = i$sampling_hash, actual_rows = s44pr_a_hash(as.character(context$d$eid)), 
        draw = i$draw_token, weights = context$draw$weights_sha256, source = i$source_hash, model = i$model_spec, 
        members = context$members, domains = context$domains, composition = composition, selected44_A_revision = "BODY_PMM23591_V2", 
        statistical_version = i$statistical_version, body_A_source = s44pr_a_sha(file.path(source_root, "R/body_group_pmm_A_v2.R")), 
        nmr23591_source = s44pr_a_sha(file.path(source_root, "R/nmr23591_logit_A_v2.R")), seed = i$seed))
}

s44pr_reuse_A <- function(cache_rds, cache_record, context, engine, methods, source_root, transfer_receipt = NULL, 
    probe_rows = 3L) {
    i <- context$identity
    am <- cache_record$A_manifest
    s44pr_a_assert(requireNamespace("jsonlite", quietly = TRUE) && requireNamespace("digest", quietly = TRUE), 
        "EXISTING_JSONLITE_DIGEST_REQUIRED_NO_INSTALL")
    s44pr_a_assert(isTRUE(cache_record$metadata_usable) && identical(cache_record$ticket$identity, i) && identical(cache_record$ticket$draw, 
        context$draw), "ORIGINAL_SERVER_METADATA_BINDING_CHANGED")
    for (key in names(context$pair_binding)) s44pr_a_assert(identical(context$pair_binding[[key]], cache_record$original_pair_binding[[key]]), 
        paste0("ORIGINAL_SERVER_PAIR_BINDING_CHANGED:", key))
    pins <- s44pr_A_source_pins()
    for (rel in names(pins)) s44pr_a_assert(identical(s44pr_a_sha(file.path(source_root, rel)), unname(pins[[rel]])), 
        paste0("ORIGINAL_A_SOURCE_CHANGED:", rel))
    s44pr_a_assert(isTRUE(am$complete) && identical(am$model_id, i$model_id) && identical(am$source_hash, i$source_hash) && 
        identical(s44pr_a_vec(am$members), context$model_ids) && identical(am$model_A_key, s44pr_A_original_key(context, 
        source_root, methods)) && identical(am$file$sha256, cache_record$A_sha256), "ORIGINAL_MODEL_A_KEY_OR_MANIFEST_CHANGED")
    if (!file.exists(cache_rds)) 
        s44pr_a_fail("ORIGINAL_A_FILE_MISSING", rebuild_eligible = TRUE)
    info <- file.info(cache_rds)
    s44pr_a_assert(!is.na(info$size) && info$size == cache_record$A_bytes, "ORIGINAL_A_FILE_SIZE_CHANGED")
    if (is.null(transfer_receipt)) {
        s44pr_a_assert(identical(s44pr_a_sha(cache_rds), am$file$sha256), "ORIGINAL_A_FILE_SHA_CHANGED")
    }
    else {
        s44pr_a_assert(isTRUE(transfer_receipt$verified) && identical(transfer_receipt$sha256, am$file$sha256) && 
            transfer_receipt$bytes == info$size && identical(normalizePath(transfer_receipt$path, mustWork = TRUE), 
            normalizePath(cache_rds, mustWork = TRUE)) && !is.null(transfer_receipt$mtime_unix) && abs(as.numeric(info$mtime) - 
            as.numeric(transfer_receipt$mtime_unix)) < 1e-06, "ONE_TIME_STREAM_HASH_AND_UNCHANGED_LOCAL_FILE_RECEIPT_REQUIRED")
    }
    A <- tryCatch(readRDS(cache_rds), error = function(e) s44pr_a_fail("ORIGINAL_A_RDS_UNREADABLE", conditionMessage(e), 
        rebuild_eligible = TRUE))
    s44pr_a_assert(is.list(A) && identical(names(A$models), context$model_ids) && identical(A$members, context$members) && 
        identical(A$conditional_order, context$model_ids), "ORIGINAL_A_MODEL_MEMBER_ORDER_OR_CONDITIONAL_LAWS_CHANGED")
    s44pr_a_assert(identical(A$body_pmm$revision, engine$body_pmm_revision()) && identical(A$body_pmm$prefix_binding, 
        engine$body_tqe_binding(context$d, context$members, context$base, context$w)) && identical(A$models$bmi$log_a_revision$revision, 
        "CORE39_BMI_LOG_FAT_LINEAR_V1") && isTRUE(A$body_pmm$raw_members_unchanged) && isTRUE(A$body_pmm$no_height) && 
        as.integer(A$body_pmm$group_uniform_slot) == 3L && identical(as.integer(A$body_pmm$unused_fixed_uniform_slots), 
        4:7), "ORIGINAL_BODY_PMM_BMI_AND_RNG_SLOT_BINDING_CHANGED")
    d <- context$d
    members <- context$members
    base <- context$base
    w <- context$w
    columns <- vapply(members, `[[`, "", "column")
    bodycols <- engine$body_tqe_columns(members)
    for (j in seq_along(A$models)) {
        f <- A$models[[j]]
        before <- columns[seq_len(j - 1L)]
        s44pr_a_assert(identical(f$member, members[[j]]), "ORIGINAL_A_MEMBER_DEFINITION_CHANGED")
        if (isTRUE(f$body_pmm_group_slot)) {
            s44pr_a_assert(as.integer(f$original_uniform_slot) == j && isTRUE(f$uniform_used) == identical(names(members)[j], 
                "weight"), "ORIGINAL_BODY_DONOR_SLOT_CHANGED")
        }
        else {
            prepared <- engine$body_tqe_affine_prepare(d, before, bodycols, reference = f$body_tqe_affine, w = w)
            expected <- engine$c39_a_design(prepared$data, base, prepared$previous)
            actual <- engine$c39_a_design(prepared$data, NULL, NULL, f$design_reference)
            s44pr_a_assert(identical(unname(f$predecessor_columns), unname(before)) && identical(engine$body_tqe_reference_copy(f$design_reference), 
                engine$body_tqe_reference_copy(expected$reference)) && identical(actual$X, expected$X), "ORIGINAL_A_FORMULA_FACTOR_LEVELS_OR_DESIGN_CHANGED")
        }
    }
    group <- A$body_pmm$group
    enc <- engine$body_tqe_encode(d, members)
    gc <- group$columns
    tuple <- cbind(T = enc$total, q = enc$fraction, e = enc$z$e, hip = d[[gc[["hip"]]]], waist = d[[gc[["waist"]]]])
    s44pr_a_assert(identical(group$tuple, tuple) && as.integer(group$donor_N) == nrow(d) && identical(group$spec, 
        engine$body_pmm_spec()), "ORIGINAL_BODY_DONOR_VALUES_OR_METHOD_SPEC_CHANGED")
    upstream <- vapply(members[c("smoking", "bmi")], `[[`, "", "column")
    gx <- engine$c39_a_design(d, base, upstream)
    ga <- engine$c39_a_design(d, NULL, NULL, group$design_reference)
    s44pr_a_assert(identical(gx$X, ga$X) && identical(engine$body_tqe_reference_copy(gx$reference), engine$body_tqe_reference_copy(group$design_reference)), 
        "ORIGINAL_BODY_MATCHING_DESIGN_CHANGED")
    coverage <- integer()
    for (leaf in group$tree$leaves) {
        rows <- as.integer(leaf$rows)
        coverage <- c(coverage, rows)
        s44pr_a_assert(length(rows) == as.integer(leaf$N) && all(rows >= 1L & rows <= nrow(d)) && length(leaf$cdf) == 
            length(rows) + 1L && leaf$cdf[1L] == 0 && tail(leaf$cdf, 1L) == 1 && all(is.finite(leaf$cdf)) && all(diff(leaf$cdf) >= 
            0), "ORIGINAL_BODY_DONOR_CDF_INVALID")
        expected_cdf <- c(0, cumsum(w[rows]/sum(w[rows])))
        expected_cdf[length(expected_cdf)] <- 1
        s44pr_a_assert(max(abs(expected_cdf - leaf$cdf)) <= 1e-12, "ORIGINAL_BODY_DONOR_WEIGHT_PROBABILITIES_CHANGED")
    }
    s44pr_a_assert(identical(sort(coverage), seq_len(nrow(d))), "ORIGINAL_BODY_DONOR_PARTITION_CHANGED")
    kind <- RNGkind()
    had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    if (had_seed) 
        seed <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
    on.exit({
        do.call(RNGkind, as.list(kind))
        if (had_seed) assign(".Random.seed", seed, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", 
            envir = .GlobalEnv)
    }, add = TRUE)
    n <- min(as.integer(probe_rows), nrow(d))
    s44pr_a_assert(n >= 1L && n <= 8L, "BOUNDED_A_PORTABILITY_PROBE_REQUIRED")
    C <- d[seq_len(n), , drop = FALSE]
    u <- matrix(0.5, n, length(A$models))
    probe <- tryCatch(lapply(0:1, function(b) {
        value <- engine$c39_draw_a(A, C, b, u)
        engine$c39_support(value, members)
        s44pr_a_assert(nrow(value) == n && all(columns %in% names(value)), "A_PROBE_SHAPE_CHANGED")
        TRUE
    }), error = function(e) s44pr_a_fail("ORIGINAL_A_PORTABILITY_PROBE_FAILED", conditionMessage(e), rebuild_eligible = TRUE))
    list(A = A, provenance = list(schema = "PROTEIN_RIDGE_ORIGINAL_A_REUSE_CERTIFICATE_V1", passed = TRUE, model_id = i$model_id, 
        b = as.integer(i$b), original_model_A_key = am$model_A_key, A_sha256 = am$file$sha256, original_ticket_sha256 = cache_record$files$ticket$sha256, 
        original_A_manifest_sha256 = cache_record$files$A_manifest$sha256, original_ticket_semantic_hash = s44pr_a_hash(cache_record$ticket), 
        original_A_manifest_semantic_hash = s44pr_a_hash(am), input_hash = i$input_hash, index_sha256 = context$draw$index_sha256, 
        weights_sha256 = context$draw$weights_sha256, sample_hash = context$pair_binding$sample_hash, family_hash = context$family_hash, 
        risk_row_mapping_hash = context$risk_row_mapping_hash, pair_binding = context$pair_binding, members = context$model_ids, 
        rng = context$rng, original_R = cache_record$original_R, current_R = R.version.string, source_pins = as.list(pins), 
        A_refitted = FALSE, A_object_modified = FALSE, Y_hash_not_an_A_key = TRUE, portability_probe_rows = n, 
        portability_probe_exposure_values = 0:1, official_MC_integrations = 0L, file_sha_verified_in_stream = !is.null(transfer_receipt)))
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        s44pr_a_fail <- function(code, detail = "", rebuild_eligible = FALSE) {
            stop(structure(list(message = paste(code, detail), call = NULL, code = code, rebuild_eligible = rebuild_eligible), 
                class = c("s44pr_a_cache_error", "error", "condition")))
        }
        s44pr_a_assert <- function(ok, code) if (!isTRUE(ok)) s44pr_a_fail(code)
        s44pr_a_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)
        s44pr_a_vec <- function(x) as.character(unlist(x, use.names = FALSE))
        s44pr_a_sha <- function(path) digest::digest(file = path, algo = "sha256")
        s44pr_A_source_pins <- function() c(`R/body_group_pmm_A_v2.R` = birth_private_setting("PRIVATE_HASH_11"), 
            `R/body_tqe_A.R` = birth_private_setting("PRIVATE_HASH_12"), `vendor/core39/R/core39_models.R` = birth_private_setting("PRIVATE_HASH_13"), 
            `vendor/core39/repairs/bmi_log_fat_linear_v1/R/mixed_a_revision.R` = birth_private_setting("PRIVATE_HASH_14"), 
            `vendor/border_fine_domain_behavior_joint_20261031_v1/f/00.core.R` = birth_private_setting("PRIVATE_HASH_15"), 
            `vendor/border_fine_domain_behavior_joint_20261031_v1/f/20.joint.R` = birth_private_setting("PRIVATE_HASH_16"), 
            `R/selected_nmr_composition_A.R` = birth_private_setting("PRIVATE_HASH_17"), `R/selected_nmr_blockspec.R` = birth_private_setting("PRIVATE_HASH_18"), 
            `vendor/NMR_BLOCKS_RAW54.source.json` = birth_private_setting("PRIVATE_HASH_19"), `R/nmr23591_logit_A_v2.R` = birth_private_setting("PRIVATE_HASH_20"))
        s44pr_A_context <- function(z, ir, ws, ticket, input_sha256, index_sha256, weights_sha256, expected_identity = NULL) {
            s44pr_a_assert(requireNamespace("digest", quietly = TRUE), "DIGEST_REQUIRED_NO_INSTALL")
            i <- ticket$identity
            draw <- ticket$draw
            s44pr_a_assert(i$model_id %in% c("S44_PROT_REF11", "S44_PROT_PLUS15") && i$line == "PROTEIN" && as.integer(i$b) %in% 
                0:100 && as.integer(i$k) == 50L, "ORIGINAL_TWO_PROTEIN_B0_100_SCOPE_REQUIRED")
            if (!is.null(expected_identity)) {
                for (key in c("model_id", "b", "line", "input_hash", "population_hash", "sampling_hash", "draw_token", 
                  "seed", "uniform_dim", "rng_spec_hash", "model_spec")) s44pr_a_assert(identical(expected_identity[[key]], 
                  i[[key]]), paste0("A_IDENTITY_CHANGED:", key))
            }
            s44pr_a_assert(identical(input_sha256, i$input_hash) && identical(index_sha256, draw$index_sha256) && 
                identical(weights_sha256, draw$weights_sha256), "ONCE_VERIFIED_ORIGINAL_INPUT_INDEX_WEIGHT_SHA_REQUIRED")
            s44pr_a_assert(identical(z$schema, "SELECTED44_LINE_FROZEN_INPUTS_V1") && z$line == "PROTEIN" && identical(s44pr_a_hash(as.character(z$data$eid)), 
                z$target_hash) && identical(z$target_hash, i$population_hash), "ORIGINAL_TARGET_AND_INPUT_BINDING_CHANGED")
            target <- if (is.null(ir$targethash)) ir$target_hash else ir$targethash
            s44pr_a_assert(identical(target, z$target_hash) && as.integer(ir$bootstrap_id) == as.integer(i$b) && 
                identical(ir$cluster_membership_hash, z$cluster_membership_hash), "ORIGINAL_BOOTSTRAP_FAMILIES_CHANGED")
            ix <- ir$indices
            s44pr_a_assert(is.numeric(ix) && !anyNA(ix) && all(ix == floor(ix) & ix >= 1 & ix <= nrow(z$data)) && 
                is.logical(ws$complete) && length(ws$complete) == length(ix) && !anyNA(ws$complete), "ORIGINAL_INDEX_OR_COMPLETENESS_MASK_INVALID")
            d <- z$data[ix[ws$complete], , drop = FALSE]
            w <- ws$w
            pb <- ws$pair_binding
            s44pr_a_assert(length(w) == nrow(d) && all(is.finite(w) & w > 0) && "cluster_id" %in% names(d) && !anyNA(d$cluster_id), 
                "ORIGINAL_A_ROWS_WEIGHTS_FAMILIES_INVALID")
            s44pr_a_assert(identical(s44pr_a_hash(as.integer(ix)), i$sampling_hash) && identical(draw$sampling_hash, 
                i$sampling_hash) && as.integer(draw$b) == as.integer(i$b) && as.integer(draw$seed) == as.integer(i$seed), 
                "ORIGINAL_DRAW_AND_RNG_SEED_CHANGED")
            s44pr_a_assert(is.list(pb) && identical(pb$target_hash, z$target_hash) && identical(pb$sample_hash, 
                s44pr_a_hash(as.character(d$eid))) && identical(pb$design_weight_hash, s44pr_a_hash(ws$design_weights)) && 
                identical(pb$selection_weight_hash, s44pr_a_hash(ws$selection_weights)) && identical(pb$sampling_hash, 
                i$sampling_hash) && identical(pb$actual_index_sha256, index_sha256) && identical(as.numeric(w), 
                as.numeric(ws$design_weights * ws$selection_weights)), "ORIGINAL_COMPLETE_MASK_AND_IPW_COMPONENT_BINDING_CHANGED")
            s44pr_a_assert(identical(draw$sample_hash, pb$sample_hash) && identical(draw$selection_weight_hash, 
                pb$selection_weight_hash) && as.integer(draw$N_fit) == nrow(d) && as.integer(draw$N_resampled) == 
                length(ix), "ORIGINAL_DRAW_POPULATION_COUNTS_OR_SELECTION_WEIGHTS_CHANGED")
            s44pr_a_assert(setequal(names(z$risks), c("IHD", "T2D")), "ORIGINAL_TWO_ENDPOINT_RISK_MAPPING_REQUIRED")
            rows <- ix[ws$complete]
            risks <- lapply(z$risks, function(r) {
                s44pr_a_assert(nrow(r) == nrow(z$data) && all(c("event", "death_event", "exit") %in% names(r)), 
                  "ORIGINAL_FULLFOLLOWUP_RISK_ROW_COUNT_CHANGED")
                value <- r[rows, , drop = FALSE]
                s44pr_a_assert(!anyNA(value[, c("event", "death_event", "exit")]) && all(value$event %in% 0:1) && 
                  all(value$death_event %in% 0:1) && all(is.finite(value$exit) & value$exit >= 0), "ORIGINAL_RISK_FLAGS_OR_FOLLOWUP_INVALID")
                value
            })
            ids <- s44pr_a_vec(if (is.null(i$model_spec$mediator_ids)) i$model_spec$members else i$model_spec$mediator_ids)
            expected_ids <- if (i$model_id == "S44_PROT_REF11") s44pr_a_vec(z$core_ids) else names(z$members)
            s44pr_a_assert(identical(ids, expected_ids) && length(ids) == if (i$model_id == "S44_PROT_REF11") 11L else 15L, 
                "ORIGINAL_MODEL_MEMBERS_AND_ORDER_CHANGED")
            members <- z$members[ids]
            cols <- vapply(members, `[[`, "", "column")
            domains <- lapply(z$domains, function(x) intersect(s44pr_a_vec(x), cols))
            domains <- Filter(length, domains)
            order <- s44pr_a_vec(i$model_spec$allocation_order)
            s44pr_a_assert(setequal(names(domains), order), "ORIGINAL_DOMAIN_ALLOCATION_CHANGED")
            domains <- domains[order]
            base <- s44pr_a_vec(z$base)
            s44pr_a_assert(identical(base, s44pr_a_vec(i$model_spec$base)) && as.integer(i$uniform_dim) == 15L, 
                "ORIGINAL_BASE_OR_PROTEIN15_COMMON_RNG_WIDTH_CHANGED")
            list(z = z, ir = ir, ws = ws, ticket = ticket, identity = i, draw = draw, d = d, w = w, risks = risks, 
                pair_binding = pb, model_ids = ids, members = members, domains = domains, base = base, family_hash = s44pr_a_hash(as.character(d$cluster_id)), 
                risk_row_mapping_hash = s44pr_a_hash(list(rows = as.integer(rows), eid = as.character(d$eid), risks = lapply(risks, 
                  function(r) r[, c("event", "death_event", "exit"), drop = FALSE]))), rng = list(kind = "L'Ecuyer-CMRG", 
                  normal = "Inversion", sample = "Rejection", seed = as.integer(i$seed), uniform_dim = 15L, independent_domain_streams = TRUE, 
                  rng_spec_hash = i$rng_spec_hash))
        }
        s44pr_A_original_key <- function(context, source_root, methods) {
            i <- context$identity
            s44pr_a_assert(is.function(methods$s44_nmr_blockspec), "FROZEN_BLOCKSPEC_FUNCTION_REQUIRED_FOR_ORIGINAL_A_KEY")
            spec <- methods$s44_nmr_blockspec(jsonlite::fromJSON(file.path(source_root, "vendor/NMR_BLOCKS_RAW54.source.json"), 
                simplifyVector = FALSE))
            composition <- s44pr_a_hash(list(spec = spec, A = s44pr_a_sha(file.path(source_root, "R/selected_nmr_composition_A.R")), 
                block_adapter = s44pr_a_sha(file.path(source_root, "R/selected_nmr_blockspec.R"))))
            s44pr_a_hash(list(input = i$input_hash, population = i$population_hash, sample = i$sampling_hash, actual_rows = s44pr_a_hash(as.character(context$d$eid)), 
                draw = i$draw_token, weights = context$draw$weights_sha256, source = i$source_hash, model = i$model_spec, 
                members = context$members, domains = context$domains, composition = composition, selected44_A_revision = "BODY_PMM23591_V2", 
                statistical_version = i$statistical_version, body_A_source = s44pr_a_sha(file.path(source_root, 
                  "R/body_group_pmm_A_v2.R")), nmr23591_source = s44pr_a_sha(file.path(source_root, "R/nmr23591_logit_A_v2.R")), 
                seed = i$seed))
        }
        s44pr_reuse_A <- function(cache_rds, cache_record, context, engine, methods, source_root, transfer_receipt = NULL, 
            probe_rows = 3L) {
            i <- context$identity
            am <- cache_record$A_manifest
            s44pr_a_assert(requireNamespace("jsonlite", quietly = TRUE) && requireNamespace("digest", quietly = TRUE), 
                "EXISTING_JSONLITE_DIGEST_REQUIRED_NO_INSTALL")
            s44pr_a_assert(isTRUE(cache_record$metadata_usable) && identical(cache_record$ticket$identity, i) && 
                identical(cache_record$ticket$draw, context$draw), "ORIGINAL_SERVER_METADATA_BINDING_CHANGED")
            for (key in names(context$pair_binding)) s44pr_a_assert(identical(context$pair_binding[[key]], cache_record$original_pair_binding[[key]]), 
                paste0("ORIGINAL_SERVER_PAIR_BINDING_CHANGED:", key))
            pins <- s44pr_A_source_pins()
            for (rel in names(pins)) s44pr_a_assert(identical(s44pr_a_sha(file.path(source_root, rel)), unname(pins[[rel]])), 
                paste0("ORIGINAL_A_SOURCE_CHANGED:", rel))
            s44pr_a_assert(isTRUE(am$complete) && identical(am$model_id, i$model_id) && identical(am$source_hash, 
                i$source_hash) && identical(s44pr_a_vec(am$members), context$model_ids) && identical(am$model_A_key, 
                s44pr_A_original_key(context, source_root, methods)) && identical(am$file$sha256, cache_record$A_sha256), 
                "ORIGINAL_MODEL_A_KEY_OR_MANIFEST_CHANGED")
            if (!file.exists(cache_rds)) s44pr_a_fail("ORIGINAL_A_FILE_MISSING", rebuild_eligible = TRUE)
            info <- file.info(cache_rds)
            s44pr_a_assert(!is.na(info$size) && info$size == cache_record$A_bytes, "ORIGINAL_A_FILE_SIZE_CHANGED")
            if (is.null(transfer_receipt)) {
                s44pr_a_assert(identical(s44pr_a_sha(cache_rds), am$file$sha256), "ORIGINAL_A_FILE_SHA_CHANGED")
            } else {
                s44pr_a_assert(isTRUE(transfer_receipt$verified) && identical(transfer_receipt$sha256, am$file$sha256) && 
                  transfer_receipt$bytes == info$size && identical(normalizePath(transfer_receipt$path, mustWork = TRUE), 
                  normalizePath(cache_rds, mustWork = TRUE)) && !is.null(transfer_receipt$mtime_unix) && abs(as.numeric(info$mtime) - 
                  as.numeric(transfer_receipt$mtime_unix)) < 1e-06, "ONE_TIME_STREAM_HASH_AND_UNCHANGED_LOCAL_FILE_RECEIPT_REQUIRED")
            }
            A <- tryCatch(readRDS(cache_rds), error = function(e) s44pr_a_fail("ORIGINAL_A_RDS_UNREADABLE", conditionMessage(e), 
                rebuild_eligible = TRUE))
            s44pr_a_assert(is.list(A) && identical(names(A$models), context$model_ids) && identical(A$members, 
                context$members) && identical(A$conditional_order, context$model_ids), "ORIGINAL_A_MODEL_MEMBER_ORDER_OR_CONDITIONAL_LAWS_CHANGED")
            s44pr_a_assert(identical(A$body_pmm$revision, engine$body_pmm_revision()) && identical(A$body_pmm$prefix_binding, 
                engine$body_tqe_binding(context$d, context$members, context$base, context$w)) && identical(A$models$bmi$log_a_revision$revision, 
                "CORE39_BMI_LOG_FAT_LINEAR_V1") && isTRUE(A$body_pmm$raw_members_unchanged) && isTRUE(A$body_pmm$no_height) && 
                as.integer(A$body_pmm$group_uniform_slot) == 3L && identical(as.integer(A$body_pmm$unused_fixed_uniform_slots), 
                4:7), "ORIGINAL_BODY_PMM_BMI_AND_RNG_SLOT_BINDING_CHANGED")
            d <- context$d
            members <- context$members
            base <- context$base
            w <- context$w
            columns <- vapply(members, `[[`, "", "column")
            bodycols <- engine$body_tqe_columns(members)
            for (j in seq_along(A$models)) {
                f <- A$models[[j]]
                before <- columns[seq_len(j - 1L)]
                s44pr_a_assert(identical(f$member, members[[j]]), "ORIGINAL_A_MEMBER_DEFINITION_CHANGED")
                if (isTRUE(f$body_pmm_group_slot)) {
                  s44pr_a_assert(as.integer(f$original_uniform_slot) == j && isTRUE(f$uniform_used) == identical(names(members)[j], 
                    "weight"), "ORIGINAL_BODY_DONOR_SLOT_CHANGED")
                } else {
                  prepared <- engine$body_tqe_affine_prepare(d, before, bodycols, reference = f$body_tqe_affine, 
                    w = w)
                  expected <- engine$c39_a_design(prepared$data, base, prepared$previous)
                  actual <- engine$c39_a_design(prepared$data, NULL, NULL, f$design_reference)
                  s44pr_a_assert(identical(unname(f$predecessor_columns), unname(before)) && identical(engine$body_tqe_reference_copy(f$design_reference), 
                    engine$body_tqe_reference_copy(expected$reference)) && identical(actual$X, expected$X), "ORIGINAL_A_FORMULA_FACTOR_LEVELS_OR_DESIGN_CHANGED")
                }
            }
            group <- A$body_pmm$group
            enc <- engine$body_tqe_encode(d, members)
            gc <- group$columns
            tuple <- cbind(T = enc$total, q = enc$fraction, e = enc$z$e, hip = d[[gc[["hip"]]]], waist = d[[gc[["waist"]]]])
            s44pr_a_assert(identical(group$tuple, tuple) && as.integer(group$donor_N) == nrow(d) && identical(group$spec, 
                engine$body_pmm_spec()), "ORIGINAL_BODY_DONOR_VALUES_OR_METHOD_SPEC_CHANGED")
            upstream <- vapply(members[c("smoking", "bmi")], `[[`, "", "column")
            gx <- engine$c39_a_design(d, base, upstream)
            ga <- engine$c39_a_design(d, NULL, NULL, group$design_reference)
            s44pr_a_assert(identical(gx$X, ga$X) && identical(engine$body_tqe_reference_copy(gx$reference), engine$body_tqe_reference_copy(group$design_reference)), 
                "ORIGINAL_BODY_MATCHING_DESIGN_CHANGED")
            coverage <- integer()
            for (leaf in group$tree$leaves) {
                rows <- as.integer(leaf$rows)
                coverage <- c(coverage, rows)
                s44pr_a_assert(length(rows) == as.integer(leaf$N) && all(rows >= 1L & rows <= nrow(d)) && length(leaf$cdf) == 
                  length(rows) + 1L && leaf$cdf[1L] == 0 && tail(leaf$cdf, 1L) == 1 && all(is.finite(leaf$cdf)) && 
                  all(diff(leaf$cdf) >= 0), "ORIGINAL_BODY_DONOR_CDF_INVALID")
                expected_cdf <- c(0, cumsum(w[rows]/sum(w[rows])))
                expected_cdf[length(expected_cdf)] <- 1
                s44pr_a_assert(max(abs(expected_cdf - leaf$cdf)) <= 1e-12, "ORIGINAL_BODY_DONOR_WEIGHT_PROBABILITIES_CHANGED")
            }
            s44pr_a_assert(identical(sort(coverage), seq_len(nrow(d))), "ORIGINAL_BODY_DONOR_PARTITION_CHANGED")
            kind <- RNGkind()
            had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
            if (had_seed) seed <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
            on.exit({
                do.call(RNGkind, as.list(kind))
                if (had_seed) assign(".Random.seed", seed, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, 
                  inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
            }, add = TRUE)
            n <- min(as.integer(probe_rows), nrow(d))
            s44pr_a_assert(n >= 1L && n <= 8L, "BOUNDED_A_PORTABILITY_PROBE_REQUIRED")
            C <- d[seq_len(n), , drop = FALSE]
            u <- matrix(0.5, n, length(A$models))
            probe <- tryCatch(lapply(0:1, function(b) {
                value <- engine$c39_draw_a(A, C, b, u)
                engine$c39_support(value, members)
                s44pr_a_assert(nrow(value) == n && all(columns %in% names(value)), "A_PROBE_SHAPE_CHANGED")
                TRUE
            }), error = function(e) s44pr_a_fail("ORIGINAL_A_PORTABILITY_PROBE_FAILED", conditionMessage(e), rebuild_eligible = TRUE))
            list(A = A, provenance = list(schema = "PROTEIN_RIDGE_ORIGINAL_A_REUSE_CERTIFICATE_V1", passed = TRUE, 
                model_id = i$model_id, b = as.integer(i$b), original_model_A_key = am$model_A_key, A_sha256 = am$file$sha256, 
                original_ticket_sha256 = cache_record$files$ticket$sha256, original_A_manifest_sha256 = cache_record$files$A_manifest$sha256, 
                original_ticket_semantic_hash = s44pr_a_hash(cache_record$ticket), original_A_manifest_semantic_hash = s44pr_a_hash(am), 
                input_hash = i$input_hash, index_sha256 = context$draw$index_sha256, weights_sha256 = context$draw$weights_sha256, 
                sample_hash = context$pair_binding$sample_hash, family_hash = context$family_hash, risk_row_mapping_hash = context$risk_row_mapping_hash, 
                pair_binding = context$pair_binding, members = context$model_ids, rng = context$rng, original_R = cache_record$original_R, 
                current_R = R.version.string, source_pins = as.list(pins), A_refitted = FALSE, A_object_modified = FALSE, 
                Y_hash_not_an_A_key = TRUE, portability_probe_rows = n, portability_probe_exposure_values = 0:1, 
                official_MC_integrations = 0L, file_sha_verified_in_stream = !is.null(transfer_receipt)))
        }
    }), envir = context)
}

