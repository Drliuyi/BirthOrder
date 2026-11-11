# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

b16_reference_load <- ff_load

b16_reference_draw <- ff_w_draw

b16_reference_a <- ff_a_scenarios

ff_load <- function(root) {
    cfg <- ff_p_cfg(root)
    h <- ff_p_host()
    oldroot <- cfg$source_fullfollowup$code_paths[[h]]
    oldrun <- cfg$source_fullfollowup$run_paths[[h]]
    ff_p_assert(identical(ff_p_sha(file.path(oldroot, "config/runtime.json")), cfg$source_fullfollowup$config_sha256), 
        "B16_SOURCE_CONFIG_CHANGED")
    before <- Sys.getenv("FF_RUN_DIR", unset = NA_character_)
    on.exit(if (is.na(before)) Sys.unsetenv("FF_RUN_DIR") else Sys.setenv(FF_RUN_DIR = before), add = TRUE)
    Sys.setenv(FF_RUN_DIR = oldrun)
    ctx <- b16_reference_load(oldroot)
    ff_p_assert(identical(as.numeric(ctx$tau_report), 11), "B16_FROZEN_TAU_NOT11")
    lock <- jsonlite::fromJSON(file.path(root, "config/METHOD_LOCK.json"), simplifyVector = FALSE)
    ff_p_assert(identical(ctx$ff_method_hash, lock$scientific_method_hash) && identical(ctx$input_lock_sha256, 
        lock$source_input_lock_sha256), "B16_FROZEN_METHOD_INPUT_CHANGED")
    selected <- read.csv(file.path(root, cfg$active_approved_path), colClasses = "character", check.names = FALSE)
    ff_p_assert(!anyDuplicated(selected$indicator_id) && all(selected$indicator_id %in% ctx$scope$indicator_id), 
        "B16_CANDIDATE_REGISTRY_MISMATCH")
    ff_p_assert(isTRUE(cfg$candidate_ready) && nrow(selected) == 9L && setequal(selected$indicator_id, c("hip", 
        "nmr__nmr_23596", "nmr__nmr_23611", "nmr__nmr_23618", "nmr__nmr_23639", "olink__fcrl6", "olink__gzma", 
        "olink__gzmh", "olink__klrd1")), "ADD9_EXACT_NINE_REQUIRED")
    ctx$scope <- selected
    ctx$ff_config <- cfg
    ctx$run <- ff_p_env("FF16_RUN_DIR", cfg$run_paths[[h]])
    ff_p_assert(!identical(ctx$run, oldrun) && grepl("/border_fullfollowup_b1000_add9_20261030_v1/run_20261030_add9_v1", 
        ctx$run, fixed = TRUE) && identical(ctx$run, cfg$run_paths[[h]]), "B16_OUTPUT_MUST_BE_ISOLATED_NEW_RUN")
    ctx$ff_run <- ctx$run
    ctx$root <- root
    ctx$run_id <- cfg$run_id
    ctx$source_fullfollowup_run <- oldrun
    decision <- jsonlite::fromJSON(file.path(root, "config/V4_METHOD_DECISION.json"), simplifyVector = FALSE)
    adoption <- jsonlite::fromJSON(file.path(root, "config/RESEARCHER_V4_ADOPTION.json"), simplifyVector = FALSE)
    ff_p_assert(isTRUE(adoption$authorized$final_V4_method_adoption) && isTRUE(decision$v4_method_adopted) && identical(decision$method_version, 
        ff_ce_method_spec()$id) && identical(decision$method_spec_sha256, ff_ce_method_spec()$sha256) && identical(decision$new_method_hash, 
        ff_ce_method_spec()$sha256) && identical(decision$old_method_hash, ctx$ff_method_hash) && identical(decision$old_source_hash, 
        ctx$source_hash) && identical(decision$input_lock_sha256, ctx$input_lock_sha256), "V4_ADOPTED_METHOD_INPUT_BINDING_REQUIRED")
    for (n in names(decision$scientific_files_sha256)) ff_p_assert(identical(ff_p_sha(file.path(root, n)), decision$scientific_files_sha256[[n]]), 
        "V4_REVIEWED_SOURCE_CHANGED")
    ff_p_assert(setequal(names(decision$endpoint_rules), as.character(unlist(cfg$endpoints))), "V4_EXACT_FIVE_CAUSE_RULES_REQUIRED")
    ctx$frozen_A_method_hash <- ctx$ff_method_hash
    ctx$frozen_A_source_hash <- ctx$source_hash
    ctx$v4_decision <- decision
    ctx$ff_method_hash <- decision$new_method_hash
    ctx$source_hash <- decision$new_source_hash
    ctx$allow_missing_a_recompute <- FALSE
    ctx$ff_draw_cache <- new.env(parent = emptyenv())
    ctx$ff_bridge_cache <- new.env(parent = emptyenv())
    ctx <- b9_a_initialize(ctx, root)
    ctx
}

ff_w_draw <- function(ctx, bid) {
    ff_w_assert(bid %in% 0:1000 && bid == as.integer(bid), "B16_FIXED_ID_RANGE")
    if (bid <= 100L) {
        actual <- b16_reference_draw(ctx, bid)
        generated <- ctx$discrete$c40_draw(ctx, bid)
        ff_w_assert(identical(as.integer(actual$ix), as.integer(generated$ix)) && identical(actual$receipt$draw_hash, 
            generated$receipt$draw_hash) && identical(actual$receipt$multiplicity_hash, generated$receipt$multiplicity_hash), 
            "ADD9_RETAINED_DRAW_AND_FROZEN_GENERATOR_DIFFER")
        if (is.environment(ctx$add9_a)) {
            assign(as.character(bid), actual$receipt, ctx$add9_a$draws)
            assign(as.character(bid), b16_a_cache_rng(), ctx$add9_a$draw_rng)
        }
        return(actual)
    }
    actual <- ctx$discrete$c40_draw(ctx, bid)
    original_run <- ctx$run
    ctx$run <- file.path(original_run, "hosts", ff_p_host(), "draw_store")
    retained <- ctx$discrete$c40_persist_draw(ctx, actual)
    ctx$run <- original_run
    actual$receipt <- c(actual$receipt, retained)
    actual$receipt$private_sha256 <- retained$actual_draw_private_sha256
    actual$receipt$cached_draw_index_sha256 <- retained$actual_draw_private_sha256
    ff_w_assert(is.character(retained$actual_draw_private_sha256) && nchar(retained$actual_draw_private_sha256) == 
        64L, "B16_DRAW_SEAL_REQUIRED")
    actual$receipt$complete <- TRUE
    if (is.environment(ctx$add9_a)) {
        assign(as.character(bid), actual$receipt, ctx$add9_a$draws)
        assign(as.character(bid), b16_a_cache_rng(), ctx$add9_a$draw_rng)
    }
    actual
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        b16_reference_load <- ff_load
        b16_reference_draw <- ff_w_draw
        b16_reference_a <- ff_a_scenarios
        ff_load <- function(root) {
            cfg <- ff_p_cfg(root)
            h <- ff_p_host()
            oldroot <- cfg$source_fullfollowup$code_paths[[h]]
            oldrun <- cfg$source_fullfollowup$run_paths[[h]]
            ff_p_assert(identical(ff_p_sha(file.path(oldroot, "config/runtime.json")), cfg$source_fullfollowup$config_sha256), 
                "B16_SOURCE_CONFIG_CHANGED")
            before <- Sys.getenv("FF_RUN_DIR", unset = NA_character_)
            on.exit(if (is.na(before)) Sys.unsetenv("FF_RUN_DIR") else Sys.setenv(FF_RUN_DIR = before), add = TRUE)
            Sys.setenv(FF_RUN_DIR = oldrun)
            ctx <- b16_reference_load(oldroot)
            ff_p_assert(identical(as.numeric(ctx$tau_report), 11), "B16_FROZEN_TAU_NOT11")
            lock <- jsonlite::fromJSON(file.path(root, "config/METHOD_LOCK.json"), simplifyVector = FALSE)
            ff_p_assert(identical(ctx$ff_method_hash, lock$scientific_method_hash) && identical(ctx$input_lock_sha256, 
                lock$source_input_lock_sha256), "B16_FROZEN_METHOD_INPUT_CHANGED")
            selected <- read.csv(file.path(root, cfg$active_approved_path), colClasses = "character", check.names = FALSE)
            ff_p_assert(!anyDuplicated(selected$indicator_id) && all(selected$indicator_id %in% ctx$scope$indicator_id), 
                "B16_CANDIDATE_REGISTRY_MISMATCH")
            ff_p_assert(isTRUE(cfg$candidate_ready) && nrow(selected) == 9L && setequal(selected$indicator_id, 
                c("hip", "nmr__nmr_23596", "nmr__nmr_23611", "nmr__nmr_23618", "nmr__nmr_23639", "olink__fcrl6", 
                  "olink__gzma", "olink__gzmh", "olink__klrd1")), "ADD9_EXACT_NINE_REQUIRED")
            ctx$scope <- selected
            ctx$ff_config <- cfg
            ctx$run <- ff_p_env("FF16_RUN_DIR", cfg$run_paths[[h]])
            ff_p_assert(!identical(ctx$run, oldrun) && grepl("/border_fullfollowup_b1000_add9_20261030_v1/run_20261030_add9_v1", 
                ctx$run, fixed = TRUE) && identical(ctx$run, cfg$run_paths[[h]]), "B16_OUTPUT_MUST_BE_ISOLATED_NEW_RUN")
            ctx$ff_run <- ctx$run
            ctx$root <- root
            ctx$run_id <- cfg$run_id
            ctx$source_fullfollowup_run <- oldrun
            decision <- jsonlite::fromJSON(file.path(root, "config/V4_METHOD_DECISION.json"), simplifyVector = FALSE)
            adoption <- jsonlite::fromJSON(file.path(root, "config/RESEARCHER_V4_ADOPTION.json"), simplifyVector = FALSE)
            ff_p_assert(isTRUE(adoption$authorized$final_V4_method_adoption) && isTRUE(decision$v4_method_adopted) && 
                identical(decision$method_version, ff_ce_method_spec()$id) && identical(decision$method_spec_sha256, 
                ff_ce_method_spec()$sha256) && identical(decision$new_method_hash, ff_ce_method_spec()$sha256) && 
                identical(decision$old_method_hash, ctx$ff_method_hash) && identical(decision$old_source_hash, 
                ctx$source_hash) && identical(decision$input_lock_sha256, ctx$input_lock_sha256), "V4_ADOPTED_METHOD_INPUT_BINDING_REQUIRED")
            for (n in names(decision$scientific_files_sha256)) ff_p_assert(identical(ff_p_sha(file.path(root, n)), 
                decision$scientific_files_sha256[[n]]), "V4_REVIEWED_SOURCE_CHANGED")
            ff_p_assert(setequal(names(decision$endpoint_rules), as.character(unlist(cfg$endpoints))), "V4_EXACT_FIVE_CAUSE_RULES_REQUIRED")
            ctx$frozen_A_method_hash <- ctx$ff_method_hash
            ctx$frozen_A_source_hash <- ctx$source_hash
            ctx$v4_decision <- decision
            ctx$ff_method_hash <- decision$new_method_hash
            ctx$source_hash <- decision$new_source_hash
            ctx$allow_missing_a_recompute <- FALSE
            ctx$ff_draw_cache <- new.env(parent = emptyenv())
            ctx$ff_bridge_cache <- new.env(parent = emptyenv())
            ctx <- b9_a_initialize(ctx, root)
            ctx
        }
        ff_w_draw <- function(ctx, bid) {
            ff_w_assert(bid %in% 0:1000 && bid == as.integer(bid), "B16_FIXED_ID_RANGE")
            if (bid <= 100L) {
                actual <- b16_reference_draw(ctx, bid)
                generated <- ctx$discrete$c40_draw(ctx, bid)
                ff_w_assert(identical(as.integer(actual$ix), as.integer(generated$ix)) && identical(actual$receipt$draw_hash, 
                  generated$receipt$draw_hash) && identical(actual$receipt$multiplicity_hash, generated$receipt$multiplicity_hash), 
                  "ADD9_RETAINED_DRAW_AND_FROZEN_GENERATOR_DIFFER")
                if (is.environment(ctx$add9_a)) {
                  assign(as.character(bid), actual$receipt, ctx$add9_a$draws)
                  assign(as.character(bid), b16_a_cache_rng(), ctx$add9_a$draw_rng)
                }
                return(actual)
            }
            actual <- ctx$discrete$c40_draw(ctx, bid)
            original_run <- ctx$run
            ctx$run <- file.path(original_run, "hosts", ff_p_host(), "draw_store")
            retained <- ctx$discrete$c40_persist_draw(ctx, actual)
            ctx$run <- original_run
            actual$receipt <- c(actual$receipt, retained)
            actual$receipt$private_sha256 <- retained$actual_draw_private_sha256
            actual$receipt$cached_draw_index_sha256 <- retained$actual_draw_private_sha256
            ff_w_assert(is.character(retained$actual_draw_private_sha256) && nchar(retained$actual_draw_private_sha256) == 
                64L, "B16_DRAW_SEAL_REQUIRED")
            actual$receipt$complete <- TRUE
            if (is.environment(ctx$add9_a)) {
                assign(as.character(bid), actual$receipt, ctx$add9_a$draws)
                assign(as.character(bid), b16_a_cache_rng(), ctx$add9_a$draw_rng)
            }
            actual
        }
    }), envir = context)
}

