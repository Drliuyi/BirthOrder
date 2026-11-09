# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

v3_tdi_cuts <- function(ctx, d) {
    p <- ctx$cfg$sources$tdi_cuts
    v2_assert(!is.null(p) && file.exists(p), "FROZEN_TDI_CUTPOINTS_MISSING")
    z <- read.csv(p)
    v2_assert(all(z$status == "OK"), "FROZEN_TDI_CUTPOINTS_INVALID")
    as.numeric(z$value)
}

v3_strata <- function(ctx) v2_strata_ses(ctx)

v2_match_prepare <- function(ctx, co) {
    v2_cp(ctx, "matched_pairs", function() {
        m <- jsonlite::fromJSON(ctx$cfg$sources$matched_manifest)
        for (k in c("data", "pair_cache")) v2_assert(identical(v2_filehash(m[[paste0(k, "_path")]]), m[[paste0(k, 
            "_sha256")]]), "ORIGINAL_MATCH_SOURCE_CHANGED")
        old <- readRDS(m$data_path)$ctx$d
        pool <- bi_m_eligible(v2_s8(co$data))
        original <- readRDS(m$pair_cache_path)
        pairs <- original[original$eid_first %in% pool$eid & original$eid_last %in% pool$eid, , drop = FALSE]
        v2_assert(!anyDuplicated(c(pairs$eid_first, pairs$eid_last)), "MATCH_MEMBER_DUPLICATE")
        ids <- c(pairs$eid_first, pairs$eid_last)
        a <- match(ids, old$eid)
        b <- match(ids, pool$eid)
        v2_assert(!anyNA(a), "ORIGINAL_PAIR_NOT_IN_ORIGINAL_POOL")
        for (k in c("border1n", "S", "birth_year", "sex", "ethnicity")) v2_assert(identical(as.character(old[[k]][a]), 
            as.character(pool[[k]][b])), paste("MATCH_SOURCE_ELIGIBILITY_CHANGED", k))
        first <- match(pairs$eid_first, pool$eid)
        last <- match(pairs$eid_last, pool$eid)
        v2_assert(all(pool$border_relative[first] == 0 & pool$border_relative[last] == 1), "MATCH_ENDPOINT_DIRECTION")
        for (k in c("S", "birth_year", "sex", "ethnicity")) v2_assert(all(as.character(pool[[k]][first]) == as.character(pool[[k]][last])), 
            paste("MATCH_EXACTNESS", k))
        v2_csv(ctx, "MATCHED_DESIGN_REUSE.csv", data.frame(original_pairs = nrow(original), S8_pairs = nrow(pairs), 
            decision = "PRESERVE_OLD_PAIRS_FILTER_COMPLETE_S8_PAIRS_ONLY", seed = ctx$cfg$matching$seed, new_independent_roster_does_not_trigger_rematching = TRUE))
        list(pairs = pairs, status = "OK", warning = "", decision = "REUSED_OLD_S8_PAIRS")
    })
}

v3_matched_siblings <- function(ctx) {
    v2_matched_siblings(ctx)
    tab <- v2_read(ctx, "results_matched_siblings")
    co <- v2_read(ctx, "cohort")
    fam <- v2_unwrap(file.path(ctx$cfg$sources$sibling, "cache/cohort.private.rds"))$data
    ids <- intersect(as.character(fam$eid), v2_s8(co$data)$eid)
    fd <- co$data[match(ids, co$data$eid), ]
    fd$family_id <- fam$family_id[match(ids, fam$eid)]
    ep <- ep_endpoints(fd)
    ps <- v2_read(ctx, "legacy_states")
    for (k in v2_codes()) ep[[paste0("y_", k)]] <- ps$states[[k]][match(ep$eid, ps$ids)]
    reg <- data.frame(test_id = paste0("sibling_endpoint_", v2_codes()), family_id = "TRUE_SIBLING_ENDPOINT_SIX", 
        in_family = TRUE)
    v2_register(ctx, "sibling_endpoint", reg)
    rows <- flows <- list()
    for (k in v2_codes()) {
        z <- ep_outcome(ep, k)
        r <- v2_conditional(ctx, z$ready, paste0("sibling_endpoint_", k), "TRUE_SIBLING_ENDPOINT_SIX", k, "family_id", 
            c("border_relative", "age", "sex"))
        r$module <- "siblings_endpoint"
        r$design <- "TRUE_SIBLING_FIRST_LAST_ENDPOINT"
        r$families <- z$flow$complete_pairs
        r$participants <- 2 * z$flow$complete_pairs
        r$discordant_pairs <- z$flow$discordant_pairs
        rows[[k]] <- r
        flows[[k]] <- z$flow
    }
    ends <- v2_adjust(bi_bind(rows), reg)
    v2_csv(ctx, "TRUE_SIBLING_ENDPOINT_FLOW.csv", bi_bind(flows))
    v2_cp(ctx, "results_sibling_endpoint", function() ends)
    tab$design <- ifelse(tab$module == "matched", "MATCHED_FIRST_LAST_ENDPOINT", "TRUE_SIBLING_RELATIVE_ALL_OBSERVED_RANKS")
    tab$participants <- tab$families <- tab$pairs <- tab$discordant_pairs <- NA_real_
    for (i in seq_len(nrow(tab))) {
        d <- v2_disease_data(ctx, tab$phecode[i])
        if (tab$module[i] == "matched") {
            d <- bi_m_pair_data(d, v2_read(ctx, "matched_pairs")$pairs)
            tab$pairs[i] <- length(unique(d$pair_id))
            tab$discordant_pairs[i] <- tab$informative_families[i]
            tab$families[i] <- NA
        }
        else {
            d <- d[d$eid %in% ids, ]
            d$family_id <- fam$family_id[match(d$eid, fam$eid)]
            d <- bi_sample(d, y ~ border_relative + age + sex, "border_relative")$data
            tab$families[i] <- length(unique(d$family_id))
        }
        tab$participants[i] <- nrow(d)
    }
    v2_cp(ctx, "results_matching_all", function() bi_bind(list(tab, ends)))
    v2_csv(ctx, "MATCHING_ALL_RESULTS.csv", bi_bind(list(tab, ends)))
    list(status = "PROCESSED", positions = nrow(tab) + nrow(ends))
}

v3_continuous <- function(ctx, module) {
    e <- v2_feature_context(ctx)
    dd <- v2_feature_dictionary(ctx$root)
    mods <- if (module == "clinical") 
        c("clinical", "blood_pressure")
    else module
    dd <- dd[dd$module %in% mods, ]
    jobs <- lapply(seq_len(nrow(dd)), function(i) list(row = as.list(dd[i]), sensitivity = "primary"))
    if (module == "clinical") {
        v2_medication_flags(ctx)
        for (i in which(dd$feature_id %in% c("glucose", "hba1c") & dd$module == "clinical")) for (s in c("no_baseline_DM", 
            "no_glucose_drug")) jobs[[length(jobs) + 1L]] <- list(row = as.list(dd[i]), sensitivity = s)
        for (i in which(dd$module == "blood_pressure")) jobs[[length(jobs) + 1L]] <- list(row = as.list(dd[i]), 
            sensitivity = "no_BP_drug")
    }
    if (module == "cmr") {
        v2_read(ctx, "incident")
        for (i in seq_len(nrow(dd))) jobs[[length(jobs) + 1L]] <- list(row = as.list(dd[i]), sensitivity = "no_prescan_IHD_T2D_351")
    }
    fam <- function(j) {
        r <- j$row
        m <- r$module
        f <- if (m == "clinical") 
            if (r$feature_id %in% c("glucose", "hba1c")) 
                "GLYCAEMIA_2"
            else "BIOCHEMISTRY_7"
        else if (m == "blood_pressure") 
            if (tolower(r$feature_id) %in% c("sbp", "dbp")) 
                "BP_MAIN_2"
            else "BP_DERIVED_2"
        else if (m == "cmr") 
            if (isTRUE(r$indexed)) 
                "CMR_INDEXED_7"
            else "CMR_PRIMARY_82"
        else toupper(paste0(m, "_", if (m == "nmr") 
            251
        else 2923))
        if (j$sensitivity != "primary") 
            f <- paste(f, j$sensitivity, sep = "_")
        f
    }
    reg <- bi_bind(lapply(jobs, function(j) data.frame(test_id = paste("continuous", j$row$mediator_id, j$sensitivity, 
        sep = "_"), family_id = fam(j), in_family = TRUE)))
    v2_register(ctx, paste0("continuous_", module), reg)
    chunks <- split(seq_along(jobs), ceiling(seq_along(jobs)/25))
    out <- v2_parallel(ctx, chunks, function(ii) bi_bind(lapply(ii, function(i) v2_continuous_one(ctx, e, jobs[[i]]$row, 
        jobs[[i]]$sensitivity))))
    tab <- v2_adjust(bi_bind(out), reg)
    tab$measurement_time <- ifelse(tab$module == "cmr", "imaging_visit_2", "baseline_visit_0")
    v2_csv(ctx, paste0(toupper(module), "_ALL_RESULTS.csv"), tab)
    v2_cp(ctx, paste0("results_", module), function() tab)
    list(status = "PROCESSED", positions = nrow(tab), valid_p = sum(is.finite(tab$p)), failed = sum(!is.finite(tab$p)))
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        v3_tdi_cuts <- function(ctx, d) {
            p <- ctx$cfg$sources$tdi_cuts
            v2_assert(!is.null(p) && file.exists(p), "FROZEN_TDI_CUTPOINTS_MISSING")
            z <- read.csv(p)
            v2_assert(all(z$status == "OK"), "FROZEN_TDI_CUTPOINTS_INVALID")
            as.numeric(z$value)
        }
        v3_strata <- function(ctx) v2_strata_ses(ctx)
        v2_match_prepare <- function(ctx, co) {
            v2_cp(ctx, "matched_pairs", function() {
                m <- jsonlite::fromJSON(ctx$cfg$sources$matched_manifest)
                for (k in c("data", "pair_cache")) v2_assert(identical(v2_filehash(m[[paste0(k, "_path")]]), m[[paste0(k, 
                  "_sha256")]]), "ORIGINAL_MATCH_SOURCE_CHANGED")
                old <- readRDS(m$data_path)$ctx$d
                pool <- bi_m_eligible(v2_s8(co$data))
                original <- readRDS(m$pair_cache_path)
                pairs <- original[original$eid_first %in% pool$eid & original$eid_last %in% pool$eid, , drop = FALSE]
                v2_assert(!anyDuplicated(c(pairs$eid_first, pairs$eid_last)), "MATCH_MEMBER_DUPLICATE")
                ids <- c(pairs$eid_first, pairs$eid_last)
                a <- match(ids, old$eid)
                b <- match(ids, pool$eid)
                v2_assert(!anyNA(a), "ORIGINAL_PAIR_NOT_IN_ORIGINAL_POOL")
                for (k in c("border1n", "S", "birth_year", "sex", "ethnicity")) v2_assert(identical(as.character(old[[k]][a]), 
                  as.character(pool[[k]][b])), paste("MATCH_SOURCE_ELIGIBILITY_CHANGED", k))
                first <- match(pairs$eid_first, pool$eid)
                last <- match(pairs$eid_last, pool$eid)
                v2_assert(all(pool$border_relative[first] == 0 & pool$border_relative[last] == 1), "MATCH_ENDPOINT_DIRECTION")
                for (k in c("S", "birth_year", "sex", "ethnicity")) v2_assert(all(as.character(pool[[k]][first]) == 
                  as.character(pool[[k]][last])), paste("MATCH_EXACTNESS", k))
                v2_csv(ctx, "MATCHED_DESIGN_REUSE.csv", data.frame(original_pairs = nrow(original), S8_pairs = nrow(pairs), 
                  decision = "PRESERVE_OLD_PAIRS_FILTER_COMPLETE_S8_PAIRS_ONLY", seed = ctx$cfg$matching$seed, 
                  new_independent_roster_does_not_trigger_rematching = TRUE))
                list(pairs = pairs, status = "OK", warning = "", decision = "REUSED_OLD_S8_PAIRS")
            })
        }
        v3_matched_siblings <- function(ctx) {
            v2_matched_siblings(ctx)
            tab <- v2_read(ctx, "results_matched_siblings")
            co <- v2_read(ctx, "cohort")
            fam <- v2_unwrap(file.path(ctx$cfg$sources$sibling, "cache/cohort.private.rds"))$data
            ids <- intersect(as.character(fam$eid), v2_s8(co$data)$eid)
            fd <- co$data[match(ids, co$data$eid), ]
            fd$family_id <- fam$family_id[match(ids, fam$eid)]
            ep <- ep_endpoints(fd)
            ps <- v2_read(ctx, "legacy_states")
            for (k in v2_codes()) ep[[paste0("y_", k)]] <- ps$states[[k]][match(ep$eid, ps$ids)]
            reg <- data.frame(test_id = paste0("sibling_endpoint_", v2_codes()), family_id = "TRUE_SIBLING_ENDPOINT_SIX", 
                in_family = TRUE)
            v2_register(ctx, "sibling_endpoint", reg)
            rows <- flows <- list()
            for (k in v2_codes()) {
                z <- ep_outcome(ep, k)
                r <- v2_conditional(ctx, z$ready, paste0("sibling_endpoint_", k), "TRUE_SIBLING_ENDPOINT_SIX", 
                  k, "family_id", c("border_relative", "age", "sex"))
                r$module <- "siblings_endpoint"
                r$design <- "TRUE_SIBLING_FIRST_LAST_ENDPOINT"
                r$families <- z$flow$complete_pairs
                r$participants <- 2 * z$flow$complete_pairs
                r$discordant_pairs <- z$flow$discordant_pairs
                rows[[k]] <- r
                flows[[k]] <- z$flow
            }
            ends <- v2_adjust(bi_bind(rows), reg)
            v2_csv(ctx, "TRUE_SIBLING_ENDPOINT_FLOW.csv", bi_bind(flows))
            v2_cp(ctx, "results_sibling_endpoint", function() ends)
            tab$design <- ifelse(tab$module == "matched", "MATCHED_FIRST_LAST_ENDPOINT", "TRUE_SIBLING_RELATIVE_ALL_OBSERVED_RANKS")
            tab$participants <- tab$families <- tab$pairs <- tab$discordant_pairs <- NA_real_
            for (i in seq_len(nrow(tab))) {
                d <- v2_disease_data(ctx, tab$phecode[i])
                if (tab$module[i] == "matched") {
                  d <- bi_m_pair_data(d, v2_read(ctx, "matched_pairs")$pairs)
                  tab$pairs[i] <- length(unique(d$pair_id))
                  tab$discordant_pairs[i] <- tab$informative_families[i]
                  tab$families[i] <- NA
                } else {
                  d <- d[d$eid %in% ids, ]
                  d$family_id <- fam$family_id[match(d$eid, fam$eid)]
                  d <- bi_sample(d, y ~ border_relative + age + sex, "border_relative")$data
                  tab$families[i] <- length(unique(d$family_id))
                }
                tab$participants[i] <- nrow(d)
            }
            v2_cp(ctx, "results_matching_all", function() bi_bind(list(tab, ends)))
            v2_csv(ctx, "MATCHING_ALL_RESULTS.csv", bi_bind(list(tab, ends)))
            list(status = "PROCESSED", positions = nrow(tab) + nrow(ends))
        }
        v3_continuous <- function(ctx, module) {
            e <- v2_feature_context(ctx)
            dd <- v2_feature_dictionary(ctx$root)
            mods <- if (module == "clinical") c("clinical", "blood_pressure") else module
            dd <- dd[dd$module %in% mods, ]
            jobs <- lapply(seq_len(nrow(dd)), function(i) list(row = as.list(dd[i]), sensitivity = "primary"))
            if (module == "clinical") {
                v2_medication_flags(ctx)
                for (i in which(dd$feature_id %in% c("glucose", "hba1c") & dd$module == "clinical")) for (s in c("no_baseline_DM", 
                  "no_glucose_drug")) jobs[[length(jobs) + 1L]] <- list(row = as.list(dd[i]), sensitivity = s)
                for (i in which(dd$module == "blood_pressure")) jobs[[length(jobs) + 1L]] <- list(row = as.list(dd[i]), 
                  sensitivity = "no_BP_drug")
            }
            if (module == "cmr") {
                v2_read(ctx, "incident")
                for (i in seq_len(nrow(dd))) jobs[[length(jobs) + 1L]] <- list(row = as.list(dd[i]), sensitivity = "no_prescan_IHD_T2D_351")
            }
            fam <- function(j) {
                r <- j$row
                m <- r$module
                f <- if (m == "clinical") if (r$feature_id %in% c("glucose", "hba1c")) "GLYCAEMIA_2" else "BIOCHEMISTRY_7" else if (m == 
                  "blood_pressure") if (tolower(r$feature_id) %in% c("sbp", "dbp")) "BP_MAIN_2" else "BP_DERIVED_2" else if (m == 
                  "cmr") if (isTRUE(r$indexed)) "CMR_INDEXED_7" else "CMR_PRIMARY_82" else toupper(paste0(m, "_", 
                  if (m == "nmr") 251 else 2923))
                if (j$sensitivity != "primary") f <- paste(f, j$sensitivity, sep = "_")
                f
            }
            reg <- bi_bind(lapply(jobs, function(j) data.frame(test_id = paste("continuous", j$row$mediator_id, 
                j$sensitivity, sep = "_"), family_id = fam(j), in_family = TRUE)))
            v2_register(ctx, paste0("continuous_", module), reg)
            chunks <- split(seq_along(jobs), ceiling(seq_along(jobs)/25))
            out <- v2_parallel(ctx, chunks, function(ii) bi_bind(lapply(ii, function(i) v2_continuous_one(ctx, 
                e, jobs[[i]]$row, jobs[[i]]$sensitivity))))
            tab <- v2_adjust(bi_bind(out), reg)
            tab$measurement_time <- ifelse(tab$module == "cmr", "imaging_visit_2", "baseline_visit_0")
            v2_csv(ctx, paste0(toupper(module), "_ALL_RESULTS.csv"), tab)
            v2_cp(ctx, paste0("results_", module), function() tab)
            list(status = "PROCESSED", positions = nrow(tab), valid_p = sum(is.finite(tab$p)), failed = sum(!is.finite(tab$p)))
        }
    }), envir = context)
}

