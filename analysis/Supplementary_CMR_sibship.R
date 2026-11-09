# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

v2_cmr_migration_shape <- function(ctx) {
    co <- v2_read(ctx, "cohort")
    results <- list()
    regs <- list()
    gate <- tryCatch({
        v2_disease_gate(ctx)
        NULL
    }, error = identity)
    jobs <- list(list(key = "migration_association", family = "MIGRATION_ASSOCIATION", kind = "association", code = ""))
    for (k in c("250.2", "411", "351")) for (s in c("interaction_log_distance", "interaction_distance_group", "base_same_sample", 
        "plus_distance")) jobs[[length(jobs) + 1]] <- list(key = paste("migration", k, s, sep = "_"), family = paste0("MIGRATION_", 
        s), kind = s, code = k)
    reg <- bi_bind(c(list(data.frame(test_id = paste0("shape_", v2_codes()), family_id = "SHAPE_SIX", in_family = TRUE)), 
        lapply(jobs, function(j) data.frame(test_id = j$key, family_id = j$family, in_family = TRUE))))
    v2_register(ctx, "shape_migration", reg)
    for (k in v2_codes()) {
        key <- paste0("shape_", k)
        regs[[key]] <- data.frame(test_id = key, family_id = "SHAPE_SIX", in_family = TRUE)
        if (inherits(gate, "error")) {
            results[[key]] <- v2_row("nonlinearity", key, "SHAPE_SIX", k, status = "BLOCKED_OBSERVATION", reason = conditionMessage(gate))
            next
        }
        a <- v2_cp(ctx, paste0("shape_detail_", k), function() {
            d <- v2_disease_data(ctx, k)
            d <- d[d$S >= 3, , drop = FALSE]
            br_shape_test(ctx, d, "border_relative", k)
        })
        r <- v2_row("nonlinearity", key, "SHAPE_SIX", k)
        for (n in intersect(names(r), names(a$test))) r[[n]] <- a$test[[n]]
        r$test_id <- key
        r$family_id <- "SHAPE_SIX"
        r$cohort <- "S3_8"
        r$reason <- a$test$error
        r$method <- a$test$method
        results[[key]] <- r
        v2_csv(ctx, paste0(key, "_curve.csv"), a$curve)
    }
    raw <- tryCatch({
        m <- ctx$cfg$migration_reuse
        v2_assert(isTRUE(m$verified), paste("BLOCKED_COMPATIBLE_DISTANCE_CACHE:", m$reason))
        v2_assert(identical(v2_filehash(m$distance_file), m$sha256), "DISTANCE_CACHE_HASH_MISMATCH")
        q <- readRDS(m$distance_file)
        v2_assert(is.data.frame(q) && all(c("eid", "distance") %in% names(q)), "DISTANCE_CACHE_SCHEMA")
        v2_ids(q$eid)
        v2_assert(all(is.na(q$distance) | (is.finite(q$distance) & q$distance >= 0)), "DISTANCE_UNITS_OR_VALUES_INVALID")
        q
    }, error = identity)
    d <- v2_s8(co$data)
    available <- !inherits(raw, "error")
    if (available) {
        d$distance <- raw$distance[match(d$eid, raw$eid)]
        d$log_distance <- log1p(d$distance)
        d$distance_group <- cut(d$distance, c(-Inf, 25, 100, Inf), right = FALSE, labels = c("<25km", "25-<100km", 
            ">=100km"))
    }
    for (j in jobs) {
        regs[[j$key]] <- data.frame(test_id = j$key, family_id = j$family, in_family = TRUE)
        if (!available) {
            results[[j$key]] <- v2_row("migration", j$key, j$family, j$code, status = if (!isTRUE(ctx$cfg$migration_reuse$verified)) 
                "SKIPPED_UNAVAILABLE"
            else "BLOCKED_DISTANCE_SOURCE", reason = conditionMessage(raw))
            next
        }
        if (j$kind == "association") {
            d$y <- d$log_distance
            results[[j$key]] <- v2_model(ctx, d, "migration", j$key, j$family, type = "gaussian")
            next
        }
        if (inherits(gate, "error")) {
            results[[j$key]] <- v2_row("migration", j$key, j$family, j$code, status = "BLOCKED_OBSERVATION", reason = conditionMessage(gate))
            next
        }
        dd <- v2_disease_data(ctx, j$code)
        ix <- match(dd$eid, d$eid)
        dd$log_distance <- d$log_distance[ix]
        dd$distance_group <- d$distance_group[ix]
        if (startsWith(j$kind, "interaction_")) {
            r <- v2_row("migration", j$key, j$family, j$code)
            a <- bi_interaction(ctx, dd, "border_relative", sub("^interaction_", "", j$kind))
            for (n in intersect(names(r), names(a))) r[[n]] <- a[[n]]
            r$reason <- a$error
            results[[j$key]] <- r
        }
        else {
            s <- bi_sample(dd, bi_formula("border_relative", "log_distance"), "border_relative")
            results[[j$key]] <- v2_model(ctx, s$data, "migration", j$key, j$family, j$code, extra = if (j$kind == 
                "plus_distance") 
                "log_distance"
            else character())
        }
    }
    tab <- v2_adjust(bi_bind(results), reg)
    v2_csv(ctx, "SHAPE_MIGRATION_RESULTS.csv", tab)
    v2_cp(ctx, "results_shape_migration", function() tab)
    list(status = if (any(grepl("BLOCKED|ERROR|NOT_ESTIMABLE", tab$status))) "PARTIAL" else "PROCESSED", positions = nrow(tab), 
        valid_p = sum(is.finite(tab$p)))
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        v2_cmr_migration_shape <- function(ctx) {
            co <- v2_read(ctx, "cohort")
            results <- list()
            regs <- list()
            gate <- tryCatch({
                v2_disease_gate(ctx)
                NULL
            }, error = identity)
            jobs <- list(list(key = "migration_association", family = "MIGRATION_ASSOCIATION", kind = "association", 
                code = ""))
            for (k in c("250.2", "411", "351")) for (s in c("interaction_log_distance", "interaction_distance_group", 
                "base_same_sample", "plus_distance")) jobs[[length(jobs) + 1]] <- list(key = paste("migration", 
                k, s, sep = "_"), family = paste0("MIGRATION_", s), kind = s, code = k)
            reg <- bi_bind(c(list(data.frame(test_id = paste0("shape_", v2_codes()), family_id = "SHAPE_SIX", in_family = TRUE)), 
                lapply(jobs, function(j) data.frame(test_id = j$key, family_id = j$family, in_family = TRUE))))
            v2_register(ctx, "shape_migration", reg)
            for (k in v2_codes()) {
                key <- paste0("shape_", k)
                regs[[key]] <- data.frame(test_id = key, family_id = "SHAPE_SIX", in_family = TRUE)
                if (inherits(gate, "error")) {
                  results[[key]] <- v2_row("nonlinearity", key, "SHAPE_SIX", k, status = "BLOCKED_OBSERVATION", 
                    reason = conditionMessage(gate))
                  next
                }
                a <- v2_cp(ctx, paste0("shape_detail_", k), function() {
                  d <- v2_disease_data(ctx, k)
                  d <- d[d$S >= 3, , drop = FALSE]
                  br_shape_test(ctx, d, "border_relative", k)
                })
                r <- v2_row("nonlinearity", key, "SHAPE_SIX", k)
                for (n in intersect(names(r), names(a$test))) r[[n]] <- a$test[[n]]
                r$test_id <- key
                r$family_id <- "SHAPE_SIX"
                r$cohort <- "S3_8"
                r$reason <- a$test$error
                r$method <- a$test$method
                results[[key]] <- r
                v2_csv(ctx, paste0(key, "_curve.csv"), a$curve)
            }
            raw <- tryCatch({
                m <- ctx$cfg$migration_reuse
                v2_assert(isTRUE(m$verified), paste("BLOCKED_COMPATIBLE_DISTANCE_CACHE:", m$reason))
                v2_assert(identical(v2_filehash(m$distance_file), m$sha256), "DISTANCE_CACHE_HASH_MISMATCH")
                q <- readRDS(m$distance_file)
                v2_assert(is.data.frame(q) && all(c("eid", "distance") %in% names(q)), "DISTANCE_CACHE_SCHEMA")
                v2_ids(q$eid)
                v2_assert(all(is.na(q$distance) | (is.finite(q$distance) & q$distance >= 0)), "DISTANCE_UNITS_OR_VALUES_INVALID")
                q
            }, error = identity)
            d <- v2_s8(co$data)
            available <- !inherits(raw, "error")
            if (available) {
                d$distance <- raw$distance[match(d$eid, raw$eid)]
                d$log_distance <- log1p(d$distance)
                d$distance_group <- cut(d$distance, c(-Inf, 25, 100, Inf), right = FALSE, labels = c("<25km", "25-<100km", 
                  ">=100km"))
            }
            for (j in jobs) {
                regs[[j$key]] <- data.frame(test_id = j$key, family_id = j$family, in_family = TRUE)
                if (!available) {
                  results[[j$key]] <- v2_row("migration", j$key, j$family, j$code, status = if (!isTRUE(ctx$cfg$migration_reuse$verified)) "SKIPPED_UNAVAILABLE" else "BLOCKED_DISTANCE_SOURCE", 
                    reason = conditionMessage(raw))
                  next
                }
                if (j$kind == "association") {
                  d$y <- d$log_distance
                  results[[j$key]] <- v2_model(ctx, d, "migration", j$key, j$family, type = "gaussian")
                  next
                }
                if (inherits(gate, "error")) {
                  results[[j$key]] <- v2_row("migration", j$key, j$family, j$code, status = "BLOCKED_OBSERVATION", 
                    reason = conditionMessage(gate))
                  next
                }
                dd <- v2_disease_data(ctx, j$code)
                ix <- match(dd$eid, d$eid)
                dd$log_distance <- d$log_distance[ix]
                dd$distance_group <- d$distance_group[ix]
                if (startsWith(j$kind, "interaction_")) {
                  r <- v2_row("migration", j$key, j$family, j$code)
                  a <- bi_interaction(ctx, dd, "border_relative", sub("^interaction_", "", j$kind))
                  for (n in intersect(names(r), names(a))) r[[n]] <- a[[n]]
                  r$reason <- a$error
                  results[[j$key]] <- r
                } else {
                  s <- bi_sample(dd, bi_formula("border_relative", "log_distance"), "border_relative")
                  results[[j$key]] <- v2_model(ctx, s$data, "migration", j$key, j$family, j$code, extra = if (j$kind == 
                    "plus_distance") "log_distance" else character())
                }
            }
            tab <- v2_adjust(bi_bind(results), reg)
            v2_csv(ctx, "SHAPE_MIGRATION_RESULTS.csv", tab)
            v2_cp(ctx, "results_shape_migration", function() tab)
            list(status = if (any(grepl("BLOCKED|ERROR|NOT_ESTIMABLE", tab$status))) "PARTIAL" else "PROCESSED", 
                positions = nrow(tab), valid_p = sum(is.finite(tab$p)))
        }
    }), envir = context)
}

