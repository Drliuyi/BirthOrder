# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

v2_assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)

v2_codes <- function() c("250", "250.2", "411", "411.3", "411.4", "351")

v2_ids <- function(x) {
    z <- as.character(x)
    v2_assert(!anyNA(z) && !anyDuplicated(z) && all(nzchar(z) & trimws(z) == z), "INVALID_PARTICIPANT_IDS")
    z
}

v2_hash <- function(x) digest::digest(x, algo = "sha256")

v2_filehash <- function(p) digest::digest(file = p, algo = "sha256")

v2_unwrap <- function(p) {
    z <- readRDS(p)
    if (is.list(z) && all(c("value", "hash", "signature") %in% names(z))) {
        h <- z$hash
        z$hash <- NULL
        v2_assert(identical(h, v2_hash(z)), paste("CACHE_CONTENT_HASH_MISMATCH", p))
        return(z$value)
    }
    z
}

v2_publish <- function(x, p) {
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(p)) {
        v2_assert(identical(readRDS(p), x), paste("REFUSE_OVERWRITE", p))
        return(invisible(p))
    }
    t <- paste0(p, ".", Sys.getpid(), ".partial")
    saveRDS(x, t)
    Sys.chmod(t, "0600")
    v2_assert(file.rename(t, p), "ATOMIC_PUBLICATION_FAILED")
    invisible(p)
}

v2_cp <- function(ctx, key, fn) {
    v2_assert(grepl("^[A-Za-z0-9_.-]+$", key), "INVALID_CHECKPOINT_KEY")
    p <- file.path(ctx$out, "private", "checkpoints", paste0(key, ".rds"))
    if (file.exists(p)) 
        return(v2_read(ctx, key))
    value <- fn()
    v2_publish(list(signature = ctx$signature, value_hash = v2_hash(value), value = value), p)
    value
}

v2_csv <- function(ctx, name, x, private = FALSE) {
    v2_assert(is.data.frame(x), "TABLE_REQUIRED")
    p <- file.path(ctx$out, if (private) 
        "private"
    else "public/raw", name)
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    if (!private) 
        x <- bi_public(x, 5L)
    t <- paste0(p, ".", Sys.getpid(), ".partial")
    write.csv(x, t, row.names = FALSE, na = "NA")
    if (file.exists(p)) {
        v2_assert(identical(v2_filehash(t), v2_filehash(p)), paste("REFUSE_CHANGED_TABLE", p))
        unlink(t)
        return(invisible(p))
    }
    v2_assert(file.rename(t, p), "TABLE_PUBLICATION_FAILED")
    if (private) 
        Sys.chmod(p, "0600")
    invisible(p)
}

v2_json <- function(ctx, name, x) {
    p <- file.path(ctx$out, name)
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(p)) {
        v2_assert(identical(jsonlite::fromJSON(p, simplifyVector = FALSE), x), "REFUSE_CHANGED_JSON")
        return()
    }
    jsonlite::write_json(x, p, auto_unbox = TRUE, pretty = TRUE, na = "null")
}

v2_read <- function(ctx, key) {
    cacheable <- key %in% c("cohort", "phenotypes", "incident", "baseline_medication")
    if (is.null(ctx$read_cache)) 
        ctx$read_cache <- new.env(parent = emptyenv())
    if (cacheable && exists(key, ctx$read_cache, inherits = FALSE)) 
        return(get(key, ctx$read_cache))
    p <- file.path(ctx$out, "private", "checkpoints", paste0(key, ".rds"))
    v2_assert(file.exists(p), paste("DEPENDENCY_NOT_READY", key))
    z <- readRDS(p)
    v2_assert(identical(z$signature, ctx$signature) && identical(v2_hash(z$value), z$value_hash), "STALE_DEPENDENCY")
    if (cacheable) 
        assign(key, z$value, ctx$read_cache)
    z$value
}

v2_workers <- function(ctx) {
    cores <- parallel::detectCores(logical = FALSE)
    if (is.na(cores)) 
        cores <- 1L
    free <- 2
    if (file.exists("/proc/meminfo")) {
        s <- readLines("/proc/meminfo")
        free <- as.numeric(gsub("[^0-9]", "", s[grepl("^MemAvailable:", s)]))/1024^2
    }
    n <- min(ctx$cfg$runtime$workers, ctx$cfg$runtime$maximum_workers, max(1, cores - 2), floor(max(0, free - ctx$cfg$runtime$minimum_free_gib)/ctx$cfg$runtime$worker_gib))
    v2_assert(n >= 1, "WAIT_RESOURCE: insufficient memory to start even one model; resume after releasing memory")
    as.integer(n)
}

v2_parallel <- function(ctx, keys, fn) {
    ans <- vector("list", length(keys))
    i <- 1L
    while (i <= length(keys)) {
        n <- v2_workers(ctx)
        ix <- seq.int(i, min(length(keys), i + n - 1L))
        z <- parallel::mclapply(keys[ix], fn, mc.cores = n, mc.preschedule = FALSE, mc.set.seed = FALSE)
        v2_assert(!any(vapply(z, inherits, logical(1), "try-error")), "WORKER_EXECUTION_ERROR")
        ans[ix] <- z
        i <- max(ix) + 1L
    }
    ans
}

v2_s8 <- function(d) d[is.finite(d$border_relative) & is.finite(d$S) & d$S >= 2 & d$S <= 8, , drop = FALSE]

v2_row <- function(module, key, family, trait = "", feature = "", status = "NOT_ESTIMABLE", reason = "") data.frame(module = module, 
    test_id = key, family_id = family, phecode = trait, feature_id = feature, exposure = "border_relative", cohort = "S2_8", 
    N = NA_real_, cases = NA_real_, events = NA_real_, beta = NA_real_, se = NA_real_, effect = NA_real_, CI_low = NA_real_, 
    CI_high = NA_real_, p = NA_real_, p_raw = NA_real_, status = status, reason = reason, method = "", scale = "", 
    formula = "", input_signature = "", result_origin = "REESTIMATED_V2", stringsAsFactors = FALSE)

v2_model <- function(ctx, d, module, key, family, trait = "", feature = "", extra = character(), omit = character(), 
    type = "binomial") {
    v2_cp(ctx, paste0("model_", key), function() {
        r <- v2_row(module, key, family, trait, feature)
        s <- bi_sample(d, bi_formula("border_relative", extra, omit), "border_relative")
        r$input_signature <- v2_hash(list(data = s$data[, unique(c("eid", all.vars(s$formula))), drop = FALSE], 
            formula = deparse(s$formula), family = type, version = ctx$signature))
        a <- bi_fit(ctx, d, "border_relative", extra, omit, family = type)
        for (k in intersect(names(r), names(a))) r[[k]] <- a[[k]]
        r$effect <- if (type == "binomial") 
            a$OR
        else a$beta
        r$method <- if (type == "binomial") 
            "ordinary_logistic_Wald"
        else "linear_regression_t"
        r$scale <- if (type == "binomial") 
            "OR_0_to_1"
        else "beta_0_to_1"
        r$reason <- a$error
        r$warning <- a$warning
        r$rank <- a$rank
        r$columns <- a$columns
        if (type == "binomial" && is.finite(r$p)) {
            suspect <- isTRUE(a$coefficient_numeric_flag) || grepl("converg|0 or 1|singular", a$warning, ignore.case = TRUE)
            r$separation_diagnostic <- if (suspect) 
                "SUSPECTED_NUMERIC_REVIEW_NOT_LP"
            else "NO_NUMERICAL_FLAG_NOT_STRICT_LP"
            if (suspect) {
                r$status <- "INFERENCE_PENDING_SEPARATION"
                r$p <- NA_real_
            }
        }
        r
    })
}

v2_adjust <- function(x, registry) {
    if (!nrow(x)) 
        return(x)
    v2_assert(!anyDuplicated(x$test_id), "DUPLICATE_RESULT_ID")
    v2_assert(all(is.na(x$p) | (is.finite(x$p) & x$p >= 0 & x$p <= 1)), "INVALID_P_VALUE")
    x$q_primary <- x$p_bonferroni <- NA_real_
    x$M <- NA_integer_
    for (f in unique(x$family_id)) {
        ix <- which(x$family_id == f)
        rr <- registry[registry$family_id == f, , drop = FALSE]
        v2_assert(nrow(rr) > 0 && !anyDuplicated(rr$test_id), "MISSING_OR_DUPLICATE_FAMILY_REGISTRY")
        v2_assert(all(x$test_id[ix] %in% rr$test_id), "UNREGISTERED_RESULT")
        M <- sum(rr$in_family)
        x$M[ix] <- M
        ok <- ix[is.finite(x$p[ix]) & x$test_id[ix] %in% rr$test_id[rr$in_family]]
        if (length(ok)) {
            x$q_primary[ok] <- p.adjust(x$p[ok], "BH", n = M)
            x$p_bonferroni[ok] <- pmin(1, x$p[ok] * M)
        }
    }
    x
}

v2_register <- function(ctx, name, x) {
    v2_assert(!anyDuplicated(x$test_id), "REGISTRY_DUPLICATE")
    v2_cp(ctx, paste0("registry_", name), function() x)
    invisible(x)
}

v2_block <- function(msg) stop(paste0("BLOCKED: ", msg), call. = FALSE)

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        v2_assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
        v2_codes <- function() c("250", "250.2", "411", "411.3", "411.4", "351")
        v2_ids <- function(x) {
            z <- as.character(x)
            v2_assert(!anyNA(z) && !anyDuplicated(z) && all(nzchar(z) & trimws(z) == z), "INVALID_PARTICIPANT_IDS")
            z
        }
        v2_hash <- function(x) digest::digest(x, algo = "sha256")
        v2_filehash <- function(p) digest::digest(file = p, algo = "sha256")
        v2_unwrap <- function(p) {
            z <- readRDS(p)
            if (is.list(z) && all(c("value", "hash", "signature") %in% names(z))) {
                h <- z$hash
                z$hash <- NULL
                v2_assert(identical(h, v2_hash(z)), paste("CACHE_CONTENT_HASH_MISMATCH", p))
                return(z$value)
            }
            z
        }
        v2_publish <- function(x, p) {
            dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
            if (file.exists(p)) {
                v2_assert(identical(readRDS(p), x), paste("REFUSE_OVERWRITE", p))
                return(invisible(p))
            }
            t <- paste0(p, ".", Sys.getpid(), ".partial")
            saveRDS(x, t)
            Sys.chmod(t, "0600")
            v2_assert(file.rename(t, p), "ATOMIC_PUBLICATION_FAILED")
            invisible(p)
        }
        v2_cp <- function(ctx, key, fn) {
            v2_assert(grepl("^[A-Za-z0-9_.-]+$", key), "INVALID_CHECKPOINT_KEY")
            p <- file.path(ctx$out, "private", "checkpoints", paste0(key, ".rds"))
            if (file.exists(p)) return(v2_read(ctx, key))
            value <- fn()
            v2_publish(list(signature = ctx$signature, value_hash = v2_hash(value), value = value), p)
            value
        }
        v2_csv <- function(ctx, name, x, private = FALSE) {
            v2_assert(is.data.frame(x), "TABLE_REQUIRED")
            p <- file.path(ctx$out, if (private) "private" else "public/raw", name)
            dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
            if (!private) x <- bi_public(x, 5L)
            t <- paste0(p, ".", Sys.getpid(), ".partial")
            write.csv(x, t, row.names = FALSE, na = "NA")
            if (file.exists(p)) {
                v2_assert(identical(v2_filehash(t), v2_filehash(p)), paste("REFUSE_CHANGED_TABLE", p))
                unlink(t)
                return(invisible(p))
            }
            v2_assert(file.rename(t, p), "TABLE_PUBLICATION_FAILED")
            if (private) Sys.chmod(p, "0600")
            invisible(p)
        }
        v2_json <- function(ctx, name, x) {
            p <- file.path(ctx$out, name)
            dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
            if (file.exists(p)) {
                v2_assert(identical(jsonlite::fromJSON(p, simplifyVector = FALSE), x), "REFUSE_CHANGED_JSON")
                return()
            }
            jsonlite::write_json(x, p, auto_unbox = TRUE, pretty = TRUE, na = "null")
        }
        v2_read <- function(ctx, key) {
            cacheable <- key %in% c("cohort", "phenotypes", "incident", "baseline_medication")
            if (is.null(ctx$read_cache)) ctx$read_cache <- new.env(parent = emptyenv())
            if (cacheable && exists(key, ctx$read_cache, inherits = FALSE)) return(get(key, ctx$read_cache))
            p <- file.path(ctx$out, "private", "checkpoints", paste0(key, ".rds"))
            v2_assert(file.exists(p), paste("DEPENDENCY_NOT_READY", key))
            z <- readRDS(p)
            v2_assert(identical(z$signature, ctx$signature) && identical(v2_hash(z$value), z$value_hash), "STALE_DEPENDENCY")
            if (cacheable) assign(key, z$value, ctx$read_cache)
            z$value
        }
        v2_workers <- function(ctx) {
            cores <- parallel::detectCores(logical = FALSE)
            if (is.na(cores)) cores <- 1L
            free <- 2
            if (file.exists("/proc/meminfo")) {
                s <- readLines("/proc/meminfo")
                free <- as.numeric(gsub("[^0-9]", "", s[grepl("^MemAvailable:", s)]))/1024^2
            }
            n <- min(ctx$cfg$runtime$workers, ctx$cfg$runtime$maximum_workers, max(1, cores - 2), floor(max(0, 
                free - ctx$cfg$runtime$minimum_free_gib)/ctx$cfg$runtime$worker_gib))
            v2_assert(n >= 1, "WAIT_RESOURCE: insufficient memory to start even one model; resume after releasing memory")
            as.integer(n)
        }
        v2_parallel <- function(ctx, keys, fn) {
            ans <- vector("list", length(keys))
            i <- 1L
            while (i <= length(keys)) {
                n <- v2_workers(ctx)
                ix <- seq.int(i, min(length(keys), i + n - 1L))
                z <- parallel::mclapply(keys[ix], fn, mc.cores = n, mc.preschedule = FALSE, mc.set.seed = FALSE)
                v2_assert(!any(vapply(z, inherits, logical(1), "try-error")), "WORKER_EXECUTION_ERROR")
                ans[ix] <- z
                i <- max(ix) + 1L
            }
            ans
        }
        v2_s8 <- function(d) d[is.finite(d$border_relative) & is.finite(d$S) & d$S >= 2 & d$S <= 8, , drop = FALSE]
        v2_row <- function(module, key, family, trait = "", feature = "", status = "NOT_ESTIMABLE", reason = "") data.frame(module = module, 
            test_id = key, family_id = family, phecode = trait, feature_id = feature, exposure = "border_relative", 
            cohort = "S2_8", N = NA_real_, cases = NA_real_, events = NA_real_, beta = NA_real_, se = NA_real_, 
            effect = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_, p_raw = NA_real_, status = status, 
            reason = reason, method = "", scale = "", formula = "", input_signature = "", result_origin = "REESTIMATED_V2", 
            stringsAsFactors = FALSE)
        v2_model <- function(ctx, d, module, key, family, trait = "", feature = "", extra = character(), omit = character(), 
            type = "binomial") {
            v2_cp(ctx, paste0("model_", key), function() {
                r <- v2_row(module, key, family, trait, feature)
                s <- bi_sample(d, bi_formula("border_relative", extra, omit), "border_relative")
                r$input_signature <- v2_hash(list(data = s$data[, unique(c("eid", all.vars(s$formula))), drop = FALSE], 
                  formula = deparse(s$formula), family = type, version = ctx$signature))
                a <- bi_fit(ctx, d, "border_relative", extra, omit, family = type)
                for (k in intersect(names(r), names(a))) r[[k]] <- a[[k]]
                r$effect <- if (type == "binomial") a$OR else a$beta
                r$method <- if (type == "binomial") "ordinary_logistic_Wald" else "linear_regression_t"
                r$scale <- if (type == "binomial") "OR_0_to_1" else "beta_0_to_1"
                r$reason <- a$error
                r$warning <- a$warning
                r$rank <- a$rank
                r$columns <- a$columns
                if (type == "binomial" && is.finite(r$p)) {
                  suspect <- isTRUE(a$coefficient_numeric_flag) || grepl("converg|0 or 1|singular", a$warning, 
                    ignore.case = TRUE)
                  r$separation_diagnostic <- if (suspect) "SUSPECTED_NUMERIC_REVIEW_NOT_LP" else "NO_NUMERICAL_FLAG_NOT_STRICT_LP"
                  if (suspect) {
                    r$status <- "INFERENCE_PENDING_SEPARATION"
                    r$p <- NA_real_
                  }
                }
                r
            })
        }
        v2_adjust <- function(x, registry) {
            if (!nrow(x)) return(x)
            v2_assert(!anyDuplicated(x$test_id), "DUPLICATE_RESULT_ID")
            v2_assert(all(is.na(x$p) | (is.finite(x$p) & x$p >= 0 & x$p <= 1)), "INVALID_P_VALUE")
            x$q_primary <- x$p_bonferroni <- NA_real_
            x$M <- NA_integer_
            for (f in unique(x$family_id)) {
                ix <- which(x$family_id == f)
                rr <- registry[registry$family_id == f, , drop = FALSE]
                v2_assert(nrow(rr) > 0 && !anyDuplicated(rr$test_id), "MISSING_OR_DUPLICATE_FAMILY_REGISTRY")
                v2_assert(all(x$test_id[ix] %in% rr$test_id), "UNREGISTERED_RESULT")
                M <- sum(rr$in_family)
                x$M[ix] <- M
                ok <- ix[is.finite(x$p[ix]) & x$test_id[ix] %in% rr$test_id[rr$in_family]]
                if (length(ok)) {
                  x$q_primary[ok] <- p.adjust(x$p[ok], "BH", n = M)
                  x$p_bonferroni[ok] <- pmin(1, x$p[ok] * M)
                }
            }
            x
        }
        v2_register <- function(ctx, name, x) {
            v2_assert(!anyDuplicated(x$test_id), "REGISTRY_DUPLICATE")
            v2_cp(ctx, paste0("registry_", name), function() x)
            invisible(x)
        }
        v2_block <- function(msg) stop(paste0("BLOCKED: ", msg), call. = FALSE)
    }), envir = context)
}

