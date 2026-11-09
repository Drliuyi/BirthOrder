# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

fd_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

fd_env <- function(k, default = "") {
    x <- Sys.getenv(k, "")
    if (nzchar(x)) 
        x
    else default
}

fd_hash <- function(x) digest::digest(x, algo = "sha256")

fd_sha <- function(p) {
    fd_assert(file.exists(p) && !dir.exists(p), paste("REQUIRED_FILE_MISSING", p))
    digest::digest(file = p, algo = "sha256")
}

fd_json <- function(x, path) {
    fd_assert(!file.exists(path), paste("PRESERVE_EXISTING", path))
    tmp <- paste0(path, ".pending_", Sys.getpid())
    fd_assert(!file.exists(tmp), "PRESERVE_PENDING_ARTIFACT")
    jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
    fd_assert(!file.exists(path) && file.rename(tmp, path), "ATOMIC_PUBLICATION_FAILED")
}

fd_read <- function(p) jsonlite::fromJSON(p, simplifyVector = FALSE)

fd_capture <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(expr, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    list(value = value, warnings = unique(warnings))
}

fd_ids <- function(model) as.character(unlist(if (!is.null(model$indicator_ids)) model$indicator_ids else model$members))

fd_model <- function(cfg, id) {
    z <- Filter(function(m) identical(m$model_id, id), cfg$models)
    fd_assert(length(z) == 1L, "MODEL_NOT_UNIQUELY_FROZEN")
    z <- z[[1]]
    ids <- fd_ids(z)
    fd_assert(length(ids) > 0 && !anyDuplicated(ids) && setequal(ids, unlist(z$generation_order)), "INVALID_FROZEN_GENERATION_ORDER")
    z
}

fd_rng <- function(seed, code) {
    kind <- RNGkind()
    had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    if (had) 
        old <- get(".Random.seed", .GlobalEnv)
    on.exit({
        do.call(RNGkind, as.list(kind))
        if (had) assign(".Random.seed", old, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", 
            envir = .GlobalEnv)
    })
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    set.seed(seed)
    force(code)
}

fd_seed <- function(b, stream = 1L) as.integer(20261031L + as.integer(b) * 1009L + as.integer(stream) * 100003L)

fd_outdir <- function(run, out) {
    run <- normalizePath(run, mustWork = TRUE)
    fd_assert(!grepl("(^|/)[.][.](/|$)", out), "UNSAFE_OUTPUT_PATH")
    fd_assert(startsWith(paste0(out, "/"), paste0(run, "/")), "OUTPUT_OUTSIDE_ISOLATED_RUN")
    if (!dir.exists(dirname(out))) 
        dir.create(dirname(out), recursive = TRUE)
    fd_assert(startsWith(paste0(normalizePath(dirname(out)), "/"), paste0(run, "/")), "OUTPUT_SYMLINK_ESCAPES_RUN")
    if (dir.exists(out)) {
        allowed <- c("process.log", "PROCESS.json", "command.json", "COMMAND.json")
        fd_assert(all(list.files(out, all.files = TRUE, no.. = TRUE) %in% allowed), "OUTPUT_ATTEMPT_CONTAINS_PRIOR_SCIENTIFIC_ARTIFACT")
    }
    else fd_assert(dir.create(out), "OUTPUT_CREATE_FAILED")
    out
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        fd_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
        fd_env <- function(k, default = "") {
            x <- Sys.getenv(k, "")
            if (nzchar(x)) x else default
        }
        fd_hash <- function(x) digest::digest(x, algo = "sha256")
        fd_sha <- function(p) {
            fd_assert(file.exists(p) && !dir.exists(p), paste("REQUIRED_FILE_MISSING", p))
            digest::digest(file = p, algo = "sha256")
        }
        fd_json <- function(x, path) {
            fd_assert(!file.exists(path), paste("PRESERVE_EXISTING", path))
            tmp <- paste0(path, ".pending_", Sys.getpid())
            fd_assert(!file.exists(tmp), "PRESERVE_PENDING_ARTIFACT")
            jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
            fd_assert(!file.exists(path) && file.rename(tmp, path), "ATOMIC_PUBLICATION_FAILED")
        }
        fd_read <- function(p) jsonlite::fromJSON(p, simplifyVector = FALSE)
        fd_capture <- function(expr) {
            warnings <- character()
            value <- withCallingHandlers(expr, warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
            list(value = value, warnings = unique(warnings))
        }
        fd_ids <- function(model) as.character(unlist(if (!is.null(model$indicator_ids)) model$indicator_ids else model$members))
        fd_model <- function(cfg, id) {
            z <- Filter(function(m) identical(m$model_id, id), cfg$models)
            fd_assert(length(z) == 1L, "MODEL_NOT_UNIQUELY_FROZEN")
            z <- z[[1]]
            ids <- fd_ids(z)
            fd_assert(length(ids) > 0 && !anyDuplicated(ids) && setequal(ids, unlist(z$generation_order)), "INVALID_FROZEN_GENERATION_ORDER")
            z
        }
        fd_rng <- function(seed, code) {
            kind <- RNGkind()
            had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
            if (had) old <- get(".Random.seed", .GlobalEnv)
            on.exit({
                do.call(RNGkind, as.list(kind))
                if (had) assign(".Random.seed", old, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", 
                  envir = .GlobalEnv)
            })
            RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
            set.seed(seed)
            force(code)
        }
        fd_seed <- function(b, stream = 1L) as.integer(20261031L + as.integer(b) * 1009L + as.integer(stream) * 
            100003L)
        fd_outdir <- function(run, out) {
            run <- normalizePath(run, mustWork = TRUE)
            fd_assert(!grepl("(^|/)[.][.](/|$)", out), "UNSAFE_OUTPUT_PATH")
            fd_assert(startsWith(paste0(out, "/"), paste0(run, "/")), "OUTPUT_OUTSIDE_ISOLATED_RUN")
            if (!dir.exists(dirname(out))) dir.create(dirname(out), recursive = TRUE)
            fd_assert(startsWith(paste0(normalizePath(dirname(out)), "/"), paste0(run, "/")), "OUTPUT_SYMLINK_ESCAPES_RUN")
            if (dir.exists(out)) {
                allowed <- c("process.log", "PROCESS.json", "command.json", "COMMAND.json")
                fd_assert(all(list.files(out, all.files = TRUE, no.. = TRUE) %in% allowed), "OUTPUT_ATTEMPT_CONTAINS_PRIOR_SCIENTIFIC_ARTIFACT")
            } else fd_assert(dir.create(out), "OUTPUT_CREATE_FAILED")
            out
        }
    }), envir = context)
}

