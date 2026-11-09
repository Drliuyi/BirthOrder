# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

bi_assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)

bi_stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

bi_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

bi_hash <- function(x, file = FALSE) if (file) digest::digest(file = x, algo = "sha256") else digest::digest(x, 
    algo = "sha256")

bi_read_table <- function(path) {
    bi_assert(file.exists(path), paste("MISSING_INPUT", path))
    if (grepl("\\.rds$", path, ignore.case = TRUE)) 
        return(as.data.frame(readRDS(path)))
    as.data.frame(data.table::fread(path, check.names = FALSE, data.table = FALSE))
}

bi_bind <- function(xs) {
    xs <- Filter(function(x) !is.null(x) && NROW(x) > 0, xs)
    if (!length(xs)) 
        return(data.frame())
    as.data.frame(data.table::rbindlist(lapply(xs, as.data.frame), fill = TRUE, use.names = TRUE))
}

bi_capture <- function(expr) {
    w <- character()
    v <- tryCatch(withCallingHandlers(expr, warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
    }), error = identity)
    list(value = v, warnings = unique(w))
}

bi_log <- function(ctx, event, ...) {
    msg <- c(list(time = bi_stamp(), stage = ctx$stage, pid = Sys.getpid(), event = event), list(...))
    line <- jsonlite::toJSON(msg, auto_unbox = TRUE, na = "null", null = "null")
    if (!is.null(ctx$out) && dir.exists(file.path(ctx$out, "logs"))) 
        cat(line, "\n", file = file.path(ctx$out, "logs", paste0(ctx$stage, "_", Sys.getpid(), ".jsonl")), append = TRUE)
    message("[border_integrated] ", line)
}

bi_public <- function(x, k = 5L) {
    bi_assert(is.data.frame(x), "PUBLIC_TABLE_REQUIRED")
    forbidden <- intersect(tolower(names(x)), c("eid", "id", "pair_id", "birth_east", "birth_north", "home_east", 
        "home_north", "latitude", "longitude", "lon.birth", "lat.birth", "lon.home", "lat.home"))
    bi_assert(!length(forbidden), paste("PARTICIPANT_OR_LOCATION_EXPORT_BLOCKED", paste(forbidden, collapse = ",")))
    countcols <- names(x)[grepl("(^N$|^n$|cases|noncases|events|count|^n_|_n$|pairs)", names(x), ignore.case = TRUE)]
    suppress <- rep(FALSE, nrow(x))
    for (n in countcols) if (is.numeric(x[[n]]) || is.integer(x[[n]])) {
        small <- is.finite(x[[n]]) & x[[n]] > 0 & x[[n]] < k
        suppress <- suppress | small
        x[[n]][small] <- NA_real_
    }
    if (any(suppress)) {
        statscols <- names(x)[grepl("(^|_)(beta|se|OR|HR|CI|p|q)(_|$)|percent|proportion|difference|delta", names(x))]
        for (n in union(countcols, statscols)) x[[n]][suppress] <- NA
        x$disclosure_suppressed <- suppress
    }
    x
}

bi_save <- function(ctx, relative_path, x) {
    bi_assert(!grepl("(^/|(^|/)\\.\\.(/|$))", relative_path), "INVALID_OUTPUT_PATH")
    path <- file.path(ctx$out, relative_path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (is.data.frame(x)) 
        x <- bi_public(x, ctx$cfg$disclosure$minimum_cell %||% 5L)
    tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
    if (is.data.frame(x)) 
        write.csv(x, tmp, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
    else writeLines(as.character(x), tmp, useBytes = TRUE)
    if (file.exists(path)) {
        if (identical(bi_hash(path, TRUE), bi_hash(tmp, TRUE))) {
            unlink(tmp)
            return(invisible(path))
        }
        unlink(tmp)
        stop("REFUSE_DIFFERENT_EXISTING_OUTPUT: ", path)
    }
    bi_assert(file.rename(tmp, path), "OUTPUT_PUBLISH_FAILED")
    invisible(path)
}

bi_checkpoint <- function(ctx, key, fn) {
    bi_assert(grepl("^[A-Za-z0-9_.-]+$", key), "INVALID_CHECKPOINT_KEY")
    path <- file.path(ctx$out, "cache", ctx$stage, paste0(key, ".rds"))
    if (file.exists(path)) {
        bi_assert(isTRUE(ctx$resume), "EXISTING_CHECKPOINT_REQUIRES_RESUME")
        z <- readRDS(path)
        h <- z$hash
        z$hash <- NULL
        bi_assert(identical(h, bi_hash(z)) && identical(z$signature, ctx$signature), "CHECKPOINT_SIGNATURE_MISMATCH")
        bi_log(ctx, "RESUMED", key = key)
        return(z$value)
    }
    bi_resource_guard(ctx)
    bi_log(ctx, "START_MODEL", key = key)
    value <- fn()
    z <- list(signature = ctx$signature, value = value)
    z$hash <- bi_hash(z)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
    saveRDS(z, tmp)
    Sys.chmod(tmp, "0600")
    bi_assert(file.rename(tmp, path), "CHECKPOINT_PUBLISH_FAILED")
    bi_log(ctx, "SAVED_MODEL", key = key, status = if (is.data.frame(value) && "status" %in% names(value)) 
        paste(unique(value$status), collapse = ";")
    else "SAVED")
    value
}

bi_resource_guard <- function(ctx) {
    if (!is.null(ctx$started) && is.finite(ctx$cfg$resources$module_timeout_hours %||% Inf)) 
        bi_assert(as.numeric(difftime(Sys.time(), ctx$started, units = "hours")) < ctx$cfg$resources$module_timeout_hours, 
            "MODULE_TIME_LIMIT_REACHED: resume with unchanged configuration")
    if (file.exists("/proc/meminfo")) {
        l <- readLines("/proc/meminfo")
        a <- as.numeric(gsub("[^0-9]", "", l[grepl("^MemAvailable:", l)]))/1024^2
        bi_assert(length(a) == 1 && a >= (ctx$cfg$resources$minimum_free_gb %||% 2), "INSUFFICIENT_FREE_MEMORY")
    }
}

bi_skip <- function(ctx, reason) {
    bi_log(ctx, "SKIPPED", reason = reason)
    list(status = "SKIPPED_UNAVAILABLE", reason = reason)
}

bi_module_config <- function(ctx, key) ctx$cfg$modules[[key]] %||% list()

bi_formula <- function(exposure, extra = character(), omit = character(), response = "y") {
    bi_assert(exposure %in% linear_exposures, "INVALID_EXPOSURE")
    reformulate(setdiff(unique(c(exposure, linear_covariates, if (exposure != "border12") "factor(S)", extra)), 
        omit), response = response)
}

bi_sample <- function(d, formula, exposure) {
    bi_assert("eid" %in% names(d) && !anyNA(d$eid) && !anyDuplicated(as.character(d$eid)), "DUPLICATE_OR_MISSING_ID")
    required <- unique(c(all.vars(formula), exposure))
    missing <- setdiff(required, names(d))
    bi_assert(!length(missing), paste("MISSING_MODEL_FIELDS", paste(missing, collapse = ",")))
    keep <- complete.cases(d[, required, drop = FALSE])
    for (k in required) if (is.numeric(d[[k]])) 
        keep <- keep & is.finite(d[[k]])
    d <- droplevels(as.data.frame(d[keep, , drop = FALSE]))
    terms <- attr(terms(formula), "term.labels")
    remove <- character()
    for (term in setdiff(terms, exposure)) {
        if (grepl(":|\\*", term)) 
            next
        vars <- all.vars(as.formula(paste("~", term)))
        if (length(vars) == 1 && length(unique(d[[vars]])) <= 1) 
            remove <- c(remove, term)
    }
    if (length(remove)) 
        formula <- reformulate(setdiff(terms, remove), response = paste(deparse(formula[[2]]), collapse = ""))
    list(data = d, formula = formula, removed_terms = remove, counts = list(N = nrow(d), cases = if ("y" %in% names(d)) sum(d$y == 
        1, na.rm = TRUE) else NA_integer_, noncases = if ("y" %in% names(d)) sum(d$y == 0, na.rm = TRUE) else NA_integer_))
}

bi_fit <- function(ctx, d, exposure, extra = character(), omit = character(), family = "binomial", formula = NULL) {
    f <- formula %||% bi_formula(exposure, extra, omit)
    r <- data.frame(exposure = exposure, N = NA_integer_, cases = NA_integer_, noncases = NA_integer_, beta = NA_real_, 
        se = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_, p_raw = NA_real_, status = "NOT_FITTED", 
        warning = "", error = "", formula = paste(deparse(f), collapse = " "), rank = NA_integer_, columns = NA_integer_, 
        sample_hash = "", removed_terms = "", separation_status = "NOT_TESTED", inference = "DESCRIPTIVE_PENDING_REVIEW", 
        stringsAsFactors = FALSE)
    tryCatch({
        if (family == "binomial") 
            bi_assert(all(is.na(d$y) | d$y %in% c(0, 1)), "UNKNOWN_OUTCOME_STATE")
        s <- bi_sample(d, f, exposure)
        dd <- s$data
        f <- s$formula
        r$N <- nrow(dd)
        r$formula <- paste(deparse(f), collapse = " ")
        r$removed_terms <- paste(s$removed_terms, collapse = ";")
        r$sample_hash <- bi_hash(sort(as.character(dd$eid), method = "radix"))
        if (family == "binomial") {
            r$cases <- sum(dd$y == 1)
            r$noncases <- sum(dd$y == 0)
            bi_assert(r$cases >= (ctx$cfg$thresholds$min_cases %||% 20) && r$noncases >= (ctx$cfg$thresholds$min_controls %||% 
                20), "INSUFFICIENT_CASE_CONTROL_SUPPORT")
        }
        else bi_assert(nrow(dd) >= (ctx$cfg$thresholds$min_continuous_n %||% 100), "INSUFFICIENT_CONTINUOUS_N")
        bi_assert(length(unique(dd[[exposure]])) > 1, "CONSTANT_EXPOSURE")
        x <- model.matrix(f, dd)
        r$rank <- qr(x)$rank
        r$columns <- ncol(x)
        bi_assert(r$rank == r$columns && all(is.finite(x)), "RANK_DEFICIENT_OR_NONFINITE_DESIGN")
        bi_assert(exposure %in% colnames(x), "TARGET_COLUMN_MISSING")
        z <- bi_capture(if (family == "binomial") 
            glm(f, data = dd, family = binomial(), control = glm.control(maxit = 25), na.action = na.fail)
        else lm(f, data = dd, na.action = na.fail))
        r$warning <- paste(z$warnings, collapse = "; ")
        if (inherits(z$value, "error")) 
            stop(conditionMessage(z$value))
        fit <- z$value
        cf <- coef(summary(fit))
        r$p_raw <- cf[exposure, ncol(cf)]
        b <- cf[exposure, 1]
        se <- cf[exposure, 2]
        r$raw_beta <- b
        r$raw_se <- se
        bi_assert(is.finite(b) && is.finite(se) && se > 0 && is.finite(r$p_raw), "INVALID_TARGET_ESTIMATE")
        if (family == "binomial") 
            bi_assert(isTRUE(fit$converged), "GLM_NOT_CONVERGED")
        r$beta <- b
        r$se <- se
        r$p <- r$p_raw
        crit <- if (family == "binomial") 
            qnorm(0.975)
        else qt(0.975, fit$df.residual)
        r$CI_low <- b - crit * se
        r$CI_high <- b + crit * se
        if (family == "binomial") {
            r$OR <- exp(b)
            r$CI_low <- exp(r$CI_low)
            r$CI_high <- exp(r$CI_high)
            bi_assert(all(is.finite(c(r$OR, r$CI_low, r$CI_high))) && r$CI_low > 0, "NONFINITE_OR_CI")
        }
        r$status <- if (length(z$warnings)) 
            "FIT_WARNING_REVIEW"
        else "OK"
        if (family != "binomial") 
            r$inference <- "ASSOCIATION_NONCAUSAL"
    }, error = function(e) {
        r$error <<- conditionMessage(e)
        r$status <<- "NOT_ESTIMABLE"
        r$beta <<- r$se <<- r$OR <<- r$CI_low <<- r$CI_high <<- r$p <<- NA_real_
    })
    r
}

bi_interaction <- function(ctx, d, exposure, modifier, extra = character(), omit = character()) {
    r <- data.frame(exposure = exposure, modifier = modifier, N = NA_integer_, cases = NA_integer_, p = NA_real_, 
        p_raw = NA_real_, df = NA_integer_, status = "NOT_ESTIMABLE", error = "", warning = "", formula = "", sample_hash = "", 
        inference = "DESCRIPTIVE_PENDING_REVIEW")
    tryCatch({
        bi_assert(modifier %in% names(d), "MISSING_MODIFIER")
        nestedS <- modifier %in% c("S_group", "family_size_group", "sibship_group")
        f0 <- bi_formula(exposure, c(extra, if (!nestedS) modifier), omit)
        d <- d[complete.cases(d[, unique(c(all.vars(f0), modifier)), drop = FALSE]), , drop = FALSE]
        s <- bi_sample(d, f0, exposure)
        d <- s$data
        f0 <- s$formula
        bi_assert(all(is.na(d$y) | d$y %in% c(0, 1)), "UNKNOWN_OUTCOME_STATE")
        r$N <- nrow(d)
        r$cases <- sum(d$y == 1)
        r$sample_hash <- bi_hash(sort(as.character(d$eid), method = "radix"))
        bi_assert(r$cases >= (ctx$cfg$thresholds$min_cases %||% 20) && sum(d$y == 0) >= (ctx$cfg$thresholds$min_controls %||% 
            20), "INSUFFICIENT_SUPPORT")
        bi_assert(length(unique(d[[modifier]])) > 1 && length(unique(d[[exposure]])) > 1, "CONSTANT_MODIFIER_OR_EXPOSURE")
        m <- if (is.factor(d[[modifier]]) || is.character(d[[modifier]])) 
            model.matrix(~z, data.frame(z = factor(d[[modifier]])))[, -1, drop = FALSE]
        else matrix(as.numeric(d[[modifier]]), ncol = 1)
        intnames <- paste0("bi_int_", seq_len(ncol(m)))
        for (j in seq_along(intnames)) d[[intnames[j]]] <- d[[exposure]] * m[, j]
        f1 <- reformulate(c(attr(terms(f0), "term.labels"), intnames), response = "y")
        x0 <- model.matrix(f0, d)
        x1 <- model.matrix(f1, d)
        bi_assert(qr(x0)$rank == ncol(x0) && qr(x1)$rank == ncol(x1), "INTERACTION_RANK_DEFICIENT")
        z <- bi_capture({
            a <- glm(f0, d, family = binomial(), na.action = na.fail)
            b <- glm(f1, d, family = binomial(), na.action = na.fail)
            list(a = a, b = b)
        })
        r$warning <- paste(z$warnings, collapse = ";")
        if (inherits(z$value, "error")) 
            stop(conditionMessage(z$value))
        a <- z$value$a
        b <- z$value$b
        bi_assert(a$converged && b$converged, "INTERACTION_NONCONVERGENCE")
        bi_assert(all(is.finite(coef(a))) && all(is.finite(coef(b))), "INTERACTION_NONFINITE_COEFFICIENT")
        r$df <- length(coef(b)) - length(coef(a))
        r$p_raw <- pchisq(max(0, 2 * (as.numeric(logLik(b)) - as.numeric(logLik(a)))), r$df, lower.tail = FALSE)
        r$p <- r$p_raw
        r$formula <- paste(deparse(f1), collapse = " ")
        r$status <- if (length(z$warnings)) 
            "FIT_WARNING_REVIEW"
        else "OK"
    }, error = function(e) {
        r$error <<- conditionMessage(e)
    })
    r
}

bi_adjust <- function(ctx, table, family_col = "family_id", p_col = "p", id_col = "test_id") {
    if (!nrow(table)) 
        return(table)
    bi_assert(all(c(family_col, p_col, id_col) %in% names(table)), "MISSING_MULTIPLICITY_REGISTRY")
    bi_assert(!anyDuplicated(paste(table[[family_col]], table[[id_col]], sep = "::")), "DUPLICATE_FAMILY_TEST")
    table$q <- NA_real_
    table$M <- NA_integer_
    regs <- list()
    for (f in unique(table[[family_col]])) {
        ix <- which(table[[family_col]] == f)
        m <- length(ix)
        p <- as.numeric(table[[p_col]][ix])
        ok <- is.finite(p) & p >= 0 & p <= 1
        table$M[ix] <- m
        if (any(ok)) 
            table$q[ix[ok]] <- p.adjust(p[ok], "BH", n = m)
        regs[[f]] <- data.frame(family_id = f, M = m, valid_p = sum(ok), missing_p = sum(!ok), correction = "BH; planned slots retained", 
            scope = ctx$stage)
    }
    ctx$registries <- c(ctx$registries %||% list(), regs)
    table
}

bi_education <- function(x) {
    tokens <- strsplit(gsub("[,; ]+", "|", as.character(x)), "|", fixed = TRUE)
    hierarchy <- c("1", "5", "6", "2", "3", "4", "-7")
    top <- vapply(tokens, function(v) {
        v <- v[nzchar(v)]
        if (!length(v) || any(!v %in% c("1", "2", "3", "4", "5", "6", "-7", "-3")) || "-3" %in% v || ("-7" %in% 
            v && length(v) > 1)) 
            return(NA_character_)
        z <- hierarchy[hierarchy %in% v]
        if (length(z)) 
            z[1]
        else NA_character_
    }, character(1))
    list(education = factor(top, levels = rev(hierarchy), labels = c("None", "CSE", "O_GCSE", "A_AS", "Professional", 
        "NVQ_HND_HNC", "Degree")), university = factor(ifelse(is.na(top), NA, as.integer(top == "1")), levels = 0:1))
}

bi_cohort <- function(ctx) {
    if (exists("dat", envir = ctx, inherits = FALSE)) 
        return(ctx$dat)
    bi_assert(isTRUE(ctx$real_authorized) || isTRUE(ctx$synthetic), "REAL_DATA_ACCESS_NOT_AUTHORIZED")
    b <- readRDS(ctx$cfg$inputs$all_rds)
    d <- linear_covars(make_birth_order(as.data.frame(b)))
    need <- c("eid", linear_exposures, linear_covariates, "S", "birth_year", "cnt_brother", "cnt_sister")
    bi_assert(all(need %in% names(d)), "CORE_COHORT_FIELDS_MISSING")
    get <- function(key) {
        col <- ctx$cfg$field_map[[key]]
        if (is.null(col) || !col %in% names(d)) 
            rep(NA, nrow(d))
        else d[[col]]
    }
    d$birth_year <- bi_num(d$birth_year)
    d$tdi <- bi_num(get("tdi"))
    d$bmi <- bi_num(get("bmi"))
    if ("edu.sco" %in% names(d)) 
        d$edu.sco <- bi_num(d$edu.sco)
    income <- bi_num(get("income"))
    income[!income %in% 1:5] <- NA
    d$income <- factor(income, levels = 1:5)
    d$income_group <- factor(ifelse(income %in% c(1, 2), "low", ifelse(income == 3, "middle", ifelse(income %in% 
        c(4, 5), "high", NA))), levels = c("low", "middle", "high"))
    ed <- bi_education(get("education"))
    d$education <- ed$education
    d$university <- ed$university
    smoking <- bi_num(get("smoking"))
    smoking[!smoking %in% 0:2] <- NA
    d$smoking <- factor(smoking, levels = 0:2)
    alcohol <- bi_num(get("alcohol"))
    alcohol[!alcohol %in% 0:2] <- NA
    d$alcohol <- factor(alcohol, levels = 0:2)
    d$family_composition <- factor(ifelse(d$cnt_brother > 0 & d$cnt_sister == 0, "Brothers_only", ifelse(d$cnt_sister > 
        0 & d$cnt_brother == 0, "Sisters_only", ifelse(d$cnt_brother > 0 & d$cnt_sister > 0, "Mixed", NA))))
    d$sibling_sex <- factor(ifelse(d$S == 2, as.character(d$family_composition), NA), levels = c("Brothers_only", 
        "Sisters_only"))
    ctx$dat <- d
    bi_log(ctx, "COHORT_READY", rows = nrow(d))
    d
}

bi_outcomes <- function(ctx) {
    if (exists("ph", envir = ctx, inherits = FALSE)) 
        return(ctx$ph)
    d <- bi_cohort(ctx)
    p <- as.data.frame(readRDS(ctx$cfg$inputs$phewas_rds))
    bi_assert("id" %in% names(p) && !anyNA(p$id) && !anyDuplicated(as.character(p$id)), "PHENOTYPE_ID_INVALID")
    codes <- unique(c(unlist(ctx$cfg$clinical$core), unlist(ctx$cfg$clinical$supplementary)))
    missing_codes <- setdiff(unlist(ctx$cfg$clinical$core), names(p))
    if (length(missing_codes)) 
        bi_log(ctx, "OUTCOME_COLUMNS_UNAVAILABLE", phecodes = paste(missing_codes, collapse = ";"), action = "Keep registered outcomes as NA; never code unknown as control")
    ix <- match(as.character(d$eid), as.character(p$id))
    out <- data.frame(eid = d$eid)
    for (k in codes) {
        y <- if (k %in% names(p)) 
            p[[k]][ix]
        else rep(NA, nrow(d))
        bi_assert(all(is.na(y) | y %in% c(0, 1)), paste("INVALID_OUTCOME", k))
        out[[k]] <- as.numeric(y)
    }
    ctx$ph <- out
    out
}

bi_data <- function(ctx, code) {
    d <- bi_cohort(ctx)
    p <- bi_outcomes(ctx)
    bi_assert(code %in% names(p), "OUTCOME_NOT_REGISTERED")
    bi_assert(identical(as.character(d$eid), as.character(p$eid)), "OUTCOME_ALIGNMENT_FAILURE")
    d$y <- p[[code]]
    d
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        `%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
        bi_assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
        bi_stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
        bi_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
        bi_hash <- function(x, file = FALSE) if (file) digest::digest(file = x, algo = "sha256") else digest::digest(x, 
            algo = "sha256")
        bi_read_table <- function(path) {
            bi_assert(file.exists(path), paste("MISSING_INPUT", path))
            if (grepl("\\.rds$", path, ignore.case = TRUE)) return(as.data.frame(readRDS(path)))
            as.data.frame(data.table::fread(path, check.names = FALSE, data.table = FALSE))
        }
        bi_bind <- function(xs) {
            xs <- Filter(function(x) !is.null(x) && NROW(x) > 0, xs)
            if (!length(xs)) return(data.frame())
            as.data.frame(data.table::rbindlist(lapply(xs, as.data.frame), fill = TRUE, use.names = TRUE))
        }
        bi_capture <- function(expr) {
            w <- character()
            v <- tryCatch(withCallingHandlers(expr, warning = function(x) {
                w <<- c(w, conditionMessage(x))
                invokeRestart("muffleWarning")
            }), error = identity)
            list(value = v, warnings = unique(w))
        }
        bi_log <- function(ctx, event, ...) {
            msg <- c(list(time = bi_stamp(), stage = ctx$stage, pid = Sys.getpid(), event = event), list(...))
            line <- jsonlite::toJSON(msg, auto_unbox = TRUE, na = "null", null = "null")
            if (!is.null(ctx$out) && dir.exists(file.path(ctx$out, "logs"))) cat(line, "\n", file = file.path(ctx$out, 
                "logs", paste0(ctx$stage, "_", Sys.getpid(), ".jsonl")), append = TRUE)
            message("[border_integrated] ", line)
        }
        bi_public <- function(x, k = 5L) {
            bi_assert(is.data.frame(x), "PUBLIC_TABLE_REQUIRED")
            forbidden <- intersect(tolower(names(x)), c("eid", "id", "pair_id", "birth_east", "birth_north", "home_east", 
                "home_north", "latitude", "longitude", "lon.birth", "lat.birth", "lon.home", "lat.home"))
            bi_assert(!length(forbidden), paste("PARTICIPANT_OR_LOCATION_EXPORT_BLOCKED", paste(forbidden, collapse = ",")))
            countcols <- names(x)[grepl("(^N$|^n$|cases|noncases|events|count|^n_|_n$|pairs)", names(x), ignore.case = TRUE)]
            suppress <- rep(FALSE, nrow(x))
            for (n in countcols) if (is.numeric(x[[n]]) || is.integer(x[[n]])) {
                small <- is.finite(x[[n]]) & x[[n]] > 0 & x[[n]] < k
                suppress <- suppress | small
                x[[n]][small] <- NA_real_
            }
            if (any(suppress)) {
                statscols <- names(x)[grepl("(^|_)(beta|se|OR|HR|CI|p|q)(_|$)|percent|proportion|difference|delta", 
                  names(x))]
                for (n in union(countcols, statscols)) x[[n]][suppress] <- NA
                x$disclosure_suppressed <- suppress
            }
            x
        }
        bi_save <- function(ctx, relative_path, x) {
            bi_assert(!grepl("(^/|(^|/)\\.\\.(/|$))", relative_path), "INVALID_OUTPUT_PATH")
            path <- file.path(ctx$out, relative_path)
            dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
            if (is.data.frame(x)) x <- bi_public(x, ctx$cfg$disclosure$minimum_cell %||% 5L)
            tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
            if (is.data.frame(x)) write.csv(x, tmp, row.names = FALSE, na = "NA", fileEncoding = "UTF-8") else writeLines(as.character(x), 
                tmp, useBytes = TRUE)
            if (file.exists(path)) {
                if (identical(bi_hash(path, TRUE), bi_hash(tmp, TRUE))) {
                  unlink(tmp)
                  return(invisible(path))
                }
                unlink(tmp)
                stop("REFUSE_DIFFERENT_EXISTING_OUTPUT: ", path)
            }
            bi_assert(file.rename(tmp, path), "OUTPUT_PUBLISH_FAILED")
            invisible(path)
        }
        bi_checkpoint <- function(ctx, key, fn) {
            bi_assert(grepl("^[A-Za-z0-9_.-]+$", key), "INVALID_CHECKPOINT_KEY")
            path <- file.path(ctx$out, "cache", ctx$stage, paste0(key, ".rds"))
            if (file.exists(path)) {
                bi_assert(isTRUE(ctx$resume), "EXISTING_CHECKPOINT_REQUIRES_RESUME")
                z <- readRDS(path)
                h <- z$hash
                z$hash <- NULL
                bi_assert(identical(h, bi_hash(z)) && identical(z$signature, ctx$signature), "CHECKPOINT_SIGNATURE_MISMATCH")
                bi_log(ctx, "RESUMED", key = key)
                return(z$value)
            }
            bi_resource_guard(ctx)
            bi_log(ctx, "START_MODEL", key = key)
            value <- fn()
            z <- list(signature = ctx$signature, value = value)
            z$hash <- bi_hash(z)
            dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
            tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
            saveRDS(z, tmp)
            Sys.chmod(tmp, "0600")
            bi_assert(file.rename(tmp, path), "CHECKPOINT_PUBLISH_FAILED")
            bi_log(ctx, "SAVED_MODEL", key = key, status = if (is.data.frame(value) && "status" %in% names(value)) paste(unique(value$status), 
                collapse = ";") else "SAVED")
            value
        }
        bi_resource_guard <- function(ctx) {
            if (!is.null(ctx$started) && is.finite(ctx$cfg$resources$module_timeout_hours %||% Inf)) bi_assert(as.numeric(difftime(Sys.time(), 
                ctx$started, units = "hours")) < ctx$cfg$resources$module_timeout_hours, "MODULE_TIME_LIMIT_REACHED: resume with unchanged configuration")
            if (file.exists("/proc/meminfo")) {
                l <- readLines("/proc/meminfo")
                a <- as.numeric(gsub("[^0-9]", "", l[grepl("^MemAvailable:", l)]))/1024^2
                bi_assert(length(a) == 1 && a >= (ctx$cfg$resources$minimum_free_gb %||% 2), "INSUFFICIENT_FREE_MEMORY")
            }
        }
        bi_skip <- function(ctx, reason) {
            bi_log(ctx, "SKIPPED", reason = reason)
            list(status = "SKIPPED_UNAVAILABLE", reason = reason)
        }
        bi_module_config <- function(ctx, key) ctx$cfg$modules[[key]] %||% list()
        bi_formula <- function(exposure, extra = character(), omit = character(), response = "y") {
            bi_assert(exposure %in% linear_exposures, "INVALID_EXPOSURE")
            reformulate(setdiff(unique(c(exposure, linear_covariates, if (exposure != "border12") "factor(S)", 
                extra)), omit), response = response)
        }
        bi_sample <- function(d, formula, exposure) {
            bi_assert("eid" %in% names(d) && !anyNA(d$eid) && !anyDuplicated(as.character(d$eid)), "DUPLICATE_OR_MISSING_ID")
            required <- unique(c(all.vars(formula), exposure))
            missing <- setdiff(required, names(d))
            bi_assert(!length(missing), paste("MISSING_MODEL_FIELDS", paste(missing, collapse = ",")))
            keep <- complete.cases(d[, required, drop = FALSE])
            for (k in required) if (is.numeric(d[[k]])) keep <- keep & is.finite(d[[k]])
            d <- droplevels(as.data.frame(d[keep, , drop = FALSE]))
            terms <- attr(terms(formula), "term.labels")
            remove <- character()
            for (term in setdiff(terms, exposure)) {
                if (grepl(":|\\*", term)) next
                vars <- all.vars(as.formula(paste("~", term)))
                if (length(vars) == 1 && length(unique(d[[vars]])) <= 1) remove <- c(remove, term)
            }
            if (length(remove)) formula <- reformulate(setdiff(terms, remove), response = paste(deparse(formula[[2]]), 
                collapse = ""))
            list(data = d, formula = formula, removed_terms = remove, counts = list(N = nrow(d), cases = if ("y" %in% 
                names(d)) sum(d$y == 1, na.rm = TRUE) else NA_integer_, noncases = if ("y" %in% names(d)) sum(d$y == 
                0, na.rm = TRUE) else NA_integer_))
        }
        bi_fit <- function(ctx, d, exposure, extra = character(), omit = character(), family = "binomial", formula = NULL) {
            f <- formula %||% bi_formula(exposure, extra, omit)
            r <- data.frame(exposure = exposure, N = NA_integer_, cases = NA_integer_, noncases = NA_integer_, 
                beta = NA_real_, se = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_, 
                p_raw = NA_real_, status = "NOT_FITTED", warning = "", error = "", formula = paste(deparse(f), 
                  collapse = " "), rank = NA_integer_, columns = NA_integer_, sample_hash = "", removed_terms = "", 
                separation_status = "NOT_TESTED", inference = "DESCRIPTIVE_PENDING_REVIEW", stringsAsFactors = FALSE)
            tryCatch({
                if (family == "binomial") bi_assert(all(is.na(d$y) | d$y %in% c(0, 1)), "UNKNOWN_OUTCOME_STATE")
                s <- bi_sample(d, f, exposure)
                dd <- s$data
                f <- s$formula
                r$N <- nrow(dd)
                r$formula <- paste(deparse(f), collapse = " ")
                r$removed_terms <- paste(s$removed_terms, collapse = ";")
                r$sample_hash <- bi_hash(sort(as.character(dd$eid), method = "radix"))
                if (family == "binomial") {
                  r$cases <- sum(dd$y == 1)
                  r$noncases <- sum(dd$y == 0)
                  bi_assert(r$cases >= (ctx$cfg$thresholds$min_cases %||% 20) && r$noncases >= (ctx$cfg$thresholds$min_controls %||% 
                    20), "INSUFFICIENT_CASE_CONTROL_SUPPORT")
                } else bi_assert(nrow(dd) >= (ctx$cfg$thresholds$min_continuous_n %||% 100), "INSUFFICIENT_CONTINUOUS_N")
                bi_assert(length(unique(dd[[exposure]])) > 1, "CONSTANT_EXPOSURE")
                x <- model.matrix(f, dd)
                r$rank <- qr(x)$rank
                r$columns <- ncol(x)
                bi_assert(r$rank == r$columns && all(is.finite(x)), "RANK_DEFICIENT_OR_NONFINITE_DESIGN")
                bi_assert(exposure %in% colnames(x), "TARGET_COLUMN_MISSING")
                z <- bi_capture(if (family == "binomial") glm(f, data = dd, family = binomial(), control = glm.control(maxit = 25), 
                  na.action = na.fail) else lm(f, data = dd, na.action = na.fail))
                r$warning <- paste(z$warnings, collapse = "; ")
                if (inherits(z$value, "error")) stop(conditionMessage(z$value))
                fit <- z$value
                cf <- coef(summary(fit))
                r$p_raw <- cf[exposure, ncol(cf)]
                b <- cf[exposure, 1]
                se <- cf[exposure, 2]
                r$raw_beta <- b
                r$raw_se <- se
                bi_assert(is.finite(b) && is.finite(se) && se > 0 && is.finite(r$p_raw), "INVALID_TARGET_ESTIMATE")
                if (family == "binomial") bi_assert(isTRUE(fit$converged), "GLM_NOT_CONVERGED")
                r$beta <- b
                r$se <- se
                r$p <- r$p_raw
                crit <- if (family == "binomial") qnorm(0.975) else qt(0.975, fit$df.residual)
                r$CI_low <- b - crit * se
                r$CI_high <- b + crit * se
                if (family == "binomial") {
                  r$OR <- exp(b)
                  r$CI_low <- exp(r$CI_low)
                  r$CI_high <- exp(r$CI_high)
                  bi_assert(all(is.finite(c(r$OR, r$CI_low, r$CI_high))) && r$CI_low > 0, "NONFINITE_OR_CI")
                }
                r$status <- if (length(z$warnings)) "FIT_WARNING_REVIEW" else "OK"
                if (family != "binomial") r$inference <- "ASSOCIATION_NONCAUSAL"
            }, error = function(e) {
                r$error <<- conditionMessage(e)
                r$status <<- "NOT_ESTIMABLE"
                r$beta <<- r$se <<- r$OR <<- r$CI_low <<- r$CI_high <<- r$p <<- NA_real_
            })
            r
        }
        bi_interaction <- function(ctx, d, exposure, modifier, extra = character(), omit = character()) {
            r <- data.frame(exposure = exposure, modifier = modifier, N = NA_integer_, cases = NA_integer_, p = NA_real_, 
                p_raw = NA_real_, df = NA_integer_, status = "NOT_ESTIMABLE", error = "", warning = "", formula = "", 
                sample_hash = "", inference = "DESCRIPTIVE_PENDING_REVIEW")
            tryCatch({
                bi_assert(modifier %in% names(d), "MISSING_MODIFIER")
                nestedS <- modifier %in% c("S_group", "family_size_group", "sibship_group")
                f0 <- bi_formula(exposure, c(extra, if (!nestedS) modifier), omit)
                d <- d[complete.cases(d[, unique(c(all.vars(f0), modifier)), drop = FALSE]), , drop = FALSE]
                s <- bi_sample(d, f0, exposure)
                d <- s$data
                f0 <- s$formula
                bi_assert(all(is.na(d$y) | d$y %in% c(0, 1)), "UNKNOWN_OUTCOME_STATE")
                r$N <- nrow(d)
                r$cases <- sum(d$y == 1)
                r$sample_hash <- bi_hash(sort(as.character(d$eid), method = "radix"))
                bi_assert(r$cases >= (ctx$cfg$thresholds$min_cases %||% 20) && sum(d$y == 0) >= (ctx$cfg$thresholds$min_controls %||% 
                  20), "INSUFFICIENT_SUPPORT")
                bi_assert(length(unique(d[[modifier]])) > 1 && length(unique(d[[exposure]])) > 1, "CONSTANT_MODIFIER_OR_EXPOSURE")
                m <- if (is.factor(d[[modifier]]) || is.character(d[[modifier]])) model.matrix(~z, data.frame(z = factor(d[[modifier]])))[, 
                  -1, drop = FALSE] else matrix(as.numeric(d[[modifier]]), ncol = 1)
                intnames <- paste0("bi_int_", seq_len(ncol(m)))
                for (j in seq_along(intnames)) d[[intnames[j]]] <- d[[exposure]] * m[, j]
                f1 <- reformulate(c(attr(terms(f0), "term.labels"), intnames), response = "y")
                x0 <- model.matrix(f0, d)
                x1 <- model.matrix(f1, d)
                bi_assert(qr(x0)$rank == ncol(x0) && qr(x1)$rank == ncol(x1), "INTERACTION_RANK_DEFICIENT")
                z <- bi_capture({
                  a <- glm(f0, d, family = binomial(), na.action = na.fail)
                  b <- glm(f1, d, family = binomial(), na.action = na.fail)
                  list(a = a, b = b)
                })
                r$warning <- paste(z$warnings, collapse = ";")
                if (inherits(z$value, "error")) stop(conditionMessage(z$value))
                a <- z$value$a
                b <- z$value$b
                bi_assert(a$converged && b$converged, "INTERACTION_NONCONVERGENCE")
                bi_assert(all(is.finite(coef(a))) && all(is.finite(coef(b))), "INTERACTION_NONFINITE_COEFFICIENT")
                r$df <- length(coef(b)) - length(coef(a))
                r$p_raw <- pchisq(max(0, 2 * (as.numeric(logLik(b)) - as.numeric(logLik(a)))), r$df, lower.tail = FALSE)
                r$p <- r$p_raw
                r$formula <- paste(deparse(f1), collapse = " ")
                r$status <- if (length(z$warnings)) "FIT_WARNING_REVIEW" else "OK"
            }, error = function(e) {
                r$error <<- conditionMessage(e)
            })
            r
        }
        bi_adjust <- function(ctx, table, family_col = "family_id", p_col = "p", id_col = "test_id") {
            if (!nrow(table)) return(table)
            bi_assert(all(c(family_col, p_col, id_col) %in% names(table)), "MISSING_MULTIPLICITY_REGISTRY")
            bi_assert(!anyDuplicated(paste(table[[family_col]], table[[id_col]], sep = "::")), "DUPLICATE_FAMILY_TEST")
            table$q <- NA_real_
            table$M <- NA_integer_
            regs <- list()
            for (f in unique(table[[family_col]])) {
                ix <- which(table[[family_col]] == f)
                m <- length(ix)
                p <- as.numeric(table[[p_col]][ix])
                ok <- is.finite(p) & p >= 0 & p <= 1
                table$M[ix] <- m
                if (any(ok)) table$q[ix[ok]] <- p.adjust(p[ok], "BH", n = m)
                regs[[f]] <- data.frame(family_id = f, M = m, valid_p = sum(ok), missing_p = sum(!ok), correction = "BH; planned slots retained", 
                  scope = ctx$stage)
            }
            ctx$registries <- c(ctx$registries %||% list(), regs)
            table
        }
        bi_education <- function(x) {
            tokens <- strsplit(gsub("[,; ]+", "|", as.character(x)), "|", fixed = TRUE)
            hierarchy <- c("1", "5", "6", "2", "3", "4", "-7")
            top <- vapply(tokens, function(v) {
                v <- v[nzchar(v)]
                if (!length(v) || any(!v %in% c("1", "2", "3", "4", "5", "6", "-7", "-3")) || "-3" %in% v || ("-7" %in% 
                  v && length(v) > 1)) return(NA_character_)
                z <- hierarchy[hierarchy %in% v]
                if (length(z)) z[1] else NA_character_
            }, character(1))
            list(education = factor(top, levels = rev(hierarchy), labels = c("None", "CSE", "O_GCSE", "A_AS", "Professional", 
                "NVQ_HND_HNC", "Degree")), university = factor(ifelse(is.na(top), NA, as.integer(top == "1")), 
                levels = 0:1))
        }
        bi_cohort <- function(ctx) {
            if (exists("dat", envir = ctx, inherits = FALSE)) return(ctx$dat)
            bi_assert(isTRUE(ctx$real_authorized) || isTRUE(ctx$synthetic), "REAL_DATA_ACCESS_NOT_AUTHORIZED")
            b <- readRDS(ctx$cfg$inputs$all_rds)
            d <- linear_covars(make_birth_order(as.data.frame(b)))
            need <- c("eid", linear_exposures, linear_covariates, "S", "birth_year", "cnt_brother", "cnt_sister")
            bi_assert(all(need %in% names(d)), "CORE_COHORT_FIELDS_MISSING")
            get <- function(key) {
                col <- ctx$cfg$field_map[[key]]
                if (is.null(col) || !col %in% names(d)) rep(NA, nrow(d)) else d[[col]]
            }
            d$birth_year <- bi_num(d$birth_year)
            d$tdi <- bi_num(get("tdi"))
            d$bmi <- bi_num(get("bmi"))
            if ("edu.sco" %in% names(d)) d$edu.sco <- bi_num(d$edu.sco)
            income <- bi_num(get("income"))
            income[!income %in% 1:5] <- NA
            d$income <- factor(income, levels = 1:5)
            d$income_group <- factor(ifelse(income %in% c(1, 2), "low", ifelse(income == 3, "middle", ifelse(income %in% 
                c(4, 5), "high", NA))), levels = c("low", "middle", "high"))
            ed <- bi_education(get("education"))
            d$education <- ed$education
            d$university <- ed$university
            smoking <- bi_num(get("smoking"))
            smoking[!smoking %in% 0:2] <- NA
            d$smoking <- factor(smoking, levels = 0:2)
            alcohol <- bi_num(get("alcohol"))
            alcohol[!alcohol %in% 0:2] <- NA
            d$alcohol <- factor(alcohol, levels = 0:2)
            d$family_composition <- factor(ifelse(d$cnt_brother > 0 & d$cnt_sister == 0, "Brothers_only", ifelse(d$cnt_sister > 
                0 & d$cnt_brother == 0, "Sisters_only", ifelse(d$cnt_brother > 0 & d$cnt_sister > 0, "Mixed", NA))))
            d$sibling_sex <- factor(ifelse(d$S == 2, as.character(d$family_composition), NA), levels = c("Brothers_only", 
                "Sisters_only"))
            ctx$dat <- d
            bi_log(ctx, "COHORT_READY", rows = nrow(d))
            d
        }
        bi_outcomes <- function(ctx) {
            if (exists("ph", envir = ctx, inherits = FALSE)) return(ctx$ph)
            d <- bi_cohort(ctx)
            p <- as.data.frame(readRDS(ctx$cfg$inputs$phewas_rds))
            bi_assert("id" %in% names(p) && !anyNA(p$id) && !anyDuplicated(as.character(p$id)), "PHENOTYPE_ID_INVALID")
            codes <- unique(c(unlist(ctx$cfg$clinical$core), unlist(ctx$cfg$clinical$supplementary)))
            missing_codes <- setdiff(unlist(ctx$cfg$clinical$core), names(p))
            if (length(missing_codes)) bi_log(ctx, "OUTCOME_COLUMNS_UNAVAILABLE", phecodes = paste(missing_codes, 
                collapse = ";"), action = "Keep registered outcomes as NA; never code unknown as control")
            ix <- match(as.character(d$eid), as.character(p$id))
            out <- data.frame(eid = d$eid)
            for (k in codes) {
                y <- if (k %in% names(p)) p[[k]][ix] else rep(NA, nrow(d))
                bi_assert(all(is.na(y) | y %in% c(0, 1)), paste("INVALID_OUTCOME", k))
                out[[k]] <- as.numeric(y)
            }
            ctx$ph <- out
            out
        }
        bi_data <- function(ctx, code) {
            d <- bi_cohort(ctx)
            p <- bi_outcomes(ctx)
            bi_assert(code %in% names(p), "OUTCOME_NOT_REGISTERED")
            bi_assert(identical(as.character(d$eid), as.character(p$eid)), "OUTCOME_ALIGNMENT_FAILURE")
            d$y <- p[[code]]
            d
        }
    }), envir = context)
}

