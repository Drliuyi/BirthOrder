# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

c39_num_assert <- function(ok, code) if (!isTRUE(ok)) stop(code, call. = FALSE)

c39_num_spec <- function() list(version = "CORE39_UNIFORM_PIECEWISE_EXPONENTIAL_DOUBLE_ONLY_V1", implementation_version = "EXACT_PURE_CAUSE_SEMIGROUP_ONE_EXP_V1", 
    pure_cause_compression = "consecutive exact-zero opposite cause only; mixed-cause barriers; overflow fallback uses all original jumps", 
    method = "UNIFORM_PIECEWISE_CONSTANT_COMPETING_HAZARD_MATRIX_EXPONENTIAL", integration_method = "UNIFORM_PIECEWISE_CONSTANT_COMPETING_HAZARD_MATRIX_EXPONENTIAL", 
    previous_approved_method = "UNIFORM_PIECEWISE_CONSTANT_COMPETING_HAZARD_MATRIX_EXPONENTIAL", comparison_base = "exponential_b100_v2/CORE39_UNIFORM_PIECEWISE_EXPONENTIAL_V1", 
    computation_precision = "DOUBLE_ONLY_BINARY64", working_mantissa_bits = 53L, precision_change = TRUE, formula_change = FALSE, 
    method_change = FALSE, method_change_relative_to = "approved exponential_b100_v2, NOT legacy discrete AJ", 
    equivalent_to_old_AJ = FALSE, epsilon = 0, all_native_and_optimized_paths_same_kernel = TRUE, nonfinite_policy = "FAIL_COMPLETE_INPUT_BEFORE_INTEGRATION", 
    negative_hazard_policy = "FAIL", clipping = FALSE, renormalization = FALSE, survival_tail_cutoff = FALSE, stable_arithmetic = paste("binary64 expm1/exp/log; compensated CIF sums;", 
        "representability-checked hazard scaling; overflow log-sum-exp;", "subnormal survival exponential half-step"), 
    mass_residual_tolerance = 1e-10, relation_residual_tolerance = 1e-10, formally_certified = FALSE, certificate = "NONE: libm exp/log are not directed interval arithmetic", 
    error_scope = "Observed floating-point mass and algebraic consistency, NOT a CIF error bound", scientific_scope = "Linear interpolation of cumulative hazards on each common event interval; distinct from discrete AJ and Cox exp-left integral")

c39_num_reject_unsafe_flags <- function(flags) {
    unsafe <- paste(c("-ffast-math", "-Ofast", "-funsafe-math-optimizations", "-fassociative-math", "-freciprocal-math", 
        "-ffinite-math-only(?:=[1-9][0-9]*)?", "-fno-signed-zeros", "-fno-rounding-math", "-fexcess-precision=fast", 
        "-ffp-contract=fast", "-mfpmath=387", "-fapprox-func", "-fno-honor-infinities", "-fno-honor-nans", "-menable-unsafe-fp-math", 
        "-fdenormal-fp-math=[^[:space:]]+", "(?:/|-)fp:fast", "-fp-model[=[:space:]]+fast[^[:space:]]*"), collapse = "|")
    c39_num_assert(!any(grepl(paste0("(^|[[:space:]\"'])\\s*(?:", unsafe, ")(?=$|[[:space:]\"'])"), flags, perl = TRUE)), 
        "UNSAFE_FLOATING_POINT_FLAGS")
    invisible(TRUE)
}

c39_num_initialize <- function(source_dir, build_dir) {
    source_dir <- normalizePath(source_dir, mustWork = TRUE)
    src <- file.path(source_dir, "piecewise_exponential.c")
    c39_num_assert(file.exists(src), "NUM_DOUBLE_ONLY_SOURCE_MISSING")
    c39_num_assert(!file.exists(build_dir), "NUM_BUILD_MUST_BE_NEW_ADDITIVE_DIRECTORY")
    c39_num_assert(requireNamespace("digest", quietly = TRUE), "DIGEST_REQUIRED_NO_INSTALL")
    rbin <- file.path(R.home("bin"), "R")
    config_names <- c("CC", "CFLAGS", "CPICFLAGS", "CPPFLAGS", "LDFLAGS", "SHLIB_CFLAGS")
    compiler_config <- setNames(lapply(config_names, function(key) {
        z <- system2(rbin, c("CMD", "config", key), stdout = TRUE, stderr = TRUE)
        c39_num_assert(is.null(attr(z, "status")) || attr(z, "status") == 0L, paste0("NUM_COMPILER_CONFIG_FAILED: ", 
            key, " ", paste(z, collapse = "\n")))
        z
    }), config_names)
    compiler_environment <- Sys.getenv(c("CC", "CFLAGS", "CPICFLAGS", "CPPFLAGS", "LDFLAGS", "SHLIB_CFLAGS", "PKG_CFLAGS", 
        "PKG_CPPFLAGS", "PKG_LIBS", "MAKEFLAGS", "R_MAKEVARS_USER", "R_MAKEVARS_SITE"), unset = "")
    c39_num_reject_unsafe_flags(c(unlist(compiler_config, use.names = FALSE), compiler_environment))
    c39_num_assert(dir.create(build_dir, recursive = TRUE), "NUM_BUILD_CREATE_FAILED")
    build_dir <- normalizePath(build_dir, mustWork = TRUE)
    dest <- file.path(build_dir, basename(src))
    c39_num_assert(file.copy(src, dest, overwrite = FALSE), "NUM_SOURCE_COPY_FAILED")
    sha <- digest::digest(file = src, algo = "sha256")
    c39_num_assert(identical(sha, digest::digest(file = dest, algo = "sha256")), "NUM_SOURCE_COPY_HASH")
    prev <- setwd(build_dir)
    on.exit(setwd(prev), add = TRUE)
    binary <- file.path(build_dir, paste0("c39_exponential", .Platform$dynlib.ext))
    log <- system2(rbin, c("CMD", "SHLIB", "-o", basename(binary), basename(src)), stdout = TRUE, stderr = TRUE)
    writeLines(log, file.path(build_dir, "COMPILER_OUTPUT.log"))
    c39_num_assert(is.null(attr(log, "status")) || attr(log, "status") == 0L, paste(log, collapse = "\n"))
    c39_num_reject_unsafe_flags(log)
    dll <- dyn.load(binary, local = TRUE)
    bits <- .Call(getNativeSymbolInfo("c39_exp_environment", PACKAGE = dll))
    c39_num_assert(identical(bits, 53L), "NUM_DOUBLE_ONLY_MANTISSA_NOT_53_BITS")
    manifest <- list(source_sha256 = sha, binary_sha256 = digest::digest(file = binary, algo = "sha256"), R = R.version.string, 
        platform = R.version$platform, compiler_output = log, compiler_config = compiler_config, compiler_environment = compiler_environment, 
        spec = c39_num_spec(), floating_point_checks = list(computation_precision = "DOUBLE_ONLY_BINARY64", binary64 = TRUE, 
            rounding = "FE_TONEAREST", gradual_underflow = TRUE, double_mantissa_bits = bits, intermediate_evaluation = "FLT_EVAL_METHOD=0", 
            unsafe_flags_rejected = TRUE, compile_time_fast_math_guard = TRUE, checked_at_every_C_entry = TRUE), 
        formally_certified = FALSE)
    saveRDS(manifest, file.path(build_dir, "BUILD_MANIFEST.rds"))
    list(symbol = getNativeSymbolInfo("c39_piecewise_exponential", PACKAGE = dll), manifest = manifest, binary = binary)
}

c39_num_run <- function(hazard, relative, backend, epsilon = 0) {
    c39_num_assert(identical(as.double(epsilon), 0), "EXP_NO_TAIL_EPSILON_ALLOWED")
    storage.mode(hazard) <- storage.mode(relative) <- "double"
    x <- .Call(backend$symbol, hazard, relative)
    p <- x[, 1:3, drop = FALSE]
    colnames(p) <- c("disease", "death", "survival")
    c39_num_assert(all(is.finite(p)) && all(p >= 0) && all(p <= 1), "EXP_INVALID_PROBABILITY_NO_CLIP")
    c39_num_assert(all(x[, 4] <= c39_num_spec()$mass_residual_tolerance) && all(x[, 5] <= c39_num_spec()$relation_residual_tolerance), 
        "EXP_MASS_OR_RELATION_RESIDUAL_EXCEEDED")
    list(probability = p, mass_residual = x[, 4], relation_residual = x[, 5], overflow_scaled_steps = x[, 6], underflow_survival = x[, 
        7] > 0, formally_certified = FALSE)
}

c39_num_quadrature <- function(hazard, relative) {
    c39_num_assert(is.matrix(hazard) && ncol(hazard) == 2L && is.matrix(relative) && ncol(relative) == 2L && nrow(relative) >= 
        1L && all(is.finite(hazard)) && all(hazard >= 0) && all(is.finite(relative)) && all(relative > 0), "EXP_REFERENCE_INPUT_INVALID")
    t(vapply(seq_len(nrow(relative)), function(i) {
        h <- sweep(hazard, 2, relative[i, ], `*`)
        q <- rowSums(h)
        c39_num_assert(all(is.finite(q)), "EXP_QUADRATURE_FINITE_PRODUCTS_REQUIRED")
        if (!length(q)) 
            return(c(0, 0, 1))
        before <- exp(-c(0, head(cumsum(q), -1L)))
        factors <- vapply(q, function(a) if (a == 0) 
            1
        else if (a <= 1) 
            stats::integrate(function(t) exp(-a * t), 0, 1, rel.tol = 1e-11, abs.tol = 1e-13)$value
        else stats::integrate(function(u) {
            z <- numeric(length(u))
            k <- u < 1
            z[k] <- exp(-u[k]/(1 - u[k]) - 2 * log1p(-u[k]))
            z
        }, 0, 1/(1 + 1/a), rel.tol = 1e-11, abs.tol = 1e-13)$value/a, numeric(1))
        c(colSums(h * (before * factors)), exp(-sum(q)))
    }, numeric(3)))
}

c39_num_install <- function(engine, backend) {
    state <- new.env(parent = emptyenv())
    state$last <- NULL
    state$calls <- 0L
    run <- function(h, r) {
        z <- c39_num_run(h, r, backend)
        state$last <- z
        state$calls <- state$calls + 1L
        z$probability
    }
    engine$ff_aj_backend <- run
    engine$ff_y_aj <- function(hazard, rr, period) {
        c39_num_assert(all(period == 1L), "EXP_ONLY_FROZEN_SINGLE_PERIOD")
        run(hazard, rr[, 1:2, drop = FALSE])
    }
    engine$ff_aj_status <- function(require = FALSE) c(list(ready = TRUE), c39_num_spec(), backend$manifest)
    engine$c39_num_last <- function() state$last
    engine$c39_num_calls <- function() state$calls
    invisible(engine)
}

c39_num_empty_audit <- function() list(backend_calls = 0, row_paths = 0, max_mass_residual = 0, max_relation_residual = 0, 
    min_probability = 1, max_probability = 0, overflow_scaled_steps = 0, underflow_survival_paths = 0, negative_or_nonfinite_inputs = 0, 
    clipping_events = 0)

c39_num_add_audit <- function(a, b) {
    for (n in c("backend_calls", "row_paths", "overflow_scaled_steps", "underflow_survival_paths", "negative_or_nonfinite_inputs", 
        "clipping_events")) a[[n]] <- a[[n]] + b[[n]]
    for (n in c("max_mass_residual", "max_relation_residual", "max_probability")) a[[n]] <- max(a[[n]], b[[n]])
    a$min_probability <- min(a$min_probability, b$min_probability)
    a
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        c39_num_assert <- function(ok, code) if (!isTRUE(ok)) stop(code, call. = FALSE)
        c39_num_spec <- function() list(version = "CORE39_UNIFORM_PIECEWISE_EXPONENTIAL_DOUBLE_ONLY_V1", implementation_version = "EXACT_PURE_CAUSE_SEMIGROUP_ONE_EXP_V1", 
            pure_cause_compression = "consecutive exact-zero opposite cause only; mixed-cause barriers; overflow fallback uses all original jumps", 
            method = "UNIFORM_PIECEWISE_CONSTANT_COMPETING_HAZARD_MATRIX_EXPONENTIAL", integration_method = "UNIFORM_PIECEWISE_CONSTANT_COMPETING_HAZARD_MATRIX_EXPONENTIAL", 
            previous_approved_method = "UNIFORM_PIECEWISE_CONSTANT_COMPETING_HAZARD_MATRIX_EXPONENTIAL", comparison_base = "exponential_b100_v2/CORE39_UNIFORM_PIECEWISE_EXPONENTIAL_V1", 
            computation_precision = "DOUBLE_ONLY_BINARY64", working_mantissa_bits = 53L, precision_change = TRUE, 
            formula_change = FALSE, method_change = FALSE, method_change_relative_to = "approved exponential_b100_v2, NOT legacy discrete AJ", 
            equivalent_to_old_AJ = FALSE, epsilon = 0, all_native_and_optimized_paths_same_kernel = TRUE, nonfinite_policy = "FAIL_COMPLETE_INPUT_BEFORE_INTEGRATION", 
            negative_hazard_policy = "FAIL", clipping = FALSE, renormalization = FALSE, survival_tail_cutoff = FALSE, 
            stable_arithmetic = paste("binary64 expm1/exp/log; compensated CIF sums;", "representability-checked hazard scaling; overflow log-sum-exp;", 
                "subnormal survival exponential half-step"), mass_residual_tolerance = 1e-10, relation_residual_tolerance = 1e-10, 
            formally_certified = FALSE, certificate = "NONE: libm exp/log are not directed interval arithmetic", 
            error_scope = "Observed floating-point mass and algebraic consistency, NOT a CIF error bound", scientific_scope = "Linear interpolation of cumulative hazards on each common event interval; distinct from discrete AJ and Cox exp-left integral")
        c39_num_reject_unsafe_flags <- function(flags) {
            unsafe <- paste(c("-ffast-math", "-Ofast", "-funsafe-math-optimizations", "-fassociative-math", "-freciprocal-math", 
                "-ffinite-math-only(?:=[1-9][0-9]*)?", "-fno-signed-zeros", "-fno-rounding-math", "-fexcess-precision=fast", 
                "-ffp-contract=fast", "-mfpmath=387", "-fapprox-func", "-fno-honor-infinities", "-fno-honor-nans", 
                "-menable-unsafe-fp-math", "-fdenormal-fp-math=[^[:space:]]+", "(?:/|-)fp:fast", "-fp-model[=[:space:]]+fast[^[:space:]]*"), 
                collapse = "|")
            c39_num_assert(!any(grepl(paste0("(^|[[:space:]\"'])\\s*(?:", unsafe, ")(?=$|[[:space:]\"'])"), flags, 
                perl = TRUE)), "UNSAFE_FLOATING_POINT_FLAGS")
            invisible(TRUE)
        }
        c39_num_initialize <- function(source_dir, build_dir) {
            source_dir <- normalizePath(source_dir, mustWork = TRUE)
            src <- file.path(source_dir, "piecewise_exponential.c")
            c39_num_assert(file.exists(src), "NUM_DOUBLE_ONLY_SOURCE_MISSING")
            c39_num_assert(!file.exists(build_dir), "NUM_BUILD_MUST_BE_NEW_ADDITIVE_DIRECTORY")
            c39_num_assert(requireNamespace("digest", quietly = TRUE), "DIGEST_REQUIRED_NO_INSTALL")
            rbin <- file.path(R.home("bin"), "R")
            config_names <- c("CC", "CFLAGS", "CPICFLAGS", "CPPFLAGS", "LDFLAGS", "SHLIB_CFLAGS")
            compiler_config <- setNames(lapply(config_names, function(key) {
                z <- system2(rbin, c("CMD", "config", key), stdout = TRUE, stderr = TRUE)
                c39_num_assert(is.null(attr(z, "status")) || attr(z, "status") == 0L, paste0("NUM_COMPILER_CONFIG_FAILED: ", 
                  key, " ", paste(z, collapse = "\n")))
                z
            }), config_names)
            compiler_environment <- Sys.getenv(c("CC", "CFLAGS", "CPICFLAGS", "CPPFLAGS", "LDFLAGS", "SHLIB_CFLAGS", 
                "PKG_CFLAGS", "PKG_CPPFLAGS", "PKG_LIBS", "MAKEFLAGS", "R_MAKEVARS_USER", "R_MAKEVARS_SITE"), unset = "")
            c39_num_reject_unsafe_flags(c(unlist(compiler_config, use.names = FALSE), compiler_environment))
            c39_num_assert(dir.create(build_dir, recursive = TRUE), "NUM_BUILD_CREATE_FAILED")
            build_dir <- normalizePath(build_dir, mustWork = TRUE)
            dest <- file.path(build_dir, basename(src))
            c39_num_assert(file.copy(src, dest, overwrite = FALSE), "NUM_SOURCE_COPY_FAILED")
            sha <- digest::digest(file = src, algo = "sha256")
            c39_num_assert(identical(sha, digest::digest(file = dest, algo = "sha256")), "NUM_SOURCE_COPY_HASH")
            prev <- setwd(build_dir)
            on.exit(setwd(prev), add = TRUE)
            binary <- file.path(build_dir, paste0("c39_exponential", .Platform$dynlib.ext))
            log <- system2(rbin, c("CMD", "SHLIB", "-o", basename(binary), basename(src)), stdout = TRUE, stderr = TRUE)
            writeLines(log, file.path(build_dir, "COMPILER_OUTPUT.log"))
            c39_num_assert(is.null(attr(log, "status")) || attr(log, "status") == 0L, paste(log, collapse = "\n"))
            c39_num_reject_unsafe_flags(log)
            dll <- dyn.load(binary, local = TRUE)
            bits <- .Call(getNativeSymbolInfo("c39_exp_environment", PACKAGE = dll))
            c39_num_assert(identical(bits, 53L), "NUM_DOUBLE_ONLY_MANTISSA_NOT_53_BITS")
            manifest <- list(source_sha256 = sha, binary_sha256 = digest::digest(file = binary, algo = "sha256"), 
                R = R.version.string, platform = R.version$platform, compiler_output = log, compiler_config = compiler_config, 
                compiler_environment = compiler_environment, spec = c39_num_spec(), floating_point_checks = list(computation_precision = "DOUBLE_ONLY_BINARY64", 
                  binary64 = TRUE, rounding = "FE_TONEAREST", gradual_underflow = TRUE, double_mantissa_bits = bits, 
                  intermediate_evaluation = "FLT_EVAL_METHOD=0", unsafe_flags_rejected = TRUE, compile_time_fast_math_guard = TRUE, 
                  checked_at_every_C_entry = TRUE), formally_certified = FALSE)
            saveRDS(manifest, file.path(build_dir, "BUILD_MANIFEST.rds"))
            list(symbol = getNativeSymbolInfo("c39_piecewise_exponential", PACKAGE = dll), manifest = manifest, 
                binary = binary)
        }
        c39_num_run <- function(hazard, relative, backend, epsilon = 0) {
            c39_num_assert(identical(as.double(epsilon), 0), "EXP_NO_TAIL_EPSILON_ALLOWED")
            storage.mode(hazard) <- storage.mode(relative) <- "double"
            x <- .Call(backend$symbol, hazard, relative)
            p <- x[, 1:3, drop = FALSE]
            colnames(p) <- c("disease", "death", "survival")
            c39_num_assert(all(is.finite(p)) && all(p >= 0) && all(p <= 1), "EXP_INVALID_PROBABILITY_NO_CLIP")
            c39_num_assert(all(x[, 4] <= c39_num_spec()$mass_residual_tolerance) && all(x[, 5] <= c39_num_spec()$relation_residual_tolerance), 
                "EXP_MASS_OR_RELATION_RESIDUAL_EXCEEDED")
            list(probability = p, mass_residual = x[, 4], relation_residual = x[, 5], overflow_scaled_steps = x[, 
                6], underflow_survival = x[, 7] > 0, formally_certified = FALSE)
        }
        c39_num_quadrature <- function(hazard, relative) {
            c39_num_assert(is.matrix(hazard) && ncol(hazard) == 2L && is.matrix(relative) && ncol(relative) == 
                2L && nrow(relative) >= 1L && all(is.finite(hazard)) && all(hazard >= 0) && all(is.finite(relative)) && 
                all(relative > 0), "EXP_REFERENCE_INPUT_INVALID")
            t(vapply(seq_len(nrow(relative)), function(i) {
                h <- sweep(hazard, 2, relative[i, ], `*`)
                q <- rowSums(h)
                c39_num_assert(all(is.finite(q)), "EXP_QUADRATURE_FINITE_PRODUCTS_REQUIRED")
                if (!length(q)) return(c(0, 0, 1))
                before <- exp(-c(0, head(cumsum(q), -1L)))
                factors <- vapply(q, function(a) if (a == 0) 1 else if (a <= 1) stats::integrate(function(t) exp(-a * 
                  t), 0, 1, rel.tol = 1e-11, abs.tol = 1e-13)$value else stats::integrate(function(u) {
                  z <- numeric(length(u))
                  k <- u < 1
                  z[k] <- exp(-u[k]/(1 - u[k]) - 2 * log1p(-u[k]))
                  z
                }, 0, 1/(1 + 1/a), rel.tol = 1e-11, abs.tol = 1e-13)$value/a, numeric(1))
                c(colSums(h * (before * factors)), exp(-sum(q)))
            }, numeric(3)))
        }
        c39_num_install <- function(engine, backend) {
            state <- new.env(parent = emptyenv())
            state$last <- NULL
            state$calls <- 0L
            run <- function(h, r) {
                z <- c39_num_run(h, r, backend)
                state$last <- z
                state$calls <- state$calls + 1L
                z$probability
            }
            engine$ff_aj_backend <- run
            engine$ff_y_aj <- function(hazard, rr, period) {
                c39_num_assert(all(period == 1L), "EXP_ONLY_FROZEN_SINGLE_PERIOD")
                run(hazard, rr[, 1:2, drop = FALSE])
            }
            engine$ff_aj_status <- function(require = FALSE) c(list(ready = TRUE), c39_num_spec(), backend$manifest)
            engine$c39_num_last <- function() state$last
            engine$c39_num_calls <- function() state$calls
            invisible(engine)
        }
        c39_num_empty_audit <- function() list(backend_calls = 0, row_paths = 0, max_mass_residual = 0, max_relation_residual = 0, 
            min_probability = 1, max_probability = 0, overflow_scaled_steps = 0, underflow_survival_paths = 0, 
            negative_or_nonfinite_inputs = 0, clipping_events = 0)
        c39_num_add_audit <- function(a, b) {
            for (n in c("backend_calls", "row_paths", "overflow_scaled_steps", "underflow_survival_paths", "negative_or_nonfinite_inputs", 
                "clipping_events")) a[[n]] <- a[[n]] + b[[n]]
            for (n in c("max_mass_residual", "max_relation_residual", "max_probability")) a[[n]] <- max(a[[n]], 
                b[[n]])
            a$min_probability <- min(a$min_probability, b$min_probability)
            a
        }
    }), envir = context)
}

