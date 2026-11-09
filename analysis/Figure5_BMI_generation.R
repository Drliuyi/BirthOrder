# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

c39_loga_assert <- function(ok, message) {
    if (!isTRUE(ok)) 
        stop(message, call. = FALSE)
}

c39_loga_ids <- function() c("bmi")

c39_loga_revision <- function() "CORE39_BMI_LOG_FAT_LINEAR_V1"

c39_loga_hash <- function(x) {
    c39_loga_assert(requireNamespace("digest", quietly = TRUE), "LOG_A_DIGEST_REQUIRED")
    digest::digest(x, algo = "sha256", serialize = TRUE, serializeVersion = 2L)
}

c39_loga_rng <- function() {
    present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    list(kind = RNGkind(), present = present, seed = if (present) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL)
}

c39_loga_assert_rng <- function(before) {
    c39_loga_assert(identical(before, c39_loga_rng()), "LOG_A_RNG_STATE_CHANGED")
    invisible(TRUE)
}

c39_loga_linear_fit <- function(X, y, w, id, rule) {
    if (!id %in% c39_loga_ids()) 
        return(.c39_loga_original_linear_fit(X, y, w, id, rule))
    c39_loga_assert(is.matrix(X) && is.numeric(X) && all(is.finite(X)) && is.numeric(y) && length(y) == nrow(X) && 
        length(y) > 0L && all(is.finite(y) & y > 0), paste("LOG_A_OBSERVATIONS_MUST_BE_FINITE_STRICTLY_POSITIVE", 
        id))
    c39_loga_assert(is.numeric(w) && length(w) == length(y) && all(is.finite(w) & w > 0) && is.finite(sum(w)), 
        "LOG_A_INVALID_EXACT_FIT_WEIGHTS")
    response <- log(y)
    c39_loga_assert(all(is.finite(response)), paste("LOG_A_NONFINITE_LOG_RESPONSE", id))
    f <- .c39_loga_original_linear_fit(X, response, w, id, rule)
    if (identical(f$kind, "continuous")) {
        c39_loga_assert(all(is.finite(f$coefficients)) && all(is.finite(f$residuals)), paste("LOG_A_NONFINITE_FIT", 
            id))
    }
    else c39_loga_assert(identical(f$kind, "constant") && is.finite(f$value), "LOG_A_UNEXPECTED_FIT_KIND")
    f$log_a_revision <- list(revision = c39_loga_revision(), indicator_id = id, response_transform = "log", inverse_transform = "exp", 
        predictor_scale = "ORIGINAL_UNITS_UNCHANGED", downstream_scale = "ORIGINAL_UNITS", empirical_residual = "ORIGINAL_WEIGHTED_EMPIRICAL_LAW_ON_LOG_RESPONSE", 
        observed_n = length(y), observed_positive_n = sum(y > 0), observed_min = min(y), observed_max = max(y), 
        offset = 0, clipping = FALSE, reseeding = FALSE)
    f
}

c39_loga_residual_draw <- function(f, X, u, lower, upper, lower_open, id, audit = NULL) {
    if (is.null(f$log_a_revision)) 
        return(.c39_loga_original_residual_draw(f, X, u, lower, upper, lower_open, id, audit))
    tag <- f$log_a_revision
    c39_loga_assert(id %in% c39_loga_ids() && identical(tag$indicator_id, id) && identical(tag$revision, c39_loga_revision()) && 
        identical(tag$response_transform, "log") && identical(tag$inverse_transform, "exp"), "LOG_A_INVALID_MODEL_MARKER")
    n <- nrow(X)
    lower <- rep_len(lower, n)
    upper <- rep_len(upper, n)
    lower_open <- rep_len(lower_open, n)
    c39_loga_assert(all(is.finite(X)) && length(u) == n && all(is.finite(u) & u >= 0 & u < 1) && !anyNA(lower) && 
        !anyNA(upper) && !anyNA(lower_open) && all(upper > 0 & upper >= lower), paste("LOG_A_INVALID_DRAW_INPUT", 
        id))
    log_lower <- rep(-Inf, n)
    positive_lower <- lower > 0
    log_lower[positive_lower] <- log(lower[positive_lower])
    value_log <- .c39_loga_original_residual_draw(f, X, u, log_lower, log(upper), lower_open & positive_lower, 
        id, audit)
    value <- exp(value_log)
    legal <- ifelse(lower_open, value > lower, value >= lower) & value <= upper
    c39_loga_assert(all(is.finite(value_log)) && all(is.finite(value) & value > 0 & legal), paste("LOG_A_EXP_DRAW_NOT_FINITE_STRICTLY_POSITIVE_NO_CLIPPING", 
        id))
    value
}

c39_loga_install <- function(engine) {
    before_rng <- c39_loga_rng()
    on.exit(c39_loga_assert_rng(before_rng), add = TRUE)
    c39_loga_assert(is.environment(engine) && !identical(engine, .GlobalEnv) && !identical(engine, baseenv()) && 
        !identical(engine, emptyenv()), "LOG_A_ISOLATED_ENGINE_REQUIRED")
    required <- c("c39_linear_fit", "c39_residual_draw", "c39_a_design", "c39_fit_a", "c39_draw_a")
    c39_loga_assert(all(vapply(required, function(n) exists(n, engine, inherits = FALSE) && is.function(get(n, 
        engine, inherits = FALSE)), logical(1))), "LOG_A_INSTALL_CORE39_MODELS_FIRST")
    c39_loga_assert(identical(environment(engine$c39_fit_a), engine) && identical(environment(engine$c39_draw_a), 
        engine), "LOG_A_ENGINE_FUNCTION_ENVIRONMENT_DRIFT")
    if (exists(".c39_loga_installed", engine, inherits = FALSE)) {
        c39_loga_assert(identical(engine$.c39_loga_installed, c39_loga_revision()) && identical(engine$c39_linear_fit, 
            engine$c39_loga_linear_fit) && identical(engine$c39_residual_draw, engine$c39_loga_residual_draw), 
            "LOG_A_INSTALLATION_DRIFT")
        return(invisible(engine))
    }
    c39_loga_assert(!exists(".c39_loga_original_linear_fit", engine, inherits = FALSE) && !exists(".c39_loga_original_residual_draw", 
        engine, inherits = FALSE), "LOG_A_RESERVED_ENGINE_NAMES_EXIST")
    engine$.c39_loga_original_linear_fit <- engine$c39_linear_fit
    engine$.c39_loga_original_residual_draw <- engine$c39_residual_draw
    source_environment <- environment(c39_loga_install)
    for (name in c("c39_loga_assert", "c39_loga_ids", "c39_loga_revision", "c39_loga_linear_fit", "c39_loga_residual_draw")) {
        f <- get(name, source_environment, inherits = FALSE)
        environment(f) <- engine
        assign(name, f, envir = engine)
    }
    engine$c39_linear_fit <- engine$c39_loga_linear_fit
    engine$c39_residual_draw <- engine$c39_loga_residual_draw
    engine$.c39_loga_installed <- c39_loga_revision()
    invisible(engine)
}

c39_mixeda_install <- function(engine) c39_loga_install(engine)

c39_mixeda_restore_bundle <- function(original, previous, weights, engine) {
    before <- c39_loga_rng()
    on.exit(c39_loga_assert_rng(before), add = TRUE)
    assert <- c39_loga_assert
    assert(is.list(original$a$models) && length(original$a$models) == 39L && identical(names(original$a$models), 
        names(previous$a$models)), "MIXED_A_SCHEMA_MISMATCH")
    for (n in setdiff(names(original), "a")) assert(identical(original[[n]], previous[[n]]), paste("MIXED_A_NON_A_CACHE_DIFFERS", 
        n))
    for (n in setdiff(names(original$a), "models")) assert(identical(original$a[[n]], previous$a[[n]]), paste("MIXED_A_METADATA_DIFFERS", 
        n))
    for (id in setdiff(names(original$a$models), c("bmi", "fat_mass"))) assert(identical(original$a$models[[id]], 
        previous$a$models[[id]]), paste("MIXED_A_OTHER_MODEL_DIFFERS", id))
    assert(identical(as.numeric(weights), as.numeric(original$weights)), "MIXED_A_STANDARDIZATION_WEIGHTS_DIFFER")
    for (ep in names(original$fits)) assert(identical(as.numeric(weights), as.numeric(original$fits[[ep]]$binding$weight)), 
        "MIXED_A_FIT_WEIGHTS_DIFFER")
    old_bmi <- original$a$models$bmi
    old_fat <- original$a$models$fat_mass
    new_bmi <- previous$a$models$bmi
    assert(is.null(old_bmi$log_a_revision) && is.null(old_fat$log_a_revision) && identical(old_fat$kind, "continuous") && 
        identical(old_fat$rule, "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL"), "MIXED_A_ORIGINAL_FAT_MASS_LINEAR_CONDITIONAL_MODEL_REQUIRED")
    assert(identical(new_bmi$log_a_revision$revision, "CORE39_LOG_A_BMI_FAT_MASS_V1") && identical(new_bmi$log_a_revision$response_transform, 
        "log") && identical(new_bmi$log_a_revision$inverse_transform, "exp"), "MIXED_A_REVIEWED_LOG_BMI_CACHE_REQUIRED")
    meta <- setdiff(names(old_bmi), c("kind", "value", "support", "coefficients", "residuals", "cumulative", "cumulative_squared_mass", 
        "rule"))
    assert(identical(old_bmi[meta], new_bmi[meta]), "MIXED_A_BMI_DESIGN_CHANGED")
    new_bmi$log_a_revision$revision <- c39_loga_revision()
    amended <- original
    amended$a$models$bmi <- new_bmi
    protected <- setdiff(names(original$a$models), "bmi")
    assert(identical(amended$a$models[protected], original$a$models[protected]), "MIXED_A_ORIGINAL38_NOT_PRESERVED")
    assert(identical(amended$fits, original$fits) && identical(amended$pop$data, original$pop$data) && identical(amended$target, 
        original$target) && identical(amended$weights, original$weights), "MIXED_A_TARGET_Y_OR_WEIGHTS_CHANGED")
    positives <- lapply(c("bmi", "fat_mass"), function(id) {
        x <- original$pop$data[[original$a$members[[id]]$column]]
        assert(all(is.finite(x) & x > 0), paste("MIXED_A_NONPOSITIVE_OBSERVATION", id))
        list(indicator_id = id, N = length(x), positive_n = sum(x > 0), min = min(x), max = max(x))
    })
    list(bundle = amended, summary = list(revision = c39_loga_revision(), A_refits = 0L, Y_refits = 0L, reused_log_response = "bmi", 
        restored_original_linear_response = "fat_mass", original_fat_mass_sha256 = c39_loga_hash(old_fat), restored_fat_mass_sha256 = c39_loga_hash(amended$a$models$fat_mass), 
        fat_mass_strict_identical = TRUE, other37_A_strict_identical = TRUE, all_Y_and_non_A_fields_strict_identical = TRUE, 
        weights_strict_identical = TRUE, observed_positive = positives, random_state_unchanged = TRUE, original_bmi_sha256 = c39_loga_hash(old_bmi), 
        reused_bmi_fit_sha256 = c39_loga_hash(new_bmi), no_offset_no_clipping_no_reseeding = TRUE))
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        c39_loga_assert <- function(ok, message) {
            if (!isTRUE(ok)) stop(message, call. = FALSE)
        }
        c39_loga_ids <- function() c("bmi")
        c39_loga_revision <- function() "CORE39_BMI_LOG_FAT_LINEAR_V1"
        c39_loga_hash <- function(x) {
            c39_loga_assert(requireNamespace("digest", quietly = TRUE), "LOG_A_DIGEST_REQUIRED")
            digest::digest(x, algo = "sha256", serialize = TRUE, serializeVersion = 2L)
        }
        c39_loga_rng <- function() {
            present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
            list(kind = RNGkind(), present = present, seed = if (present) get(".Random.seed", envir = .GlobalEnv, 
                inherits = FALSE) else NULL)
        }
        c39_loga_assert_rng <- function(before) {
            c39_loga_assert(identical(before, c39_loga_rng()), "LOG_A_RNG_STATE_CHANGED")
            invisible(TRUE)
        }
        c39_loga_linear_fit <- function(X, y, w, id, rule) {
            if (!id %in% c39_loga_ids()) return(.c39_loga_original_linear_fit(X, y, w, id, rule))
            c39_loga_assert(is.matrix(X) && is.numeric(X) && all(is.finite(X)) && is.numeric(y) && length(y) == 
                nrow(X) && length(y) > 0L && all(is.finite(y) & y > 0), paste("LOG_A_OBSERVATIONS_MUST_BE_FINITE_STRICTLY_POSITIVE", 
                id))
            c39_loga_assert(is.numeric(w) && length(w) == length(y) && all(is.finite(w) & w > 0) && is.finite(sum(w)), 
                "LOG_A_INVALID_EXACT_FIT_WEIGHTS")
            response <- log(y)
            c39_loga_assert(all(is.finite(response)), paste("LOG_A_NONFINITE_LOG_RESPONSE", id))
            f <- .c39_loga_original_linear_fit(X, response, w, id, rule)
            if (identical(f$kind, "continuous")) {
                c39_loga_assert(all(is.finite(f$coefficients)) && all(is.finite(f$residuals)), paste("LOG_A_NONFINITE_FIT", 
                  id))
            } else c39_loga_assert(identical(f$kind, "constant") && is.finite(f$value), "LOG_A_UNEXPECTED_FIT_KIND")
            f$log_a_revision <- list(revision = c39_loga_revision(), indicator_id = id, response_transform = "log", 
                inverse_transform = "exp", predictor_scale = "ORIGINAL_UNITS_UNCHANGED", downstream_scale = "ORIGINAL_UNITS", 
                empirical_residual = "ORIGINAL_WEIGHTED_EMPIRICAL_LAW_ON_LOG_RESPONSE", observed_n = length(y), 
                observed_positive_n = sum(y > 0), observed_min = min(y), observed_max = max(y), offset = 0, clipping = FALSE, 
                reseeding = FALSE)
            f
        }
        c39_loga_residual_draw <- function(f, X, u, lower, upper, lower_open, id, audit = NULL) {
            if (is.null(f$log_a_revision)) return(.c39_loga_original_residual_draw(f, X, u, lower, upper, lower_open, 
                id, audit))
            tag <- f$log_a_revision
            c39_loga_assert(id %in% c39_loga_ids() && identical(tag$indicator_id, id) && identical(tag$revision, 
                c39_loga_revision()) && identical(tag$response_transform, "log") && identical(tag$inverse_transform, 
                "exp"), "LOG_A_INVALID_MODEL_MARKER")
            n <- nrow(X)
            lower <- rep_len(lower, n)
            upper <- rep_len(upper, n)
            lower_open <- rep_len(lower_open, n)
            c39_loga_assert(all(is.finite(X)) && length(u) == n && all(is.finite(u) & u >= 0 & u < 1) && !anyNA(lower) && 
                !anyNA(upper) && !anyNA(lower_open) && all(upper > 0 & upper >= lower), paste("LOG_A_INVALID_DRAW_INPUT", 
                id))
            log_lower <- rep(-Inf, n)
            positive_lower <- lower > 0
            log_lower[positive_lower] <- log(lower[positive_lower])
            value_log <- .c39_loga_original_residual_draw(f, X, u, log_lower, log(upper), lower_open & positive_lower, 
                id, audit)
            value <- exp(value_log)
            legal <- ifelse(lower_open, value > lower, value >= lower) & value <= upper
            c39_loga_assert(all(is.finite(value_log)) && all(is.finite(value) & value > 0 & legal), paste("LOG_A_EXP_DRAW_NOT_FINITE_STRICTLY_POSITIVE_NO_CLIPPING", 
                id))
            value
        }
        c39_loga_install <- function(engine) {
            before_rng <- c39_loga_rng()
            on.exit(c39_loga_assert_rng(before_rng), add = TRUE)
            c39_loga_assert(is.environment(engine) && !identical(engine, .GlobalEnv) && !identical(engine, baseenv()) && 
                !identical(engine, emptyenv()), "LOG_A_ISOLATED_ENGINE_REQUIRED")
            required <- c("c39_linear_fit", "c39_residual_draw", "c39_a_design", "c39_fit_a", "c39_draw_a")
            c39_loga_assert(all(vapply(required, function(n) exists(n, engine, inherits = FALSE) && is.function(get(n, 
                engine, inherits = FALSE)), logical(1))), "LOG_A_INSTALL_CORE39_MODELS_FIRST")
            c39_loga_assert(identical(environment(engine$c39_fit_a), engine) && identical(environment(engine$c39_draw_a), 
                engine), "LOG_A_ENGINE_FUNCTION_ENVIRONMENT_DRIFT")
            if (exists(".c39_loga_installed", engine, inherits = FALSE)) {
                c39_loga_assert(identical(engine$.c39_loga_installed, c39_loga_revision()) && identical(engine$c39_linear_fit, 
                  engine$c39_loga_linear_fit) && identical(engine$c39_residual_draw, engine$c39_loga_residual_draw), 
                  "LOG_A_INSTALLATION_DRIFT")
                return(invisible(engine))
            }
            c39_loga_assert(!exists(".c39_loga_original_linear_fit", engine, inherits = FALSE) && !exists(".c39_loga_original_residual_draw", 
                engine, inherits = FALSE), "LOG_A_RESERVED_ENGINE_NAMES_EXIST")
            engine$.c39_loga_original_linear_fit <- engine$c39_linear_fit
            engine$.c39_loga_original_residual_draw <- engine$c39_residual_draw
            source_environment <- environment(c39_loga_install)
            for (name in c("c39_loga_assert", "c39_loga_ids", "c39_loga_revision", "c39_loga_linear_fit", "c39_loga_residual_draw")) {
                f <- get(name, source_environment, inherits = FALSE)
                environment(f) <- engine
                assign(name, f, envir = engine)
            }
            engine$c39_linear_fit <- engine$c39_loga_linear_fit
            engine$c39_residual_draw <- engine$c39_loga_residual_draw
            engine$.c39_loga_installed <- c39_loga_revision()
            invisible(engine)
        }
        c39_mixeda_install <- function(engine) c39_loga_install(engine)
        c39_mixeda_restore_bundle <- function(original, previous, weights, engine) {
            before <- c39_loga_rng()
            on.exit(c39_loga_assert_rng(before), add = TRUE)
            assert <- c39_loga_assert
            assert(is.list(original$a$models) && length(original$a$models) == 39L && identical(names(original$a$models), 
                names(previous$a$models)), "MIXED_A_SCHEMA_MISMATCH")
            for (n in setdiff(names(original), "a")) assert(identical(original[[n]], previous[[n]]), paste("MIXED_A_NON_A_CACHE_DIFFERS", 
                n))
            for (n in setdiff(names(original$a), "models")) assert(identical(original$a[[n]], previous$a[[n]]), 
                paste("MIXED_A_METADATA_DIFFERS", n))
            for (id in setdiff(names(original$a$models), c("bmi", "fat_mass"))) assert(identical(original$a$models[[id]], 
                previous$a$models[[id]]), paste("MIXED_A_OTHER_MODEL_DIFFERS", id))
            assert(identical(as.numeric(weights), as.numeric(original$weights)), "MIXED_A_STANDARDIZATION_WEIGHTS_DIFFER")
            for (ep in names(original$fits)) assert(identical(as.numeric(weights), as.numeric(original$fits[[ep]]$binding$weight)), 
                "MIXED_A_FIT_WEIGHTS_DIFFER")
            old_bmi <- original$a$models$bmi
            old_fat <- original$a$models$fat_mass
            new_bmi <- previous$a$models$bmi
            assert(is.null(old_bmi$log_a_revision) && is.null(old_fat$log_a_revision) && identical(old_fat$kind, 
                "continuous") && identical(old_fat$rule, "EMPIRICAL_RESIDUAL_PHYSICAL_SUPPORT_CONDITIONAL"), "MIXED_A_ORIGINAL_FAT_MASS_LINEAR_CONDITIONAL_MODEL_REQUIRED")
            assert(identical(new_bmi$log_a_revision$revision, "CORE39_LOG_A_BMI_FAT_MASS_V1") && identical(new_bmi$log_a_revision$response_transform, 
                "log") && identical(new_bmi$log_a_revision$inverse_transform, "exp"), "MIXED_A_REVIEWED_LOG_BMI_CACHE_REQUIRED")
            meta <- setdiff(names(old_bmi), c("kind", "value", "support", "coefficients", "residuals", "cumulative", 
                "cumulative_squared_mass", "rule"))
            assert(identical(old_bmi[meta], new_bmi[meta]), "MIXED_A_BMI_DESIGN_CHANGED")
            new_bmi$log_a_revision$revision <- c39_loga_revision()
            amended <- original
            amended$a$models$bmi <- new_bmi
            protected <- setdiff(names(original$a$models), "bmi")
            assert(identical(amended$a$models[protected], original$a$models[protected]), "MIXED_A_ORIGINAL38_NOT_PRESERVED")
            assert(identical(amended$fits, original$fits) && identical(amended$pop$data, original$pop$data) && 
                identical(amended$target, original$target) && identical(amended$weights, original$weights), "MIXED_A_TARGET_Y_OR_WEIGHTS_CHANGED")
            positives <- lapply(c("bmi", "fat_mass"), function(id) {
                x <- original$pop$data[[original$a$members[[id]]$column]]
                assert(all(is.finite(x) & x > 0), paste("MIXED_A_NONPOSITIVE_OBSERVATION", id))
                list(indicator_id = id, N = length(x), positive_n = sum(x > 0), min = min(x), max = max(x))
            })
            list(bundle = amended, summary = list(revision = c39_loga_revision(), A_refits = 0L, Y_refits = 0L, 
                reused_log_response = "bmi", restored_original_linear_response = "fat_mass", original_fat_mass_sha256 = c39_loga_hash(old_fat), 
                restored_fat_mass_sha256 = c39_loga_hash(amended$a$models$fat_mass), fat_mass_strict_identical = TRUE, 
                other37_A_strict_identical = TRUE, all_Y_and_non_A_fields_strict_identical = TRUE, weights_strict_identical = TRUE, 
                observed_positive = positives, random_state_unchanged = TRUE, original_bmi_sha256 = c39_loga_hash(old_bmi), 
                reused_bmi_fit_sha256 = c39_loga_hash(new_bmi), no_offset_no_clipping_no_reseeding = TRUE))
        }
    }), envir = context)
}

