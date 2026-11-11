# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

ff_aj_status <- function(require = getOption("ff.y.require_compiled_aj", FALSE)) {
    out <- .ff_aj_backend$metadata
    if (is.null(out)) 
        out <- list(state = .ff_aj_backend$state, ready = FALSE)
    if (require && !isTRUE(.ff_aj_backend$ready)) 
        ff_y_stop("FF_AJ_COMPILED_BACKEND_REQUIRED_NOT_READY", .ff_aj_backend$state)
    out
}

ff_aj_backend <- function(hazard, relative) {
    ff_aj_status(require = TRUE)
    ff_y_assert(is.matrix(hazard) && ncol(hazard) == 2L && is.matrix(relative) && ncol(relative) == 2L && nrow(relative) > 
        0L, "FF_AJ_BACKEND_MATRIX_DIMENSION_INVALID")
    storage.mode(hazard) <- storage.mode(relative) <- "double"
    out <- .Call(.ff_aj_backend$symbol, hazard, relative)
    colnames(out) <- c("disease", "death", "survival")
    ff_y_probability_check(out)
    out
}

ff_aj_compare_fixture <- function(hazard, relative, name) {
    reference <- ff_y_aj_original(hazard, cbind(relative, relative), rep(1L, nrow(hazard)))
    compiled <- ff_aj_backend(hazard, relative)
    difference <- max(abs(reference - compiled), 0)
    data.frame(fixture = name, n = nrow(relative), jumps = nrow(hazard), disease_max_abs = max(abs(reference[, 
        1L] - compiled[, 1L]), 0), death_max_abs = max(abs(reference[, 2L] - compiled[, 2L]), 0), survival_max_abs = max(abs(reference[, 
        3L] - compiled[, 3L]), 0), tolerance = 1e-12, passed = is.finite(difference) && difference <= 1e-12)
}

ff_aj_acceptance <- function() {
    j <- seq_len(127L)
    i <- seq_len(29L)
    h <- cbind(abs(sin(j))/500, abs(cos(j))/700)
    rr <- cbind(exp(seq(-2, 2, length.out = length(i))), exp(seq(1, -1, length.out = length(i))))
    checks <- rbind(ff_aj_compare_fixture(h, rr, "varying_relative_hazards_fractional_weight_like"), ff_aj_compare_fixture(h, 
        matrix(c(1.3, 0.7), 1L), "singleton"), ff_aj_compare_fixture(matrix(numeric(), 0L, 2L), rr, "zero_jumps"), 
        ff_aj_compare_fixture(h/1000, rr * 100, "high_relative_hazards_small_increments"))
    bad <- matrix(c(0.7, 0.4), 1L)
    native_error <- tryCatch({
        ff_aj_backend(bad, matrix(1, 2L, 2L))
        "NO_ERROR"
    }, error = conditionMessage)
    r_error <- tryCatch({
        ff_y_aj_original(bad, matrix(1, 2L, 4L), 1L)
        "NO_ERROR"
    }, error = conditionMessage)
    invalid_pass <- grepl("B_AJ_JUMP_PROBABILITY_INVALID", native_error, fixed = TRUE) && grepl("B_AJ_JUMP_PROBABILITY_INVALID", 
        r_error, fixed = TRUE)
    list(passed = all(checks$passed) && invalid_pass, checks = checks, invalid_jump = list(compiled_error = native_error, 
        reference_error = r_error, passed = invalid_pass), tolerance_basis = "EXACT_FROZEN_R_RECURRENCE_1E_12_ABSOLUTE_EACH_STATE_NO_CLIP_NO_EPSILON", 
        unit_RR_and_weighted_native_validation_still_required = TRUE)
}

ff_aj_initialize <- function(source_root = Sys.getenv("FF_SOURCE_ROOT", ""), work_root = Sys.getenv("FF_WORK", 
    ""), require_compiled = TRUE) {
    ff_y_assert(isTRUE(require_compiled), "FF_AJ_INITIALIZATION_REQUIRES_COMPILED_BACKEND")
    options(ff.y.require_compiled_aj = TRUE)
    if (isTRUE(.ff_aj_backend$ready)) {
        current <- ff_aj_status(TRUE)
        ff_y_assert(nzchar(source_root) && nzchar(work_root) && identical(normalizePath(work_root, mustWork = TRUE), 
            current$work_root) && identical(digest::digest(file = file.path(source_root, "src", "aj.c"), algo = "sha256"), 
            current$source_sha256), "FF_AJ_FIXED_BACKEND_REINITIALIZATION_BINDING_CHANGED")
        return(current)
    }
    ff_y_assert(identical(.ff_aj_backend$state, "NOT_INITIALIZED"), "FF_AJ_FAILED_INITIALIZATION_NOT_SILENTLY_REPLACED_START_NEW_PROCESS")
    .ff_aj_backend$state <- "INITIALIZING"
    tryCatch({
        ff_y_assert(requireNamespace("digest", quietly = TRUE), "FF_AJ_DIGEST_MISSING")
        ff_y_assert(nzchar(source_root) && nzchar(work_root), "FF_AJ_EXPLICIT_SOURCE_ROOT_AND_FF_WORK_REQUIRED")
        source_root <- normalizePath(source_root, mustWork = TRUE)
        ff_y_assert(dir.exists(work_root), "FF_AJ_FF_WORK_DIRECTORY_MISSING")
        work_root <- normalizePath(work_root, mustWork = TRUE)
        ff_y_assert(!identical(work_root, source_root), "FF_AJ_SOFTWARE_MUST_NOT_REPLACE_SOURCE_TREE")
        source <- file.path(source_root, "src", "aj.c")
        ff_y_assert(file.exists(source), "FF_AJ_C_SOURCE_MISSING")
        source_sha <- digest::digest(file = source, algo = "sha256")
        runtime <- list(R = R.version$version.string, platform = R.version$platform, arch = R.version$arch, pointer_size = .Machine$sizeof.pointer, 
            extension = .Platform$dynlib.ext, isolated_build_contract = "R_CMD_SHLIB_V2_EXPLICIT_RUNTIME_COMPILER_PATH", 
            compiler_paths = Sys.which(c("cc", "clang", "gcc", "arm64-apple-darwin20.0.0-clang")))
        runtime_sha <- digest::digest(runtime, algo = "sha256")
        software <- file.path(work_root, "software")
        if (!dir.exists(software)) 
            ff_y_assert(dir.create(software), "FF_AJ_SOFTWARE_DIRECTORY_CREATE_FAILED")
        build <- file.path(software, paste0("aj_", substr(source_sha, 1L, 16L), "_", substr(runtime_sha, 1L, 12L)))
        manifest <- file.path(build, "BUILD_MANIFEST.rds")
        binary <- file.path(build, paste0("ff_aj", .Platform$dynlib.ext))
        copied <- file.path(build, "aj.c")
        if (!dir.exists(build)) {
            ff_y_assert(dir.create(build), "FF_AJ_BUILD_DIRECTORY_CLAIM_FAILED")
            ff_y_assert(file.copy(source, copied, overwrite = FALSE), "FF_AJ_ADDITIVE_SOURCE_COPY_FAILED")
            ff_y_assert(identical(digest::digest(file = copied, algo = "sha256"), source_sha), "FF_AJ_COPIED_SOURCE_HASH_CHANGED")
            previous <- setwd(build)
            on.exit(setwd(previous), add = TRUE)
            executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") 
                "R.exe"
            else "R")
            result <- system2(executable, c("CMD", "SHLIB", "-o", basename(binary), "aj.c"), stdout = TRUE, stderr = TRUE)
            status <- attr(result, "status")
            ff_y_assert((is.null(status) || identical(status, 0L)) && file.exists(binary), "FF_AJ_COMPILATION_FAILED", 
                paste(result, collapse = "\n"))
            record <- list(source_sha256 = source_sha, binary_sha256 = digest::digest(file = binary, algo = "sha256"), 
                runtime_sha256 = runtime_sha, runtime = runtime, compiler_output = result, source_path = source, 
                work_root = work_root, build_directory = build, binary_path = binary, build_scope = "OWN_FF_WORK_SOFTWARE_ONLY_NO_PACKAGE_INSTALL")
            ff_y_assert(!file.exists(manifest), "FF_AJ_BUILD_MANIFEST_ALREADY_EXISTS_REFUSE_OVERWRITE")
            saveRDS(record, manifest)
        }
        else {
            ff_y_assert(file.exists(manifest) && file.exists(binary) && file.exists(copied), "FF_AJ_EXISTING_BUILD_INCOMPLETE_PREPARE_ONCE_PER_HOST_NO_OVERWRITE")
            record <- readRDS(manifest)
            ff_y_assert(identical(record$source_sha256, source_sha) && identical(record$runtime_sha256, runtime_sha) && 
                identical(digest::digest(file = copied, algo = "sha256"), source_sha) && identical(digest::digest(file = binary, 
                algo = "sha256"), record$binary_sha256), "FF_AJ_EXISTING_BUILD_HASH_OR_RUNTIME_MISMATCH")
        }
        record$work_root <- work_root
        dll <- dyn.load(binary, local = TRUE)
        .ff_aj_backend$symbol <- getNativeSymbolInfo("ff_aj_native", PACKAGE = dll, withRegistrationInfo = TRUE)
        .ff_aj_backend$ready <- TRUE
        .ff_aj_backend$state <- "VALIDATING"
        .ff_aj_backend$metadata <- c(record, list(state = "VALIDATING", ready = TRUE))
        acceptance <- ff_aj_acceptance()
        ff_y_assert(isTRUE(acceptance$passed), "FF_AJ_COMPILED_VS_FROZEN_R_VALIDATION_FAILED")
        .ff_aj_backend$state <- "READY_FIXED_FOR_PROCESS"
        .ff_aj_backend$metadata <- c(record, list(state = .ff_aj_backend$state, ready = TRUE, acceptance = acceptance, 
            reference_kernel = ff_y_kernel_provenance(), no_midrun_fallback = TRUE, estimator_change = FALSE))
        ff_aj_status(TRUE)
    }, error = function(e) {
        .ff_aj_backend$ready <- FALSE
        .ff_aj_backend$state <- "FAILED_NO_FALLBACK"
        stop(e)
    })
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        .ff_aj_backend <- new.env(parent = emptyenv())
        .ff_aj_backend$state <- "NOT_INITIALIZED"
        .ff_aj_backend$ready <- FALSE
        ff_aj_status <- function(require = getOption("ff.y.require_compiled_aj", FALSE)) {
            out <- .ff_aj_backend$metadata
            if (is.null(out)) out <- list(state = .ff_aj_backend$state, ready = FALSE)
            if (require && !isTRUE(.ff_aj_backend$ready)) ff_y_stop("FF_AJ_COMPILED_BACKEND_REQUIRED_NOT_READY", 
                .ff_aj_backend$state)
            out
        }
        ff_aj_backend <- function(hazard, relative) {
            ff_aj_status(require = TRUE)
            ff_y_assert(is.matrix(hazard) && ncol(hazard) == 2L && is.matrix(relative) && ncol(relative) == 2L && 
                nrow(relative) > 0L, "FF_AJ_BACKEND_MATRIX_DIMENSION_INVALID")
            storage.mode(hazard) <- storage.mode(relative) <- "double"
            out <- .Call(.ff_aj_backend$symbol, hazard, relative)
            colnames(out) <- c("disease", "death", "survival")
            ff_y_probability_check(out)
            out
        }
        ff_aj_compare_fixture <- function(hazard, relative, name) {
            reference <- ff_y_aj_original(hazard, cbind(relative, relative), rep(1L, nrow(hazard)))
            compiled <- ff_aj_backend(hazard, relative)
            difference <- max(abs(reference - compiled), 0)
            data.frame(fixture = name, n = nrow(relative), jumps = nrow(hazard), disease_max_abs = max(abs(reference[, 
                1L] - compiled[, 1L]), 0), death_max_abs = max(abs(reference[, 2L] - compiled[, 2L]), 0), survival_max_abs = max(abs(reference[, 
                3L] - compiled[, 3L]), 0), tolerance = 1e-12, passed = is.finite(difference) && difference <= 1e-12)
        }
        ff_aj_acceptance <- function() {
            j <- seq_len(127L)
            i <- seq_len(29L)
            h <- cbind(abs(sin(j))/500, abs(cos(j))/700)
            rr <- cbind(exp(seq(-2, 2, length.out = length(i))), exp(seq(1, -1, length.out = length(i))))
            checks <- rbind(ff_aj_compare_fixture(h, rr, "varying_relative_hazards_fractional_weight_like"), ff_aj_compare_fixture(h, 
                matrix(c(1.3, 0.7), 1L), "singleton"), ff_aj_compare_fixture(matrix(numeric(), 0L, 2L), rr, "zero_jumps"), 
                ff_aj_compare_fixture(h/1000, rr * 100, "high_relative_hazards_small_increments"))
            bad <- matrix(c(0.7, 0.4), 1L)
            native_error <- tryCatch({
                ff_aj_backend(bad, matrix(1, 2L, 2L))
                "NO_ERROR"
            }, error = conditionMessage)
            r_error <- tryCatch({
                ff_y_aj_original(bad, matrix(1, 2L, 4L), 1L)
                "NO_ERROR"
            }, error = conditionMessage)
            invalid_pass <- grepl("B_AJ_JUMP_PROBABILITY_INVALID", native_error, fixed = TRUE) && grepl("B_AJ_JUMP_PROBABILITY_INVALID", 
                r_error, fixed = TRUE)
            list(passed = all(checks$passed) && invalid_pass, checks = checks, invalid_jump = list(compiled_error = native_error, 
                reference_error = r_error, passed = invalid_pass), tolerance_basis = "EXACT_FROZEN_R_RECURRENCE_1E_12_ABSOLUTE_EACH_STATE_NO_CLIP_NO_EPSILON", 
                unit_RR_and_weighted_native_validation_still_required = TRUE)
        }
        ff_aj_initialize <- function(source_root = Sys.getenv("FF_SOURCE_ROOT", ""), work_root = Sys.getenv("FF_WORK", 
            ""), require_compiled = TRUE) {
            ff_y_assert(isTRUE(require_compiled), "FF_AJ_INITIALIZATION_REQUIRES_COMPILED_BACKEND")
            options(ff.y.require_compiled_aj = TRUE)
            if (isTRUE(.ff_aj_backend$ready)) {
                current <- ff_aj_status(TRUE)
                ff_y_assert(nzchar(source_root) && nzchar(work_root) && identical(normalizePath(work_root, mustWork = TRUE), 
                  current$work_root) && identical(digest::digest(file = file.path(source_root, "src", "aj.c"), 
                  algo = "sha256"), current$source_sha256), "FF_AJ_FIXED_BACKEND_REINITIALIZATION_BINDING_CHANGED")
                return(current)
            }
            ff_y_assert(identical(.ff_aj_backend$state, "NOT_INITIALIZED"), "FF_AJ_FAILED_INITIALIZATION_NOT_SILENTLY_REPLACED_START_NEW_PROCESS")
            .ff_aj_backend$state <- "INITIALIZING"
            tryCatch({
                ff_y_assert(requireNamespace("digest", quietly = TRUE), "FF_AJ_DIGEST_MISSING")
                ff_y_assert(nzchar(source_root) && nzchar(work_root), "FF_AJ_EXPLICIT_SOURCE_ROOT_AND_FF_WORK_REQUIRED")
                source_root <- normalizePath(source_root, mustWork = TRUE)
                ff_y_assert(dir.exists(work_root), "FF_AJ_FF_WORK_DIRECTORY_MISSING")
                work_root <- normalizePath(work_root, mustWork = TRUE)
                ff_y_assert(!identical(work_root, source_root), "FF_AJ_SOFTWARE_MUST_NOT_REPLACE_SOURCE_TREE")
                source <- file.path(source_root, "src", "aj.c")
                ff_y_assert(file.exists(source), "FF_AJ_C_SOURCE_MISSING")
                source_sha <- digest::digest(file = source, algo = "sha256")
                runtime <- list(R = R.version$version.string, platform = R.version$platform, arch = R.version$arch, 
                  pointer_size = .Machine$sizeof.pointer, extension = .Platform$dynlib.ext, isolated_build_contract = "R_CMD_SHLIB_V2_EXPLICIT_RUNTIME_COMPILER_PATH", 
                  compiler_paths = Sys.which(c("cc", "clang", "gcc", "arm64-apple-darwin20.0.0-clang")))
                runtime_sha <- digest::digest(runtime, algo = "sha256")
                software <- file.path(work_root, "software")
                if (!dir.exists(software)) ff_y_assert(dir.create(software), "FF_AJ_SOFTWARE_DIRECTORY_CREATE_FAILED")
                build <- file.path(software, paste0("aj_", substr(source_sha, 1L, 16L), "_", substr(runtime_sha, 
                  1L, 12L)))
                manifest <- file.path(build, "BUILD_MANIFEST.rds")
                binary <- file.path(build, paste0("ff_aj", .Platform$dynlib.ext))
                copied <- file.path(build, "aj.c")
                if (!dir.exists(build)) {
                  ff_y_assert(dir.create(build), "FF_AJ_BUILD_DIRECTORY_CLAIM_FAILED")
                  ff_y_assert(file.copy(source, copied, overwrite = FALSE), "FF_AJ_ADDITIVE_SOURCE_COPY_FAILED")
                  ff_y_assert(identical(digest::digest(file = copied, algo = "sha256"), source_sha), "FF_AJ_COPIED_SOURCE_HASH_CHANGED")
                  previous <- setwd(build)
                  on.exit(setwd(previous), add = TRUE)
                  executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")
                  result <- system2(executable, c("CMD", "SHLIB", "-o", basename(binary), "aj.c"), stdout = TRUE, 
                    stderr = TRUE)
                  status <- attr(result, "status")
                  ff_y_assert((is.null(status) || identical(status, 0L)) && file.exists(binary), "FF_AJ_COMPILATION_FAILED", 
                    paste(result, collapse = "\n"))
                  record <- list(source_sha256 = source_sha, binary_sha256 = digest::digest(file = binary, algo = "sha256"), 
                    runtime_sha256 = runtime_sha, runtime = runtime, compiler_output = result, source_path = source, 
                    work_root = work_root, build_directory = build, binary_path = binary, build_scope = "OWN_FF_WORK_SOFTWARE_ONLY_NO_PACKAGE_INSTALL")
                  ff_y_assert(!file.exists(manifest), "FF_AJ_BUILD_MANIFEST_ALREADY_EXISTS_REFUSE_OVERWRITE")
                  saveRDS(record, manifest)
                } else {
                  ff_y_assert(file.exists(manifest) && file.exists(binary) && file.exists(copied), "FF_AJ_EXISTING_BUILD_INCOMPLETE_PREPARE_ONCE_PER_HOST_NO_OVERWRITE")
                  record <- readRDS(manifest)
                  ff_y_assert(identical(record$source_sha256, source_sha) && identical(record$runtime_sha256, runtime_sha) && 
                    identical(digest::digest(file = copied, algo = "sha256"), source_sha) && identical(digest::digest(file = binary, 
                    algo = "sha256"), record$binary_sha256), "FF_AJ_EXISTING_BUILD_HASH_OR_RUNTIME_MISMATCH")
                }
                record$work_root <- work_root
                dll <- dyn.load(binary, local = TRUE)
                .ff_aj_backend$symbol <- getNativeSymbolInfo("ff_aj_native", PACKAGE = dll, withRegistrationInfo = TRUE)
                .ff_aj_backend$ready <- TRUE
                .ff_aj_backend$state <- "VALIDATING"
                .ff_aj_backend$metadata <- c(record, list(state = "VALIDATING", ready = TRUE))
                acceptance <- ff_aj_acceptance()
                ff_y_assert(isTRUE(acceptance$passed), "FF_AJ_COMPILED_VS_FROZEN_R_VALIDATION_FAILED")
                .ff_aj_backend$state <- "READY_FIXED_FOR_PROCESS"
                .ff_aj_backend$metadata <- c(record, list(state = .ff_aj_backend$state, ready = TRUE, acceptance = acceptance, 
                  reference_kernel = ff_y_kernel_provenance(), no_midrun_fallback = TRUE, estimator_change = FALSE))
                ff_aj_status(TRUE)
            }, error = function(e) {
                .ff_aj_backend$ready <- FALSE
                .ff_aj_backend$state <- "FAILED_NO_FALLBACK"
                stop(e)
            })
        }
    }), envir = context)
}

