# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

required_packages <- c("dplyr", "tidyr", "ggplot2", "ggrepel", "digest", "jsonlite", "detectseparation")

compat <- NULL

input_snapshot <- function() {
    z <- file.info(paths)
    data.frame(input = names(paths), path = paths, bytes = z$size, mtime = as.numeric(z$mtime), row.names = NULL)
}

stored_version <- function(name) {
    x <- attr(phenos, name, exact = TRUE)
    if (is.character(x) && length(x) == 1L) 
        x
    else NA_character_
}

write_new_csv <- function(x, name) {
    p <- file.path(out, "raw", name)
    if (file.exists(p)) 
        stop("Refuse existing output: ", p)
    utils::write.csv(x, p, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}

state <- function(stage, done = 0L, total = length(codes), M = NA_integer_) {
    jsonlite::write_json(list(state = "RUNNING", stage = stage, branch = branch, completed = done, total = total, 
        M = M, updated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), file.path(out, "progress.json"), auto_unbox = TRUE, 
        pretty = TRUE)
    cat("[PheWAS] ", branch, " ", stage, " ", done, "/", total, "\n", sep = "")
    flush.console()
}

run_registered_models <- function() {
    executor <- create_linear_executor(dat, package, execution_workers, branch)
    on.exit(executor$close(), add = TRUE)
    jsonlite::write_json(executor$info(), file.path(out, "raw", "execution_start.json"), auto_unbox = TRUE, pretty = TRUE)
    for (first in seq.int(1L, length(codes), by = execution_workers)) {
        batch <- seq.int(first, min(length(codes), first + execution_workers - 1L))
        state(paste0("DISPATCH ", paste(codes[batch], collapse = ",")), first - 1L, length(codes), M)
        tasks <- lapply(batch, function(i) {
            fit_allowed <- registry$planned[i] && meta$metadata_verified[i]
            list(index = i, y = phenos[[codes[i]]][shared], exposure = branch, phecode = codes[i], sex_rule = meta$sex_rule[i], 
                min_cases = if (fit_allowed) linear_config$min_cases else Inf, min_controls = if (fit_allowed) linear_config$min_controls else Inf)
        })
        got <- executor$run(tasks)
        if (length(got) != length(batch)) 
            stop("Parallel batch lost rows; no final q or completion marker")
        for (j in seq_along(batch)) {
            i <- batch[j]
            r <- got[[j]]
            y <- tasks[[j]]$y
            if (!identical(as.character(r$phecode), codes[i])) 
                stop("Parallel task order/key conflict: ", codes[i])
            if (!(registry$planned[i] && meta$metadata_verified[i])) {
                r$reason <- paste("Original case/control threshold", linear_config$min_cases, linear_config$min_controls)
                if (registry$cases[i] >= linear_config$min_cases && registry$noncases[i] >= linear_config$min_controls) {
                  r$status <- "NO_EXPOSURE_VARIATION"
                  r$reason <- "Birth-order target is constant in the disease-complete sample"
                }
                if (!branch_ok) {
                  r$status <- "BRANCH_INSUFFICIENT_SHARED_DATA"
                  r$reason <- "Original branch N>=500 or exposure variation gate not met"
                }
                else if (!meta$metadata_verified[i]) {
                  r$status <- "PHECODE_METADATA_UNVERIFIED"
                  r$reason <- "Exact name/category/sex rule not found in pinned package; row and otherwise eligible planned M retained"
                }
            }
            if (r$N != registry$N[i] || r$cases != registry$cases[i] || r$noncases != registry$noncases[i] || r$actual_formula != 
                registry$actual_formula[i]) 
                stop("Prefit registry/sample drift: ", codes[i])
            r <- cbind(r, meta[i, setdiff(names(meta), "phecode"), drop = FALSE])
            r$planned <- registry$planned[i]
            r$run_id <- owner$run_id
            r$model_id <- "Primary_linear_baseline_age_v1"
            r$exposure <- r$exposure_id
            r$phenotype <- r$phecode
            r$snp <- r$exposure_id
            r$type <- "logistic"
            r$n_total <- r$N
            r$n_cases <- r$cases
            r$n_controls <- r$noncases
            r$p_raw <- r$p
            r$execution_backend <- execution_backend
            r$execution_workers <- execution_workers
            r$sample_id_digest <- digest::digest(sort(as.character(linear_sample(dat, y, branch, meta$sex_rule[i])$data$eid)), 
                algo = "sha256")
            results[[i]] <<- r
        }
        partial <- adjust_branch(dplyr::bind_rows(results[seq_len(max(batch))]), M, complete = FALSE)
        utils::write.csv(partial, file.path(out, "raw", "phewas_partial.csv"), row.names = FALSE, na = "NA")
        state(paste0("BATCH_COMPLETED through ", codes[max(batch)]), max(batch), length(codes), M)
    }
    last_info <- executor$info()
    executor$close()
    jsonlite::write_json(list(before_close = last_info, after_close = executor$info()), file.path(out, "raw", "execution_end.json"), 
        auto_unbox = TRUE, pretty = TRUE)
    invisible(NULL)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        args <- commandArgs(TRUE)
        if (length(args) < 2L) stop("Expected MODE PACKAGE [BRANCH OUTPUT]")
        mode <- args[1]
        package <- normalizePath(args[2], mustWork = TRUE)
        if (!mode %in% c("check", "inspect", "run")) stop("Only check/inspect/run are supported")
        source(file.path(package, "config.R"))
        required_packages <- c("dplyr", "tidyr", "ggplot2", "ggrepel", "digest", "jsonlite", "detectseparation")
        for (p in required_packages) {
            if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p, "; run the isolated dependency setup, not a shared-library update")
            cat("PACKAGE ", p, " ", as.character(utils::packageVersion(p)), "\n", sep = "")
        }
        source(file.path(package, "code", "linear_functions.R"))
        source(file.path(package, "code", "parallel_executor.R"))
        worker_text <- Sys.getenv("BORDER_LINEAR_WORKERS", as.character(linear_config$execution_default_workers))
        if (!grepl("^[1-8]$", worker_text)) stop("Execution workers must be an integer from 1 to 8; no automatic core overcommit")
        execution_workers <- as.integer(worker_text)
        execution_backend <- if (execution_workers > 1L) Sys.getenv("BORDER_LINEAR_EXECUTOR", "PSOCK") else "serial"
        if (!execution_backend %in% c("PSOCK", "FORK_PIPE", "serial")) stop("Unreviewed execution backend")
        cat("EXECUTION_BACKEND=", execution_backend, "; WORKERS=", execution_workers, "; native threads per worker=1\n", 
            sep = "")
        indir <- Sys.getenv("PHEDIR", linear_config$phenotype_root_default)
        script_root <- Sys.getenv("SCRIPT_ROOT", birth_private_setting("PRIVATE_SCRIPTS"))
        source_root <- Sys.getenv("BORDER_SOURCE_ROOT", file.path(script_root, "border"))
        helper_dir <- Sys.getenv("HELPER_DIR", file.path(script_root, "0f"))
        cat("INPUT_ROOT=", indir, "\nSCRIPT_ROOT=", script_root, "\nSOURCE_ROOT=", source_root, "\nHELPER_DIR=", 
            helper_dir, "\n", sep = "")
        paths <- c(all = Sys.getenv("ALL_RDS", file.path(indir, "Rdata", "all.rds")), phenotypes = Sys.getenv("PHEWAS_RDS", 
            file.path(indir, "Rdata", "PheWAS.rds")))
        cat("ALL_RDS=", paths[["all"]], "\nPHEWAS_RDS=", paths[["phenotypes"]], "\n", sep = "")
        if (any(!file.exists(paths))) stop("Missing exact input paths: ", paste(paths[!file.exists(paths)], collapse = "; "))
        compat <- NULL
        compat_path <- Sys.getenv("BORDER_ALL_COMPAT_RECORD", "")
        if (nzchar(compat_path)) {
            if (!file.exists(compat_path) || digest::digest(file = compat_path, algo = "sha256") != linear_config$adopted_compat_record_sha256) stop("Known full adoption-compatible input record missing/changed")
            compat <- jsonlite::read_json(compat_path, simplifyVector = TRUE)
            info <- file.info(paths[["all"]])
            if (!isTRUE(compat$original_columns_identical) || compat$N != 501937 || compat$source_field != "p1767_i0" || 
                compat$source_instance != 0 || compat$added_column != "adopted" || compat$unmatched_N != 1 || info$size != 
                linear_config$adopted_compat_all_bytes || floor(as.numeric(info$mtime)) != linear_config$adopted_compat_all_mtime) stop("Full compatible input provenance/size/mtime conflicts; no rebuilding or fill-zero")
            cat("FULL_ALL_INPUT=Huang authenticated adoption-only compatible copy; original full population, no education/TDI filtering\n")
        }
        pins <- linear_config$source_sha256
        for (n in names(pins)) {
            f <- file.path(source_root, "f", n)
            if (!file.exists(f) || digest::digest(file = f, algo = "sha256") != pins[[n]]) stop("Huang source pin missing/changed: ", 
                f)
        }
        if (digest::digest(file = file.path(package, "code", "source0common.R"), algo = "sha256") != pins[["0.common.R"]]) stop("Immutable helper SHA conflict")
        plot_helper <- file.path(helper_dir, "plotting.R")
        if (!file.exists(plot_helper) || digest::digest(file = plot_helper, algo = "sha256") != linear_config$plotting_sha256) stop("Shared plotting helper missing/changed")
        for (n in names(linear_config$metadata_sha256)) {
            f <- file.path(package, "metadata", n)
            if (!file.exists(f) || digest::digest(file = f, algo = "sha256") != linear_config$metadata_sha256[[n]]) stop("Pinned original PheWAS metadata snapshot conflict: ", 
                f)
        }
        dictionary <- utils::read.csv(file.path(package, "metadata", "pheinfo.csv"), colClasses = "character", 
            check.names = FALSE)
        restrictions <- utils::read.csv(file.path(package, "metadata", "sex_restriction.csv"), colClasses = "character", 
            check.names = FALSE)
        for (n in c("male_only", "female_only")) {
            if (any(!tolower(restrictions[[n]]) %in% c("true", "false"))) stop("Invalid pinned restriction encoding")
            restrictions[[n]] <- tolower(restrictions[[n]]) == "true"
        }
        cat("MAPPING_VERSION=", linear_config$mapping_version, "; original WinPC package tables, unchanged\n", 
            sep = "")
        if (mode == "check") {
            files <- c(file.path(package, "config.R"), list.files(file.path(package, "code"), pattern = "[.]R$", 
                full.names = TRUE), list.files(file.path(package, "tests"), pattern = "[.]R$", full.names = TRUE))
            for (f in files) invisible(parse(f))
            print(vapply(linear_exposures, function(x) paste(deparse(linear_formulas(x)), collapse = " "), character(1)))
            cat("CHECK_OK; no RDS read, no participant model fitted\n")
            quit(status = 0)
        }
        if (mode == "run") {
            if (length(args) != 4L) stop("run requires exactly BRANCH OUTPUT")
            branch <- args[3]
            out <- normalizePath(args[4], mustWork = TRUE)
            if (!branch %in% linear_exposures || !file.exists(file.path(out, "owner.json"))) stop("Explicit branch/owned output missing")
            owner <- jsonlite::read_json(file.path(out, "owner.json"), simplifyVector = TRUE)
            if (owner$package != linear_config$version || owner$branch != branch || owner$output != out) stop("Wrong output ownership")
            if (file.exists(file.path(out, "complete.json"))) stop("Completed result cannot be overwritten")
            if (isTRUE(owner$synthetic_only)) {
                fixture <- file.path(dirname(out), "inputs", "Rdata")
                if (!identical(normalizePath(paths[["all"]], mustWork = TRUE), normalizePath(file.path(fixture, 
                  "all.rds"), mustWork = TRUE)) || !identical(normalizePath(paths[["phenotypes"]], mustWork = TRUE), 
                  normalizePath(file.path(fixture, "PheWAS.rds"), mustWork = TRUE))) stop("Synthetic owner must use its generated fixture, never real UKB inputs")
            }
        }
        input_snapshot <- function() {
            z <- file.info(paths)
            data.frame(input = names(paths), path = paths, bytes = z$size, mtime = as.numeric(z$mtime), row.names = NULL)
        }
        before <- input_snapshot()
        all <- readRDS(paths[["all"]])
        if (!is.data.frame(all) || anyDuplicated(names(all))) stop("all.rds must be a data frame with unique columns")
        if (!is.null(compat) && nrow(all) != compat$N) stop("Full compatible baseline population changed")
        map <- c(eid = linear_config$input_columns[["id"]], age = linear_config$input_columns[["age"]], sex = linear_config$input_columns[["sex"]], 
            ethnicity = linear_config$input_columns[["ethnicity"]], center = linear_config$input_columns[["center"]], 
            date_attend = linear_config$input_columns[["baseline_date"]], cnt_older = linear_config$count_columns[["older"]], 
            cnt_brother = linear_config$count_columns[["brothers"]], cnt_sister = linear_config$count_columns[["sisters"]], 
            adopted = linear_config$count_columns[["adopted"]])
        if (anyDuplicated(unname(map))) stop("Input field mapping reuses a source column")
        if (any(!map %in% names(all))) stop("Required fields missing: ", paste(map[!map %in% names(all)], collapse = "; "))
        if (linear_config$ethnicity_encoding != "UKB21000_Coding1001_numeric" || linear_config$sex_encoding != 
            "UKB31_0_female_1_male") stop("Unreviewed input encoding declared")
        raw <- as.data.frame(all[, unname(map), drop = FALSE])
        names(raw) <- names(map)
        field_map <- data.frame(model_column = names(map), source_column = unname(map), source = paths[["all"]], 
            class = vapply(raw, function(x) paste(class(x), collapse = ";"), character(1)), stringsAsFactors = FALSE)
        field_map$compatibility_record <- if (nzchar(compat_path)) compat_path else NA_character_
        field_map$adoption_source_field <- ifelse(field_map$model_column == "adopted" & !is.null(compat), "p1767_i0", 
            "")
        rm(all)
        invisible(gc())
        dat <- linear_covars(make_birth_order(raw))
        if (all(is.na(dat$ethnicity)) || all(is.na(dat$sex)) || all(is.na(dat$assessment_year))) stop("A required categorical adapter produced no valid values; encoding/date must be reviewed, not omitted")
        audit <- data.frame(metric = c("raw_all_rows", "known_adopted_excluded", "unknown_adoption_retained_if_otherwise_eligible", 
            "inconsistent_rank", "age_missing", "sex_missing_or_unmapped", "ethnicity_missing_or_unmapped", "center_missing", 
            "baseline_date_missing_or_invalid", paste0(linear_exposures, "_eligible")), value = c(nrow(dat), sum(dat$adopted_num == 
            1, na.rm = TRUE), sum(is.na(dat$adopted_num) & dat$eligible_border), sum(dat$impossible_rank, na.rm = TRUE), 
            sum(is.na(dat$age)), sum(is.na(dat$sex)), sum(is.na(dat$ethnicity)), sum(is.na(dat$center)), sum(is.na(dat$assessment_year)), 
            vapply(linear_exposures, function(x) sum(!is.na(dat[[x]])), integer(1))))
        rm(raw)
        invisible(gc())
        phenos <- readRDS(paths[["phenotypes"]])
        if (!is.data.frame(phenos) || !("id" %in% names(phenos)) || anyDuplicated(names(phenos))) stop("PheWAS matrix must have id and unique character column keys")
        phenos$id <- as.character(phenos$id)
        if (anyNA(phenos$id) || any(!nzchar(phenos$id)) || anyDuplicated(phenos$id)) stop("PheWAS IDs missing/duplicated; no deduplication")
        codes <- setdiff(names(phenos), "id")
        if (!length(codes)) stop("Empty PheWAS disease spectrum")
        if (any(!vapply(phenos[codes], function(v) all(is.na(v) | v %in% c(0, 1)), logical(1)))) stop("PheWAS states must be 0/1/NA; no automatic recoding of unknown values")
        if (anyDuplicated(as.character(dictionary$phecode)) || anyDuplicated(as.character(restrictions$phecode))) stop("Ambiguous package PheCode metadata")
        ix <- match(codes, as.character(dictionary$phecode))
        sx <- match(codes, as.character(restrictions$phecode))
        meta <- data.frame(phecode = codes, description = dictionary$description[ix], category = dictionary$group[ix], 
            groupnum = dictionary$groupnum[ix], color = dictionary$color[ix], metadata_verified = !is.na(ix) & 
                !is.na(sx), sex_rule = "both", stringsAsFactors = FALSE)
        male <- !is.na(sx) & restrictions$male_only[sx] %in% TRUE
        female <- !is.na(sx) & restrictions$female_only[sx] %in% TRUE
        if (any(male & female)) stop("Conflicting official sex restriction")
        meta$sex_rule[male] <- "male"
        meta$sex_rule[female] <- "female"
        meta$mapping_version <- linear_config$mapping_version
        meta$definition_version <- linear_config$definition_version
        stored_version <- function(name) {
            x <- attr(phenos, name, exact = TRUE)
            if (is.character(x) && length(x) == 1L) x else NA_character_
        }
        meta$stored_definition_version <- stored_version("definition_version")
        meta$stored_mapping_version <- stored_version("mapping_version")
        meta$original_source <- paste0(paths[["phenotypes"]], "::", meta$phecode)
        if (mode == "inspect") {
            print(field_map)
            print(audit)
            cat("PHENOTYPE_ROWS=", nrow(phenos), "\nPHENOTYPE_CODES=", length(codes), "\nUNMAPPED_METADATA_CODES=", 
                sum(!meta$metadata_verified), "\n", sep = "")
            print(meta[!meta$metadata_verified, c("phecode", "description", "metadata_verified")])
            if (!identical(before, input_snapshot())) stop("Source changed during inspection")
            cat("INSPECT_OK; aggregate fields/counts only; no participant model fitted\n")
            quit(status = 0)
        }
        write_new_csv <- function(x, name) {
            p <- file.path(out, "raw", name)
            if (file.exists(p)) stop("Refuse existing output: ", p)
            utils::write.csv(x, p, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
        }
        state <- function(stage, done = 0L, total = length(codes), M = NA_integer_) {
            jsonlite::write_json(list(state = "RUNNING", stage = stage, branch = branch, completed = done, total = total, 
                M = M, updated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), file.path(out, "progress.json"), auto_unbox = TRUE, 
                pretty = TRUE)
            cat("[PheWAS] ", branch, " ", stage, " ", done, "/", total, "\n", sep = "")
            flush.console()
        }
        write_new_csv(field_map, "field_mapping.csv")
        write_new_csv(audit, "cohort_qc.csv")
        write_new_csv(before, "input_snapshot.csv")
        write_new_csv(meta, "phecode_metadata.csv")
        writeLines(capture.output(sessionInfo()), file.path(out, "raw", "sessionInfo.txt"))
        shared <- match(as.character(dat$eid), phenos$id)
        has_shared <- !is.na(shared)
        unmatched <- sum(!has_shared & !is.na(dat[[branch]]))
        dat <- dat[has_shared, , drop = FALSE]
        shared <- shared[has_shared]
        required <- setdiff(all.vars(linear_formulas(branch)), "y")
        base_complete <- !is.na(dat[[branch]]) & complete.cases(dat[, required, drop = FALSE])
        branch_N <- sum(base_complete)
        branch_ok <- branch_N >= linear_config$branch_min_shared_predictor_complete && length(unique(dat[[branch]][base_complete])) >= 
            2L
        write_new_csv(data.frame(branch = branch, N_predictor_complete_shared = branch_N, exposure_eligible_unmatched_phenotype_ID = unmatched, 
            branch_minimum = linear_config$branch_min_shared_predictor_complete, branch_eligible = branch_ok), 
            "branch_sample.csv")
        plan <- vector("list", length(codes))
        state("REGISTERING_ALL_PHECODES")
        for (i in seq_along(codes)) {
            s <- linear_sample(dat, phenos[[codes[i]]][shared], branch, meta$sex_rule[i])
            eligible <- branch_ok && s$counts$cases >= linear_config$min_cases && s$counts$noncases >= linear_config$min_controls && 
                length(unique(s$data[[branch]])) >= 2L
            plan[[i]] <- data.frame(meta[i, ], planned = eligible, N = s$counts$N, cases = s$counts$cases, noncases = s$counts$noncases, 
                actual_formula = paste(deparse(s$formula), collapse = " "), removed_terms = paste(s$removed_terms, 
                  collapse = ";"))
        }
        registry <- dplyr::bind_rows(plan)
        M <- sum(registry$planned)
        write_new_csv(registry, "planned_tests.csv")
        rm(plan)
        invisible(gc())
        results <- vector("list", length(codes))
        run_registered_models <- function() {
            executor <- create_linear_executor(dat, package, execution_workers, branch)
            on.exit(executor$close(), add = TRUE)
            jsonlite::write_json(executor$info(), file.path(out, "raw", "execution_start.json"), auto_unbox = TRUE, 
                pretty = TRUE)
            for (first in seq.int(1L, length(codes), by = execution_workers)) {
                batch <- seq.int(first, min(length(codes), first + execution_workers - 1L))
                state(paste0("DISPATCH ", paste(codes[batch], collapse = ",")), first - 1L, length(codes), M)
                tasks <- lapply(batch, function(i) {
                  fit_allowed <- registry$planned[i] && meta$metadata_verified[i]
                  list(index = i, y = phenos[[codes[i]]][shared], exposure = branch, phecode = codes[i], sex_rule = meta$sex_rule[i], 
                    min_cases = if (fit_allowed) linear_config$min_cases else Inf, min_controls = if (fit_allowed) linear_config$min_controls else Inf)
                })
                got <- executor$run(tasks)
                if (length(got) != length(batch)) stop("Parallel batch lost rows; no final q or completion marker")
                for (j in seq_along(batch)) {
                  i <- batch[j]
                  r <- got[[j]]
                  y <- tasks[[j]]$y
                  if (!identical(as.character(r$phecode), codes[i])) stop("Parallel task order/key conflict: ", 
                    codes[i])
                  if (!(registry$planned[i] && meta$metadata_verified[i])) {
                    r$reason <- paste("Original case/control threshold", linear_config$min_cases, linear_config$min_controls)
                    if (registry$cases[i] >= linear_config$min_cases && registry$noncases[i] >= linear_config$min_controls) {
                      r$status <- "NO_EXPOSURE_VARIATION"
                      r$reason <- "Birth-order target is constant in the disease-complete sample"
                    }
                    if (!branch_ok) {
                      r$status <- "BRANCH_INSUFFICIENT_SHARED_DATA"
                      r$reason <- "Original branch N>=500 or exposure variation gate not met"
                    } else if (!meta$metadata_verified[i]) {
                      r$status <- "PHECODE_METADATA_UNVERIFIED"
                      r$reason <- "Exact name/category/sex rule not found in pinned package; row and otherwise eligible planned M retained"
                    }
                  }
                  if (r$N != registry$N[i] || r$cases != registry$cases[i] || r$noncases != registry$noncases[i] || 
                    r$actual_formula != registry$actual_formula[i]) stop("Prefit registry/sample drift: ", codes[i])
                  r <- cbind(r, meta[i, setdiff(names(meta), "phecode"), drop = FALSE])
                  r$planned <- registry$planned[i]
                  r$run_id <- owner$run_id
                  r$model_id <- "Primary_linear_baseline_age_v1"
                  r$exposure <- r$exposure_id
                  r$phenotype <- r$phecode
                  r$snp <- r$exposure_id
                  r$type <- "logistic"
                  r$n_total <- r$N
                  r$n_cases <- r$cases
                  r$n_controls <- r$noncases
                  r$p_raw <- r$p
                  r$execution_backend <- execution_backend
                  r$execution_workers <- execution_workers
                  r$sample_id_digest <- digest::digest(sort(as.character(linear_sample(dat, y, branch, meta$sex_rule[i])$data$eid)), 
                    algo = "sha256")
                  results[[i]] <<- r
                }
                partial <- adjust_branch(dplyr::bind_rows(results[seq_len(max(batch))]), M, complete = FALSE)
                utils::write.csv(partial, file.path(out, "raw", "phewas_partial.csv"), row.names = FALSE, na = "NA")
                state(paste0("BATCH_COMPLETED through ", codes[max(batch)]), max(batch), length(codes), M)
            }
            last_info <- executor$info()
            executor$close()
            jsonlite::write_json(list(before_close = last_info, after_close = executor$info()), file.path(out, 
                "raw", "execution_end.json"), auto_unbox = TRUE, pretty = TRUE)
            invisible(NULL)
        }
        run_registered_models()
        res <- dplyr::bind_rows(results)
        if (nrow(res) != length(codes) || anyDuplicated(res$phecode) || !identical(res$phecode, codes) || sum(res$planned) != 
            M) stop("Incomplete/duplicated result coverage; no final q")
        if (any(is.finite(res$p) & !res$planned)) stop("An unregistered model produced P; no final q")
        if (!identical(before, input_snapshot())) stop("Input size/mtime changed; no final q and partial results preserved")
        res <- adjust_branch(res, M, complete = TRUE)
        res$FDR <- res$q
        res$bonferroni_p <- res$bonf
        res$BH_q <- res$q
        res$bonferroni_threshold <- if (M > 0) 0.05/M else NA_real_
        res$p_underflow <- !is.na(res$p) & res$p == 0
        res$logp <- -log10(pmax(res$p, .Machine$double.xmin))
        res$multiplicity_scope <- "separate branch; not four-branch joint correction"
        write_new_csv(res, "phewas.csv")
        saveRDS(res, file.path(out, "raw", "phewas.rds"))
        state("RENDER_SAVED_RESULTS", length(codes), length(codes), M)
        source(file.path(package, "code", "phewas_plot.R"))
        plot_status <- tryCatch({
            as.character(render_linear_phewas(res, dat, base_complete, branch, helper_dir, out, dictionary))
        }, error = function(e) paste("PLOT_FAILED_RESULT_TABLE_RETAINED:", conditionMessage(e)))
        writeLines(plot_status, file.path(out, "raw", "plot_status.txt"))
        jsonlite::write_json(list(state = "PHEWAS_DONE", branch = branch, registered_PheCodes = length(codes), 
            M = M, estimable = sum(is.finite(res$p)), planned_failed = sum(res$planned & !is.finite(res$p)), all_codes_retained = TRUE, 
            input_size_mtime_unchanged = TRUE, plot_status = plot_status, real_scope = "hospital ever-recorded disease; OR, not incident HR", 
            definition_limit = linear_config$phenotype_limit, finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), 
            file.path(out, "complete.json"), auto_unbox = TRUE, pretty = TRUE)
        cat("PHEWAS_DONE ", branch, "; M=", M, "; saved all ", length(codes), " rows; no Cox/omics executed\n", 
            sep = "")
    }), envir = context)
}

