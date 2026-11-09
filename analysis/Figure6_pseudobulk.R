# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

env_or <- function(name, default) {
    value <- Sys.getenv(name, unset = "")
    if (nzchar(value)) 
        value
    else default
}

resolve_path <- function(path, must_exist = FALSE) {
    path <- path.expand(path)
    if (!grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)) 
        path <- file.path(getwd(), path)
    normalizePath(path, winslash = "/", mustWork = must_exist)
}

fail <- function(...) stop(paste0(...), call. = FALSE)

assert <- function(condition, ...) if (!isTRUE(condition)) fail(...)

path_exists <- function(path) {
    link <- Sys.readlink(path)
    file.exists(path) || dir.exists(path) || (!is.na(link) && nzchar(link))
}

write_csv <- function(x, path) {
    assert(!path_exists(path), "Refusing to overwrite: ", path)
    write.csv(x, path, row.names = FALSE, na = "")
}

write_text <- function(x, path) {
    assert(!path_exists(path), "Refusing to overwrite: ", path)
    writeLines(x, path, useBytes = TRUE)
}

bind_rows <- function(x) {
    x <- Filter(function(z) !is.null(z) && nrow(z) > 0L, x)
    if (!length(x)) 
        return(data.frame())
    do.call(rbind, x)
}

bh_fixed <- function(p, family_n) {
    out <- rep(NA_real_, length(p))
    valid <- is.finite(p) & p >= 0 & p <= 1
    if (any(valid)) 
        out[valid] <- p.adjust(p[valid], method = "BH", n = family_n)
    out
}

needed <- c("Matrix", "edgeR", "locfit", "statmod")

candidates <- c("FCRL6", "GZMA", "GZMH", "KLRD1")

primary_types <- c("NK", "CD8_T")

min_cells <- 20L

min_donors <- 3L

input_names <- c("pseudobulk_counts.mtx.gz", "pseudobulk_genes.tsv", "pseudobulk_samples.tsv")

empty_all <- function() data.frame(dataset = character(), cell_type = character(), contrast = character(), gene_id = character(), 
    gene = character(), log2FC = numeric(), SE = numeric(), CI_low = numeric(), CI_high = numeric(), LR = numeric(), 
    P = numeric(), q_genome = numeric(), logCPM = numeric(), dispersion = numeric(), n_Control = integer(), n_Disease = integer(), 
    model_formula = character(), covariate_status = character(), uncertainty_status = character())

build_design <- function(sm, disease) {
    group <- factor(sm$group, levels = c("Control", disease))
    design <- model.matrix(~group)
    model_formula <- "~ group"
    covariate_status <- "T2D_group_only_prespecified"
    if (disease == "CAD") {
        assert("sex" %in% names(sm), "CAD samples table requires the public sex column.")
        sex_complete <- !anyNA(sm$sex) && all(nzchar(trimws(as.character(sm$sex))))
        if (sex_complete && length(unique(sm$sex)) > 1L) {
            sex <- factor(sm$sex)
            proposed <- model.matrix(~sex + group)
            if (qr(proposed)$rank == ncol(proposed) && nrow(proposed) > ncol(proposed)) {
                design <- proposed
                model_formula <- "~ sex + group"
                covariate_status <- "CAD_public_sex_adjusted_full_rank"
            }
            else covariate_status <- "CAD_group_only_fallback_sex_design_rank_deficient_or_no_residual_df"
        }
        else if (sex_complete) 
            covariate_status <- "CAD_group_only_sex_constant_within_cell_type"
        else covariate_status <- "CAD_group_only_fallback_public_sex_incomplete"
    }
    assert(qr(design)$rank == ncol(design) && nrow(design) > ncol(design), "Non-estimable group model.")
    disease_column <- match(paste0("group", disease), colnames(design))
    assert(!is.na(disease_column), "Disease coefficient absent from design matrix.")
    contrast <- rep(0, ncol(design))
    contrast[disease_column] <- 1
    list(design = design, contrast = contrast, group = group, model_formula = model_formula, covariate_status = covariate_status)
}

analyse_dataset <- function(dataset) {
    ds_out <- file.path(output_root, dataset)
    assert(dir.create(ds_out), "Could not create new dataset output: ", ds_out)
    in_dir <- file.path(input_root, dataset)
    genes <- read.delim(file.path(in_dir, "pseudobulk_genes.tsv"), sep = "\t", quote = "", comment.char = "", check.names = FALSE, 
        na.strings = c("", "NA"))
    samples <- read.delim(file.path(in_dir, "pseudobulk_samples.tsv"), sep = "\t", quote = "", comment.char = "", 
        check.names = FALSE, na.strings = c("", "NA"))
    assert(all(c("gene_id", "gene") %in% names(genes)), dataset, ": genes table requires gene_id,gene.")
    assert(all(c("pb_id", "sample", "group", "cell_type", "n_cells") %in% names(samples)), dataset, ": samples table requires pb_id,sample,group,cell_type,n_cells.")
    genes$gene_id <- as.character(genes$gene_id)
    genes$gene <- as.character(genes$gene)
    for (column in c("pb_id", "sample", "group", "cell_type")) {
        samples[[column]] <- as.character(samples[[column]])
        assert(!anyNA(samples[[column]]) && all(nzchar(samples[[column]])), dataset, ": missing ", column, ".")
    }
    assert(nrow(genes) > 0L && !anyNA(genes$gene_id) && all(nzchar(genes$gene_id)) && !anyDuplicated(genes$gene_id), 
        dataset, ": gene_id must be nonempty and unique.")
    assert(!anyDuplicated(samples$pb_id), dataset, ": pb_id must be unique.")
    assert(!anyDuplicated(samples[c("sample", "cell_type")]), dataset, ": more than one library per nominal donor/cell type; aggregate upstream first.")
    assert(all(samples$group %in% c("Control", "CAD", "T2D", "Treatment")), dataset, ": unexpected group label.")
    group_by_donor <- split(samples$group, samples$sample)
    assert(all(vapply(group_by_donor, function(x) length(unique(x)) == 1L, logical(1))), dataset, ": inconsistent group within sample/donor.")
    samples$n_cells <- suppressWarnings(as.numeric(samples$n_cells))
    assert(all(is.finite(samples$n_cells)) && all(samples$n_cells >= 0) && all(abs(samples$n_cells - round(samples$n_cells)) < 
        1e-08), dataset, ": n_cells must be nonnegative integers.")
    disease <- intersect(c("CAD", "T2D"), unique(samples$group))
    assert(length(disease) == 1L, dataset, ": expected exactly one of CAD/T2D; no pooling or post-hoc contrast selection.")
    if (disease == "CAD") 
        assert("sex" %in% names(samples), dataset, ": CAD samples require public sex metadata.")
    contrast_name <- paste0(disease, "_vs_Control")
    handle <- gzfile(file.path(in_dir, "pseudobulk_counts.mtx.gz"), open = "rt")
    counts <- tryCatch(Matrix::readMM(handle), finally = close(handle))
    counts <- as(counts, "CsparseMatrix")
    assert(identical(dim(counts), as.integer(c(nrow(genes), nrow(samples)))), dataset, ": matrix dimensions do not match ordered genes/samples TSV rows.")
    assert(all(is.finite(counts@x)) && all(counts@x >= 0) && all(abs(counts@x - round(counts@x)) < 1e-08), dataset, 
        ": counts must be finite nonnegative raw integers.")
    rownames(counts) <- genes$gene_id
    colnames(counts) <- samples$pb_id
    samples$total_counts <- as.numeric(Matrix::colSums(counts))
    samples$primary_group <- samples$group %in% c("Control", disease)
    samples$eligible <- samples$primary_group & samples$n_cells >= min_cells & samples$total_counts > 0
    samples$exclusion_reason <- ifelse(!samples$primary_group, "Treatment_excluded", ifelse(samples$n_cells < min_cells, 
        "fewer_than_20_cells", ifelse(samples$total_counts <= 0, "zero_library", "eligible")))
    write_csv(samples, file.path(ds_out, "LIBRARY_ELIGIBILITY.csv"))
    cell_types <- unique(c(primary_types, sort(unique(samples$cell_type))))
    grid <- expand.grid(gene = candidates, cell_type = cell_types, stringsAsFactors = FALSE)
    grid$dataset <- dataset
    grid$contrast <- contrast_name
    grid$gene_id <- NA_character_
    grid$primary_family <- grid$cell_type %in% primary_types
    grid$status <- NA_character_
    grid$status_detail <- ""
    grid$n_Control <- 0L
    grid$n_Disease <- 0L
    grid$model_formula <- "not_fit"
    grid$covariate_status <- "not_fit"
    grid$total_counts_Control <- grid$total_counts_Disease <- NA_real_
    for (field in c("log2FC", "SE", "CI_low", "CI_high", "P", "q_genome", "q_primary8", "q_fixed_all")) grid[[field]] <- NA_real_
    qc_list <- list()
    norm_list <- list()
    for (ct in cell_types) {
        cat("Dataset ", dataset, "; cell type ", ct, "\n", sep = "")
        gi <- which(grid$cell_type == ct)
        ct_input <- which(samples$cell_type == ct & samples$primary_group)
        sel <- which(samples$cell_type == ct & samples$eligible)
        sm <- samples[sel, , drop = FALSE]
        n_control <- length(unique(sm$sample[sm$group == "Control"]))
        n_disease <- length(unique(sm$sample[sm$group == disease]))
        grid$n_Control[gi] <- n_control
        grid$n_Disease[gi] <- n_disease
        qc <- data.frame(dataset = dataset, cell_type = ct, contrast = contrast_name, n_input_Control = sum(samples$group[ct_input] == 
            "Control"), n_input_Disease = sum(samples$group[ct_input] == disease), n_Control = n_control, n_Disease = n_disease, 
            n_tested_genes = 0L, model_formula = "not_fit", covariate_status = "not_fit", status = "", warnings = "", 
            error = "")
        for (i in gi) {
            hit <- which(!is.na(genes$gene) & genes$gene == grid$gene[i])
            if (!length(hit)) 
                grid$status[i] <- "not_covered"
            else if (length(hit) > 1L) {
                grid$status[i] <- "ambiguous_multiple_gene_rows"
                grid$status_detail[i] <- paste(genes$gene_id[hit], collapse = ";")
            }
            else {
                grid$gene_id[i] <- genes$gene_id[hit]
                grid$total_counts_Control[i] <- sum(counts[hit, sel[sm$group == "Control"], drop = FALSE])
                grid$total_counts_Disease[i] <- sum(counts[hit, sel[sm$group == disease], drop = FALSE])
            }
        }
        available <- gi[is.na(grid$status[gi])]
        all_results <- empty_all()
        if (!length(ct_input)) {
            qc$status <- "cell_type_not_present_in_primary_groups"
            grid$status[available] <- qc$status
        }
        else if (n_control < min_donors || n_disease < min_donors) {
            qc$status <- "fewer_than_3_nominal_donors_per_group"
            grid$status[available] <- qc$status
        }
        else {
            warning_messages <- character()
            fitted <- tryCatch(withCallingHandlers({
                cm <- as.matrix(counts[, sel, drop = FALSE])
                model <- build_design(sm, disease)
                design <- model$design
                contrast <- model$contrast
                y <- edgeR::DGEList(counts = cm, genes = genes)
                y <- edgeR::calcNormFactors(y, method = "TMM")
                keep <- edgeR::filterByExpr(y, design = design, min.count = 10, min.total.count = 15, large.n = 10, 
                  min.prop = 0.7)
                if (!any(keep)) 
                  return_value <- list(keep = keep, table = empty_all(), norm = NULL, model = model)
                else {
                  y <- y[keep, , keep.lib.sizes = FALSE]
                  y <- edgeR::calcNormFactors(y, method = "TMM")
                  y <- edgeR::estimateDisp(y, design, robust = TRUE)
                  fit <- edgeR::glmFit(y, design, prior.count = 0)
                  lrt <- edgeR::glmLRT(fit, contrast = contrast)
                  tab <- lrt$table
                  mu <- fit$fitted.values
                  phi <- rep_len(as.numeric(fit$dispersion), nrow(mu))
                  weight <- mu/(1 + mu * phi)
                  se <- vapply(seq_len(nrow(mu)), function(i) {
                    information <- crossprod(design, design * weight[i, ])
                    tryCatch(sqrt(drop(crossprod(contrast, solve(information, contrast))))/log(2), error = function(e) NA_real_)
                  }, numeric(1))
                  out <- data.frame(dataset = dataset, cell_type = ct, contrast = contrast_name, gene_id = y$genes$gene_id, 
                    gene = y$genes$gene, log2FC = tab$logFC, SE = se, CI_low = tab$logFC - qnorm(0.975) * se, CI_high = tab$logFC + 
                      qnorm(0.975) * se, LR = tab$LR, P = tab$PValue, q_genome = bh_fixed(tab$PValue, nrow(tab)), 
                    logCPM = tab$logCPM, dispersion = phi, n_Control = n_control, n_Disease = n_disease, model_formula = model$model_formula, 
                    covariate_status = model$covariate_status, uncertainty_status = ifelse(is.finite(se) & is.finite(tab$logFC), 
                      "conditional_dispersion_normal_approximation", "nonfinite_estimate_or_Fisher_information"))
                  norm <- data.frame(dataset = dataset, cell_type = ct, pb_id = sm$pb_id, sample = sm$sample, group = sm$group, 
                    n_cells = sm$n_cells, sex = if ("sex" %in% names(sm)) 
                      as.character(sm$sex)
                    else NA_character_, lib_size_after_gene_filter = y$samples$lib.size, TMM_factor = y$samples$norm.factors, 
                    effective_lib_size = y$samples$lib.size * y$samples$norm.factors)
                  return_value <- list(keep = keep, table = out, norm = norm, model = model)
                }
                return_value
            }, warning = function(w) {
                warning_messages <<- c(warning_messages, conditionMessage(w))
                invokeRestart("muffleWarning")
            }), error = function(e) e)
            qc$warnings <- paste(unique(warning_messages), collapse = " | ")
            if (inherits(fitted, "error")) {
                qc$status <- "model_error"
                qc$error <- conditionMessage(fitted)
                grid$status[available] <- "model_error"
                grid$status_detail[available] <- qc$error
            }
            else {
                all_results <- fitted$table
                qc$model_formula <- fitted$model$model_formula
                qc$covariate_status <- fitted$model$covariate_status
                grid$model_formula[gi] <- fitted$model$model_formula
                grid$covariate_status[gi] <- fitted$model$covariate_status
                qc$n_tested_genes <- nrow(all_results)
                qc$status <- if (nrow(all_results)) 
                  "fitted"
                else "no_genes_pass_filterByExpr"
                if (!is.null(fitted$norm)) 
                  norm_list[[ct]] <- fitted$norm
                design_table <- data.frame(pb_id = sm$pb_id, sample = sm$sample, fitted$model$design, check.names = FALSE)
                write_csv(design_table, file.path(ds_out, sprintf("DESIGN_%02d.csv", match(ct, cell_types))))
                for (i in available) {
                  hit <- match(grid$gene_id[i], all_results$gene_id)
                  if (is.na(hit)) 
                    grid$status[i] <- "low_counts_filterByExpr"
                  else {
                    for (field in c("log2FC", "SE", "CI_low", "CI_high", "P", "q_genome")) {
                      grid[i, field] <- all_results[hit, field]
                    }
                    grid$status[i] <- if (is.finite(grid$P[i])) 
                      "tested"
                    else "nonfinite_LRT_result"
                    grid$status_detail[i] <- all_results$uncertainty_status[hit]
                  }
                }
            }
        }
        ct_index <- match(ct, cell_types)
        ct_file <- sprintf("ALL_TESTED_%02d.csv", ct_index)
        write_csv(all_results, file.path(ds_out, ct_file))
        qc$all_tested_file <- ct_file
        qc_list[[ct]] <- qc
    }
    primary <- which(grid$primary_family)
    assert(length(primary) == 8L, "Internal error: primary family is not exactly 8 slots.")
    grid$q_primary8[primary] <- bh_fixed(grid$P[primary], 8L)
    grid$q_fixed_all <- bh_fixed(grid$P, nrow(grid))
    grid$q_primary8_family_n <- 8L
    grid$q_fixed_all_family_n <- nrow(grid)
    grid$q_fixed_all_scope <- "exploratory_fixed_candidates_all_cell_types_NA_slots_retained"
    write_csv(grid, file.path(ds_out, "FIXED4_ALL_CELL_TYPES.csv"))
    write_csv(grid[primary, , drop = FALSE], file.path(ds_out, "PRIMARY8_FAMILY.csv"))
    qc <- bind_rows(qc_list)
    write_csv(qc, file.path(ds_out, "CELL_TYPE_MODEL_AUDIT.csv"))
    norms <- bind_rows(norm_list)
    if (nrow(norms)) 
        write_csv(norms, file.path(ds_out, "MODEL_LIBRARY_NORMALIZATION.csv"))
    list(candidates = grid, qc = qc, error = "")
}

dataset_errors <- list()

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        options(stringsAsFactors = FALSE, warn = 1)
        env_or <- function(name, default) {
            value <- Sys.getenv(name, unset = "")
            if (nzchar(value)) value else default
        }
        resolve_path <- function(path, must_exist = FALSE) {
            path <- path.expand(path)
            if (!grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)) path <- file.path(getwd(), path)
            normalizePath(path, winslash = "/", mustWork = must_exist)
        }
        fail <- function(...) stop(paste0(...), call. = FALSE)
        assert <- function(condition, ...) if (!isTRUE(condition)) fail(...)
        path_exists <- function(path) {
            link <- Sys.readlink(path)
            file.exists(path) || dir.exists(path) || (!is.na(link) && nzchar(link))
        }
        write_csv <- function(x, path) {
            assert(!path_exists(path), "Refusing to overwrite: ", path)
            write.csv(x, path, row.names = FALSE, na = "")
        }
        write_text <- function(x, path) {
            assert(!path_exists(path), "Refusing to overwrite: ", path)
            writeLines(x, path, useBytes = TRUE)
        }
        bind_rows <- function(x) {
            x <- Filter(function(z) !is.null(z) && nrow(z) > 0L, x)
            if (!length(x)) return(data.frame())
            do.call(rbind, x)
        }
        bh_fixed <- function(p, family_n) {
            out <- rep(NA_real_, length(p))
            valid <- is.finite(p) & p >= 0 & p <= 1
            if (any(valid)) out[valid] <- p.adjust(p[valid], method = "BH", n = family_n)
            out
        }
        file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
        script_default <- if (length(file_arg)) {
            dirname(resolve_path(sub("^--file=", "", file_arg[[1L]]), TRUE))
        } else file.path(getwd(), "scripts")
        script_root <- resolve_path(env_or("SC_SCRIPT_ROOT", script_default), TRUE)
        task_root <- resolve_path(env_or("SC_TASK_ROOT", dirname(script_root)), TRUE)
        input_root <- resolve_path(env_or("SC_INPUT_ROOT", file.path(task_root, "results")), TRUE)
        output_root <- resolve_path(env_or("SC_OUTPUT_ROOT", file.path(task_root, "results", "pseudobulk_v1")))
        r_lib <- resolve_path(env_or("SC_R_LIB", file.path(dirname(task_root), "birthorder_geo_transcriptome_validation_20261101_v2", 
            "runtime", "Rlib")), TRUE)
        cat("Resolved task root: ", task_root, "\nResolved script root: ", script_root, "\nResolved input root: ", 
            input_root, "\nResolved output root: ", output_root, "\nRead-only supplemental R library: ", r_lib, 
            "\n", sep = "")
        assert(!path_exists(output_root), "Output root already exists; select a new versioned SC_OUTPUT_ROOT: ", 
            output_root)
        assert(output_root != input_root && output_root != script_root, "Output root must differ from input and script roots.")
        .libPaths(unique(c(r_lib, .libPaths())))
        needed <- c("Matrix", "edgeR", "locfit", "statmod")
        for (pkg in needed) assert(requireNamespace(pkg, quietly = TRUE), "Missing R package ", pkg, " in resolved libraries; no package installation is performed.")
        candidates <- c("FCRL6", "GZMA", "GZMH", "KLRD1")
        primary_types <- c("NK", "CD8_T")
        min_cells <- 20L
        min_donors <- 3L
        dataset_args <- commandArgs(trailingOnly = TRUE)
        datasets_env <- Sys.getenv("SC_DATASETS", unset = "")
        datasets <- if (length(dataset_args)) dataset_args else if (nzchar(datasets_env)) {
            trimws(strsplit(datasets_env, ",", fixed = TRUE)[[1L]])
        } else {
            dirs <- list.dirs(input_root, full.names = FALSE, recursive = FALSE)
            dirs[file.exists(file.path(input_root, dirs, "pseudobulk_samples.tsv"))]
        }
        datasets <- unique(datasets[nzchar(datasets)])
        assert(length(datasets) > 0L, "No dataset inputs found under ", input_root, ". Supply dataset names as arguments or comma-separated SC_DATASETS.")
        assert(all(grepl("^[A-Za-z0-9_.-]+$", datasets)) && !any(datasets %in% c(".", "..")), "Dataset names must be simple directory names, without path separators.")
        input_names <- c("pseudobulk_counts.mtx.gz", "pseudobulk_genes.tsv", "pseudobulk_samples.tsv")
        input_files <- do.call(rbind, lapply(datasets, function(dataset) {
            data.frame(dataset = dataset, file = file.path(input_root, dataset, input_names))
        }))
        assert(all(file.exists(input_files$file)), "Missing required input file(s): ", paste(input_files$file[!file.exists(input_files$file)], 
            collapse = "; "))
        for (dataset in datasets) {
            dataset_path <- resolve_path(file.path(input_root, dataset), TRUE)
            assert(output_root != dataset_path && !startsWith(output_root, paste0(dataset_path, "/")), "Outputs must be outside each dataset input directory: ", 
                dataset_path)
        }
        input_files$bytes <- file.info(input_files$file)$size
        input_files$md5 <- unname(tools::md5sum(input_files$file))
        if (!dir.exists(dirname(output_root))) dir.create(dirname(output_root), recursive = TRUE)
        assert(dir.create(output_root, recursive = FALSE), "Could not create new output root: ", output_root)
        write_csv(input_files, file.path(output_root, "INPUT_MANIFEST.csv"))
        write_text(c(paste("UTC start:", format(Sys.time(), tz = "UTC", usetz = TRUE)), paste("Datasets:", paste(datasets, 
            collapse = ", ")), paste("Candidates:", paste(candidates, collapse = ", ")), "Primary family per dataset: 4 fixed genes x NK/CD8_T = 8, including untestable NA slots.", 
            "Exploratory candidate family: 4 fixed genes x all observed cell types plus NK/CD8_T.", "Minimum cells per nominal donor/cell type = 20; minimum eligible nominal donors/group = 3.", 
            "Groups compared independently within each dataset; Treatment is excluded.", "Model: edgeR TMM, filterByExpr, robust estimateDisp, glmFit(prior.count=0), glmLRT.", 
            "CAD: sex + group if the public sex covariate is complete and design is full rank; otherwise documented group-only fallback.", 
            "T2D: group only; no inferred sex or additional age model.", "filterByExpr: min.count=10, min.total.count=15, large.n=10, min.prop=0.7.", 
            "Direction: positive log2FC means disease (CAD or T2D) above Control.", "95% CI: inverse expected NB Fisher information, weights mu/(1+phi*mu), conditional on EB tagwise dispersion.", 
            "CI uses normal approximation; P is likelihood-ratio test, not a Wald test; dispersion uncertainty is not propagated.", 
            "No per-cell disease tests; sample IDs are nominal donors whose biological independence must be audited upstream.", 
            "No birth-order covariate is available here; disease expression is not birth-order mediation validation.", 
            "Cell-level detection percentages and compositions are supplied by upstream preprocessing, without added tests."), 
            file.path(output_root, "METHOD_AND_RUN.txt"))
        empty_all <- function() data.frame(dataset = character(), cell_type = character(), contrast = character(), 
            gene_id = character(), gene = character(), log2FC = numeric(), SE = numeric(), CI_low = numeric(), 
            CI_high = numeric(), LR = numeric(), P = numeric(), q_genome = numeric(), logCPM = numeric(), dispersion = numeric(), 
            n_Control = integer(), n_Disease = integer(), model_formula = character(), covariate_status = character(), 
            uncertainty_status = character())
        build_design <- function(sm, disease) {
            group <- factor(sm$group, levels = c("Control", disease))
            design <- model.matrix(~group)
            model_formula <- "~ group"
            covariate_status <- "T2D_group_only_prespecified"
            if (disease == "CAD") {
                assert("sex" %in% names(sm), "CAD samples table requires the public sex column.")
                sex_complete <- !anyNA(sm$sex) && all(nzchar(trimws(as.character(sm$sex))))
                if (sex_complete && length(unique(sm$sex)) > 1L) {
                  sex <- factor(sm$sex)
                  proposed <- model.matrix(~sex + group)
                  if (qr(proposed)$rank == ncol(proposed) && nrow(proposed) > ncol(proposed)) {
                    design <- proposed
                    model_formula <- "~ sex + group"
                    covariate_status <- "CAD_public_sex_adjusted_full_rank"
                  } else covariate_status <- "CAD_group_only_fallback_sex_design_rank_deficient_or_no_residual_df"
                } else if (sex_complete) covariate_status <- "CAD_group_only_sex_constant_within_cell_type" else covariate_status <- "CAD_group_only_fallback_public_sex_incomplete"
            }
            assert(qr(design)$rank == ncol(design) && nrow(design) > ncol(design), "Non-estimable group model.")
            disease_column <- match(paste0("group", disease), colnames(design))
            assert(!is.na(disease_column), "Disease coefficient absent from design matrix.")
            contrast <- rep(0, ncol(design))
            contrast[disease_column] <- 1
            list(design = design, contrast = contrast, group = group, model_formula = model_formula, covariate_status = covariate_status)
        }
        analyse_dataset <- function(dataset) {
            ds_out <- file.path(output_root, dataset)
            assert(dir.create(ds_out), "Could not create new dataset output: ", ds_out)
            in_dir <- file.path(input_root, dataset)
            genes <- read.delim(file.path(in_dir, "pseudobulk_genes.tsv"), sep = "\t", quote = "", comment.char = "", 
                check.names = FALSE, na.strings = c("", "NA"))
            samples <- read.delim(file.path(in_dir, "pseudobulk_samples.tsv"), sep = "\t", quote = "", comment.char = "", 
                check.names = FALSE, na.strings = c("", "NA"))
            assert(all(c("gene_id", "gene") %in% names(genes)), dataset, ": genes table requires gene_id,gene.")
            assert(all(c("pb_id", "sample", "group", "cell_type", "n_cells") %in% names(samples)), dataset, ": samples table requires pb_id,sample,group,cell_type,n_cells.")
            genes$gene_id <- as.character(genes$gene_id)
            genes$gene <- as.character(genes$gene)
            for (column in c("pb_id", "sample", "group", "cell_type")) {
                samples[[column]] <- as.character(samples[[column]])
                assert(!anyNA(samples[[column]]) && all(nzchar(samples[[column]])), dataset, ": missing ", column, 
                  ".")
            }
            assert(nrow(genes) > 0L && !anyNA(genes$gene_id) && all(nzchar(genes$gene_id)) && !anyDuplicated(genes$gene_id), 
                dataset, ": gene_id must be nonempty and unique.")
            assert(!anyDuplicated(samples$pb_id), dataset, ": pb_id must be unique.")
            assert(!anyDuplicated(samples[c("sample", "cell_type")]), dataset, ": more than one library per nominal donor/cell type; aggregate upstream first.")
            assert(all(samples$group %in% c("Control", "CAD", "T2D", "Treatment")), dataset, ": unexpected group label.")
            group_by_donor <- split(samples$group, samples$sample)
            assert(all(vapply(group_by_donor, function(x) length(unique(x)) == 1L, logical(1))), dataset, ": inconsistent group within sample/donor.")
            samples$n_cells <- suppressWarnings(as.numeric(samples$n_cells))
            assert(all(is.finite(samples$n_cells)) && all(samples$n_cells >= 0) && all(abs(samples$n_cells - round(samples$n_cells)) < 
                1e-08), dataset, ": n_cells must be nonnegative integers.")
            disease <- intersect(c("CAD", "T2D"), unique(samples$group))
            assert(length(disease) == 1L, dataset, ": expected exactly one of CAD/T2D; no pooling or post-hoc contrast selection.")
            if (disease == "CAD") assert("sex" %in% names(samples), dataset, ": CAD samples require public sex metadata.")
            contrast_name <- paste0(disease, "_vs_Control")
            handle <- gzfile(file.path(in_dir, "pseudobulk_counts.mtx.gz"), open = "rt")
            counts <- tryCatch(Matrix::readMM(handle), finally = close(handle))
            counts <- as(counts, "CsparseMatrix")
            assert(identical(dim(counts), as.integer(c(nrow(genes), nrow(samples)))), dataset, ": matrix dimensions do not match ordered genes/samples TSV rows.")
            assert(all(is.finite(counts@x)) && all(counts@x >= 0) && all(abs(counts@x - round(counts@x)) < 1e-08), 
                dataset, ": counts must be finite nonnegative raw integers.")
            rownames(counts) <- genes$gene_id
            colnames(counts) <- samples$pb_id
            samples$total_counts <- as.numeric(Matrix::colSums(counts))
            samples$primary_group <- samples$group %in% c("Control", disease)
            samples$eligible <- samples$primary_group & samples$n_cells >= min_cells & samples$total_counts > 0
            samples$exclusion_reason <- ifelse(!samples$primary_group, "Treatment_excluded", ifelse(samples$n_cells < 
                min_cells, "fewer_than_20_cells", ifelse(samples$total_counts <= 0, "zero_library", "eligible")))
            write_csv(samples, file.path(ds_out, "LIBRARY_ELIGIBILITY.csv"))
            cell_types <- unique(c(primary_types, sort(unique(samples$cell_type))))
            grid <- expand.grid(gene = candidates, cell_type = cell_types, stringsAsFactors = FALSE)
            grid$dataset <- dataset
            grid$contrast <- contrast_name
            grid$gene_id <- NA_character_
            grid$primary_family <- grid$cell_type %in% primary_types
            grid$status <- NA_character_
            grid$status_detail <- ""
            grid$n_Control <- 0L
            grid$n_Disease <- 0L
            grid$model_formula <- "not_fit"
            grid$covariate_status <- "not_fit"
            grid$total_counts_Control <- grid$total_counts_Disease <- NA_real_
            for (field in c("log2FC", "SE", "CI_low", "CI_high", "P", "q_genome", "q_primary8", "q_fixed_all")) grid[[field]] <- NA_real_
            qc_list <- list()
            norm_list <- list()
            for (ct in cell_types) {
                cat("Dataset ", dataset, "; cell type ", ct, "\n", sep = "")
                gi <- which(grid$cell_type == ct)
                ct_input <- which(samples$cell_type == ct & samples$primary_group)
                sel <- which(samples$cell_type == ct & samples$eligible)
                sm <- samples[sel, , drop = FALSE]
                n_control <- length(unique(sm$sample[sm$group == "Control"]))
                n_disease <- length(unique(sm$sample[sm$group == disease]))
                grid$n_Control[gi] <- n_control
                grid$n_Disease[gi] <- n_disease
                qc <- data.frame(dataset = dataset, cell_type = ct, contrast = contrast_name, n_input_Control = sum(samples$group[ct_input] == 
                  "Control"), n_input_Disease = sum(samples$group[ct_input] == disease), n_Control = n_control, 
                  n_Disease = n_disease, n_tested_genes = 0L, model_formula = "not_fit", covariate_status = "not_fit", 
                  status = "", warnings = "", error = "")
                for (i in gi) {
                  hit <- which(!is.na(genes$gene) & genes$gene == grid$gene[i])
                  if (!length(hit)) grid$status[i] <- "not_covered" else if (length(hit) > 1L) {
                    grid$status[i] <- "ambiguous_multiple_gene_rows"
                    grid$status_detail[i] <- paste(genes$gene_id[hit], collapse = ";")
                  } else {
                    grid$gene_id[i] <- genes$gene_id[hit]
                    grid$total_counts_Control[i] <- sum(counts[hit, sel[sm$group == "Control"], drop = FALSE])
                    grid$total_counts_Disease[i] <- sum(counts[hit, sel[sm$group == disease], drop = FALSE])
                  }
                }
                available <- gi[is.na(grid$status[gi])]
                all_results <- empty_all()
                if (!length(ct_input)) {
                  qc$status <- "cell_type_not_present_in_primary_groups"
                  grid$status[available] <- qc$status
                } else if (n_control < min_donors || n_disease < min_donors) {
                  qc$status <- "fewer_than_3_nominal_donors_per_group"
                  grid$status[available] <- qc$status
                } else {
                  warning_messages <- character()
                  fitted <- tryCatch(withCallingHandlers({
                    cm <- as.matrix(counts[, sel, drop = FALSE])
                    model <- build_design(sm, disease)
                    design <- model$design
                    contrast <- model$contrast
                    y <- edgeR::DGEList(counts = cm, genes = genes)
                    y <- edgeR::calcNormFactors(y, method = "TMM")
                    keep <- edgeR::filterByExpr(y, design = design, min.count = 10, min.total.count = 15, large.n = 10, 
                      min.prop = 0.7)
                    if (!any(keep)) return_value <- list(keep = keep, table = empty_all(), norm = NULL, model = model) else {
                      y <- y[keep, , keep.lib.sizes = FALSE]
                      y <- edgeR::calcNormFactors(y, method = "TMM")
                      y <- edgeR::estimateDisp(y, design, robust = TRUE)
                      fit <- edgeR::glmFit(y, design, prior.count = 0)
                      lrt <- edgeR::glmLRT(fit, contrast = contrast)
                      tab <- lrt$table
                      mu <- fit$fitted.values
                      phi <- rep_len(as.numeric(fit$dispersion), nrow(mu))
                      weight <- mu/(1 + mu * phi)
                      se <- vapply(seq_len(nrow(mu)), function(i) {
                        information <- crossprod(design, design * weight[i, ])
                        tryCatch(sqrt(drop(crossprod(contrast, solve(information, contrast))))/log(2), error = function(e) NA_real_)
                      }, numeric(1))
                      out <- data.frame(dataset = dataset, cell_type = ct, contrast = contrast_name, gene_id = y$genes$gene_id, 
                        gene = y$genes$gene, log2FC = tab$logFC, SE = se, CI_low = tab$logFC - qnorm(0.975) * se, 
                        CI_high = tab$logFC + qnorm(0.975) * se, LR = tab$LR, P = tab$PValue, q_genome = bh_fixed(tab$PValue, 
                          nrow(tab)), logCPM = tab$logCPM, dispersion = phi, n_Control = n_control, n_Disease = n_disease, 
                        model_formula = model$model_formula, covariate_status = model$covariate_status, uncertainty_status = ifelse(is.finite(se) & 
                          is.finite(tab$logFC), "conditional_dispersion_normal_approximation", "nonfinite_estimate_or_Fisher_information"))
                      norm <- data.frame(dataset = dataset, cell_type = ct, pb_id = sm$pb_id, sample = sm$sample, 
                        group = sm$group, n_cells = sm$n_cells, sex = if ("sex" %in% names(sm)) as.character(sm$sex) else NA_character_, 
                        lib_size_after_gene_filter = y$samples$lib.size, TMM_factor = y$samples$norm.factors, effective_lib_size = y$samples$lib.size * 
                          y$samples$norm.factors)
                      return_value <- list(keep = keep, table = out, norm = norm, model = model)
                    }
                    return_value
                  }, warning = function(w) {
                    warning_messages <<- c(warning_messages, conditionMessage(w))
                    invokeRestart("muffleWarning")
                  }), error = function(e) e)
                  qc$warnings <- paste(unique(warning_messages), collapse = " | ")
                  if (inherits(fitted, "error")) {
                    qc$status <- "model_error"
                    qc$error <- conditionMessage(fitted)
                    grid$status[available] <- "model_error"
                    grid$status_detail[available] <- qc$error
                  } else {
                    all_results <- fitted$table
                    qc$model_formula <- fitted$model$model_formula
                    qc$covariate_status <- fitted$model$covariate_status
                    grid$model_formula[gi] <- fitted$model$model_formula
                    grid$covariate_status[gi] <- fitted$model$covariate_status
                    qc$n_tested_genes <- nrow(all_results)
                    qc$status <- if (nrow(all_results)) "fitted" else "no_genes_pass_filterByExpr"
                    if (!is.null(fitted$norm)) norm_list[[ct]] <- fitted$norm
                    design_table <- data.frame(pb_id = sm$pb_id, sample = sm$sample, fitted$model$design, check.names = FALSE)
                    write_csv(design_table, file.path(ds_out, sprintf("DESIGN_%02d.csv", match(ct, cell_types))))
                    for (i in available) {
                      hit <- match(grid$gene_id[i], all_results$gene_id)
                      if (is.na(hit)) grid$status[i] <- "low_counts_filterByExpr" else {
                        for (field in c("log2FC", "SE", "CI_low", "CI_high", "P", "q_genome")) {
                          grid[i, field] <- all_results[hit, field]
                        }
                        grid$status[i] <- if (is.finite(grid$P[i])) "tested" else "nonfinite_LRT_result"
                        grid$status_detail[i] <- all_results$uncertainty_status[hit]
                      }
                    }
                  }
                }
                ct_index <- match(ct, cell_types)
                ct_file <- sprintf("ALL_TESTED_%02d.csv", ct_index)
                write_csv(all_results, file.path(ds_out, ct_file))
                qc$all_tested_file <- ct_file
                qc_list[[ct]] <- qc
            }
            primary <- which(grid$primary_family)
            assert(length(primary) == 8L, "Internal error: primary family is not exactly 8 slots.")
            grid$q_primary8[primary] <- bh_fixed(grid$P[primary], 8L)
            grid$q_fixed_all <- bh_fixed(grid$P, nrow(grid))
            grid$q_primary8_family_n <- 8L
            grid$q_fixed_all_family_n <- nrow(grid)
            grid$q_fixed_all_scope <- "exploratory_fixed_candidates_all_cell_types_NA_slots_retained"
            write_csv(grid, file.path(ds_out, "FIXED4_ALL_CELL_TYPES.csv"))
            write_csv(grid[primary, , drop = FALSE], file.path(ds_out, "PRIMARY8_FAMILY.csv"))
            qc <- bind_rows(qc_list)
            write_csv(qc, file.path(ds_out, "CELL_TYPE_MODEL_AUDIT.csv"))
            norms <- bind_rows(norm_list)
            if (nrow(norms)) write_csv(norms, file.path(ds_out, "MODEL_LIBRARY_NORMALIZATION.csv"))
            list(candidates = grid, qc = qc, error = "")
        }
        dataset_results <- list()
        dataset_errors <- list()
        for (dataset in datasets) {
            result <- tryCatch(analyse_dataset(dataset), error = function(e) {
                data.frame(dataset = dataset, status = "dataset_input_or_execution_error", error = conditionMessage(e))
            })
            if (is.data.frame(result)) {
                dataset_errors[[dataset]] <- result
                cat("Dataset failed: ", dataset, ": ", result$error, "\n", sep = "")
            } else dataset_results[[dataset]] <- result
        }
        if (length(dataset_results)) {
            write_csv(bind_rows(lapply(dataset_results, `[[`, "candidates")), file.path(output_root, "FIXED4_ALL_DATASETS.csv"))
            write_csv(bind_rows(lapply(dataset_results, `[[`, "qc")), file.path(output_root, "MODEL_AUDIT_ALL_DATASETS.csv"))
        }
        if (length(dataset_errors)) write_csv(bind_rows(dataset_errors), file.path(output_root, "DATASET_ERRORS.csv"))
        write_text(capture.output(sessionInfo()), file.path(output_root, "SESSION_INFO.txt"))
        write_text(c(paste("UTC finish:", format(Sys.time(), tz = "UTC", usetz = TRUE)), paste("Datasets completed:", 
            length(dataset_results)), paste("Datasets with errors:", length(dataset_errors)), paste("Output root:", 
            output_root)), file.path(output_root, "RUN_FINISHED.txt"))
        if (!length(dataset_results)) fail("No dataset completed. Inspect DATASET_ERRORS.csv in ", output_root)
        cat("Pseudobulk outputs written to ", output_root, "\n", sep = "")
    }), envir = context)
}

