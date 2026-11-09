# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

sha <- function(p) digest(file = p, algo = "sha256")

hash <- function(x) digest(x, algo = "sha256", serializeVersion = 2L)

readj <- function(p) fromJSON(p, simplifyVector = FALSE)

vec <- function(x) as.character(unlist(x, use.names = FALSE))

writej <- function(x, p) {
    stopifnot(!file.exists(p))
    tmp <- paste0(p, ".pending_", Sys.getpid())
    stopifnot(!file.exists(tmp))
    write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
    stopifnot(!file.exists(p), file.rename(tmp, p))
}

writer <- function(x, p) {
    stopifnot(!file.exists(p))
    tmp <- paste0(p, ".pending_", Sys.getpid())
    stopifnot(!file.exists(tmp))
    saveRDS(x, tmp, version = 2L, compress = FALSE)
    stopifnot(!file.exists(p), file.rename(tmp, p))
}

selection_for <- function(d, columns, base) {
    z <- se$b3_missing(d, columns, base)
    stopifnot(identical(z$complete, complete.cases(d[, unique(c(base, columns)), drop = FALSE])), length(z$w) == 
        sum(z$complete), all(is.finite(z$w) & z$w > 0))
    z$design_weights <- rep(1, length(z$w))
    z$selection_weights <- as.numeric(z$w)
    z$selection_coefficients <- if (is.null(z$fit)) 
        NULL
    else stats::coef(z$fit)
    z$selection_formula <- if (is.null(z$fit)) 
        NULL
    else paste(deparse(stats::formula(z$fit)), collapse = " ")
    z$fit <- NULL
    z
}

preserve_rng <- function(code) {
    kind <- RNGkind()
    had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    if (had) 
        before <- get(".Random.seed", .GlobalEnv)
    on.exit({
        do.call(RNGkind, as.list(kind))
        if (had) assign(".Random.seed", before, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", 
            envir = .GlobalEnv)
    })
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    force(code)
}

endpoint_counts <- function(r, keep) list(N_target = nrow(r), N_fit = sum(keep), target_events_full = sum(r$event), 
    target_events_11y = sum(r$event == 1 & r$exit <= 11), target_competing_deaths_full = sum(r$death_event), target_competing_deaths_11y = sum(r$death_event == 
        1 & r$exit <= 11), fit_events_full = sum(r$event[keep]), fit_events_11y = sum(r$event[keep] == 1 & r$exit[keep] <= 
        11), fit_competing_deaths_full = sum(r$death_event[keep]), fit_competing_deaths_11y = sum(r$death_event[keep] == 
        1 & r$exit[keep] <= 11), followup_max = max(r$exit), source_followup_not_capped_at_11 = TRUE)

prepared_lines <- list()

line_context <- list()

line_models <- list(BASE = "S44_BASE11", NMR = c("S44_NMR_REF11", "S44_NMR_PLUS40"), PROTEIN = c("S44_PROT_REF11", 
    "S44_PROT_PLUS15"))

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        args <- commandArgs(TRUE)
        stopifnot(length(args) == 2L, Sys.info()[["sysname"]] != "Darwin", Sys.getenv("FF_HOST") == "winpc", Sys.getenv("FD_REAL_EXECUTION_AUTHORIZED") == 
            "YES_FOR_THIS_RUN")
        root <- normalizePath(args[1], mustWork = TRUE)
        out <- args[2]
        library(jsonlite)
        library(digest)
        sha <- function(p) digest(file = p, algo = "sha256")
        hash <- function(x) digest(x, algo = "sha256", serializeVersion = 2L)
        readj <- function(p) fromJSON(p, simplifyVector = FALSE)
        vec <- function(x) as.character(unlist(x, use.names = FALSE))
        writej <- function(x, p) {
            stopifnot(!file.exists(p))
            tmp <- paste0(p, ".pending_", Sys.getpid())
            stopifnot(!file.exists(tmp))
            write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
            stopifnot(!file.exists(p), file.rename(tmp, p))
        }
        writer <- function(x, p) {
            stopifnot(!file.exists(p))
            tmp <- paste0(p, ".pending_", Sys.getpid())
            stopifnot(!file.exists(tmp))
            saveRDS(x, tmp, version = 2L, compress = FALSE)
            stopifnot(!file.exists(p), file.rename(tmp, p))
        }
        stopifnot(!grepl("(^|/)[.][.](/|$)", out), grepl("^/mnt/g/UKB_results/birth/border_selected44_stratified_20261103_v1/[^/]+$", 
            out), !file.exists(out), !dir.exists(out))
        input_root <- Sys.getenv("SELECTED44_INPUT_ROOT", birth_private_setting("PRIVATE_BIRTH"))
        core_run <- Sys.getenv("SELECTED44_CORE_INPUT_RUN", Sys.getenv("OMICS_CORE39_RUN", file.path(input_root, 
            "border_jointdomain_decomposition_20261101_v2", "run_20261101_CORE39_V21_B100")))
        v5_run <- Sys.getenv("SELECTED44_V5_INPUT_RUN", Sys.getenv("OMICS_V5_RUN", file.path(input_root, "border_165_typewise_b100_v5", 
            "run_v5_b0572ca19997")))
        cat(toJSON(list(input_root = input_root, script_root = root, output_root = out, shared_incident_input_run = core_run, 
            frozen_feature_input_run = v5_run), auto_unbox = TRUE), "\n")
        cp <- file.path(core_run, "private/CORE39_INPUTS.private.rds")
        vp <- file.path(v5_run, "private/V5_INPUTS.private.rds")
        clp <- file.path(core_run, "private/INPUT_MANIFEST.json")
        vlp <- file.path(v5_run, "private/V5_INPUT_LOCK.json")
        required <- c(cp, vp, clp, vlp, file.path(root, "config/MODELS.json"))
        if (!all(file.exists(required))) stop("MISSING_INPUT_PATH: ", paste(required[!file.exists(required)], collapse = "; "))
        cl <- readj(clp)
        vl <- readj(vlp)
        stopifnot(identical(sha(cp), cl$input_hash), identical(sha(vp), vl$prepared_inputs_sha256))
        core <- readRDS(cp)
        v5 <- readRDS(vp)
        registry <- readj(file.path(root, "config/MODELS.json"))
        base_ids <- vec(registry$base11_members)
        stopifnot(identical(registry$schema, "BORDER_SELECTED44_STRATIFIED_REGISTRY_V1_20261103"), length(base_ids) == 
            11L, length(vec(registry$nmr29_members)) == 29L, length(vec(registry$protein4_members)) == 4L, length(unique(vec(registry$selected44_members))) == 
            44L, setequal(names(registry$models), c("S44_BASE11", "S44_NMR_REF11", "S44_NMR_PLUS40", "S44_PROT_REF11", 
            "S44_PROT_PLUS15")), identical(core$schema, "CORE39_FROZEN_ACTUAL_INPUTS_V21"), !anyDuplicated(core$data$eid), 
            !anyDuplicated(v5$d_full$eid), nrow(core$data) == core$counts$N_target, identical(hash(core$selection_functions), 
                core$selection_source_hash), identical(vec(registry$base), as.character(core$model$base)))
        j <- match(as.character(core$data$eid), as.character(v5$d_full$eid))
        stopifnot(!anyNA(j), length(core$parent_rows) == nrow(core$data), identical(as.integer(j), as.integer(core$parent_rows)))
        se <- new.env(parent = globalenv())
        for (n in names(core$selection_functions)) se[[n]] <- eval(parse(text = core$selection_functions[[n]]), 
            envir = se)
        selection_aux_cols <- grep("^aux_", names(core$data), value = TRUE)
        selection_aux_required <- unique(unlist(regmatches(core$selection_functions$b3_missing, gregexpr("aux_[A-Za-z0-9_]+", 
            core$selection_functions$b3_missing))))
        if (!length(selection_aux_required) || !all(selection_aux_required %in% selection_aux_cols)) stop("FROZEN_SELECTION_AUXILIARY_MISSING: ", 
            paste(setdiff(selection_aux_required, selection_aux_cols), collapse = "; "))
        selection_for <- function(d, columns, base) {
            z <- se$b3_missing(d, columns, base)
            stopifnot(identical(z$complete, complete.cases(d[, unique(c(base, columns)), drop = FALSE])), length(z$w) == 
                sum(z$complete), all(is.finite(z$w) & z$w > 0))
            z$design_weights <- rep(1, length(z$w))
            z$selection_weights <- as.numeric(z$w)
            z$selection_coefficients <- if (is.null(z$fit)) NULL else stats::coef(z$fit)
            z$selection_formula <- if (is.null(z$fit)) NULL else paste(deparse(stats::formula(z$fit)), collapse = " ")
            z$fit <- NULL
            z
        }
        preserve_rng <- function(code) {
            kind <- RNGkind()
            had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
            if (had) before <- get(".Random.seed", .GlobalEnv)
            on.exit({
                do.call(RNGkind, as.list(kind))
                if (had) assign(".Random.seed", before, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, 
                  inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
            })
            RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
            force(code)
        }
        endpoint_counts <- function(r, keep) list(N_target = nrow(r), N_fit = sum(keep), target_events_full = sum(r$event), 
            target_events_11y = sum(r$event == 1 & r$exit <= 11), target_competing_deaths_full = sum(r$death_event), 
            target_competing_deaths_11y = sum(r$death_event == 1 & r$exit <= 11), fit_events_full = sum(r$event[keep]), 
            fit_events_11y = sum(r$event[keep] == 1 & r$exit[keep] <= 11), fit_competing_deaths_full = sum(r$death_event[keep]), 
            fit_competing_deaths_11y = sum(r$death_event[keep] == 1 & r$exit[keep] <= 11), followup_max = max(r$exit), 
            source_followup_not_capped_at_11 = TRUE)
        dir.create(out, recursive = TRUE)
        counts <- list()
        prepared_lines <- list()
        line_context <- list()
        line_models <- list(BASE = "S44_BASE11", NMR = c("S44_NMR_REF11", "S44_NMR_PLUS40"), PROTEIN = c("S44_PROT_REF11", 
            "S44_PROT_PLUS15"))
        for (line in c("BASE", "NMR", "PROTEIN")) {
            model_ids <- line_models[[line]]
            ms <- registry$models[[tail(model_ids, 1L)]]
            all_ids <- vec(ms$members)
            addon_ids <- setdiff(all_ids, base_ids)
            expected <- switch(line, BASE = 11L, NMR = 40L, PROTEIN = 15L)
            stopifnot(length(all_ids) == expected, identical(all_ids, vec(registry$line_members[[line]])), all(vapply(registry$models[model_ids], 
                function(m) !isTRUE(m$total_only) && isTRUE(m$compute_domain_decomposition), FALSE)))
            base <- vec(registry$base)
            d <- core$data[, unique(c("eid", "cluster_id", base, selection_aux_cols)), drop = FALSE]
            members <- setNames(vector("list", length(all_ids)), all_ids)
            member_source <- list()
            for (id in all_ids) {
                if (id %in% names(core$members)) {
                  m <- core$members[[id]]
                  stopifnot(m$column %in% names(core$data))
                  d[[m$column]] <- core$data[[m$column]]
                  member_source[[id]] <- "HASH_LOCKED_SHARED_INCIDENT_INPUT_SELECTED_COLUMN_ONLY"
                } else {
                  f <- v5$features[[id]]
                  vr <- v5$variables[[id]]
                  stopifnot(is.list(f), identical(vr$status, "AVAILABLE_FROZEN"), f$column %in% names(v5$d_full))
                  stopifnot(id == "weight" || id %in% c(vec(registry$nmr29_members), vec(registry$protein4_members)))
                  x <- v5$d_full[[f$column]]
                  type <- "continuous"
                  if (type == "continuous") stopifnot(is.numeric(x), identical(f$prepared_transform, "identity"))
                  d[[f$column]] <- x[j]
                  m <- list(indicator_id = id, column = f$column, type = type, levels = if (type != "continuous") levels(x) else NULL, 
                    lower = if (is.null(f$physical_lower_bound)) -Inf else f$physical_lower_bound, upper = if (is.null(f$physical_upper_bound)) Inf else f$physical_upper_bound, 
                    prepared_transform = f$prepared_transform, prepared_unit = f$prepared_unit, raw_unit = f$raw_unit, 
                    inverse_transform = f$inverse_transform, technical_covariates = f$technical_covariates)
                  member_source[[id]] <- "HASH_LOCKED_V5_SELECTED_COLUMN_EID_ALIGNED"
                }
                members[[id]] <- m
                for (tc in m$technical_covariates) if (!tc %in% names(d)) {
                  stopifnot(tc %in% names(v5$d_full))
                  d[[tc]] <- v5$d_full[[tc]][j]
                }
            }
            base <- unique(c(base, unlist(lapply(members, `[[`, "technical_covariates"))))
            columns <- vapply(members, `[[`, "", "column")
            stopifnot(!anyDuplicated(columns), all(c(base, columns) %in% names(d)), !any(selection_aux_cols %in% 
                c(base, columns)))
            addon_cols <- vapply(members[addon_ids], `[[`, "", "column")
            measured <- if (line == "BASE") rep(TRUE, nrow(d)) else Reduce(`|`, lapply(d[addon_cols], function(x) is.finite(x)))
            rows <- which(measured)
            stopifnot(length(rows) > 0L)
            d <- d[rows, , drop = FALSE]
            rownames(d) <- NULL
            risks <- core$risks[c("411", "250.2")]
            names(risks) <- c("IHD", "T2D")
            risks <- lapply(risks, function(r) r[rows, , drop = FALSE])
            stopifnot(all(vapply(risks, function(r) identical(as.character(r$eid), as.character(d$eid)) && all(is.finite(r$exit) & 
                r$exit > 0 & r$entry == 0 & r$event %in% 0:1 & r$death_event %in% 0:1 & r$event + r$death_event <= 
                1), FALSE)), !anyNA(d$cluster_id), !anyDuplicated(d$eid), all(as.numeric(as.character(d$S)) %in% 
                2:8), all(is.finite(d$B) & d$B >= 0 & d$B <= 1))
            selection <- selection_for(d, columns, base)
            domains <- lapply(ms$domains, function(ids) vapply(members[vec(ids)], `[[`, "", "column"))
            stopifnot(length(domains) == if (line == "BASE") 4L else 5L, setequal(unlist(domains, use.names = FALSE), 
                columns))
            folder <- file.path(out, line)
            dir.create(folder)
            dir.create(file.path(folder, "draws"))
            population_rule <- registry$population_policy[[line]]
            object <- list(schema = "SELECTED44_LINE_FROZEN_INPUTS_V1", line = line, data = d, risks = risks, members = members, 
                core_ids = base_ids, base11_ids = base_ids, base = base, domains = domains, continuous_rule = registry$continuous_rule, 
                count_rule = registry$count_rule, selection_functions = core$selection_functions, selection_source_hash = core$selection_source_hash, 
                target_hash = hash(as.character(d$eid)), cluster_membership_hash = hash(paste(d$eid, d$cluster_id, 
                  sep = ":")), source_core_sha256 = sha(cp), source_v5_sha256 = sha(vp), source_core_input_lock_sha256 = sha(clp), 
                source_v5_input_lock_sha256 = sha(vlp), registry_sha256 = sha(file.path(root, "config/MODELS.json")), 
                population_rule = population_rule, unused39_or121_complete_mask_applied = FALSE, original_eligibility = "Frozen B/C and followup eligible nonadopted S2..8; BOTH baseline IHD and T2D disease-free; original per-disease censoring unchanged", 
                member_source = member_source, actual_A_generation_columns = columns, selection_auxiliary_columns = selection_aux_cols, 
                selection_auxiliary_required_by_frozen_IPW = selection_aux_required, selection_auxiliary_in_A_Y_or_completeness = FALSE, 
                uniform_dim = as.integer(registry$line_uniform_dim[[line]]), model_specifications = registry$models[model_ids])
            ip <- file.path(folder, "LINE_INPUTS.private.rds")
            writer(object, ip)
            base_missing <- !complete.cases(d[, base, drop = FALSE])
            nonomics_missing <- !complete.cases(d[, columns[base_ids], drop = FALSE])
            omics_missing <- if (line == "BASE") rep(FALSE, nrow(d)) else !complete.cases(d[, columns[addon_ids], 
                drop = FALSE])
            ct <- list(line = line, models = model_ids, registered_line_members = expected, source_shared_incident_target_N = nrow(core$data), 
                N_target = nrow(d), N_fit = sum(selection$complete), N_standardization = sum(selection$complete), 
                N_target_families = length(unique(d$cluster_id)), N_fit_families = length(unique(d$cluster_id[selection$complete])), 
                assay_positive_observed_N = length(rows), assay_no_selected_value_N = sum(!measured), N_missing_base = sum(base_missing), 
                N_missing_selected11 = sum(nonomics_missing), N_missing_selected_omics = sum(omics_missing), ordinary_missing_loss_N = sum(!selection$complete), 
                missing_fraction = mean(!selection$complete), missing_by_selected_member = setNames(lapply(columns, 
                  function(cn) sum(is.na(d[[cn]]))), names(members)), completeness_columns = unique(c(base, columns)), 
                selection_auxiliary_columns = selection_aux_cols, selection_auxiliary_in_A_Y_or_completeness = FALSE, 
                unused39_or121_complete_mask_applied = FALSE, target_hash = object$target_hash, sample_hash = hash(as.character(d$eid[selection$complete])), 
                weights = list(min = min(selection$w), median = median(selection$w), max = max(selection$w), sum = sum(selection$w), 
                  ess = sum(selection$w)^2/sum(selection$w^2), selection_formula = selection$selection_formula, 
                  weight_hash = hash(as.numeric(selection$w)), design_weight_hash = hash(selection$design_weights), 
                  selection_weight_hash = hash(selection$selection_weights), coefficient_hash = hash(selection$selection_coefficients), 
                  selection_source_hash = core$selection_source_hash), endpoints = lapply(risks, endpoint_counts, 
                  keep = selection$complete), pair_shared_rows_weights_and_draws = TRUE, source_39_fit_mask_ignored = TRUE, 
                original_eligibility = object$original_eligibility, scientific_scope = population_rule)
            writej(ct, file.path(folder, "POPULATION_MANIFEST.json"))
            counts[[line]] <- ct
            message("SELECTED44_COUNTS_READY ", toJSON(ct, auto_unbox = TRUE, digits = 17))
            flush.console()
            line_context[[line]] <- list(object = object, selection = selection, columns = columns, input_path = ip, 
                folder = folder, model_ids = model_ids)
        }
        writej(counts, file.path(out, "POPULATION_MANIFEST.json"))
        message("SELECTED44_ALL_THREE_POPULATIONS_READY ", file.path(out, "POPULATION_MANIFEST.json"))
        flush.console()
        for (line in c("BASE", "NMR", "PROTEIN")) {
            lc <- line_context[[line]]
            object <- lc$object
            d <- object$data
            selection <- lc$selection
            columns <- lc$columns
            base <- object$base
            ip <- lc$input_path
            folder <- lc$folder
            model_ids <- lc$model_ids
            ct <- counts[[line]]
            clusters <- unique(as.character(d$cluster_id))
            ci <- match(as.character(d$cluster_id), clusters)
            sizes <- tabulate(ci, length(clusters))
            descriptors <- preserve_rng(lapply(0:100, function(b) {
                ix <- se$b3_boot_index(d, as.integer(b), 20261023L)
                stopifnot(all(is.finite(ix) & ix == floor(ix) & ix >= 1L & ix <= nrow(d)))
                if (b == 0L) stopifnot(identical(as.integer(ix), seq_len(nrow(d))))
                mult <- tabulate(ci[ix], length(clusters))/sizes
                stopifnot(all(is.finite(mult) & mult >= 0 & mult == floor(mult)), sum(mult) == length(clusters))
                actual <- list(bootstrap_id = as.integer(b), indices = as.integer(ix), cluster_multiplicity = as.integer(mult), 
                  target_hash = object$target_hash, cluster_membership_hash = object$cluster_membership_hash, bootstrap_seed = as.integer(20261023L + 
                    1009L * b), RNG_kind = RNGkind())
                wp <- file.path(folder, "draws", sprintf("weights_%04d.private.rds", b))
                dp <- file.path(folder, "draws", sprintf("bootstrap_%04d.private.rds", b))
                writer(actual, dp)
                mi <- if (b == 0L) selection else selection_for(d[ix, , drop = FALSE], columns, base)
                mi$pair_binding <- list(target_hash = object$target_hash, sample_hash = hash(as.character(d$eid[ix[mi$complete]])), 
                  design_weight_hash = hash(mi$design_weights), selection_weight_hash = hash(mi$selection_weights), 
                  sampling_hash = hash(actual$indices), actual_index_sha256 = sha(dp))
                writer(mi, wp)
                message("SELECTED44_DRAW_WEIGHTS_FROZEN ", line, " B", b, " Nfit=", sum(mi$complete))
                flush.console()
                list(b = as.integer(b), seed = as.integer(20261103L + 1009L * b + switch(line, BASE = 0L, NMR = 200000L, 
                  PROTEIN = 100000L)), bootstrap_seed = actual$bootstrap_seed, sampling_hash = hash(actual$indices), 
                  index_path = dp, index_sha256 = sha(dp), weights_path = wp, weights_sha256 = sha(wp), N_resampled = length(ix), 
                  N_fit = sum(mi$complete), sample_hash = mi$pair_binding$sample_hash, selection_weight_hash = mi$pair_binding$selection_weight_hash)
            }))
            writej(list(schema = "SELECTED44_ACTUAL_LINE_DRAWS_V1", line = line, fixed_attempts = 100L, unit = "Frozen cluster_id family", 
                reuse_old_indices = FALSE, fixed_ids = as.list(0:100), target_hash = object$target_hash, cluster_membership_hash = object$cluster_membership_hash, 
                uniform_dim = object$uniform_dim, draws = descriptors), file.path(folder, "DRAW_MANIFEST.json"))
            prep <- list(schema = "SELECTED44_PREPARATION_V1", line = line, models = model_ids, input_path = ip, 
                input_hash = sha(ip), population_hash = object$target_hash, draw_manifest = file.path(folder, "DRAW_MANIFEST.json"), 
                draw_manifest_hash = sha(file.path(folder, "DRAW_MANIFEST.json")), inputs_prepared = TRUE, A_Y_models_fitted = 0L, 
                production_started = FALSE, scientific_representation_ready = TRUE, representation_ready_is_not_B0_fitted_support_or_rank_validation = TRUE, 
                uniform_dim = object$uniform_dim, N = ct)
            writej(prep, file.path(folder, "PREPARATION.json"))
            prepared_lines[[line]] <- prep
        }
        writej(list(schema = "SELECTED44_ALL_LINES_PREPARATION_V1", inputs_prepared = TRUE, original44_unique_members = 44L, 
            models = names(registry$models), lines = prepared_lines, original_membership_sha256 = registry$sources$original_panel44$sha256, 
            A_Y_models_fitted = 0L, production_started = FALSE), file.path(out, "PREPARATION.json"))
        cat(toJSON(list(prepared = TRUE, output = out, counts = counts, A_Y_models_fitted = 0L), auto_unbox = TRUE, 
            digits = 17), "\n")
    }), envir = context)
}

