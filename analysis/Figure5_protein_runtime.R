# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

pr_assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)

pr_sha <- function(p) digest::digest(file = p, algo = "sha256")

pr_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)

pr_vec <- function(x) unlist(x, use.names = FALSE)

pr_json <- function(x, p) {
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(p, ".", Sys.getpid(), ".tmp")
    jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, digits = 17, null = "null", na = "null")
    pr_assert(file.rename(tmp, p), "ATOMIC_JSON_COMMIT_FAILED")
}

pr_rds <- function(x, p) {
    pr_assert(!file.exists(p), paste("PRESERVE_EXISTING_RDS", p))
    tmp <- paste0(p, ".", Sys.getpid(), ".tmp")
    saveRDS(x, tmp, version = 2)
    pr_assert(file.rename(tmp, p), "ATOMIC_RDS_COMMIT_FAILED")
    list(path = basename(p), sha256 = pr_sha(p), bytes = file.info(p)$size)
}

pr_read <- function(p) jsonlite::fromJSON(p, simplifyVector = FALSE)

pr_runtime <- function(config, initialize_native = FALSE) {
    pr_assert(identical(Sys.info()[["nodename"]], birth_private_setting("PRIVATE_HOST")) && Sys.info()[["sysname"]] == 
        "Darwin", "APPROVED_MAC_MINI_ONLY")
    for (pkg in c("survival", "splines", "nnet", "digest", "jsonlite")) suppressPackageStartupMessages(library(pkg, 
        character.only = TRUE))
    root <- config$source_root
    old <- file.path(config$y_source_root, "vendor/core39")
    methods <- new.env(parent = globalenv())
    for (p in c("vendor/prior121/R/product_domain_kernel.R", "vendor/prior121/R/y_prediction_adapter.R", "R/product_reverse_isolation.R", 
        "R/product_kernel_runtime.R", "R/selected_nmr_blockspec.R")) sys.source(file.path(root, p), methods)
    sys.source(file.path(root, "R/selected_nmr_composition_A.R"), methods)
    methods$nmr_install_product_runtime(methods)
    num <- new.env(parent = globalenv())
    nd <- file.path(old, "repairs/exponential_fast_b100_v1/numerics")
    sys.source(file.path(nd, "bounded_aj.R"), num)
    if (initialize_native) 
        backend <- num$c39_num_initialize(nd, config$native_build)
    else {
        manifest <- readRDS(file.path(config$native_build, "BUILD_MANIFEST.rds"))
        binary <- file.path(config$native_build, paste0("c39_exponential", .Platform$dynlib.ext))
        pr_assert(identical(pr_sha(binary), manifest$binary_sha256), "NEW_NATIVE_BINARY_HASH_MISMATCH")
        pr_assert(identical(pr_sha(file.path(nd, "piecewise_exponential.c")), manifest$source_sha256), "FROZEN_DOUBLE_SOURCE_CHANGED")
        dll <- dyn.load(binary, local = TRUE)
        backend <- list(symbol = getNativeSymbolInfo("c39_piecewise_exponential", PACKAGE = dll), manifest = manifest, 
            binary = binary)
    }
    engine <- methods$new121_y_engine(old, source_parent = dirname(old), numerical_backend = backend)
    lock <- pr_read(file.path(old, "config/SOURCE_LOCK.json"))
    for (suffix in c("f/00.core.R", "f/20.joint.R")) {
        entry <- Filter(function(x) endsWith(x$path, suffix), lock$R_sources)
        pr_assert(length(entry) == 1L, "PINNED_A_SOURCE_NOT_UNIQUE")
        path <- file.path(dirname(old), entry[[1]]$path)
        pr_assert(identical(pr_sha(path), entry[[1]]$sha256), "PINNED_A_SOURCE_CHANGED")
        sys.source(path, engine)
    }
    engine$c39_install_models(engine)
    sys.source(file.path(old, "repairs/bmi_log_fat_linear_v1/R/mixed_a_revision.R"), engine)
    engine$c39_mixeda_install(engine)
    sys.source(file.path(root, "R/body_tqe_A.R"), engine)
    engine$body_tqe_install(engine)
    sys.source(file.path(root, "R/body_group_pmm_A_v2.R"), engine)
    engine$body_pmm_install(engine)
    sys.source(file.path(config$code_root, "R/ridge_fit.R"), methods)
    sys.source(file.path(config$code_root, "R/ridge_predict.R"), methods)
    sys.source(file.path(config$code_root, "R/a_cache_reuse.R"), methods)
    list(engine = engine, methods = methods, config = config)
}

pr_data <- function(job, rt) {
    config <- rt$config
    ticket <- job$ticket
    id <- ticket$identity
    draw <- ticket$draw
    pr_assert(id$model_id %in% c("S44_PROT_REF11", "S44_PROT_PLUS15") && id$line == "PROTEIN" && id$b >= 0 && id$b <= 
        100 && id$k == 50 && id$uniform_dim == 15, "FROZEN_PROTEIN_SCOPE_CHANGED")
    checked <- function(p, h) {
        pr_assert(file.exists(p) && identical(pr_sha(p), h), paste("FROZEN_INPUT_HASH_MISMATCH", p))
        p
    }
    z <- readRDS(checked(file.path(config$input_root, "LINE_INPUTS.private.rds"), id$input_hash))
    b <- as.integer(id$b)
    ir <- readRDS(checked(file.path(config$input_root, "draws", sprintf("bootstrap_%04d.private.rds", b)), draw$index_sha256))
    ws <- readRDS(checked(file.path(config$input_root, "draws", sprintf("weights_%04d.private.rds", b)), draw$weights_sha256))
    target <- if (is.null(ir$targethash)) 
        ir$target_hash
    else ir$targethash
    ix <- ir$indices
    pr_assert(z$schema == "SELECTED44_LINE_FROZEN_INPUTS_V1" && z$line == "PROTEIN" && identical(pr_hash(as.character(z$data$eid)), 
        z$target_hash) && identical(target, z$target_hash) && identical(z$target_hash, id$population_hash) && as.integer(ir$bootstrap_id) == 
        b && identical(ir$cluster_membership_hash, z$cluster_membership_hash) && is.numeric(ix) && !anyNA(ix) && 
        all(ix == floor(ix) & ix >= 1 & ix <= nrow(z$data)) && identical(pr_hash(as.integer(ix)), id$sampling_hash) && 
        is.logical(ws$complete) && length(ws$complete) == length(ix) && !anyNA(ws$complete), "STRICT_ORIGINAL_SAMPLE_BINDING_CHANGED")
    rows <- ix[ws$complete]
    d <- z$data[rows, , drop = FALSE]
    w <- ws$w
    pb <- ws$pair_binding
    pr_assert(length(w) == nrow(d) && all(is.finite(w) & w > 0) && identical(pb$target_hash, z$target_hash) && 
        identical(pb$sample_hash, pr_hash(as.character(d$eid))) && identical(pb$design_weight_hash, pr_hash(ws$design_weights)) && 
        identical(pb$selection_weight_hash, pr_hash(ws$selection_weights)) && identical(pb$sampling_hash, id$sampling_hash) && 
        identical(pb$actual_index_sha256, draw$index_sha256) && identical(as.numeric(w), as.numeric(ws$design_weights * 
        ws$selection_weights)), "STRICT_ORIGINAL_WEIGHT_FAMILY_BINDING_CHANGED")
    ids <- pr_vec(id$model_spec$mediator_ids)
    members <- z$members[ids]
    pr_assert(length(members) == if (id$model_id == "S44_PROT_REF11") 
        11L
    else 15L, "FROZEN_MEMBERS_CHANGED")
    cols <- vapply(members, `[[`, "", "column")
    base <- pr_vec(z$base)
    domains <- Filter(length, lapply(z$domains, function(x) intersect(pr_vec(x), cols)))
    domains <- domains[pr_vec(id$model_spec$allocation_order)]
    design <- rt$methods$new121_y_design_factory(d, members, base, domains, rt$engine, cache_entries = 2L)
    context <- rt$methods$s44pr_A_context(z, ir, ws, ticket, id$input_hash, draw$index_sha256, draw$weights_sha256)
    list(d = d, w = w, risks = lapply(z$risks, function(r) r[rows, , drop = FALSE]), members = members, columns = cols, 
        A_context = context, base = base, domains = domains, design = design, X = attr(design, "fullX"), identity = id, 
        pair_binding = pb, binding = list(input_sha256 = id$input_hash, index_sha256 = draw$index_sha256, weights_sha256 = draw$weights_sha256, 
            sampling_hash = id$sampling_hash, eid_hash = pr_hash(as.character(d$eid)), family_hash = pr_hash(as.character(d$cluster_id)), 
            weights_hash = pr_hash(as.numeric(w)), risk_hash = pr_hash(lapply(z$risks, function(r) r[rows, , drop = FALSE])), 
            X_hash = pr_hash(attr(design, "fullX")), spline_reference_hash = pr_hash(attr(design, "reference")), 
            A_source_hash = id$source_hash, seed = id$seed, pair_binding = pb))
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        pr_assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
        pr_sha <- function(p) digest::digest(file = p, algo = "sha256")
        pr_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)
        pr_vec <- function(x) unlist(x, use.names = FALSE)
        pr_json <- function(x, p) {
            dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
            tmp <- paste0(p, ".", Sys.getpid(), ".tmp")
            jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, digits = 17, null = "null", na = "null")
            pr_assert(file.rename(tmp, p), "ATOMIC_JSON_COMMIT_FAILED")
        }
        pr_rds <- function(x, p) {
            pr_assert(!file.exists(p), paste("PRESERVE_EXISTING_RDS", p))
            tmp <- paste0(p, ".", Sys.getpid(), ".tmp")
            saveRDS(x, tmp, version = 2)
            pr_assert(file.rename(tmp, p), "ATOMIC_RDS_COMMIT_FAILED")
            list(path = basename(p), sha256 = pr_sha(p), bytes = file.info(p)$size)
        }
        pr_read <- function(p) jsonlite::fromJSON(p, simplifyVector = FALSE)
        pr_runtime <- function(config, initialize_native = FALSE) {
            pr_assert(identical(Sys.info()[["nodename"]], birth_private_setting("PRIVATE_HOST")) && Sys.info()[["sysname"]] == 
                "Darwin", "APPROVED_MAC_MINI_ONLY")
            for (pkg in c("survival", "splines", "nnet", "digest", "jsonlite")) suppressPackageStartupMessages(library(pkg, 
                character.only = TRUE))
            root <- config$source_root
            old <- file.path(config$y_source_root, "vendor/core39")
            methods <- new.env(parent = globalenv())
            for (p in c("vendor/prior121/R/product_domain_kernel.R", "vendor/prior121/R/y_prediction_adapter.R", 
                "R/product_reverse_isolation.R", "R/product_kernel_runtime.R", "R/selected_nmr_blockspec.R")) sys.source(file.path(root, 
                p), methods)
            sys.source(file.path(root, "R/selected_nmr_composition_A.R"), methods)
            methods$nmr_install_product_runtime(methods)
            num <- new.env(parent = globalenv())
            nd <- file.path(old, "repairs/exponential_fast_b100_v1/numerics")
            sys.source(file.path(nd, "bounded_aj.R"), num)
            if (initialize_native) backend <- num$c39_num_initialize(nd, config$native_build) else {
                manifest <- readRDS(file.path(config$native_build, "BUILD_MANIFEST.rds"))
                binary <- file.path(config$native_build, paste0("c39_exponential", .Platform$dynlib.ext))
                pr_assert(identical(pr_sha(binary), manifest$binary_sha256), "NEW_NATIVE_BINARY_HASH_MISMATCH")
                pr_assert(identical(pr_sha(file.path(nd, "piecewise_exponential.c")), manifest$source_sha256), 
                  "FROZEN_DOUBLE_SOURCE_CHANGED")
                dll <- dyn.load(binary, local = TRUE)
                backend <- list(symbol = getNativeSymbolInfo("c39_piecewise_exponential", PACKAGE = dll), manifest = manifest, 
                  binary = binary)
            }
            engine <- methods$new121_y_engine(old, source_parent = dirname(old), numerical_backend = backend)
            lock <- pr_read(file.path(old, "config/SOURCE_LOCK.json"))
            for (suffix in c("f/00.core.R", "f/20.joint.R")) {
                entry <- Filter(function(x) endsWith(x$path, suffix), lock$R_sources)
                pr_assert(length(entry) == 1L, "PINNED_A_SOURCE_NOT_UNIQUE")
                path <- file.path(dirname(old), entry[[1]]$path)
                pr_assert(identical(pr_sha(path), entry[[1]]$sha256), "PINNED_A_SOURCE_CHANGED")
                sys.source(path, engine)
            }
            engine$c39_install_models(engine)
            sys.source(file.path(old, "repairs/bmi_log_fat_linear_v1/R/mixed_a_revision.R"), engine)
            engine$c39_mixeda_install(engine)
            sys.source(file.path(root, "R/body_tqe_A.R"), engine)
            engine$body_tqe_install(engine)
            sys.source(file.path(root, "R/body_group_pmm_A_v2.R"), engine)
            engine$body_pmm_install(engine)
            sys.source(file.path(config$code_root, "R/ridge_fit.R"), methods)
            sys.source(file.path(config$code_root, "R/ridge_predict.R"), methods)
            sys.source(file.path(config$code_root, "R/a_cache_reuse.R"), methods)
            list(engine = engine, methods = methods, config = config)
        }
        pr_data <- function(job, rt) {
            config <- rt$config
            ticket <- job$ticket
            id <- ticket$identity
            draw <- ticket$draw
            pr_assert(id$model_id %in% c("S44_PROT_REF11", "S44_PROT_PLUS15") && id$line == "PROTEIN" && id$b >= 
                0 && id$b <= 100 && id$k == 50 && id$uniform_dim == 15, "FROZEN_PROTEIN_SCOPE_CHANGED")
            checked <- function(p, h) {
                pr_assert(file.exists(p) && identical(pr_sha(p), h), paste("FROZEN_INPUT_HASH_MISMATCH", p))
                p
            }
            z <- readRDS(checked(file.path(config$input_root, "LINE_INPUTS.private.rds"), id$input_hash))
            b <- as.integer(id$b)
            ir <- readRDS(checked(file.path(config$input_root, "draws", sprintf("bootstrap_%04d.private.rds", b)), 
                draw$index_sha256))
            ws <- readRDS(checked(file.path(config$input_root, "draws", sprintf("weights_%04d.private.rds", b)), 
                draw$weights_sha256))
            target <- if (is.null(ir$targethash)) ir$target_hash else ir$targethash
            ix <- ir$indices
            pr_assert(z$schema == "SELECTED44_LINE_FROZEN_INPUTS_V1" && z$line == "PROTEIN" && identical(pr_hash(as.character(z$data$eid)), 
                z$target_hash) && identical(target, z$target_hash) && identical(z$target_hash, id$population_hash) && 
                as.integer(ir$bootstrap_id) == b && identical(ir$cluster_membership_hash, z$cluster_membership_hash) && 
                is.numeric(ix) && !anyNA(ix) && all(ix == floor(ix) & ix >= 1 & ix <= nrow(z$data)) && identical(pr_hash(as.integer(ix)), 
                id$sampling_hash) && is.logical(ws$complete) && length(ws$complete) == length(ix) && !anyNA(ws$complete), 
                "STRICT_ORIGINAL_SAMPLE_BINDING_CHANGED")
            rows <- ix[ws$complete]
            d <- z$data[rows, , drop = FALSE]
            w <- ws$w
            pb <- ws$pair_binding
            pr_assert(length(w) == nrow(d) && all(is.finite(w) & w > 0) && identical(pb$target_hash, z$target_hash) && 
                identical(pb$sample_hash, pr_hash(as.character(d$eid))) && identical(pb$design_weight_hash, pr_hash(ws$design_weights)) && 
                identical(pb$selection_weight_hash, pr_hash(ws$selection_weights)) && identical(pb$sampling_hash, 
                id$sampling_hash) && identical(pb$actual_index_sha256, draw$index_sha256) && identical(as.numeric(w), 
                as.numeric(ws$design_weights * ws$selection_weights)), "STRICT_ORIGINAL_WEIGHT_FAMILY_BINDING_CHANGED")
            ids <- pr_vec(id$model_spec$mediator_ids)
            members <- z$members[ids]
            pr_assert(length(members) == if (id$model_id == "S44_PROT_REF11") 11L else 15L, "FROZEN_MEMBERS_CHANGED")
            cols <- vapply(members, `[[`, "", "column")
            base <- pr_vec(z$base)
            domains <- Filter(length, lapply(z$domains, function(x) intersect(pr_vec(x), cols)))
            domains <- domains[pr_vec(id$model_spec$allocation_order)]
            design <- rt$methods$new121_y_design_factory(d, members, base, domains, rt$engine, cache_entries = 2L)
            context <- rt$methods$s44pr_A_context(z, ir, ws, ticket, id$input_hash, draw$index_sha256, draw$weights_sha256)
            list(d = d, w = w, risks = lapply(z$risks, function(r) r[rows, , drop = FALSE]), members = members, 
                columns = cols, A_context = context, base = base, domains = domains, design = design, X = attr(design, 
                  "fullX"), identity = id, pair_binding = pb, binding = list(input_sha256 = id$input_hash, index_sha256 = draw$index_sha256, 
                  weights_sha256 = draw$weights_sha256, sampling_hash = id$sampling_hash, eid_hash = pr_hash(as.character(d$eid)), 
                  family_hash = pr_hash(as.character(d$cluster_id)), weights_hash = pr_hash(as.numeric(w)), risk_hash = pr_hash(lapply(z$risks, 
                    function(r) r[rows, , drop = FALSE])), X_hash = pr_hash(attr(design, "fullX")), spline_reference_hash = pr_hash(attr(design, 
                    "reference")), A_source_hash = id$source_hash, seed = id$seed, pair_binding = pb))
        }
    }), envir = context)
}

