# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

prepare <- function() {
    if (file.exists(file.path(out, "cache/core.rds"))) 
        return(invisible(cp_load(file.path(out, "cache/core.rds"))))
    src <- paths()
    assert(all(file.exists(src)), "REQUIRED_INPUT_MISSING")
    manifest_path <- file.path(out, "INPUT_SHA256.csv")
    if (!file.exists(manifest_path)) {
        ss <- file.info(src)
        h <- vapply(seq_along(src), function(i) {
            log_event("HASH_SOURCE", role = names(src)[i])
            digest(file = src[i], algo = "sha256")
        }, "")
        after <- file.info(src)
        assert(identical(ss$size, after$size) && identical(ss$mtime, after$mtime), "SOURCE_CHANGED_DURING_HASH")
        names(h) <- names(src)
        for (k in names(cfg$expected_sha256)) assert(identical(h[[k]], cfg$expected_sha256[[k]]), paste("FROZEN_HASH_MISMATCH", 
            k))
        publish(data.table(role = names(src), path = unname(src), sha256 = unname(h), bytes = ss$size, mtime = as.numeric(ss$mtime)), 
            manifest_path, "csv")
    }
    else {
        a <- fread(manifest_path)
        ss <- file.info(a$path)
        assert(all(ss$size == a$bytes & as.numeric(ss$mtime) == a$mtime), "SOURCE_CHANGED_SINCE_PREPARE")
    }
    inc <- cp_load(src[["incident"]], strict = FALSE)
    zz <- readRDS(src[["incident"]])
    m <- read_json(file.path(results, cfg$incident_root, "check_complete.json"))
    assert(identical(zz$signature, m$signature), "INCIDENT_SIGNATURE_MISMATCH")
    all <- as.data.frame(readRDS(src[["all"]]))
    d <- linear_covars(make_birth_order(all))
    d <- d[rowSums(!is.na(d[, unlist(cfg$exposures)])) > 0, , drop = FALSE]
    for (k in c("birth_year", "cnt_older", "cnt_brother", "cnt_sister")) assert(k %in% names(d), paste("MISSING_FIELD", 
        k))
    d$eid <- as.character(d$eid)
    ix <- match(d$eid, as.character(inc$data$eid))
    assert(!anyNA(ix), "INCIDENT_CACHE_DOES_NOT_COVER_FULL_EXPOSURE_COHORT")
    for (k in c("S", unlist(cfg$exposures), "age", "sex", "ethnicity", "center", "assessment_year")) assert(isTRUE(all.equal(as.character(d[[k]]), 
        as.character(inc$data[[k]][ix]))), paste("FROZEN_ADAPTER_MISMATCH", k))
    inc$events <- lapply(inc$events, function(v) v[ix])
    inc$baseline <- lapply(inc$baseline, function(v) v[ix])
    inc$data <- inc$data[ix, , drop = FALSE]
    d$date_death_uncapped <- inc$data$date_death_uncapped
    d$date_lost_registry <- inc$data$date_lost_registry
    d <- d[, unique(c("eid", "birth_year", "date_attend", "date_death_uncapped", "date_lost_registry", "age", "sex", 
        "ethnicity", "center", "assessment_year", "S", unlist(cfg$exposures))), drop = FALSE]
    p <- as.data.frame(readRDS(src[["phewas"]]))
    assert("id" %in% names(p) && !anyDuplicated(p$id), "PHEWAS_INVALID_ID")
    j <- match(d$eid, as.character(p$id))
    history <- lapply(unlist(cfg$phecodes), function(k) {
        assert(k %in% names(p), paste("PHECODE_MISSING", k))
        v <- p[[k]][j]
        assert(all(is.na(v) | v %in% c(0, 1)), "UNKNOWN_PHEWAS_STATUS")
        num(v)
    })
    names(history) <- unlist(cfg$phecodes)
    rm(p, all)
    gc(FALSE)
    kin <- bs_kin(src[["kinship"]])
    clusters <- bs_components(kin[kinship >= 0.0442])
    d$.cluster <- paste0("U", d$eid)
    j <- match(d$eid, clusters$eid)
    d$.cluster[!is.na(j)] <- paste0("R", clusters$component[j[!is.na(j)]])
    dd <- dict()
    qmap <- read_json(file.path(root, "config/nmr_qc_map.json"), simplifyVector = TRUE)
    spec <- list(metadata = unique(c("p53_i0", "p53_i2", "p21003_i2", "p54_i2", "p50_i2", "p21002_i2", "p74_i0", 
        "p20282_i0", "p23650_i0", "p30901_i0", "p26227", "p26285", "p22000", paste0("p22009_a", 1:40), "p4080_i0_a0", 
        "p4080_i0_a1", "p4079_i0_a0", "p4079_i0_a1")), clinical = unique(c(dd[module == "clinical", column], dd[module == 
        "clinical", qc_column], paste0("p", num(sub("_i0", "", sub("p", "", dd[module == "clinical", column]))) + 
        1, "_i0"))), nmr = unique(c(dd[module == "nmr", column], unlist(qmap), paste0("p", 23652:23655, "_i0"))), 
        cmr = dd[module == "cmr", column])
    ids <- file.path(out, "cache/eligible_ids.private.txt")
    if (!file.exists(ids)) 
        publish(d$eid, ids, "txt")
    sp <- file.path(out, "cache/projection_spec.json")
    if (!file.exists(sp)) 
        publish(toJSON(spec, auto_unbox = FALSE), sp, "txt")
    projected <- file.path(out, "cache/projected")
    cmd <- c(shQuote(file.path(root, "f/project.py")), shQuote(src[["pheno"]]), shQuote(sp), shQuote(ids), shQuote(projected))
    code <- system2(Sys.getenv("PYTHON", "python3"), cmd)
    assert(code == 0L, "PROJECTION_FAILED: partials retained")
    proof <- read_json(file.path(projected, "complete.json"), simplifyVector = TRUE)
    for (i in seq_len(nrow(proof$files))) assert(identical(digest(file = file.path(projected, paste0(proof$files$module[i], 
        ".csv")), algo = "sha256"), proof$files$sha256[i]), "PROJECTION_HASH_MISMATCH")
    meta <- fread(file.path(projected, "metadata.csv"), data.table = FALSE, colClasses = "character")
    j <- match(d$eid, meta$eid)
    assert(!anyNA(j), "PROJECTED_LINK_FAILED")
    meta <- meta[j, , drop = FALSE]
    assert("p53_i0" %in% names(meta) && all(as.Date(meta$p53_i0) == as.Date(d$date_attend), na.rm = TRUE), "BASELINE_DATES_MISMATCH")
    prs <- list(T2D = "p26285", CAD = "p26227")
    for (k in names(prs)) d[[paste0("PRS_", k)]] <- if (prs[[k]] %in% names(meta)) 
        num(meta[[prs[[k]]]])
    else NA_real_
    for (i in 1:40) d[[paste0("PC", i)]] <- if (paste0("p22009_a", i) %in% names(meta)) 
        num(meta[[paste0("p22009_a", i)]])
    else NA_real_
    batch <- if ("p22000" %in% names(meta)) 
        num(meta$p22000)
    else rep(NA, nrow(d))
    batch[!batch %in% c(-11:-1, 1:95)] <- NA
    d$gen_batch <- factor(batch)
    core <- list(data = d, events = inc$events, baseline = inc$baseline, history = history, metadata = meta, kin = kin, 
        missing = proof$missing)
    if (!file.exists(file.path(out, "cache/olink_chunks_complete.rds"))) {
        op <- fread(src[["olink"]], data.table = FALSE)
        assert("eid" %in% names(op) && !anyDuplicated(op$eid), "OLINK_ID_INVALID")
        oo <- dd[module == "olink"]
        chunks <- split(seq_len(nrow(oo)), ceiling(seq_len(nrow(oo))/25L))
        for (nm in names(chunks)) {
            p <- file.path(out, "cache/olink", paste0("chunk_", nm, ".rds"))
            cols <- intersect(c("eid", oo$column[chunks[[nm]]]), names(op))
            if (!file.exists(p)) 
                cp_save(op[as.character(op$eid) %in% d$eid, cols, drop = FALSE], p)
            else cp_load(p)
        }
        cp_save(list(columns = names(op), missing = setdiff(oo$column, names(op))), file.path(out, "cache/olink_chunks_complete.rds"))
        rm(op)
        gc(FALSE)
    }
    cp_save(core, file.path(out, "cache/core.rds"))
    flow <- rbindlist(lapply(unlist(cfg$cohorts), function(c) rbindlist(lapply(unlist(cfg$exposures), function(e) data.table(cohort = c, 
        exposure = e, N = sum(is.finite(d[[e]]) & if (c == "S8") d$S <= 8 else TRUE))))))
    publish(flow, file.path(out, "EXPOSURE_COHORT_FLOW.csv"), "csv")
    log_event("PREPARED_FULL_S", N = nrow(d), S_gt8 = sum(d$S > 8))
}

load_feature <- function(core, row) {
    module <- row$module
    col <- row$column
    d <- core$data
    meta <- core$metadata
    extra <- character()
    source_qc <- "FROZEN_QC"
    cachekey <- module
    if (module == "olink") 
        cachekey <- paste0("olink_", ceiling(match(col, dict()[module == "olink", column])/25L))
    if (!exists(cachekey, envir = feature_cache, inherits = FALSE)) {
        old <- ls(feature_cache, pattern = "^olink_")
        if (length(old)) 
            rm(list = old, envir = feature_cache)
        if (module == "olink") 
            z <- cp_load(file.path(out, "cache/olink", paste0("chunk_", sub("olink_", "", cachekey), ".rds")))
        else if (module == "blood_pressure") 
            z <- meta
        else z <- fread(file.path(out, "cache/projected", paste0(module, ".csv")), data.table = FALSE)
        assert(!anyDuplicated(z$eid), "FEATURE_DUPLICATE_ID")
        assign(cachekey, z, envir = feature_cache)
    }
    z <- get(cachekey, envir = feature_cache)
    j <- match(d$eid, as.character(z$eid))
    if (module == "blood_pressure") {
        sc <- c("p4080_i0_a0", "p4080_i0_a1")
        dc <- c("p4079_i0_a0", "p4079_i0_a1")
        assert(all(c(sc, dc) %in% names(z)), "BP_FIELDS_MISSING")
        s <- sapply(z[sc], num)
        b <- sapply(z[dc], num)
        ok <- is.finite(s) & is.finite(b) & s > b & b > 0
        s[!ok] <- b[!ok] <- NA
        sb <- rowMeans(s, na.rm = TRUE)
        db <- rowMeans(b, na.rm = TRUE)
        v <- switch(col, SBP = sb, DBP = db, PP = sb - db, MAP = db + (sb - db)/3)[j]
    }
    else {
        assert(col %in% names(z), paste("FEATURE_COLUMN_MISSING", col))
        v <- num(z[[col]][j])
    }
    d$measurement_date <- as.Date(meta$p53_i0)
    if (module == "clinical") {
        qc <- row$qc_column
        assert(qc %in% names(z), "REPORTABILITY_FIELD_MISSING")
        q <- num(z[[qc]][j])
        v[is.na(q) | q != 1 | !is.finite(v) | v <= 0] <- NA
        dtcol <- paste0("p", num(sub("_i0", "", sub("p", "", col))) + 1, "_i0")
        assert(dtcol %in% names(z), "ASSAY_DATE_MISSING")
        d$assay_month <- factor(format(as.Date(z[[dtcol]][j]), "%Y-%m"))
        d$fasting <- num(meta$p74_i0)
        d$fasting[d$fasting < 0] <- NA
        extra <- c("assay_month", "fasting")
    }
    if (module == "nmr") {
        qcmap <- read_json(file.path(root, "config/nmr_qc_map.json"), simplifyVector = TRUE)
        qc <- c(unlist(qcmap[[col]]), paste0("p", 23652:23655, "_i0"))
        assert(all(qc %in% names(z)), "NMR_QC_FIELD_MISSING")
        bad <- rep(FALSE, nrow(z))
        for (k in qc) bad <- bad | (!is.na(z[[k]]) & nzchar(trimws(as.character(z[[k]]))))
        v[bad[j]] <- NA
        d$nmr_batch <- factor(meta$p20282_i0)
        d$nmr_phase <- factor(meta$p23650_i0)
        d$fasting <- num(meta$p74_i0)
        d$fasting[d$fasting < 0] <- NA
        extra <- c("nmr_batch", "nmr_phase", "fasting")
    }
    if (module == "cmr") {
        assert(all(c("p53_i2", "p21003_i2", "p54_i2") %in% names(meta)), "IMAGING_METADATA_MISSING")
        d$measurement_date <- as.Date(meta$p53_i2)
        d$age <- num(meta$p21003_i2)
        d$center <- factor(meta$p54_i2)
        d$assessment_year <- factor(format(d$measurement_date, "%Y"))
        if (isTRUE(row$indexed)) {
            assert(all(c("p50_i2", "p21002_i2") %in% names(meta)), "SAME_VISIT_BSA_MISSING")
            h <- num(meta$p50_i2)
            w <- num(meta$p21002_i2)
            bsa <- sqrt(h * w/3600)
            bsa[!is.finite(bsa) | h <= 0 | w <= 0] <- NA
            v <- v/bsa
        }
    }
    if (module == "olink") {
        assert("p30901_i0" %in% names(meta), "PLATE_ID_MISSING")
        mp <- read.delim(file.path(root, "config/olink_batch_number.dat"), colClasses = "character")
        canon <- function(x) sub("^0+", "", trimws(as.character(x)))
        assert(!anyDuplicated(canon(mp$PlateID)), "BATCH_MAP_DUPLICATE")
        d$olink_batch <- factor(mp$Batch[match(canon(meta$p30901_i0), canon(mp$PlateID))])
        extra <- "olink_batch"
        source_qc <- "PROVISIONAL_SOURCE_QC_UNVERIFIED"
    }
    v[!is.finite(v)] <- NA
    d$marker_original <- v
    if (row$transform == "log") {
        v[v <= 0] <- NA
        v <- log(v)
    }
    d$marker <- v
    list(data = d, extra = extra, source_qc = source_qc)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        prepare <- function() {
            if (file.exists(file.path(out, "cache/core.rds"))) return(invisible(cp_load(file.path(out, "cache/core.rds"))))
            src <- paths()
            assert(all(file.exists(src)), "REQUIRED_INPUT_MISSING")
            manifest_path <- file.path(out, "INPUT_SHA256.csv")
            if (!file.exists(manifest_path)) {
                ss <- file.info(src)
                h <- vapply(seq_along(src), function(i) {
                  log_event("HASH_SOURCE", role = names(src)[i])
                  digest(file = src[i], algo = "sha256")
                }, "")
                after <- file.info(src)
                assert(identical(ss$size, after$size) && identical(ss$mtime, after$mtime), "SOURCE_CHANGED_DURING_HASH")
                names(h) <- names(src)
                for (k in names(cfg$expected_sha256)) assert(identical(h[[k]], cfg$expected_sha256[[k]]), paste("FROZEN_HASH_MISMATCH", 
                  k))
                publish(data.table(role = names(src), path = unname(src), sha256 = unname(h), bytes = ss$size, 
                  mtime = as.numeric(ss$mtime)), manifest_path, "csv")
            } else {
                a <- fread(manifest_path)
                ss <- file.info(a$path)
                assert(all(ss$size == a$bytes & as.numeric(ss$mtime) == a$mtime), "SOURCE_CHANGED_SINCE_PREPARE")
            }
            inc <- cp_load(src[["incident"]], strict = FALSE)
            zz <- readRDS(src[["incident"]])
            m <- read_json(file.path(results, cfg$incident_root, "check_complete.json"))
            assert(identical(zz$signature, m$signature), "INCIDENT_SIGNATURE_MISMATCH")
            all <- as.data.frame(readRDS(src[["all"]]))
            d <- linear_covars(make_birth_order(all))
            d <- d[rowSums(!is.na(d[, unlist(cfg$exposures)])) > 0, , drop = FALSE]
            for (k in c("birth_year", "cnt_older", "cnt_brother", "cnt_sister")) assert(k %in% names(d), paste("MISSING_FIELD", 
                k))
            d$eid <- as.character(d$eid)
            ix <- match(d$eid, as.character(inc$data$eid))
            assert(!anyNA(ix), "INCIDENT_CACHE_DOES_NOT_COVER_FULL_EXPOSURE_COHORT")
            for (k in c("S", unlist(cfg$exposures), "age", "sex", "ethnicity", "center", "assessment_year")) assert(isTRUE(all.equal(as.character(d[[k]]), 
                as.character(inc$data[[k]][ix]))), paste("FROZEN_ADAPTER_MISMATCH", k))
            inc$events <- lapply(inc$events, function(v) v[ix])
            inc$baseline <- lapply(inc$baseline, function(v) v[ix])
            inc$data <- inc$data[ix, , drop = FALSE]
            d$date_death_uncapped <- inc$data$date_death_uncapped
            d$date_lost_registry <- inc$data$date_lost_registry
            d <- d[, unique(c("eid", "birth_year", "date_attend", "date_death_uncapped", "date_lost_registry", 
                "age", "sex", "ethnicity", "center", "assessment_year", "S", unlist(cfg$exposures))), drop = FALSE]
            p <- as.data.frame(readRDS(src[["phewas"]]))
            assert("id" %in% names(p) && !anyDuplicated(p$id), "PHEWAS_INVALID_ID")
            j <- match(d$eid, as.character(p$id))
            history <- lapply(unlist(cfg$phecodes), function(k) {
                assert(k %in% names(p), paste("PHECODE_MISSING", k))
                v <- p[[k]][j]
                assert(all(is.na(v) | v %in% c(0, 1)), "UNKNOWN_PHEWAS_STATUS")
                num(v)
            })
            names(history) <- unlist(cfg$phecodes)
            rm(p, all)
            gc(FALSE)
            kin <- bs_kin(src[["kinship"]])
            clusters <- bs_components(kin[kinship >= 0.0442])
            d$.cluster <- paste0("U", d$eid)
            j <- match(d$eid, clusters$eid)
            d$.cluster[!is.na(j)] <- paste0("R", clusters$component[j[!is.na(j)]])
            dd <- dict()
            qmap <- read_json(file.path(root, "config/nmr_qc_map.json"), simplifyVector = TRUE)
            spec <- list(metadata = unique(c("p53_i0", "p53_i2", "p21003_i2", "p54_i2", "p50_i2", "p21002_i2", 
                "p74_i0", "p20282_i0", "p23650_i0", "p30901_i0", "p26227", "p26285", "p22000", paste0("p22009_a", 
                  1:40), "p4080_i0_a0", "p4080_i0_a1", "p4079_i0_a0", "p4079_i0_a1")), clinical = unique(c(dd[module == 
                "clinical", column], dd[module == "clinical", qc_column], paste0("p", num(sub("_i0", "", sub("p", 
                "", dd[module == "clinical", column]))) + 1, "_i0"))), nmr = unique(c(dd[module == "nmr", column], 
                unlist(qmap), paste0("p", 23652:23655, "_i0"))), cmr = dd[module == "cmr", column])
            ids <- file.path(out, "cache/eligible_ids.private.txt")
            if (!file.exists(ids)) publish(d$eid, ids, "txt")
            sp <- file.path(out, "cache/projection_spec.json")
            if (!file.exists(sp)) publish(toJSON(spec, auto_unbox = FALSE), sp, "txt")
            projected <- file.path(out, "cache/projected")
            cmd <- c(shQuote(file.path(root, "f/project.py")), shQuote(src[["pheno"]]), shQuote(sp), shQuote(ids), 
                shQuote(projected))
            code <- system2(Sys.getenv("PYTHON", "python3"), cmd)
            assert(code == 0L, "PROJECTION_FAILED: partials retained")
            proof <- read_json(file.path(projected, "complete.json"), simplifyVector = TRUE)
            for (i in seq_len(nrow(proof$files))) assert(identical(digest(file = file.path(projected, paste0(proof$files$module[i], 
                ".csv")), algo = "sha256"), proof$files$sha256[i]), "PROJECTION_HASH_MISMATCH")
            meta <- fread(file.path(projected, "metadata.csv"), data.table = FALSE, colClasses = "character")
            j <- match(d$eid, meta$eid)
            assert(!anyNA(j), "PROJECTED_LINK_FAILED")
            meta <- meta[j, , drop = FALSE]
            assert("p53_i0" %in% names(meta) && all(as.Date(meta$p53_i0) == as.Date(d$date_attend), na.rm = TRUE), 
                "BASELINE_DATES_MISMATCH")
            prs <- list(T2D = "p26285", CAD = "p26227")
            for (k in names(prs)) d[[paste0("PRS_", k)]] <- if (prs[[k]] %in% names(meta)) num(meta[[prs[[k]]]]) else NA_real_
            for (i in 1:40) d[[paste0("PC", i)]] <- if (paste0("p22009_a", i) %in% names(meta)) num(meta[[paste0("p22009_a", 
                i)]]) else NA_real_
            batch <- if ("p22000" %in% names(meta)) num(meta$p22000) else rep(NA, nrow(d))
            batch[!batch %in% c(-11:-1, 1:95)] <- NA
            d$gen_batch <- factor(batch)
            core <- list(data = d, events = inc$events, baseline = inc$baseline, history = history, metadata = meta, 
                kin = kin, missing = proof$missing)
            if (!file.exists(file.path(out, "cache/olink_chunks_complete.rds"))) {
                op <- fread(src[["olink"]], data.table = FALSE)
                assert("eid" %in% names(op) && !anyDuplicated(op$eid), "OLINK_ID_INVALID")
                oo <- dd[module == "olink"]
                chunks <- split(seq_len(nrow(oo)), ceiling(seq_len(nrow(oo))/25L))
                for (nm in names(chunks)) {
                  p <- file.path(out, "cache/olink", paste0("chunk_", nm, ".rds"))
                  cols <- intersect(c("eid", oo$column[chunks[[nm]]]), names(op))
                  if (!file.exists(p)) cp_save(op[as.character(op$eid) %in% d$eid, cols, drop = FALSE], p) else cp_load(p)
                }
                cp_save(list(columns = names(op), missing = setdiff(oo$column, names(op))), file.path(out, "cache/olink_chunks_complete.rds"))
                rm(op)
                gc(FALSE)
            }
            cp_save(core, file.path(out, "cache/core.rds"))
            flow <- rbindlist(lapply(unlist(cfg$cohorts), function(c) rbindlist(lapply(unlist(cfg$exposures), function(e) data.table(cohort = c, 
                exposure = e, N = sum(is.finite(d[[e]]) & if (c == "S8") d$S <= 8 else TRUE))))))
            publish(flow, file.path(out, "EXPOSURE_COHORT_FLOW.csv"), "csv")
            log_event("PREPARED_FULL_S", N = nrow(d), S_gt8 = sum(d$S > 8))
        }
        load_feature <- function(core, row) {
            module <- row$module
            col <- row$column
            d <- core$data
            meta <- core$metadata
            extra <- character()
            source_qc <- "FROZEN_QC"
            cachekey <- module
            if (module == "olink") cachekey <- paste0("olink_", ceiling(match(col, dict()[module == "olink", column])/25L))
            if (!exists(cachekey, envir = feature_cache, inherits = FALSE)) {
                old <- ls(feature_cache, pattern = "^olink_")
                if (length(old)) rm(list = old, envir = feature_cache)
                if (module == "olink") z <- cp_load(file.path(out, "cache/olink", paste0("chunk_", sub("olink_", 
                  "", cachekey), ".rds"))) else if (module == "blood_pressure") z <- meta else z <- fread(file.path(out, 
                  "cache/projected", paste0(module, ".csv")), data.table = FALSE)
                assert(!anyDuplicated(z$eid), "FEATURE_DUPLICATE_ID")
                assign(cachekey, z, envir = feature_cache)
            }
            z <- get(cachekey, envir = feature_cache)
            j <- match(d$eid, as.character(z$eid))
            if (module == "blood_pressure") {
                sc <- c("p4080_i0_a0", "p4080_i0_a1")
                dc <- c("p4079_i0_a0", "p4079_i0_a1")
                assert(all(c(sc, dc) %in% names(z)), "BP_FIELDS_MISSING")
                s <- sapply(z[sc], num)
                b <- sapply(z[dc], num)
                ok <- is.finite(s) & is.finite(b) & s > b & b > 0
                s[!ok] <- b[!ok] <- NA
                sb <- rowMeans(s, na.rm = TRUE)
                db <- rowMeans(b, na.rm = TRUE)
                v <- switch(col, SBP = sb, DBP = db, PP = sb - db, MAP = db + (sb - db)/3)[j]
            } else {
                assert(col %in% names(z), paste("FEATURE_COLUMN_MISSING", col))
                v <- num(z[[col]][j])
            }
            d$measurement_date <- as.Date(meta$p53_i0)
            if (module == "clinical") {
                qc <- row$qc_column
                assert(qc %in% names(z), "REPORTABILITY_FIELD_MISSING")
                q <- num(z[[qc]][j])
                v[is.na(q) | q != 1 | !is.finite(v) | v <= 0] <- NA
                dtcol <- paste0("p", num(sub("_i0", "", sub("p", "", col))) + 1, "_i0")
                assert(dtcol %in% names(z), "ASSAY_DATE_MISSING")
                d$assay_month <- factor(format(as.Date(z[[dtcol]][j]), "%Y-%m"))
                d$fasting <- num(meta$p74_i0)
                d$fasting[d$fasting < 0] <- NA
                extra <- c("assay_month", "fasting")
            }
            if (module == "nmr") {
                qcmap <- read_json(file.path(root, "config/nmr_qc_map.json"), simplifyVector = TRUE)
                qc <- c(unlist(qcmap[[col]]), paste0("p", 23652:23655, "_i0"))
                assert(all(qc %in% names(z)), "NMR_QC_FIELD_MISSING")
                bad <- rep(FALSE, nrow(z))
                for (k in qc) bad <- bad | (!is.na(z[[k]]) & nzchar(trimws(as.character(z[[k]]))))
                v[bad[j]] <- NA
                d$nmr_batch <- factor(meta$p20282_i0)
                d$nmr_phase <- factor(meta$p23650_i0)
                d$fasting <- num(meta$p74_i0)
                d$fasting[d$fasting < 0] <- NA
                extra <- c("nmr_batch", "nmr_phase", "fasting")
            }
            if (module == "cmr") {
                assert(all(c("p53_i2", "p21003_i2", "p54_i2") %in% names(meta)), "IMAGING_METADATA_MISSING")
                d$measurement_date <- as.Date(meta$p53_i2)
                d$age <- num(meta$p21003_i2)
                d$center <- factor(meta$p54_i2)
                d$assessment_year <- factor(format(d$measurement_date, "%Y"))
                if (isTRUE(row$indexed)) {
                  assert(all(c("p50_i2", "p21002_i2") %in% names(meta)), "SAME_VISIT_BSA_MISSING")
                  h <- num(meta$p50_i2)
                  w <- num(meta$p21002_i2)
                  bsa <- sqrt(h * w/3600)
                  bsa[!is.finite(bsa) | h <= 0 | w <= 0] <- NA
                  v <- v/bsa
                }
            }
            if (module == "olink") {
                assert("p30901_i0" %in% names(meta), "PLATE_ID_MISSING")
                mp <- read.delim(file.path(root, "config/olink_batch_number.dat"), colClasses = "character")
                canon <- function(x) sub("^0+", "", trimws(as.character(x)))
                assert(!anyDuplicated(canon(mp$PlateID)), "BATCH_MAP_DUPLICATE")
                d$olink_batch <- factor(mp$Batch[match(canon(meta$p30901_i0), canon(mp$PlateID))])
                extra <- "olink_batch"
                source_qc <- "PROVISIONAL_SOURCE_QC_UNVERIFIED"
            }
            v[!is.finite(v)] <- NA
            d$marker_original <- v
            if (row$transform == "log") {
                v[v <= 0] <- NA
                v <- log(v)
            }
            d$marker <- v
            list(data = d, extra = extra, source_qc = source_qc)
        }
    }), envir = context)
}

