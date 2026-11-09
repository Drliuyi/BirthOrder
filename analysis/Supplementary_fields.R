# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

v2_supplemental_fields <- function(ctx) {
    v2_cp(ctx, "supplemental_fields", function() {
        cfg <- ctx$cfg$supplemental_fields
        p <- ctx$cfg$inputs$raw_pheno
        v2_assert(ctx$real_authorized && Sys.info()[["sysname"]] == "Linux", "WINPC_RAW_PROJECTION_ONLY")
        before <- file.info(p)
        v2_assert(identical(v2_filehash(p), cfg$source_sha256), "SUPPLEMENTAL_RAW_SOURCE_CHANGED")
        dest <- file.path(ctx$out, "private", paste0("supplemental_", Sys.getpid(), ".csv.partial"))
        audit <- paste0(dest, ".audit.json")
        log <- paste0(dest, ".log")
        v2_assert(!file.exists(dest), "EXISTING_PROJECTION_PARTIAL_REQUIRES_REVIEW")
        args <- c(shQuote(file.path(ctx$root, "vendor/streaming_tsv_select.py")), "--input", shQuote(p), "--columns-json", 
            shQuote(jsonlite::toJSON(unlist(cfg$columns), auto_unbox = FALSE)), "--dialect", cfg$dialect, "--audit-json", 
            shQuote(audit))
        rc <- system2("python3", args, stdout = dest, stderr = log)
        v2_assert(rc == 0, paste("SUPPLEMENTAL_PROJECTION_FAILED", log))
        evidence <- jsonlite::fromJSON(audit)
        v2_assert(isTRUE(evidence$full_scan) && evidence$status == "PASS", "INCOMPLETE_PROJECTION")
        after <- file.info(p)
        v2_assert(identical(before$size, after$size) && identical(before$mtime, after$mtime), "SOURCE_CHANGED_DURING_PROJECTION")
        z <- data.table::fread(dest, data.table = FALSE, colClasses = "character", na.strings = c("", "NA"))
        v2_ids(z$eid)
        v2_assert(setequal(names(z), unlist(cfg$columns)), "SUPPLEMENTAL_FIELD_SCHEMA")
        v2_csv(ctx, "SUPPLEMENTAL_SOURCE_AUDIT.csv", data.frame(source = p, sha256 = cfg$source_sha256, projection_sha256 = v2_filehash(dest), 
            rows = nrow(z), columns = ncol(z), dialect = cfg$dialect, full_scan = TRUE))
        Sys.chmod(dest, "0600")
        z
    })
}

v2_restriction_specs <- function(raw) {
    clean <- function(x) {
        v <- suppressWarnings(as.numeric(x))
        v2_assert(all(is.na(x) | (!is.na(v) & v %in% c(0, 1, -1, -3))), "BASELINE_100349_CODE_INVALID")
        !is.na(v) & v == 0
    }
    a <- clean(raw$p1767_i0)
    m <- clean(raw$p1777_i0)
    list(known_not_adopted = a, known_not_multiple_birth = m, known_not_adopted_and_not_multiple_birth = a & m)
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        v2_supplemental_fields <- function(ctx) {
            v2_cp(ctx, "supplemental_fields", function() {
                cfg <- ctx$cfg$supplemental_fields
                p <- ctx$cfg$inputs$raw_pheno
                v2_assert(ctx$real_authorized && Sys.info()[["sysname"]] == "Linux", "WINPC_RAW_PROJECTION_ONLY")
                before <- file.info(p)
                v2_assert(identical(v2_filehash(p), cfg$source_sha256), "SUPPLEMENTAL_RAW_SOURCE_CHANGED")
                dest <- file.path(ctx$out, "private", paste0("supplemental_", Sys.getpid(), ".csv.partial"))
                audit <- paste0(dest, ".audit.json")
                log <- paste0(dest, ".log")
                v2_assert(!file.exists(dest), "EXISTING_PROJECTION_PARTIAL_REQUIRES_REVIEW")
                args <- c(shQuote(file.path(ctx$root, "vendor/streaming_tsv_select.py")), "--input", shQuote(p), 
                  "--columns-json", shQuote(jsonlite::toJSON(unlist(cfg$columns), auto_unbox = FALSE)), "--dialect", 
                  cfg$dialect, "--audit-json", shQuote(audit))
                rc <- system2("python3", args, stdout = dest, stderr = log)
                v2_assert(rc == 0, paste("SUPPLEMENTAL_PROJECTION_FAILED", log))
                evidence <- jsonlite::fromJSON(audit)
                v2_assert(isTRUE(evidence$full_scan) && evidence$status == "PASS", "INCOMPLETE_PROJECTION")
                after <- file.info(p)
                v2_assert(identical(before$size, after$size) && identical(before$mtime, after$mtime), "SOURCE_CHANGED_DURING_PROJECTION")
                z <- data.table::fread(dest, data.table = FALSE, colClasses = "character", na.strings = c("", "NA"))
                v2_ids(z$eid)
                v2_assert(setequal(names(z), unlist(cfg$columns)), "SUPPLEMENTAL_FIELD_SCHEMA")
                v2_csv(ctx, "SUPPLEMENTAL_SOURCE_AUDIT.csv", data.frame(source = p, sha256 = cfg$source_sha256, 
                  projection_sha256 = v2_filehash(dest), rows = nrow(z), columns = ncol(z), dialect = cfg$dialect, 
                  full_scan = TRUE))
                Sys.chmod(dest, "0600")
                z
            })
        }
        v2_restriction_specs <- function(raw) {
            clean <- function(x) {
                v <- suppressWarnings(as.numeric(x))
                v2_assert(all(is.na(x) | (!is.na(v) & v %in% c(0, 1, -1, -3))), "BASELINE_100349_CODE_INVALID")
                !is.na(v) & v == 0
            }
            a <- clean(raw$p1767_i0)
            m <- clean(raw$p1777_i0)
            list(known_not_adopted = a, known_not_multiple_birth = m, known_not_adopted_and_not_multiple_birth = a & 
                m)
        }
    }), envir = context)
}

