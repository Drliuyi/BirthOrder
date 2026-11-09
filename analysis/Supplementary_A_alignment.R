# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

num <- function(x) suppressWarnings(as.numeric(x))

blank <- function(x) is.na(x) | x %in% c("", "NA")

close_relative <- function(a, b, tolerance = 1e-12) {
    length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)) && all(abs(a - b) <= tolerance * pmax(abs(a), 
        abs(b), .Machine$double.xmin))
}

read_character_csv <- function(path) read.csv(path, colClasses = "character", check.names = FALSE, na.strings = NULL, 
    stringsAsFactors = FALSE)

formula_clean <- function(x) gsub("[[:space:]]+", "", x)

sha <- function(p) digest::digest(file = p, algo = "sha256")

read_text <- function(p) paste(readLines(p, warn = FALSE), collapse = "\n")

checks <- list()

check <- function(name, ok, evidence) {
    assert(ok, paste("FAILED", name, evidence))
    checks[[length(checks) + 1L]] <<- data.frame(check = name, pass = TRUE, evidence = evidence)
}

family <- c("packyears", "age_education", "qualification", "income", "townsend", "mobility")

p_restored <- p_original

fmt <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))

original_all <- a

write_new <- function(d, name) {
    p <- file.path(OUTPUT_ROOT, name)
    assert(!file.exists(p), paste("Refusing overwrite:", p))
    write.csv(d, p, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        options(stringsAsFactors = FALSE, digits = 17)
        assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
        num <- function(x) suppressWarnings(as.numeric(x))
        blank <- function(x) is.na(x) | x %in% c("", "NA")
        close_relative <- function(a, b, tolerance = 1e-12) {
            length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)) && all(abs(a - b) <= tolerance * pmax(abs(a), 
                abs(b), .Machine$double.xmin))
        }
        read_character_csv <- function(path) read.csv(path, colClasses = "character", check.names = FALSE, na.strings = NULL, 
            stringsAsFactors = FALSE)
        formula_clean <- function(x) gsub("[[:space:]]+", "", x)
        expected_formula <- "M~border_relative+age+sex+ethnicity+center+assessment_year+factor(S)"
        DIR0 <- Sys.getenv("DIR0", if (.Platform$OS.type == "windows") birth_private_setting("PRIVATE_INPUT") else "/mnt/d")
        SCRIPT_ROOT <- Sys.getenv("APATH_SCRIPT_ROOT", file.path(DIR0, "scripts"))
        INPUT_ROOT <- Sys.getenv("APATH_INPUT_ROOT", file.path(DIR0, "analysis", "ukb", "birth"))
        OUTPUT_ROOT <- Sys.getenv("APATH_OUTPUT_ROOT", file.path(INPUT_ROOT, "border_apath_M6_restored_20261022_v1"))
        cat("Resolved script root:", SCRIPT_ROOT, "\nResolved aggregate input root:", INPUT_ROOT, "\nResolved output root:", 
            OUTPUT_ROOT, "\n")
        assert(dir.exists(SCRIPT_ROOT), paste("Missing script root:", SCRIPT_ROOT))
        assert(dir.exists(INPUT_ROOT), paste("Missing aggregate input root:", INPUT_ROOT))
        assert(!file.exists(OUTPUT_ROOT), paste("Refusing existing output:", OUTPUT_ROOT))
        assert(dir.exists(dirname(OUTPUT_ROOT)), paste("Missing output parent:", dirname(OUTPUT_ROOT)))
        for (package in c("digest", "jsonlite")) {
            assert(requireNamespace(package, quietly = TRUE), paste("Missing required package:", package))
        }
        supplement <- file.path(INPUT_ROOT, "border_pathway_supplement_phase1_20261021_122704")
        recovery <- file.path(INPUT_ROOT, "border_mobility_recovery_S8_v1_20261021")
        paths <- c(original_A6 = file.path(supplement, "A_NEW_ALIGNMENT_6.csv"), recovered_mobility_A = file.path(recovery, 
            "MOBILITY_A_PATH.csv"), original_config = file.path(supplement, "code", "config", "analysis_config.json"), 
            original_variables = file.path(supplement, "code", "config", "variables.json"), original_prepare = file.path(supplement, 
                "code", "f", "1.prepare.R"), original_merge = file.path(supplement, "code", "f", "4.merge.R"), 
            original_sample_flow = file.path(supplement, "SAMPLE_FLOW.csv"), recovered_derivation = file.path(recovery, 
                "DISTANCE_DERIVATION_MANIFEST.json"), recovery_prepare = file.path(SCRIPT_ROOT, "birth", "border_mobility_recovery_S8_v1", 
                "f", "1.prepare.R"), recovery_config = file.path(SCRIPT_ROOT, "birth", "border_mobility_recovery_S8_v1", 
                "config", "analysis_config.json"), this_script = file.path(SCRIPT_ROOT, "birth", "border_apath_M6_restored_20261022_v1", 
                "restore_M6.R"))
        assert(all(file.exists(paths)), paste("Missing frozen aggregate/provenance input:", paste(paths[!file.exists(paths)], 
            collapse = "; ")))
        paths <- vapply(paths, normalizePath, "", mustWork = TRUE)
        sha <- function(p) digest::digest(file = p, algo = "sha256")
        hash_before <- vapply(paths, sha, "")
        a <- read_character_csv(paths[["original_A6"]])
        g <- read_character_csv(paths[["recovered_mobility_A"]])
        cfg <- jsonlite::fromJSON(paths[["original_config"]])
        vars <- jsonlite::fromJSON(paths[["original_variables"]])
        deriv <- jsonlite::fromJSON(paths[["recovered_derivation"]])
        gcfg <- jsonlite::fromJSON(paths[["recovery_config"]])
        flow <- read_character_csv(paths[["original_sample_flow"]])
        read_text <- function(p) paste(readLines(p, warn = FALSE), collapse = "\n")
        prepare_text <- read_text(paths[["original_prepare"]])
        recovery_prepare_text <- read_text(paths[["recovery_prepare"]])
        merge_text <- read_text(paths[["original_merge"]])
        checks <- list()
        check <- function(name, ok, evidence) {
            assert(ok, paste("FAILED", name, evidence))
            checks[[length(checks) + 1L]] <<- data.frame(check = name, pass = TRUE, evidence = evidence)
        }
        family <- c("packyears", "age_education", "qualification", "income", "townsend", "mobility")
        check("original_family_exact_six", nrow(a) == 6L && !anyDuplicated(a$marker) && setequal(a$marker, family), 
            paste(family, collapse = ";"))
        check("original_preregistered_A6", cfg$M_alignment_A == 6L && setequal(vars$id[vars$A_action == "NEW_ALIGNMENT_6"], 
            family), "Saved config and variable registry")
        check("original_main_A_only", all(a$analysis == "A") && all(blank(a$phecode)) && all(num(a$M_alignment) == 
            6L) && all(num(a$M_reported) == 6L) && all(a$q_family == "SUPPLEMENT_NEW_ALIGNMENT_A_M6"), "Independent main A; original family M6")
        check("original_S2_8_exposure_roster", grepl("d$S%in%2:8", prepare_text, fixed = TRUE) && grepl("d$border_relative>=0&d$border_relative<=1", 
            prepare_text, fixed = TRUE) && num(flow$N[flow$step == "S2_8_frozen_roster"]) == 143736L, "Saved producer eligibility and aggregate roster N143736; no participant read")
        check("original_fixed_M6_BH_rule", grepl("a6$q_alignment <- psm_fixed_bh(a6$p, psm_reliable(a6), 6L)", 
            merge_text, fixed = TRUE) && grepl("p.adjust(p[ok], method = \"BH\", n = M)", merge_text, fixed = TRUE), 
            "Frozen producer fixed-M6 reliable finite-P BH rule")
        mi <- match("mobility", a$marker)
        oi <- setdiff(seq_len(nrow(a)), mi)
        check("five_original_reliable_A_fits", all(a$status[oi] == "SUCCESS") && all(is.finite(num(a$p[oi]))) && 
            close_relative(num(a$p[oi]), num(a$p_reliable[oi])) && all(num(a$rank[oi]) == num(a$columns[oi])) && 
            all(num(a$rank[oi]) == 25L) && all(blank(a$removed_terms[oi])) && all(blank(a$warning[oi])) && all(blank(a$reason[oi])), 
            "Five SUCCESS, reliable P, rank=columns=25, no removed terms/warnings/reasons")
        check("original_five_same_formula", all(formula_clean(a$formula[oi]) == expected_formula) && all(formula_clean(a$a_formula[oi]) == 
            expected_formula), expected_formula)
        check("old_mobility_was_blocked_not_negative", a$status[mi] == "GEO_SOURCE_BLOCKED" && blank(a$p[mi]) && 
            blank(a$q_alignment_M6[mi]), "Original slot remains untouched in source")
        check("current_mobility_formal_reliable_main_A", nrow(g) == 1L && g$marker == "mobility" && g$analysis == 
            "A" && blank(g$phecode) && g$status == "SUCCESS" && g$branch_status == "FORMAL_BASELINE_I0_REDERIVED" && 
            g$source_approved == "TRUE" && close_relative(num(g$p), num(g$p_reliable)) && num(g$rank) == 25L && 
            num(g$columns) == 25L && blank(g$removed_terms) && blank(g$warning) && blank(g$reason), "Approved recovery SUCCESS; reliable P; full rank25; no warnings/removed terms")
        check("all_current_models_formula_exposure_match", all(a$exposure_scale == "border_relative_0_to_1_not_Z") && 
            g$exposure_scale == "border_relative_0_to_1_not_Z" && formula_clean(g$formula) == expected_formula && 
            formula_clean(g$a_formula) == expected_formula, "Exposure border_relative 0-to-1, not standardized; frozen linear-age adjustment")
        check("current_measurement_and_exact_sample", g$measurement_version == "BASELINE_I0_GB_NET_DISTANCE_KM_20261021" && 
            deriv$measurement_version == g$measurement_version && deriv$instance == 0L && isTRUE(deriv$source_approved) && 
            deriv$source_status == g$branch_status && num(g$N) == 123611L && deriv$valid_base_complete_N == num(g$N) && 
            deriv$prepared_N == 143736L && gcfg$measurement_instance == 0L && grepl("all(d$S%in%2:8)", recovery_prepare_text, 
            fixed = TRUE), "Formal baseline-i0 GB net distance; S2-8 inherited frozen roster N143736; exact main A N123611")
        check("current_mobility_scale_CI_valid", num(g$marker_sd) > 0 && close_relative(num(g$beta_raw)/num(g$marker_sd), 
            num(g$beta_SD)) && close_relative(num(g$CI_low_raw)/num(g$marker_sd), num(g$CI_low_SD)) && close_relative(num(g$CI_high_raw)/num(g$marker_sd), 
            num(g$CI_high_SD)) && num(g$CI_low_SD) < num(g$beta_SD) && num(g$CI_high_SD) > num(g$beta_SD), "Saved same-main-A outcome SD conversion; no model or scale change")
        p_original <- num(a$p)
        p_restored <- p_original
        p_restored[mi] <- num(g$p)
        check("restored_complete_six_valid_P", all(is.finite(p_restored)) && all(p_restored >= 0 & p_restored <= 
            1), "Six reliable saved primary-A P values")
        q_restored <- p.adjust(p_restored, method = "BH", n = 6L)
        ord <- order(p_restored, a$marker)
        q_sorted <- pmin(1, rev(cummin(rev(6 * p_restored[ord]/seq_len(6L)))))
        q_explicit <- numeric(6L)
        q_explicit[ord] <- q_sorted
        check("R_BH_matches_independent_sorted_formula", close_relative(q_restored, q_explicit, 1e-14), "p.adjust(BH,n=6) vs min_{j>=i}(6*P_(j)/j); no estimator used")
        check("old_saved_q_reproduces_old_fixed_M6", close_relative(p.adjust(p_original[oi], method = "BH", n = 6L), 
            num(a$q_alignment_M6[oi])), "Original five saved q reproduced under original missing-slot n6 rule")
        check("other_five_P_unchanged", identical(p_restored[oi], p_original[oi]), "Only restored mobility P differs from original vector")
        restored_source <- rep(paths[["original_A6"]], 6L)
        restored_source[mi] <- paths[["recovered_mobility_A"]]
        restored_source_hash <- rep(hash_before[["original_A6"]], 6L)
        restored_source_hash[mi] <- hash_before[["recovered_mobility_A"]]
        N_restored <- a$N
        N_restored[mi] <- g$N
        formula_restored <- a$formula
        formula_restored[mi] <- g$formula
        sample_hash_restored <- a$sample_hash
        sample_hash_restored[mi] <- g$sample_hash
        rank_restored <- a$rank
        rank_restored[mi] <- g$rank
        column_count_restored <- a$columns
        column_count_restored[mi] <- g$columns
        p_restored_text <- a$p
        p_restored_text[mi] <- g$p
        fmt <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))
        q_old <- num(a$q_alignment_M6)
        primary <- data.frame(marker = a$marker, p_original = a$p, q_original = a$q_alignment_M6, p_restored = p_restored_text, 
            q_BH_M6_restored = fmt(q_restored), q_family_M6_restored = "PRIMARY_A_M6_RESTORED_20261022", family_M_restored = 6L, 
            N_restored = N_restored, source_file = restored_source, source_sha256 = restored_source_hash, original_source_file = paths[["original_A6"]], 
            original_source_sha256 = hash_before[["original_A6"]], q_BH_M6_explicit_check = fmt(q_explicit), BH_absolute_check_difference = fmt(q_restored - 
                q_explicit), q_restored_minus_original = fmt(q_restored - q_old), q_restored_over_original = fmt(q_restored/q_old), 
            q_change_status = rep("", 6L), significant_q05_restored = q_restored < 0.05, formula_restored = formula_restored, 
            exposure_scale_restored = "border_relative_0_to_1_not_Z", S_range_restored = "2-8", measurement_version_restored = ifelse(a$marker == 
                "mobility", g$measurement_version, "baseline_instance_0"), sample_hash_restored = sample_hash_restored, 
            rank_restored = rank_restored, columns_restored = column_count_restored, status_restored = "SUCCESS", 
            p_source_operation = ifelse(a$marker == "mobility", "FORMAL_BASELINE_I0_MAIN_A_RECOVERY", "ORIGINAL_FIVE_MAIN_A_P_UNCHANGED"), 
            correction_operation = "WHOLE_ORIGINAL_SIX_TEST_FAMILY_BH_N6_NO_REFIT")
        primary$q_change_status <- vapply(seq_len(6L), function(i) {
            if (is.na(q_old[i])) "RESTORED_PREVIOUSLY_MISSING" else if (close_relative(q_restored[i], q_old[i])) "NUMERICALLY_UNCHANGED" else "UPDATED_WHOLE_FAMILY_AFTER_RECOVERY"
        }, "")
        original_all <- a
        names(original_all) <- paste0("original__", names(a))
        primary <- cbind(primary, original_all)
        hash_after <- vapply(paths, sha, "")
        check("all_frozen_inputs_unchanged", identical(hash_before, hash_after), "SHA256 before and after exact equality")
        manifest <- data.frame(input_role = names(paths), source_file = unname(paths), sha256_before = unname(hash_before), 
            sha256_after = unname(hash_after), unchanged = TRUE, bytes = file.info(paths)$size)
        validation <- do.call(rbind, checks)
        assert(!file.exists(OUTPUT_ROOT), paste("Refusing existing output:", OUTPUT_ROOT))
        assert(dir.create(OUTPUT_ROOT), paste("Could not create output:", OUTPUT_ROOT))
        assert(dir.create(file.path(OUTPUT_ROOT, "data")), "Could not create isolated data output")
        write_new <- function(d, name) {
            p <- file.path(OUTPUT_ROOT, name)
            assert(!file.exists(p), paste("Refusing overwrite:", p))
            write.csv(d, p, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
        }
        write_new(primary, "data/PRIMARY_A_M6_RESTORED.csv")
        write_new(manifest, "INPUT_HASH_MANIFEST.csv")
        write_new(validation, "VALIDATION.csv")
        write_new(primary[, c("marker", "p_original", "q_original", "p_restored", "q_BH_M6_restored", "q_restored_minus_original", 
            "q_restored_over_original", "q_change_status", "significant_q05_restored")], "M6_OLD_TO_RESTORED_DIFFERENCES.csv")
        writeLines(capture.output(sessionInfo()), file.path(OUTPUT_ROOT, "SOFTWARE_SESSION.txt"))
        summary <- list(version = "PRIMARY_A_M6_RESTORED_20261022", family_M = 6L, markers = a$marker, old_valid_P = 5L, 
            restored_valid_P = 6L, significant_q05_restored = sum(q_restored < 0.05), tests_passed = nrow(validation), 
            formula = expected_formula, exposure = "border_relative 0 to 1, not Z standardized", S_range = "2-8", 
            inputs_unchanged = TRUE, model_fits = 0L, bootstrap_replicates = 0L, other_families_modified = FALSE, 
            original_tables_modified = FALSE, script_root = SCRIPT_ROOT, input_root = INPUT_ROOT, output_root = OUTPUT_ROOT, 
            primary_output = file.path(OUTPUT_ROOT, "data", "PRIMARY_A_M6_RESTORED.csv"))
        jsonlite::write_json(summary, file.path(OUTPUT_ROOT, "COMPLETION.json"), pretty = TRUE, auto_unbox = TRUE)
        cat("Completed whole-family M6 recovery without model fitting:\n")
        print(primary[, c("marker", "q_original", "q_BH_M6_restored", "q_change_status")], row.names = FALSE)
        cat("Primary output:", summary$primary_output, "\n")
    }), envir = context)
}

