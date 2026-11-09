# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        args <- commandArgs(TRUE)
        if (length(args) < 2L) stop("Usage: run.R resolved_config.json stage [branch]")
        cfg <- jsonlite::fromJSON(args[1], simplifyVector = TRUE)
        stage <- args[2]
        out <- cfg$out
        options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.poly"), warn = 1)
        set.seed(cfg$seed)
        for (k in c("jsonlite", "digest", "dplyr", "logistf", "ggplot2")) if (!requireNamespace(k, quietly = TRUE)) stop("DEPENDENCY_MISSING: ", 
            k)
        if (as.character(packageVersion("logistf")) != cfg$logistf_version) stop("LOGISTF_VERSION_MISMATCH")
        source(file.path(cfg$code_root, "frozen", "phewas.main.helpers.R"))
        source(file.path(cfg$code_root, "frozen", "method.R"))
        for (n in c("00.common.R", "20.formal.R", "30.shape.R", "40.report.R")) source(file.path(cfg$code_root, 
            "f", n))
        dir.create(file.path(out, "logs"), showWarnings = FALSE, recursive = TRUE)
        capture.output(sessionInfo(), file = file.path(out, "logs", paste0("software_", Sys.getpid(), ".txt")))
        log_event("STAGE_START", stage_name = stage, output = out, method = cfg$method_id)
        if (stage == "check") {
            reg <- main_registry(cfg)
            write_new_csv(data.frame(branch = linear_exposures, registered = 1582, planned = vapply(linear_exposures, 
                function(b) sum(reg$planned[reg$branch == b]), integer(1))), file.path(out, "audit", "registry_check.csv"))
            log_event("CHECK_PASS_NO_PARTICIPANT_MODELS", registered = nrow(reg))
        } else if (stage == "test") {
            for (p in list.files(file.path(cfg$code_root, "tests"), pattern = "^test_.*[.]R$", full.names = TRUE)) source(p)
            for (n in c("test_common", "test_formal", "test_shape", "test_report")) {
                if (!exists(n, mode = "function")) stop("TEST_ENTRY_MISSING: ", n)
                log_event("SYNTHETIC_TEST", test = n)
                a <- get(n)()
                if (is.data.frame(a)) {
                  print(a)
                  write_new_csv(a, file.path(out, "logs", paste0(n, "_", Sys.getpid(), ".csv")))
                  if ("status" %in% names(a) && any(a$status != "PASS")) stop("SYNTHETIC_TEST_FAILED: ", n)
                }
            }
            log_event("ALL_SYNTHETIC_TESTS_PASS")
        } else if (stage == "audit56") {
            old_firth_audit(cfg, out)
        } else if (stage == "benchmark") {
            old_firth_audit(cfg, out)
            run_benchmark(cfg, out)
        } else if (stage == "phewas") {
            old_firth_audit(cfg, out)
            ans <- run_formal(cfg, out, args[3])
            if (identical(ans$status, "RESOURCE_HOLD")) quit(status = 78L)
        } else if (stage == "shape") {
            ans <- run_shape(cfg, out)
            if (identical(ans$status, "RESOURCE_HOLD")) quit(status = 78L)
        } else if (stage == "report") {
            run_report(cfg, out)
        } else stop("UNSUPPORTED_STAGE")
        log_event("STAGE_FINISH", stage_name = stage)
    }), envir = context)
}

