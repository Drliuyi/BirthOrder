# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

bi_stage_nmr <- function(ctx) bi_feature_stage(ctx, "nmr")

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        bi_stage_nmr <- function(ctx) bi_feature_stage(ctx, "nmr")
    }), envir = context)
}

