# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

nmr_install_product_runtime <- function(methods) {
    stopifnot(is.environment(methods), is.function(methods$omics_product_stream))
    f <- methods$omics_product_stream
    replacements <- 0L
    declarations <- 0L
    domain_definition <- quote(draw_domain <- function(b, u, g, label) {
        M <- deterministic(function() draw_joint(b = b, C = cc, uniforms = u, requested_domain = g))
        cols <- domains[[g]]
        if (!is.data.frame(M) || nrow(M) != nrow(cc) || !all(cols %in% names(M))) new121_abort("domain-only draw returned invalid requested columns")
        requested <- M[, cols, drop = FALSE]
        if (anyNA(requested) || any(!vapply(requested, function(x) if (is.numeric(x)) all(is.finite(x)) else TRUE, 
            logical(1)))) new121_abort("domain-only requested values are not finite")
        if (!isTRUE(deterministic(function() validate_support(requested, cc, paste(context, label, paste0("requested_domain=", 
            g), sep = ";"))))) new121_abort("domain-only requested support failure")
        M
    })
    walk <- function(x) {
        if (!is.call(x)) 
            return(x)
        if (identical(x[[1L]], as.name("<-")) && identical(x[[2L]], as.name("draw")) && is.call(x[[3L]]) && identical(x[[3L]][[1L]], 
            as.name("function"))) {
            declarations <<- declarations + 1L
            return(as.call(list(as.name("{"), x, domain_definition)))
        }
        if (identical(x[[1L]], as.name("draw")) && length(x) >= 3L && identical(x[[3L]], as.name("ug"))) {
            replacements <<- replacements + 1L
            return(as.call(c(list(as.name("draw_domain"), x[[2L]], x[[3L]], as.name("g")), as.list(x)[-seq_len(3L)])))
        }
        as.call(lapply(as.list(x), walk))
    }
    body(f) <- walk(body(f))
    stopifnot(declarations == 1L, replacements == 2L)
    environment(f) <- methods
    methods$nmr_product_stream <- f
    invisible(methods)
}

nmr_ref_width105 <- function(original) {
    stopifnot(is.function(original), identical(formals(original)$uniform_dim, 46L))
    f <- original
    formals(f)$uniform_dim <- 105L
    walk <- function(x) {
        if (identical(x, 46L)) 
            return(105L)
        if (is.character(x)) 
            return(gsub("WIDTH46", "WIDTH105", x, fixed = TRUE))
        if (is.call(x)) 
            return(as.call(lapply(as.list(x), walk)))
        x
    }
    body(f) <- walk(body(f))
    f
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        nmr_install_product_runtime <- function(methods) {
            stopifnot(is.environment(methods), is.function(methods$omics_product_stream))
            f <- methods$omics_product_stream
            replacements <- 0L
            declarations <- 0L
            domain_definition <- quote(draw_domain <- function(b, u, g, label) {
                M <- deterministic(function() draw_joint(b = b, C = cc, uniforms = u, requested_domain = g))
                cols <- domains[[g]]
                if (!is.data.frame(M) || nrow(M) != nrow(cc) || !all(cols %in% names(M))) new121_abort("domain-only draw returned invalid requested columns")
                requested <- M[, cols, drop = FALSE]
                if (anyNA(requested) || any(!vapply(requested, function(x) if (is.numeric(x)) all(is.finite(x)) else TRUE, 
                  logical(1)))) new121_abort("domain-only requested values are not finite")
                if (!isTRUE(deterministic(function() validate_support(requested, cc, paste(context, label, paste0("requested_domain=", 
                  g), sep = ";"))))) new121_abort("domain-only requested support failure")
                M
            })
            walk <- function(x) {
                if (!is.call(x)) return(x)
                if (identical(x[[1L]], as.name("<-")) && identical(x[[2L]], as.name("draw")) && is.call(x[[3L]]) && 
                  identical(x[[3L]][[1L]], as.name("function"))) {
                  declarations <<- declarations + 1L
                  return(as.call(list(as.name("{"), x, domain_definition)))
                }
                if (identical(x[[1L]], as.name("draw")) && length(x) >= 3L && identical(x[[3L]], as.name("ug"))) {
                  replacements <<- replacements + 1L
                  return(as.call(c(list(as.name("draw_domain"), x[[2L]], x[[3L]], as.name("g")), as.list(x)[-seq_len(3L)])))
                }
                as.call(lapply(as.list(x), walk))
            }
            body(f) <- walk(body(f))
            stopifnot(declarations == 1L, replacements == 2L)
            environment(f) <- methods
            methods$nmr_product_stream <- f
            invisible(methods)
        }
        nmr_ref_width105 <- function(original) {
            stopifnot(is.function(original), identical(formals(original)$uniform_dim, 46L))
            f <- original
            formals(f)$uniform_dim <- 105L
            walk <- function(x) {
                if (identical(x, 46L)) return(105L)
                if (is.character(x)) return(gsub("WIDTH46", "WIDTH105", x, fixed = TRUE))
                if (is.call(x)) return(as.call(lapply(as.list(x), walk)))
                x
            }
            body(f) <- walk(body(f))
            f
        }
    }), envir = context)
}

