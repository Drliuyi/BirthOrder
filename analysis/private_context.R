# Path and lock values are supplied by the data custodian, never by GitHub.
birth_private_setting <- function(key) {
  stopifnot(is.character(key),length(key)==1L,grepl("^[A-Z0-9_]+$",key))
  value <- Sys.getenv(paste0("BIRTH_",key),unset="")
  if (!nzchar(value)) stop("Missing private setting: BIRTH_",key,call.=FALSE)
  value
}

# Each original script has its own environment: do not merge names across the
# single-screen, single-refinement or joint estimators.
birth_method_module <- function(code_root, filename, context=new.env(parent=globalenv())) {
  if (!is.environment(context) || !grepl("^[A-Za-z0-9_]+[.]R$",filename))
    stop("Explicit method filename and context required",call.=FALSE)
  context$birth_private_setting <- birth_private_setting
  sys.source(file.path(code_root,"analysis",filename),envir=context)
  context
}
