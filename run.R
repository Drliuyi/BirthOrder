#!/usr/bin/env Rscript
# No action is implicit. 'check' and 'methods' never open a research input.
arg_file <- grep("^--file=", commandArgs(), value=TRUE)[1]
root <- dirname(normalizePath(sub("^--file=", "", arg_file), mustWork=TRUE))
source(file.path(root, "config.R"))
for (name in c("common.R", "data.R", "association.R", "mediation.R", "joint.R", "Figure2.R", "SupplementaryFigure2.R", "Table1.R"))
  sys.source(file.path(root, "analysis", name), envir=.GlobalEnv)
args <- commandArgs(TRUE)
if (!length(args) || args[1] %in% c("-h", "--help")) {
  cat("Birth Order | analysis code only\n",
      "  bash run.sh check\n",
      "  bash run.sh methods Figure4\n",
      "  bash run.sh analyse Figure2 NEW_OUTPUT\n",
      "  bash run.sh analyse SupplementaryFigure2 NEW_OUTPUT\n",
      "  bash run.sh analyse Table1 NEW_OUTPUT\n",
      "  bash run.sh draw Figure4 PRIVATE_PLOT_WORKBOOK NEW_OUTPUT\n",
      "Mediation and expression stages are supplied as function libraries; see README.\n",
      sep="")
  quit(status=0)
}
config <- birth_config()
config$code_root <- root
action <- args[1]
if (action == "check") {
  files <- list.files(root, pattern="[.]R$", recursive=TRUE, full.names=TRUE)
  for (f in files) invisible(parse(f))
  cat("Parsed",length(files),"R source files; no private inputs read or models run.\n")
} else if (action == "methods") {
  assert(length(args)==2L, "Specify Figure2-5, SupplementaryFigure2, Table1 or Supplementary")
  require_packages("jsonlite")
  modules <- jsonlite::fromJSON(file.path(root,"reference/modules.json"), simplifyVector=FALSE)
  selected <- Filter(function(x) identical(x$scope,args[2]), modules)
  for (x in selected) cat(x$file,"|",length(x$functions),"functions\n")
  assert(length(selected)>0L,"No extracted method module for this target; see README")
} else if (action == "analyse") {
  assert(length(args)==3L && args[2] %in% c("Figure2","SupplementaryFigure2","Table1"),
    "Input analyses: Figure2, SupplementaryFigure2 or Table1. Other stages use the method libraries.")
  require_packages(c("survival","digest"))
  out <- new_output(args[3],config,root)
  result <- get(paste0("analyse_",args[2]),mode="function")(config,root)
  for (name in names(result)) write_table(result[[name]],out,paste0(name,".csv"))
} else if (action == "draw") {
  assert(length(args)==4L && args[2] %in% c(paste0("Figure",1:5),"SupplementaryFigure2"),
    "Specify Figure1-5 or SupplementaryFigure2, input workbook and new output")
  sys.source(file.path(root,"figures/style.R"),envir=.GlobalEnv)
  sys.source(file.path(root,"figures",paste0(args[2],".R")),envir=.GlobalEnv)
  out <- new_output(args[4],config,root)
  get(paste0("draw_",args[2]),mode="function")(args[3],out,config)
} else stop("Unknown action. No model or Bootstrap started.",call.=FALSE)
