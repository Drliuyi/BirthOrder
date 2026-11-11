# Figure 4b: domain effects, proportion mediated and joint-effect intervals.
draw_Figure4 <- function(workbook, out, config) {
  domains <- read_source(workbook,"heatmap")
  forest <- read_source(workbook,"forest")
  models <- c("S44_BASE11","S44_NMR_PLUS40","S44_PROT_PLUS15")
  ids <- c("SMOKING","BODY_COMPOSITION","CLINICAL_LIPIDS","GLYCAEMIA",
    "NMR_ALL","CANDIDATE_PROTEINS")
  domains <- domains[domains$model_id %in% models & domains$domain_id %in% ids,,drop=FALSE]
  # Accept either one saved row per cell or the wide source-table schema.
  if (all(c("outcome","point_per_1000","q_BH") %in% names(domains))) {
    heatmap <- domains
  } else {
    heatmap <- do.call(rbind,lapply(c("IHD","T2D","Delta"),function(effect) {
      cols <- paste0(effect,c("_component_point_per_1000","_q_BH"))
      if (!all(cols %in% names(domains))) stop("Missing saved domain estimate/q columns")
      data.frame(model_id=domains$model_id,domain_id=domains$domain_id,
        outcome=effect,point_per_1000=domains[[cols[1]]],q_BH=domains[[cols[2]]])
    }))
  }
  forest <- forest[forest$model_id %in% models,,drop=FALSE]
  if (all(c("point_per_1000","CI95_low_per_1000","CI95_high_per_1000",
      "mediated_proportion_percent") %in% names(forest))) {
    joint <- forest
  } else {
    cols <- c("point_B0_per_1000","CI_type7_low_per_1000","CI_type7_high_per_1000","PM")
    if (!all(cols %in% names(forest))) stop("Missing saved joint-effect/proportion columns")
    joint <- data.frame(model_id=forest$model_id,point_per_1000=forest[[cols[1]]],
      CI95_low_per_1000=forest[[cols[2]]],CI95_high_per_1000=forest[[cols[3]]],
      mediated_proportion_percent=forest[[cols[4]]])
  }
  render_saved_python("Figure4_heatmap.py",list(heatmap=heatmap,joint=joint),out,
    function(p)c("--heatmap",file.path(p,"heatmap.csv"),"--joint",file.path(p,"joint.csv")),config)
}
