# Figure 3d: six paired disease contrasts from saved estimates.
draw_Figure3 <- function(workbook, out, config) {
  values <- read_source(workbook, "contrasts")
  directions <- read_source(workbook, "a_directions")
  required <- c("indicator_id","indicator_name_en","domain","mediator_order",
    "contrast","Delta_IE_per_1000","Delta_TE_per_1000","Pcontrast_percent")
  if (!all(required %in% names(values)) ||
      !all(c("indicator_id","sign","layer") %in% names(directions)))
    stop("Provide contrasts and a_directions sheets as documented by Figure3_contrasts.py")
  render_saved_python("Figure3_contrasts.py", list(values=values,directions=directions), out,
    function(p) c("--values",file.path(p,"values.csv"),
      "--a-directions",file.path(p,"directions.csv")),config)
}
