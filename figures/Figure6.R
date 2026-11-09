# Optional summary of existing expression estimates. Original UMAP/boxplot
# render methods are separate Figure6_* libraries; data and artwork are external.
draw_Figure6 <- function(workbook,out,config) {
  require_packages(c("ggplot2","svglite"))
  d <- read_source(workbook,"bulk")
  d <- d[d$gene %in% c("FCRL6","GZMA","GZMH","KLRD1"),,drop=FALSE]
  d <- plot_numbers(d,c("effect","CI_low","CI_high","P","FDR_genome"))
  p <- ggplot2::ggplot(d,ggplot2::aes(effect,gene,colour=contrast)) +
    ggplot2::geom_vline(xintercept=0,linewidth=.25,colour="grey65") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin=CI_low,xmax=CI_high),orientation="y",width=.15,
      position=ggplot2::position_dodge(.6),linewidth=.35) +
    ggplot2::geom_point(position=ggplot2::position_dodge(.6),size=1.7) +
    ggplot2::labs(x="Whole-blood expression contrast",y=NULL,colour=NULL) + paper_theme()
  save_panel(p,out,"Figure6_Bulk_summary",height=95)
}
