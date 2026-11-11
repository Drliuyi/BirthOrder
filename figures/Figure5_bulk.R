# Render whole-blood expression using externally supplied values and annotations.
figure5_transcriptome_theme <- function(base_size = 7.6) {
  base_size <- base_size * 1.15
  ggplot2::theme_classic(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold", colour = "#202020"),
      axis.text = ggplot2::element_text(colour = "#202020"),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "#202020"),
      axis.ticks.length = grid::unit(1.15, "mm"),
      panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.975, colour = "#202020"),
      strip.background = ggplot2::element_rect(fill = "#F1F1F1", colour = "#202020", linewidth = 0.45),
      strip.text = ggplot2::element_text(face = "bold", colour = "#202020"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.height = grid::unit(3.2, "mm"),
      legend.spacing.x = grid::unit(1.4, "mm"),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(size = ggplot2::rel(.86), hjust = 0, colour = "#4A4A4A"),
      plot.margin = ggplot2::margin(3.2, 3.2, 3.2, 3.2)
    )
}

draw_Figure5_bulk <- function(input_root, out, formats = c("png", "pdf", "svg")) {
  packages <- c("ggplot2", "data.table", "scales")
  if ("png" %in% formats) packages <- c(packages, "ragg")
  if ("svg" %in% formats) packages <- c(packages, "svglite")
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg)
  }
  if (!all(formats %in% c("png", "pdf", "svg")) || anyDuplicated(formats)) {
    stop("Export formats must be unique values from png, pdf and svg")
  }
  suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
  })
  required <- file.path(input_root, c("PLOT_SOURCE.csv", "ALL_RAW_P_AND_STARS.csv",
                                     "DISPLAYED_STAR_BRACKETS.csv"))
  if (!all(file.exists(required))) stop("Missing external plot inputs")
  if (dir.exists(out) || file.exists(out)) stop("Output directory must be new")
  before <- tools::md5sum(required)
  pal <- c(CTRL='#9B7BC6',CAD='#D96C80',T2D='#5B9FC8')
  d_all <- fread(required[1]); ann_all <- fread(required[2]); sig_all <- fread(required[3])
  stopifnot(nrow(d_all)==115,uniqueN(d_all$gsm)==23,nrow(ann_all)==15,nrow(sig_all)==9)
  genes <- c('FCRL6','GZMA','GZMH','KLRD1')
  d <- d_all[gene %chin% genes]; ann <- ann_all[gene %chin% genes]; sig <- sig_all[gene %chin% genes]
  stopifnot(nrow(d)==92,uniqueN(d$gsm)==23,all(d$dataset=='GSE23561'),nrow(ann)==12)
  d[,group:=factor(group,levels=c('CTRL','CAD','T2D'))]
  d[,x_center:=match(gene,genes)]
  d[,x:=x_center+c(CTRL=-.24,CAD=0,T2D=.24)[as.character(group)]]
  d[,box_id:=interaction(gene,group,drop=TRUE)]
  span <- diff(range(d_all$expression_z))
  stopifnot(nrow(sig)==7,all(sig$stars %in% c('*','**','***')),
            all(sig$P == ann$P[match(paste(sig$gene,sig$contrast),paste(ann$gene,ann$contrast))]))
  strip_top <- min(d_all$expression_z)-.030*span
  strip_bottom <- strip_top-.075*span
  upper <- max(c(max(d_all$expression_z),sig_all$bracket_y+.045*span))+.055*span
  width <- (13.4*1.92/(1.92+.68))/2
  height <- 12.52*.90/(.12+.90+.90+.90+.88)
  p <- ggplot(d,aes(x=x,y=expression_z,fill=group,colour=group,group=box_id))+
    geom_boxplot(width=.088,outlier.shape=NA,linewidth=.34,position='identity')+
    geom_segment(data=sig,aes(x=x_left,xend=x_right,y=bracket_y,yend=bracket_y),inherit.aes=FALSE,colour='#111111',linewidth=.25)+
    geom_segment(data=sig,aes(x=x_left,xend=x_left,y=bracket_y-.014*span,yend=bracket_y),inherit.aes=FALSE,colour='#111111',linewidth=.25)+
    geom_segment(data=sig,aes(x=x_right,xend=x_right,y=bracket_y-.014*span,yend=bracket_y),inherit.aes=FALSE,colour='#111111',linewidth=.25)+
    geom_text(data=sig,aes(x=(x_left+x_right)/2,y=bracket_y+.028*span,label=stars),inherit.aes=FALSE,size=4.4,fontface='bold',colour='#000000')+
    annotate('rect',xmin=.5,xmax=4.5,ymin=strip_bottom,ymax=strip_top,fill='#4F7D70',colour=NA,alpha=.92)+
    annotate('text',x=2.5,y=(strip_bottom+strip_top)/2,label='GSE23561',colour='white',size=3.0,fontface='bold')+
    scale_x_continuous(breaks=1:4,labels=genes,limits=c(.5,4.5),expand=c(0,0))+
    scale_y_continuous(limits=c(strip_bottom,upper),expand=c(0,0))+
    scale_fill_manual(values=setNames(scales::alpha(pal,.22),names(pal)))+
    scale_colour_manual(values=pal)+
    labs(x=NULL,y='Gene expression (z score)',fill=NULL,colour=NULL)+
    figure5_transcriptome_theme(8.6)+
    theme(legend.position='top',legend.direction='horizontal',
          legend.text=element_text(size=9.5,colour='black'),
          axis.title.y=element_text(size=10.3,colour='black'),
          axis.text.x=element_text(angle=58,hjust=1,vjust=1,size=8.0,lineheight=.95,colour='black'),
          axis.text.y=element_text(size=8.8,colour='black'),plot.margin=margin(4,4,4,4))+
    guides(fill=guide_legend(nrow=1,override.aes=list(fill=unname(scales::alpha(pal,.22)),colour=unname(pal))),colour='none')

  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  for (ext in formats) {
    path <- file.path(out, paste0("Figure5_bulk.", ext))
    if (ext == "png") ggsave(path, p, width = width, height = height, dpi = 600,
                            device = ragg::agg_png, bg = "white")
    if (ext == "pdf") ggsave(path, p, width = width, height = height,
                            device = "pdf", family = "Helvetica", useDingbats = FALSE, bg = "white")
    if (ext == "svg") ggsave(path, p, width = width, height = height,
                            device = svglite::svglite, bg = "white")
  }
  stopifnot(identical(before, tools::md5sum(required)))
  invisible(p)
}
