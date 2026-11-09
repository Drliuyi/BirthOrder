# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

keep <- c("transcriptome_theme", "transcriptome_program_palette")

pal <- c(CTRL = "#7651A6", CAD = "#B65F64", T2D = "#5683A8")

genes <- c("FCRL6", "GZMA", "GZMH", "KLRD1", "TOP1")

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        suppressPackageStartupMessages({
            library(ggplot2)
            library(data.table)
            library(jsonlite)
        })
        args <- commandArgs(FALSE)
        script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
        project <- Sys.getenv("PROJECT_ROOT")
        if (!nzchar(project) || !dir.exists(project)) stop("Set PROJECT_ROOT to the existing Birth Order project root")
        root <- Sys.getenv("OUTPUT_ROOT", file.path(project, "reports", "birth", "birthorder_expression_style_20261105_v1", 
            "bulk"))
        input <- Sys.getenv("INPUT_ROOT", file.path(project, "reports", "birth", "birthorder_bulk_GSE23561_YYstyle_clean_halfwidth_20261101_v3", 
            "data"))
        script_root <- Sys.getenv("SCRIPT_ROOT", file.path(project, "scripts"))
        reference <- Sys.getenv("YY_FIGURE_REFERENCE", file.path(script_root, "yy", "yystage", "R", "transcriptome_figure.R"))
        cat("INPUT_ROOT:", input, "\nSCRIPT_ROOT:", script_root, "\nOUTPUT_ROOT:", root, "\nYY_FIGURE_REFERENCE:", 
            reference, "\n")
        required <- c(file.path(input, "PLOT_SOURCE.csv"), file.path(input, "ALL_RAW_P_AND_STARS.csv"), file.path(input, 
            "DISPLAYED_STAR_BRACKETS.csv"), reference)
        if (!all(file.exists(required))) stop("Missing inputs: ", paste(required[!file.exists(required)], collapse = ", "))
        dir.create(root, recursive = TRUE, showWarnings = FALSE)
        for (child in c("figures", "metadata", "code")) dir.create(file.path(root, child), showWarnings = FALSE)
        before <- tools::md5sum(required)
        formats <- strsplit(Sys.getenv("EXPORT_FORMATS", "png,pdf,svg"), ",", fixed = TRUE)[[1]]
        stopifnot(all(formats %in% c("png", "pdf", "svg")))
        for (ext in formats) if (file.exists(file.path(root, "figures", paste0("GSE23561_YYstyle_clean.", ext)))) stop("Figure exists; use a new OUTPUT_ROOT")
        reference_ast <- parse(reference)
        ref_env <- new.env(parent = globalenv())
        keep <- c("transcriptome_theme", "transcriptome_program_palette")
        for (expr in reference_ast) if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) && 
            as.character(expr[[2]]) %in% keep) eval(expr, envir = ref_env)
        stopifnot(all(vapply(keep, exists, logical(1), envir = ref_env, inherits = FALSE)))
        pal <- c(CTRL = "#7651A6", CAD = "#B65F64", T2D = "#5683A8")
        d <- fread(required[1])
        ann <- fread(required[2])
        sig <- fread(required[3])
        genes <- c("FCRL6", "GZMA", "GZMH", "KLRD1", "TOP1")
        stopifnot(nrow(d) == 115, uniqueN(d$gsm) == 23, all(d$dataset == "GSE23561"), nrow(ann) == 15)
        d[, `:=`(group, factor(group, levels = c("CTRL", "CAD", "T2D")))]
        d[, `:=`(x_center, match(gene, genes))]
        d[, `:=`(x, x_center + c(CTRL = -0.24, CAD = 0, T2D = 0.24)[as.character(group)])]
        d[, `:=`(box_id, interaction(gene, group, drop = TRUE))]
        span <- diff(range(d$expression_z))
        stopifnot(nrow(sig) == 9, all(sig$stars %in% c("*", "**", "***")), all(sig$P == ann$P[match(paste(sig$gene, 
            sig$contrast), paste(ann$gene, ann$contrast))]))
        strip_top <- min(d$expression_z) - 0.03 * span
        strip_bottom <- strip_top - 0.075 * span
        upper <- max(c(max(d$expression_z), sig$bracket_y + 0.045 * span)) + 0.055 * span
        width <- (13.4 * 1.92/(1.92 + 0.68))/2
        height <- 12.52 * 0.9/(0.12 + 0.9 + 0.9 + 0.9 + 0.88)
        p <- ggplot(d, aes(x = x, y = expression_z, fill = group, colour = group, group = box_id)) + geom_boxplot(width = 0.088, 
            outlier.shape = NA, linewidth = 0.34, position = "identity") + geom_segment(data = sig, aes(x = x_left, 
            xend = x_right, y = bracket_y, yend = bracket_y), inherit.aes = FALSE, colour = "#111111", linewidth = 0.25) + 
            geom_segment(data = sig, aes(x = x_left, xend = x_left, y = bracket_y - 0.014 * span, yend = bracket_y), 
                inherit.aes = FALSE, colour = "#111111", linewidth = 0.25) + geom_segment(data = sig, aes(x = x_right, 
            xend = x_right, y = bracket_y - 0.014 * span, yend = bracket_y), inherit.aes = FALSE, colour = "#111111", 
            linewidth = 0.25) + geom_text(data = sig, aes(x = (x_left + x_right)/2, y = bracket_y + 0.028 * span, 
            label = stars), inherit.aes = FALSE, size = 4.4, fontface = "bold", colour = "#000000") + annotate("rect", 
            xmin = 0.5, xmax = 5.5, ymin = strip_bottom, ymax = strip_top, fill = "#4F7D70", colour = NA, alpha = 0.92) + 
            annotate("text", x = 3, y = (strip_bottom + strip_top)/2, label = "GSE23561", colour = "white", size = 3, 
                fontface = "bold") + scale_x_continuous(breaks = 1:5, labels = genes, limits = c(0.5, 5.5), expand = c(0, 
            0)) + scale_y_continuous(limits = c(strip_bottom, upper), expand = c(0, 0)) + scale_fill_manual(values = setNames(scales::alpha(pal, 
            c(0.22, 0.67, 0.67)), names(pal))) + scale_colour_manual(values = pal) + labs(x = NULL, y = "Gene expression (z score)", 
            fill = NULL, colour = NULL) + ref_env$transcriptome_theme(8.6) + theme(legend.position = "top", legend.direction = "horizontal", 
            legend.text = element_text(size = 9.5, colour = "black"), axis.title.y = element_text(size = 10.3, 
                colour = "black"), axis.text.x = element_text(angle = 58, hjust = 1, vjust = 1, size = 8, lineheight = 0.95, 
                colour = "black"), axis.text.y = element_text(size = 8.8, colour = "black"), plot.margin = margin(4, 
                4, 4, 4)) + guides(fill = guide_legend(nrow = 1, override.aes = list(fill = unname(pal))), colour = "none")
        if ("png" %in% formats) ggsave(file.path(root, "figures/GSE23561_YYstyle_clean.png"), p, width = width, 
            height = height, dpi = 600, device = ragg::agg_png, bg = "white")
        if ("pdf" %in% formats) ggsave(file.path(root, "figures/GSE23561_YYstyle_clean.pdf"), p, width = width, 
            height = height, device = "pdf", family = "Helvetica", useDingbats = FALSE, bg = "white")
        if ("svg" %in% formats) ggsave(file.path(root, "figures/GSE23561_YYstyle_clean.svg"), p, width = width, 
            height = height, device = svglite::svglite, bg = "white")
        stopifnot(identical(before, tools::md5sum(required)))
        fwrite(sig, file.path(root, "metadata/DISPLAYED_STAR_BRACKETS.csv"))
        write_json(list(reference = reference, theme = "original evaluated reference theme, enlarged typography only", 
            source_files = as.list(required), source_md5 = as.list(before), source_unchanged = TRUE, palette = as.list(pal), 
            box_width = 0.088, group_offsets = c(-0.24, 0, 0.24), frame_linewidth = 0.975, figure_width_inches = width, 
            figure_height_inches = height, horizontal_scale = 0.5, zero_reference_line = FALSE, reference_geometry = "reference bulk-panel width halved; height unchanged", 
            title = NULL, caption = NULL, subtitle = NULL, scatter_points = FALSE, bottom_strip = "GSE23561", datasets = "GSE23561", 
            coordinate_panels = nrow(ggplot_build(p)$layout$layout), preserved_values = 115, preserved_tests = 15, 
            star_comparisons = 9, model_refits = 0, shared_control_not_duplicated = TRUE, statistics_recomputed = 0, 
            typography = list(axis_title = 10.3, axis_tick = 8.8, gene = 8, legend = 9.5, star_mm = 4.4, dataset_mm = 3), 
            star_rule = "exact nine saved unadjusted-P annotations; no significance recalculation"), file.path(root, 
            "metadata/STYLE_AND_AUDIT.json"), pretty = TRUE, auto_unbox = TRUE)
        if (!file.copy(script_file, file.path(root, "code", "draw_bulk_style.R"), overwrite = FALSE)) stop("Cannot archive exact renderer")
        cat("CLEAN_YY_REFERENCE_PLOT_COMPLETE\n")
    }), envir = context)
}

