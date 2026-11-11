# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

core <- c("FCRL6", "GZMA", "GZMH", "KLRD1", "KIR2DL2", "KIR2DL3", "TOP1")

write_new <- function(x, f) {
    p <- file.path(out, f)
    if (file.exists(p)) 
        stop("Output already exists: ", p)
    fwrite(x, p, sep = if (endsWith(f, ".tsv")) 
        intToUtf8(9)
    else ",", na = "NA")
}

read_matrix <- function(a) {
    l <- readLines(gzfile(file.path(input, a, paste0(a, "_series_matrix.txt.gz"))))
    i <- grep("^!series_matrix_table_begin", l)
    j <- grep("^!series_matrix_table_end", l)
    d <- fread(text = paste(l[(i + 1):(j - 1)], collapse = "\n"))
    m <- as.matrix(d[, -1])
    rownames(m) <- as.character(d[[1]])
    storage.mode(m) <- "double"
    m
}

meta <- function(a) as.data.frame(fread(file.path(out, "metadata", paste0(a, "_SAMPLES.tsv"))))

normalize_symbols <- function(s) {
    s <- trimws(s)
    ans <- ifelse(s %in% official, s, NA_character_)
    u <- unique(s[is.na(ans) & nzchar(s) & !grepl("///|;|\\|", s)])
    if (length(u)) {
        valid <- intersect(u, keys(org.Hs.eg.db, keytype = "ALIAS"))
        if (length(valid)) {
            z <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, valid, c("SYMBOL"), "ALIAS"))
            sp <- split(z$SYMBOL, z$ALIAS)
            sp <- lapply(sp, unique)
            one <- sp[lengths(sp) == 1]
            mp <- setNames(unlist(one), names(one))
            ix <- which(is.na(ans))
            ans[ix] <- mp[s[ix]]
        }
    }
    ans
}

annotation <- function(a, gpl, m) {
    z <- as.data.frame(fread(file.path(out, "metadata", paste0(gpl, "_ANNOTATION.tsv")), colClasses = "character"))
    z <- z[match(rownames(m), z$probe), ]
    z$gene <- normalize_symbols(z$original_symbol)
    z$gene[z$control == "TRUE" & !is.na(z$control)] <- NA
    z$dataset <- a
    z$mapping_status <- ifelse(is.na(z$gene), "unmapped_or_ambiguous", ifelse(z$gene == z$original_symbol, "official_symbol", 
        "unique_alias"))
    z$mapping_status[z$gene %in% c("KIR2DL2", "KIR2DL3")] <- "KIR_homology_unresolved"
    z$representative <- FALSE
    good <- which(!is.na(z$gene) & z$mapping_status != "KIR_homology_unresolved" & rowSums(is.finite(m)) == ncol(m) & 
        apply(m, 1, var) > 0)
    ord <- if (gpl == "GPL10775") 
        order(as.numeric(z$probe[good]))
    else order(z$probe[good])
    g <- good[ord]
    g <- g[!duplicated(z$gene[g])]
    z$representative[g] <- TRUE
    list(map = z, y = m[g, , drop = FALSE], gene = z$gene[g], probe = z$probe[g])
}

fit <- function(y, design, contrasts, trend = FALSE) {
    stopifnot(qr(design)$rank == ncol(design), nrow(design) == ncol(y), nrow(design) > ncol(design))
    f <- lmFit(y, design)
    eBayes(contrasts.fit(f, contrasts), trend = trend, robust = TRUE)
}

tab <- function(f, genes, probes, a, model) {
    rbindlist(lapply(seq_len(ncol(f$coefficients)), function(k) {
        se <- f$stdev.unscaled[, k] * sqrt(f$s2.post)
        q <- qt(0.975, f$df.total)
        data.table(dataset = a, model = model, contrast = colnames(f$coefficients)[k], gene = genes, probe = probes, 
            effect = f$coefficients[, k], SE = se, CI_low = f$coefficients[, k] - q * se, CI_high = f$coefficients[, 
                k] + q * se, P = f$p.value[, k], FDR_genome = p.adjust(f$p.value[, k], "BH"), t = f$t[, k], AveExpr = f$Amean, 
            n_tested = nrow(f$coefficients))
    }))
}

coretab <- function(t, genes = core) {
    z <- rbindlist(lapply(unique(t$contrast), function(c) {
        merge(data.table(gene = genes, contrast = c), t[contrast == c], all.x = TRUE, by = c("gene", "contrast"), 
            sort = FALSE)
    }), fill = TRUE)
    z$dataset <- unique(t$dataset)
    z$model <- unique(t$model)
    z$measurement <- ifelse(is.na(z$effect), "not_reliably_covered_or_filtered", "annotation_resolved_no_detection_p")
    z$measurement[z$gene %in% c("KIR2DL2", "KIR2DL3")] <- "KIR_homology_unresolved"
    z[, `:=`(FDR_core7, p.adjust(P, "BH", n = 7)), by = contrast]
    z
}

qc <- function(y, s, a) {
    p <- prcomp(t(y), center = TRUE, scale. = FALSE)
    pct <- 100 * p$sdev^2/sum(p$sdev^2)
    z <- data.table(gsm = s$gsm, group = s$group, PC1 = p$x[, 1], PC2 = p$x[, 2])
    write_new(z, paste0("results/", a, "_PCA.csv"))
    pdf(file.path(out, "figures", paste0("QC_", a, ".pdf")), width = 11, height = 5)
    par(mfrow = c(1, 2), mar = c(6, 4, 3, 1))
    cc <- as.integer(factor(s$group))
    boxplot(y, col = cc + 1, outline = FALSE, names = FALSE, main = paste(a, "expression distribution"), ylab = "Analysis expression scale")
    plot(p$x[, 1:2], col = cc + 1, pch = 19, xlab = sprintf("PC1 (%.1f%%)", pct[1]), ylab = sprintf("PC2 (%.1f%%)", 
        pct[2]), main = paste(a, "PCA"))
    legend("topright", levels(factor(s$group)), pch = 19, col = seq_along(unique(s$group)) + 1, cex = 0.65)
    dev.off()
}

audit <- list()

maps <- list()

allcore <- list()

cache <- list()

A <- "GSE23561"

A <- "GSE90074"

A <- "GSE184050"

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        suppressPackageStartupMessages({
            library(limma)
            library(data.table)
            library(AnnotationDbi)
            library(org.Hs.eg.db)
        })
        args <- commandArgs(FALSE)
        sf <- sub("^--file=", "", args[grep("^--file=", args)][1])
        task <- normalizePath(file.path(dirname(sf), ".."))
        input <- Sys.getenv("INPUT_ROOT", file.path(task, "data"))
        out <- Sys.getenv("OUTPUT_ROOT", task)
        scriptroot <- Sys.getenv("SCRIPT_ROOT", dirname(sf))
        cat("INPUT_ROOT=", input, "\nOUTPUT_ROOT=", out, "\nSCRIPT_ROOT=", scriptroot, "\n")
        stopifnot(dir.exists(input), dir.exists(out))
        core <- c("FCRL6", "GZMA", "GZMH", "KLRD1", "KIR2DL2", "KIR2DL3", "TOP1")
        write_new <- function(x, f) {
            p <- file.path(out, f)
            if (file.exists(p)) stop("Output already exists: ", p)
            fwrite(x, p, sep = if (endsWith(f, ".tsv")) intToUtf8(9) else ",", na = "NA")
        }
        read_matrix <- function(a) {
            l <- readLines(gzfile(file.path(input, a, paste0(a, "_series_matrix.txt.gz"))))
            i <- grep("^!series_matrix_table_begin", l)
            j <- grep("^!series_matrix_table_end", l)
            d <- fread(text = paste(l[(i + 1):(j - 1)], collapse = "\n"))
            m <- as.matrix(d[, -1])
            rownames(m) <- as.character(d[[1]])
            storage.mode(m) <- "double"
            m
        }
        meta <- function(a) as.data.frame(fread(file.path(out, "metadata", paste0(a, "_SAMPLES.tsv"))))
        official <- keys(org.Hs.eg.db, keytype = "SYMBOL")
        normalize_symbols <- function(s) {
            s <- trimws(s)
            ans <- ifelse(s %in% official, s, NA_character_)
            u <- unique(s[is.na(ans) & nzchar(s) & !grepl("///|;|\\|", s)])
            if (length(u)) {
                valid <- intersect(u, keys(org.Hs.eg.db, keytype = "ALIAS"))
                if (length(valid)) {
                  z <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, valid, c("SYMBOL"), "ALIAS"))
                  sp <- split(z$SYMBOL, z$ALIAS)
                  sp <- lapply(sp, unique)
                  one <- sp[lengths(sp) == 1]
                  mp <- setNames(unlist(one), names(one))
                  ix <- which(is.na(ans))
                  ans[ix] <- mp[s[ix]]
                }
            }
            ans
        }
        annotation <- function(a, gpl, m) {
            z <- as.data.frame(fread(file.path(out, "metadata", paste0(gpl, "_ANNOTATION.tsv")), colClasses = "character"))
            z <- z[match(rownames(m), z$probe), ]
            z$gene <- normalize_symbols(z$original_symbol)
            z$gene[z$control == "TRUE" & !is.na(z$control)] <- NA
            z$dataset <- a
            z$mapping_status <- ifelse(is.na(z$gene), "unmapped_or_ambiguous", ifelse(z$gene == z$original_symbol, 
                "official_symbol", "unique_alias"))
            z$mapping_status[z$gene %in% c("KIR2DL2", "KIR2DL3")] <- "KIR_homology_unresolved"
            z$representative <- FALSE
            good <- which(!is.na(z$gene) & z$mapping_status != "KIR_homology_unresolved" & rowSums(is.finite(m)) == 
                ncol(m) & apply(m, 1, var) > 0)
            ord <- if (gpl == "GPL10775") order(as.numeric(z$probe[good])) else order(z$probe[good])
            g <- good[ord]
            g <- g[!duplicated(z$gene[g])]
            z$representative[g] <- TRUE
            list(map = z, y = m[g, , drop = FALSE], gene = z$gene[g], probe = z$probe[g])
        }
        fit <- function(y, design, contrasts, trend = FALSE) {
            stopifnot(qr(design)$rank == ncol(design), nrow(design) == ncol(y), nrow(design) > ncol(design))
            f <- lmFit(y, design)
            eBayes(contrasts.fit(f, contrasts), trend = trend, robust = TRUE)
        }
        tab <- function(f, genes, probes, a, model) {
            rbindlist(lapply(seq_len(ncol(f$coefficients)), function(k) {
                se <- f$stdev.unscaled[, k] * sqrt(f$s2.post)
                q <- qt(0.975, f$df.total)
                data.table(dataset = a, model = model, contrast = colnames(f$coefficients)[k], gene = genes, probe = probes, 
                  effect = f$coefficients[, k], SE = se, CI_low = f$coefficients[, k] - q * se, CI_high = f$coefficients[, 
                    k] + q * se, P = f$p.value[, k], FDR_genome = p.adjust(f$p.value[, k], "BH"), t = f$t[, k], 
                  AveExpr = f$Amean, n_tested = nrow(f$coefficients))
            }))
        }
        coretab <- function(t, genes = core) {
            z <- rbindlist(lapply(unique(t$contrast), function(c) {
                merge(data.table(gene = genes, contrast = c), t[contrast == c], all.x = TRUE, by = c("gene", "contrast"), 
                  sort = FALSE)
            }), fill = TRUE)
            z$dataset <- unique(t$dataset)
            z$model <- unique(t$model)
            z$measurement <- ifelse(is.na(z$effect), "not_reliably_covered_or_filtered", "annotation_resolved_no_detection_p")
            z$measurement[z$gene %in% c("KIR2DL2", "KIR2DL3")] <- "KIR_homology_unresolved"
            z[, `:=`(FDR_core7, p.adjust(P, "BH", n = 7)), by = contrast]
            z
        }
        qc <- function(y, s, a) {
            p <- prcomp(t(y), center = TRUE, scale. = FALSE)
            pct <- 100 * p$sdev^2/sum(p$sdev^2)
            z <- data.table(gsm = s$gsm, group = s$group, PC1 = p$x[, 1], PC2 = p$x[, 2])
            write_new(z, paste0("results/", a, "_PCA.csv"))
            pdf(file.path(out, "figures", paste0("QC_", a, ".pdf")), width = 11, height = 5)
            par(mfrow = c(1, 2), mar = c(6, 4, 3, 1))
            cc <- as.integer(factor(s$group))
            boxplot(y, col = cc + 1, outline = FALSE, names = FALSE, main = paste(a, "expression distribution"), 
                ylab = "Analysis expression scale")
            plot(p$x[, 1:2], col = cc + 1, pch = 19, xlab = sprintf("PC1 (%.1f%%)", pct[1]), ylab = sprintf("PC2 (%.1f%%)", 
                pct[2]), main = paste(a, "PCA"))
            legend("topright", levels(factor(s$group)), pch = 19, col = seq_along(unique(s$group)) + 1, cex = 0.65)
            dev.off()
        }
        audit <- list()
        maps <- list()
        allcore <- list()
        cache <- list()
        A <- "GSE23561"
        m <- read_matrix(A)
        s <- meta(A)
        stopifnot(identical(colnames(m), s$gsm), ncol(m) == 35, all(is.finite(m)), min(m) >= 0)
        s$group <- c(Control = "CTRL", `Coronary Artery Disease` = "CAD", `Type 2 Diabetes` = "T2D")[s$disease_state]
        keep <- !is.na(s$group)
        s <- s[keep, ]
        m <- m[, keep]
        s$group <- factor(s$group, levels = c("CTRL", "CAD", "T2D"))
        s$age <- as.numeric(s[["age_(yr)"]])
        s$sex <- factor(s$gender)
        stopifnot(identical(as.integer(table(s$group)), c(9L, 6L, 8L)))
        floorv <- min(m[m > 0])/2
        y <- log2(pmax(m, floorv))
        an <- annotation(A, "GPL10775", y)
        maps[[A]] <- an$map
        rownames(an$y) <- an$gene
        d <- model.matrix(~0 + group, s)
        colnames(d) <- levels(s$group)
        ct <- makeContrasts(CAD_vs_CTRL = CAD - CTRL, T2D_vs_CTRL = T2D - CTRL, CAD_vs_T2D = CAD - T2D, levels = d)
        f <- fit(an$y, d, ct, TRUE)
        t <- tab(f, an$gene, an$probe, A, "primary_group")
        write_new(t, "results/GSE23561_ALL_GENES.csv")
        c1 <- coretab(t)
        write_new(c1, "results/GSE23561_CORE7.csv")
        allcore[[A]] <- c1
        ix <- which(an$map$gene %in% core)
        fp <- fit(y[ix, , drop = FALSE], d, ct, TRUE)
        write_new(tab(fp, an$map$gene[ix], an$map$probe[ix], A, "all_core_probes_diagnostic_BH_within_probes"), 
            "results/GSE23561_CORE7_ALL_PROBES_DIAGNOSTIC.csv")
        complete <- is.finite(s$age) & !is.na(s$sex)
        da <- model.matrix(~0 + group + age + sex, s[complete, ])
        colnames(da) <- sub("^group", "", colnames(da))
        ca <- matrix(0, ncol(da), 3, dimnames = list(colnames(da), colnames(ct)))
        ca[rownames(ct), ] <- ct
        ta <- tab(fit(an$y[, complete, drop = FALSE], da, ca, TRUE), an$gene, an$probe, A, "age_sex_complete21_sensitivity")
        write_new(ta, "results/GSE23561_AGE_SEX_ALL_GENES.csv")
        qc(an$y, s, A)
        cache[[A]] <- list(y = an$y, s = s, results = t, adjusted = ta)
        audit[[A]] <- data.table(dataset = A, n_downloaded = 35, n_analyzed = 23, n_donors = 23, n_genes = nrow(an$y), 
            scale = "log2(sum10000 median; zeros floored)", zero_floor = floorv, n_input_features = nrow(m), n_model_columns = ncol(d), 
            n_residual_df = nrow(d) - ncol(d), note = "RA and metabolic syndrome excluded before transform/model; no second normalization; age+sex sensitivity")
        cat("PRIORITY1_COMPLETE\n")
        print(c1[gene %in% core[1:4], .(gene, contrast, effect, P, FDR_genome)])
        flush.console()
        saveRDS(cache[[A]], file.path(out, "results/GSE23561_analysis_cache.rds"))
        A <- "GSE90074"
        m <- read_matrix(A)
        s <- meta(A)
        stopifnot(identical(colnames(m), s$gsm), ncol(m) == 143, all(is.finite(m)))
        s$group <- factor(ifelse(s$obstructive_cad == "Y", "Obstructive", "Nonobstructive"), levels = c("Nonobstructive", 
            "Obstructive"))
        s$sex <- factor(s$gender)
        s$ancestry <- factor(s$ancestry)
        s$diabetes <- factor(s$diabetes)
        s$bmi <- as.numeric(s$bmi)
        s$batch <- factor(s$batchid)
        an <- annotation(A, "GPL6480", m)
        maps[[A]] <- an$map
        rownames(an$y) <- an$gene
        d <- model.matrix(~0 + group, s)
        colnames(d) <- levels(s$group)
        ct <- makeContrasts(Obstructive_vs_Nonobstructive = Obstructive - Nonobstructive, levels = d)
        t <- tab(fit(an$y, d, ct, TRUE), an$gene, an$probe, A, "primary_group")
        write_new(t, "results/GSE90074_ALL_GENES.csv")
        allcore[[A]] <- coretab(t)
        da <- model.matrix(~0 + group + sex + ancestry + diabetes + bmi + batch, s)
        colnames(da) <- sub("^group", "", colnames(da))
        ca <- matrix(0, ncol(da), 1, dimnames = list(colnames(da), colnames(ct)))
        ca[rownames(ct), ] <- ct
        ta <- tab(fit(an$y, da, ca, TRUE), an$gene, an$probe, A, "sex_ancestry_diabetes_BMI_batch_sensitivity")
        write_new(ta, "results/GSE90074_ADJUSTED_ALL_GENES.csv")
        qc(an$y, s, A)
        cache[[A]] <- list(y = an$y, s = s, results = t, adjusted = ta)
        audit[[A]] <- data.table(dataset = A, n_downloaded = 143, n_analyzed = 143, n_donors = 143, n_genes = nrow(an$y), 
            scale = "submitted normalized log2 test/reference", zero_floor = NA_real_, n_input_features = nrow(m), 
            n_model_columns = ncol(d), n_residual_df = nrow(d) - ncol(d), note = "No repeat log/normalization; age and medication absent; nonobstructive is not healthy")
        cat("PRIORITY2_COMPLETE\n")
        flush.console()
        A <- "GSE184050"
        s <- meta(A)
        m <- as.matrix(read.delim(gzfile(file.path(input, A, "GSE184050_CCHC_T2D_DEseq2_normalized_count.txt.gz")), 
            row.names = 1, check.names = FALSE))
        storage.mode(m) <- "double"
        s$expression_column <- ifelse(s$disease_state == "T2D case", sub("T2D case sample ", "case_sample_", s$title), 
            sub("Healthy control sample ", "control_sample_", s$title))
        stopifnot(setequal(colnames(m), s$expression_column))
        m <- m[, s$expression_column]
        colnames(m) <- s$gsm
        s$donor <- factor(s$individual)
        s$time <- ifelse(s$timepoint == "follow-up", 1, 0)
        s$case <- as.integer(s$disease_state == "T2D case")
        s$group <- ifelse(s$case == 1, "Transition", "Control")
        stopifnot(nlevels(s$donor) == 58, all(table(s$donor) == 2), all(table(s$donor, s$time) == 1), all(is.finite(m)), 
            min(m) >= 1)
        genes <- normalize_symbols(rownames(m))
        map <- data.frame(dataset = A, probe = rownames(m), original_symbol = rownames(m), gene = genes, mapping_status = ifelse(is.na(genes), 
            "unmapped_or_ambiguous", ifelse(genes %in% c("KIR2DL2", "KIR2DL3"), "KIR_homology_unresolved", "official_or_unique_alias")), 
            representative = FALSE)
        good <- !is.na(genes) & !(genes %in% c("KIR2DL2", "KIR2DL3")) & rowMeans(m >= 10) >= 0.2 & apply(m, 1, 
            var) > 0
        ix <- which(good)
        ix <- ix[order(rownames(m)[ix])]
        ix <- ix[!duplicated(genes[ix])]
        map$representative[ix] <- TRUE
        maps[[A]] <- map
        y <- log2(m[ix, , drop = FALSE])
        rownames(y) <- genes[ix]
        s$time_control <- s$time * (1 - s$case)
        s$time_transition <- s$time * s$case
        d <- model.matrix(~0 + donor + time_control + time_transition, s)
        ct <- matrix(0, ncol(d), 3, dimnames = list(colnames(d), c("Transition_change", "Control_change", "Transition_minus_Control_change")))
        ct["time_transition", c(1, 3)] <- 1
        ct["time_control", c(2, 3)] <- c(1, -1)
        t <- tab(fit(y, d, ct, TRUE), rownames(y), rownames(m)[ix], A, "donor_fixed_time_interaction")
        write_new(t, "results/GSE184050_ALL_GENES.csv")
        allcore[[A]] <- coretab(t)
        qc(y, s, A)
        cache[[A]] <- list(y = y, s = s, results = t, normalized_counts = m[ix, , drop = FALSE])
        audit[[A]] <- data.table(dataset = A, n_downloaded = 116, n_analyzed = 116, n_donors = 58, n_genes = nrow(y), 
            scale = "log2(DESeq2 normalized count clipped >=1 by submitter)", zero_floor = 1, n_input_features = nrow(m), 
            n_model_columns = ncol(d), n_residual_df = nrow(d) - ncol(d), note = "25 transition and 33 control donors; log2 normalized expression limma-trend; >=10 in >=20%; longitudinal difference-in-differences")
        write_new(rbindlist(maps, fill = TRUE), "metadata/GENE_MAPPING.tsv")
        write_new(rbindlist(audit, fill = TRUE), "metadata/DATASET_MANIFEST.tsv")
        cc <- rbindlist(allcore, fill = TRUE)
        write_new(cc, "results/CORE7_EXPRESSION_SUMMARY.csv")
        write_new(cc[gene %in% core[1:4]], "results/CORE4_EXPRESSION_SUMMARY.csv")
        saveRDS(cache, file.path(out, "results/analysis_cache.rds"))
        writeLines(capture.output(sessionInfo()), file.path(out, "logs/R_SESSIONINFO.txt"))
        cat("ALL_THREE_COMPLETE\n")
    }), envir = context)
}

