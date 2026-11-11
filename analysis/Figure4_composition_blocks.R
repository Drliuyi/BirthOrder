# Function library extracted without evaluating the original script.
# Source-time jobs, datasets, machine paths and sample-specific locks are not distributed.
# Supply the private context explicitly; this file does not launch an analysis.

s44_nmr_fields <- function() c("23406", "23418", "23435", "23587", "23588", "23591", "23594", "23595", "23596", 
    "23597", "23598", "23599", "23600", "23602", "23603", "23605", "23606", "23608", "23610", "23611", "23613", 
    "23617", "23618", "23622", "23623", "23628", "23639", "23640", "23648")

s44_nmr_blockspec <- function(source_spec) {
    stopifnot(identical(source_spec$schema, "NMR_COMPOSITION_BLOCKSPEC_V1"))
    composition <- function(id, parts, raw, fixed = TRUE, aux = list(), percent = TRUE) list(id = id, type = "composition", 
        parts = parts, raw = raw, aux = aux, raw_fields = names(raw), conditional_raw_fields = character(), total_scale = if (fixed) "fixed" else "variable", 
        fixed_total = if (fixed) 100 else NULL, component_unit = if (percent) "reported_percent" else "mmol/l", 
        strict_positive_parts_required = FALSE, zero_policy = "INHERITED_EXPLICIT_MASK_NO_PSEUDOCOUNT_NO_ROW_EXCLUSION")
    partial <- function(id, first, second, first_name, second_name) {
        parts <- setNames(list(paste0("F", first), paste0("F", second), paste0("100-F", first, "-F", second)), 
            c(first_name, second_name, "OTHER"))
        composition(id, parts, setNames(list(first_name, second_name), c(first, second)))
    }
    cholesterol <- function(id, C, CE, TG) composition(id, list(CE = paste0("F", CE), FC = paste0("F", C, "-F", 
        CE), TG = paste0("F", TG), OTHER = paste0("100-F", C, "-F", TG)), setNames(list("CE+FC", "CE", "TG"), c(C, 
        CE, TG)))
    scalar <- function(field) list(id = paste0("UNTIED_", field), type = "scalar", field = field, raw_fields = field, 
        conditional_raw_fields = character(), unchanged_conditional_family = TRUE)
    blocks <- list(composition("HDL_CHOLESTEROL", list(CE = "F23418", FC = "F23406-F23418"), list(`23406` = "TOTAL", 
        `23418` = "CE"), fixed = FALSE, percent = FALSE), partial("XL_VLDL_FC_TG_PERCENT", "23587", "23588", "FC", 
        "TG"), scalar("23591"), composition("M_VLDL_REPORTED_PERCENT", list(PL = "F23594", CE = "F23596", FC = "F23597", 
        TG = "F23598"), list(`23594` = "PL", `23595` = "CE+FC+E_C", `23596` = "CE", `23597` = "FC", `23598` = "TG"), 
        fixed = FALSE, aux = list(E_C = "F23595-F23596-F23597")), composition("S_VLDL_REPORTED_PERCENT", list(PL = "F23599", 
        CE = "F23600-F23602", FC = "F23602", TG = "F23603"), list(`23599` = "PL", `23600` = "CE+FC", `23602` = "FC", 
        `23603` = "TG"), fixed = FALSE), cholesterol("XS_VLDL_PARTIAL_PERCENT", "23605", "23606", "23608"), cholesterol("IDL_PARTIAL_PERCENT", 
        "23610", "23611", "23613"), partial("L_LDL_FC_TG_PERCENT", "23617", "23618", "FC", "TG"), partial("M_LDL_FC_TG_PERCENT", 
        "23622", "23623", "FC", "TG"), scalar("23628"), partial("M_HDL_PL_C_PERCENT", "23639", "23640", "PL", "C"), 
        scalar("23648"), scalar("23435"))
    covered <- unlist(lapply(blocks, nmrA_block_fields), use.names = FALSE)
    stopifnot(length(covered) == 29L, !anyDuplicated(covered), setequal(covered, s44_nmr_fields()), all(s44_nmr_fields() %in% 
        names(source_spec$raw_support)))
    for (block in blocks) if (identical(block$type, "composition")) {
        expr <- c(block$parts, block$raw, block$aux)
        symbols <- unique(unlist(lapply(expr, function(x) all.names(parse(text = x)[[1L]], functions = FALSE, unique = TRUE))))
        stopifnot(all(symbols %in% c(paste0("F", s44_nmr_fields()), names(block$parts), names(block$aux), "TOTAL")))
    }
    list(schema = "NMR_COMPOSITION_BLOCKSPEC_V1", version = "SELECTED44_RAW29_V1", information_members = 29L, coordinate_count_positive_interior = 29L, 
        expected_uniform_dim_with_base11 = 49L, reserved_zero_mask_uniform_slots = 9L, selected_raw_fields = s44_nmr_fields(), 
        blocks = blocks, raw_support = source_spec$raw_support[s44_nmr_fields()], sources = source_spec$sources, 
        raw_Y_columns_preserved = TRUE, all_recorded_discrepancies_preserved = TRUE, no_unselected_NMR_inputs_required = TRUE, 
        no_CORE39_completeness_required = TRUE, physical_definition = "Disjoint measured parts and aggregate remainders within each original lipoprotein denominator; never a whole29 simplex", 
        rounding_policy = "Retain measured M_VLDL cholesterol discrepancy and variable full-component reported-percent closure; never discard recorded variation", 
        partial_percent_policy = "Only selected disjoint fractions plus100-sum aggregate remainder; no unavailable component imputation", 
        zero_policy = "Inherited exact-zero conditional mask; no epsilon, clipping, row exclusion or retry")
}

run_in_private_context <- function(context) {
    stopifnot(is.environment(context))
    eval(quote({
        s44_nmr_fields <- function() c("23406", "23418", "23435", "23587", "23588", "23591", "23594", "23595", 
            "23596", "23597", "23598", "23599", "23600", "23602", "23603", "23605", "23606", "23608", "23610", 
            "23611", "23613", "23617", "23618", "23622", "23623", "23628", "23639", "23640", "23648")
        s44_nmr_blockspec <- function(source_spec) {
            stopifnot(identical(source_spec$schema, "NMR_COMPOSITION_BLOCKSPEC_V1"))
            composition <- function(id, parts, raw, fixed = TRUE, aux = list(), percent = TRUE) list(id = id, type = "composition", 
                parts = parts, raw = raw, aux = aux, raw_fields = names(raw), conditional_raw_fields = character(), 
                total_scale = if (fixed) "fixed" else "variable", fixed_total = if (fixed) 100 else NULL, component_unit = if (percent) "reported_percent" else "mmol/l", 
                strict_positive_parts_required = FALSE, zero_policy = "INHERITED_EXPLICIT_MASK_NO_PSEUDOCOUNT_NO_ROW_EXCLUSION")
            partial <- function(id, first, second, first_name, second_name) {
                parts <- setNames(list(paste0("F", first), paste0("F", second), paste0("100-F", first, "-F", second)), 
                  c(first_name, second_name, "OTHER"))
                composition(id, parts, setNames(list(first_name, second_name), c(first, second)))
            }
            cholesterol <- function(id, C, CE, TG) composition(id, list(CE = paste0("F", CE), FC = paste0("F", 
                C, "-F", CE), TG = paste0("F", TG), OTHER = paste0("100-F", C, "-F", TG)), setNames(list("CE+FC", 
                "CE", "TG"), c(C, CE, TG)))
            scalar <- function(field) list(id = paste0("UNTIED_", field), type = "scalar", field = field, raw_fields = field, 
                conditional_raw_fields = character(), unchanged_conditional_family = TRUE)
            blocks <- list(composition("HDL_CHOLESTEROL", list(CE = "F23418", FC = "F23406-F23418"), list(`23406` = "TOTAL", 
                `23418` = "CE"), fixed = FALSE, percent = FALSE), partial("XL_VLDL_FC_TG_PERCENT", "23587", "23588", 
                "FC", "TG"), scalar("23591"), composition("M_VLDL_REPORTED_PERCENT", list(PL = "F23594", CE = "F23596", 
                FC = "F23597", TG = "F23598"), list(`23594` = "PL", `23595` = "CE+FC+E_C", `23596` = "CE", `23597` = "FC", 
                `23598` = "TG"), fixed = FALSE, aux = list(E_C = "F23595-F23596-F23597")), composition("S_VLDL_REPORTED_PERCENT", 
                list(PL = "F23599", CE = "F23600-F23602", FC = "F23602", TG = "F23603"), list(`23599` = "PL", `23600` = "CE+FC", 
                  `23602` = "FC", `23603` = "TG"), fixed = FALSE), cholesterol("XS_VLDL_PARTIAL_PERCENT", "23605", 
                "23606", "23608"), cholesterol("IDL_PARTIAL_PERCENT", "23610", "23611", "23613"), partial("L_LDL_FC_TG_PERCENT", 
                "23617", "23618", "FC", "TG"), partial("M_LDL_FC_TG_PERCENT", "23622", "23623", "FC", "TG"), scalar("23628"), 
                partial("M_HDL_PL_C_PERCENT", "23639", "23640", "PL", "C"), scalar("23648"), scalar("23435"))
            covered <- unlist(lapply(blocks, nmrA_block_fields), use.names = FALSE)
            stopifnot(length(covered) == 29L, !anyDuplicated(covered), setequal(covered, s44_nmr_fields()), all(s44_nmr_fields() %in% 
                names(source_spec$raw_support)))
            for (block in blocks) if (identical(block$type, "composition")) {
                expr <- c(block$parts, block$raw, block$aux)
                symbols <- unique(unlist(lapply(expr, function(x) all.names(parse(text = x)[[1L]], functions = FALSE, 
                  unique = TRUE))))
                stopifnot(all(symbols %in% c(paste0("F", s44_nmr_fields()), names(block$parts), names(block$aux), 
                  "TOTAL")))
            }
            list(schema = "NMR_COMPOSITION_BLOCKSPEC_V1", version = "SELECTED44_RAW29_V1", information_members = 29L, 
                coordinate_count_positive_interior = 29L, expected_uniform_dim_with_base11 = 49L, reserved_zero_mask_uniform_slots = 9L, 
                selected_raw_fields = s44_nmr_fields(), blocks = blocks, raw_support = source_spec$raw_support[s44_nmr_fields()], 
                sources = source_spec$sources, raw_Y_columns_preserved = TRUE, all_recorded_discrepancies_preserved = TRUE, 
                no_unselected_NMR_inputs_required = TRUE, no_CORE39_completeness_required = TRUE, physical_definition = "Disjoint measured parts and aggregate remainders within each original lipoprotein denominator; never a whole29 simplex", 
                rounding_policy = "Retain measured M_VLDL cholesterol discrepancy and variable full-component reported-percent closure; never discard recorded variation", 
                partial_percent_policy = "Only selected disjoint fractions plus100-sum aggregate remainder; no unavailable component imputation", 
                zero_policy = "Inherited exact-zero conditional mask; no epsilon, clipping, row exclusion or retry")
        }
    }), envir = context)
}

