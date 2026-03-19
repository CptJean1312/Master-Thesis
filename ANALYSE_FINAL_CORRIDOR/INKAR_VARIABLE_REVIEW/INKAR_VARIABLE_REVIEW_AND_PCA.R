#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
})

options(stringsAsFactors = FALSE)

paths <- list(
  inventory = "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/tables/corridor_inkar_indicator_inventory_dimensioned.csv",
  full_corridor_csv = "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/tables/corridor_full_inkar_original_codes.csv",
  output_dir = "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs"
)

dir.create(file.path(paths$output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(paths$output_dir, "logs", "inkar_variable_review_log.txt")
log_lines <- character()

log_message <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  log_lines <<- c(log_lines, line)
  message(line)
}

impute_median <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  x[is.na(x)] <- stats::median(x, na.rm = TRUE)
  x
}

rescale_01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / diff(rng)
}

prepare_matrix <- function(df, vars) {
  out <- df %>%
    select(all_of(vars)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

  keep_rows <- rowSums(!is.na(out)) > 0
  out <- out[keep_rows, , drop = FALSE]
  out <- out %>% mutate(across(everything(), impute_median))
  out
}

top_loading_table <- function(pca, vars, n_pc = 8, top_n = 8) {
  loadings <- as.data.frame(pca$rotation)
  loadings$variable <- rownames(loadings)
  pcs <- paste0("PC", seq_len(min(n_pc, ncol(loadings) - 1)))

  bind_rows(lapply(pcs, function(pc) {
    loadings %>%
      transmute(
        PC = pc,
        variable = variable,
        loading = .data[[pc]],
        abs_loading = abs(.data[[pc]])
      ) %>%
      arrange(desc(abs_loading)) %>%
      slice_head(n = top_n)
  }))
}

top_corr_pairs <- function(cor_mat, top_n = 40) {
  idx <- upper.tri(cor_mat, diag = FALSE)
  pairs <- data.frame(
    var1 = rownames(cor_mat)[row(cor_mat)[idx]],
    var2 = colnames(cor_mat)[col(cor_mat)[idx]],
    correlation = cor_mat[idx],
    abs_correlation = abs(cor_mat[idx]),
    stringsAsFactors = FALSE
  )
  pairs %>%
    arrange(desc(abs_correlation), var1, var2) %>%
    slice_head(n = top_n)
}

plot_corr_heatmap <- function(cor_mat, title, out_path, base_size = 7, hide_labels = FALSE) {
  hc <- hclust(as.dist(1 - cor_mat), method = "average")
  ord <- hc$order
  cor_ord <- cor_mat[ord, ord]

  long <- as.data.frame(as.table(cor_ord), stringsAsFactors = FALSE)
  names(long) <- c("var1", "var2", "correlation")
  long$var1 <- factor(long$var1, levels = rownames(cor_ord))
  long$var2 <- factor(long$var2, levels = colnames(cor_ord))

  p <- ggplot(long, aes(var1, var2, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "#8e1b1b", mid = "white", high = "#0b5f8a", midpoint = 0, limits = c(-1, 1)) +
    coord_equal() +
    labs(title = title, x = NULL, y = NULL, fill = "r") +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = if (hide_labels) element_blank() else element_text(angle = 90, vjust = 0.5, hjust = 1, size = base_size * 0.8),
      axis.text.y = if (hide_labels) element_blank() else element_text(size = base_size * 0.8)
    )

  ggsave(out_path, plot = p, width = 12, height = 11, dpi = 300)
}

run_pca_bundle <- function(df, vars, prefix, title) {
  data_imputed <- prepare_matrix(df, vars)
  scaled <- scale(data_imputed)
  pca <- prcomp(scaled)

  eigenvalues <- pca$sdev^2
  variance <- eigenvalues / sum(eigenvalues)
  scree <- tibble(
    PC = seq_along(eigenvalues),
    eigenvalue = eigenvalues,
    variance = variance,
    cumulative = cumsum(variance)
  )

  scores <- as.data.frame(pca$x)
  names(scores) <- paste0("PC", seq_len(ncol(scores)))

  loadings_top <- top_loading_table(pca, vars)
  cor_mat <- cor(data_imputed)
  corr_pairs <- top_corr_pairs(cor_mat)

  write_csv(scree, file.path(paths$output_dir, "tables", paste0(prefix, "_scree.csv")))
  write_csv(loadings_top, file.path(paths$output_dir, "tables", paste0(prefix, "_top_loadings.csv")))
  write_csv(as.data.frame(cor_mat) %>% tibble::rownames_to_column("variable"), file.path(paths$output_dir, "tables", paste0(prefix, "_correlation_matrix.csv")))
  write_csv(corr_pairs, file.path(paths$output_dir, "tables", paste0(prefix, "_top_correlation_pairs.csv")))

  plot_corr_heatmap(
    cor_mat = cor_mat,
    title = paste0(title, " correlation matrix"),
    out_path = file.path(paths$output_dir, "plots", paste0(prefix, "_correlation_heatmap.png")),
    base_size = if (length(vars) > 60) 3.5 else 6,
    hide_labels = length(vars) > 120
  )

  scree_plot <- ggplot(scree, aes(PC, variance)) +
    geom_line(color = "#0b5f8a", linewidth = 0.6) +
    geom_point(color = "#0b5f8a", size = 1.5) +
    geom_hline(yintercept = 1 / ncol(data_imputed), linetype = "dashed", color = "#8e1b1b") +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(
      title = paste0(title, " scree plot"),
      x = "Principal component",
      y = "Explained variance share"
    ) +
    theme_minimal(base_size = 11)

  ggsave(
    file.path(paths$output_dir, "plots", paste0(prefix, "_scree.png")),
    plot = scree_plot,
    width = 8.5,
    height = 5.5,
    dpi = 300
  )

  list(
    data_imputed = data_imputed,
    pca = pca,
    scree = scree,
    loadings_top = loadings_top,
    cor_mat = cor_mat,
    corr_pairs = corr_pairs,
    scores = scores
  )
}

log_message("Loading corridor INKAR inventory and full table.")
inventory <- read_csv(paths$inventory, show_col_types = FALSE)
full_corridor <- read_csv(paths$full_corridor_csv, show_col_types = FALSE)

all_codes <- inventory$Kuerzel

old51_codes <- c(
  "a_ALGII_SGBII", "a_BG1P", "a_BG5um", "a_BGKind", "a_Unterkunft_SGBII",
  "a_aloLang", "q_alo_u25_einw", "q_alo_ü55_einw",
  "q_einkst_bev", "q_kaufkraft", "q_gewst_bev", "d_steuereinnahme",
  "a_bev0003", "a_bev0306", "a_bev0618", "a_bev1825", "a_bev2530", "a_bev3050",
  "a_bev5065", "a_bev6575", "a_bev65um", "a_bev75um",
  "i_saldo_nat", "i_wans",
  "q_abhg_alt", "q_abhg_jung", "r_ewf_jungalt",
  "q_HH1", "a_hh_kind",
  "a_hheink_hoch", "a_hheink_mittel", "a_hheink_niedrig",
  "q_stud", "q_stud_1825", "q_stud_fh",
  "q_allgemeinärzte_bev", "q_hausarzt_bev", "q_internist_bev", "q_kinderarzt_kinder", "q_ärzte_bev",
  "a_bb_1000Mbits", "a_bb_100Mbits", "a_bb_50Mbits", "a_bb_4G",
  "m_G02_SUP_DIST", "m_Q01_APO_DIST", "m_Q07_HA_DIST", "m_OEV20_DIST", "m_P01_PRIM_DIST",
  "q_bev_fl", "q_bevsva_qkm"
)

curated_codes <- c(
  "a_ALGII_SGBII",
  "a_aloLang",
  "q_alo_u25_einw",
  "a_Minijobs",
  "q_kaufkraft",
  "q_einkst_bev",
  "a_hheink_niedrig",
  "a_bev65um",
  "q_abhg_alt",
  "q_HH1",
  "a_hh_kind",
  "m_G02_SUP_DIST",
  "m_Q01_APO_DIST",
  "m_Q07_HA_DIST",
  "m_OEV20_DIST",
  "q_ärzte_bev",
  "a_bb_100Mbits"
)

if (length(setdiff(old51_codes, names(full_corridor))) > 0) {
  stop("Some original wide-PCA variables are missing in the full corridor table.")
}
if (length(setdiff(curated_codes, names(full_corridor))) > 0) {
  stop("Some curated variables are missing in the full corridor table.")
}

subgroup_redundant <- c(
  "a_aloLang_f", "a_aloLang_m", "a_alo_ausländer", "a_alo_ausländer_f", "a_alo_ausländer_m",
  "a_alo_f", "a_alo_m", "a_alo_u25", "a_alo_u25_f", "a_alo_u25_m", "a_alo_ü55", "a_alo_ü55_f", "a_alo_ü55_m",
  "a_gb_Frauen", "a_gb_Männer", "a_bev6575_f", "a_bev65um_f", "a_bev75um_f", "a_bev1825_f", "a_bev2530_f",
  "a_ewfBG_55um_", "a_ewfBG_f", "a_ewfBG_u25_"
)

age_composition <- c(
  "a_bev0003", "a_bev0306", "a_bev0618", "a_bev1825", "a_bev2530", "a_bev3050", "a_bev5065", "a_bev6575",
  "a_bev75um", "q_abhg_jung", "r_ewf_jungalt"
)

income_complements <- c("a_hheink_hoch", "a_hheink_mittel", "a_Unterkunft_SGBII", "q_gewst_bev", "d_steuereinnahme")
digital_duplicates <- c("a_bb_1000Mbits", "a_bb_50Mbits", "a_bb_4G")
health_duplicates <- c("q_allgemeinärzte_bev", "q_hausarzt_bev", "q_internist_bev", "q_kinderarzt_kinder")
education_context <- c("q_stud", "q_stud_1825", "q_stud_fh")
education_low_coverage <- c("a_stud_1", "a_stud_a", "a_stud_m", "a_stud_w")
size_controls <- c("TN23-kataster_qkm", "alo", "xbev", "xbevf", "xbevm", "bev_korr", "q_HH", "ewf_1565_ges", "q_bev_fl", "q_bevsva_qkm", "sva", "svw")
access_alternatives <- c("a_G02_SUP_ANT", "a_OEV20_ANT", "a_P01_PRIM_ANT", "a_Q01_APO_ANT", "a_Q07_HA_ANT", "m_P01_PRIM_DIST")
sgbii_review <- c("a_BG1P", "a_BG5um", "a_BGKind", "a_ewfBG", "a_ewfBG_allein")
tourism_context <- c("m_übern", "q_schlafg_bev", "q_übern_bev")
public_finance <- c("q_sach", "q_investZ", "q_schlüsselzuw")
landuse_context <- c("a_freifläche", "a_landwirtschaft", "a_naturnah", "a_wald", "a_wasser", "a_suv_fl", "a_wo_r12", "a_wo_r5um", "a_wg_wo12", "a_wg_wo3um", "a_wo_wg12", "a_wo_wg3um", "a_fert_wg12", "a_fert_wo12", "a_fert_wohn", "a_gen_wo12", "a_gen_wo3um", "a_gest_bev")
dynamics_context <- c("i_saldo_nat", "i_wans", "a_geb_bev", "a_gest_bev")

review_table <- inventory %>%
  mutate(
    in_old51 = Kuerzel %in% old51_codes,
    in_curated = Kuerzel %in% curated_codes,
    proposed_status = case_when(
      in_curated ~ "core_curated_pca",
      Kuerzel == "a_Unterkunft_SGBII" ~ "drop_overlap_with_algii",
      Kuerzel == "a_BG1P" ~ "drop_definition_mismatch",
      Kuerzel %in% education_context ~ "review_education_context",
      Kuerzel %in% education_low_coverage ~ "drop_low_coverage",
      coverage_pct < 90 ~ "drop_low_coverage",
      Kuerzel %in% subgroup_redundant ~ "drop_redundant_subgroup",
      Kuerzel %in% age_composition ~ "drop_nested_age_block",
      Kuerzel %in% income_complements ~ "drop_complement_or_context",
      Kuerzel %in% digital_duplicates ~ "drop_duplicate_digital_measure",
      Kuerzel %in% health_duplicates ~ "drop_duplicate_health_measure",
      Kuerzel %in% size_controls ~ "context_or_control_only",
      Kuerzel %in% access_alternatives ~ "review_redundant_access_measure",
      Kuerzel %in% sgbii_review ~ "review_sgbii_specific_household_measure",
      Kuerzel %in% tourism_context ~ "drop_tourism_context",
      Kuerzel %in% public_finance ~ "drop_public_finance_low_coverage",
      Kuerzel %in% landuse_context ~ "context_land_use",
      Kuerzel %in% dynamics_context ~ "context_demographic_dynamics",
      TRUE ~ "review_not_in_curated_set"
    ),
    handling_note = case_when(
      in_curated ~ "Keep as core candidate in the curated vulnerability PCA.",
      Kuerzel == "a_Unterkunft_SGBII" ~ "Conceptually overlaps with a_ALGII_SGBII and is difficult to interpret as an independent vulnerability signal.",
      Kuerzel == "a_BG1P" ~ "Officially this is Einpersonen-Bedarfsgemeinschaften, not single-parent households; the old label was misleading.",
      Kuerzel %in% education_context ~ "Student concentration may be analytically interesting, but it mainly captures university-center structure rather than core social vulnerability.",
      Kuerzel %in% education_low_coverage | coverage_pct < 90 ~ "Too sparse in the corridor sample for a stable main PCA.",
      Kuerzel %in% subgroup_redundant ~ "Specific subgroup breakdown of a broader indicator; likely to duplicate the same latent dimension.",
      Kuerzel %in% age_composition ~ "Nested age-composition block; retaining many of these together overweights the same demographic structure.",
      Kuerzel %in% income_complements ~ "Either compositional counterpart or broader fiscal context rather than a distinct vulnerability mechanism.",
      Kuerzel %in% digital_duplicates ~ "Measures the same digital access construct at multiple thresholds; keep only one representative.",
      Kuerzel %in% health_duplicates ~ "Specialist-specific supply measure; keep total physician availability instead.",
      Kuerzel %in% size_controls ~ "Useful as descriptive context or control, but not as part of the main vulnerability index.",
      Kuerzel %in% access_alternatives ~ "Access concept already captured by distance-based measures; use as optional robustness check only.",
      Kuerzel %in% sgbii_review ~ "Potentially informative, but highly specific to SGB II household composition and hard to generalize.",
      Kuerzel %in% tourism_context ~ "Not part of the social-vulnerability concept for the main index.",
      Kuerzel %in% public_finance ~ "Outside the core vulnerability concept and additionally incomplete across municipalities.",
      Kuerzel %in% landuse_context ~ "Land-use and housing context variable, not a direct social-vulnerability indicator.",
      Kuerzel %in% dynamics_context ~ "Demographic dynamics may be reported descriptively, but are not core vulnerability dimensions here.",
      TRUE ~ "Review individually; not obviously required for the core vulnerability index."
    )
  ) %>%
  arrange(broad_dimension, desc(in_curated), desc(in_old51), Kuerzel)

write_csv(review_table, file.path(paths$output_dir, "tables", "corridor_variable_review_table.csv"))
write_csv(review_table %>% filter(in_old51), file.path(paths$output_dir, "tables", "original_wide51_variable_list.csv"))
write_csv(review_table %>% filter(in_curated), file.path(paths$output_dir, "tables", "curated_variable_list.csv"))

log_message("Running original 51-variable PCA review bundle.")
old51_bundle <- run_pca_bundle(full_corridor, old51_codes, "old51", "Original 51-variable set")

log_message("Running all-176-variable PCA review bundle.")
all176_bundle <- run_pca_bundle(full_corridor, all_codes, "all176", "All 176 corridor indicators")

log_message("Running curated PCA review bundle.")
curated_bundle <- run_pca_bundle(full_corridor, curated_codes, "curated", "Curated vulnerability set")

comparison <- tibble(
  model = c("original_51", "all_176", "curated"),
  variables = c(length(old51_codes), length(all_codes), length(curated_codes)),
  municipalities_used = c(nrow(old51_bundle$data_imputed), nrow(all176_bundle$data_imputed), nrow(curated_bundle$data_imputed)),
  pc1_variance = c(old51_bundle$scree$variance[1], all176_bundle$scree$variance[1], curated_bundle$scree$variance[1]),
  pc1_to_pc4_cumulative = c(old51_bundle$scree$cumulative[4], all176_bundle$scree$cumulative[4], curated_bundle$scree$cumulative[4]),
  pc1_to_pc8_cumulative = c(old51_bundle$scree$cumulative[8], all176_bundle$scree$cumulative[8], curated_bundle$scree$cumulative[min(8, nrow(curated_bundle$scree))])
)
write_csv(comparison, file.path(paths$output_dir, "tables", "pca_model_comparison.csv"))

pair_check <- read_csv(
  "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/tables/corridor_sgb2_pair_check.csv",
  show_col_types = FALSE
)

selected_review_codes <- c(
  curated_codes,
  "a_ALGII_SGBII", "a_Unterkunft_SGBII", "a_BG1P", "a_BGKind", "a_BG5um", "a_ewfBG_allein",
  "a_bev65um", "a_bev75um", "q_abhg_alt", "q_HH1", "a_hh_kind", "a_hheink_niedrig",
  "a_bb_100Mbits", "a_bb_50Mbits", "a_bb_1000Mbits", "q_ärzte_bev", "q_hausarzt_bev",
  "m_G02_SUP_DIST", "m_Q01_APO_DIST", "m_Q07_HA_DIST", "m_OEV20_DIST", "q_bev_fl", "q_bevsva_qkm"
)

selected_review <- review_table %>%
  filter(Kuerzel %in% selected_review_codes) %>%
  arrange(match(Kuerzel, selected_review_codes))

format_pct <- function(x) sprintf("%.2f%%", x * 100)

comparison_lines <- apply(comparison, 1, function(row) {
  sprintf(
    "- `%s`: %s variables, %s municipalities, PC1 variance `%s`, cumulative PC1-PC4 `%s`, cumulative PC1-PC8 `%s`",
    row[["model"]],
    row[["variables"]],
    row[["municipalities_used"]],
    format_pct(as.numeric(row[["pc1_variance"]])),
    format_pct(as.numeric(row[["pc1_to_pc4_cumulative"]])),
    format_pct(as.numeric(row[["pc1_to_pc8_cumulative"]]))
  )
})

dimension_lines <- lapply(split(review_table, review_table$broad_dimension), function(df) {
  c(
    sprintf("### %s", unique(df$broad_dimension)),
    "",
    sprintf("- Indicators in corridor inventory: `%s`", nrow(df)),
    sprintf("- Median coverage in this block: `%s%%`", format(round(median(df$coverage_pct, na.rm = TRUE), 2), nsmall = 2)),
    "",
    vapply(seq_len(nrow(df)), function(i) {
      sprintf(
        "- `%s` — %s | coverage `%s%%` | status `%s` | %s",
        df$Kuerzel[i],
        df$Indikator[i],
        format(df$coverage_pct[i], nsmall = 2),
        df$proposed_status[i],
        df$handling_note[i]
      )
    }, character(1)),
    ""
  )
})

report_lines <- c(
  "# INKAR Corridor Variable Review and PCA Note",
  "",
  sprintf("Created: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 1. Why this note exists",
  "",
  "This note does two things at once:",
  "- it prepares a concrete answer to the supervisor email about potentially counterintuitive INKAR variables;",
  "- it documents a cleaner decision path from the old exploratory wide PCA towards a curated, literature-guided vulnerability PCA for the EFAS corridor sample.",
  "",
  "## 2. Main conclusion",
  "",
  "- A wide PCA with many variables is still useful as an exploratory diagnostic tool.",
  "- It should not automatically be the final vulnerability index if the variable set contains compositional counterparts, nested age blocks, duplicated access measures, or variables with counterintuitive definitions.",
  "- The original 51-variable set was a broad exploratory screening set. It was not yet fully cleaned with respect to definition conflicts and duplicated constructs.",
  "- For the final thesis, the stronger strategy is: keep the old wide PCA as an exploratory/sensitivity result, but build the main interpretable index from a curated set of variables with clear theoretical roles.",
  "",
  "## 3. Direct answer to the email concern",
  "",
  sprintf("- In the current corridor raw-INKAR build, `a_ALGII_SGBII` and `a_Unterkunft_SGBII` have `%s` complete municipality pairs and a correlation of `%.3f`.", pair_check$n_complete_pairs[1], pair_check$correlation[1]),
  "- That is a strong negative relationship, but it is not a perfect `-1.0` complement relationship in the current corridor data.",
  "- So the safe interpretation is: they are conceptually overlapping and partly compositional, but they should not be treated as a clean mathematical complement pair without checking the exact official denominator definitions.",
  "- For PCA interpretation, both variables should not remain together in the final curated index, because they would overweight one SGB-II-related dimension.",
  "- My recommendation is to keep `a_ALGII_SGBII` as the more directly interpretable deprivation signal and to drop `a_Unterkunft_SGBII` from the curated main PCA.",
  "",
  "## 4. Corridor full-INKAR layer",
  "",
  "- Full original-code corridor layer:",
  "  `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/gpkg/corridor_full_inkar_original_codes.gpkg`",
  "- One-row-per-municipality CSV:",
  "  `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/tables/corridor_full_inkar_original_codes.csv`",
  sprintf("- Corridor municipalities in geometry layer: `%s`", 835),
  sprintf("- Municipalities with any INKAR row in the latest-build table: `%s`", 834),
  sprintf("- Original INKAR indicators available in the corridor build: `%s`", length(all_codes)),
  "- One municipality (`16076094`, Berga-Wünschendorf) remains missing in the raw latest-build table.",
  "",
  "## 5. What the old 51-variable PCA was actually doing",
  "",
  "- The old 51-variable set was a manually assembled, broad exploratory screening set.",
  "- It mixed poverty/welfare, detailed age composition, household structure, income composition, health-care supply, digital infrastructure, accessibility, and density in one block.",
  "- That is useful to let the data show broad covariance structure.",
  "- But it also means that some concepts appear several times, for example:",
  "  - age structure through multiple nested age shares plus dependency ratios;",
  "  - household income through low, medium, and high shares together;",
  "  - health access through both total physicians and specialist-specific sub-indicators;",
  "  - digital access through several broadband thresholds;",
  "  - SGB-II deprivation through several related but definition-sensitive indicators.",
  "- So the old 51-variable PCA is defensible as an exploratory starting point, but not ideal as the final, most interpretable vulnerability index.",
  "",
  "## 6. Initial correlation diagnosis",
  "",
  "- Original 51-variable correlation heatmap:",
  "  `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/old51_correlation_heatmap.png`",
  "- All-176-variable correlation heatmap:",
  "  `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/all176_correlation_heatmap.png`",
  "- Curated-variable correlation heatmap:",
  "  `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/curated_correlation_heatmap.png`",
  "",
  "The correlation step should come first because it makes duplication visible before PCA loadings are interpreted.",
  "",
  "Top correlation pairs from the old 51-variable set are saved here:",
  "- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/old51_top_correlation_pairs.csv`",
  "",
  "## 7. PCA comparison",
  "",
  "Comparison across model versions:",
  comparison_lines,
  "",
  "Interpretation:",
  "- The all-176 PCA is useful as a broad structural diagnostic, not as the main final index.",
  "- The original 51-variable PCA is a cleaner exploratory subset than all 176, but it still contains several duplicated constructs.",
  "- The curated PCA sacrifices breadth for interpretability and should be the strongest candidate for the main thesis index.",
  "",
  "## 8. Proposed curated variable set for the main vulnerability PCA",
  "",
  "These are the proposed core variables for the main interpretable PCA:",
  vapply(seq_len(nrow(selected_review %>% filter(in_curated))), function(i) {
    df <- selected_review %>% filter(in_curated)
    sprintf(
      "- `%s` — %s | coverage `%s%%` | %s",
      df$Kuerzel[i],
      df$Indikator[i],
      format(df$coverage_pct[i], nsmall = 2),
      df$handling_note[i]
    )
  }, character(1)),
  "",
  "Why this curated set is stronger:",
  "- it keeps one representative for each main vulnerability mechanism instead of several near-duplicates;",
  "- it avoids explicit definition conflicts like `a_BG1P` and the SGB-II housing-cost share problem;",
  "- it reduces compositional over-weighting from nested age, income, physician, and broadband blocks;",
  "- it stays close to the literature-based logic of deprivation, demographic sensitivity, household/social structure, and accessibility/adaptive capacity.",
  "",
  "## 9. Selected variables that need explicit caution",
  "",
  vapply(seq_len(nrow(selected_review)), function(i) {
    sprintf(
      "- `%s` — %s | coverage `%s%%` | status `%s` | %s",
      selected_review$Kuerzel[i],
      selected_review$Indikator[i],
      format(selected_review$coverage_pct[i], nsmall = 2),
      selected_review$proposed_status[i],
      selected_review$handling_note[i]
    )
  }, character(1)),
  "",
  "## 10. Full corridor variable inventory",
  "",
  "The complete review table is saved here:",
  "- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/corridor_variable_review_table.csv`",
  "",
  "Below, all corridor variables are listed by broad socio-economic dimension with a proposed handling note.",
  "",
  unlist(dimension_lines, use.names = FALSE),
  "## 11. Recommended next step",
  "",
  "- Keep the original 51-variable PCA as an exploratory comparison / appendix result.",
  "- Use the all-176 PCA only as a diagnostic exercise, not as the main thesis index.",
  "- Build the main vulnerability index from the curated set, then inspect its loadings and signs carefully.",
  "- If needed, we can next turn this directly into a keep/review/drop decision meeting note or a reply draft to the supervisors."
)

report_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/INKAR_VARIABLE_REVIEW_AND_PCA.md"
writeLines(report_lines, report_path)

writeLines(log_lines, log_path)
log_message("Wrote report: ", report_path)
log_message("Done.")
