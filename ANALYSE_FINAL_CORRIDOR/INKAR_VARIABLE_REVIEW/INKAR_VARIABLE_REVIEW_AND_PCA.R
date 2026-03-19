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

thesis_candidate_codes <- c(
  "a_ALGII_SGBII",
  "q_newfBGu15_bev",
  "a_aloLang",
  "q_alo_u25_einw",
  "q_alo_ü55_einw",
  "a_Minijobs",
  "q_svw",
  "q_kaufkraft",
  "q_einkst_bev",
  "a_hheink_niedrig",
  "a_bev65um",
  "a_bev_0006",
  "q_HH1",
  "a_hh_kind",
  "a_ewfBG_allein",
  "m_G02_SUP_DIST",
  "m_Q01_APO_DIST",
  "m_Q07_HA_DIST",
  "m_OEV20_DIST",
  "m_P01_PRIM_DIST",
  "q_ärzte_bev",
  "a_bb_50Mbits",
  "a_bb_4G"
)

student_sensitivity_codes <- c(
  thesis_candidate_codes,
  "q_stud",
  "q_stud_1825",
  "q_stud_fh"
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
if (length(setdiff(thesis_candidate_codes, names(full_corridor))) > 0) {
  stop("Some thesis-candidate variables are missing in the full corridor table.")
}
if (length(setdiff(student_sensitivity_codes, names(full_corridor))) > 0) {
  stop("Some student-sensitivity variables are missing in the full corridor table.")
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
landuse_context <- setdiff(landuse_context, "a_gest_bev")
dynamics_context <- c("i_saldo_nat", "i_wans", "a_geb_bev", "a_gest_bev")

inventory <- inventory %>%
  mutate(
    Indikator = if_else(Kuerzel == "a_bb_4G", "4G-Mobilfunkverfügbarkeit", Indikator)
  )

review_table <- inventory %>%
  mutate(
    in_old51 = Kuerzel %in% old51_codes,
    in_thesis_candidate = Kuerzel %in% thesis_candidate_codes,
    in_student_sensitivity = Kuerzel %in% student_sensitivity_codes,
    in_curated = Kuerzel %in% curated_codes,
    proposed_status = case_when(
      in_curated ~ "core_curated_pca",
      in_thesis_candidate ~ "core_thesis_candidate",
      Kuerzel %in% education_context ~ "student_sensitivity_only",
      Kuerzel == "a_Unterkunft_SGBII" ~ "drop_overlap_with_algii",
      Kuerzel == "a_BG1P" ~ "drop_definition_mismatch",
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
      in_thesis_candidate ~ "Keep in the wider thesis-candidate PCA set; broad enough for a wide PCA, but cleaner than the original 51-variable block.",
      Kuerzel %in% education_context ~ "Use only in a student-sensitivity PCA. These variables mainly capture university-town structure rather than core municipality-level social vulnerability.",
      Kuerzel == "a_Unterkunft_SGBII" ~ "Conceptually overlaps with a_ALGII_SGBII and is difficult to interpret as an independent vulnerability signal.",
      Kuerzel == "a_BG1P" ~ "Officially this is Einpersonen-Bedarfsgemeinschaften, not single-parent households; the old label was misleading.",
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
write_csv(review_table %>% filter(in_thesis_candidate), file.path(paths$output_dir, "tables", "thesis_candidate_variable_list.csv"))
write_csv(review_table %>% filter(in_student_sensitivity & Kuerzel %in% student_sensitivity_codes), file.path(paths$output_dir, "tables", "student_sensitivity_variable_list.csv"))
write_csv(review_table %>% filter(in_curated), file.path(paths$output_dir, "tables", "curated_variable_list.csv"))

log_message("Running original 51-variable PCA review bundle.")
old51_bundle <- run_pca_bundle(full_corridor, old51_codes, "old51", "Original 51-variable set")

log_message("Running all-176-variable PCA review bundle.")
all176_bundle <- run_pca_bundle(full_corridor, all_codes, "all176", "All 176 corridor indicators")

log_message("Running thesis-candidate PCA review bundle.")
thesis_candidate_bundle <- run_pca_bundle(full_corridor, thesis_candidate_codes, "thesis_candidate", "Thesis-candidate set")

log_message("Running student-sensitivity PCA review bundle.")
student_sensitivity_bundle <- run_pca_bundle(full_corridor, student_sensitivity_codes, "student_sensitivity", "Thesis candidate plus student sensitivity set")

log_message("Running curated PCA review bundle.")
curated_bundle <- run_pca_bundle(full_corridor, curated_codes, "curated", "Curated vulnerability set")

comparison <- tibble(
  model = c("original_51", "all_176", "thesis_candidate_23", "student_sensitivity_26", "curated_17"),
  variables = c(length(old51_codes), length(all_codes), length(thesis_candidate_codes), length(student_sensitivity_codes), length(curated_codes)),
  municipalities_used = c(nrow(old51_bundle$data_imputed), nrow(all176_bundle$data_imputed), nrow(thesis_candidate_bundle$data_imputed), nrow(student_sensitivity_bundle$data_imputed), nrow(curated_bundle$data_imputed)),
  pc1_variance = c(old51_bundle$scree$variance[1], all176_bundle$scree$variance[1], thesis_candidate_bundle$scree$variance[1], student_sensitivity_bundle$scree$variance[1], curated_bundle$scree$variance[1]),
  pc1_to_pc4_cumulative = c(old51_bundle$scree$cumulative[4], all176_bundle$scree$cumulative[4], thesis_candidate_bundle$scree$cumulative[4], student_sensitivity_bundle$scree$cumulative[4], curated_bundle$scree$cumulative[4]),
  pc1_to_pc8_cumulative = c(
    old51_bundle$scree$cumulative[8],
    all176_bundle$scree$cumulative[8],
    thesis_candidate_bundle$scree$cumulative[8],
    student_sensitivity_bundle$scree$cumulative[8],
    curated_bundle$scree$cumulative[min(8, nrow(curated_bundle$scree))]
  )
)
write_csv(comparison, file.path(paths$output_dir, "tables", "pca_model_comparison.csv"))

pair_check <- read_csv(
  "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/tables/corridor_sgb2_pair_check.csv",
  show_col_types = FALSE
)

selected_review_codes <- c(
  thesis_candidate_codes,
  student_sensitivity_codes,
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
  "- The thesis-candidate set is the pragmatic middle ground: broader than the strict curated core, but already cleaned from the most obvious definition conflicts and duplicate blocks.",
  "- The student-sensitivity set adds the student concentration indicators only as a robustness run, not as part of the main thesis PCA.",
  "- The curated PCA sacrifices breadth for interpretability and remains the strictest option.",
  "",
  "## 8. Thesis-candidate variable set from all available corridor indicators",
  "",
  "This is the proposed working set for the thesis right now: broader than the strict 17-variable core, but already cleaned enough to remove the most problematic overlaps.",
  vapply(seq_len(nrow(selected_review %>% filter(Kuerzel %in% thesis_candidate_codes) %>% arrange(match(Kuerzel, thesis_candidate_codes)))), function(i) {
    df <- selected_review %>% filter(Kuerzel %in% thesis_candidate_codes) %>% arrange(match(Kuerzel, thesis_candidate_codes))
    sprintf(
      "- `%s` — %s | coverage `%s%%` | %s",
      df$Kuerzel[i],
      df$Indikator[i],
      format(df$coverage_pct[i], nsmall = 2),
      df$handling_note[i]
    )
  }, character(1)),
  "",
  "Working logic of the thesis-candidate set:",
  "- it keeps more nuance than the strict curated 17, especially within labour market strain and accessibility;",
  "- it still removes the clearest duplicates and definition problems;",
  "- it is broad enough for a meaningful PCA, but much easier to interpret than the old 51-variable block;",
  "- it avoids keeping several variables that say almost the same thing at the same time.",
  "",
  "Thesis-candidate correlation heatmap:",
  "- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/thesis_candidate_correlation_heatmap.png`",
  "",
  "Top correlation pairs from the thesis-candidate set:",
  "- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/thesis_candidate_top_correlation_pairs.csv`",
  "",
  "## 9. Student sensitivity block",
  "",
  "Because student populations can be socio-economically precarious but are also strongly tied to urban university locations, student indicators are not included in the main thesis-candidate PCA.",
  "- Instead, they are added only in a dedicated sensitivity run.",
  "- Student sensitivity variable list: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/student_sensitivity_variable_list.csv`",
  "- Student sensitivity correlation heatmap: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/student_sensitivity_correlation_heatmap.png`",
  "- Student sensitivity top correlation pairs: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/student_sensitivity_top_correlation_pairs.csv`",
  "- Student sensitivity scree table: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/student_sensitivity_scree.csv`",
  "",
  "## 10. Proposed strict curated variable set for the main vulnerability PCA",
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
  "## 11. Selected variables that need explicit caution",
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
  "## 12. Full corridor variable inventory",
  "",
  "The complete review table is saved here:",
  "- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/corridor_variable_review_table.csv`",
  "",
  "Below, all corridor variables are listed by broad socio-economic dimension with a proposed handling note.",
  "",
  unlist(dimension_lines, use.names = FALSE),
  "## 13. Recommended next step",
  "",
  "- Keep the original 51-variable PCA as an exploratory comparison / appendix result.",
  "- Use the all-176 PCA only as a diagnostic exercise, not as the main thesis index.",
  "- Use the thesis-candidate set as the current working set for substantive thesis analyses.",
  "- Use the student-sensitivity run only as a robustness check for how strongly university-center structure changes the PCA.",
  "- Keep the strict curated 17-variable set as a robustness / interpretability check.",
  "- If needed, we can next turn this directly into a keep/review/drop decision meeting note or a reply draft to the supervisors."
)

report_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/INKAR_VARIABLE_REVIEW_AND_PCA.md"
writeLines(report_lines, report_path)

all_inventory_lines <- lapply(split(review_table, review_table$broad_dimension), function(df) {
  c(
    sprintf("### %s", unique(df$broad_dimension)),
    "",
    sprintf("- Indicators in this block: `%s`", nrow(df)),
    sprintf("- Median coverage: `%s%%`", format(round(median(df$coverage_pct, na.rm = TRUE), 2), nsmall = 2)),
    "",
    vapply(seq_len(nrow(df)), function(i) {
      sprintf(
        "- `%s` — %s | coverage `%s%%` | status `%s`",
        df$Kuerzel[i],
        df$Indikator[i],
        format(df$coverage_pct[i], nsmall = 2),
        df$proposed_status[i]
      )
    }, character(1)),
    ""
  )
})

overview_lines <- c(
  "# INKAR HQ500 OVERVIEW",
  "",
  sprintf("Created: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 1. Study area",
  "",
  "The study area is defined as the HQ500 Elbe flood corridor.",
  "- Municipalities were retained if they intersect the RP500 flood extent derived from the prepared EFAS flood rasters.",
  "- This replaces the older basin-wide definition and focuses the analysis on municipalities that are actually flood-relevant under an extreme return period.",
  "- Corridor municipalities in the current geometry layer: `835`.",
  "- Municipalities with a latest available raw-INKAR record after corridor filtering: `834`.",
  "- One municipality (`16076094`, Berga-Wünschendorf) remains missing in the current raw-INKAR latest-build table and should be documented explicitly.",
  "",
  "Core corridor files:",
  "- Geometry + all original INKAR codes: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/gpkg/corridor_full_inkar_original_codes.gpkg`",
  "- Flat table + all original INKAR codes: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/tables/corridor_full_inkar_original_codes.csv`",
  "",
  "## 2. All available INKAR variables in the corridor",
  "",
  "- Number of original INKAR indicators available in the corridor latest-build table: `176`.",
  "- Inventory with official names, coverage, and broad dimensions: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/tables/corridor_inkar_indicator_inventory_dimensioned.csv`",
  "- Short summary by broad dimension: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_FULL_CORRIDOR/outputs/tables/corridor_inkar_dimension_summary.csv`",
  "",
  "All available corridor indicators grouped by broad dimension:",
  "",
  unlist(all_inventory_lines, use.names = FALSE),
  "## 3. Correlation matrix of all available corridor variables",
  "",
  "The full 176-variable correlation matrix is useful as a diagnostic step, not as a final model input on its own.",
  "- Correlation heatmap: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/all176_correlation_heatmap.png`",
  "- Top absolute correlation pairs: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/all176_top_correlation_pairs.csv`",
  "",
  "Why this step matters:",
  "- it shows where several indicators measure the same underlying construct;",
  "- it makes compositional counterparts visible;",
  "- it helps identify definition-sensitive indicators before PCA loadings are interpreted.",
  "",
  "## 4. Original 51-variable wide PCA set",
  "",
  "The original 51-variable set came from the earlier INKAR preprocessing script and was designed as a broad exploratory screening block.",
  "- Source logic: one manually assembled selection spanning poverty/welfare, unemployment, demography, dependency, household structure, income groups, health care, digital infrastructure, accessibility, and density.",
  "- This was a reasonable exploratory starting point because it allowed the covariance structure of a broad vulnerability field to emerge from the data.",
  "- At the same time, the block was not yet fully cleaned for duplicated constructs or definition conflicts.",
  "",
  "Original 51-variable files:",
  "- Variable list: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/original_wide51_variable_list.csv`",
  "- Correlation heatmap: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/old51_correlation_heatmap.png`",
  "- Top correlation pairs: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/old51_top_correlation_pairs.csv`",
  "- Scree table: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/old51_scree.csv`",
  "",
  "Main methodological issue with the old 51-variable set:",
  "- several variables partly say the same thing;",
  "- nested age blocks over-weight demographic composition;",
  "- multiple broadband thresholds duplicate one digital-access construct;",
  "- specialist-specific physician indicators overlap strongly with total physician supply;",
  "- SGB-II-related indicators include at least one definition-sensitive pair and one naming mismatch (`a_BG1P`).",
  "",
  "## 5. Working thesis candidate set from all available variables",
  "",
  "This is the recommended working set for the thesis at the current stage: broader than the strict 17-variable core, but cleaned enough to avoid the clearest double-counting problems.",
  "",
  vapply(seq_len(length(thesis_candidate_codes)), function(i) {
    code <- thesis_candidate_codes[i]
    row <- review_table %>% filter(Kuerzel == code) %>% slice(1)
    sprintf(
      "- `%s` — %s | coverage `%s%%` | %s",
      row$Kuerzel,
      row$Indikator,
      format(row$coverage_pct, nsmall = 2),
      row$handling_note
    )
  }, character(1)),
  "",
  "Selection logic of the thesis candidate set:",
  "- keep a broad but interpretable coverage of deprivation, labour-market strain, demographic sensitivity, household/social structure, accessibility, health access, and digital infrastructure;",
  "- remove direct duplicates or near-duplicates wherever possible;",
  "- remove variables with clearly problematic definitions for interpretation;",
  "- keep the set broad enough for PCA while avoiding obvious over-weighting of the same latent construct.",
  "",
  "Thesis-candidate files:",
  "- Variable list: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/thesis_candidate_variable_list.csv`",
  "- Correlation heatmap: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/thesis_candidate_correlation_heatmap.png`",
  "- Top correlation pairs: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/thesis_candidate_top_correlation_pairs.csv`",
  "- Scree table: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/thesis_candidate_scree.csv`",
  "- Top loadings: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/thesis_candidate_top_loadings.csv`",
  "",
  "Student variables are intentionally excluded from this main set and are handled only in a separate sensitivity PCA.",
  "",
  "## 6. Student sensitivity block",
  "",
  "- Student sensitivity variable list: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/student_sensitivity_variable_list.csv`",
  "- Student sensitivity correlation heatmap: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/plots/student_sensitivity_correlation_heatmap.png`",
  "- Student sensitivity top correlation pairs: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/student_sensitivity_top_correlation_pairs.csv`",
  "- Student sensitivity scree table: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/outputs/tables/student_sensitivity_scree.csv`",
  "",
  "Rationale:",
  "- student populations can be economically precarious;",
  "- but at municipality scale, student indicators also capture university-center structure very strongly;",
  "- therefore they are better treated as a sensitivity block than as part of the main thesis index.",
  "",
  "## 7. PCA comparison across all sets",
  "",
  comparison_lines,
  "",
  "Interpretation:",
  "- `all 176` is too broad to be the main thesis PCA, but very useful as a structural diagnostic.",
  "- `original 51` remains useful as an exploratory benchmark and historical reference.",
  "- `thesis candidate 23` is currently the most balanced working set for the thesis.",
  "- `student sensitivity 26` shows what changes once student concentration variables are allowed into the PCA.",
  "- `curated 17` remains the stricter fallback / robustness set.",
  "",
  "## 8. Immediate methodological takeaway",
  "",
  "- Yes, the old wide PCA made sense as an exploratory first step.",
  "- No, it should probably not remain the only final vulnerability index without additional variable curation.",
  "- The safer thesis strategy is to document three layers clearly:",
  "  - `all 176` for inventory and structural diagnosis,",
  "  - `original 51` for exploratory comparison,",
  "  - `thesis candidate 23` as the main working PCA set,",
  "  - `student sensitivity 26` as an additional robustness run.",
  "",
  "## 9. Related note",
  "",
  "A more detailed review note with additional handling commentary is saved here:",
  "- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_VARIABLE_REVIEW/INKAR_VARIABLE_REVIEW_AND_PCA.md`"
)

overview_path <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/INKAR_HQ500_OVERVIEW.md"
writeLines(overview_lines, overview_path)

writeLines(log_lines, log_path)
log_message("Wrote report: ", report_path)
log_message("Wrote overview: ", overview_path)
log_message("Done.")
