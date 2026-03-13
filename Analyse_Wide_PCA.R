#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(ggspatial)
})

options(scipen = 999)
set.seed(42)

municipality_file <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/ANALYSIS.nosync/gemeinden_elbe_final_full.gpkg"
elbe_file <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PHYSISCH.nosync/RIGHT PROJECTION/Elbe.gpkg"

output_root <- Sys.getenv("ANALYSE_WIDE_PCA_OUT_DIR", unset = "outputs_wide_pca")
plot_dir <- file.path(output_root, "plots")
table_dir <- file.path(output_root, "tables")
map_dir <- file.path(output_root, "maps")

for (dir_path in c(output_root, plot_dir, table_dir, map_dir)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

save_plot <- function(plot_obj, filename, width = 8, height = 5, subdir = "plots") {
  ggsave(
    filename = file.path(output_root, subdir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
}

save_table <- function(data_obj, filename) {
  write_csv(data_obj, file.path(table_dir, filename))
}

map_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      plot.caption = element_text(size = base_size - 3, hjust = 1)
    )
}

map_annotations <- function() {
  list(
    annotation_north_arrow(
      location = "tl",
      which_north = "true",
      style = north_arrow_fancy_orienteering,
      height = unit(1.2, "cm"),
      width = unit(1.2, "cm")
    ),
    annotation_scale(
      location = "bl",
      width_hint = 0.25,
      line_width = 0.6,
      text_cex = 0.7
    )
  )
}

impute_median <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

build_weighted_index <- function(spatial_data, n_components, variance_table) {
  component_weights <- variance_table$Variance[seq_len(n_components)]
  names(component_weights) <- paste0("PC", seq_len(n_components))

  component_matrix <- spatial_data %>%
    st_drop_geometry() %>%
    select(all_of(names(component_weights))) %>%
    as.matrix()

  raw_index <- as.numeric(component_matrix %*% component_weights)
  z_index <- as.numeric(scale(raw_index))

  list(raw = raw_index, z = z_index)
}

find_current_script <- function() {
  script_arg <- commandArgs(trailingOnly = FALSE)
  script_arg <- script_arg[grepl("^--file=", script_arg)]
  if (length(script_arg) == 0) return(NA_character_)
  sub("^--file=", "", script_arg[1])
}

export_code_pdf <- function(script_path, pdf_path, lines_per_page = 58, cex = 0.58) {
  if (is.na(script_path) || !file.exists(script_path)) return(invisible(NULL))

  source_lines <- readLines(script_path, warn = FALSE, encoding = "UTF-8")
  pdf_lines <- iconv(source_lines, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "?")
  n_pages <- ceiling(length(pdf_lines) / lines_per_page)

  grDevices::pdf(pdf_path, width = 8.27, height = 11.69, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (page_id in seq_len(n_pages)) {
    line_start <- ((page_id - 1) * lines_per_page) + 1
    line_end <- min(page_id * lines_per_page, length(pdf_lines))
    page_lines <- pdf_lines[line_start:line_end]
    line_numbers <- seq.int(line_start, line_end)
    y_positions <- seq(0.96, 0.05, length.out = length(page_lines))

    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    title(
      main = basename(script_path),
      sub = sprintf("Page %d of %d", page_id, n_pages),
      line = -1,
      adj = 0
    )

    text(
      x = rep(0.02, length(page_lines)),
      y = y_positions,
      labels = sprintf("%3d  %s", line_numbers, page_lines),
      adj = c(0, 1),
      cex = cex,
      family = "mono"
    )
  }

  invisible(pdf_path)
}

municipalities <- st_read(municipality_file, quiet = TRUE)
elbe_river <- st_read(elbe_file, quiet = TRUE)
elbe_river <- st_zm(elbe_river, drop = TRUE, what = "ZM")
elbe_river <- st_make_valid(elbe_river)

if (is.na(st_crs(elbe_river))) {
  stop("The Elbe layer has no CRS.")
}
if (st_crs(elbe_river) != st_crs(municipalities)) {
  elbe_river <- st_transform(elbe_river, st_crs(municipalities))
}
elbe_river <- suppressWarnings(st_cast(elbe_river, "LINESTRING"))

pca_columns <- c(
  "exposure_share_alg2_sgb2",
  "exposure_share_bg_single_parent",
  "exposure_share_bg_5plus",
  "exposure_share_bg_with_children",
  "exposure_share_sgb2_with_housing_costs",
  "exposure_share_longterm_unemp",
  "exposure_unemp_u25_per_1000",
  "exposure_unemp_55plus_per_1000",
  "exposure_income_tax_per_capita",
  "exposure_purchasing_power",
  "exposure_trade_tax_per_capita",
  "exposure_tax_revenue_total",
  "exposure_share_age_0_3",
  "exposure_share_age_3_6",
  "exposure_share_age_6_18",
  "exposure_share_age_18_25",
  "exposure_share_age_25_30",
  "exposure_share_age_30_50",
  "exposure_share_age_50_65",
  "exposure_share_age_65_75",
  "exposure_share_age_65plus",
  "exposure_share_age_75plus",
  "exposure_natural_pop_change",
  "exposure_migration_balance",
  "exposure_old_age_dependency",
  "exposure_youth_dependency",
  "exposure_young_old_ratio",
  "exposure_share_single_households",
  "exposure_share_households_with_children",
  "exposure_share_hh_income_high",
  "exposure_share_hh_income_medium",
  "exposure_share_hh_income_low",
  "exposure_students_total_per_1000",
  "exposure_students_18_25_per_1000",
  "exposure_students_fh_per_1000",
  "exposure_gp_general_per_1000",
  "exposure_gp_primary_per_1000",
  "exposure_internists_per_1000",
  "exposure_pediatricians_per_1000_children",
  "exposure_doctors_total_per_1000",
  "exposure_share_bb_1000mbit",
  "exposure_share_bb_100mbit",
  "exposure_share_bb_50mbit",
  "exposure_share_4g",
  "exposure_dist_supermarket_m",
  "exposure_dist_pharmacy_m",
  "exposure_dist_gp_m",
  "exposure_dist_public_transport_m",
  "exposure_dist_primary_school_m",
  "exposure_pop_density_per_km2",
  "exposure_employment_density_per_km2"
)

missing_pca_columns <- setdiff(pca_columns, names(municipalities))
if (length(missing_pca_columns) > 0) {
  stop("Missing PCA columns: ", paste(missing_pca_columns, collapse = ", "))
}

wide_pca_data <- municipalities %>%
  st_drop_geometry() %>%
  select(all_of(pca_columns)) %>%
  mutate(across(everything(), as.numeric)) %>%
  rename_with(~ str_remove(.x, "^exposure_"))

wide_pca_imputed <- wide_pca_data %>%
  mutate(across(everything(), impute_median))

wide_pca_scaled <- scale(wide_pca_imputed)
wide_pca <- prcomp(wide_pca_scaled)

component_scores <- as.data.frame(wide_pca$x)
names(component_scores) <- paste0("PC", seq_len(ncol(component_scores)))
municipalities <- bind_cols(municipalities, component_scores)

eigenvalues <- wide_pca$sdev^2
variance_share <- eigenvalues / sum(eigenvalues)

scree_table <- tibble(
  PC = seq_along(eigenvalues),
  Eigenvalue = eigenvalues,
  Variance = variance_share,
  Cumulative = cumsum(variance_share)
)

main_index <- build_weighted_index(municipalities, n_components = 8, variance_table = scree_table)
municipalities$vuln_index_main <- main_index$raw
municipalities$vuln_index_main_z <- main_index$z

sensitivity_index <- build_weighted_index(municipalities, n_components = 12, variance_table = scree_table)
municipalities$vuln_index_sens12 <- sensitivity_index$raw
municipalities$vuln_index_sens12_z <- sensitivity_index$z

anchor_correlation <- suppressWarnings(
  cor(
    municipalities$vuln_index_main_z,
    municipalities$exposure_share_alg2_sgb2,
    use = "complete.obs"
  )
)

if (!is.na(anchor_correlation) && anchor_correlation < 0) {
  municipalities$vuln_index_main <- -municipalities$vuln_index_main
  municipalities$vuln_index_main_z <- -municipalities$vuln_index_main_z
  municipalities$vuln_index_sens12 <- -municipalities$vuln_index_sens12
  municipalities$vuln_index_sens12_z <- -municipalities$vuln_index_sens12_z
  message("Index flipped: higher values now indicate higher vulnerability.")
} else {
  message("Index direction OK.")
}

scree_plot <- ggplot(scree_table, aes(x = PC, y = Eigenvalue)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Scree Plot (Kaiser)", x = "PC", y = "Eigenvalue") +
  theme_classic()

cumulative_variance_plot <- ggplot(scree_table, aes(x = PC, y = Cumulative)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.7, linetype = "dashed") +
  labs(title = "Cumulative Variance", x = "PC", y = "Cumulative") +
  theme_classic()

save_plot(scree_plot, "scree_kaiser.png")
save_plot(cumulative_variance_plot, "cumulative_variance.png")
save_table(scree_table, "scree_table.csv")

loading_table <- as.data.frame(wide_pca$rotation)
loading_table$variable <- rownames(loading_table)

top_loadings <- loading_table %>%
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "loading") %>%
  mutate(abs_loading = abs(loading)) %>%
  group_by(PC) %>%
  arrange(desc(abs_loading)) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  arrange(as.integer(str_remove(PC, "PC")), desc(abs_loading))

save_table(top_loadings, "pca_top_loadings_top8_per_pc.csv")

map_vulnerability <- ggplot(municipalities) +
  geom_sf(aes(fill = vuln_index_main_z), color = NA) +
  geom_sf(data = elbe_river, inherit.aes = FALSE, color = "black", linewidth = 0.35) +
  coord_sf(crs = 25832) +
  map_annotations() +
  labs(
    title = "Socio-economic vulnerability index (PCA-weighted)",
    subtitle = "Variance-weighted index based on the first eight principal components",
    fill = "Index (z)",
    caption = "Source: INKAR (BBSR). Own processing."
  ) +
  scale_fill_viridis_c(option = "C") +
  map_theme()

save_plot(map_vulnerability, "map_vulnerability_index.png", width = 9, height = 7, subdir = "maps")

pc_titles <- c(
  PC1 = "PC1: Socio-economic disadvantage and welfare dependency",
  PC2 = "PC2: Urbanisation, density and accessibility",
  PC3 = "PC3: Demographic ageing and dependency",
  PC4 = "PC4: Education, students and human capital"
)

pc_subtitles <- c(
  PC1 = "High loadings: ALG II, low income, unemployment, single households",
  PC2 = "High loadings: population density, transport access, broadband",
  PC3 = "High loadings: 65+, old-age dependency, low youth share",
  PC4 = "High loadings: students, higher income, education indicators"
)

for (pc_name in paste0("PC", 1:4)) {
  pc_map <- ggplot(municipalities) +
    geom_sf(aes(fill = .data[[pc_name]]), color = NA) +
    geom_sf(data = elbe_river, inherit.aes = FALSE, color = "black", linewidth = 0.3) +
    coord_sf(crs = 25832) +
    map_annotations() +
    labs(
      title = pc_titles[[pc_name]],
      subtitle = pc_subtitles[[pc_name]],
      fill = pc_name,
      caption = "Source: INKAR (BBSR). PCA-based component scores. Own processing."
    ) +
    scale_fill_viridis_c(option = "C") +
    map_theme()

  save_plot(pc_map, paste0("map_", pc_name, ".png"), width = 9, height = 7, subdir = "maps")
}

save_table(tibble(variable = names(wide_pca_data)), "wide_pca_variables.csv")
saveRDS(municipalities, file.path(output_root, "FULL_with_PCs_index.rds"))

if (!identical(Sys.getenv("ANALYSE_WIDE_PCA_SKIP_CODE_PDF", unset = "FALSE"), "TRUE")) {
  script_path <- find_current_script()
  export_code_pdf(
    script_path = script_path,
    pdf_path = file.path(output_root, "Analyse_Wide_PCA_code.pdf")
  )
}

message("DONE. Outputs saved to: ", normalizePath(output_root))
