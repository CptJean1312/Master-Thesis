#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

# ============================================================
# Reproducible INKAR availability check for Germany
# ------------------------------------------------------------
# Goal:
# - filter the original INKAR 2025 raw CSV to municipality level
# - keep the latest available value per municipality x indicator
# - check which indicators are present in all 16 federal states
# - check which indicators have full AGS-level coverage
# - summarise results by official and broad socio-economic dimensions
# ============================================================

raw_file <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/SOCIOECONOMIC.nosync/Downloads/inkar_2025/inkar_2025.csv"
meta_file <- "/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/SOCIOECONOMIC.nosync/INKAR Uebersicht der Indikatoren.xlsx"
out_dir <- file.path(getwd(), "outputs_inkar_overview")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(raw_file)) {
  stop("Raw INKAR CSV not found: ", raw_file)
}

if (!file.exists(meta_file)) {
  stop("INKAR metadata workbook not found: ", meta_file)
}

all_states <- sprintf("%02d", 1:16)
broad_order <- c(
  "Labour market and employment",
  "Demography and household structure",
  "Health, services and accessibility",
  "Housing and land-use context",
  "Income, deprivation and social transfers",
  "Basic population and structural counts",
  "Education",
  "Economic structure and tourism",
  "Other / uncategorised"
)

save_csv <- function(dt, filename) {
  fwrite(dt, file.path(out_dir, filename))
}

safe_empty_to_na <- function(x) {
  if (!is.character(x)) return(x)
  y <- trimws(x)
  y[y == ""] <- NA_character_
  y
}

message("1) Reading INKAR metadata workbook ...")

meta_raw <- read_excel(
  meta_file,
  sheet = "Raumbeobachtung DE",
  skip = 1,
  col_names = c(
    "Kurzname", "Name", "Algorithmus", "M_ID", "Kuerzel",
    "Anmerkungen", "Grundlagen", "Gemeinden", "Kreise"
  )
)

meta_raw[] <- lapply(meta_raw, safe_empty_to_na)
meta_raw$row_id <- seq_len(nrow(meta_raw))

current_main <- NA_character_
current_sub <- NA_character_
main_vec <- character(nrow(meta_raw))
sub_vec <- character(nrow(meta_raw))

for (i in seq_len(nrow(meta_raw))) {
  short_label <- meta_raw$Kurzname[i]
  kuerzel <- meta_raw$Kuerzel[i]

  if (is.na(kuerzel) && !is.na(short_label) && short_label != "Kurzname") {
    if (grepl("\\s[-\u2013]\\s", short_label)) {
      current_sub <- short_label
    } else {
      current_main <- short_label
      current_sub <- NA_character_
    }
  }

  main_vec[i] <- current_main
  sub_vec[i] <- current_sub
}

meta <- as.data.table(meta_raw)[
  !is.na(Kuerzel) & Kuerzel != "K\u00fcrzel",
  .(
    Kuerzel,
    official_short = Kurzname,
    official_name = Name,
    main_dimension = main_vec[row_id],
    sub_dimension = sub_vec[row_id],
    period_gemeinden = Gemeinden,
    period_kreise = Kreise
  )
]

# The metadata overview sheet does not assign a clean label to a_bb_4G.
meta[Kuerzel == "a_bb_4G", `:=`(
  official_short = fifelse(
    is.na(official_short),
    "4G coverage (metadata label missing in overview sheet)",
    official_short
  ),
  official_name = fifelse(
    is.na(official_name),
    "Share of households with 4G mobile-network coverage",
    official_name
  ),
  main_dimension = fifelse(
    is.na(main_dimension),
    "Erreichbarkeit",
    main_dimension
  ),
  sub_dimension = fifelse(
    is.na(sub_dimension),
    "Digitale Erreichbarkeit",
    sub_dimension
  )
)]

message("2) Reading raw INKAR CSV and filtering Germany-wide municipality data ...")

inkar <- fread(
  raw_file,
  sep = ";",
  quote = "",
  select = c("Bereich", "Raumbezug", "Kennziffer", "Kuerzel", "Wert", "Zeitbezug"),
  encoding = "UTF-8",
  showProgress = FALSE
)

inkar <- inkar[
  Bereich == "LRB" &
    Raumbezug == "Gemeinden" &
    substr(Kennziffer, 1, 2) %in% all_states
]

inkar[, state := substr(Kennziffer, 1, 2)]

message("3) Keeping latest available observation per municipality x indicator ...")

setorder(inkar, Kennziffer, Kuerzel, -Zeitbezug)
latest <- inkar[, .SD[1], by = .(Kennziffer, state, Kuerzel)]

latest[, value_chr := trimws(as.character(Wert))]
latest[, has_value := !is.na(Wert) & !(value_chr %in% c("", ".", "..", "-"))]

state_totals <- latest[, .(municipalities_total = uniqueN(Kennziffer)), by = state]

message("4) Calculating indicator availability by state ...")

coverage_state <- latest[
  has_value == TRUE,
  .(municipalities_with_value = uniqueN(Kennziffer)),
  by = .(Kuerzel, state)
]

coverage_state <- merge(
  CJ(Kuerzel = sort(unique(latest$Kuerzel)), state = sort(unique(latest$state))),
  coverage_state,
  by = c("Kuerzel", "state"),
  all.x = TRUE
)

coverage_state <- merge(coverage_state, state_totals, by = "state", all.x = TRUE)
coverage_state[is.na(municipalities_with_value), municipalities_with_value := 0L]
coverage_state[, state_coverage := municipalities_with_value / municipalities_total]

coverage_summary <- coverage_state[
  ,
  .(
    states_with_data = sum(municipalities_with_value > 0),
    all_states_present = all(municipalities_with_value > 0),
    full_coverage_all_states = all(state_coverage == 1),
    min_state_coverage = min(state_coverage),
    mean_state_coverage = mean(state_coverage),
    max_state_coverage = max(state_coverage),
    states_below_95 = sum(state_coverage < 0.95),
    states_below_100 = sum(state_coverage < 1),
    municipalities_with_value_total = sum(municipalities_with_value),
    municipalities_total = sum(municipalities_total)
  ),
  by = Kuerzel
]

coverage_summary <- merge(coverage_summary, meta, by = "Kuerzel", all.x = TRUE)

coverage_summary[Kuerzel == "a_bb_4G", `:=`(
  official_short = fifelse(
    is.na(official_short) | official_short == "",
    "Bandbreitenverfugbarkeit mindestens 4G",
    official_short
  ),
  official_name = fifelse(
    is.na(official_name) | official_name == "",
    "Anteil der Haushalte mit einer 4G-Mobilfunkversorgung in %",
    official_name
  ),
  main_dimension = fifelse(
    is.na(main_dimension) | main_dimension == "",
    "Erreichbarkeit",
    main_dimension
  ),
  sub_dimension = fifelse(
    is.na(sub_dimension) | sub_dimension == "",
    "Digitale Erreichbarkeit",
    sub_dimension
  )
)]

coverage_summary[, main_dimension := fifelse(
  is.na(main_dimension),
  "Other / uncategorised in metadata",
  main_dimension
)]

coverage_summary[, sub_dimension := fifelse(
  is.na(sub_dimension),
  "(no sub-dimension)",
  sub_dimension
)]

coverage_summary[, broad_dimension := fifelse(
  main_dimension %in% c("Bev\u00f6lkerung", "Siedlungsstruktur"),
  "Demography and household structure",
  fifelse(
    main_dimension %in% c("Arbeitslosigkeit", "Besch\u00e4ftigung und Erwerbst\u00e4tigkeit"),
    "Labour market and employment",
    fifelse(
      main_dimension %in% c("Sozialleistungen", "Privateinkommen und private Schulden", "\u00d6ffentliche Finanzen"),
      "Income, deprivation and social transfers",
      fifelse(
        main_dimension %in% c("Medizinische Versorgung", "Erreichbarkeit"),
        "Health, services and accessibility",
        fifelse(
          main_dimension %in% c("Bauen und Wohnen", "Fl\u00e4chennutzung"),
          "Housing and land-use context",
          fifelse(
            main_dimension == "Bildung",
            "Education",
            fifelse(
              main_dimension == "Wirtschaft",
              "Economic structure and tourism",
              fifelse(
                main_dimension == "Absolutzahlen",
                "Basic population and structural counts",
                "Other / uncategorised"
              )
            )
          )
        )
      )
    )
  )
)]

all_states_present_tbl <- coverage_summary[all_states_present == TRUE]
full_coverage_tbl <- coverage_summary[full_coverage_all_states == TRUE]
missing_states_tbl <- coverage_summary[all_states_present == FALSE]

summary_official <- all_states_present_tbl[
  ,
  .(
    n_indicators_present_all_states = .N,
    n_full_coverage_all_states = sum(full_coverage_all_states),
    mean_min_state_coverage = mean(min_state_coverage),
    mean_mean_state_coverage = mean(mean_state_coverage)
  ),
  by = .(main_dimension, sub_dimension)
][order(main_dimension, sub_dimension)]

summary_broad <- all_states_present_tbl[
  ,
  .(
    n_indicators = .N,
    n_full_coverage = sum(full_coverage_all_states)
  ),
  by = broad_dimension
]

summary_broad[, broad_dimension := factor(broad_dimension, levels = broad_order)]
setorder(summary_broad, broad_dimension)
summary_broad[, broad_dimension := as.character(broad_dimension)]

message("5) Writing output tables ...")

save_csv(
  coverage_summary[order(main_dimension, sub_dimension, official_short)],
  "inkar_indicator_availability_all_states.csv"
)

save_csv(
  all_states_present_tbl[order(main_dimension, sub_dimension, official_short)],
  "inkar_indicators_present_in_all_16_states.csv"
)

save_csv(
  full_coverage_tbl[order(main_dimension, sub_dimension, official_short)],
  "inkar_indicators_full_coverage_all_states.csv"
)

save_csv(
  missing_states_tbl[order(main_dimension, sub_dimension, official_short)],
  "inkar_indicators_missing_some_states.csv"
)

save_csv(summary_official, "inkar_availability_summary_by_official_dimension.csv")
save_csv(summary_broad, "inkar_availability_summary_by_broad_dimension.csv")

headline <- data.table(
  municipality_indicators_total = uniqueN(coverage_summary$Kuerzel),
  indicators_present_all_16_states = nrow(all_states_present_tbl),
  indicators_full_coverage_all_16_states = nrow(full_coverage_tbl),
  indicators_not_present_all_16_states = nrow(missing_states_tbl),
  municipalities_germany = uniqueN(latest$Kennziffer)
)

save_csv(headline, "inkar_germany_availability_headline.csv")

message("")
message("Headline result")
message("- municipality-level indicators in total: ", uniqueN(coverage_summary$Kuerzel))
message("- indicators present in all 16 federal states: ", nrow(all_states_present_tbl))
message("- indicators with full municipality coverage in all 16 federal states: ", nrow(full_coverage_tbl))
message("- indicators not present in all 16 federal states: ", nrow(missing_states_tbl))
message("")
message("Indicators not present in all 16 federal states:")
print(
  missing_states_tbl[
    ,
    .(Kuerzel, official_short, official_name, main_dimension)
  ]
)
message("")
message("Broad socio-economic dimensions:")
print(summary_broad)
message("")
message("Outputs written to: ", out_dir)
