suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(TAM)
})

source("dif_script_anchor.R")

DATA_PATH <- "simulated_rasch_dif_responses.csv"
TRUTH_PATH <- "simulated_rasch_dif_truth.csv"
LATENT_PATH <- "simulated_rasch_dif_latent.csv"

OUTDIR <- "sample size effects"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

LOOKUP_SHEET <- "SIMULATED_RASCH_DIF"

TARGET_SAMPLE_SIZES <- c(
  seq(50, 300, by = 25),
  seq(500, 1500, by = 250),
  seq(2000, 5750, by = 250),
  7000,
  9000
)

SEED <- 67
CONVERGENCE_THRESHOLD <- 0.02
N_ANCHORS <- 10L
ANCHOR_ADJ <- 0.5

ATTAINMENT_LABELS <- c(
  "simulated_rasch" = "Simulated Rasch"
)

build_simulated_item_response_object <- function(data_path) {
  raw_df <- readr::read_csv(data_path, show_col_types = FALSE)

  if (!all(c("student_id", "group") %in% names(raw_df))) {
    stop("Simulated response file must contain columns: student_id, group.")
  }

  item_cols <- setdiff(names(raw_df), c("student_id", "group"))

  if (length(item_cols) == 0) {
    stop("No item columns found in simulated response file.")
  }

  gender_vec <- ifelse(raw_df$group == 1, "F", "M")

  lookup_map <- tibble(
    item_col = item_cols,
    item_id = item_cols,
    split_type = "SIMULATED",
    fit_group = "simulated_rasch",
    item_type = "SIMULATED",
    output_order = seq_along(item_cols)
  )

  list(
    raw = as.data.frame(raw_df),
    items = raw_df %>% dplyr::select(dplyr::all_of(item_cols)),
    gender = gender_vec,
    retained = raw_df %>% dplyr::select(student_id, group),
    lookup_map = lookup_map
  )
}

subset_irt_object <- function(irt_obj, idx) {
  out <- irt_obj
  out$raw <- irt_obj$raw[idx, , drop = FALSE]
  out$items <- irt_obj$items[idx, , drop = FALSE]
  out$gender <- irt_obj$gender[idx]

  if (!is.null(irt_obj$retained)) {
    out$retained <- irt_obj$retained[idx, , drop = FALSE]
  }

  out
}

sample_stratified_gender_idx <- function(gender_vec, n_total) {
  idx_m <- which(gender_vec == "M")
  idx_f <- which(gender_vec == "F")

  p_m <- length(idx_m) / (length(idx_m) + length(idx_f))
  n_m <- floor(n_total * p_m)
  n_f <- n_total - n_m

  if (n_m == 0) {
    n_m <- 1
    n_f <- n_total - 1
  }
  if (n_f == 0) {
    n_f <- 1
    n_m <- n_total - 1
  }

  sampled_m <- sample(idx_m, size = n_m, replace = FALSE)
  sampled_f <- sample(idx_f, size = n_f, replace = FALSE)

  sample(c(sampled_m, sampled_f), replace = FALSE)
}

get_wright_source_data_anchor <- function(irt_obj, n_students, rep, attainment_labels = NULL, anchor_df = NULL, anchor_adj = 0.5) {
  valid_groups <- irt_obj$lookup_map %>%
    filter(!is.na(fit_group)) %>%
    distinct(fit_group) %>%
    pull(fit_group)

  out_list <- vector("list", length(valid_groups))
  names(out_list) <- valid_groups

  for (i in seq_along(valid_groups)) {
    sp <- valid_groups[i]

    cols_sp <- irt_obj$lookup_map %>%
      filter(fit_group == sp) %>%
      arrange(output_order) %>%
      pull(item_col)

    items_sp <- irt_obj$items %>% dplyr::select(dplyr::all_of(cols_sp))
    prep <- prepare_resp_for_rasch(items_sp, irt_obj$gender)

    anchor_sp <- NULL
    if (!is.null(anchor_df)) {
      anchor_sp <- anchor_df %>% filter(item_id %in% colnames(prep$resp))
    }

    mod <- fit_rasch_model(prep$resp, anchor_df = anchor_sp, adj = anchor_adj, bias = FALSE)
    attainment_lab <- if (is.null(attainment_labels)) sp else unname(attainment_labels[sp])

    ability_df <- tibble(
      n_students = n_students,
      rep = rep,
      fit_group = sp,
      attainment = attainment_lab,
      curve = "ability",
      subgroup = ifelse(prep$gender == "M", "Boys", "Girls"),
      value = as.numeric(mod$theta)
    )

    item_df <- tibble(
      n_students = n_students,
      rep = rep,
      fit_group = sp,
      attainment = attainment_lab,
      curve = "item_difficulty",
      subgroup = "All",
      value = as.numeric(mod$item$xsi.item)
    )

    out_list[[i]] <- bind_rows(ability_df, item_df)
  }

  bind_rows(out_list)
}

get_n_reps <- function(n_now, n_total_full) {
  if (n_now == n_total_full) return(1L)
  if (n_now <= 300) return(200L)
  if (n_now <= 1500) return(100L)
  if (n_now <= 4000) return(50L)
  25L
}

set.seed(SEED)

irt_obj_full <- build_simulated_item_response_object(DATA_PATH)

usable_idx <- which(!is.na(irt_obj_full$gender) & irt_obj_full$gender %in% c("M", "F"))
irt_obj_full$raw <- irt_obj_full$raw[usable_idx, , drop = FALSE]
irt_obj_full$items <- irt_obj_full$items[usable_idx, , drop = FALSE]
irt_obj_full$gender <- irt_obj_full$gender[usable_idx]

if (!is.null(irt_obj_full$retained)) {
  irt_obj_full$retained <- irt_obj_full$retained[usable_idx, , drop = FALSE]
}

truth_df <- readr::read_csv(TRUTH_PATH, show_col_types = FALSE)

if (!all(c("item_id", "delta_group1", "dif_flag", "dif_type") %in% names(truth_df))) {
  stop("Truth file must contain item_id, delta_group1, dif_flag, dif_type.")
}

true_difficulty_col <- find_true_difficulty_col(truth_df)
anchor_df <- choose_anchor_items(truth_df, n_anchors = N_ANCHORS, difficulty_col = true_difficulty_col)

missing_anchor_items <- setdiff(anchor_df$item_id, colnames(irt_obj_full$items))
if (length(missing_anchor_items) > 0) {
  stop(sprintf(
    "These anchor items are not present in the response data: %s",
    paste(missing_anchor_items, collapse = ", ")
  ))
}

n_total_full <- length(irt_obj_full$gender)
sample_sizes <- sort(unique(c(TARGET_SAMPLE_SIZES, n_total_full)))
sample_sizes <- sample_sizes[sample_sizes <= n_total_full]

dif_results_list <- list()
item_level_dif_list <- list()
effect_type_list <- list()
wright_source_list <- list()

for (n_now in sample_sizes) {
  reps_now <- get_n_reps(n_now, n_total_full)

  for (rep_i in seq_len(reps_now)) {
    idx <- if (n_now == n_total_full) {
      seq_len(n_total_full)
    } else {
      sample_stratified_gender_idx(irt_obj_full$gender, n_now)
    }

    irt_sub <- subset_irt_object(irt_obj_full, idx)

    dif_out <- run_dif_from_item_response_object(
      irt_obj = irt_sub,
      attainment_labels = ATTAINMENT_LABELS,
      outdir = NULL,
      file_stub = NULL,
      digits_difficulty = 2,
      digits_facility = 4,
      digits_contrast = 2,
      write_split_csvs = FALSE,
      write_combined_csv = FALSE,
      anchor_df = anchor_df,
      anchor_adj = ANCHOR_ADJ
    )

    item_level_dif_list[[length(item_level_dif_list) + 1L]] <- dif_out %>%
      left_join(
        truth_df %>% select(item_id, delta_group1, dif_flag, dif_type),
        by = "item_id"
      ) %>%
      mutate(
        contrast_error = contrast - delta_group1,
        abs_contrast_error = abs(contrast_error),
        n_students = n_now,
        rep = rep_i
      )

    dif_results_list[[length(dif_results_list) + 1L]] <- summarise_effect_error(
      dif_df = dif_out,
      truth_df = truth_df
    ) %>%
      mutate(
        n_students = n_now,
        rep = rep_i
      ) %>%
      select(
        n_students, rep, n_dif_items, mean_contrast, mean_true_delta,
        bias, mae, rmse, max_abs_error,
        prop_within_0_01, prop_within_0_02, prop_within_0_05
      )

    effect_type_list[[length(effect_type_list) + 1L]] <- summarise_effect_error_by_type(
      dif_df = dif_out,
      truth_df = truth_df
    ) %>%
      mutate(
        n_students = n_now,
        rep = rep_i
      ) %>%
      select(
        n_students, rep, dif_type, n_dif_items, mean_contrast, mean_true_delta,
        bias, mae, rmse, max_abs_error,
        prop_within_0_01, prop_within_0_02, prop_within_0_05
      )

    wright_source_list[[length(wright_source_list) + 1L]] <- get_wright_source_data_anchor(
      irt_obj = irt_sub,
      n_students = n_now,
      rep = rep_i,
      attainment_labels = ATTAINMENT_LABELS,
      anchor_df = anchor_df,
      anchor_adj = ANCHOR_ADJ
    )
  }
}

dif_results_long <- bind_rows(dif_results_list)
effect_type_long <- bind_rows(effect_type_list)
item_level_dif <- bind_rows(item_level_dif_list)
wright_source_long <- bind_rows(wright_source_list)

write_csv(item_level_dif, file.path(OUTDIR, "dif_item_level_results.csv"))
write_csv(effect_type_long, file.path(OUTDIR, "dif_effect_convergence_by_type_raw.csv"))

dif_summary <- dif_results_long %>%
  group_by(n_students) %>%
  summarise(
    n_reps = n(),
    mean_bias = mean(bias, na.rm = TRUE),
    sd_bias = ifelse(n() > 1, sd(bias, na.rm = TRUE), 0),
    mean_mae = mean(mae, na.rm = TRUE),
    sd_mae = ifelse(n() > 1, sd(mae, na.rm = TRUE), 0),
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = ifelse(n() > 1, sd(rmse, na.rm = TRUE), 0),
    mean_max_abs_error = mean(max_abs_error, na.rm = TRUE),
    mean_prop_within_0_01 = mean(prop_within_0_01, na.rm = TRUE),
    mean_prop_within_0_02 = mean(prop_within_0_02, na.rm = TRUE),
    mean_prop_within_0_05 = mean(prop_within_0_05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ymin_mae = pmax(0, mean_mae - sd_mae),
    ymax_mae = mean_mae + sd_mae,
    ymin_rmse = pmax(0, mean_rmse - sd_rmse),
    ymax_rmse = mean_rmse + sd_rmse,
    log10_n_students = log10(n_students),
    converged_mae_0_01 = mean_mae <= CONVERGENCE_THRESHOLD,
    converged_rmse_0_01 = mean_rmse <= CONVERGENCE_THRESHOLD
  ) %>%
  arrange(n_students)

write_csv(dif_results_long, file.path(OUTDIR, "dif_effect_convergence_raw.csv"))
write_csv(dif_summary, file.path(OUTDIR, "dif_effect_convergence_summary.csv"))
write_csv(anchor_df, file.path(OUTDIR, "anchor_items_used.csv"))
write_csv(wright_source_long, file.path(OUTDIR, "wright_plot_source_values.csv"))

analysis_metadata <- tibble(
  data_path = DATA_PATH,
  truth_path = TRUTH_PATH,
  lookup_sheet = LOOKUP_SHEET,
  n_total_full = n_total_full,
  n_true_dif_items = sum(truth_df$dif_flag, na.rm = TRUE),
  seed = SEED,
  convergence_threshold = CONVERGENCE_THRESHOLD,
  n_anchors = N_ANCHORS,
  anchor_adj = ANCHOR_ADJ
)

write_csv(analysis_metadata, file.path(OUTDIR, "sample_size_effect_analysis_metadata.csv"))

print(anchor_df)
print(dif_summary)
cat("\nSaved analysis outputs to:", normalizePath(OUTDIR), "\n")
