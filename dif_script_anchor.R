suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(TAM)
})

prepare_resp_for_rasch <- function(resp, gender_vec) {
  keep_gender <- !is.na(gender_vec) & gender_vec %in% c("M", "F")
  resp <- resp[keep_gender, , drop = FALSE]
  gender_vec <- gender_vec[keep_gender]

  keep_persons <- rowSums(!is.na(resp)) > 0
  resp <- resp[keep_persons, , drop = FALSE]
  gender_vec <- gender_vec[keep_persons]

  keep_items <- apply(resp, 2, function(x) {
    x <- x[!is.na(x)]
    length(x) > 0 && length(unique(x)) > 1 && sum(x > 0, na.rm = TRUE) > 0
  })
  resp <- resp[, keep_items, drop = FALSE]

  keep_persons2 <- rowSums(!is.na(resp)) > 0
  resp <- resp[keep_persons2, , drop = FALSE]
  gender_vec <- gender_vec[keep_persons2]

  list(resp = as.data.frame(resp), gender = gender_vec)
}

find_true_difficulty_col <- function(truth_df) {
  candidates <- c("true_difficulty", "b_base")
  hit <- candidates[candidates %in% names(truth_df)]

  if (length(hit) == 0) {
    stop(
      paste0(
        "Truth file must contain a true item difficulty column for anchoring. Accepted names: ",
        paste(candidates, collapse = ", "),
        "."
      )
    )
  }

  hit[[1]]
}

choose_anchor_items <- function(truth_df, n_anchors = 10L, difficulty_col = NULL) {
  if (is.null(difficulty_col)) {
    difficulty_col <- find_true_difficulty_col(truth_df)
  }

  pool <- truth_df %>%
    filter(!is.na(item_id), !is.na(.data[[difficulty_col]]), !dif_flag) %>%
    arrange(.data[[difficulty_col]], item_id)

  if (nrow(pool) < n_anchors) {
    stop(sprintf(
      "Requested %d anchors but only %d eligible non-DIF items are available.",
      n_anchors,
      nrow(pool)
    ))
  }

  pick_idx <- unique(round(seq(1, nrow(pool), length.out = n_anchors)))
  anchors <- pool[pick_idx, , drop = FALSE]

  if (nrow(anchors) < n_anchors) {
    remainder <- pool %>% anti_join(anchors, by = "item_id")
    anchors <- bind_rows(anchors, head(remainder, n_anchors - nrow(anchors)))
  }

  anchors %>%
    transmute(
      item_id = item_id,
      anchor_value = .data[[difficulty_col]]
    ) %>%
    distinct(item_id, .keep_all = TRUE)
}

build_xsi_fixed <- function(resp_df, anchor_df) {
  pos <- match(anchor_df$item_id, colnames(resp_df))

  if (any(is.na(pos))) {
    stop(sprintf(
      "These anchor items are missing from the response matrix: %s",
      paste(anchor_df$item_id[is.na(pos)], collapse = ", ")
    ))
  }

  cbind(pos, anchor_df$anchor_value)
}

fit_rasch_model <- function(resp, anchor_df = NULL, adj = 0.5, bias = FALSE) {
  if (is.null(anchor_df) || nrow(anchor_df) == 0) {
    return(TAM::tam.jml(resp = resp, bias = bias))
  }

  TAM::tam.jml(
    resp = resp,
    adj = adj,
    bias = bias,
    xsi.fixed = build_xsi_fixed(resp, anchor_df),
    constraint = "cases"
  )
}

compute_dif_table <- function(
  resp,
  gender_vec,
  mod,
  lookup_map,
  attainment_lab,
  digits_difficulty = 2,
  digits_facility = 4,
  digits_contrast = 2
) {
  dff <- as.data.frame(resp)
  item_names <- colnames(dff)
  i <- as.numeric(mod$item$xsi.item)
  names(i) <- item_names
  p <- as.numeric(mod$theta)

  probs <- matrix(0, nrow = nrow(dff), ncol = ncol(dff), dimnames = list(NULL, item_names))
  for (m in seq_len(ncol(dff))) {
    for (n in seq_len(nrow(dff))) {
      probs[n, m] <- 1 / (1 + exp(i[m] - p[n]))
    }
  }

  group <- gender_vec
  M <- dff[group == "M", , drop = FALSE]
  Ff <- dff[group == "F", , drop = FALSE]
  Mthet <- p[group == "M"]
  Fthet <- p[group == "F"]

  Mprobs <- probs[group == "M", , drop = FALSE]
  Fprobs <- probs[group == "F", , drop = FALSE]

  sums <- rep(NA_real_, ncol(dff))
  Fsums <- rep(NA_real_, ncol(dff))

  for (ij in seq_len(ncol(dff))) {
    Mm <- M[[ij]]
    mask <- !is.na(Mm)
    Mp <- Mprobs[mask, ij]
    denom <- sum(Mp * (1 - Mp))
    sums[ij] <- if (denom > 0) sum(Mm[mask] - Mp) / denom else NA_real_
  }

  for (ij in seq_len(ncol(dff))) {
    Fm <- Ff[[ij]]
    mask <- !is.na(Fm)
    Fp <- Fprobs[mask, ij]
    denom <- sum(Fp * (1 - Fp))
    Fsums[ij] <- if (denom > 0) sum(Fm[mask] - Fp) / denom else NA_real_
  }

  jj <- i - sums
  ff <- i - Fsums

  newprobsmale <- matrix(NA_real_, nrow = nrow(M), ncol = ncol(dff))
  for (m in seq_len(ncol(dff))) {
    for (n in seq_len(nrow(M))) {
      if (!is.na(M[n, m])) {
        newprobsmale[n, m] <- 1 / (1 + exp(jj[m] - Mthet[n]))
      }
    }
  }

  newprobsfemale <- matrix(NA_real_, nrow = nrow(Ff), ncol = ncol(dff))
  for (m in seq_len(ncol(dff))) {
    for (n in seq_len(nrow(Ff))) {
      if (!is.na(Ff[n, m])) {
        newprobsfemale[n, m] <- 1 / (1 + exp(ff[m] - Fthet[n]))
      }
    }
  }

  se_nongender <- 1 / sqrt(colSums(probs * (1 - probs), na.rm = TRUE))
  se_M <- 1 / sqrt(colSums(newprobsmale * (1 - newprobsmale), na.rm = TRUE))
  se_F <- 1 / sqrt(colSums(newprobsfemale * (1 - newprobsfemale), na.rm = TRUE))
  tdist <- (Fsums - sums) / sqrt(se_M^2 + se_F^2)

  n_m <- colSums(!is.na(M))
  n_f <- colSums(!is.na(Ff))

  denom <- (
    (1 + 4 / (n_m + n_f - 2)) *
      ((se_M^4 / (n_m + 1)) + (se_F^4 / (n_f + 1)))
  )
  num <- (se_M^2 + se_F^2)^2
  dfs <- num / denom

  Prob <- rep(NA_real_, ncol(dff))
  favours <- rep(NA_character_, ncol(dff))
  stat_signif <- rep(NA_character_, ncol(dff))
  contrast_signif <- rep(NA_character_, ncol(dff))

  contrast_vals <- ff - jj
  abs_contrast <- abs(contrast_vals)

  for (k in seq_len(ncol(dff))) {
    if (is.finite(tdist[k]) && is.finite(dfs[k]) && dfs[k] > 0) {
      Prob[k] <- pt(abs(tdist[k]), as.numeric(dfs[k]), lower.tail = FALSE) * 2
    }

    if (!is.na(Prob[k]) && Prob[k] <= 0.05) {
      favours[k] <- if (jj[k] < ff[k]) "Boys" else "Girls"
      stat_signif[k] <- "*"
    }
    if (!is.na(Prob[k]) && Prob[k] <= 0.01) stat_signif[k] <- "**"
    if (!is.na(Prob[k]) && Prob[k] <= 0.001) stat_signif[k] <- "***"

    if (!is.na(abs_contrast[k])) {
      if (abs_contrast[k] < 0.43) contrast_signif[k] <- "*"
      if (abs_contrast[k] >= 0.43 && abs_contrast[k] < 0.64) contrast_signif[k] <- "**"
      if (abs_contrast[k] >= 0.64) contrast_signif[k] <- "***"
    }
  }

  all_b <- c(i, jj, ff)
  theta_min <- min(all_b, na.rm = TRUE)
  theta_max <- max(all_b, na.rm = TRUE)
  scale_range <- theta_max - theta_min

  make_facility <- function(x) {
    if (!is.finite(scale_range) || scale_range <= 0) {
      return(rep(NA_real_, length(x)))
    }
    100 - (((x - theta_min) / scale_range) * 100)
  }

  meta <- lookup_map %>%
    transmute(
      item_id = item_id,
      item_type = item_type,
      split_type = split_type,
      fit_group = fit_group
    )

  tibble(
    item_id = item_names,
    attainment = attainment_lab,
    difficulty = round(i, digits_difficulty),
    difficulty_M = round(jj, digits_difficulty),
    difficulty_F = round(ff, digits_difficulty),
    Facility = round(make_facility(i), digits_facility),
    Facility_M = round(make_facility(jj), digits_facility),
    Facility_F = round(make_facility(ff), digits_facility),
    contrast = round(contrast_vals, digits_contrast),
    se_nongender = se_nongender,
    se_M = se_M,
    se_F = se_F,
    t_value = tdist,
    df = dfs,
    Prob = Prob,
    favours = favours,
    `Stat Signif` = stat_signif,
    `Contrast Signif` = contrast_signif
  ) %>%
    left_join(meta, by = "item_id") %>%
    relocate(
      fit_group,
      attainment,
      item_id,
      item_type,
      split_type,
      difficulty,
      difficulty_M,
      difficulty_F,
      Facility,
      Facility_M,
      Facility_F,
      contrast,
      `Contrast Signif`,
      se_nongender,
      se_M,
      se_F,
      t_value,
      df,
      Prob,
      `Stat Signif`,
      favours
    )
}

run_dif_from_item_response_object <- function(
  irt_obj,
  attainment_labels = NULL,
  outdir = NULL,
  file_stub = NULL,
  digits_difficulty = 2,
  digits_facility = 4,
  digits_contrast = 2,
  write_split_csvs = FALSE,
  write_combined_csv = FALSE,
  anchor_df = NULL,
  anchor_adj = 0.5
) {
  valid_groups <- irt_obj$lookup_map %>%
    filter(!is.na(fit_group)) %>%
    distinct(fit_group) %>%
    pull(fit_group)

  out_list <- vector("list", length(valid_groups))
  names(out_list) <- valid_groups

  for (ii in seq_along(valid_groups)) {
    sp <- valid_groups[ii]

    cols_sp <- irt_obj$lookup_map %>%
      filter(fit_group == sp) %>%
      arrange(output_order) %>%
      pull(item_col)

    items_sp <- irt_obj$items %>% dplyr::select(dplyr::all_of(cols_sp))
    prep <- prepare_resp_for_rasch(items_sp, irt_obj$gender)

    anchor_sp <- NULL
    if (!is.null(anchor_df)) {
      anchor_sp <- anchor_df %>% filter(item_id %in% colnames(prep$resp))
      if (nrow(anchor_sp) == 0) {
        stop(sprintf("No anchors available in fit group '%s' after preprocessing.", sp))
      }
    }

    mod <- fit_rasch_model(prep$resp, anchor_df = anchor_sp, adj = anchor_adj, bias = FALSE)
    attainment_lab <- if (is.null(attainment_labels)) sp else unname(attainment_labels[sp])

    out_list[[ii]] <- compute_dif_table(
      resp = prep$resp,
      gender_vec = prep$gender,
      mod = mod,
      lookup_map = irt_obj$lookup_map %>% filter(fit_group == sp),
      attainment_lab = attainment_lab,
      digits_difficulty = digits_difficulty,
      digits_facility = digits_facility,
      digits_contrast = digits_contrast
    )
  }

  bind_rows(out_list)
}

summarise_effect_error <- function(dif_df, truth_df) {
  dif_df %>%
    inner_join(
      truth_df %>%
        filter(dif_flag) %>%
        select(item_id, delta_group1, dif_type),
      by = "item_id"
    ) %>%
    mutate(
      error = contrast - delta_group1,
      abs_error = abs(error),
      sq_error = error^2
    ) %>%
    summarise(
      n_dif_items = n(),
      mean_contrast = mean(contrast, na.rm = TRUE),
      mean_true_delta = mean(delta_group1, na.rm = TRUE),
      bias = mean(error, na.rm = TRUE),
      mae = mean(abs_error, na.rm = TRUE),
      rmse = sqrt(mean(sq_error, na.rm = TRUE)),
      max_abs_error = max(abs_error, na.rm = TRUE),
      prop_within_0_01 = mean(abs_error <= 0.01, na.rm = TRUE),
      prop_within_0_02 = mean(abs_error <= 0.02, na.rm = TRUE),
      prop_within_0_05 = mean(abs_error <= 0.05, na.rm = TRUE)
    )
}

summarise_effect_error_by_type <- function(dif_df, truth_df) {
  dif_df %>%
    inner_join(
      truth_df %>%
        filter(dif_flag) %>%
        select(item_id, delta_group1, dif_type),
      by = "item_id"
    ) %>%
    mutate(
      error = contrast - delta_group1,
      abs_error = abs(error),
      sq_error = error^2
    ) %>%
    group_by(dif_type) %>%
    summarise(
      n_dif_items = n(),
      mean_contrast = mean(contrast, na.rm = TRUE),
      mean_true_delta = mean(delta_group1, na.rm = TRUE),
      bias = mean(error, na.rm = TRUE),
      mae = mean(abs_error, na.rm = TRUE),
      rmse = sqrt(mean(sq_error, na.rm = TRUE)),
      max_abs_error = max(abs_error, na.rm = TRUE),
      prop_within_0_01 = mean(abs_error <= 0.01, na.rm = TRUE),
      prop_within_0_02 = mean(abs_error <= 0.02, na.rm = TRUE),
      prop_within_0_05 = mean(abs_error <= 0.05, na.rm = TRUE),
      .groups = "drop"
    )
}
