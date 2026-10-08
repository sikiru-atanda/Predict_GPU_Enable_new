#' Summarize Hybrid Prediction Outputs
#'
#' @param predicted_object Hybrid prediction table returned by a hybrid model path.
#' @param response Response column name.
#' @param female_parent Female parent column name.
#' @param male_parent Male parent column name.
#' @param mode Hybrid mode label such as `hybrid_asreml`, `hybrid_bayes`,
#'   `hybrid_gp`, `hybrid_ml`, or `hybrid_dl`.
#' @param model_type Model name.
#'
#' @return A list of hybrid summary data frames.
#' @export
summary_statistics_hybrid <- function(predicted_object,
                                      response,
                                      female_parent,
                                      male_parent,
                                      mode,
                                      model_type) {
  pred_df <- as.data.frame(predicted_object, stringsAsFactors = FALSE)
  if (!nrow(pred_df)) {
    return(list(
      summary_statistics = data.frame(stat = character(), summary = character(), stringsAsFactors = FALSE),
      hybrid_component_summary = data.frame(),
      hybrid_train_test_summary = data.frame()
    ))
  }

  train_idx <- pred_df$Train_Test_Label != "Test"
  test_idx <- pred_df$Train_Test_Label == "Test"

  component_basis <- "none"
  component_cols <- NULL
  component_labels <- NULL
  if (all(c("Female_GCA", "Male_GCA", "SCA_effect") %in% names(pred_df))) {
    component_basis <- "model_native_gca_sca"
    component_cols <- c("Female_GCA", "Male_GCA", "SCA_effect")
    component_labels <- c("Female_GCA", "Male_GCA", "SCA_effect")
  } else if (all(c(
    "Female_additive_contribution",
    "Male_additive_contribution",
    "Hybrid_interaction_contribution"
  ) %in% names(pred_df))) {
    component_basis <- unique(pred_df$Predictive_decomposition_basis)
    component_basis <- component_basis[!is.na(component_basis) & nzchar(component_basis)]
    component_basis <- if (length(component_basis)) paste(component_basis, collapse = "; ") else "predictive_decomposition"
    component_cols <- c(
      "Female_additive_contribution",
      "Male_additive_contribution",
      "Hybrid_interaction_contribution"
    )
    component_labels <- c(
      "Female_additive_contribution",
      "Male_additive_contribution",
      "Hybrid_interaction_contribution"
    )
  }

  summary_rows <- data.frame(
    stat = c(
      "Mode",
      "Response_Family",
      "Model_Type",
      "Response",
      "Observed_Training_Hybrids",
      "Predicted_Test_Hybrids",
      "Total_Hybrids",
      "Female_Parent_Column",
      "Male_Parent_Column",
      "Female_Parent_Count",
      "Male_Parent_Count",
      "Hybrid_Component_Basis",
      "Train_Mean_Observed",
      "Train_Mean_Predicted",
      "Test_Mean_Predicted"
    ),
    summary = c(
      mode,
      "gaussian",
      model_type,
      response,
      sum(train_idx, na.rm = TRUE),
      sum(test_idx, na.rm = TRUE),
      nrow(pred_df),
      female_parent,
      male_parent,
      length(unique(as.character(pred_df[[female_parent]]))),
      length(unique(as.character(pred_df[[male_parent]]))),
      component_basis,
      sprintf("%.6f", mean(pred_df$Observed_value[train_idx], na.rm = TRUE)),
      sprintf("%.6f", mean(pred_df$Predicted_value[train_idx], na.rm = TRUE)),
      sprintf("%.6f", mean(pred_df$Predicted_value[test_idx], na.rm = TRUE))
    ),
    stringsAsFactors = FALSE
  )

  if (all(c("Observed_value", "Predicted_value") %in% names(pred_df))) {
    hybrid_train_test_summary <- do.call(
      rbind,
      lapply(c("Train", "Test"), function(lbl) {
        idx <- pred_df$Train_Test_Label == lbl
        data.frame(
          Train_Test_Label = lbl,
          n_hybrids = sum(idx, na.rm = TRUE),
          mean_observed = mean(pred_df$Observed_value[idx], na.rm = TRUE),
          mean_predicted = mean(pred_df$Predicted_value[idx], na.rm = TRUE),
          mean_absolute_error = mean(abs(pred_df$Observed_value[idx] - pred_df$Predicted_value[idx]), na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      })
    )
  } else {
    hybrid_train_test_summary <- data.frame()
  }

  if (!is.null(component_cols)) {
    hybrid_component_summary <- do.call(
      rbind,
      lapply(seq_along(component_cols), function(i) {
        vals <- as.numeric(pred_df[[component_cols[[i]]]])
        data.frame(
          component = component_labels[[i]],
          mean = mean(vals, na.rm = TRUE),
          mean_abs = mean(abs(vals), na.rm = TRUE),
          sd = stats::sd(vals, na.rm = TRUE),
          min = min(vals, na.rm = TRUE),
          max = max(vals, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      })
    )
  } else {
    hybrid_component_summary <- data.frame()
  }

  list(
    summary_statistics = summary_rows,
    hybrid_component_summary = hybrid_component_summary,
    hybrid_train_test_summary = hybrid_train_test_summary
  )
}
