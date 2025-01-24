compute_model_agreement <- function(
    preds_matrix, 
    disagreement_threshold = NULL, 
    threshold_quantile = 0.90
) {
  # preds_matrix: matrix of predictions from different models
  #     rows = test instances, columns = models
  # disagreement_threshold: if NULL, define by threshold quantile of disagreement
  # threshold_quantile: e.g., 0.90 means if disagreement is above the 90th percentile, flag untrustworthy
  
  # 1. Compute standard deviation across models for each test instance
  disagreement <- apply(preds_matrix, 1, sd)
  
  # 2. If no threshold is specified, pick one based on a quantile
  if (is.null(disagreement_threshold)) {
    disagreement_threshold <- quantile(disagreement, threshold_quantile)
  }
  
  # 3. Flag test instances as untrusted if disagreement > threshold
  untrusted_flags <- as.numeric(disagreement > disagreement_threshold)
  
  return(list(
    disagreement           = disagreement,
    disagreement_threshold = disagreement_threshold,
    untrusted_flags        = untrusted_flags
  ))
}

### Example usage:
# Suppose you have 3 models, each providing predictions for N test samples.
# preds_matrix might look like:
#     pred_model1 <- c( ... )  # length N
#     pred_model2 <- c( ... )
#     pred_model3 <- c( ... )
# preds_matrix <- cbind(pred_model1, pred_model2, pred_model3)
#
# agreement_result <- compute_model_agreement(preds_matrix, threshold_quantile = 0.90)
# disagreement_vals  <- agreement_result$disagreement
# disagree_threshold <- agreement_result$disagreement_threshold
# untrusted_flags    <- agreement_result$untrusted_flags  # 1 = "too high disagreement"
