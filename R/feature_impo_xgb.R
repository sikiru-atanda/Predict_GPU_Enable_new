### XGBoost Feature Importance

#' Extract top-ranked XGBoost feature importances
#'
#' Returns the top `N_feature_impo` features from a fitted XGBoost model ranked
#' by importance / gain, used for downstream feature selection.
#'
#' @param xgb_fit A trained XGBoost model object.
#' @param X_train Numeric matrix of genomic / omics features used to train
#'   `xgb_fit` (column names are the feature labels).
#' @param N_feature_impo Number of top features to return ranked by importance
#'   / weight.
#' @param xgb_booster The XGBoost booster type used (`"gbtree"`, `"gblinear"` or
#'   `"dart"`); selects the appropriate importance branch.
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return A data frame of the top features and their importance scores.
#' @export
#'
#' @examples
feature_impo_xgb <- function(xgb_fit = NULL,
                             X_train = NULL,
                             N_feature_impo = NULL,
                             xgb_booster = NULL,
                             ...){

  if (xgb_booster == "gblinear"){
importance_matrix <- xgboost::xgb.importance(
  feature_names = colnames(X_train),
  model = xgb_fit
)
#importance_matrix

#names(importance_matrix) <- c("Feature", "Importance")
# Sort the features based on absolute coefficient values in descending order
importance_matrix <- importance_matrix[order(-importance_matrix$Weight), ]
if(!is.null(N_feature_impo)){
  importance_matrix <- importance_matrix[1:N_feature_impo, ]

}

  } else {

    if (xgb_booster == "gbtree"){
      importance_matrix <-  xgboost::xgb.importance(feature_names = colnames(X_train), model = xgb_fit)

      importance_matrix <- importance_matrix[order(-importance_matrix$Gain), ]
      if(!is.null(N_feature_impo)){
        importance_matrix <- importance_matrix[1:N_feature_impo, ]

      }

    }
}

#output <- list(feature_weight = importance_matrix)

return(importance_matrix)

}

#xgboost::xgb.plot.importance(importance_matrix)


# if(nrow(importance_matrix)>N_feature_impo){
#
#    ggplot2::ggplot(importance_matrix[1:N_feature_impo, ],  ggplot2::aes(x=reorder(Feature, Importance),
#                                                          y=Importance,
#                                                          color=as.factor(Feature)))+
#      ggplot2::geom_point() +
#      ggplot2::geom_segment( ggplot2::aes(x=Feature,
#                       xend=Feature,
#                       y=0,
#                       yend=Importance)) +
#      ggplot2::scale_color_discrete(name="VariableName") +
#      ggplot2::ylab("Importance") +
#      ggplot2::xlab("VariableName") +
#      ggplot2::coord_flip()
#
#
#
#
# } else {
#
#  ggplot2::ggplot(importance_matrix,  ggplot2::aes(x=reorder(Feature, Importance),
#                                                   y=Importance,
#                                                   color=as.factor(Feature))) +
#
#      ggplot2::geom_point() +
#      ggplot2::geom_segment(ggplot2::aes(x=Feature,
#                       xend=Feature,
#                       y=0,
#                       yend=Importance)) +
#      ggplot2::scale_color_discrete(name="VariableName") +
#      ggplot2::ylab("Importance") +
#      ggplot2::xlab("VariableName") +
#      ggplot2::coord_flip()
#
#
# }
#
# img <-   ggplot2::ggsave("feature_impo_xgb.jpeg",width = 10, height = 6, units = "in", dpi = 300)
#
#
# output <- list(importance_matrix, img)
#
# names(output) <- c('feature_weight', 'plot_top10_Impo_feature')
#
#  return(output)
#
#
# }
