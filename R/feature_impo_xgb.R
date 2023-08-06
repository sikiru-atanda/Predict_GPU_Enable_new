### XGBoost Feature Importance

feature_impo_xgb <- function(xgb_fit,
                             X_train,
                             ...){

importance_matrix <- xgboost::xgb.importance(
  feature_names = colnames(X_train),
  model = xgb_fit
)
importance_matrix

names(importance_matrix) <- c("Feature", "Importance")

#xgboost::xgb.plot.importance(importance_matrix)

if(nrow(importance_matrix)>10){

  (p <-   ggplot2::ggplot(importance_matrix[1:10, ],  ggplot2::aes(x=reorder(Feature, Importance),
                                                         y=Importance,
                                                         color=as.factor(Feature)))+
     ggplot2::geom_point() +
     ggplot2::geom_segment( ggplot2::aes(x=Feature,
                      xend=Feature,
                      y=0,
                      yend=Importance)) +
     ggplot2::scale_color_discrete(name="VariableName") +
     ggplot2::ylab("Importance") +
     ggplot2::xlab("VariableName") +
     ggplot2::coord_flip()

  )


} else {

  (p <-    ggplot2::ggplot(importance_matrix,  ggplot2::aes(x=reorder(Feature, Importance),
                                                  y=Importance,
                                                  color=as.factor(Feature))) +

     ggplot2::geom_point() +
     ggplot2::geom_segment(ggplot2::aes(x=Feature,
                      xend=Feature,
                      y=0,
                      yend=Importance)) +
     ggplot2::scale_color_discrete(name="VariableName") +
     ggplot2::ylab("Importance") +
     ggplot2::xlab("VariableName") +
     ggplot2::coord_flip()

  )

}

img <-   ggplot2::ggsave("feature_impo_xgb.jpeg",width = 10, height = 6, units = "in", dpi = 300)


output <- list(importance_matrix, img)

names(output) <- c('feature_weight', 'plot_top10_Impo_feature')

 return(output)


}
