bayes_variance_componentsnew <-  function(var_u_mean_omics_list,
                                       se_var_u_omics_list,
                                       var_u_total,
                                       var_residual,
                                       se_var_residual){

  genomic_h2 <- (var_u_total)/(var_u_total+var_residual)
  se_genomic_h2 <- standard_deviation(genomic_h2)
  genomic_h2 <- mean(genomic_h2)
  var_residual <- mean(var_residual)
  se_var_u_total = standard_deviation(var_u_total)
  if(length(var_u_mean_omics_list)>1){
variance_components <- data.frame(Components = c(unlist(var_u_mean_omics_list), mean(var_u_total), var_residual, genomic_h2),
                                  Standard_error = c(unlist(se_var_u_omics_list), se_var_u_total, se_var_residual, se_genomic_h2),
                                  row.names = c(paste("genetic_variance", seq_along(var_u_mean_omics_list), sep = "_"), "total_genetic_variance",
                                                "residual_variance", "heritability"),
                                  stringsAsFactors = FALSE
)
  }

  if(length(var_u_mean_omics_list)==1){
    variance_components <- data.frame(Components = c(unlist(var_u_mean_omics_list), var_residual, genomic_h2),
                                      Standard_error = c(unlist(se_var_u_omics_list), se_var_residual, se_genomic_h2),
                                      row.names = c("genetic_variance",
                                                    "residual_variance", "heritability"),
                                      stringsAsFactors = FALSE
    )
  }

return(variance_components)

}
