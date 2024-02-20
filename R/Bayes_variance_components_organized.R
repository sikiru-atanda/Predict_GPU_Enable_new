

#' Title
#'
#' @param Var_U_1_Se
#' @param Var_U_2_Se
#' @param Var_U_3_Se
#' @param Var_U_4_Se
#' @param Var_U_Se
#' @param Var_E
#' @param genomic_h2
#' @param genomic_h2_Se
#' @param Var_U_1
#' @param Var_U_2
#' @param Var_U_3
#' @param Var_U_4
#' @param Var_U
#' @param Var_E_Se
#'
#' @return
#' @export
#'
#' @examples
#'
#'


Bayes_variance_components <- function(
                                Var_U_1 = NULL,
                                Var_U_2 = NULL,
                                Var_U_3 = NULL,
                                Var_U_4 = NULL,
                                Var_U = NULL,
                                Var_E = NULL,
                                genomic_h2= NULL,
                                Var_U_1_Se = NULL,
                                Var_U_2_Se = NULL,
                                Var_U_3_Se= NULL,
                                Var_U_4_Se = NULL,
                                Var_U_Se = NULL,
                                Var_E_Se = NULL,
                                genomic_h2_Se = NULL){

# browser()
#
#   cat(Var_U_1)
#   cat(Var_U_2)
#   cat(Var_U_3)
#   cat(Var_U_4)
#   cat(Var_U)
#   cat(Var_E)
#   cat(genomic_h2)
#   cat(Var_U_1_Se)
#   cat(Var_U_2_Se)
#   cat(Var_U_3_Se)
#   cat(Var_U_4_Se)
#   cat(Var_U_Se)
#   cat(Var_E_Se)
#   cat(genomic_h2_Se)

# if(isTRUE(all(!sapply(list("Var_U", "Var_E", "genomic_h2","Var_U_Se","Var_E_Se", "genomic_h2_Se"), is.null)))
#     & isFALSE(all(sapply(list("Var_U_1", "Var_U_2", "Var_U_3", "Var_U_4",
#                              "Var_U_1_Se", "Var_U_2_Se", "Var_U_3_Se", "Var_U_4_Se"), is.null)))){

  if(((((is.null(Var_U_1)& is.null(Var_U_2))&is.null(Var_U_3)) & ((is.null(Var_U_4)& !is.null(Var_U))& !is.null(genomic_h2))) &
      (((is.null(Var_U_1_Se)& is.null(Var_U_2_Se))&is.null(Var_U_3_Se)) & ((is.null(Var_U_4_Se)& !is.null(Var_U_Se))& !is.null(genomic_h2_Se))))){

  Variance_components <- data.frame(Components = c(Var_U, Var_E, genomic_h2),
                                    Standard_error = c(Var_U_Se, Var_E_Se, genomic_h2_Se),
                                    row.names = c("genetic_variance", "residual_variance", "heritability"),
                                    stringsAsFactors = FALSE
  )

  }

  if(((((!is.null(Var_U_1)& !is.null(Var_U_2))&is.null(Var_U_3)) & ((is.null(Var_U_4)& !is.null(Var_U))& !is.null(genomic_h2))) &
      (((!is.null(Var_U_1_Se)& !is.null(Var_U_2_Se))&is.null(Var_U_3_Se)) & ((is.null(Var_U_4_Se)& !is.null(Var_U_Se))& !is.null(genomic_h2_Se))))){


    Variance_components <- data.frame(Components = c(Var_U_1,
                                                     Var_U_2,
                                                     Var_U,
                                                     Var_E,
                                                     genomic_h2),
                                      Standard_error = c(Var_U_1_Se,
                                                         Var_U_2_Se,
                                                         Var_U_Se,
                                                         Var_E_Se,
                                                         genomic_h2_Se),
                                      row.names = c("genetic_variance1","genetic_variance2",
                                                    "total_genetic_variance",
                                                    "residual_variance", "heritability"),
                                      stringsAsFactors = FALSE
    )

  }

  #####
  if(((((!is.null(Var_U_1)& !is.null(Var_U_2))&!is.null(Var_U_3)) & ((is.null(Var_U_4)& !is.null(Var_U))& !is.null(genomic_h2))) &
      (((!is.null(Var_U_1_Se)& !is.null(Var_U_2_Se))&!is.null(Var_U_3_Se)) & ((is.null(Var_U_4_Se)& !is.null(Var_U_Se))& !is.null(genomic_h2_Se))))){


    Variance_components <- data.frame(Components = c(Var_U_1,
                                                     Var_U_2,
                                                     Var_U_3,
                                                     Var_U,
                                                     Var_E,
                                                     genomic_h2),
                                      Standard_error = c(Var_U_1_Se,
                                                         Var_U_2_Se,
                                                         Var_U_3_Se,
                                                         Var_U_Se,
                                                         Var_E_Se,
                                                         genomic_h2_Se),
                                      row.names = c("genetic_variance1","genetic_variance2",
                                                    "genetic_variance3","total_genetic_variance",
                                                    "residual_variance", "heritability"),
                                      stringsAsFactors = FALSE
    )

  }


  #####
  # if(((((!is.null(Var_U_1)& !is.null(Var_U_2))&!is.null(Var_U_3)) & ((!is.null(Var_U_4)& !is.null(Var_U))& !is.null(genomic_h2))) &
  #     (((!is.null(Var_U_1_Se)& !is.null(Var_U_2_Se))&!is.null(Var_U_3_Se)) & ((!is.null(Var_U_4_Se)& !is.null(Var_U_Se))& !is.null(genomic_h2_Se))))){
  if(((((!is.null(Var_U_1)& !is.null(Var_U_2))&!is.null(Var_U_3)) & ((!is.null(Var_U_4)& !is.null(Var_U))& !is.null(genomic_h2))) &
      (((!is.null(Var_U_1_Se)& !is.null(Var_U_2_Se))&!is.null(Var_U_3_Se)) & ((!is.null(Var_U_4_Se)& !is.null(Var_U_Se))& !is.null(genomic_h2_Se))))){


    Variance_components <- data.frame(Components = c(Var_U_1,
                                                     Var_U_2,
                                                     Var_U_3,
                                                     Var_U_4,
                                                     Var_U,
                                                     Var_E,
                                                     genomic_h2),
                                      Standard_error = c(Var_U_1_Se,
                                                         Var_U_2_Se,
                                                         Var_U_3_Se,
                                                         Var_U_4_Se,
                                                         Var_U_Se,
                                                         Var_E_Se,
                                                         genomic_h2_Se),
                                      row.names = c("genetic_variance1","genetic_variance2",
                                                    "genetic_variance3", "genetic_variance4",
                                                    "total_genetic_variance",
                                                    "residual_variance", "heritability"),
                                      stringsAsFactors = FALSE
    )

  }

  return(Variance_components)

}

