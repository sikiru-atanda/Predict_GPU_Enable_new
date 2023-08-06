
predict_proR <- function(object_coeff1 = NULL,
                      object_coeff2 = NULL,
                      object_coeff3 = NULL,
                      object_coeff4 = NULL,
                      mu = NULL,
                      object_geno = NULL,
                      object_omic1 = NULL,
                      object_omic2 = NULL,
                      object_omic3 = NULL,
                      message = TRUE,
                      ...) {
  if((!is.null('object_coeff1') & (is.null('object_coeff2') & (is.null('object_coeff3') & is.null('object_coeff4'))))){

    if((!is.null('object_geno') & (is.null('object_omic1') & (is.null('object_omic2') & is.null('object_omic3'))))){

      pred <- object_geno %*%object_coeff1 + mu
    }

  }


  if((is.null('object_coeff1') & (!is.null('object_coeff2') & (is.null('object_coeff3') & is.null('object_coeff4'))))){

    if((!is.null('object_geno') & (!is.null('object_omic1') & (is.null('object_omic2') & is.null('object_omic3'))))){

      pred <- object_omic1 %*%object_coeff2  + mu
    }

  }

  if((is.null('object_coeff1') & (is.null('object_coeff2') & (!is.null('object_coeff3') & is.null('object_coeff4'))))){

    if((!is.null('object_geno') & (is.null('object_omic1') & (!is.null('object_omic2') & is.null('object_omic3'))))){

      pred <- object_omic2 %*%object_coeff3 + mu

    }

  }


  if((is.null('object_coeff1') & (is.null('object_coeff2') & (is.null('object_coeff3') & !is.null('object_coeff4'))))){

    if((!is.null('object_geno') & (is.null('object_omic1') & (is.null('object_omic2') & !is.null('object_omic3'))))){

      pred <- object_omic3 %*%object_coeff4 + mu
    }

  }

  if((is.null('object_coeff1') & (is.null('object_coeff2') & (is.null('object_coeff3') & !is.null('object_coeff4'))))){

    if((!is.null('object_geno') & (is.null('object_omic1') & (is.null('object_omic2') & !is.null('object_omic3'))))){

      pred <- object_omic3 %*%object_coeff4  + mu
    }

  }

  if((!is.null('object_coeff1') & (!is.null('object_coeff2') & (is.null('object_coeff3') & is.null('object_coeff4'))))){

    if((!is.null('object_geno') & (!is.null('object_omic1') & (is.null('object_omic2') & is.null('object_omic3'))))){

      pred <- object_geno %*%object_coeff1  + object_omic1 %*%object_coeff2 + mu



      }

  }


  if((!is.null('object_coeff1') & (is.null('object_coeff2') & (!is.null('object_coeff3') & is.null('object_coeff4'))))){

    if((!is.null('object_geno') & (is.null('object_omic1') & (!is.null('object_omic2') & is.null('object_omic3'))))){

      pred <- object_geno %*%object_coeff1  + object_omic2 %*%object_coeff3 + mu



    }

  }


  if((!is.null('object_coeff1') & (is.null('object_coeff2') & (is.null('object_coeff3') & !is.null('object_coeff4'))))){

    if((!is.null('object_geno') & (is.null('object_omic1') & (is.null('object_omic2') & !is.null('object_omic3'))))){

      pred <- object_geno %*%object_coeff1  + object_omic3 %*%object_coeff4 + mu


    }

  }


  if((is.null('object_coeff1') & (!is.null('object_coeff2') & (!is.null('object_coeff3') & is.null('object_coeff4'))))){

    if((is.null('object_geno') & (!is.null('object_omic1') & (!is.null('object_omic2') & !is.null('object_omic3'))))){

      pred <- object_omic1 %*%object_coeff2  + object_omic2 %*%object_coeff3 + mu


    }

  }


  if((is.null('object_coeff1') & (!is.null('object_coeff2') & (is.null('object_coeff3') & !is.null('object_coeff4'))))){

    if((is.null('object_geno') & (!is.null('object_omic1') & (is.null('object_omic2') & !is.null('object_omic3'))))){

      pred <- object_omic1 %*%object_coeff2  + object_omic3 %*%object_coeff4 + mu


    }

  }


  if((is.null('object_coeff1') & (is.null('object_coeff2') & (!is.null('object_coeff3') & !is.null('object_coeff4'))))){

    if((is.null('object_geno') & (is.null('object_omic1') & (!is.null('object_omic2') & !is.null('object_omic3'))))){

      pred <- object_omic2 %*%object_coeff3  + object_omic3 %*%object_coeff4 + mu


    }

  }


  if((is.null('object_coeff1') & (!is.null('object_coeff2') & (!is.null('object_coeff3') & !is.null('object_coeff4'))))){

    if((is.null('object_geno') & (!is.null('object_omic1') & (!is.null('object_omic2') & !is.null('object_omic3'))))){

      pred <- object_omic1 %*%object_coeff2  + object_omic2 %*%object_coeff3 + object_omic3 %*%object_coeff4 + mu


    }

  }


  if((!is.null('object_coeff1') & (!is.null('object_coeff2') & (!is.null('object_coeff3') & !is.null('object_coeff4'))))){

    if((is.null('object_geno') & (!is.null('object_omic1') & (!is.null('object_omic2') & !is.null('object_omic3'))))){

      pred <- object_geno%*%object_coeff1 + object_omic1 %*%object_coeff2  + object_omic2 %*%object_coeff3 + object_omic3 %*%object_coeff4 + mu


    }

  }

 return( pred)
}
