#' Title
#'
#' @param response
#' @param cova
#' @param fixed
#' @param random
#' @param heter_resid
#' @param heter_groups
#' @param VarCov_str
#' h@param weights
#' @param GS_model
#' @param core
#' @param weights
#' @param pheno_data
#' @param gmatrix
#' @param gkernel
#' @param omic1_kernel
#' @param omic2_kernel
#' @param omic3_kernel
#' @param workspace
#' @param maxit
#' @param gen_name
#' @param ...
#'
#' @return
#' @examples
#' @importFrom foreach %dopar%

asreml_utilisOLDD <- function(
    fixed = NULL,
    random = NULL,
    cova=NULL,
    GS_model = NULL,
    response = NULL,
    pheno_data = NULL,
    gmatrix = NULL,
    gkernel = NULL,
    omic1_kernel = NULL,
    omic2_kernel = NULL,
    omic3_kernel = NULL,
    gen_name = NULL,
    heter_groups = NULL,
    heter_resid = FALSE,
    VarCov_str = NULL,
    weights = NULL,
    core = NULL,
    workspace=1e08,
    #pworkspace= 1e06,
    maxit = 50,
    ...
) {

  msg <- sprintf("==================================================\n")


  ######

  if (!is.null(gmatrix)){
    G_inv <-  chol2inv(chol(gmatrix))

    rownames(G_inv) <- rownames(gmatrix)
    colnames(G_inv) <- colnames(gmatrix)
    # asreml requires the Gmatrix to be inverted in sparse format
    # G_inv <- solve(G_bend)
    # #check the attribute rowNames
    attr(gmatrix, "rowNames") <- rownames(gmatrix)
    attr(gmatrix, "colNames") <- colnames(gmatrix)
    attr(G_inv, "rowNames") <- rownames(gmatrix)
    attr(G_inv, "colNames") <- colnames(gmatrix)
    attr(G_inv, "INVERSE") <- TRUE

  }

  if (!is.null(gkernel)){
    GK_inv <-  chol2inv(chol(gkernel))

    rownames(GK_inv) <- rownames(gkernel)
    colnames(GK_inv) <- colnames(gkernel)
    attr(gkernel, "rowNames") <- rownames(gkernel)
    attr(gkernel, "colNames") <- colnames(gkernel)
    attr(GK_inv, "rowNames") <- rownames(gkernel)
    attr(GK_inv, "colNames") <- colnames(gkernel)
    attr(GK_inv, "INVERSE") <- TRUE

  }

  ####
  if (!is.null(omic1_kernel)){
    omic1_inv <<-  chol2inv(chol(omic1_kernel))

    rownames(omic1_inv) <- rownames(omic1_kernel)
    colnames(omic1_inv) <- colnames(omic1_kernel)
    attr(omic1_kernel, "rowNames") <- rownames(omic1_kernel)
    attr(omic1_kernel, "colNames") <- colnames(omic1_kernel)
    attr(omic1_inv, "rowNames") <- rownames(omic1_kernel)
    attr(omic1_inv, "colNames") <- colnames(omic1_kernel)
    attr(omic1_inv, "INVERSE") <- TRUE

  }
  #############
  if (!is.null(omic2_kernel)){
    omic2_inv <<-  chol2inv(chol(omic2_kernel))

    rownames(omic2_inv) <- rownames(omic2_kernel)
    colnames(omic2_inv) <- colnames(omic2_kernel)
    attr(omic2_kernel, "rowNames") <- rownames(omic2_kernel)
    attr(omic2_kernel, "colNames") <- colnames(omic2_kernel)
    attr(omic2_inv, "rowNames") <- rownames(omic2_kernel)
    attr(omic2_inv, "colNames") <- colnames(omic2_kernel)
    attr(omic2_inv, "INVERSE") <- TRUE

  }
  ####
  if (!is.null(omic3_kernel)){
    omic3_inv <<-  chol2inv(chol(omic3_kernel))

    rownames(omic3_inv) <- rownames(omic3_kernel)
    colnames(omic3_inv) <- colnames(omic3_kernel)
    attr(omic3_kernel, "rowNames") <- rownames(omic3_kernel)
    attr(omic3_kernel, "colNames") <- colnames(omic3_kernel)
    attr(omic3_inv, "rowNames") <- rownames(omic3_kernel)
    attr(omic3_inv, "colNames") <- colnames(omic3_kernel)
    attr(omic3_inv, "INVERSE") <- TRUE

  }
  #######

  if ((exists("G_inv") | exists("GK_inv")) & ((!exists("omic1_inv") & !exists("omic3_inv")) & !exists("omic3_inv"))){

    #G_copy = 1
    if(exists("G_inv")){

      G_list <<- list('G')

    } else {

      if(exists("GK_inv")){

        #G_list <<- list('GK')
        G_list <<- list('G')

      }

    }


  }

  if ((!exists("G_inv") | !exists("GK_inv")) & ((exists("omic1_inv") & !exists("omic3_inv")) & !exists("omic3_inv"))){

    #G_copy = 1

    G_list <<- list('omic1')
  }


  if ((!exists("G_inv") | !exists("GK_inv")) & ((!exists("omic1_inv") & exists("omic2_inv")) & !exists("omic3_inv"))){

    #G_copy = 1

    G_list <<- list('omic2')
  }

  if ((!exists("G_inv") | !exists("GK_inv")) & ((!exists("omic1_inv") & !exists("omic2_inv")) & exists("omic3_inv"))){

    #G_copy = 1

    G_list <<- list('omic3')
  }

  #####
  if ((exists("G_inv") | exists("GK_inv")) & ((exists("omic1_inv") & !exists("omic2_inv")) & !exists("omic3_inv"))){

    #G_copy = 2
    if(exists("G_inv")){
      G_list <<- list('G', 'omic1')

    } else {
      if(exists("GK_inv")){
        G_list <<- list('G', 'omic1')
        #G_list <<- list('GK', 'omic1')

      }

    }
  }

  if ((exists("G_inv") | exists("GK_inv")) & ((!exists("omic1_inv") & exists("omic2_inv")) & !exists("omic3_inv"))){

    #G_copy = 2

    if(exists("G_inv")){
      G_list <<- list('G', 'omic2')

    } else {
      if(exists("GK_inv")){
        #G_list <<- list('GK', 'omic2')
        G_list <<- list('G', 'omic2')
      }

    }

  }


  if ((exists("G_inv") | exists("GK_inv")) & ((!exists("omic1_inv") & !exists("omic2_inv")) & exists("omic3_inv"))){

    #G_copy = 2

    if(exists("G_inv")){
      G_list <<- list('G', 'omic3')

    } else {
      if(exists("GK_inv")){
        G_list <<- list('G', 'omic3')

        #G_list <<- list('GK', 'omic3')

      }

    }
  }

  if ((!exists("G_inv") | !exists("GK_inv")) & ((exists("omic1_inv") & exists("omic2_inv")) & !exists("omic3_inv"))){

    #G_copy = 2

    G_list <<- list('omic1', 'omic2')
  }

  if ((!exists("G_inv") | !exists("GK_inv")) & ((exists("omic1_inv") & !exists("omic2_inv")) & exists("omic3_inv"))){

    #G_copy = 2

    G_list <<- list('omic1', 'omic3')
  }

  if ((!exists("G_inv") | !exists("GK_inv")) & ((!exists("omic1_inv") & exists("omic2_inv")) & exists("omic3_inv"))){

    #G_copy = 2

    G_list <<- list('omic2', 'omic3')
  }

  if ((!exists("G_inv") | !exists("GK_inv")) & ((exists("omic1_inv") & exists("omic2_inv")) & exists("omic3_inv"))){

    #G_copy = 3

    G_list <<- list('omic1', 'omic2', 'omic3')
  }

  if ((exists("G_inv") | exists("GK_inv")) & ((exists("omic1_inv") & exists("omic2_inv")) & !exists("omic3_inv"))){

    #G_copy = 3
    if(exists("G_inv")){
      G_list <<- list('G', 'omic1', 'omic2')

    } else {
      if(exists("GK_inv")){
        G_list <<- list('G', 'omic1', 'omic2')
        #G_list <<- list('GK', 'omic1', 'omic2')

      }

    }


  }

  if ((exists("G_inv") | exists("GK_inv")) & ((exists("omic1_inv") & !exists("omic2_inv")) & exists("omic3_inv"))){

    #G_copy = 3
    if(exists("G_inv")){
      G_list <<- list('G', 'omic1', 'omic3')

    } else {
      if(exists("GK_inv")){
        G_list <<- list('G', 'omic1', 'omic3')
        #G_list <<- list('GK', 'omic1', 'omic3')

      }

    }


  }

  if ((exists("G_inv") | exists("GK_inv")) & ((!exists("omic1_inv") & exists("omic2_inv")) & exists("omic3_inv"))){

    #G_copy = 3
    if(exists("G_inv")){
      G_list <<- list('G', 'omic2', 'omic3')

    } else {
      if(exists("GK_inv")){
        G_list <<- list('G', 'omic2', 'omic3')
        #G_list <<- list('GK', 'omic2', 'omic3')
      }

    }

  }

  if ((exists("G_inv") | exists("GK_inv")) & ((exists("omic1_inv") & exists("omic2_inv")) & exists("omic3_inv"))){

    #G_copy = 4
    if(exists("G_inv")){
      G_list <<- list('G', 'omic1', 'omic2', 'omic3')

    } else {
      if(exists("GK_inv")){
        G_list <<- list('G', 'omic1', 'omic2', 'omic3')
        #G_list <<- list('GK', 'omic1', 'omic2', 'omic3')
      }

    }


  }


  ## Here pworkspace is not included because predict function is not done here
  asreml::asreml.options(trace=FALSE, workspace = workspace,
                         #pworkspace = pworkspace,
                         maxit = maxit)


  #### When the gen_name are present in more than one environment/location
  if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
    if (!is.null(heter_resid)) {
      if (is.null(heter_groups)) {
        stop(print(paste(msg,'No column of heterogeneous groups provided.')), call. = FALSE)
      } else {

        if(!heter_groups%in%colnames(pheno_data)) {stop(print(paste(msg,'heterogenous group provided did not match column names in pheno_data')), call. = FALSE)}

        if(!all(sapply(heter_groups, function(x, pheno_data) is.factor(pheno_data[,x]),  pheno_data))) {
          pheno_data[,heter_groups] <- as.factor(pheno_data[,heter_groups])
        }
      }

    }

  }

  # if (is.null(core)){
  #   cl = parallel::detectCores()
  #
  #   if (cl> 4){
  #     # Try in parallel
  #     cl <<- parallel::makeCluster(4)
  #   } else{
  #     cl <<- parallel::makeCluster(2)
  #   }
  #
  # } else {
  #   if(!is.null(core)){
  #
  #     cl <<- parallel::makeCluster(core)
  #   }
  # }
  #
  # doParallel::registerDoParallel(cl)

  #if(length(response)>1){
  #a= 0


  if (exists("G_inv")){.GlobalEnv$G_inv <- G_inv}
  if (exists("GK_inv")){.GlobalEnv$GK_inv <- GK_inv}
  if (exists("omic1_inv")){.GlobalEnv$omic1_inv <- omic1_inv}
  if (exists("omic2_inv")){.GlobalEnv$omic2_inv <- omic2_inv}
  if (exists("omic3_inv")){.GlobalEnv$omic3_inv <- omic3_inv}

  #a <- a + 1
  # Code Strings for all factors y= XB + UZ + e
  code.asr <- as.character()

  # colnames(pheno_data)[colnames(pheno_data) == trait] <-
  #   deparse(substitute(trait))

  code.asr[1] <- paste0(paste('asreml::asreml(fixed=', 'trait'),  '~1')
  code.asr[2] <- 'random=~'
  code.asr[3] <- 'residual=~'

  # Adding covariates (fixed)
  if (!is.null(cova)) {

    cova_term <- strsplit(as.character(cova[2]), split = "[+]")[[1]]

    if (length(cova_term )>1) {

      for (c in 1:length(cova_term )) {
        code.asr[1] <- paste(code.asr[1], cova_term[c], sep='+')
      }

    } else {
      code.asr[1] <- paste(code.asr[1], cova_term, sep='+')

      #code.asr[1] <- paste(code.asr[1], cova)

    }

  }

  # Adding fixed factors
  if (!is.null(fixed)) {

    fixed_term <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]

    if (length(all.vars(fixed))>1) {

      for (v in 1:length(all.vars(fixed))) {
        code.asr[1] <- paste(code.asr[1], all.vars(fixed)[v], sep='+')
      }

    } else {
      code.asr[1] <- paste(code.asr[1], all.vars(fixed), sep='+')

    }

  } else{

    if (is.null(fixed)) {
      fixed_term = NULL
      Check_heter.grp.Fixed = NULL
    }
  }


  # Adding random factors

  if (is.null(random)){ stop(print('provide random term'))}
  if (!is.null(random)){
    #if (length(all.vars(random))>1) {
    rand_term <- strsplit(as.character(random[2]), split = "[+]")[[1]] # random parts
    rand_term = gsub(" ", "", rand_term)

    # if (is.null(heter_groups) & is.null(VarCov_str)){
    #
    #   if(exists("G_inv")){
    #     if (dim(pheno_data)[1]> dim(G_inv)[1]){
    #       stop(print('If gen_name is the only random term, number of rows for the pheno_datatypic data must equal number of row in the GRM matrix'),
    #            call. = FALSE)
    #     }
    #
    #   }
    #
    #   if(exists("GK_inv")){
    #     if (dim(pheno_data)[1]> dim(GK_inv)[1]){
    #       stop(print('If gen_name is the only random term, number of rows for the pheno_datatypic data must equal number of row in the GRM matrix'),
    #            call. = FALSE)
    #     }
    #
    #   }
    # }

    #rand_InterPresent = grep(":", rand_term)
    ## No interaction term
    rand_term_No_Inter = rand_term[!grepl(":", rand_term)]

    ## Interaction term
    Check_rand_Inter = rand_term[grepl(":", rand_term)]

    ## Extract other side expect the gen_name
    rand_term_No_Inter_No_Gen = rand_term_No_Inter[!rand_term_No_Inter%in%gen_name]
    ## Get the position of other terms in the random that is not gen_name
    Non_Gen_pos = match(rand_term_No_Inter_No_Gen, rand_term)
    if(anyNA(Non_Gen_pos)){Non_Gen_pos = NULL}

    if(length(Check_rand_Inter)!=0){

      TestPresentofGeno = grep(gen_name, Check_rand_Inter, value = TRUE)
      if(anyNA(TestPresentofGeno)) {TestPresentofGeno = NULL}

      if(length(TestPresentofGeno)!=0){

        Inter_Gen_pos = match(TestPresentofGeno, rand_term)
        if(anyNA(Inter_Gen_pos)){Inter_Gen_pos = NULL}

        Non_Gen_Inter_Test =  Check_rand_Inter[!Check_rand_Inter%in%TestPresentofGeno]
        if(anyNA(Non_Gen_Inter_Test)){Non_Gen_Inter_Test = NULL}

      } else {

        Non_Gen_Inter_Test = NULL

        Inter_Gen_pos = NULL
      }

      if(length(Non_Gen_Inter_Test)!=0){

        Inter_Non_Gen_pos = match(Non_Gen_Inter_Test, rand_term)

      }

      if(is.null(Inter_Gen_pos)){warning(paste(gen_name, 'missing in the interaction term. This implies you cannot estimate GxE'))}

    } else {

      if(length(Check_rand_Inter)==0){

        Inter_Gen_pos = NULL
      }
    }

    Gen_pos=  match(gen_name, rand_term)
    if(anyNA(Gen_pos)){Gen_pos = NULL}

    ## Eg when GID:Env with no variance structure is the only term in the
    ## random effect provided by the user
    if(length(rand_term)==length(Check_rand_Inter) & is.null(VarCov_str)){
      stop(print(paste(msg, "provide variance-covariance structure")), call. = FALSE)

    }


    if (length(Check_rand_Inter)>=1){

      if (is.null(heter_groups)) {stop(print(paste(msg, "hetero.groups cannot be NULL")), call. = FALSE)}
    }


    if(!is.null(VarCov_str)){

      if(is.null(heter_groups)){stop(print(paste(msg, "Provide heter_groups to model specified variance-covariance structure")), call. = FALSE)}

      ### Check if the number of hetero.Grp is greater 5 or greater than 5
      NN = nlevels(pheno_data[, heter_groups])
      if (NN >=5 & isFALSE(grepl("fa", VarCov_str))){

        msg <- sprintf("\r==================================================\n")

        warning(paste(msg, "The number of", heter_groups, " is", NN,  "consider using factor analytic model"), immediate. = TRUE, call. =FALSE)
      }

      if (!is.null(fixed_term)){
        Check_heter.grp.Fixed  = match(heter_groups, fixed_term)
        if(anyNA(Check_heter.grp.Fixed)) {Check_heter.grp.Fixed = NULL}
      } else {
        Check_heter.grp.Fixed = NULL
      }


      if (!is.null(rand_term)){
        Check_heter.grp.Rand  = match(heter_groups, rand_term)
        if(anyNA(Check_heter.grp.Rand)) {Check_heter.grp.Rand = NULL}
      } else {
        Check_heter.grp.Rand = NULL
      }

      #if(VarCov_str=="us" |VarCov_str=="corgh" | VarCov_str=="corgv" | VarCov_str=="corh" | VarCov_str=="corv"){

      #### When user provide only the Interaction term was provided by the user
      #if(!is.null(Inter_Gen_pos) & is.null(Gen_pos)){
      if(!is.null(Inter_Gen_pos)){
        if (length(Check_heter.grp.Rand)==0 & length(Check_heter.grp.Fixed)==0){

          ### if Inter_Gen_pos is greater than 1 (Multiple kernel) but Gen_pos is null
          #if (length(Inter_Gen_pos)>1){

          #if (exists('G_list')){ Inter_Gen_pos_use = length(G_list) }


          for (i in 1:length(G_list)){

            if(VarCov_str=="us" |VarCov_str=="corgh" | VarCov_str=="corgv" | VarCov_str=="corh" | VarCov_str=="corv"){
              #if(exists("G_inv") | exists("GK_inv")){
              random= stats::update(random,
                                    paste(paste("~ . +", (paste(paste0(VarCov_str, paste0("(",heter_groups,")")),
                                                                paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
                                                                sep = ":"))), "+", heter_groups))



            }else{

              if(isTRUE(grepl("fa", VarCov_str))) {
                N_fa = substr(VarCov_str, 3, 100)
                #if(exists("GK_inv") | exists("G_inv")){
                random= stats::update(random,
                                      paste(paste("~ . +", (paste(paste0("fa", paste0("(",paste0(heter_groups, ",", N_fa),")")),
                                                                  paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
                                                                  sep = ":"))), "+", heter_groups))


              }


              if(isTRUE(grepl("rr", VarCov_str))) {
                N_rr = substr(VarCov_str, 3, 100)
                #if(exists("GK_inv") | exists("G_inv")){
                random= stats::update(random,
                                      paste(paste("~ . +", (paste(paste0("rr", paste0("(",paste0(heter_groups, ",", N_rr),")")),
                                                                  paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
                                                                  sep = ":"))), "+", heter_groups))


              }

            }


          } ### End



        } ### when Env is missing in both fixed and random terms


        ## If user provide ENV in the fixed term and missing in the random term
        if (length(Check_heter.grp.Rand)==0 & length(Check_heter.grp.Fixed)==1){


          for (i in 1:length(G_list)){


            if(VarCov_str=="us" |VarCov_str=="corgh" | VarCov_str=="corgv" | VarCov_str=="corh" | VarCov_str=="corv"){
              random =stats::update(random,
                                    paste("~ . +",paste(paste0(VarCov_str, paste0("(",heter_groups,")")),
                                                        paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))

            } else{

              if(isTRUE(grepl("fa", VarCov_str))) {
                N_fa = substr(VarCov_str, 3, 100)

                random =stats::update(random,
                                      paste("~ . +",paste(paste0("fa", paste0("(",paste0(heter_groups, ",", N_fa),")")),
                                                          paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))



              }


              if(isTRUE(grepl("rr", VarCov_str))) {
                N_rr = substr(VarCov_str, 3, 100)

                random =stats::update(random,
                                      paste("~ . +",paste(paste0("rr", paste0("(",paste0(heter_groups, ",", N_rr),")")),
                                                          paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))


              }


            }



          } ### End


        } ### End  If user provide ENV in the fixed term and missing in the random term

        ## If user provide ENV in the random term and missing in the fixed term
        if (length(Check_heter.grp.Rand)==1 & length(Check_heter.grp.Fixed)==0){


          for (i in 1:length(G_list)){

            if(VarCov_str=="us" |VarCov_str=="corgh" | VarCov_str=="corgv" | VarCov_str=="corh" | VarCov_str=="corv"){

              random= stats::update(random,
                                    paste("~ . +",paste(paste0(VarCov_str, paste0("(",heter_groups,")")),
                                                        paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))

            } else{

              if(isTRUE(grepl("fa", VarCov_str))) {
                N_fa = substr(VarCov_str, 3, 100)

                random= stats::update(random,
                                      paste("~ . +",paste(paste0("fa",paste0("(",paste0(heter_groups, ",", N_fa),")")),
                                                          paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))


              }


              if(isTRUE(grepl("rr", VarCov_str))) {
                N_rr = substr(VarCov_str, 3, 100)

                random= stats::update(random,
                                      paste("~ . +",paste(paste0("rr",paste0("(",paste0(heter_groups, ",", N_rr),")")),
                                                          paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))


              }
            }


          }
        }

      }
      ### If user provide GID:Env
      ### use this step to drop the orginal GID:Env
      if(!is.null(Inter_Gen_pos) & length(Gen_pos)==0){
        rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]

        rand_termCopy = gsub(" ", "", rand_termCopy)

        #Inter_Gen_pos_copy = match(rand_term[Inter_Gen_pos], rand_termCopy)
        if(!is.null(heter_groups)){
          Inter_Gen_pos_copy = match(paste(gen_name,heter_groups, sep = ":"),rand_termCopy)
          if(anyNA(Inter_Gen_pos_copy)){
            Inter_Gen_pos_copy= match(paste(heter_groups,gen_name, sep = ":"),rand_termCopy)
          }

        }

        random = stats::formula(stats::drop.terms(stats::terms(random),Inter_Gen_pos_copy, keep.response = F))

      } ### End

      ### If user provide GID, GID:Env
      if(!is.null(Inter_Gen_pos) & !is.null(Gen_pos)){
        ### use this step to drop the orginal GID and GID:Env
        rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]

        rand_termCopy = gsub(" ", "", rand_termCopy)

        #Inter_Gen_pos_copy = match(rand_term[Inter_Gen_pos], rand_termCopy)
        if(!is.null(heter_groups)){
          Inter_Gen_pos_copy = match(paste(gen_name,heter_groups, sep = ":"),rand_termCopy)
          if(anyNA(Inter_Gen_pos_copy)){
            Inter_Gen_pos_copy= match(paste(heter_groups,gen_name, sep = ":"),rand_termCopy)
          }

        }
        Gen_pos_copy = match(rand_term[Gen_pos], rand_termCopy)

        random = stats::formula(stats::drop.terms(stats::terms(random), c(Gen_pos_copy,Inter_Gen_pos_copy), keep.response = F))

      }

      #}

      #}### End of when only Gen_pos and Inter_Gen_pos are provided
      ## That is user provide GID, GID:Env

      #####
      # #### When Interaction term and Gen were provided by the user that is GID + GID:Env
      # if(!is.null(Inter_Gen_pos) & !is.null(Gen_pos)){
      #   if (length(Check_heter.grp.Rand)==0 & length(Check_heter.grp.Fixed)==0){
      #
      #     ### if Inter_Gen_pos is greater than 1 (Multiple kernel) but Gen_pos is null
      #     #if (length(Inter_Gen_pos)>1){
      #
      #     #if (exists('G_list')){ Inter_Gen_pos_use = length(G_list) }
      #
      #
      #     for (i in 1:length(G_list)){
      #
      #       if(VarCov_str=="us" |VarCov_str=="corgh" | VarCov_str=="corgv" | VarCov_str=="corh" | VarCov_str=="corv"){
      #         #if(exists("G_inv") | exists("GK_inv")){
      #         random= stats::update(random,
      #                               paste(paste(paste("~ . +", (paste(paste0(VarCov_str, paste0("(",heter_groups,")")),
      #                                                                 paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
      #                                                                 sep = ":"))), "+", heter_groups), "+",
      #                                     paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))
      #
      #
      #
      #       }else{
      #
      #         if(isTRUE(grepl("fa", VarCov_str))) {
      #           N_fa = substr(VarCov_str, 3, 100)
      #           #if(exists("GK_inv") | exists("G_inv")){
      #           random= stats::update(random,
      #                                 paste(paste(paste("~ . +", (paste(paste0("fa", paste0("(",paste0(heter_groups, ",", N_fa),")")),
      #                                                                   paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
      #                                                                   sep = ":"))), "+", heter_groups), "+",
      #                                       paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))
      #
      #
      #         }
      #
      #
      #         if(isTRUE(grepl("rr", VarCov_str))) {
      #           N_rr = substr(VarCov_str, 3, 100)
      #           #if(exists("GK_inv") | exists("G_inv")){
      #           random= stats::update(random,
      #                                 paste(paste(paste("~ . +", (paste(paste0("rr", paste0("(",paste0(heter_groups, ",", N_rr),")")),
      #                                                                   paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
      #                                                                   sep = ":"))), "+",  heter_groups), "+",
      #                                       paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))
      #
      #         }
      #
      #       }
      #
      #
      #     } ### End
      #
      #
      #
      #   } ### when Env is missing in both fixed and random terms
      #
      #
      #   ## If user provide ENV in the fixed term and missing in the random term
      #   if (length(Check_heter.grp.Rand)==0 & length(Check_heter.grp.Fixed)==1){
      #
      #
      #     for (i in 1:length(G_list)){
      #
      #
      #       if(VarCov_str=="us" |VarCov_str=="corgh" | VarCov_str=="corgv" | VarCov_str=="corh" | VarCov_str=="corv"){
      #         random =stats::update(random,
      #                               paste(paste("~ . +",paste(paste0(VarCov_str, paste0("(",heter_groups,")")),
      #                                                         paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")),"+",
      #                                     paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))
      #
      #
      #       } else{
      #
      #         if(isTRUE(grepl("fa", VarCov_str))) {
      #           N_fa = substr(VarCov_str, 3, 100)
      #
      #           random =stats::update(random,
      #                                 paste(paste("~ . +",paste(paste0("fa", paste0("(",paste0(heter_groups, ",", N_fa),")")),
      #                                                           paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")),"+",
      #                                       paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))
      #
      #
      #
      #         }
      #
      #
      #         if(isTRUE(grepl("rr", VarCov_str))) {
      #           N_rr = substr(VarCov_str, 3, 100)
      #
      #           random =stats::update(random,
      #                                 paste(paste("~ . +",paste(paste0("rr", paste0("(",paste0(heter_groups, ",", N_rr),")")),
      #                                                           paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")), "+",
      #                                       paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))
      #
      #
      #         }
      #
      #
      #       }
      #
      #
      #
      #     } ### End
      #
      #
      #   } ### End  If user provide ENV in the fixed term and missing in the random term
      #
      #   ## If user provide ENV in the random term and missing in the fixed term
      #   if (length(Check_heter.grp.Rand)==1 & length(Check_heter.grp.Fixed)==0){
      #
      #
      #     for (i in 1:length(G_list)){
      #
      #       if(VarCov_str=="us" |VarCov_str=="corgh" | VarCov_str=="corgv" | VarCov_str=="corh" | VarCov_str=="corv"){
      #
      #         random=  stats::update(random,
      #                                paste(paste("~ . +",paste(paste0(VarCov_str, paste0("(",heter_groups,")")),
      #                                                          paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")), "+",
      #                                      paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))
      #       } else{
      #
      #         if(isTRUE(grepl("fa", VarCov_str))) {
      #           N_fa = substr(VarCov_str, 3, 100)
      #
      #           random= stats::update(random,
      #                                 paste(paste("~ . +",paste(paste0("fa",paste0("(",paste0(heter_groups, ",", N_fa),")")),
      #                                                           paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")), "+",
      #                                       paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))
      #
      #
      #         }
      #
      #
      #         if(isTRUE(grepl("rr", VarCov_str))) {
      #           N_rr = substr(VarCov_str, 3, 100)
      #
      #           random= stats::update(random,
      #                                 paste(paste("~ . +",paste(paste0("rr",paste0("(",paste0(heter_groups, ",", N_rr),")")),
      #                                                           paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")), "+",
      #                                       paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))
      #
      #
      #         }
      #       }
      #
      #
      #     }
      #   }
      #
      #   # ### If user provide GID:Env
      #   # ### use this step to drop the orginal GID:Env
      #   # if(!is.null(Inter_Gen_pos) & length(Gen_pos)==0){
      #   #   rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]
      #   #
      #   #   rand_termCopy = gsub(" ", "", rand_termCopy)
      #   #
      #   #   #Inter_Gen_pos_copy = match(rand_term[Inter_Gen_pos], rand_termCopy)
      #   #   if(!is.null(heter_groups)){
      #   #     Inter_Gen_pos_copy = match(paste(gen_name,heter_groups, sep = ":"),rand_termCopy)
      #   #     if(anyNA(Inter_Gen_pos_copy)){
      #   #       Inter_Gen_pos_copy= match(paste(heter_groups,gen_name, sep = ":"),rand_termCopy)
      #   #     }
      #   #
      #   #   }
      #   #
      #   #   random = stats::formula(stats::drop.terms(stats::terms(random),Inter_Gen_pos_copy, keep.response = F))
      #   #
      #   # } ### End
      #
      #   ### If user provide GID, GID:Env
      #   #if(!is.null(Inter_Gen_pos) & !is.null(Gen_pos)){
      #   ### use this step to drop the orginal GID and GID:Env
      #   rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]
      #
      #   rand_termCopy = gsub(" ", "", rand_termCopy)
      #
      #   #Inter_Gen_pos_copy = match(rand_term[Inter_Gen_pos], rand_termCopy)
      #   if(!is.null(heter_groups)){
      #     Inter_Gen_pos_copy = match(paste(gen_name,heter_groups, sep = ":"),rand_termCopy)
      #     if(anyNA(Inter_Gen_pos_copy)){
      #       Inter_Gen_pos_copy= match(paste(heter_groups,gen_name, sep = ":"),rand_termCopy)
      #     }
      #
      #   }
      #   Gen_pos_copy = match(rand_term[Gen_pos], rand_termCopy)
      #
      #   random = stats::formula(stats::drop.terms(stats::terms(random), c(Gen_pos_copy,Inter_Gen_pos_copy), keep.response = F))
      #
      #   #}
      #
      #   #}
      #
      # }### End of when only Gen_pos and Inter_Gen_pos are provided
      # ## That is user provide GID, GID:Env

    }

    ##################################################
    ### If user provide only the gen_name in the random term as gen_name effect
    ######################################################

    if(is.null(Inter_Gen_pos) & length(Gen_pos)==1){
      #if (length(Gen_pos)==1 & length(Inter_Gen_pos)==0){
      ### stats::formula must be at least  length 1 so since it just one length. The
      ## adjusted random term was added and the old one was removed
      if (exists('G_list')){
        for (i in 1:length(G_list)) {


          random =stats::update(random, paste("~ . +",paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))

        }
      }
      ### use this step to drop the orginal GID and GID:Env
      rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]

      rand_termCopy = gsub(" ", "", rand_termCopy)

      Gen_pos_copy = match(rand_term[Gen_pos], rand_termCopy)

      random = stats::formula(stats::drop.terms(stats::terms(random), Gen_pos_copy, keep.response = F))

    }
    ######################################################################
    #### When Variance-Covariance Structure is missing. Compound Symmetry
    ####################################################################
    if(is.null(VarCov_str) & (length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name])))){

      #if(is.null(VarCov_str) & !is.null(heter_groups)){
      #if(dim(pheno_data)[1]!=dim(Geno_data)[1] | dim(pheno_data)[1]!=dim(Gmatrix)[1]){
      # if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
      #   warning(paste(msg, "Variance-Covariance structure is missing."), immediate. = TRUE, call. =FALSE)
      #   if(is.null(heter_groups)){
      #     stop(print('Provide heter_groups, since gen_name is present in more than one location/environment.'), call. = FALSE)
      #   }
      #
      # }

      ###############################
      if (!is.null(fixed_term)){
        Check_heter.grp.Fixed  = match(heter_groups, fixed_term)
      } else {
        Check_heter.grp.Fixed = NULL
      }


      if (!is.null(rand_term)){
        Check_heter.grp.Rand  = match(heter_groups, rand_term)
      } else {
        Check_heter.grp.Rand = NULL
      }

      #if(VarCov_str=="us" |VarCov_str=="corgh" | VarCov_str=="corgv" | VarCov_str=="corh" | VarCov_str=="corv"){

      if(!is.null(Inter_Gen_pos)){
        if (length(Check_heter.grp.Rand)==0 & length(Check_heter.grp.Fixed)==0){

          ### if Inter_Gen_pos is greater than 1 (Multiple kernel) but Gen_pos is null
          #if (length(Inter_Gen_pos)>1){

          for (i in 1:length(G_list)){

            # random =stats::update(random, paste(paste("~ . +",paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))),
            #                                     "+", heter_groups))

            random= stats::update(random,
                                  paste(paste("~ . +", (paste(paste0("idv", paste0("(",heter_groups,")")),
                                                              paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
                                                              sep = ":"))), "+", heter_groups))


          }

        }### when Env is missing in both fixed and random terms


        ## If user provide ENV in the fixed term and missing in the random term
        if (length(Check_heter.grp.Fixed)==1 & length(Check_heter.grp.Rand)==0){

          for (i in 1:length(G_list)){

            #random =stats::update(random, paste("~ . +",paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))

            random= stats::update(random,
                                  paste("~ . +",paste(paste0("idv", paste0("(",heter_groups,")")),
                                                      paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))


          }

        }### End  If user provide ENV in the fixed term and missing in the random term

        ## If user provide ENV in the random term and missing in the fixed term
        if (length(Check_heter.grp.Fixed)==0 & length(Check_heter.grp.Rand)==1){

          if (exists('G_list')){Inter_Gen_pos_use = length(G_list)}
          for (i in 1:length(G_list)){

            #random =stats::update(random, paste("~ . +",paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))

            random= stats::update(random,
                                  paste("~ . +",paste(paste0("idv", paste0("(",heter_groups,")")),
                                                      paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))


          }

        }


      }

      ### use this step to drop the orginal GID and GID:Env
      if(!is.null(Gen_pos) & is.null(Inter_Gen_pos)){
        rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]

        rand_termCopy = gsub(" ", "", rand_termCopy)

        Gen_pos_copy = match(rand_term[Gen_pos], rand_termCopy)

        random = stats::formula(stats::drop.terms(stats::terms(random), Gen_pos_copy, keep.response = F))
      }

      if(!is.null(heter_groups)){
        Inter_Gen_pos_copy = match(paste(gen_name,heter_groups, sep = ":"),rand_termCopy)
        if(anyNA(Inter_Gen_pos_copy)){
          Inter_Gen_pos_copy= match(paste(heter_groups,gen_name, sep = ":"),rand_termCopy)
        }

      }
      Gen_pos_copy = match(rand_term[Gen_pos], rand_termCopy)

      random = stats::formula(stats::drop.terms(stats::terms(random), c(Gen_pos_copy,Inter_Gen_pos_copy), keep.response = F))


    } ## End of when variance-covariance str is not provided by user

    #####################################################
    ####
    # Extract each random component/terms
    ####################################################

    ranTerms <- strsplit(as.character(random[2]), split = "[+]")[[1]]
    ranTerms = gsub(" ", "", ranTerms)

    for (r in 1:length(ranTerms)) {

      if(r == 1) {
        code.asr[2] <- paste(code.asr[2], ranTerms[r], sep='')

      } else {

        code.asr[2] <- paste(code.asr[2], ranTerms[r], sep='+')
      }
    }


  } ## End of fixing the random terms

  # Heterogeneous errors
  if (!is.null(heter_groups)&isTRUE(heter_resid)) {

    if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
      code.asr[3] <- paste(code.asr[3], paste0('dsum(~id(units)|', paste0(heter_groups, ')')), sep='')

    } else {
      warning(paste(msg,'Heterogenous residual is not possible with one environment/location. We fix it for you.'),
              call. = FALSE)
      if(length(pheno_data[,gen_name])==length(unique(pheno_data[,gen_name]))){
        code.asr[3] <- paste(code.asr[3], 'id(units)', sep='')
      }
    }
    #code.asr[3] <- paste(code.asr[3], paste0('dsum(~idv(units)|', paste0(heter_groups, ')')), sep='')

    #code.asr[3] <- paste(code.asr[3], 'dsum(~idv(units|', 'heter_groups)', sep='')
  } else {
    #code.asr[3] <- paste(code.asr[3], 'idv(units)', sep='')
    code.asr[3] <- paste(code.asr[3], 'id(units)', sep='')
  }
  # Adding individual/gen_name (random)
  #code.asr[2] <- paste(code.asr[2], 'vm(indiv,ainv)', sep='+')

  Univariate = c()
  for (trait in 1: length(response)) {

    if (trait==1){
      code.asr[1] <-  gsub("trait", response[trait], code.asr[1])
    } else {

      code.asr[1] <-  gsub(response[trait-1], response[trait], code.asr[1])
    }


    if(is.null(weights)){
      code.asr[4] <- 'na.action=list(x="include",y="include"),data=pheno_data)'
    } else {
      if(!is.null(weights)){
        code.asr[4] <- 'na.action=list(x="include",y="include"), weights = weights, family = asr_gaussian(dispersion = 1), data=pheno_data)'
      }

      if(length(unique(response[trait]))<=10 & is.null(weights)){
        code.asr[4] <- 'na.action=list(x="include",y="include"), family = asr_multinomial(),  data=pheno_data)'
      }

      if(length(unique(response[trait]))<=10 & !is.null(weights)){
        code.asr[4] <- 'na.action=list(x="include",y="include"), weights = weights, family = asr_multinomial(dispersion = 1),  data=pheno_data)'
      }

      if(length(unique(response[trait]))==2 & !is.null(weights)){
        code.asr[4] <- 'na.action=list(x="include",y="include"), weights = weights, family = asr_binomial(dispersion = 1),  data=pheno_data)'
      }

      if(length(unique(response[trait]))==2 & is.null(weights)){
        code.asr[4] <- 'na.action=list(x="include",y="include"), family = asr_binomial(),  data=pheno_data)'
      }

    }



    ####
    code.asr[1] <- paste('mod<-', code.asr[1], sep='')
    str.mod <- paste(code.asr[1],code.asr[2],code.asr[3],code.asr[4],sep=',')
    ## Calls the current environment for evaluation
    eval(parse(text=str.mod), envir=environment())
    if (!mod$converge) { eval(parse(text='mod<-asreml::update.asreml(mod)')) }

    ###################################################
    ##### Start the process of processing the results
    ################################################



    #IND = paste(paste0('vm(', gen_name), collapse = ',', 'G_inv)')
    # Initiate Process to extract BLUPs/BV from the model
    # BLUP <- summary(mod, coef=TRUE)$coef.random
    #
    # #Heter.Grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))
    #
    # if (!is.null(heter_groups)){
    #   Heter.Grp <- as.character(unique(pheno_data[, heter_groups]))
    # }
    # #ENV_Ids = as.character(unique(pheno_data[, heter_groups]))
    # #### Extract Breeding values/genetic effect estimate for all omics
    # BV_All = vector(mode = 'list', length = length(G_list))
    # for (b in 1:length(G_list)) {
    #
    #   if(b ==1){
    #     BV_All[[b]] <- BLUP[grep(paste(paste(G_list[[b]], '_inv', sep = ""),"\\)", sep = ""),rownames(BLUP)),]
    #   }
    #
    #   if(b>1){
    #
    #     BV_All[[b]] <- BLUP[grep(paste(paste(paste(G_list[[b]],b, sep=""), '_inv', sep = ""),"\\)", sep = ""),rownames(BLUP)),]
    #   }
    # }
    #
    #
    # #BV <- BLUP[grep(paste(paste0('vm\\(', gen_name), sep = ',', 'G_inv\\)'),rownames(BLUP)),]
    # #BV <- BLUP[grep('G_inv\\)',rownames(BLUP)),]
    #
    # if(!is.null(VarCov_str)){
    #
    #   if(isTRUE(grepl("fa", VarCov_str))){
    #     ## Extract the number of factors
    #     #N_fa = substr(VarCov_str, 3, 100)
    #
    #     for (bb in 1:length(G_list)) {
    #       BV_All[[bb]] <- BV_All[[bb]][!rownames(BV_All[[bb]])%in%rownames(BV_All[[bb]][grep('Comp',rownames(BV_All[[bb]])),]), ]
    #
    #       BV_All[[bb]] <- as.data.frame(BV_All[[bb]])
    #
    #       BV_All[[bb]][, "GID"]<-as.character(stringr::str_split_fixed(rownames(BV_All[[bb]]), "\\)_", 3)[,3])
    #
    #     }
    #
    #   } else {
    #
    #     if(VarCov_str=="us" |VarCov_str=="corgh" | VarCov_str=="corgv" | VarCov_str=="corh" | VarCov_str=="corv") {
    #
    #
    #       for (bb in 1:length(G_list)) {
    #
    #         BV_All[[bb]] <- as.data.frame(BV_All[[bb]])
    #
    #         BV_All[[bb]][, "GID"] <- as.character(stringr::str_split_fixed(rownames(BV_All[[bb]]), "\\)_", 3)[,2])
    #
    #       }
    #
    #       # BV <- as.data.frame(BV)
    #       #
    #       # BV[, "GID"]<-as.character(stringr::str_split_fixed(rownames(BV), "\\)_", 3)[,2])
    #     }
    #
    #
    #   } ## End
    #
    #   #################################
    #   for (bb in 1:length(G_list)) {
    #
    #     BV_All[[bb]] <- BV_All[[bb]][, c(4, 1:2)]
    #
    #     BV_All[[bb]][, heter_groups] <- rep(Heter.Grp, each=length(unique(BV_All[[bb]][, "GID"])))
    #
    #     BV_All[[bb]] <- BV_All[[bb]][, c(1, 4, 2:3)]
    #
    #     colnames(BV_All[[bb]])[1:3] <- c(gen_name, heter_groups, "BLUP")
    #
    #     BV_All[[bb]][, "PEV"] <-  BV_All[[bb]][, "std.error"]^2
    #
    #   }
    #
    #   #######################################
    #
    #   # BV = BV[, c(4, 1:2)]
    #   #
    #   # BV[, heter_groups] = rep(Heter.Grp, each=length(unique(BV$GID)))
    #   #
    #   # BV = BV[, c(1, 4, 2:3)]
    #   #
    #   # colnames(BV)[1:3] <- c(gen_name, heter_groups, "BLUP")
    #
    #
    #   Res= VarCov_Herr(model= mod,
    #                    heter_groups=heter_groups,
    #                    VarCov_str= VarCov_str,
    #                    heter_resid=heter_resid,
    #                    N_Omics = G_list)
    #
    #
    #   #BV$PEV <- BV$std.error^2
    #   #VA = Res$Genetic_Var
    #   Heter.Grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))
    #
    #   for (bb in 1:length(G_list)){
    #
    #     VA = Res$Genetic_Var[bb, ]
    #
    #     if(length(VA)< length(Heter.Grp)){
    #
    #       print(paste("Reliability cannot be estimated. Not all varaince components for", heter_groups, "are postive definitive"))
    #     } else{
    #
    #       if(length(VA) == length(Heter.Grp)){
    #
    #         BV_All[[bb]][, "Reliability"] = NA
    #         #BV$Reliability = NA
    #         for (i in 1:length(VA)) {
    #
    #           BV_All[[bb]][, "Reliability"] <- ifelse(BV_All[[bb]][, heter_groups]%in% Heter.Grp[i],
    #                                                   round(1 - BV_All[[bb]][, "PEV"]/VA[i],6), BV_All[[bb]][, "Reliability"])
    #
    #         }
    #       }
    #
    #
    #     }
    #
    #   }
    #
    #
    # } else {
    #   ## Problem
    #   if(is.null(VarCov_str)){
    #
    #     for (bb in 1:length(G_list)){
    #       BV_All[[bb]] <- as.data.frame(BV_All[[bb]])
    #       BV_All[[bb]][, "GID"] <- as.character(stringr::str_split_fixed(rownames(BV_All[[bb]]), "\\)_", 3)[,2])
    #
    #       BV_All[[bb]] =  BV_All[[bb]][, c(4, 1:2)]
    #       colnames(BV_All[[bb]])[1:2] <- c(gen_name, "BLUP")
    #     }
    #
    #     vc <- summary(mod)$varcomp
    #     #VAR_check_Pos <- which(vc$bound=='F' | vc$bound=='U' | vc$bound=="?" | vc$bound=="S")
    #     VAR_check_Pos <- which(vc$bound=="?" | vc$bound=="S")
    #     if(length(VAR_check_Pos)>1) {
    #
    #       stop(print( paste(paste("variance component for", as.character(rownames(vc)[VAR_check_Pos]), collapse =" and" ),
    #                         " is unstable, refix the model")), call. = FALSE)
    #
    #     } else {
    #
    #       if(length(VAR_check_Pos)==1) {
    #
    #         stop(print(paste(paste("variance component for", as.character(rownames(vc)[VAR_check_Pos]), collapse =" " ),
    #                          " is unstable, refix the model")), call. = FALSE)
    #
    #       }
    #
    #     }
    #
    #     VarG_All = vector("list", length = length(G_list))
    #     VarG <- vc[grep(paste("^vm\\(", gen_name, sep = ""), rownames(vc)), drop=FALSE, ]
    #     for (g in 1:length(G_list)) {
    #
    #       if(g ==1){
    #         VarG_All[[g]] <-VarG[grep(paste(G_list[[g]],"_inv", sep = ""), rownames(VarG)), ]
    #         CheckR <- VarG_All[[g]][grep("!R", rownames(VarG_All[[g]])), "component"]
    #         if (length(CheckR)!=0){
    #           VarG_All[[g]] <-  VarG_All[[g]][-grep("!R", rownames(VarG_All[[g]])), "component"]
    #         }else{
    #           VarG_All[[g]] <-  VarG_All[[g]][, "component"]
    #
    #         }
    #       }
    #
    #       if(g>1){
    #
    #         VarG_All[[g]] <- VarG[grep(paste(paste(G_list[[g]], g, sep=""),"_inv", sep = ""), rownames(VarG)), ]
    #         CheckR <- VarG_All[[g]][grep("!R", rownames(VarG_All[[g]])), "component"]
    #         if (length(CheckR)!=0){
    #           VarG_All[[g]] <-  VarG_All[[g]][-grep("!R", rownames(VarG_All[[g]])), "component"]
    #         }else{
    #           VarG_All[[g]] <-  VarG_All[[g]][, "component"]
    #
    #         }
    #       }
    #     }
    #
    #     #VA <- vc[match(paste0("vm(", gen_name, paste(",","G_inv)", sep = " ")), rownames(vc)),"component"]
    #     VE <- vc[grep("!R", row.names(vc)), "component"]
    #
    #     varG_matrix = unlist(VarG_All)
    #
    #     H = matrix(NA, nrow = 1, ncol = length(VE))
    #
    #
    #     #BV$PEV <- BV$std.error^2
    #     BV_All$PEV <- BV_All$std.error^2
    #     for (bb in 1:length(G_list)){
    #
    #       BV_All[[bb]][, "PEV"] <- BV_All[[bb]][, "std.error"]^2
    #
    #       BV_All[[bb]][, "Reliability"] <- round(1 - BV_All[[bb]][, "PEV"]/mean(varG_matrix[bb]),6)
    #
    #       if(length(VE)>1){
    #
    #         for (i in 1:length(VE)) {
    #           H[, i] <- sum(varG_matrix)/(sum(varG_matrix)+VE[i])
    #
    #         }
    #
    #       } else {
    #
    #         if (length(VE)==1){
    #
    #           #BV$Reliability <- round(1 - BV$PEV/VA,6)
    #
    #           H <- sum(varG_matrix)/(sum(varG_matrix)+VE)
    #         }
    #       }
    #
    #     }
    #
    #   } ## End
    #
    # } #### End var_Covar
    #
    # ### Predicted Values
    # pred_value <- asreml::predict.asreml(mod, classify=gen_name, sed=FALSE)$pvals
    #
    # if (is.null(heter_groups) & is.null(VarCov_str)){Inter_Gen_pos= NULL}
    # if(length(Gen_pos) == length(rand_term)){Inter_Gen_pos= NULL}
    # if(!is.null(Inter_Gen_pos)){
    #
    #   pred_heter_groups <- asreml::predict.asreml(mod, classify= rand_term[[Inter_Gen_pos]], sed=FALSE)$pvals
    #
    # }
    # ### Results
    #
    # # if (!is.null(VarCov_str)){
    # #
    # #   Result = list(call=str.mod, mod=mod, ebv=BV,pred_value =pred_value,
    # #                 BIC=summary(mod)$bic, AIC=summary(mod)$aic, H2 =Res$Heritability,
    # #                 Va = Res$Genetic_Var, Ve = Res$Residual_Var, CORR= Res$Correlation, COV=Res$Covariance)
    # #
    # #} else {
    #
    # if (!is.null(VarCov_str) & length(Inter_Gen_pos)!=0){
    #
    #   Result = list(call=str.mod, mod=mod, ebv=BV_All,pred_value =pred_value, pred_heter_groups = pred_heter_groups,
    #                 BIC=summary(mod)$bic, AIC=summary(mod)$aic, H2 =Res$Heritability,
    #                 Va = Res$Genetic_Var, Ve = Res$Residual_Var, CORR= Res$Correlation, COV=Res$Covariance)
    #
    #
    # } else {
    #
    #   #if (is.null(VarCov_str) & length(rand_inter_Pos)==0){
    #   if (is.null(VarCov_str) & length(Inter_Gen_pos)==0){
    #
    #     Result = list(call=str.mod, mod=mod, ebv=BV_All,pred_value =pred_value,
    #                   BIC=summary(mod)$bic, AIC=summary(mod)$aic, H2 =H,
    #                   Va = varG_matrix, Ve = VE)
    #
    #
    #   }
    #
    #   #if (is.null(VarCov_str) & length(rand_inter_Pos)!=0){
    #   if (is.null(VarCov_str) & length(Inter_Gen_pos)!=0){
    #
    #     Result = list(call=str.mod, mod=mod, ebv=BV_All,pred_value =pred_value, pred_heter_groups = pred_heter_groups,
    #                   BIC=summary(mod)$bic, AIC=summary(mod)$aic, H2 =H,
    #                   Va = varG_matrix, Ve = VE)
    #
    #
    #
    #   }
    #
    # } ## Result Ends

    #Univariate <-  Result
    #Univariate[[trait]] = Result

    Univariate[[trait]] = mod
  } ## End loop for multiple response variables

  names(Univariate) <- response

  #names(Univariate)[trait] <- response[trait]
  #parallel::stopCluster(cl)

  #doParallel::stopImplicitCluster()

  output <- list(Univariate,
                 str.mod,
                 G_list,
                 Gen_pos,
                 Inter_Gen_pos,
                 rand_term)

  names(output) <- c("model",
                     "str.mod",
                     "G_list",
                     "Gen_pos",
                     "Inter_Gen_pos",
                     "rand_term")
  #return(c(Univariate, G_list))
  return(output)

}




