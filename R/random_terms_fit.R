
random_terms_fit <- function(random = NULL,
                              fixed = NULL,
                              fixed_term = NULL,
                              heter_groups = NULL,
                              heter_resid = NULL,
                              VarCov_str = NULL,
                              code_asr = NULL,
                              G_list = NULL,
                              gen_name = NULL,
                             pheno_data= NULL,
                              ...){
  if (!is.null(random)){
    #if (length(all.vars(random))>1) {
    rand_term <- strsplit(as.character(random[2]), split = "[+]")[[1]] # random parts
    rand_term = gsub(" ", "", rand_term)


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

        Inter_Non_Gen_pos <- match(Non_Gen_Inter_Test, rand_term)

      }

      if(is.null(Inter_Gen_pos)){warning(paste(gen_name, 'missing in the interaction term. This implies you cannot estimate GxE'))}

    } else {

      if(length(Check_rand_Inter)==0){

        Inter_Gen_pos = NULL
      }
    }

    Gen_pos=  match(gen_name, rand_term)
    if(anyNA(Gen_pos)){Gen_pos = NULL}
    if(length(Gen_pos)>1){
      stop(message(paste(msg, "Genotype main effect cannot be present more than one time in the model")), call. = FALSE)
    }

    ## Eg when GID:Env with no variance structure is the only term in the
    ## random effect provided by the user
    if(length(rand_term)==length(Check_rand_Inter) & is.null(VarCov_str)){
      stop(message(paste(msg, "provide variance-covariance structure")), call. = FALSE)

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

      ### This check if heterogeneous group/environment/location is present in the fixed term
      if (!is.null(fixed_term)){
        check_heter_grp_fixed  = match(heter_groups, fixed_term)
        if(anyNA(check_heter_grp_fixed)) {check_heter_grp_fixed = NULL}
      } else {
        check_heter_grp_fixed = NULL
      }


      if (!is.null(rand_term)){
        check_heter_grp_rand  = match(heter_groups, rand_term)
        if(anyNA(check_heter_grp_rand)) {check_heter_grp_rand = NULL}
      } else {
        check_heter_grp_rand = NULL
      }

      #if(VarCov_str=="us" |VarCov_str=="corgh" | VarCov_str=="corgv" | VarCov_str=="corh" | VarCov_str=="corv"){

      #### When user provide only the Interaction term was provided by the user
      #if(!is.null(Inter_Gen_pos) & is.null(Gen_pos)){
      if(!is.null(Inter_Gen_pos)){
        if (length(check_heter_grp_rand)==0 & length(check_heter_grp_fixed)==0){

          #### Add the environment to the fixed term by default if not provided by the user
          ## Providing it in fixed term allow for optimal model fit compared to the random term

          # Adding fixed factors
          if (!is.null(fixed_term)) {


            ### Update the fixed term if not null
            fixed =  stats::update(fixed,
                                   paste("~ . +", heter_groups))

            fixed_term <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]


            for (v in 1:length(all.vars(fixed))) {
              code_asr[1] <- paste(code_asr[1], all.vars(fixed)[v], sep='+')
            }



          } else {
            #fixed = heter_groups
            fixed_term <- heter_groups
            code_asr[1] <- paste(code_asr[1], fixed_term, sep='+')

          }

          for (i in 1:length(G_list)){

            if (VarCov_str %in% c("us", "corgh", "corgv", "corh", "corv")) {
              #if(exists("G_inv") | exists("GK_inv")){
              # ### This add the environment to the random term even when when Env is missing in both fixed and random terms
              # random= stats::update(random,
              #                       paste(paste("~ . +", (paste(paste0(VarCov_str, paste0("(",heter_groups,")")),
              #                                                   paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
              #                                                   sep = ":"))), "+", heter_groups))


              ### This add the environment to the fixed term when interaction term was only provided
              ## by the user. That is Env is missing in both fixed and random terms
              random = stats::update(random,
                                     paste("~ . +", (paste(paste0(VarCov_str, paste0("(",heter_groups,")")),
                                                           paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
                                                           sep = ":"))))



            }else{

              if(isTRUE(grepl("fa", VarCov_str))) {
                N_fa = substr(VarCov_str, 3, 100)
                #if(exists("GK_inv") | exists("G_inv")){
                # random= stats::update(random,
                #                       paste(paste("~ . +", (paste(paste0("fa", paste0("(",paste0(heter_groups, ",", N_fa),")")),
                #                                                   paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
                #                                                   sep = ":"))), "+", heter_groups))

                random=  stats::update(random,
                                       paste("~ . +", (paste(paste0("fa", paste0("(",paste0(heter_groups, ",", N_fa),")")),
                                                             paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
                                                             sep = ":"))))


              }


              if(isTRUE(grepl("rr", VarCov_str))) {
                N_rr = substr(VarCov_str, 3, 100)
                #if(exists("GK_inv") | exists("G_inv")){
                # random= stats::update(random,
                #                       paste(paste("~ . +", (paste(paste0("rr", paste0("(",paste0(heter_groups, ",", N_rr),")")),
                #                                                   paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
                #                                                   sep = ":"))), "+", heter_groups))

                random= stats::update(random,
                                      paste("~ . +", (paste(paste0("rr", paste0("(",paste0(heter_groups, ",", N_rr),")")),
                                                            paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
                                                            sep = ":"))))


              }

            }


          } ### End



        } ### when Env is missing in both fixed and random terms


        ## If user provide ENV in the fixed term and missing in the random term
        if (length(check_heter_grp_rand)==0 & length(check_heter_grp_fixed)==1){


          for (i in 1:length(G_list)){


            if (VarCov_str %in% c("us", "corgh", "corgv", "corh", "corv")) {
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
        if (length(check_heter_grp_rand)==1 & length(check_heter_grp_fixed)==0){


          for (i in 1:length(G_list)){

            if (VarCov_str %in% c("us", "corgh", "corgv", "corh", "corv")) {

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
        ### use this step to drop the original GID and GID:Env
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

      ###############################
      if (!is.null(fixed_term)){
        check_heter_grp_fixed  = match(heter_groups, fixed_term)
      } else {
        check_heter_grp_fixed = NULL
      }


      if (!is.null(rand_term)){
        check_heter_grp_rand  = match(heter_groups, rand_term)
      } else {
        check_heter_grp_rand = NULL
      }

      #if(VarCov_str=="us" |VarCov_str=="corgh" | VarCov_str=="corgv" | VarCov_str=="corh" | VarCov_str=="corv"){

      if(!is.null(Inter_Gen_pos)){
        if (length(check_heter_grp_rand)==0 & length(check_heter_grp_fixed)==0){

          #### Add the environment to the fixed term by default if not provided by the user
          ## Providing it in fixed term allow for optimal model fit compared to the random term

          # Adding fixed factors
          if (!is.null(fixed_term )) {

            ### Update the fixed term if not null
            fixed =  stats::update(fixed,
                                   paste("~ . +", heter_groups))

            fixed_term <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]


            for (v in 1:length(all.vars(fixed))) {
              code_asr[1] <- paste(code_asr[1], all.vars(fixed)[v], sep='+')
            }

          } else {

            fixed_term <- heter_groups

            code_asr[1] <- paste(code_asr[1], fixed_term, sep='+')

          }

          ###

          for (i in 1:length(G_list)){

            # random =stats::update(random, paste(paste("~ . +",paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))),
            #                                     "+", heter_groups))

            # random= stats::update(random,
            #                       paste(paste("~ . +", (paste(paste0("idv", paste0("(",heter_groups,")")),
            #                                                   paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
            #                                                   sep = ":"))), "+", heter_groups))

            random= stats::update(random,
                                  paste("~ . +", (paste(paste0("idv", paste0("(",heter_groups,")")),
                                                        paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
                                                        sep = ":"))))

          }

        }### when Env is missing in both fixed and random terms


        ## If user provide ENV in the fixed term and missing in the random term
        if (length(check_heter_grp_fixed)==1 & length(check_heter_grp_rand)==0){

          for (i in 1:length(G_list)){

            #random =stats::update(random, paste("~ . +",paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))

            random= stats::update(random,
                                  paste("~ . +",paste(paste0("idv", paste0("(",heter_groups,")")),
                                                      paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))


          }

        }### End  If user provide ENV in the fixed term and missing in the random term

        ## If user provide ENV in the random term and missing in the fixed term
        if (length(check_heter_grp_fixed)==0 & length(check_heter_grp_rand)==1){

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
      # if(!is.null(Gen_pos) & is.null(Inter_Gen_pos)){
      # rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]
      #
      # rand_termCopy = gsub(" ", "", rand_termCopy)
      #
      # Gen_pos_copy = match(rand_term[Gen_pos], rand_termCopy)
      #
      # random = stats::formula(stats::drop.terms(stats::terms(random), Gen_pos_copy, keep.response = F))
      # }

      ### If user provide GID, GID:Env
      if(!is.null(Inter_Gen_pos) & !is.null(Gen_pos)){
        ### use this step to drop the original GID and GID:Env
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

    } ## End of when variance-covariance str is not provided by user

    #####################################################
    ####
    # Extract each random component/terms
    ####################################################

    ranTerms <- strsplit(as.character(random[2]), split = "[+]")[[1]]
    ranTerms = gsub(" ", "", ranTerms)

    for (r in 1:length(ranTerms)) {

      if(r == 1) {
        code_asr[2] <- paste(code_asr[2], ranTerms[r], sep='')

      } else {

        code_asr[2] <- paste(code_asr[2], ranTerms[r], sep='+')
      }
    }


  }

  out <- list(code_asr = code_asr,
              Gen_pos = Gen_pos,
              Inter_Gen_pos = Inter_Gen_pos,
              rand_term = rand_term
              )

  return(out)
}
