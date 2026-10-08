
#' Title
#'
#' @param random
#' @param fixed
#' @param fixed_term
#' @param heter_groups
#' @param heter_resid
#' @param var_cov_str
#' @param code_asr
#' @param names_in_inv_list
#' @param gen_name
#' @param pheno_data
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
random_terms_fit_new <- function(random = NULL,
                                 fixed = NULL,
                                 fixed_term = NULL,
                                 heter_groups = NULL,
                                 heter_resid = NULL,
                                 var_cov_str = NULL,
                                 code_asr = NULL,
                                 names_in_inv_list = NULL,
                                 gen_name = NULL,
                                 pheno_data= NULL,
                                 ...){
  gp_reject_obsolete_asreml_structure(var_cov_str)
  msg <- ""
  heter_control <- gp_normalize_single_environment_heter_controls(
    pheno_data = pheno_data,
    gen_name = gen_name,
    heter_groups = heter_groups,
    heter_resid = heter_resid,
    var_cov_str = var_cov_str
  )
  heter_groups <- heter_control$heter_groups
  heter_resid <- heter_control$heter_resid
  var_cov_str <- heter_control$var_cov_str
  build_structured_vm_term <- function(i) {
    vm_term <- paste0("vm(", gen_name, ",", names_in_inv_list[i], ")")

    if (var_cov_str %in% c("us", "corgh", "corh", "corv")) {
      return(paste0(var_cov_str, "(", heter_groups, "):", vm_term))
    }

    if (isTRUE(grepl("fa", var_cov_str))) {
      N_fa <- substr(var_cov_str, 3, 100)
      vc_prefix <- paste0("fa(", heter_groups, ",", N_fa, ")")
      return(paste0(vc_prefix, ":", vm_term))
    }

    if (isTRUE(grepl("rr", var_cov_str))) {
      N_rr <- substr(var_cov_str, 3, 100)
      vc_prefix <- paste0("rr(", heter_groups, ",", N_rr, ")")
      return(paste0(vc_prefix, ":", vm_term))
    }

    vm_term
  }

  if (!is.null(random)){
    #if (length(all.vars(random))>1) {
    rand_term <- strsplit(as.character(random[2]), split = "[+]")[[1]] # random parts
    rand_term <- gsub(" ", "", rand_term)


    #rand_InterPresent = grep(":", rand_term)
    ## No interaction term
    rand_term_no_inter <- rand_term[!grepl(":", rand_term)]

    ## Interaction term
    check_rand_inter <- rand_term[grepl(":", rand_term)]

    ## Extract other side expect the gen_name
    rand_term_no_inter_no_gen <-  rand_term_no_inter[!rand_term_no_inter%in%gen_name]
    ## Get the position of other terms in the random that is not gen_name
    non_gen_pos <- match(rand_term_no_inter_no_gen, rand_term)
    if(anyNA(non_gen_pos)){non_gen_pos <- NULL}

    if(length(check_rand_inter)!=0){

      test_present_of_geno <-  grep(gen_name, check_rand_inter, value = TRUE)
      if(anyNA(test_present_of_geno)) {test_present_of_geno <- NULL}

      if(length(test_present_of_geno)!=0){

        inter_gen_pos <-  match(test_present_of_geno, rand_term)
        if(anyNA(inter_gen_pos)){inter_gen_pos <-  NULL}

        non_gen_inter_test <- check_rand_inter[!check_rand_inter%in%test_present_of_geno]
        if(anyNA(non_gen_inter_test)){non_gen_inter_test <-  NULL}

      } else {

        non_gen_inter_test <-  NULL

        inter_gen_pos <-  NULL
      }

      if(length(non_gen_inter_test)!=0){

        inter_non_gen_pos <- match(non_gen_inter_test, rand_term)

      }

      if(is.null(inter_gen_pos)){warning(paste(gen_name, 'missing in the interaction term. This implies you cannot estimate GxE'))}

    } else {

      if(length(check_rand_inter)==0){

        inter_gen_pos <-  NULL
      }
    }

    gen_pos <- match(gen_name, rand_term)
    if(anyNA(gen_pos)){gen_pos <-  NULL}
    if(length(gen_pos)>1){
      stop(message(paste(msg, "Genotype main effect cannot be present more than one time in the model")), call. = FALSE)
    }

    ## Eg when GID:Env with no variance structure is the only term in the
    ## random effect provided by the user
    if(length(rand_term)==length(check_rand_inter) & is.null(var_cov_str)){
      stop(message(paste(msg, "provide variance-covariance structure")), call. = FALSE)

    }


    if (length(check_rand_inter)>=1){
      if (length(pheno_data[[gen_name]]) ==length(unique(pheno_data[[gen_name]]))){
        stop(paste(msg, "Phenotypic data contain single environment but you specify multi-environment analysis."), call. = FALSE)
      }
      if (is.null(heter_groups)) {stop(paste(msg, heter_groups,"cannot be NULL"), call. = FALSE)}
    }


    if(!is.null(var_cov_str)){

      if(is.null(heter_groups)){stop(paste(msg, "Provide heter_groups to model specified variance-covariance structure"), call. = FALSE)}

      ### Check if the number of hetero.Grp is greater 5 or greater than 5
      NN <-  nlevels(pheno_data[[heter_groups]])
      if (NN >=5 & isFALSE(grepl("fa", var_cov_str))){

        #msg <- ""

        warning(paste(msg, "The number of", heter_groups, " is", NN,  "consider using factor analytic model"), immediate. = TRUE, call. =FALSE)
      }

      ### This check if heterogeneous group/environment/location is present in the fixed term
      if (!is.null(fixed_term)){
        check_heter_grp_fixed  <-  match(heter_groups, fixed_term)
        if(anyNA(check_heter_grp_fixed)) {check_heter_grp_fixed <-  NULL}
      } else {
        check_heter_grp_fixed <-  NULL
      }


      if (!is.null(rand_term)){
        check_heter_grp_rand  <-  match(heter_groups, rand_term)
        if(anyNA(check_heter_grp_rand)) {check_heter_grp_rand <- NULL}
      } else {
        check_heter_grp_rand <-  NULL
      }

      #if(var_cov_str=="us" |var_cov_str=="corgh" | var_cov_str=="corh" | var_cov_str=="corv"){

      #### When user provide only the Interaction term was provided by the user
      #if(!is.null(inter_gen_pos) & is.null(gen_pos)){
      if(!is.null(inter_gen_pos)){
        if (length(check_heter_grp_rand)==0 & length(check_heter_grp_fixed)==0){

          #### Add the environment to the fixed term by default if not provided by the user
          ## Providing it in fixed term allow for optimal model fit compared to the random term

          # Adding fixed factors
          if (!is.null(fixed_term)) {


            ### Update the fixed term if not null
            fixed <-   stats::update(fixed,
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

          for (i in 1:length(names_in_inv_list)){
            random <- stats::update(random,
                                    paste("~ . +", build_structured_vm_term(i)))

          } ### End

        } ### when Env is missing in both fixed and random terms


        ## If user provide ENV in the fixed term and missing in the random term
        if (length(check_heter_grp_rand)==0 & length(check_heter_grp_fixed)==1){

          for (i in 1:length(names_in_inv_list)){
            random <- stats::update(random,
                                    paste("~ . +", build_structured_vm_term(i)))

          } ### End

        } ### End  If user provide ENV in the fixed term and missing in the random term

        ## If user provide ENV in the random term and missing in the fixed term
        if (length(check_heter_grp_rand)==1 & length(check_heter_grp_fixed)==0){
          for (i in 1:length(names_in_inv_list)){
            random <- stats::update(random,
                                    paste("~ . +", build_structured_vm_term(i)))

          }
        }

      }
      ### If user provide GID:Env
      ### use this step to drop the orginal GID:Env
      if(!is.null(inter_gen_pos) & length(gen_pos)==0){
        rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]
        rand_termCopy <-  gsub(" ", "", rand_termCopy)
        #inter_gen_pos_copy = match(rand_term[inter_gen_pos], rand_termCopy)
        if(!is.null(heter_groups)){
          inter_gen_pos_copy <-  match(paste(gen_name,heter_groups, sep = ":"),rand_termCopy)
          if(anyNA(inter_gen_pos_copy)){
            inter_gen_pos_copy <-  match(paste(heter_groups,gen_name, sep = ":"),rand_termCopy)
          }

        }
        random <-  stats::formula(stats::drop.terms(stats::terms(random),inter_gen_pos_copy, keep.response = F))
      } ### End

      ### If user provide GID, GID:Env
      if(!is.null(inter_gen_pos) & !is.null(gen_pos)){
        ### use this step to drop the original GID and GID:Env
        rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]
        rand_termCopy <-  gsub(" ", "", rand_termCopy)
        #inter_gen_pos_copy = match(rand_term[inter_gen_pos], rand_termCopy)
        if(!is.null(heter_groups)){
          inter_gen_pos_copy <-  match(paste(gen_name,heter_groups, sep = ":"),rand_termCopy)
          if(anyNA(inter_gen_pos_copy)){
            inter_gen_pos_copy <-  match(paste(heter_groups,gen_name, sep = ":"),rand_termCopy)
          }
        }
        gen_pos_copy <-  match(rand_term[gen_pos], rand_termCopy)
        random <-  stats::formula(stats::drop.terms(stats::terms(random), c(gen_pos_copy,inter_gen_pos_copy), keep.response = F))
      }

      # End of when only gen_pos and inter_gen_pos are provided
      ## That is user provide GID, GID:Env

    }

    ##################################################
    ### If user provide only the gen_name in the random term as gen_name effect
    ######################################################
    if(is.null(inter_gen_pos) & length(gen_pos)==1){
      #if (length(gen_pos)==1 & length(inter_gen_pos)==0){
      ### stats::formula must be at least  length 1 so since it just one length. The
      ## adjusted random term was added and the old one was removed
      if (exists('names_in_inv_list')){
        for (i in 1:length(names_in_inv_list)) {
          random <- stats::update(random, paste("~ . +",paste(paste0('vm(', gen_name), sep = ',', paste(names_in_inv_list[i], ")", sep = ""))))

        }
      }
      ### use this step to drop the orginal GID and GID:Env
      rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]
      rand_termCopy <-  gsub(" ", "", rand_termCopy)
      gen_pos_copy <-  match(rand_term[gen_pos], rand_termCopy)
      random <-  stats::formula(stats::drop.terms(stats::terms(random), gen_pos_copy, keep.response = F))
    }
    ######################################################################
    #### When Variance-Covariance Structure is missing. Compound Symmetry
    ####################################################################
    if(is.null(var_cov_str) & (length(pheno_data[[gen_name]])>length(unique(pheno_data[[gen_name]])))){

      ###############################
      if (!is.null(fixed_term)){
        check_heter_grp_fixed  <-  match(heter_groups, fixed_term)
        if(anyNA(check_heter_grp_fixed)) {check_heter_grp_fixed <- NULL}
      } else {
        check_heter_grp_fixed <-  NULL
      }


      if (!is.null(rand_term)){
        check_heter_grp_rand  <-  match(heter_groups, rand_term)
        if(anyNA(check_heter_grp_rand)) {check_heter_grp_rand <- NULL}
      } else {
        check_heter_grp_rand <-  NULL
      }

      #if(var_cov_str=="us" |var_cov_str=="corgh" | var_cov_str=="corh" | var_cov_str=="corv"){

      if(!is.null(inter_gen_pos)){
        if (length(check_heter_grp_rand)==0 & length(check_heter_grp_fixed)==0){
          #### Add the environment to the fixed term by default if not provided by the user
          ## Providing it in fixed term allow for optimal model fit compared to the random term

          # Adding fixed factors
          if (!is.null(fixed_term )) {
            ### Update the fixed term if not null
            fixed <-   stats::update(fixed,
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

          for (i in 1:length(names_in_inv_list)){

            random <-  stats::update(random,
                                  paste("~ . +", (paste(paste0("idv", paste0("(",heter_groups,")")),
                                                        paste(paste0('vm(', gen_name), sep = ',', paste(names_in_inv_list[i], ")", sep = "")),
                                                        sep = ":"))))

          }

        }### when Env is missing in both fixed and random terms

        ## If user provide ENV in the fixed term and missing in the random term
        if (length(check_heter_grp_fixed)==1 & length(check_heter_grp_rand)==0){

          for (i in 1:length(names_in_inv_list)){

            random <- stats::update(random,
                                  paste("~ . +",paste(paste0("idv", paste0("(",heter_groups,")")),
                                                      paste(paste0('vm(', gen_name), sep = ',', paste(names_in_inv_list[i], ")", sep = "")), sep = ":")))


          }

        }### End  If user provide ENV in the fixed term and missing in the random term

        ## If user provide ENV in the random term and missing in the fixed term
        if (length(check_heter_grp_fixed)==0 & length(check_heter_grp_rand)==1){

          if (exists('names_in_inv_list')){inter_gen_pos_use = length(names_in_inv_list)}
          for (i in 1:length(names_in_inv_list)){

            random <-  stats::update(random,
                                  paste("~ . +",paste(paste0("idv", paste0("(",heter_groups,")")),
                                                      paste(paste0('vm(', gen_name), sep = ',', paste(names_in_inv_list[i], ")", sep = "")), sep = ":")))

          }

        }

      }

      ### use this step to drop the orginal GID and GID:Env

      ### If user provide GID, GID:Env
      if(!is.null(inter_gen_pos) & !is.null(gen_pos)){
        ### use this step to drop the original GID and GID:Env
        rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]

        rand_termCopy <-  gsub(" ", "", rand_termCopy)

        #inter_gen_pos_copy = match(rand_term[inter_gen_pos], rand_termCopy)
        if(!is.null(heter_groups)){
          inter_gen_pos_copy <-  match(paste(gen_name,heter_groups, sep = ":"),rand_termCopy)
          if(anyNA(inter_gen_pos_copy)){
            inter_gen_pos_copy <-  match(paste(heter_groups,gen_name, sep = ":"),rand_termCopy)
          }

        }
        gen_pos_copy <-  match(rand_term[gen_pos], rand_termCopy)

        random <-  stats::formula(stats::drop.terms(stats::terms(random), c(gen_pos_copy,inter_gen_pos_copy), keep.response = F))

      }

    } ## End of when variance-covariance str is not provided by user

    #####################################################
    ####
    # Extract each random component/terms
    ####################################################

    ranTerms <- strsplit(as.character(random[2]), split = "[+]")[[1]]
    ranTerms <-  gsub(" ", "", ranTerms)

    for (r in 1:length(ranTerms)) {

      if(r == 1) {
        code_asr[2] <- paste(code_asr[2], ranTerms[r], sep='')

      } else {

        code_asr[2] <- paste(code_asr[2], ranTerms[r], sep='+')
      }
    }


  }

  return(list(code_asr = code_asr,
              gen_pos = gen_pos,
              inter_gen_pos = inter_gen_pos,
              rand_term = rand_term))


}
