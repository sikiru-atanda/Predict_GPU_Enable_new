#' Title
#' Overall, this serve as gateway between phenotype-precheck function and readiness of
#' the phenotypic data for model fitting
#' The objective of this function is to do the following:
#' 1. Check the output from phenotype-precheck function for pheno_data before declaring it for model fit.
#'    NA is allowed for response variable with the assumption that those individuals are the testing set.
#' 2. Check output from phenotype-precheck function for pheno_data_train and
#'    check for NA. If NA present it will not be declare for model fit because
#'    NA is not expected in training set. If all good it will be declared for model fit
#' 3. Check the output from phenotype-precheck function for pheno_data_test before declaring it for model fit
#'    NA is allowed for response variable here.
#' 4. If pheno_data_train and pheno_data_test were present and pass through the pre-check process,
#'    it will be processed that is.
#'    1) it check that column name for both data match/the same
#'    2) Combined the dataset for model fit and prediction.
#'    It is assumed here that pheno_data is missing/not provided by the user.
#' 5. The output will be pheno_data declared for model fit and test_set if pheno_data_test was provided.
#'
#' @param pheno_data phenotypic object NA is allowed
#' @param pheno_data_train phenotypic object for the training set. NA is allowed.
#' This is expected if pheno_data was not provided
#' @param pheno_data_test phenotypic object for the testing set. NA is allowed
#' @param response y variables/lables
#' @param gen_name column name containing individuals/genotypes
#' @param response_family Response distribution family
#'   (`"gaussian"`, `"binomial"`, `"multinomial"`, `"ordinal"`).
#' @param heter_groups Optional column carrying the heterogeneous-group label
#'   (typically environment) for MET runs.
#' @param random Optional one-sided formula for random terms (e.g. `~ GID`).
#' @param fixed Optional one-sided formula for additional fixed terms.
#' @param type_pheno Optional phenotype-type hint
#'   (e.g. `"continuous"`, `"binary"`).
#' @param train_set Optional integer / character vector of training-set row
#'   indices or genotype IDs.
#' @param test_set Optional integer / character vector of test-set row indices
#'   or genotype IDs (rows with `NA` response are inferred as test if not
#'   supplied).
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return The phenotype frame declared for model fitting (and prediction),
#'   plus the resolved `test_set` indices.
#' @export
#'
#' @examples

phenotype_to_model <- function(
                              pheno_data = NULL,
                              pheno_data_train = NULL,
                              pheno_data_test = NULL,
                              train_set = NULL,
                              test_set = NULL,
                              response=NULL,
                              response_family = "gaussian",
                              gen_name=NULL,
                              heter_groups = NULL,
                              random = NULL,
                              fixed = NULL,
                              type_pheno = NULL,
                              ...
                          ) {

  pheno_data_train_ <-  NULL

  pheno_data_test_ <- NULL
  #test_set <-  NULL

msg <- ""
normalize_test_id_vector <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.data.frame(x) || is.matrix(x)) {
    x <- x[, 1]
  } else if (is.list(x)) {
    stop(paste(msg,'The testing set cannot be a list. Should be either dataframe, matrix or a vector.'), call. = FALSE)
  }
  x <- unique(as.character(x))
  x[!is.na(x) & nzchar(x)]
}

build_test_set_by_trait <- function(pheno_df, response_cols, id_col) {
  out <- setNames(vector("list", length(response_cols)), response_cols)
  for (col in response_cols) {
    out[[col]] <- unique(as.character(pheno_df[[id_col]][is.na(pheno_df[[col]])]))
  }
  out
}

explicit_test_set_supplied <- !is.null(test_set)
explicit_train_set_supplied <- !is.null(train_set)
test_set_source <- NULL
test_set_by_trait <- NULL
test_set_by_trait_source <- NULL

  ## Check availability of pheno_datatypic data (training and testing set). This
  ## accommodate missing value with the assumption that testing will have NA
  ## Check for NA is not test here for response variables because testing set might have NA
  # thus, NA is expected in pheno_data.

  if ((!is.null(pheno_data) & (is.null(pheno_data_train) & is.null(pheno_data_test)))){

    pheno_data <- phenotype_precheck(pheno_data= pheno_data,
                                     gen_name = gen_name,
                                     response = response,
                                     response_family = response_family,
                                     heter_groups = heter_groups,
                                     random = random,
                                     fixed = fixed)

    if (isTRUE(explicit_test_set_supplied)) {
      test_set <- normalize_test_id_vector(test_set)
      test_set_source <- "explicit"
    } else if (isTRUE(explicit_train_set_supplied)) {
      train_set <- normalize_test_id_vector(train_set)
      non_train_rows <- !as.character(pheno_data[[gen_name]]) %in% as.character(train_set)
      test_set <- unique(as.character(pheno_data[[gen_name]][non_train_rows]))
      test_set_source <- "train_set"

      if(length(test_set)==0){

        message(paste( insight::print_color("WARNINGS\n", "blue"),
                       insight::print_color(paste(msg,paste("The training set size is the same size as the unique genotypes in the pheno_data.")), "blue")))

      }
    } else {
      # Find the rows with NA in each response column if the user has NA as testing set
      na_rows <- lapply(response, function(col) which(is.na(pheno_data[[col]])))
      has_missing_response <- any(vapply(na_rows, length, integer(1)) > 0L)

      if (isTRUE(has_missing_response)) {
        same_na_pattern <- all(sapply(na_rows, function(x) identical(x, na_rows[[1]])))
        test_set_by_trait <- build_test_set_by_trait(pheno_data, response, gen_name)
        test_set_by_trait_source <- "inferred_missing_response"

        if (!isTRUE(same_na_pattern)) {
          warning(
            paste(
              msg,
              sprintf(
                "Rows containing NA did not match across response columns: %s. Proceeding with trait-specific inferred testing sets.",
                paste(response, collapse = ", ")
              )
            ),
            call. = FALSE
          )
        }
      } else {
        test_set <- NULL
      }
    }

    if(!is.null(test_set) && length(test_set) > 0L){

      # pheno_data[[response]] <- ifelse(pheno_data[[gen_name]]%in%test_set, NA,
      #                                  pheno_data[[response]])
      # Loop through each column name in response and apply the ifelse function

      if (dplyr::is_grouped_df(pheno_data)) {
        pheno_data <- dplyr::ungroup(pheno_data)
      }

      test_rows <- as.character(pheno_data[[gen_name]]) %in% test_set
      for (col in response) {
        pheno_data[[col]][test_rows] <- NA
      }

      # test_set_ = test_set
      #
      # rm(test_set)
    }

    ### The result object has to pass the test attribute before it can be stored/
    ## pass through for the next step of check and declared good for model fit
    if(attr(pheno_data, "cleared")=="pass") {

    # Assign appropriate class.
    #class(pheno_data) <- c("data.frame", "phenotype")

    attr(pheno_data, "cleared") <- "for_model_fit"

    output <- list(pheno_clean_data = pheno_data)
    if (!is.null(test_set_by_trait)) {
      output[["test_set_by_trait"]] <- test_set_by_trait
    }
    if (!is.null(test_set_by_trait_source)) {
      output[["test_set_by_trait_source"]] <- test_set_by_trait_source
    }

    #names(output) <- c("pheno_data")

    #rm(pheno_data)

    } else {

      stop(paste(msg, "Data is not class phenotype."), call. = FALSE)
    }

  } else {

    ### If pheno_data is missing, pheno_data traning should be available for traning and it is expected
    ## pheno_data testing is also available
    ## However, pheno_data testing might be missing if the user objective is for cross-validation
    ## In general, missing value is not expected in the pheno_data training.
    if (!is.null(pheno_data_train)){

      pheno_data_train_ <- phenotype_precheck(pheno_data= pheno_data_train,
                                              gen_name = gen_name,
                                              response = response,
                                              response_family = response_family,
                                              heter_groups = heter_groups,
                                              random = random,
                                              fixed = fixed)

      if(attr(pheno_data_train_, "cleared")!="pass") {

        stop(paste(msg, 'pheno_data_train is not object phenotype.'), call. = FALSE)



      }

      ## Though pass the pre-check test but NA is not expected in the pheno_data training.
      if(anyNA(pheno_data_train_)){

        stop(paste(msg, "Missing value in not accepted in training set."), call. = FALSE)

        }

      rm(pheno_data_train)
    }

    ### pheno_data test can be present or absent. It is expected this will contain
    ## missing value for the response variable

    if (!is.null(pheno_data_test)){

      pheno_data_test_ <- phenotype_precheck(pheno_data= pheno_data_test,
                                             gen_name = gen_name,
                                             response = response,
                                             response_family = response_family,
                                             type_pheno = type_pheno)

      #if(attr(pheno_data_test_, "cleared")!="pass" && all(class(pheno_data_test_)!=c("data.frame", "phenotype"))) {
      if(attr(pheno_data_test_, "cleared")!="pass") {
        stop(paste(msg, 'pheno_data_test is not object phenotype.'), call. = FALSE)

      }

      rm(pheno_data_test)
    }

    #### if both pheno_data_train and pheno_data_test are provided

    if (!is.null(pheno_data_train_) & !is.null(pheno_data_test_)){

      if (!identical(colnames(pheno_data_train_), colnames(pheno_data_test_))){
        stop(paste(msg,'Columns name in the pheno_data_train not the same as pheno_data_test.'), call. = FALSE)

      } else{

        pheno_data <- rbind(pheno_data_train_, pheno_data_test_)

        test_set <- data.frame(name = as.character(unique(pheno_data_test_[, gen_name])), stringsAsFactors = FALSE)
        names(test_set) <- gen_name
        test_set_by_trait <- setNames(
          replicate(length(response), unique(as.character(pheno_data_test_[, gen_name])), simplify = FALSE),
          response
        )
        test_set_source <- "pheno_data_test"
        test_set_by_trait_source <- "pheno_data_test"


        ### The result object has to pass the test attribute before it can be stored

          # Assign appropriate class.
          #class(pheno_data) <- c("data.frame", "phenotype")

          attr(pheno_data, "cleared") <- "for_model_fit"



      }

    } else {
      if (!is.null(pheno_data_train_) & is.null(pheno_data_test_)){
        pheno_data <- pheno_data_train_
        attr(pheno_data, "cleared") <- "for_model_fit"

      }
    }



  }



  #if(!is.null(pheno_data) & (!is.null(test_set) && !is.null(train_set))){

  if(!is.null(pheno_data) & !is.null(test_set)){

    output <- list(pheno_clean_data = pheno_data,
                   test_set = test_set)
    if (!is.null(test_set_source)) {
      output[["test_set_source"]] <- test_set_source
    }
    if (exists("test_set_by_trait") && !is.null(test_set_by_trait)) {
      output[["test_set_by_trait"]] <- test_set_by_trait
    }
    if (exists("test_set_by_trait_source") && !is.null(test_set_by_trait_source)) {
      output[["test_set_by_trait_source"]] <- test_set_by_trait_source
    }

    #rm(pheno_data, test_set)

    }else {

      if(!is.null(pheno_data) & is.null(test_set)){
        ## here it is expected only training set is provided so no testing set, thus no NA
        ## Recheck to be sure no NA

        na_rows <- lapply(response, function(col) which(is.na(pheno_data[[col]])))
        if (!(exists("test_set_by_trait") && !is.null(test_set_by_trait)) && length(na_rows[[1]]) > 0) {
          stop(paste(msg, "Missing value is not expected in the training set."), call. = FALSE)
        }
        output <- list(pheno_clean_data = pheno_data)
        if (exists("test_set_by_trait") && !is.null(test_set_by_trait)) {
          output[["test_set_by_trait"]] <- test_set_by_trait
        }
        if (exists("test_set_by_trait_source") && !is.null(test_set_by_trait_source)) {
          output[["test_set_by_trait_source"]] <- test_set_by_trait_source
        }
        #rm(pheno_data, test_set)

      }

  }

  return(output)
}


