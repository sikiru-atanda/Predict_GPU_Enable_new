

genetic_space_recommendation <- function(recommendation) {

  msg <- ""

  if (recommendation=="All metrics indicate low coverage. Consider expanding the training set.") {
    message(paste(msg,"All the metrics recommendation for determining the power of the training data to predict the testing set is not optimal."))
    message(paste(msg, "Press [Enter] to continue the analysis anyway, or type [stop] to terminate."))

    # Capture user input
    if (interactive()) {
      user_input <- readline(prompt = "Your choice: ")
      if (tolower(user_input) == "stop") {

        stop(paste(msg, "Analysis terminated by the user."), call. = FALSE)
      } else {
        message(paste(msg, "Continuing with the analysis..."))
      }
    } else {
      message(paste(msg, "Non-interactive session: automatically continuing analysis..."))

    }
  } else if (recommendation=="All coverage metrics are good. The training set is likely adequate."){
    message(paste(msg, "All coverage metrics are good. The training set is likely adequate."))

  } else{
    message(paste(msg, recommendation))
  }


}




