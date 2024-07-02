
# Function to print messages with optional warning, error, or info types and stopping behavior, with updated ordering
print_custom_message <- function(msg, type = "info", should_stop = FALSE) {

  # Define a separator line
  separator <- rep("=", 50)  # Create a line of 50 equal signs
  separator_line <- paste(separator, collapse="")

  # Define message color based on the type
  message_color <- switch(type,
                          "info" = crayon::blue,
                          "warning" = crayon::yellow,
                          "error" = crayon::red,
                          crayon::blue) # Default to blue if type is unspecified

  # Optionally stop the execution if should_stop is TRUE
  if (should_stop && type == "error") {
    stop_message <- "Stopping execution due to critical error."
    cat(crayon::red(stop_message), "\n", separator_line, "\n")
    stop(crayon::red(msg), call. = FALSE)
  } else {
    # Print the separator, then the message in the specified color
    cat(separator_line, "\n", message_color(msg), "\n", separator_line, "\n")
  }
}

# # Example usage of the function
# print_custom_message("This is an informational message.", "info")
# print_custom_message("This is a warning message. Please take note.", "warning")
# print_custom_message("This is an error message and will stop execution.", "error", should_stop = TRUE)


# # Function to print messages in specified colors with an integrated separator
# print_colored_message <- function(msg, color_name) {
#
#   # Define a separator
#   separator <- rep("=", 50)  # Create a line of 50 equal signs
#   separator_line <- paste(separator, collapse="")
#
#   # Create a list of available color functions from crayon using explicit reference
#   colors <- list(
#     blue = crayon::blue,
#     red = crayon::red,
#     green = crayon::green,
#     yellow = crayon::yellow,
#     magenta = crayon::magenta,
#     cyan = crayon::cyan
#   )
#
#   # Check if the specified color is available
#   if (!color_name %in% names(colors)) {
#     stop("Specified color is not supported.")
#   }
#
#   # Print the message with separators
#   cat(crayon::yellow(separator_line), "\n", colors[[color_name]](msg), "\n", crayon::yellow(separator_line), "\n")
# }
#
# # Example usage of the function
# #print_colored_message("Number of iteration is missing. Default value of 26000 was assigned.\nCheck if this is appropriate for your data.", "blue")
