# .onAttach <- function(libname, pkgname) {
#   # Path to the logo file
#   logo_path <- system.file("NDSUlogo.png", package = pkgname)
#
#   # Display the logo and custom message
#   packageStartupMessage("This is a product of sik
# ----------------------------------------
# ")
#   if (file.exists(logo_path)) {
#     logo <- png::readPNG(logo_path)
#     grid::grid.raster(logo)
#   }
#
#   packageStartupMessage("
# Welcome to North Dakota State University R Package!
# ----------------------------------------
# ")
# }


.onAttach <- function(libname, pkgname) {
  # ANSI escape codes for colors
  green <- "\033[32m"
  yellow <- "\033[33m"
  reset <- "\033[0m"

  # ASCII art of the NDSU logo with colors
  logo_ascii <- paste0(
    green, " _   _  ____  ____  _    _  ", reset, "\n",
    green, "| \\ | ||  _ \\|  _ \\| |  | | ", reset, "\n",
    green, "|  \\| || | | | | | | |  | | ", reset, "\n",
    green, "| . ` || |_| | |_| | |__| | ", reset, "\n",
    green, "|_|\\_||____/|____/ \\____/  ", reset, "\n",
    yellow, "  North Dakota State University", reset, "\n"
  )

  # Display the ASCII logo and custom message
  packageStartupMessage("
PredictProR Product of NDSU!
----------------------------------------
")
  packageStartupMessage(logo_ascii)

}
