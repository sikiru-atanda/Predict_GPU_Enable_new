# Function to revert scaling of labels
revert_scaling <- function(scaled_values, scaler) {
  center <- scaler$mean
  scale <- scaler$std
  original_values <- scaled_values * scale + center
  return(original_values)
}
