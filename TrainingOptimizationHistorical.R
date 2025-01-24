library(caret)
library(ipred)
library(xgboost)
library(randomForest)
#####
##Data prep
rm(list = ls()); ls()
if(!is.null(dev.list())) dev.off()
cat("\014")
gc()

library(data.table)


setwd("D:/FILES PILE UP/2017_2018_Data")

rm(list = ls()); ls()
if(!is.null(dev.list())) dev.off()
cat("\014")
gc()

library(data.table)


#setwd("/home/asa259/2017_2018AVG")

#setwd("F:/FILES PILE UP/2017_2018_Data")

#setwd("D:/PredictProR")
ref = read.csv("ref2017_2018DanMatchPerTrial.csv", header = T, as.is = T, sep = ",", stringsAsFactors = F)


ref = subset(ref, select = c("GID", "GY_BLUE", "popClass", "Year"))

ref$GID = as.character(ref$GID)

ref = ref[order(ref$GID),]

table(ref$popClass)

#setwd("D:/FILES PILE UP/2017_2018_Data")
mark = read.csv("2017_2018Kenya_Dan0.05.csv", header = T, as.is = T,
                stringsAsFactors = F, row.names = 1)

size = as.data.frame(table(ref$popClass))

colnames(size) = c("popClass", "size")

setdiff(rownames(mark), ref$GID)

monomorphic_markers <- which(apply(mark, 2, function(x) length(table(x)) <= 1))
if (length(monomorphic_markers) > 0) {

  mark_clean <- mark[, -monomorphic_markers]
}

#mark2 = detect_genomic_coding(as.matrix(mark))



ref2017 = ref[ref$Year==2017, ]

table(ref2017$popClass)

ref2018 = ref[ref$Year==2018, ]

#table(ref2018$popClass)

ref2017_pop5 = ref2017[ref2017$popClass=="pop5@2017", "GID"]


ref2017_pop3 = ref2017[ref2017$popClass=="pop3@2017", "GID"]

ref2017_pop2 = ref2017[ref2017$popClass=="pop2@2017", "GID"]

ref2017_pop4 = ref2017[ref2017$popClass=="pop4@2017", "GID"]

ref2017_pop13 = ref2017[ref2017$popClass=="pop13@2017", "GID"]

ref2017_pop12 = ref2017[ref2017$popClass=="pop12@2017", "GID"]

ref2017_pop5 = ref2017[ref2017$popClass=="pop5@2017", "GID"]

ref2017_pop8 = ref2017[ref2017$popClass=="pop8@2017", "GID"]

ref2017_pop9 = ref2017[ref2017$popClass=="pop9@2017", "GID"]
######
ref2018_pop20 = ref2018[ref2018$popClass=="pop20@2018", "GID"]

ref2018_pop3 = ref2018[ref2018$popClass=="pop3@2018", "GID"]

ref2018_pop5 = ref2018[ref2018$popClass=="pop5@2018", "GID"]

ref2018_pop6 = ref2018[ref2018$popClass=="pop6@2018", "GID"]

ref2018_pop41 = ref2018[ref2018$popClass=="pop41@2018", "GID"]

ref2018_pop7 = ref2018[ref2018$popClass=="pop7@2018", "GID"]

ref2018_pop39 = ref2018[ref2018$popClass=="pop39@2018", "GID"]

ref2018_pop10 = ref2018[ref2018$popClass=="pop10@2018", "GID"]
###
tst = c(ref2017_pop3) #ref2017_pop3

trn = c(ref2017_pop4, ref2017_pop5,
        ref2017_pop13, ref2017_pop12,
        ref2017_pop8, ref2017_pop9,
        ref2018_pop20, ref2018_pop3,
        ref2018_pop5, ref2018_pop6,
        ref2018_pop41,ref2018_pop7,
        ref2018_pop39, ref2018_pop10)

trn_tst = c(ref2017_pop2, ref2017_pop3,
            ref2017_pop4, ref2017_pop5,
            ref2018_pop20, ref2018_pop3,
            ref2018_pop5)

y_trn = ref[ref$GID%in%trn, "GY_BLUE"]

geno_trn = mark[rownames(mark)%in%trn, ]


y_tst = ref[ref$GID%in%tst, "GY_BLUE"]

geno_tst = mark[rownames(mark)%in%tst, ]

trn_ref= ref[ref$GID%in%trn, ]

trn_ref_opt= trn_ref[train_indices, ]
###### Data prep end
opt = ref[ref$GID%in%rownames(X_train_selected), ]

opt2 = ref[ref$GID%in%fTS, ]
table(opt$popClass)

historical_data <- cbind(y=y_trn, geno_trn)

# Initialize training set with a small subset of the historical data
initial_indices <- sample(1:nrow(historical_data), size = 50)  # Start with 50 samples
train_indices <- initial_indices
remaining_indices <- setdiff(1:nrow(historical_data), train_indices)

# Set the expected reliability threshold and maximum number of iterations
expected_reliability_threshold <- 0.80  # Define your reliability threshold
max_iterations <- 20  # Maximum number of iterations
max_train_size <- 200  # Maximum size of the training set

train_predict_xgboost <- function(data_label_geno, indices, test_geno, params, nrounds) {
  train_data <- data_label_geno[, -1]
  y_train <- data_label_geno[, 1]
  # Subset the data
  dtrain_boot <- xgboost::xgb.DMatrix(data = as.matrix(train_data[indices,]), label = y_train[indices])

  # Train the model
  model <- xgboost::xgboost(data = dtrain_boot, params = params, nrounds = nrounds, verbose = 0)

  # Predict on the original data

  pred <- stats::predict(model, as.matrix(test_geno), reshape = TRUE)


  return(pred)
}

# Get the best-tuned parameters
params <- list(booster = "dart",
               objective = "reg:squarederror",
               sample_type = "uniform",
               normalize_type = "tree",
               rate_drop = 0.1,
               skip_drop = 0.5,
               eta = 0.1,
               max_depth = 6,
               min_child_weight = 1,
               subsample = 0.8,
               colsample_bytree = 0.8,
               lambda = 1,
               alpha = 0)
# Function to train models and estimate uncertainty
train_and_estimate_uncertainty <- function(train_data, test_data,
                                           params = NULL,
                                           nrounds = 100, verbose = 0) {
  # Train bagging model

  boot_results <- boot::boot(
    data = train_data,
    statistic = train_predict_xgboost,
    R = 100,  # Number of bootstrap samples
    #sim = "ordinary",
    #dtrain = geno_omic_object_train,
    test_geno = test_data,
    params = params,
    nrounds = nrounds
  )

  pred_variances <- apply(boot_results$t, 2, var)

  pred_SE <- apply(boot_results$t, 2, sd)

  predmean <- apply(boot_results$t, 2, mean)

  lower_bound <- apply(boot_results$t, 2, quantile, probs = 0.05)
  upper_bound <- apply(boot_results$t, 2, quantile, probs = 0.95)
  interval_width <- upper_bound - lower_bound

  # Calculate reliability (e.g., 1 - (mean interval width / mean prediction))
  #interval_width <- prediction_intervals$upper_bound - prediction_intervals$lower_bound
  reliability <- 1 - (pred_variances / var(predmean))


  return(list(
    prediction_intervals = interval_width,
    sd_predictions = pred_SE,
    reliability = mean(reliability),
    predmean = predmean,
    rel = reliability
  ))
}

# Iterative process to dynamically update the training set
for (iteration in 1:max_iterations) {
  # Create the current training and validation sets
  train_data <- historical_data[train_indices, ]
  #test_data <- geno_tst  # Assuming you have a validation set

  # Train models and estimate uncertainty
  results <- train_and_estimate_uncertainty(train_data = train_data,
                                            test_data = geno_tst, params = params,
                                            nrounds = 100, verbose = 0)

  # Check if reliability meets the expected threshold
  if (results$reliability >= expected_reliability_threshold) {
    cat("Expected reliability threshold met. Stopping iterations.\n")
    break
  }

  # Select additional data points with the highest uncertainty
  uncertainty <- results$sd_predictions
  uncertain_indices <- order(uncertainty, decreasing = TRUE)[1:10]  # Select top 10 uncertain points

  # Update training and remaining indices
  new_train_indices <- remaining_indices[uncertain_indices]
  train_indices <- c(train_indices, new_train_indices)
  remaining_indices <- setdiff(remaining_indices, new_train_indices)

  # Optionally, drop the least informative data points from the training set
  # Here, we remove the 10 points that contribute the least to reducing uncertainty if train size exceeds max size
  if (length(train_indices) > max_train_size) {
    # Identify least informative points (you can use different criteria)
    train_uncertainty <- results$sd_predictions[train_indices]
    least_informative_indices <- train_indices[order(train_uncertainty, decreasing = FALSE)[1:10]]
    train_indices <- setdiff(train_indices, least_informative_indices)
  }

  # Check stopping criteria (e.g., if no more data to add or desired accuracy achieved)
  if (length(remaining_indices) == 0) break
}

# Final training set
final_train_data <- historical_data[train_indices, ]

# Final model training on the selected optimal training set
final_results <- train_and_estimate_uncertainty(train_data = final_train_data,
                                                test_data = geno_tst,
                                                params = params)

# Display final prediction intervals and reliability
print(final_results$prediction_intervals)
cat("Final reliability:", final_results$reliability, "\n")

# Check if the final reliability meets the expected threshold
if (final_results$reliability < expected_reliability_threshold) {
  cat("The historical data did not contain individuals with high reliability for the test set.\n")
  cat("The test set might need to be planted for better prediction accuracy.\n")
}

cor(final_results$predmean,y_tst)
View(final_results$prediction_intervals)

View(data.frame(y_tst))
