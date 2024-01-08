# trains a random forest model on a given set of rows and predicts on a disjunct set of rows
# predictors: array of predictor column names
train_test_by_fold_scaled <- function(df,
                               idx_train,
                               idx_val,
                               predictors = c("elv", "mat", "map", "ndep", "mai", "Species")
                               ){

  # dataframe with scaled training predictors
  # needs to be computed on each fold seperately to avoid data leakage
  pred_train_scaled <- scale(df[idx_train, predictors[-6]])                # omit Species column (cannot be scaled)
  pred_train_scaled <- cbind(pred_train_scaled, df[idx_train, "Species"])  # add unscaled Species column

  # dataframe with scaled validation predictors
  pred_val_scaled <- scale(df[idx_val, predictors[-6]])
  pred_val_scaled <- cbind(pred_val_scaled, df[idx_val, "Species"])

  # build a model on the spatial train folds
  mod <- ranger::ranger(
    x =  pred_train_scaled,  # data frame with columns corresponding to predictors
    y =  df[idx_train,]$leafN,       # a vector of the target leafN values

    # use same hyperparameters as above
    mtry = 3,
    min.node.size = 12,
    splitrule = "variance", # default
    seed = 42)

  # df with only validation data
  val_fold <- df[idx_val,]

  # add predictions to the train fold
  val_fold$fitted <- predict(mod,                         # the fitted model object
                             data = pred_val_scaled       # predictor data frame (validation data)
                             )$predictions                # extract predicted values

  # get evaluation metrics for validation set
  metrics <- val_fold |>
    yardstick::metrics(leafN, fitted)

  # extract r squared and rmse from metrics table
  rsq <- metrics |>
    filter(.metric == "rsq") |>
    pull(.estimate)
  rmse <- metrics |>
    filter(.metric == "rmse") |>
    pull(.estimate)

  return(tibble(rsq = rsq, rmse = rmse))
}
