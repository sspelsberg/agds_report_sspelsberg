# run model and compare to true values
# returns the RMSE
rmse_gdd <- function(par, data) {

  # split out data
  drivers <- data$drivers
  validation <- data$validation

  # calculate phenology predictions
  # and put in a data frame
  predictions <- drivers |>
    group_by(year) |>
    summarise(
      predictions = gdd_model(
        temp = tmean,
        par = par
      )
    )

  predictions <- left_join(predictions, validation, by = "year")

  rmse <- predictions |>
    summarise(
      rmse = sqrt(mean((predictions - doy)^2, na.rm = TRUE))
    ) |>
    pull(rmse)

  # return rmse value
  return(rmse)
}
