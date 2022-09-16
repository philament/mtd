
#' Plots a Modified Target Diagram (MTD)
#'
#' @description  This function creates a plot of a Modified Target Diagram according to Yatkin et al. (2022). Such diagrams are useful for comparing, for example,
#' a low-cost air quality sensor against reference instrumentation, particularly with regard to the Data Quality Objectives of the European air quality directive.
#'
#' @references Yatkin, S., Gerboles, M., Borowiak, A., Davila, S., Spinelle, L., Bartonova, A., Dauge, F., Schneider, P., Van Poppel, M., Peters, J.,
#' Matheeussen, C., & Signorini, M. (2022). Modified Target Diagram to check compliance of low-cost sensors with the Data Quality Objectives of the
#' European air quality directive. Atmospheric Environment, 273, 118967. \href{https://doi.org/10.1016/j.atmosenv.2022.118967}{https://doi.org/10.1016/j.atmosenv.2022.118967}
#'
#' @param pollutant Required. String giving pollutant (one of "CO","CO2","CO2","NO","NO2","O3","PM1","PM10","PM2.5")
#' @param dates Required. Vector of datetimes of the measurements
#' @param reference_data Required. Numeric vector of observations from reference
#' @param sensor_data Required. Numeric vector of observations from sensor
#' @param ubs_rm Required. Random standard uncertainty of the results of the reference method, given as a constant value for all reference values (in same units as measurements)
#' @param ubss Required. Between standard uncertainty of the sensor, given as a constant value for all yis sensor values (in same units as measurements)
#' @param unit_ref Optional. Character string prodiving the unit of the reference values (one of "ppm", "ppb", "ug/m3", "mg/m3"). Default is 'ug/m3'.
#' @param unit_sensor Optional. Character string prodiving the unit of the sensor values (one of "ppm", "ppb", "ug/m3", "mg/m3"). Default is 'ug/m3'.
#' @param x_axis_label,y_axis_label Optional. Character string. Labels of the x and y axes.
#' @param x_lim,y_lim Optional. Numeric vectors of two values min and max values. Limits of x and y axis, default values is NA. [DOES NOT WORK CURRENTLY]
#' @param sensor_name Optional. Character string. Name of the sensor to be written in front of the calibration equation. If NULL (default), do not print sensor name.
#'
#' @return A plot
#' @export
#'
#' @examples
#' data(pm25)
#' mtd(pollutant = 'PM2.5',
#'     dates = pm25$datetime,
#'     reference_data = pm25$reference_pm25,
#'     sensor_data = pm25$sensor_pm25,
#'     ubs_rm = 0.5,
#'     ubss = 5)

mtd <- function(pollutant, dates, reference_data, sensor_data, ubs_rm, ubss, unit_ref = 'ug/m3', unit_sensor = 'ug/m3', x_axis_label = NULL, y_axis_label = NULL,
                x_lim = NA, y_lim = NA, sensor_name = NULL)
{
    pollutant = match.arg(pollutant,
                          choices = c("CO","CO2","CO2","NO","NO2","O3","PM1","PM10","PM2.5","RH","RH_int","Temp","Temp_int","Press"))

    unit_ref = match.arg(unit_ref, choices = c('ug/m3', 'ppm', 'ppb', 'mg/m3'))
    unit_sensor = match.arg(unit_sensor, choices = c('ug/m3', 'ppm', 'ppb', 'mg/m3'))

    stopifnot(length(dates) == length(sensor_data))
    stopifnot(length(reference_data) == length(sensor_data))
    n_obs = length(reference_data)

    # Create a data.frame for U_orth_DF
    df = data.frame(case = 1:n_obs,
                    date = dates,
                    xis = reference_data,
                    yis = sensor_data,
                    ubs_rm = ubs_rm,
                    ubs_s = ubss)

    # Run orthogonal regression
    U.orth.List = U_orth_DF(df, Regression = 'OLS', Add.ubss = FALSE)

    # Get the DQO thresholds
    DQO = get.DQO(name.gas = pollutant, unit.ref = unit_ref)

    # Run the target diagram function
    Target.Diagram(Sensor_name = sensor_name,
                   Mat = U.orth.List[["Mat"]],
                   ubsRM = ubs_rm,
                   ubss = ubss,
                   b0 = U.orth.List$b0,
                   b1 = U.orth.List$b1,
                   Unit.Ref = unit_ref,
                   Unit.sensor = unit_sensor,
                   DQO.1 = DQO$DQO.1/DQO$LV,
                   DQO.2 = DQO$DQO.2/DQO$LV,
                   DQO.3 = DQO$DQO.3/DQO$LV,
                   LAT = DQO$LAT,
                   UAT = DQO$UAT,
                   LV = DQO$LV,
                   AT = DQO$AT,
                   sdm_sdo = U.orth.List[["sdm"]] > U.orth.List[["sdo"]],
                   BeginEnd = as.character(dates),
                   with.ubss = FALSE,
                   xAxisLabel = x_axis_label,
                   yAxisLabel = y_axis_label,
                   Xlim = x_lim,
                   Ylim = y_lim
                   )



}
