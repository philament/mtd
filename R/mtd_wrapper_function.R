
#' Plots a Modified Target Diagram (MTD)
#'
#' @description  This function creates a plot of a Modified Target Diagram according to Yatkin et al. (2022). Such diagrams are useful for comparing, for example,
#' a low-cost air quality sensor against reference instrumentation, particularly with regard to the Data Quality Objectives of the European air quality directive.
#' 
#' The relative expanded measurement uncertainty (Ur) is given by the distance between the origin (0, 0) and any point on the bold coloured line read on concentric circles.
#' The color scale indicates the reference data (xis).
#' The extent of random error (Rel.RSS) is read on the x-axis, the extent of relative bias (Rel.bias) is read on the y-axis.
#' Within Rel.bias, the contribution of the slope (b1, thin vertical coloured line) and intercept (b0, thin oblique coloured line) are read on the x-axis.

#' The primary use of the Modified Target Diagram is to check whether the Data Quality Objective (DQO) of indicative methods set in the European Directive (EC/50/2008) are satisfied.
#' DQOs are satisfied when the Limit Value (LV) blue star found on the coloured bold line falls within the 1st target circle representing the DQO.
#' The Modified Target Diagram gives additional information such as: 
#'  . The bold coloured line being above the x-axis (1st or 4th quadrant) indicates an overestimation of the sensor compared to reference data. Conversely, the bold coloured line being below the x-axis (2nd or 3rd quadrant) indicates an underestimation of the sensor compared to reference data.
#'  . Additionally, the bold coloured line being at the right of y-axis (1st or 2nd quadrant) indicates a higher sensitivity of the sensor compared to reference data. Conversely, the bold coloured line being at left of the y-axis (3rd or 4th quadrant) indicates a lower sensitivity of the sensor compared to reference data.
#'  . Variation of contributors to Rel.bias: high contribution from b0 would indicate an offset between the sensor and reference data, possibly correctable by offset subtraction, while high contribution from b1 may indicate an erroneous slope of calibration, possibly correctable by re-calibration or readjustment. Significant contribution from b0 and/or b1 are evidenced when corresponding thin coloured line(s) are far from y-axis.
#'  . The comparison of Rel.RSS and Rel.bias: overwhelming value of Rel.RSS compared to Rel.bias shows that Ur is dominated by the random errors likely resulting from several parameters including the electronic noise of sensor or by mis-calibration as missing covariates ...
#'  . Improvement by adjustment of b0 and/or b1 in order to set the Rel.bias to zero: the adjustment of b0 and b1 could allow to set the Rel.bias to zero for the entire range reference data with a rotation of the bold coloured line.
#'  However, this adjustment of b0 and b1 is only significant if Ur is not dominated by rel.RSS. In this case, setting the Rel.bias to zero by adjustment of b0 and b1 would not allow any significant decrease of Ur.
#' 
#' There are cases when the effects of b0 and b1 never cancel each other and thus Rel.bias nevercrosses the x-axis. Such cases might happen due to either:
#'  . the contributions of b0 and b1 are located on the same side of the y-axis with the combinations of b0 negative and b1 < 1, or, b0 positive and b1 > 1;
#'  . the contributions of b0 and b1 are located on different sides of the y-axis but one of these contributions is overwhelming with the b0 effect being higher than the b1 effect.
#'
#' @references Yatkin, S., Gerboles, M., Borowiak, A., Davila, S., Spinelle, L., Bartonova, A., Dauge, F., Schneider, P., Van Poppel, M., Peters, J.,
#' Matheeussen, C., & Signorini, M. (2022). Modified Target Diagram to check compliance of low-cost sensors with the Data Quality Objectives of the
#' European air quality directive. Atmospheric Environment, 273, 118967. \href{https://doi.org/10.1016/j.atmosenv.2022.118967}{https://doi.org/10.1016/j.atmosenv.2022.118967}
#'
#' @param pollutant Required. String giving pollutant (one of "CH4", CO","CO2","NO","NO2","O3","PM1","PM10","PM2.5")
#' @param dates Required. Vector of datetimes of the measurements. It shall be a POSIXct or date.
#' @param sensor_data Required. Numeric vector of observations from sensor. Unit shall be the same as for reference_data
#' @param Add.ubss Optional, logical, default is FALSE. If TRUE ubss is added to the Random Error (see Yatlin et al.). If FALSE ubss is not added to the Random Error.
#' @param variable.ubss Optional, logical, default is FALSE. If FALSE, ubss is used as constant random standard uncertainties for all sensor_data. If TRUE ubss shall be a numeric vector of same length as sensor_data
#' and is used for each sensor_data in the same order
#' @param ubss Optional, default is NULL, required if Add.ubss is TRUE. Numeric value or numeric vector. Between standard uncertainty of the sensor_data ubss is either a constant value used for all sensor_data or
#' a numeric vector with one value used for each sensor_data. In this case, ubss shall be in the same order as sensor_data.
#' @param unit_sensor Optional. Character string prodiving the unit of the sensor values (one of "ppm", "ppb", "ug/m3", "mg/m3"). Default is 'ug/m3'.
#' @param reference_data Required. Numeric vector of observations from reference analyser with same length as sensor_data. Unit shall be the same as for sensor_data
#' @param variable.ubsRM Optional, logical, default is FALSE. If FALSE, ubsRM is used as constant random standard uncertainties for all xis reference values. 
#' If TRUE ubsRM shall be a numeric vector with one value for each reference-data in the same order.
#' @param ubsRM Required, numeric value or numeric vector. Random standard uncertainty of reference_data. ubsRM is either a constant value used for all reference_data or
#' a numeric vector (if variable.ubsRM is TRUE) with one value used for each reference_data. In this case, ubsRM shall be in the same order as reference_data. ubsRM shall be in same units as reference_data.
#' @param unit_ref Optional. Character string prodiving the unit of the reference values (one of "ppm", "ppb", "ug/m3", "mg/m3"). Default is 'ug/m3'.
#' @param Fitted.RS Optional, logical, default is FALSE. If TRUE the square residuals (RSi) are fitted using a General Additive Model (GAM), provided that the null hypothesis of no correlation between xis and RS is rejected when the probability is lower than 0.05, (p < 0.05)
#' @param Forced.Fitted.RS Optional, logical, default is FALSE. If TRUE even if the variance of residuals is constant (the null hypothesis), RSi are GAM fitted.
#' @param Regression character, default is "OLS", possible values are "OLS" ,"Deming" and "Orthogonal". For Orthogonal regression Delta is 1 and for "Deming" Delta is ubss^2/ubsRM^2. See https://en.wikipedia.org/wiki/Deming_regression

#' @param sensor_name Optional, character string, default is NULL. Name of the sensor to be written in front of the calibration equation. If NULL (default), do not print sensor name.
#' @param Max.percent Optional, numeric in percent, default is NULL. Maximum extent of the x and y axis of the Target Diagram. This is respected provided that there exists Ur smaller than Max.percent. If NULL, DQO.3 is Max.percent .
#' @param x_axis_label,y_axis_label Optional. Character string. Labels of the x and y axes.
#' @param Show.Diag.Ur Optional, logical, default is FALSE. If TRUE a diagonal is printed between the origin and farest point, indicating relative expanded measurement uncertainty. A that the point the contribution of Bias and Relative Error is also plotted.


#' @param Verbose Optional, logical, default is FALSE. If TRUE messages are displayed during execution.
#'
#' @return A plot
#' @export
#'
#' @examples
#' data(pm25)
#' # Simple example for a PM2.5 sensor. The between reference standard uncertainty, ubsRM, is assumed to be constant along the PM2.5 range of concentrations. The relative expanded measurement uncertainty, Ur, is computed without adding any between sensor standard uncertainty, ubss. Ordinary Lest Square is used to compute the regression line. The default unit for sensor and reference data: 'ug/m3' is used.
#' mtd(pollutant = 'PM2.5',
#'     dates = pm25$datetime,
#'     reference_data = pm25$reference_pm25,
#'     sensor_data = pm25$sensor_pm25,
#'     ubsRM = 0.5)
#' # Example for a PM2.5 sensor with ubsRM being fitted along the PM2.5 range of concentrations. The relative expanded measurement uncertainty, Ur, is computed adding a constant between sensor standard uncertainty, ubss, of 2 µg/m³. Ordinary Lest Square is used to compute the regression line. The default unit for sensor and reference data: 'ug/m3' is used.
#' mtd(pollutant = 'PM2.5',
#'     dates = pm25$datetime,
#'     reference_data = pm25$reference_pm25,
#'     sensor_data = pm25$sensor_pm25, 
#'     Add.ubss = T, variable.ubss = T, ubss = rep(2, length(pm25$sensor_pm25)),
#'     variable.ubsRM = T, Fitted.RS = T, ubsRM = rep(0.2, length(pm25$datetime)), Verbose = T)

mtd <- function(pollutant, dates, sensor_data, Add.ubss = FALSE, variable.ubss = FALSE, ubss = NULL, unit_sensor = 'ug/m3', 
                reference_data, variable.ubsRM = FALSE, ubsRM = NULL, unit_ref = 'ug/m3', Fitted.RS = FALSE, Forced.Fitted.RS = FALSE, Regression = 'OLS',
                sensor_name = NULL, Max.percent = NULL, x_axis_label = NULL, y_axis_label = NULL,
                x_lim = NA, y_lim = NA, Show.Diag.Ur = FALSE, Verbose = FALSE)
{
    pollutant = match.arg(pollutant,
                          choices = c("CH4","CO","CO2","CO2", "H2O","NO","NO2","O3","PM1","PM10","PM2.5","RH","RH_int","Temp","Temp_int","Press"))
    
    unit_ref    = match.arg(unit_ref, choices = c('ug/m3', 'ppm', 'ppb', 'mg/m3', 'Percent'))
    unit_sensor = match.arg(unit_sensor, choices = c('ug/m3', 'ppm', 'ppb', 'mg/m3', 'Percent'))
    
    stopifnot(length(dates) == length(sensor_data))
    stopifnot(length(reference_data) == length(sensor_data))
    n_obs = length(reference_data)
    
    # Create a data.table for U_orth_DF
    Mat <- data.table::data.table(case   = 1:n_obs,
                                  date   = dates,
                                  xis    = reference_data,
                                  yis    = sensor_data)
    
    
    # Setting ubsRM and ubss in df
    if (!variable.ubsRM) {
        stopifnot(!is.null(ubsRM))
        data.table::set(Mat,  j = "ubsRM", value = rep(ubsRM, nrow(Mat))) 
    } else {
        stopifnot(length(ubsRM) == length(reference_data))
        data.table::set(Mat,  j = "ubsRM", value = ubsRM)}
    if (Add.ubss){
        if (!variable.ubss) {
            stopifnot(!is.null(ubss))
            data.table::set(Mat,  j = "ubss", value = rep(ubss, nrow(Mat)))
        } else {
            stopifnot(length(ubss) == length(sensor_data))
            data.table::set(Mat,  j = "ubss", value = ubss)}}
    
    # Change to data.table to be faster and check availability of data
    Mat <- Mat[is.finite(rowSums(Mat[, c("xis", "yis")]))]
    stopifnot(nrow(Mat) > 0)
    
    # Run orthogonal regression
    U.orth.List = U_orth_DF(Mat = Mat, Add.ubss = Add.ubss, variable.ubss = variable.ubss, variable.ubsRM = variable.ubsRM,
                            Fitted.RS = Fitted.RS, Forced.Fitted.RS = Forced.Fitted.RS, Regression = Regression, Verbose = Verbose)
    
    # Get the DQO thresholds
    DQO = get.DQO(name.gas = pollutant, unit.ref = unit_ref)
    
    # Run the target diagram function
    Target.Diagram(Sensor_name  = sensor_name,
                   Mat           = U.orth.List[["Mat"]],
                   variable.ubsRM = variable.ubsRM,
                   ubsRM         = ubsRM,
                   with.ubss     = Add.ubss,
                   ubss          = ubss,
                   b0            = U.orth.List$b0,
                   b1            = U.orth.List$b1,
                   Unit.Ref      = unit_ref,
                   Unit.sensor   = unit_sensor,
                   DQO.1         = DQO$DQO.1/DQO$LV,
                   DQO.2         = DQO$DQO.2/DQO$LV,
                   DQO.3         = DQO$DQO.3/DQO$LV,
                   LAT           = DQO$LAT,
                   UAT           = DQO$UAT,
                   LV            = DQO$LV,
                   AT            = DQO$AT,
                   sdm_sdo       = U.orth.List[["sdm"]] > U.orth.List[["sdo"]],
                   BeginEnd      = as.Date(range(dates, na.rm = T)),
                   xAxisLabel    = x_axis_label,
                   yAxisLabel    = y_axis_label,
                   Xlim          = NA,
                   Ylim          = NA,
                   Show.Diag.Ur  = Show.Diag.Ur)
}
