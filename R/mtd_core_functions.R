# Functions required for creating the Modified Target Diagram (MTD)

#================================================================CR
### U_orth_DF: Function Orthogonal regression without plotting (Vs 180505) ====
#================================================================CR
#' Function Orthogonal regression without plotting
#'
#' @param Mat : Data.table or DataFrame of data including Case number, Date, x, y + optional ubsRM and/or ubss if ubsRM and/or ubss are not constant for all xi. The columns shall be in the order: "case", "Date", "xis", "yis","ubsRM", "ubss" with whatever column names.
#' @param ubsRM : numeric (default = NULL ), random standard uncertainty of reference measurements, xis, given as a constant value for all xis reference values.
#' @param variable.ubsRM : logical, default is FALSE. If FALSE, ubsRM is used as constant random standard uncertainties for all xis reference values. If TRUE ubsRM given in Mat and is used for each reference value
#' @param perc.ubsRM : numeric default value 0.03. Use to compute ubsRm in case variable.ubsRM is TRUE as Mat$ubsRM = perc.ubsRM * Mat[["xis"]] 
#' @param ubss : numeric (default = NULL ), random standard uncertainty of sensor measurements, yis, given as a constant value for all yis sensor values
#' @param variable.ubss : logical, default is FALSE. If FALSE, ubss is used as constant random standard uncertainties for all yis sensor values. If TRUE ubss given in Mat and is used for each sensor value
#' @param perc.ubss : numeric default value 0.03. Use to compute ubss in case variable.ubss is TRUE as Mat$ubss = perc.ubss * Mat[["yis"]] 
#' @param Fitted.RS : logical, default is FALSE. If TRUE the square residuals (RS) can be fitted using a General Additive Model, according to the result of provided that the probability that the correlation between xis and RS is null is lower than 0.01, (p < 0.01)
#' @param Forced.Fitted.RS : logical, default is FALSE. If TRUE even if the variance of residuals is constant, RS is Gam fitted.
#' @param Regression character, default is "Orthogonal", possible values are "OLS" ,"Deming" and "Orthogonal". For Orthogonal Delta  is 1 and for "Deming" Delta is ubss^2/ubsRM^2. See https://en.wikipedia.org/wiki/Deming_regression
#' @param Add.ubss : logical, default is TRUE If TRUE ubss is added to  Mat$Rel.RSS. If FALSE ubss is not added to  Mat$Rel.RSS.
#' @param Verbose : logical, default is FALSE. If TRUE messages are displayed during execution.
#' @param Versus character, default is NUL. If not NULL name of the column in data.table Mat which is used with the gam fitting to fit RSi. If NULL, RSi will befitted versus reference data (xis). 

#' @return a list with parameters: "mo","sdo", "mm","sdm", "b1", "ub1", "b0", "ub0", "RSS","rmse", "mbe", "Correlation", "nb", "CRMSE", "NMSD" and a data.table called "Mat" with columns: "case", "Date", "xis", "yis","ubsRM", "RS", "Ur", "U", "Rel.bias", "Rel.RSS"
#' @details: Homogeneity of variance of residuals is tested For the calculation of RSi adding Ur In a new column of Mat
#' returning a list with slope (b and ub), intercept (a and ua), the sum of square of residuals (RSS),
#' the root means square of error (RMSE), the mean bias error (mbe), the coefficeint of correlation (Correlation),
#' the number of valid measurements (nb), the centered root mean square of error (CRMSE), the normalised mean standard deviation (NMSD) and Mat with relative expanded uncertainty.
#' Negative Mat$RS - Mat$ubsRM^2) are set to 0

#' @import data.table
#' @import car
#' @import futile.logger
#' @import lmtest
#' @import MethComp
#' @import plotrix
#' @import shape
#' @import skedastic
#'
#' @examples: empty
#'
#' @noRd
#'
U_orth_DF <- function(Mat, ubsRM = NULL, variable.ubsRM = FALSE, perc.ubsRM = 0.03, ubss = NULL, variable.ubss = FALSE, perc.ubss = 0.03, Fitted.RS = FALSE, Forced.Fitted.RS = FALSE, Regression = "Orthogonal",
                      Verbose = FALSE, Add.ubss = TRUE, Versus = NULL) {
    #checking that the Mat  is not empty
    if (exists("Mat") && !is.null(Mat) && nrow(Mat)>0) {
        
        # Setting Versus with xis is not NULL
        if (is.null(Versus)) Versus <- "xis"
        data.table::setkeyv(Mat, Versus)
        
        # Setting column names and adding ubsRM and ubss as constant values
        colnames(Mat) <- c("case", "Date", Versus, "yis","ubsRM", "ubss")[1:length(colnames(Mat))]
        
        # Convert Mat to data.table if needed
        if (!"data.table" %in% class(Mat)) Mat <- data.table(Mat, key = "Date")
        # Filtering for the validation data only
        Mat <- Mat[is.finite(rowSums(Mat[, c(Versus,"yis"), with = FALSE]))]
        if (!nrow(Mat) > 0) return(futile.logger::flog.error("[U_orth_DF] Mat does not contains any complete rows with xis and yis"))
        
        # Setting ubsRM and ubss in Mat
        if (!variable.ubsRM) {
            if (!is.null(ubsRM)) data.table::set(Mat,  j = "ubsRM", value = rep(ubsRM, nrow(Mat))) 
        } else if (!"ubsRM" %in% names(Mat)) {
            return(futile.logger::flog.error("[U_orth_DF] u(bs,RM) not given in data.table Mat. Stopping the function."))
        }  else Mat$ubsRM = perc.ubsRM * Mat[[Versus]]
        if (!variable.ubss) {
            if (!is.null(ubss)) data.table::set(Mat,  j = "ubss", value = rep(ubss, nrow(Mat)))
        } else if (!"ubss" %in% names(Mat)) {
            return(futile.logger::flog.error("[U_orth_DF] u(bs,s) not given in data.table Mat. Stopping the function."))
        } else Mat$ubss = perc.ubss * Mat[["yis"]]
        
        # Determination of delta for Deming or orthogonal regression
        if (Regression == "Deming") {
            if (variable.ubsRM || variable.ubss) {
                return(futile.logger::flog.error("[U_orth_DF] regression cannot be set to \"Deming\" with variable u(bs,RM) and/or variable u(bs,s)"))
            } else if (is.null(ubsRM) || is.null(ubsRM)) {
                return(futile.logger::flog.error("[U_orth_DF] with \"Deming\" regression both u(bs,RM) and u(bs,s) cannot be null"))
            } else {
                Delta <- (ubss/ubsRM)^2
                futile.logger::flog.info(paste0("[U_orth_DF] regression type: \"",Regression, "\" with Delta = ", Delta))}
        } else if (Regression == "Orthogonal") {
            Delta <- 1
            futile.logger::flog.info(paste0("[U_orth_DF] regression type: \"",Regression, "\" (Delta = ", Delta,")"))
        } else if (Regression != "OLS") return(futile.logger::flog.error("[U_orth_DF] unknown regression type. Only \"OLS\",  \"Deming\" (Delta is ubss^2/ubsRM^2) or \"orthogonal\" (Delta is 1) regressions can be used"))
        
        # Common parameters
        nb  <- nrow(Mat)
        mo  <- mean(Mat[[Versus]])
        mm  <- mean(Mat[["yis"]])
        sdo <- sd(Mat[[Versus]])
        sdm <- sd(Mat[["yis"]])
        Sxx <- sum((Mat[[Versus]] - mo)^2)
        Syy <- sum((Mat[["yis"]] - mm)^2)
        Sxy <- sum((Mat[[Versus]] - mo) * (Mat[["yis"]] - mm))
        # Linear Regression Orthogonal regression ()
        if (Regression == "Orthogonal") {
            # as in annex b of Guide for The Demonstration of Equivalence
            m2  <- MethComp::Deming(x = Mat[[Versus]], y = Mat[["yis"]], vr = Delta)
            b1  <- (Syy - Sxx + sqrt((Syy- Sxx)^2 + 4*Sxy^2))/(2*Sxy)
            b0  <- mm - b1 * mo
            ub1 <- sqrt((Syy - (Sxy^2/Sxx))/((nb-2)*Sxx))
            ub0 <- sqrt(ub1^2 * sum(Mat[[Versus]]^2)/nb)
        } else if (Regression == "Deming"){
            #tls::tls(yis ~ xis, data = Mat[, .(xis, yis)])
            m2  <- MethComp::Deming(x = Mat[[Versus]], y = Mat[["yis"]], vr = Delta)
            b0  <- m2[1] 
            b1  <- m2[2]
        } else if (Regression == "OLS") {
            # b1 = Sxy/Sxx
            # b0 = mm - b1 * mo
            # Syx = sqrt(sum((Mat[["yis"]] - (b0 + b1 * Mat[[Versus]]))^2)/(nrow(Mat)-2))
            # ub = Syx/sqrt(Sxx)
            # ua = Syx * sqrt(sum((Mat[[Versus]]^2)/(nrow(Mat) * Sxx)))
            Formula <- as.formula(paste0("yis ~ ", Versus))
            m2 <- lm(Formula, data = Mat)
            b0  <- coef(m2)[1] 
            b1  <- coef(m2)[2]
            st.errors <- sqrt(diag(vcov(m2)))
            ub0 <- st.errors[1]
            ub1 <- st.errors[2]}
        
        # Regression statistics for Target Diagram (see delta tool user guide)
        rmse  <- sqrt((sum((Mat[["yis"]] - (b0 + b1 * Mat[[Versus]]))^2))/(nb-2))
        mbe   <- mean(Mat[["yis"]] - Mat[[Versus]])
        mae   <- mean(abs(Mat[["yis"]] - Mat[[Versus]]))
        CRMSE <- sqrt(mean(((Mat[["yis"]] - mm) - (Mat[[Versus]] - mo))^2))
        NMSD  <- (sd(Mat[["yis"]]) - sd(Mat[[Versus]])) / sd(Mat[[Versus]])
        Correlation <- cor(Mat[[Versus]],Mat[["yis"]])
        # Squares of Residuals and bias (vector of values)
        Mat[, fitted := (b0 + b1 * Mat[[Versus]])]
        Mat[,   bias := (b0 + (b1 - 1) * Mat[[Versus]])] # Bias from identity line x
        Mat[, residuals := (Mat[["yis"]] - Mat[["fitted"]])]
        Mat[,     RS := (Mat[["yis"]] - Mat[["fitted"]])^2]
        # # creating an OLS model to apply Breusch Pagan test and to compute ub0 and ub1
        Formula <- as.formula(paste0("yis ~ ", Versus))
        OLS <- lm(Formula, data = Mat)
        OLS$coefficients  <- c(b0,b1)
        OLS$residuals     <- Mat$residuals
        OLS$fitted.values <- Mat$fitted
        if (Regression == "Deming") {
            # https://bookdown.org/ccolonescu/RPoE4/heteroskedasticity.html
            ub0 <- lmtest::coeftest(OLS, vcov. = car::hccm(OLS, type = "hc1"))[1, "Std. Error"]
            ub1 <- lmtest::coeftest(OLS, vcov. = car::hccm(OLS, type = "hc1"))[2, "Std. Error"]}
        
        # testing for heteroscedacity with Breusch Pagan test
        # https://www.r-bloggers.com/2016/01/how-to-detect-heteroscedasticity-and-rectify-it/
        Breusch.Pagan     <- lmtest::bptest(formula = OLS)
        skedastic::breusch_pagan(OLS)
        # testing significance of correlation between s and square of absolute residuals - The calculation does not work only possibility the constrant RSS
        rtest <- Breusch.Pagan
        
        # if fitting the residuals is needed, ordering before
        futile.logger::flog.info("[U_orth_DF] Breusch-Pagan: null hypothesis means constant variance of residuals along x axis (homoscedacity). The hyposthesis is rejected if p-value < 0.05 (heteroscedacity)")
        futile.logger::flog.info(paste0("[U_orth_DF] and the residuals would be heteroscedastic with 0.95 level of statistical confidence. Finally, p-value =  ", format(rtest$p.value, digits = 4)))
        if (!is.na(rtest$p.value) && Fitted.RS && (rtest$p.value < 0.05 || Forced.Fitted.RS)) {
            futile.logger::flog.info("[U_orth_DF] The variance of residuals is not constant. RSS is calculated after applying a General Additive Model fitting.")
            # Fitting with gam Vs Versus ("Xi")
            # if any y value is zero getting Warning: Error in eval: non-positive values not allowed for the 'gamma' family (we had 0.5 % of min(xis) to avoid this
            
            Formula <- as.formula(paste0("sqrt(RS) ~ s(", Versus, ")"))
            z <- mgcv::gam(Formula, data = Mat,family=Gamma(link=log) )
            
            # # see https://stats.stackexchange.com/questions/270124/how-to-choose-the-type-of-gam-parameters
            # z <- mgcv::gam(Formula, data = Mat, method = "REML", select = TRUE)
            
            Mat[, RS := fitted(z)^2]
            # Sum of squares of Residuals (one constant value)
            RSS     <- sum(Mat$RS)
            
            if (Verbose) print(summary(z))
        } else {
            futile.logger::flog.info("[U_orth_DF] The variance of residuals is constant. Constant RSS is calculated.")
            if (Fitted.RS) {
                futile.logger::flog.info(paste0("[U_orth_DF] Argument \"Fitted.RS\" in U_orth_DF(): ", Fitted.RS, ". If FALSE the square residuals are not fitted and constant RSS is computed."))
                futile.logger::flog.info(paste0("[U_orth_DF] if Breusch-Pagan lagrange multiplier statistic p-value > 0.05, residuals are homoscedastic with 0.95 level of statistical confidence. p -value = ",
                                                format(rtest$p.value, digits = 4)))}
            futile.logger::flog.info("[U_orth_DF] RSS is calculated with equation for constant residuals.")
            # Sum of squares of Residuals (one constant value)
            RSS     <- sum(Mat$RS)
            futile.logger::flog.info(paste0("[U_orth_DF] RSS is the square root of sum of squares of Residuals divided by n - 2 = ", format(sqrt(RSS/(nb-2)), digit = 3)))
            # No need to fit a line in this case
            Mat[, RS := rep(RSS/(nb-2), times = .N)]}
        
        # Plotting RS
        if (Verbose) {
            # See https://stackoverflow.com/questions/17093935/r-scatter-plot-symbol-color-represents-number-of-overlapping-points
            
            ## Use densCols() output to get density at each point
            x <- grDevices::densCols(Mat[[Versus]],sqrt(Mat$residuals^2), colramp=colorRampPalette(c("black", "white")))
            Mat$dens <- col2rgb(x)[1,] + 1L
            ## Map densities to colors
            cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                                                 "#FCFF00", "#FF9400", "#FF3100"))(256)
                                                 
            Mat$col <- cols[Mat$dens]
            
            plot(sqrt(Mat$residuals^2) ~ get(Versus), data = Mat, type = "p", col = col, xlab = Versus)
            lines(Mat[[Versus]],sqrt(Mat[["RS"]]), col = "red", xlab = Versus); grid()}
        
        # Checking if RSS^2 - Mat$ubsRM^2 < 0 that results in an error using sqrt(RSS^2 - Mat$ubsRM^2) of the rs.RSS. Replacing with 0
        neg.RSS <- which(Mat$RS - Mat[["ubsRM"]]^2 < 0)
        if (length(neg.RSS) > 0) {
            futile.logger::flog.warn("[U_orth_DF] Some \"RSS/(nb - 2) - ubsRM^2\" are negative and square roots cannot be calculated.")
            futile.logger::flog.info("[U_orth_DF] ubsRM maybe too high and could be modified.")
            futile.logger::flog.info("[U_orth_DF] The \"RSS/(nb - 2) - ubsRM^2\" that are negative will be set to 0 when computing uncertainties.")
        }  else {
            # mat$RS are not changed and they are already calculated
            futile.logger::flog.info("[U_orth_DF] All \"RSS/(nb - 2) - ubsRM^2\" are positives. ubsRM makes sence.")}
        #### Calculating parameters for modified Target diagram Rel.bias and Rel.RSS
        Mat[, Rel.bias := 2 * (b0/Mat[[Versus]] + (b1 - 1))]
        # (Mat$RS - Mat$ubsRM^2) cannot be negative, adding or not ubss
        if (Add.ubss){
            if (length(neg.RSS) == 0) {
                Mat[, Rel.RSS := 2 * sqrt(Mat$ubss^2 + Mat$RS - Mat$ubsRM^2) / Mat[[Versus]]]
            } else {
                Positives <- setdiff(1:nrow(Mat),neg.RSS)
                if (length(Positives) > 0) Mat[Positives, Rel.RSS := 2 * sqrt(Mat[Positives,ubss]^2 + Mat[Positives,RS] - Mat[Positives,ubsRM]^2) / Mat[Positives][[Versus]]]
                Mat[neg.RSS,  Rel.RSS := 0]}
        } else {
            if (length(neg.RSS) == 0) {
                Mat[, Rel.RSS := 2 * sqrt(Mat$RS - Mat$ubsRM^2) / Mat[[Versus]]]
            } else {
                Positives <- setdiff(1:nrow(Mat),neg.RSS)
                if (length(Positives) > 0) Mat[Positives, Rel.RSS := 2 * sqrt(Mat[Positives,RS] - Mat[Positives,ubsRM]^2) / Mat[Positives][[Versus]]]
                Mat[neg.RSS,  Rel.RSS := 0]}} 
        #### Calculating uncertainty
        Mat[, Ur := sqrt(Mat$Rel.bias^2 + Mat$Rel.RSS^2) * 100]
        Mat[, U  := Mat$Ur / 100 * Mat[[Versus]]]
        # Indicators
        Mat[, Max.ubsRM := sqrt(RSS/(nb-2) + Mat$bias^2)]
        Mat[, Max.RSD   := sqrt(RSS/(nb-2) + Mat$bias^2)/Mat[[Versus]]]
        # Printing
        if (Verbose) {
            cat("--------------------------------\n")
            cat(sprintf("mean of x   : %.1g +/- %.1g",mo,sdo),"\n")
            cat(sprintf("Intercept b0: %.4g +/- %.4g",mm,sdm), "\n")
            cat(sprintf("u of b1     : %.4g +/- %.4g",b1,ub1),"\n")
            cat(sprintf("u of b0     : %.4g +/- %.4g",b0,ub0), "\n")
            cat(sprintf("R2: %.4g",Correlation^2), "\n")
            if (!is.na(rtest$p.value) && Fitted.RS && (rtest$p.value < 0.05 || Forced.Fitted.RS)) {
                cat("The residuals are not constant. RS are fitted with a general Additive model (k=5) see in returned matrix. \n")
            } else {
                cat("The residuals are constant. RSS is calculated with equation for constant residuals:")
                cat(sprintf("RSS: %.4g ",Mat$RSS[1]), "\n")}
            cat(sprintf("RMSE : %.4g ",rmse), "\n")
            cat(sprintf("mbe  : %.4g ",mbe), "\n")
            cat(sprintf("CRMSE: %.4g ",CRMSE), "\n")
            cat(sprintf("NMSD : %.4g ",NMSD), "\n")
            cat(sprintf("n    : %.4g ",nb), "\n")}
        calib <- list(mo,sdo, mm,sdm, b1, ub1, b0, ub0, RSS, rmse, mbe, Correlation, nb, CRMSE, NMSD, Mat, Fitted.RS, Regression)
    } else {
        cat("Mat is empty. Returning NAs.")
        calib <- list(mo = NA,sdo = NA, mm = NA,sdm = NA, b1 = NA, ub1 = NA, b0 = NA, ub0 = NA, RSS = NA,rmse = NA, mbe = NA, Correlation = NA, nb = NA, CRMSE = NA, NMSD = NA, Mat = NA, Regression = Regression)}
    names(calib) <- c("mo","sdo", "mm","sdm", "b1", "ub1", "b0", "ub0", "RSS","rmse", "mbe", "Correlation", "nb", "CRMSE", "NMSD", "Mat", "RS.Fitted", "Regression")
    
    # if (Verbose) print(data.frame(x              = Mat[[Versus]],
    #                               ubsRM          = Mat$ubsRM,
    #                               Max.ubsRM      = sqrt(RSS/(nb-2) + Mat$bias^2),
    #                               Max.RSD        = sqrt(RSS/(nb-2) + Mat$bias^2)/Mat[[Versus]],
    #                               Decrease.ubsRM = (Mat$ubsRM - sqrt(RSS/(nb-2) + Mat$bias^2)) > 0))
    
    return(calib)
}


#' @description  Get Critical Value, Information Threshold, Alert Threshold, Limit Value, Upper Assessment Threshold, Lower Assessment Threshold and Data Quality Objectives from name of the molecule, the molecule or the sensor name according to what is set into the European Air Quality Directive (2008). ---
#' @param name.gas,name.sensor,gas.sensor character strings, default is NULL. One of them shall not be NULL, Molecule or pollutant symbol (CO, CO2, NO, NO2, O3, PM1, PM1, PM10, PM2.5),
#'  brand name of sensor (CO_A4_P1, D300, NO_B4_P1, NO2_B43F_P1, OX_A431_P1, 5301CAT, 5301CST, OPCN3PM1, 5310CAT, 5310CST, OPCN3PM10, 5325CAT, 5325CST, OPCN3PM25) and
#'   molecule name("Carbon_monoxide", "Nitrogen_dioxide", "Nitric_oxide", "Ozone", "Sulphur_dioxide", "Benzene", "Particulate_Matter_10", "PM10_PMSCal", "Particulate_Matter_25", "PM25_PMSCal", "Carbon_dioxide")
#' @param Averaging.Period  character, Averaging periods as defined in the European Air Quality Directive. 
#' The Averaging Period changes the LV of pollutants. Values can be: "1hour", "8hour" and "1year". Default is "1hour"
#' @param unit.ref  character, units in which parameters are returned. It can be: "ug/m3", "mg/m3", "ppb" or "ppm". Default is "ug/m3".
#' @param Candidate  character, Data Quality Objectives as stated in the European Air Quality. It can be "Sensor" or "Model". Default is "Sensor"
#' According to Candidate the percentage of DQO changes. For "Sensor", the DQO of the CENT TC264/WG42 are used. For "Model" the DQO of the Air Quality Directive are used.
#' @return a list with , CL, IT, AT, LV, LAT, UAT, DQO.1, DQO.2, DQO.3 with NA for undefined parameters. They are returned in unit given by unit.ref.
#'
#' @import data.table
#' @import car
#' @import futile.logger
#' @import lmtest
#' @import MethComp
#' @import plotrix
#' @import shape
#' @import skedastic
#'
#' @examples NA
#'
#' @noRd
#'
get.DQO <- function(gas.sensor = NULL, name.sensor = NULL, name.gas = NULL, Averaging.Period = "1hour", unit.ref = "ppb", Candidate = "Sensor") {
    # Determining gas.sensor
    DT.gas <- data.table::data.table(name.gas    = c(   "H2O"          ,        "CH4","CO"             ,            "CO2",            "CO2",           "NO",              "NO2",         "O3",        "PM1",        "PM1",                  "PM1",        "PM10",        "PM10",                  "PM10",       "PM2.5",       "PM2.5",                 "PM2.5",                "RH",                "RH_int",        "Temp",        "Temp_int",                "Press"),
                                     name.sensor = c("MPLCPC"          ,     "LPLCPC","CO_A4_P1"       ,         "SPLCPC",           "D300",     "NO_B4_P1",      "NO2_B43F_P1", "OX_A431_P1",    "5301CAT",    "5301CST",             "OPCN3PM1",     "5310CAT",     "5310CST",             "OPCN3PM10",     "5325CAT",     "5325CST",             "OPCN3PM25",           "SHT31HE",               "SHT31HI",     "SHT31TE",         "SHT31TI",               "BMP280"),
                                     gas.sensor  = c("K96_Water_vapour","K96_Methane","Carbon_monoxide", "Carbon_dioxide", "Carbon_dioxide", "Nitric_oxide", "Nitrogen_dioxide",      "Ozone", "PM1_PMSCal", "PM1_PMSraw", "Particulate_Matter_1", "PM10_PMSCal", "PM10_PMSraw", "Particulate_Matter_10", "PM25_PMSCal", "PM25_PMSraw", "Particulate_Matter_25", "Relative_humidity", "Relative_humidity_int", "Temperature", "Temperature_int", "Atmospheric_pressure"))
    
    if (is.null(name.gas) && is.null(name.sensor) && is.null(gas.sensor)){
        stop(futile.logger::flog.error("[get.DQO] either name.gas, name.sensor, gas.sensor shall be given"))
    } else if (!is.null(name.sensor)){
        row.DT     <- which(DT.gas$name.sensor == name.sensor)
        if (length(row.DT) > 0){
            name.gas   <- DT.gas[row.DT, name.gas][1]
            gas.sensor <- DT.gas[row.DT, gas.sensor][1]
        } else stop(futile.logger::flog.error("[get.DQO] unknown name.sensor"))
    } else if (!is.null(gas.sensor)){
        row.DT      <- which(DT.gas$gas.sensor == gas.sensor)
        if (length(row.DT) > 0){
            name.gas    <- DT.gas[row.DT, name.gas][1]
            name.sensor <- DT.gas[row.DT, name.sensor][1]
        } else stop(futile.logger::flog.error("[get.DQO] unknown gas.sensor"))
    } else if (!is.null(name.gas)){
        row.DT      <- which(DT.gas$name.gas == name.gas)
        if (length(row.DT) > 0){
            gas.sensor <- DT.gas[row.DT, gas.sensor][1]
            name.sensor <- DT.gas[row.DT, name.sensor][1]
        } else stop(futile.logger::flog.error("[get.DQO] unknown name.gas"))
    }
    # Checking if gas sensors or name.gas is given
    # Checking consistency of arguments
    if (!gas.sensor %in% c("K96_Water_vapour","K96_Methane","Carbon_monoxide", "Carbon_dioxide", "Nitrogen_dioxide", "Nitric_oxide", "Ozone", "Sulphur_dioxide", "Benzene",
                           "Particulate_Matter_10", "PM10_PMSCal", "PM10_PMSraw", "Particulate_Matter_25", "Particulate_Matter_1", "PM25_PMSCal", "PM25_PMSraw", 
                           "PM1_PMSCal", "PM1_PMSraw","Carbon_dioxide", "Relative_humidity", "Relative_humidity_int", "Temperature", "Temperature_int", "Atmospheric_pressure")) {
        return(futile.logger::flog.error("[get.DQO] unknown compound, no DQO"))}
    if (!Averaging.Period %in% c("1hour", "8hour", "24hour","1year")) return(futile.logger::flog.error("[get.DQO] unknown Averaging.Period"))
    if (!unit.ref %in% c("ug/m3", "mg/m3", "ppb", "ppm", "percent", "Celsius", "hPa")) return(futile.logger::flog.error("[get.DQO] unknown unit.ref"))
    if (!Candidate %in% c("Sensor", "Model")) return(futile.logger::flog.error("[get.DQO] unknown Candidate"))
    # Limit Value, DQOs, UAT, LAT
    IT = NA
    AT = NA
    if (name.gas == "H2O") {
        LV = 2000
        DQO.1 = 0.25 * LV
        DQO.2 = 0.75 * LV
        DQO.3 = 2.00 * LV
        LAT   = 0.50 * LV
        UAT   = 0.70 * LV
    } else if (name.gas == "CH4") {
        LV = 2000
        DQO.1 = 0.25 * LV
        DQO.2 = 0.75 * LV
        DQO.3 = 2.00 * LV
        LAT   = 0.50 * LV
        UAT   = 0.70 * LV
    } else if (name.gas == "CO") {
        if (unit.ref == "ppm") {
            LV = 10/1.34
        } else if (unit.ref == "ppb") {
            LV = 10000/1.34
        } else if (unit.ref == "mg/m3") {
            LV = 10
        } else if (unit.ref == "ug/m3") {
            LV = 10000}
        DQO.1 = 0.25 * LV
        DQO.2 = 0.75 * LV
        DQO.3 = 2.00 * LV
        LAT   = 0.50 * LV
        UAT   = 0.70 * LV
    } else if (name.gas == "NO2") {
        if (unit.ref == "ppb") {
            LV = 200/1.91
            AT = 400/1.91
        } else if (unit.ref == "ug/m3") {
            LV = 200
            AT = 400}
        DQO.1 = 0.25 * LV
        DQO.2 = 0.75 * LV
        DQO.3 = 2.00 * LV
        LAT   = 0.50 * LV
        UAT   = 0.70 * LV
    } else if (name.gas == "NO") {
        if (unit.ref == "ppb") {
            LV = 200/1.25
            AT = 400/1.25
        } else if (unit.ref == "ug/m3") {
            LV = 200
            AT = 400
        } else cat(paste0("Wrong unit for ",gas.sensor, "\n"))
        DQO.1 = 0.25 * LV
        DQO.2 = 0.75 * LV
        DQO.3 = 2.00 * LV
        LAT   = 0.50 * LV
        UAT   = 0.70 * LV
    } else if (name.gas == "O3") {
        if (unit.ref == "ppb") {
            LV = 120/2.05
            IT = 180/2.05
            AT = 240/2.05
        } else if (unit.ref == "ug/m3") {
            LV = 120
            IT = 180
            AT = 240
        } else cat(paste0("Wrong unit for ",gas.sensor, "\n"))
        DQO.1 = 0.30 * LV
        DQO.2 = 0.75 * LV
        DQO.3 = 2.00 * LV
        LAT   = NA
        UAT   = NA
    } else if (name.gas == "SO2") {
        # Using LV for 1 year time average
        if (unit.ref == "ppb") {
            LV = 120/2.05
            AT = 500/2.05
        } else if (unit.ref == "ug/m3") {LV = 350; IT = NA; AT = 500}
        DQO.1 = 0.25 * LV
        DQO.2 = 0.70 * LV
        DQO.3 = 2.00 * LV
        LAT   = 0.40 * LV
        UAT   = 0.60 * LV
    } else if (name.gas == "Benzene") {
        # Using LV for 1 year time average
        if (unit.ref == "ppb") {
            LV = 5/2.05
        } else if (unit.ref == "ug/m3") {LV = 5; IT = NA; AT = NA}
        DQO.1 = 0.30 * LV
        DQO.2 = 1.00 * LV
        DQO.3 = 2.00 * LV
        LAT   = 0.40 * LV
        UAT   = 0.70 * LV
    } else if (name.gas == "PM10") {
        # Using LV for 24 hours time average
        LV    = 50
        DQO.1 = 0.50 * LV
        DQO.2 = 1.00 * LV
        DQO.3 = 2.00 * LV
        LAT   = 0.50 * LV
        UAT   = 0.70 * LV
    }  else if (name.gas == "PM2.5") {
        # Using LV for 24 hours time average
        LV    = 25
        DQO.1 = 0.50 * LV
        DQO.2 = 1.00 * LV
        DQO.3 = 2.00 * LV
        LAT   = 0.50 * LV
        UAT   = 0.70 * LV
    }  else if (name.gas == "PM1") {
        # Fake DQO only for plotting
        LV    = 20
        IT    = NA
        AT    = NA
        DQO.1 = 0.50 * LV
        DQO.2 = 1.00 * LV
        DQO.3 = 2.00 * LV
        LAT   = 0.50 * LV
        UAT   = 0.70 * LV
    } else if (name.gas == "CO2") {
        # Using LV for 24 hours time average
        LV    = 500
        DQO.1 = 0.30 * LV
        DQO.2 = 0.50 * LV
        DQO.3 = 1.00 * LV
        LAT   = 0.50 * LV
        UAT   = 0.70 * LV
    } else if (name.gas %in% c("RH", "RH_int")) {
        # Using LV for 24 hours time average
        LV    = 50
        DQO.1 = 0.30 * LV
        DQO.2 = 0.50 * LV
        DQO.3 = 1.00 * LV
        LAT   = NA
        UAT   = NA
    } else if (name.gas %in% c("Temp","Temp_int")) {
        # Using LV for 24 hours time average
        LV    = 15
        DQO.1 = 0.30 * LV
        DQO.2 = 0.50 * LV
        DQO.3 = 1.00 * LV
        LAT   = NA
        UAT   = NA
    } else if (name.gas == "Press") {
        # Using LV for 24 hours time average
        LV    = 1013
        DQO.1 = 0.30 * LV
        DQO.2 = 0.50 * LV
        DQO.3 = 1.00 * LV
        LAT   = NA
        UAT   = NA
    } else {
        LV    = NA
        IT    = NA
        AT    = NA
        DQO.1 = 0.05 * LV
        DQO.2 = 0.20 * LV
        DQO.3 = 0.30 * LV
        LAT   = 0.50 * LV
        UAT   = 0.70 * LV}
    return(list(IT = IT, AT = AT, LV = LV, UAT = UAT, LAT = LAT, DQO.1 = DQO.1, DQO.2 = DQO.2, DQO.3 = DQO.3))
}



#' Plot a Modified Target Diagram ---
#'
#' @description
#' Base on script by Aaron Albin
#' Use this script when the following is true:
#' Inside Mat there are four sets of numeric data stored as three separate columns (xis, yis, Rel.bias, Rel.RSS).
#' Slope and intercept of an OLS, Deming or orthognal regression fitted to Mat$yis vs Mat$yis are given.
#' xis, yis, Rel.bias, Rel.RSS, Slope and intercept are computed with function U_orth_DF().
#' You wish to plot Rel.bias, Rel.RSS onto a 'scatterplot', where the horizontal axis corresponds to Rel.RSS and the vertical axis corresponds to Rel.bias.
#' Each individual pair of numbers Rel.bias and Rel.RSS, then, is represented as a point in this two-dimensional space.
#' @param Sensor_name name of the sensor to be written in front of the calibration equation. If NULL, do not print sensor name.
#' @param Mat  Datatable or dataframe of data including Case number, Date, xis, yis, [ ubss and ubsRM if not constant], Rel.bias, Rel.RSS. Rel.bias, Rel.RSS, xis must be included into dataFrame Mat. It is easier to get it from function U.orth.List()
#' @param ubsRM numeric (default = NULL ). Random standard uncertainty of the results of the reference method, xis, given as a constant value for all xis reference values
#' @param ubss numeric (default = NULL ): Between standard uncertainty of the sensor, yis, given as a constant value for all yis sensor values
#' @param b0, numeric, intercept of the orthogonal regression, default: NULL. If not NULL b0/xi is plotted
#' @param b1 numeric, slope of the orthogonal regression, default is NULL. If not NULL b1 - 1 is plotted
#' @param Unit.Ref Character, default is NULL. Unit of reference values, can "ppm", "ppb", "ug/m3", "mg/m3" ...
#' @param Unit.sensor character vector, unit for the expanded uncertainty of yis and yis. Default is Unit.Ref.
#' @param Xlabel,Ylabel label of the x and axis
#' @param Xlim,Ylim limits of x and y axis, default values is NA, vectors of two values min and max values. Xlim and Ylim overruled Max.Percent
#' @param Max.percent: numeric in percent, maximum extent of the x and y axis of the Target Diagram. This is reespected provided that Ur smaller than Max.percent exist.
#' @param MainTitle character, title to appear On the top of the scatter plot of x And y values - NOT USED ANYMORE
#' @param DQO.1,DQO.2,DQO.3 numeric, data quality objectives for Indicative measurements, Modelling and objective estimation. The DQOs are expressed in percentage. Defaul is NA, if NA no DQO target circle is plotted. Use function get DQO()
#' @param LAT,UAT,LV,AT,CL numeric, lower and upper assessment threshold, limit value, Alert threshold and Critical level of the European Air Quality Directive for Mat$xis, same unit as Mat$xis, default value = NA, used for color scale and target circles
#' @param Disk,WD,Dir Character vectors where you put the graphic files (Disk, working directory, directory), It is sufficient if only Dir is supplied
#' @param variable.ubsRM logical (default = FALSE ), if FALSE ubsRM is used as constant random standard uncertainties for all xis reference values. If TRUE ubsRM given in Mat and is used for each reference values
#' @param f_coef1,f_coef2,f_R2 numeric, number of digit for intercept, slope and R2 using sprintf syntax. f_coef1 is used both for intercept and s(Res), f_coef2 is used for all parameters of the model apart from the intercept
#' @param nameModel character, name of model to be used to save uncertainty plots, character, default NULL
#' @param sdm_sdo logical, default is NULL. It shall be TRUE if the standard deviation of yis is lower than the one of xi and conversely.
#' If Null, sdm_sdo is determined using sd(Mat[Index.Good, yis]) and sd(Mat[Index.Good, yis])
#' @param SavePlot logical, default is TRUE if TRUE uncertainty polts are saved
#' @param Model.used character, default is NULL. Name of calibration model used to compute yis. Only used if MainTitle is null.
#' @param BeginEnd character vector representing the starting and ending date, default is null. Only used if MainTitle is null.
#' @param Time.average: numeric, default is null. Tme average in minutes. . Only used if MainTitle is null.
#' @param with.ubss: logical, default is TRUE. If FALSE, ubss is subtracted to Rel.RSS. In this case Rel.RSS = 2 * (sqrt((Mat$RS - ubsRM^2) / Mat$xis).
#'                   If TRUE Rel.RSS is 2 * (sqrt(ubss^2 + Mat$RS - ubsRM^2)/ Mat$xis)
#' @param Ref_Analysers name of reference analyser, default is NULL. If not null the name is added in Target_diagram
#' @return Plot a Target diagram that is saved as pdf file provided that SavePlot is TRUE and unless error message is returned.
#'
#' @import data.table
#' @import car
#' @import futile.logger
#' @import lmtest
#' @import MethComp
#' @import plotrix
#' @import shape
#' @import skedastic
#'
#' @examples TBD
#'
#' @noRd
#'
Target.Diagram <- function(Sensor_name, Mat, ubsRM = NULL, ubss = NULL, b0 = NULL, b1 = NULL,  Unit.Ref = NULL, Unit.sensor = Unit.Ref,
                           xAxisLabel = NULL, yAxisLabel = NULL, Xlim = NA, Ylim = NA, MainTitle = NULL,
                           DQO.1 = NA, DQO.2 = NA, DQO.3 = NA, CL = NA, IT = NA, AT = NA, LV = NA, LAT = NA, UAT = NA,
                           Disk = NA, WD = NA, Dir = NA, sdm_sdo = NULL, SavePlot = TRUE,
                           Model.used = NULL, BeginEnd = NULL, Time.average = NULL, Max.percent = NULL, with.ubss = TRUE, Ref_Analysers = NULL,
                           Regression = "OLS", Fitted.RS  = FALSE, variable.ubsRM = FALSE) {
    #=====[Consistency checks]=====
    # Ordering Mat according to xis to be able to create a color pallete e
    Mat <- Mat[order(Mat$xis),]
    #checking that Mat is dataFrame
    if (!is.data.frame(Mat) || !data.table::is.data.table(Mat)) return(futile.logger::flog.error("[Target.Diagram] Mat is not a data.table or dataFrame."))
    # convert to data.table if needed
    if (!data.table::is.data.table(Mat)) Mat <- data.table::data.table(Mat, key = "date")
    #checking that the mat dataFrame is not empty, considering only the complete cases with Xis and yis > 0
    Mat <- Mat[complete.cases(Mat) & xis > 0 & yis > 0,]
    if (!nrow(Mat) > 0) return(futile.logger::flog.error("[Target.Diagram] Mat is empty or not complete. Returning NAs.\n"))
    # checking if Mat includes "Rel.bias", "Rel.RSS", "xis"
    if (!all(c("Rel.bias", "Rel.RSS", "xis") %in% colnames(Mat)) || any(is.null(c("b0","b1")))) {
        return(futile.logger::flog.error("[Target.Diagram] Mat does not contain Rel.bias, Rel.RSS or xis or b0 and/or b1 are null.\n"))}
    # plotting the Target Diagram
    # Create X-Y scatterplot: Coded by color by Aaron Albin
    #=====[Set XData, YData, b0 and b1 variables]=====
    # Rel.RSS is on the x axis
    # Rel.Bias is on the y axis
    # Values are multiplied by 100 to plot in percentage
    if (!with.ubss) {
        Sign.RS <- Mat$RS - Mat$ubsRM^2
        Negative.RS <- which(Sign.RS < 0)
        if (length(Negative.RS) > 0) {
            set(Mat, i = Negative.RS, j = "Rel.RSS", value = rep(0, times = length(Negative.RS)))
            Positive.RS <- setdiff(seq(Sign.RS), Negative.RS)
            if (length(Positive.RS) > 0) set(Mat, i = Positive.RS, j = "Rel.RSS", value = 2 * sqrt(Mat$RS[Positive.RS] - Mat$ubsRM[Positive.RS]^2) / Mat$xis[Positive.RS])
        } else set(Mat, j = "Rel.RSS", value = 2 * sqrt(Mat$RS - Mat$ubsRM^2) / Mat$xis)}
    Mat[, XData :=  Rel.RSS * 100]
    Mat[, YData :=  Rel.bias * 100]
    # re-computing Ur in case RSS was changed
    set(Mat, j = "Ur", value = sqrt(Mat$XData^2 + Mat$YData^2))
    # For the decomposition of Rel.Bias
    # Contribution of b0 to Rel.Bias
    Mat[, b0 := 2 * b0 / xis * 100]
    # Contribution of b1 to Rel.Bias
    Mat[, b1 := 2 * rep((b1-1) * 100, times = nrow(Mat))]
    #=====[Set Index.Good, Xlim and Ylim]=====
    # To adjust xlim and ylim, set the arguments Xlim and Ylim when calling Target.Diagram()
    # Set Xlim or Ylim equal to 'c( ___, ___ )', where the two '___'s indicate the lower and upper bound (respectively) of the range for the relevant axis.
    # If you want the axis to be backwards/flip-flopped, add the upper bound first, followed by the lower bound, e.g. rather than c(5,100) say c(100,5).
    # When leaving the values Xlim and Ylim as NULL, R determines the limits (i.e. range) of the x (horizontal) axis and y (vertical) axis for you.
    if (is.null(Max.percent)) Max.percent <- DQO.3 * 100
    if (all(is.na(Xlim))) {
        if (any(abs(Mat$XData) <= Max.percent)) {
            # There are XData < Max.percent
            Index.Good.X <- which(abs(Mat$XData) <= Max.percent)
        } else {
            # All XData > Max.percent, creating a WRONG INDEX
            Index.Good.X <- NULL}}
    if (all(is.na(Ylim))) {
        if (any(abs(Mat$YData) < Max.percent)) {
            # There are YData < Max.percent
            Index.Good.Y <- which(abs(Mat$YData) <= Max.percent)
        } else {
            # All YData > Max.percent, creating a WRONG INDEX nothing to plot
            Index.Good.Y <- NULL}}
    # index of intersection of Index.Good.X and Index.Good.Y
    if (!is.null(Index.Good.X)) {
        if (!is.null(Index.Good.Y)) {
            Index.Good <- intersect(x = Index.Good.X, y = Index.Good.Y)
        } else Index.Good <- NULL
    } else Index.Good <- NULL
    if (all(is.na(Xlim))) {
        if (!is.null(Index.Good) && length(Index.Good) > 0) {
            Xlim <- c(max(-Max.percent,min(Mat$b1[Index.Good],Mat$b0[Index.Good],Mat$XData[Index.Good])), min(Max.percent, max(Mat$XData[Index.Good])))
        } else Xlim <- c(0, min(Max.percent, DQO.3 * 100))}
    if (all(is.na(Ylim))) {
        if (!is.null(Index.Good) && length(Index.Good) > 0) {
            Ylim <- c(max(-Max.percent,min(0, min(Mat$YData[Index.Good]))), min(Max.percent, max(DQO.3, max(Mat$YData[Index.Good]))))
        } else Ylim <- c(0, min(Max.percent, DQO.3 * 100))}
    if (is.null(Index.Good))     return(futile.logger::flog.warn(paste0("[Target.Diagram] ", Sensor_name, " : no data with Ur < ", Max.percent, "%")))
    if (!length(Index.Good) > 0) return(futile.logger::flog.warn(paste0("[Target.Diagram] ", Sensor_name, " : no data with bias or RSS < ", Max.percent, "%")))
    # Flipping Xlim if standard deviation of yis is lower than the one of yis using sdm_sdo
    if (is.null(sdm_sdo)) sdm_sdo <- sd(Mat[Index.Good, yis]) <= sd(Mat[Index.Good, xis])
    if (sdm_sdo) Xlim <- rev(Xlim)
    #=====[MainTitle]=====
    # Set the name for the main title (centered along the top of the plot) here.
    if (is.null(MainTitle)) {
        if (!is.null(Model.used)) {
            MainTitle = gsub("__", "_", sub('.rdata', '', basename(Model.used)))
        } else {
            if (!is.null(Sensor_name)) MainTitle = paste0(Sensor_name, " - Target Diagram - Relative expanded uncertainty") else MainTitle = "Target Diagram - Relative expanded uncertainty"}}
    #=====[axis labels]=====
    # In the argument 'xAxisLabel', type the words(s) you want to see displayed beneath the x (horizontal) axis.
    # In the argument 'yAxisLabel', type the word(s) you want to see displayed to the left of the y (vertical) axis.
    # If you don't want to include an axis label at all for either of these, just leave the double-quotes empty, i.e. "", for that one.
    if (is.null(xAxisLabel)) {
        if (with.ubss) {
            xAxisLabel = bquote(paste("Relative random effect in bold: RR = 2 ", sqrt("u"[bs_s]^2*" + RMSE"^2*" - u"[bs_RM]^2), "/X"[i]*"; oblique line: 2b"[0]*"/X"[i]*"; vertical line: 2(b"[1]*"- 1) in %"))
        } else {
            xAxisLabel = bquote(paste("Relative random effect in bold: RR = 2 ", sqrt("RMSE"^2*" - u"[bs_RM]^2), "/X"[i]*"; oblique line: 2b"[0]*"/X"[i]*"; vertical line: 2(b"[1]*"- 1) in %"))}}
    if (is.null(yAxisLabel)) yAxisLabel = bquote("Relative bias: RB = 2 (b"[0]*"/X"[i]*" + (b"[1]*" - 1)) in %")
    #=====[LogAxis]=====
    # In some applications, the interesting patterns in your numbers are happening at small numbers (e.g. under 10), but you have a few outliers (e.g. over 100) that are forcing R to zoom out far enough to encompass all the data, so far that you can't see the interesting patterns in your data in the small-number range.
    # To overcome this problem, you can have R plot one or both of the axes 'logarithmically'.
    # What this means is that, visually, smaller numbers will be 'given more space', whereas the higher numbers will be 'scrunched together' more.
    # Below, indicate which axes
    # PlotWhichAxesLogarithmically = "" means you don't want either axis to be plotted logarithmically
    # PlotWhichAxesLogarithmically = "x" means you want the horizontal axis to be plotted logarithmically
    # PlotWhichAxesLogarithmically = "y" means you want the vertical axis to be plotted logarithmically
    # PlotWhichAxesLogarithmically = "xy" means you want *both* axes to be plotted logarithmically
    # Note: If you choose to plot an axis logarithmically, any zeroes and any negative numbers will be omitted from the plot. (R will give you a warning telling you how many total points this was.) This is necessary because of how the underlying logarithmic transformation itself works.
    PlotWhichAxesLogarithmically = ""
    #=====[PointSymbol]=====
    # Choose here the symbol that you want to use for all the points that get plotted in the scatter plot. Each symbol is represented by a special 'code number':
    #  1 = unfilled circle
    #  0 = unfilled square
    #  2 = unfilled triangle
    #  5 = unfilled diamond
    # 16 =   filled circle
    # 15 =   filled square
    # 17 =   filled triangle
    # 18 =   filled diamond
    #  3 = a plus sign (i.e. '+')
    #  4 = an 'x' symbol
    #  8 = a star composed of '+' and 'x' superimposed
    # You should *not* surround any numbers like these with quotes.
    # Indicate which one you like by typing the code number below after the equals sign (without quotes).
    PointSymbol = 16
    #=====[LegendLabels]=====
    # The labels plotted are the ones in factor.Color that does not include the legislative levels CL, IT, AT, LV, LAT, UAT
    # Levels include Which Targets are not NULL
    Levels       <- c(CL, IT, AT, LV, LAT, UAT)[which(!is.na(c(CL, IT, AT, LV, LAT, UAT)))]
    Levels       <- Levels[Levels >= min(Mat[Index.Good,xis]) & Levels <= max(Mat[Index.Good,xis])]
    factor.Color <- pretty(Mat[Index.Good,xis], n = (10 - length(levels)))
    factor.Color <- round(seq(from = min(Mat[Index.Good,xis]), to = max(Mat[Index.Good,xis]), length.out = 10), digits = 1)
    LegendLabels <- round(unique(sort(c(factor.Color, Levels))), digits = 0)
    #=====[colorRampPalette]=====
    # Here, you can choose the colors that you would like R to use when drawing all of your points on the scatterplot.
    # Be sure to surround the name of the colors with double-quotes, e.g. "black".
    # Type 'colors()' (without the single-quotes) at the R console and hit ENTER to see the list of all the color names R accepts.
    # Hint: Colors like red, green, blue, cyan, magenta, yellow, and many others can be easily darkened by putting a number after them.
    #       Taking green for instance, "green" is the same thing as "green1".
    #       To darken the color, you can make it "green2", "green3" or "green4".
    #       This stops at 4, however - there is no "green5" or anything beyond that.
    #       Grey colors are different. There is "grey0" (black) to "grey100" (white) and everything in-between.
    #          This gives you lots of flexibility in specifying greyscale colors. (Note: The spelling 'gray' also works.)
    # Below, inside the 'combine' function 'c()', list which color you would like associated with each character by giving the color name for each in quotes.
    # The ordering will match the order of labels you just saw in the last step - in other words, the first thing there will be the first color here, the second will be the second, and so on.
    # It's perfectly OK to have more colors specified here than you will actually use. Think of this list here as the full default hierarchy that is referred to in making each plot.
    #         [1]     [2]     [3]       [4]       [5]        [6]          [7]       [8]       [9]        [10]        [11]   [12]      [13]
    coul <- c( "black", "red", "green3", "brown4", "magenta", "magenta4" , "purple", "purple4", "orange", "orange4", "gray", "gray33", "gray66", topo.colors(12),rainbow(12))
    coul <- rainbow(length(factor.Color))
    PointColors<-coul[1:length(factor.Color)]
    # https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
    # functions that interpolate a set of given colors to create new color palettes
    colfunc <- colorRampPalette(PointColors)
    Mat[, col := colfunc(nrow(Mat))]
    Mat[Index.Good, col := plotrix::color.scale(Mat[Index.Good, xis],c(0,1,1),c(1,1,0),0)]
    #=====[PointSizeModifier]=====
    # If you would like your points to be bigger or smaller, the following lets you adjust their size.
    # '1' represents 100% of the default size, so if you want to make the points larger you could type, for example, 1.5 (i.e. 150%). Similarly, if you want to make them smaller you could type 0.5 (i.e. 50%).
    PointSizeModifier = 1.2
    #=====[GridlineType]=====
    # SPECIFY PARAMETERS FOR GRID #
    # This section lets you add a 'grid' to the plot, running across the main plotting region to help you keep track of where things align with the axes.
    # R will draw lines wherever the 'tick marks' are on both axes. If you wish to adjust these, change 'Xlim' and 'Ylim' above.
    # Choose the type of line you would like to use for your gridlines. You have the following six choices:
    #  "solid"    = -----------------------------------------------             CR
    #  "dashed"   = ----    ----    ----    ----    ----    ----                CR
    #  "dotted"   = -   -   -   -   -   -   -   -   -   -   -   -               CR
    #  "dotdash"  = -   ----   -   ----   -   ----   -   ----   -               CR
    #  "longdash" = -------   -------   -------   -------   -------             CR
    #  "twodash"  = --  ------  --  ------  --  ------  --  ------              CR
    # If you do not want gridlines on your plot, choose "blank" instead.
    GridlineType = "dashed"
    #=====[GridlineColor]=====
    # Choose the color for your grid lines. See step 12 for help.
    GridlineColor = "lightgrey"
    #=====[GridlineWidthModifier]=====
    # If you would like your gridlines to be thicker, the following lets you increase its width.
    # '1' represents 100% of the default thickness, so for example, you could type 2 to make it 200% that thickness.
    GridlineWidthModifier = 1
    #=====[Target Diagram]=====
    op <- par(no.readonly = TRUE)
    # Restoring graphical parameters on exit of function
    par(mar = c(2.6, 2.6, 2.5, 3.5),  # mar=c(margin below, lines at left, lines at top, lines at right)
        mgp = c(1.35, 0.0, 0),         # mgp=c(label, tick mar label, tick mark)
        cex.axis = 0.8)

    on.exit(par(op))
    plot.default(x           = Mat$XData[Index.Good], y    = Mat$YData[Index.Good],
                 xlim        = Xlim                 , ylim = Ylim,
                 xlab        = xAxisLabel           , ylab = yAxisLabel,
                 cex.lab     = 0.8,
                 pch         = PointSymbol,
                 col         = Mat$col[Index.Good],
                 cex         = PointSizeModifier,
                 type        = "p",
                 ann         = TRUE,
                 axes        = TRUE,
                 frame.plot  = TRUE,
                 panel.first = grid(nx  = NULL,
                                    ny  = NULL,
                                    lty = GridlineType,
                                    col = GridlineColor,
                                    lwd = GridlineWidthModifier),
                 panel.last  = NULL,
                 asp         = 1, tck = 0)
    title(main = MainTitle, line = 1.5, font.main = 1, cex.main = 0.9)
    # get the limits of x and y axis
    usr <- par('usr')
    #=====[highest standard deviation]====
    if (sdm_sdo) {
        label.sigma <- " < "
    } else label.sigma <- " > "
    if (abs(Ylim[2]) > abs(Ylim[1])) label.Bias <-  "Bias > 0" else label.Bias <-  "Bias < 0"
    #=====[Text highest standard deviation and bias]=====
    # see https://stackoverflow.com/questions/4973898/combining-paste-and-expression-functions-in-plot-labels
    text(x      = if (sdm_sdo) usr[1]/2 else usr[2]/2,
         y      = if (abs(Ylim[2]) > abs(Ylim[1])) usr[4] else usr[3],
         pos    = if (abs(Ylim[2]) > abs(Ylim[1])) 1 else 3,
         labels = bquote(sigma[Sensor] ~ .(label.sigma) ~ sigma[Reference] ~ " and " ~ .(label.Bias)), cex = 0.8)
    #=====[2nd line Main Title]=====
    if (!is.null(BeginEnd)) {
        Text.Dates <- paste0(BeginEnd[1], " to ", BeginEnd[2])
    } else Text.Dates <- paste0(format(min(Mat$Date, na.rm = T),"%Y%m%d"), " to ", format(max(Mat$Date, na.rm = T),"%Y%m%d"))
    Ref_Analysers <- ifelse(!is.null(Ref_Analysers), paste0(" (",Ref_Analysers,") "),"")
    # Formatting Decimal places in R: https://stackoverflow.com/questions/3443687/formatting-decimal-places-in-r
    if (variable.ubsRM) {
        Text <- substitute(paste(Text.Dates,", Y"["i"]," = b"[0]," + b"[1]," X"["i"],", with b"[0]," = ", b0.Digits, " and b"[1]," = ", b1.Digits, " (", Type.Regression,")",
                                 Ref_Analysers),
                           list(Text.Dates = Text.Dates, b0.Digits = format(round(b0,2), digits = 2), b1.Digits = format(round(b1, digits = 2), nsmall = 2), Type.Regression = Regression, Ref_Analysers = Ref_Analysers))
    } else {
        Text <- substitute(paste(Text.Dates,", Y"["i"]," = b"[0]," + b"[1]," X"["i"],", with b"[0]," = ", b0.Digits, " and b"[1]," = ", b1.Digits, " (", Type.Regression,")",
                                 ", u"[bs_RM]," = ", ubs_RM.Digits, " ", Unit.Ref, Ref_Analysers),
                           list(Text.Dates = Text.Dates, b0.Digits = format(round(b0,2), digits = 2), b1.Digits = format(round(b1, digits = 2), nsmall = 2), Type.Regression = Regression,ubs_RM.Digits = format(round(ubsRM, 2), digits = 2),
                                Unit.Ref = Unit.Ref, Ref_Analysers = Ref_Analysers))
    }
    if (with.ubss) Text <- paste0(Text, substitute(paste(" and u"[bs_s]," = ", ubss.Digits, " ", Unit.sensor), list(ubss.Digits = format(round(ubss,2), digits = 2), Unit.sensor = Unit.sensor)))
    if (Unit.Ref == "ppm" || Unit.Ref == "mg/m3") {
        mtext(text = Text, side = 3, line = 0.1, cex  = 0.65)
    } else {
        if (is.null(Time.average)) {
            mtext(text = Text, side = 3, line = 0.1, cex  = 0.65)
        } else mtext(text = paste0(Text, format(Time.average, digits = 0), " min average"), side = 3, line = 0.1, cex  = 1.0)}
    #=====[Contribution of b1 : (b1 - 1)]=====
    points(x = Mat$b1[Index.Good],
           y = Mat$YData[Index.Good],
           col = Mat$col[Index.Good], pch  = "-", cex  = 0.5)
    #=====[Contribution of b0 : b0/x]=====
    points(x = Mat$b0[Index.Good],
           y = Mat$YData[Index.Good],
           col = Mat$col[Index.Good], pch  = "-", cex  = 0.5)
    #=====[draw axis]=====
    abline(h = 0)
    abline(v = 0)
    #=====[colorscale]=====
    # Add colorscale
    testcol<-plotrix::color.gradient(c(0,1,1),c(1,1,0),0,nslices=100)
    plotrix::color.legend(usr[2] + 0.005 *(usr[2]-usr[1]), usr[3], usr[2] + 0.03 *(usr[2]-usr[1]),usr[4] - 0.05 * (usr[4]-usr[3]),
                          factor.Color,testcol, cex=0.6, align="rb",gradient="y")
    LegendTitle <- substitute(paste("X"[i], ", ", Unit,sep=""),list(Unit=Unit.Ref))
    #=====[target circles]=====
    DQO.step <- min(c(DQO.1, DQO.2 - DQO.1, DQO.3 - DQO.2)) * 100
    Quadrants <- c(sum(usr[c(2,4)]^2),sum(usr[c(1,4)]^2),sum(usr[c(1,3)]^2),sum(usr[2:3]^2))
    #Which.Quadrants <- which.max(Quadrants)
    #usr.Quadrant    <- switch (Which.Quadrants, c(2,4), c(1,4), c(1,3), c(2,3))
    Dist.Quadrants  <- sqrt(max(Quadrants))
    if (DQO.1 * 100 <= Dist.Quadrants && diff(Xlim) != 0 && diff(Ylim) != 0) {
        if (Xlim[2] > Xlim[1]) x = seq(Xlim[1], Xlim[2]) else x = seq(Xlim[2], Xlim[1])
        if (Ylim[2] > Ylim[1]) y = seq(Ylim[1], Ylim[2]) else y = seq(Ylim[2], Ylim[1])
        UR <- function(x,y) {z <- sqrt(x^2+y^2); return(z)}
        z <- outer(x, y, UR)
        Steps <- seq(from = DQO.1 * 100, to = Dist.Quadrants, by = DQO.step)
        Steps <- Steps[-which(as.character(Steps) %in% as.character(c(DQO.1 * 100 , DQO.2  * 100, DQO.3  * 100)))]

        # Plotting intermediary levels
        if (length(Steps) > 0) {
            contour(x = x, y = y, z = z, col = "grey", add = TRUE, method = "edge", levels = Steps,
                    vfont = NULL, axes = FALSE, frame.plot = axes, lty = 'dotdash', lwd = 1, labcex = 0.7)
            # for (Step.level in Steps) {
            #     cercle <- rbind(cbind(seq(-Step.level, Step.level, by =  Step.level/1000),  sqrt(Step.level^2 - seq(-Step.level,Step.level, by = Step.level/1000)^2)),
            #                     cbind(seq( Step.level,-Step.level, by = -Step.level/1000), -sqrt(Step.level^2 - seq(-Step.level,Step.level, by = Step.level/1000)^2)))
            #     lines(cercle,type = "l", lty = GridlineType, col = GridlineColor)
            #     text(x = Step.level/Dist.Quadrants * usr[usr.Quadrant[1]],
            #          y = Step.level/Dist.Quadrants * usr[usr.Quadrant[2]],
            #          paste0(Step.level,"%"),
            #          pos = ifelse(sdm_sdo,2,4),
            #          cex = 0.8)}
        }

        # Plotting DQOs
        contour(x = x, y = y, z = z, col = "black", add = TRUE, method = "edge", levels = c(DQO.1*100, DQO.2*100, DQO.3*100),
                vfont = NULL, axes = axes, frame.plot = axes, lty = 'solid', lwd  = 2, labcex = 0.8)
        # for (i in c(DQO.1*100, DQO.2*100, DQO.3*100)) {
        #     if (!is.na(i)) {
        #         cercle <- rbind(cbind(seq(-i, i, by =  i/1000),  sqrt(i^2 - seq(-i,i, by = i/1000)^2)),
        #                         cbind(seq( i,-i, by = -i/1000), -sqrt(i^2 - seq(-i,i, by = i/1000)^2)))
        #         lines(cercle,type = "l")
        #         text(x = i/Dist.Quadrants * usr[usr.Quadrant[1]], y = i/Dist.Quadrants * usr[usr.Quadrant[2]], paste0(i,"%"), pos = ifelse(sdm_sdo,2,4), cex = 0.8)}}
    }
    #=====[Arrows]=====
    # Plotting segments for relative measurement uncertainty
    # checking that data is visible
    # Select last quartile of these 90 % the biggest Ur within Index.Good
    Index.Big.Ur <- which(Mat[Index.Good,Ur] < quantile(Mat[Index.Good,Ur], probs = 0.9) & Mat[Index.Good,Ur] > quantile(Mat[Index.Good,Ur], probs = 0.75))
    # Absolute difference between Angle of Ur and diagonal 45 and -45 degress
    Mat[, Angle45  := abs(atan(YData/XData) * 180 / pi  -   45 )]
    Mat[, Angle_45 := abs(atan(YData/XData) * 180 / pi  - (-45))]
    # Selecting best point for plotting arrows
    if (min( Mat[Index.Good, Angle45][Index.Big.Ur], na.rm = T) < min(Mat[Index.Good, Angle_45][Index.Big.Ur], na.rm = T)) {
        Index.med.UR <- which(Mat[["Angle45"]] == min(Mat[["Angle45"]][Index.Good][Index.Big.Ur] ))[1]
    } else Index.med.UR <- which(Mat[["Angle_45"]] == min(Mat[["Angle_45"]][Index.Good][Index.Big.Ur]))[1]
    # Index.med.UR based on biggest visible uncertainty
    Index.med.UR <- which(Mat$Ur == Mat$Ur[Index.Good][which.max(Mat[Index.Good][["Ur"]])])[1]
    # plotting Arrows and text Rel.Bias, rel.Rc, Ur
    shape::Arrows(x0 = Mat$XData[Index.med.UR], y0 = 0,
                  x1 = Mat$XData[Index.med.UR], y1 = Mat$YData[Index.med.UR],
                  col = "black",
                  lty = 1, lwd = 1,
                  arr.type = "curved", arr.length = 0.1, code = 3, arr.adj = 1)
    text(x = Mat$XData[Index.med.UR], y = Mat$YData[Index.med.UR]/2,
         labels = c("Bias in %"), pos = 4, srt = 90, cex = 0.8)
    shape::Arrows(x0 = 0                  , y0 = Mat$YData[Index.med.UR],
                  x1 = Mat$XData[Index.med.UR], y1 = Mat$YData[Index.med.UR],
                  col = "black", lty = 1, lwd = 1, arr.type = "curved", arr.length = 0.1,
                  code = 3, arr.adj = 1)
    text(x = Mat$XData[Index.med.UR]/2, y = Mat$YData[Index.med.UR],
         labels = c("Random effect"), pos = 3, adj = 0.5, cex = 0.8)
    shape::Arrows(x0 = 0                  , y0 = 0,
                  x1 = Mat$XData[Index.med.UR], y1 = Mat$YData[Index.med.UR],
                  col = "black", lwd = 2, arr.type = "curved", arr.length = 0.15, code = 3, arr.adj = 1)
    Srt <- atan(Mat$YData[Index.med.UR]/Mat$XData[Index.med.UR]) * 180 /pi
    if (sdm_sdo) {
        if (Srt >= 0) Srt <- - Srt else Srt <- - Srt}
    text(x = Mat$XData[Index.med.UR]/2, y = Mat$YData[Index.med.UR]/2,
         labels = c("Relative expanded uncertainty in %"),
         srt = Srt, pos = 3)
    #=====[Text and line b1 - 1]=====
    text(x = Mat$b1[1]/2,
         y = ifelse(label.Bias ==  "Bias > 0",max(Mat$YData[Index.Good]),min(Mat$YData[Index.Good])),
         labels = bquote("2(b"[1]*"-1)"),
         cex = 0.8, pos = 3) # pos label above the coordinate
    text(x = Mat$b1[1]/2,
         y = ifelse(label.Bias ==  "Bias > 0",max(Mat$YData[Index.Good]),min(Mat$YData[Index.Good])),
         labels = paste0(format(2*(b1 - 1)*100, digits = 1),"%"),
         cex = 0.8, pos = 1) # pos label above the coordinate
    shape::Arrows(x0 = 0         , y0 = ifelse(label.Bias ==  "Bias > 0", max(Mat$YData[Index.Good]), min(Mat$YData[Index.Good])),
                  x1 =  Mat$b1[1], y1 = ifelse(label.Bias ==  "Bias > 0",max(Mat$YData[Index.Good]),min(Mat$YData[Index.Good])),
                  col = "black", lty = 1, lwd = 1, arr.type = "curved", arr.length = 0.1, code = 3, arr.adj = 1)
    #=====[Text and point b0]=====
    Label = bquote("2 b"[0]*"/X"[i]*" in %")
    if (b0 > 0) {
        if (any(Mat$b0 <= ifelse(sdm_sdo,Xlim[1],Xlim[2]))) {
            # There are Mat$b0 <= Xlim[2]
            Index.Good.b0 <- which(Mat[Index.Good, b0] <= ifelse(sdm_sdo,Xlim[1],Xlim[2]))
            x.b0 <- median(Mat[Index.Good, b0][Index.Good.b0], na.rm = FALSE)
            Lower.x.b0 <- which(Mat$b0 <= x.b0)
            if (length(Lower.x.b0) > 0) {
                y.b0 <- Mat[Lower.x.b0][which.max(Mat$b0[Lower.x.b0]), Rel.bias] * 100
                #y.b0 <- Mat[Mat$b0 == x.b0, Rel.bias][1] * 100
                Srt = atan(y.b0 /(x.b0 + Mat$b1[1])) * 180 / pi
                if (sdm_sdo) if (Srt >= 0) Srt <- - Srt else Srt <- - Srt
                text(x = x.b0, y = y.b0, labels = Label, srt = Srt, pos = 3)}}
    } else {
        if (any(Mat$b0 >= ifelse(sdm_sdo,Xlim[2],Xlim[1]))) {
            # There are Mat$b0 >= Xlim[1]
            Index.Good.b0 <- which(Mat[Index.Good, b0] >= ifelse(sdm_sdo,Xlim[2],Xlim[1]))
            x.b0 <- median(Mat[Index.Good, b0][Index.Good.b0], na.rm = FALSE)
            Lower.x.b0 <- which(Mat$b0 <= x.b0)
            if (length(Lower.x.b0) > 0) {
                y.b0 <- Mat[Lower.x.b0][which.max(Mat$b0[Lower.x.b0]), Rel.bias] * 100
                Srt = atan(y.b0 /(x.b0 + Mat$b1[1])) * 180 / pi
                if (sdm_sdo) if (Srt >= 0) Srt <- - Srt else Srt <- - Srt
                text(x = x.b0, y = y.b0, labels = Label, srt = Srt, pos = 4)}}}
    points(x = Mat$b0[Index.med.UR], y = Mat$YData[Index.med.UR], type = "p", col = "black")
    segments(x0 = Mat$b0[Index.med.UR], y0 = 0,
             x1 = Mat$b0[Index.med.UR], y1 = Mat$YData[Index.med.UR],
             col = "black", lty = 2, lwd = 1)
    Label <- substitute(paste("2 b"[0], "/X"[i], ":\n", b0.Digits, "%"), list(b0.Digits = format(Mat$b0[Index.med.UR],digits =0)))
    text(x = Mat$b0[Index.med.UR], y = 0, labels = Label, pos = 1, cex = 0.7)
    #=====[Limit Values identifier]=====
    Limits2plots        <- c(CL, IT, AT, LV, UAT, LAT)
    names(Limits2plots) <- c("CL", "IT", "AT", "LV", "UAT", "LAT")
    for (Limit in names(Limits2plots)) {
        # Checking if Limit is not NA
        if (!is.na(Limits2plots[Limit])) {
            if (Limits2plots[Limit] < max(Mat[Index.Good,xis], na.rm = T)) {
                # checking is Limit is witin plotted x values not necessary, plot Limit is will apper only if within xlim an ylim
                if (any(Mat[Index.Good,xis] == Limits2plots[Limit])) {
                    Index.Limit <- which(Mat[Index.Good,xis] == Limits2plots[Limit])[1]
                } else Index.Limit <- which(abs(Mat[Index.Good,xis] - Limits2plots[Limit]) == min(abs(Mat[Index.Good,xis] - Limits2plots[Limit]), na.rm = T))[1]
                # PLotting a point for Limit
                points(x = Mat$XData[Index.Good][Index.Limit],
                       y = Mat$YData[Index.Good][Index.Limit],
                       type = "p",
                       col  = "blue",
                       cex  = 1.5,
                       pch  = 8,
                       lwd  = 2) #"+")
                # plotting label Limit
                text(x = Mat$XData[Index.Good][Index.Limit],
                     y = Mat$YData[Index.Good][Index.Limit],
                     labels = Limit,
                     pos = 1)}
            if (!exists("Name.Limits2plots")) {
                Name.Limits2plots <- paste0(Limit, " = ", round(Limits2plots[Limit], digits = 1))
            } else Name.Limits2plots <- paste0(Name.Limits2plots, ", ",paste0(Limit, " = ", round(Limits2plots[Limit], digits = 1)))}}
    Blank.space <- 0.01 * (usr[3] - usr[4])
    if (exists("Name.Limits2plots") && length(Name.Limits2plots) > 0 ) text(x      = if (sdm_sdo) (usr[1] + usr[2])/2 else (usr[2] + usr[1])/2,
                                                                            y      = if (abs(Ylim[2]) > abs(Ylim[1])) usr[3]-Blank.space else usr[4]+Blank.space,
                                                                            adj    = c(0.5, ifelse(abs(Ylim[2]) > abs(Ylim[1]), 0, 1)),
                                                                            #pos    = ifelse(abs(Ylim[2]) > abs(Ylim[1]), 3, 1),
                                                                            labels = paste0(Name.Limits2plots," ", Unit.sensor), cex = 0.8)
    # setting all margin accessible with mtext
    par(xpd = NA)
    text(x = usr[2] + 0.005 * (usr[2]-usr[1]), y = usr[4] - 0.005 * (usr[4]-usr[3]), labels  = LegendTitle, adj = c(0,1), cex = 0.8)
}
