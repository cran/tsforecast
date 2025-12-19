##### Generic Forecast Function #####
#' @rdname tsforecast
#' @param object a time series or time series model for which forecasts are required.
#' @param ... additional arguments affecting the forecasts produced.
#' @usage ## S3 method for default
#' tsforecast(object, n.ahead = 1, 
#'   alpha = 0.05, 
#'   show.plot = TRUE, 
#'   forecast.incl = c("all", "forecast", "predict"), 
#'   nobs.incl = NULL, 
#'   log = NULL, 
#'   newxreg = NULL, 
#'   x.name = NULL, pred.name = "Forecasts", ...)
#' @export
tsforecast <- function(object, n.ahead = 1, alpha = 0.05, show.plot = TRUE,
                       forecast.incl = c("all", "forecast", "predict"), nobs.incl = NULL,
                       log = NULL, newxreg = NULL, x.name = NULL, pred.name = "Forecasts", ...)
{
    UseMethod("tsforecast")
}

##### Forecast Time Series Models with Plot #####
#' Forecast Time Series based on Fitted Models
#' @description The generic function `\code{tsforecast}` is used for forecasting time series or time series models.
#' @name tsforecast
#' @rdname tsforecast
#' @param object a time series or time series model for which forecasts are required.
#' @param n.ahead number of forecasting periods. Default is \code{1}.
#' @param alpha significance level. (1 - \code{alpha}) indicates is the confidence level of the prediction intervals. Default is \code{0.05}.
#' @param show.plot logical. If \code{TRUE}, forecasting plot will be displayed directly. Default is \code{TRUE}.
#' @param forecast.incl character string giving which part of the series should be predicted or forecasted. Default value is "\code{all}", which indicates that predicted values of the entire series, together with the forecasted values of the number of future periods specified with \code{n.ahead}, will be included. Specify with "\code{forecast}" if only forecasts of future periods should be included, and "\code{predict}" is used if only fitted values of past periods should be included in the output.
#' @param nobs.incl number of past observations for which the predicted values should be included in the output. If \code{forecast.incl} is set as "\code{forecast}", the value of this parameter will be ignored. Default is \code{NULL},
#' @param log optional. A logical value indicating whether the forecasted values are log-transformed and should be inverted back to the original series scale. If the object is an \code{tsarima} model and this parameter is omitted, the value will be taken over by the settings of the model given in object. Default is \code{NULL} here.
#' @param newxreg new values of the regressors. Only necessary if ARIMA model is built with independent variables.
#' @param x.name name of the series. If omitted here, the series name found by \code{tsname()} will be taken over here. If \code{tsname()} is \code{NULL}, the variable name will be used instead. Default is \code{NULL}.
#' @param pred.name name of the forecasted series here. Default is `\code{Forecasts}`.
#' @param ... additional arguments affecting the forecasts produced.
#' @returns The function \code{print} is used to obtain and print the forecasting results, while the function \code{plot} produces a plot of the forecasts and prediction intervals.
#' @returns An object of class "\code{tsforecast}" is a list usually containing the following elements: 
#' @return \item{x}{original series data. Note that for future periods, no data can be displayed, and the corresponding cells are therefore blank.}
#' @return \item{x.time}{list of time in which the series values were observed.} 
#' @return \item{x.timegap}{time gap between the series and forecasted values.}
#' @return \item{x.name}{name of the time series for which forecasts was requested.}
#' @return \item{pred}{predicted past values and forecasted future values.}
#' @return \item{pred.time}{list of time in which the predictions/forecasts were estimated.} 
#' @return \item{pred.name}{name of the series containing the predicted/forecasted values.}
#' @return \item{se}{standard errors of the forecasted values. Note that standard errors are only calculated for future periods. The corresponding cells for past observations are therefore blank.}
#' @return \item{cil, ciu}{lower and upper limits of the prediction interval. Note that confidence (prediction) intervals are only calculated for future periods. The corresponding cells for past observations are therefore blank.}
#' @return \item{n.ahead}{number of forecasting periods.}
#' @return \item{forecast.incl}{indication of the series part that should be predicted or forecasted.}
#' @return \item{log}{logical. Indicates whether series values are log-transformed for model fitting or not.}
#' @return \item{alpha}{significance level.}
#' @author Ka Yui Karl Wu
#' @references Box, G. E. P., & Jenkins, G. M. (1970). Time series analysis: Forecasting and control. Holden-Day.
#' @references Hyndman, R. J., & Athanasopoulos, G. (2021). Forecasting: Principles and practice (3rd ed.). OTexts. \cr \url{https://otexts.com/fpp3/} 
#' @example inst/examples/tsforecast/tsforecast.R
#' @usage NULL
get_forecast <- function(object, n.ahead = 1, alpha = 0.05, show.plot = TRUE,
                       forecast.incl = c("all", "forecast", "predict"), nobs.incl = NULL,
                       log = NULL, newxreg = NULL, x.name = NULL, pred.name = "Forecasts", ...)
{
    orgmatch <- match.call()
    arglist <- as.list(orgmatch)[-1]
    argnam <- names(arglist)
    x.name <- if (is.null(x.name) & is.null(tsname(object$x))) {object$series} else if (is.null(x.name) & !is.null(tsname(object$x))) {tsname(object$x)} else {x.name}
    forecast.incl <- tolower(match.arg(forecast.incl))
    log <- if (is.null(log) & class(object)[1] == "tsarima") {object$log} else {log <- FALSE}
    if (!is.null(object$xreg))
    {
        ncxreg <- ncol(object$xreg)
        if (ncxreg > 0 & is.null(newxreg))
        {
            colxreg <- colnames(object$xreg)
            newxreg <- NULL
            for (j in 1:ncxreg)
            {
                regj <- tsattrcopy(x = object$xreg[, j], x.orig = object$x.used)
                xreg.arima <- tsarima(regj, order = c(object$arma[1], object$arma[6], object$arma[2]), seasonal = list(order = c(object$arma[3], object$arma[7], object$arma[4]), period = object$arma[5]), include.const = object$include.const)
                xreg.arima.pred <- predict.tsarima(xreg.arima, n.ahead = n.ahead)$pred
                newxreg <- cbind(newxreg, xreg.arima.pred)
            }
            colnames(newxreg) <- colxreg
        }
    }
    fc_arglist <- list(object = object, newxreg = newxreg, n.ahead = n.ahead, alpha = alpha, log = log)
    fc <- if ("tsarima" %in% class(object)) {do.call(predict.tsarima, fc_arglist)} else if ("tsesm" %in% class(object)) {do.call(predict.tsesm, fc_arglist)}
    if (!is.null(nobs.incl) && (nobs.incl == 0 && forecast.incl == "all")) {forecast.incl <- "forecast"}
    if (forecast.incl == "forecast")
    {
        nobs.incl <- 0
        plot.posixdate <- fc$pred.time
        plot.date <- plot.posixdate
        plot.xorig <- plot.xdate <- NULL
        plot.xpred <- fc$pred
        plot.se <- fc$se
        plot.cil <- fc$cil
        plot.ciu <- fc$ciu
    }
    else if (forecast.incl == "predict")
    {
        nobs.incl <- if(is.null(nobs.incl)) {length(object$x)} else if (nobs.incl == 0) {stop("You must at least include some observations in this plot.")} else {nobs.incl}
        n.ahead <- 0
        plot.posixdate <- tail(object$x.time.used, nobs.incl)
        plot.date <- plot.xdate <- plot.posixdate
        plot.xpred <- tail(object[[if (log) {"exp.fitted"} else {"fitted"}]], nobs.incl)
        plot.xorig <- tail(if (log) {exp(object$x.used)} else {object$x.used}, nobs.incl)
        plot.se <- plot.cil <- plot.ciu <- NULL
    }
    else if (forecast.incl == "all")
    {
        nobs.incl <- if(is.null(nobs.incl)) {length(object$x.used)} else {nobs.incl}
        n.plot <- nobs.incl + n.ahead
        cdate <- c(object$x.time.used, fc$pred.time)
        cpred <- c(object[[if (log) {"exp.fitted"} else {"fitted"}]], fc$pred)
        plot.posixdate <- tail(cdate, n.plot)
        plot.date <- plot.posixdate
        plot.xdate <- tail(object$x.time.used, nobs.incl)
        plot.xpred <- tail(cpred, n.plot)
        plot.xorig <- tail(if (log) {exp(object$x.used)} else {object$x.used}, nobs.incl)
        plot.se <- c(rep(NA, nobs.incl), fc$se)
        plot.cil <- c(rep(NA, nobs.incl), fc$cil)
        plot.ciu <- c(rep(NA, nobs.incl), fc$ciu)
    }
    plot.xpred <- tsconvert(x = plot.xpred, t = plot.posixdate)
    if (!forecast.incl == "predict")
    {
        plot.se <- tsconvert(x = plot.se, t = plot.posixdate)
        plot.cil <- tsconvert(x = plot.cil, t = plot.posixdate)
        plot.ciu <- tsconvert(x = plot.ciu, t = plot.posixdate)
    }
    if (!forecast.incl == "forecast")
    {
        plot.xorig <- tsconvert(x = plot.xorig, t = plot.posixdate)
        plot.mainx <- plot.xorig
        plot.maint <- plot.xdate
        plot.predx <- plot.xpred
        plot.predt <- plot.date
    }
    else
    {
        plot.mainx <- NULL
        plot.maint <- NULL
        plot.predx <- plot.xpred
        plot.predt <- plot.date
    }
    tgap <- if (n.ahead == 1) {tsfreq(object$x.used)} else {tstimegap(plot.date)}
    out <- c(list(x = object$x, x.used = plot.mainx, x.time = plot.maint, x.timegap = tgap, x.name = x.name, pred = plot.predx, pred.time = plot.predt, pred.name = pred.name, se = plot.se, cil = plot.cil, ciu = plot.ciu, n.ahead = n.ahead, forecast.incl = forecast.incl, log = log, alpha = alpha))
    if (show.plot)
    {
        plot_arg <- c("ylim", "title", "x.lwidth", "pred.lwidth", "x.col", "pred.col", "ci.col")
        user_plot_arglist <- arglist[plot_arg[plot_arg %in% argnam]]
        plot_arglist <- c(list(x = out), user_plot_arglist)
        do.call(plot.tsforecast, args = plot_arglist)
    }
    return(structure(out, class = "tsforecast"))
}

##### Forecast ARIMA Models with Plot #####
#' @rdname tsforecast
#' @exportS3Method 
tsforecast.tsarima <- function(object, ...)
{
    args <- as.list(match.call())
    args$object <- object
    do.call(get_forecast, args[-1L])
}

##### Forecast ESM with Plot #####
#' @rdname tsforecast
#' @exportS3Method 
tsforecast.tsesm <- function(object, ...)
{
    args <- as.list(match.call())
    args$object <- object
    do.call(get_forecast, args[-1L])
}

##### Print ARIMA Forecasts #####
#' @rdname tsforecast
#' @param x an object of class `\code{tsforecast}`.
#' @exportS3Method 
print.tsforecast <- function(x, ...)
{
    args <- as.list(match.call())
    args$object <- x
    do.call(tsfprint, args[-1L])
}

##### Plot Forecasted Values #####
#' @rdname tsforecast
#' @param ylim value limit of the y-axis. The values should be specified in \code{c(lower_limit, upper_limit)}, where \code{lower_limit} and \code{upper_limit} are the values of the smallest and largest number of the y-axis, respectively. Default is \code{NULL}.
#' @param title title of the forecasting chart. Default is \code{NULL}.
#' @param x.lwidth line width of the original series in the output plot. Default is \code{0.7}.
#' @param pred.lwidth line width of the predicted past values or forecasted future values in the output plot. Default is \code{0.7}.
#' @param x.col colour of the original series in the output plot. Default is `\code{darkgrey}`.
#' @param pred.col colour of the predicted past values or forecasted future values in the output plot. Default is `\code{red}`.
#' @param ci.col colour of the forecasted prediction interval in the output plot. Default is `\code{royalblue}`.
#' @exportS3Method 
plot.tsforecast <- function(x, ylim = NULL, title = NULL, x.lwidth = 0.7, pred.lwidth = 0.7, 
                            x.col = "darkgrey", pred.col = "red", ci.col = "royalblue", ...)
{
    if (!x$forecast.incl == "forecast")
    {
        plot.colx <- x.col
        plot.colpred <- pred.col
        plot.widthx <- x.lwidth
        plot.widthpred <- pred.lwidth
        plot.mainx <- x$x.used
        plot.maint <- x$x.time
        plot.predx <- x$pred
        plot.predt <- x$pred.time
    }
    else
    {
        plot.colx <- pred.col
        plot.colpred <- pred.col
        plot.widthx <- pred.lwidth
        plot.widthpred <- pred.lwidth
        plot.mainx <- x$pred
        plot.maint <- x$pred.time
        plot.predx <- NULL
        plot.predt <- x$pred.time
    }
    plot.cit <- x$pred.time[-(1:length(x$x.time))]
    plot.cil <- x$cil
    plot.ciu <- x$ciu
    x.name <- x$x.name
    pred.name <- x$pred.name
    tslineplot(x = plot.mainx, t = plot.maint, pred = plot.predx, pred.t = plot.predt, cil = plot.cil, ciu = plot.ciu, 
               ci.t = plot.predt, x.name = x.name, pred.name = if (is.null(pred.name)) {"Forecasts"} else {pred.name}, 
               title = if (is.null(title)) {paste0("Forecasts of ", x.name)} else {title}, ylim = ylim, 
               x.lwidth = plot.widthx, pred.lwidth = plot.widthpred, x.col = plot.colx, pred.col = plot.colpred, ci.col = ci.col)
}

##### Print Forecasts of TS Models #####
tsfprint <- function(object, digits = max(3L, getOption("digits") - 3L), se = TRUE, head = NULL, tail = NULL, print.index = NULL, ...)
{
    cilevel <- paste0("CI(", (1 - object$alpha) * 100, "%)")
    #x.col <- if (is.null(object$x.used)) {NULL} else {if (object$log) {exp(object$x.used)} else {object$x.used}}
    x.col <- if (is.null(object$x.used)) {NULL} else {object$x.used}
    if ((!is.null(object$x.time)) & (!is.null(object$pred.time)))
    {
        if (length(object$pred.time) - length(object$x.time) != 0 & (!is.null(x.col)))
        {
            x.col <- c(x.col, rep(NA, length(object$pred.time) - length(object$x.time)))
        }
    }
    x.col <- if (!is.null(x.col)) {list(x.col)}
    se.col <- if (se && !is.null(object$se)) {list(object$se)} else {NULL}
    cil.col <- if (is.null(object$cil)) {NULL} else {list(object$cil)}
    ciu.col <- if (is.null(object$ciu)) {NULL} else {list(object$ciu)}
    flist <- c(x = x.col, if (!is.null(object$pred)) {pred = list(object$pred)}, se = se.col, cil = cil.col, ciu = ciu.col)
    ftab <- data.frame(flist)
    tabtime <- if (is.null(object$pred.time)) {object$x.time} else {object$pred.time}
    fdate <- tstimeformat(tabtime, object$x.timegap)
    rownames(ftab) <- fdate
    colnames(ftab) <- c(if (length(object$x.used) > 0 && !is.null(x.col)) {"Original"}, "Forecast", if (se && !is.null(object$se)) {"Std. Error"}, if (!is.null(object$cil)) {paste(cilevel, "Lower")}, if (!is.null(object$ciu)) {paste(cilevel, "Upper")})
    ptab <- round(if (!is.null(tail)) {tail(ftab, tail)} else if (!is.null(head)) {head(ftab, head)} else if (!is.null(print.index)) {ftab[print.index, ]} else {ftab}, digits = digits)
    ptab[is.na(ptab)] <- ""
    return(print(ptab, digits = digits, print.gap = 3))
}

##### Generate Forecasting Periods #####
fperiod <- function(x, n.ahead = 1)
{
    speriod <- tstime(x)
    pperiod <- seq.POSIXt(from = tail(speriod$time, 1L), length.out = n.ahead + 1L, by = speriod$frequency)
    lperiod <- pperiod[-1L]
    return(lperiod)
}

##### Check whether an Observation is an Outlier #####
#' Outlier Identification
#' @description The function `\code{is.outlier}` checks whether any time series observations are outliers based on the interquartile range (IQR) rule.
#' @param x a time series or any other R data type.
#' @param method method based on which the outliers are identified. Available options are `\code{iqr}`, `\code{sigma}`, and `\code{zscore}`.
#' @param param parameter value for setting specific boundary criteria. Default is \code{NULL}.
#' @author Ka Yui Karl Wu
#' @details With \code{method = "iqr"}, the interquartile range rule for outlier identification is applied. An observation \eqn{x_i} will be identified as outlier if one of the following conditions fulfils:
#' @details \deqn{x_i < q_1 - m \cdot (q_3 - q_1)} 
#' @details \deqn{x_i > q_3 + m \cdot (q_3 - q_1)}
#' @details where \eqn{q_1} and \eqn{q_3} are the 1st and 3rd quartiles of the time series \code{x}, respectively. \code{m} is the value specified by \code{param}. If omitted, it will be set as 1.5.
#' @details By using \code{method = "sigma"}, the following criteria for outlier identification, known as the 3-sigma rule, are applied:
#' @details \deqn{x_i < \mu(x) - m \cdot \sigma(x)} 
#' @details \deqn{x_i > \mu(x) + m \cdot \sigma(x)}
#' @details where \eqn{\mu(x)} and \eqn{\sigma(x)} are the mean and standard deviation of the time series \code{x}, respectively. \code{m} is the value specified by \code{param}. If omitted, it will be set as 3.
#' @details The z-score rule, specified by \code{method = "zscore"}, compares the standardised observation values to a specific threshold:
#' @details \deqn{\left|\dfrac{(x_i - \mu(x))}{\sigma(x)}\right| > m} 
#' @details where \eqn{\mu(x)} and \eqn{\sigma(x)} are the mean and standard deviation of the time series \code{x}, respectively. \code{m} is the value specified by \code{param}. If omitted, it will be set as 2. Note that 2 is the threshold for mild outliers. If checking for extreme outliers is required, the value should be set as 3.
#' @return A vector indicating whether the values in \code{x} are outliers (\code{TRUE}) or not (\code{FALSE}).
#' @example inst/examples/is.outlier/is.outlier.R
#' @export
is.outlier <- function(x, method = c("iqr", "sigma", "zscore"), param = NULL)
{
    method <- tolower(match.arg(method))
    if (method == "iqr")
    {
        if (is.null(param)) {param <- 1.5}
        x.q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
        bound <- x.q + c(-param, param) * diff(x.q)
    }
    else if (method == "sigma")
    {
        if (is.null(param)) {param <- 3}
        x.sigma <- sd(x, na.rm = TRUE)
        bound <- mean(x, na.rm = TRUE) + c(-param, param) * x.sigma
    }
    else if (method == "zscore")
    {
        attrx <- attributes(x)
        if (is.null(param)) {param <- 2}
        x <- (x[!is.na(x)] - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
        attributes(x) <- attrx
        bound <- c(-param, param)
    }
    return(x <= bound[1L] | x >= bound[2L])
}

##### Calculate the Optimal Time Breaks for Plots #####
optimTimeBreaks <- function(len, frequency, numbreaks)
{
    if (len > frequency & len > numbreaks)
    {
        if (len / frequency > frequency / 2L)
        {
            base <- len / numbreaks
            stepwidth <- round(base / frequency, 0L)
            if (stepwidth == 0L) stepwidth <- 1L
            step <- stepwidth * frequency
        }
        else
        {
            step <- ceiling(len / frequency)
        }
    }
    else
    {
        step <- 1L
    }
    return(seq(from = 1L, to = len, by = step))
}

##### Fit ARIMA Models #####
#' Fitting ARIMA Models
#' @description The `\code{tsarima}` function is used to fit an ARIMA model to a univariate time series.
#' @name tsarima
#' @rdname tsarima
#' @param x a univariate time series or an `\code{tsarima}` object.
#' @param order a specification of the non-seasonal part of the ARIMA model: the three integer components \eqn{(p, d, q)} are the AR order, the degree of differencing, and the MA order.
#' @param seasonal a specification of the seasonal part of the ARIMA model \eqn{(P, D, Q)}, the seasonal AR order, the degree of seasonal differencing, and the seasonal MA order, plus the period (which defaults to \code{frequency(x)}). This should be a list with components \code{order} and \code{period}, but a specification of just a numeric vector of length 3 will be turned into a suitable list with the specification as the order. 
#' @param xreg optional. A vector or matrix of external regressors, which must have the same number of rows as \code{x}.
#' @param include.const logical. Indicates if the ARMA model should include a mean/intercept term. The default is \code{TRUE} for non-differenced series. For ARIMA models with differencing, it may fail to estimate the standard errors.
#' @param log optional. A logical value indicating whether the forecasted values are log-transformed and should be inverted back to the original series scale. If the object is an \code{tsarima} model and this parameter is omitted, the value will be taken over by the settings of the model given in object. Default is \code{NULL} here.
#' @param train.prop a numerical value specifying the proportion of training data in the series. The value must be between 0 and 1. Default is \code{1}.
#' @param arch.test optional. A logical value indicating whether the ARCH effect in the residuals should be tested by the McLeod-Li test of not. Default is \code{FALSE}.
#' @param transform.pars logical. If \code{TRUE}, the AR parameters are transformed to ensure that they remain in the region of stationarity. Not used for \code{method = "CSS"}. For \code{method = "ML"}, it has been advantageous to set \code{transform.pars = FALSE} in some cases, see also \code{fixed}.
#' @param fixed optional. Numeric vector of the same length as the total number of coefficients to be estimated. It should be of the form 
#' \deqn{(\phi_1,\ldots,\phi_p,\theta_1,\ldots,\theta_q,\Phi_1,\ldots,\Phi_P,\Theta_1,\ldots,\Theta_Q,\mu)} 
#' where \eqn{\phi_i} are the AR coefficients, \eqn{\theta_i} are the MA coefficients, \eqn{\Phi_i} are the seasonal AR coefficients, \eqn{\Theta_i} are the seasonal MA coefficients and \eqn{\mu} is the intercept term. 
#' \cr Note that the \eqn{\mu} entry is required if and only if \code{include.const} is \code{TRUE}. In particular it should not be present if the model is an ARIMA model with differencing.
#' \cr The entries of the \code{fixed} vector should consist of the values at which the user wishes to `fix` the corresponding coefficient, or \code{NA} if that coefficient should \emph{not} be fixed, but estimated.
#' \cr The argument \code{transform.pars} will be set to \code{FALSE} if any AR parameters are fixed. A warning will be given if \code{transform.pars} is set to (or left at its default) \code{TRUE}. It may be wise to set \code{transform.pars = FALSE} even when fixing MA parameters, especially at values that cause the model to be nearly non-invertible.
#' @param init optional. Numeric vector of initial parameter values. Missing values will be filled in, by zeroes except for regression coefficients. Values already specified in \code{fixed} will be ignored.
#' @param method fitting method. Maximum likelihood or minimize conditional sum-of-squares. The default (unless there are missing values) is to use conditional-sum-of-squares to find starting values, then maximum likelihood. Can be abbreviated.
#' @param SSinit a string specifying the algorithm to compute the state-space initialization of the likelihood; see \code{\link{KalmanLike}} for details. Can be abbreviated.
#' @param optim.method The value passed as the method argument to \code{\link{optim}}.
#' @param optim.control List of control parameters for \code{\link{optim}}.
#' @param kappa the prior variance (as a multiple of the innovations variance) for the past observations in a differenced model. Do not reduce this.
#' @author Ka Yui Karl Wu
#' @details Different definitions of ARMA models have different signs for the AR and/or MA coefficients. The definition used here is the original Box & Jenkins (1970) formulation:
#' @details \eqn{x_t=\phi_1 x_{t-1}+\ldots+\phi_p x_{t-p}+\varepsilon_t-\theta_1\varepsilon_{t-1}-\ldots-\theta_p\varepsilon_{t-q}}
#' @details and so the MA coefficients differ in sign from the output of \code{stats::arima}. Further, if \code{include.const} is \code{TRUE} (the default for an ARMA model), this formula applies to \eqn{x_t-\mu} rather than \eqn{x_t}. For ARIMA models with differencing, the differenced series usually follows a zero-mean ARMA model, but \code{include.const} is still available in case a constant term is required for the model. However, the estimation of the coefficients' standard error may not be successful. If an \code{xreg} term is included, a linear regression (with a constant term if \code{include.mean} is \code{TRUE} and there is no differencing) is fitted with an ARMA model for the error term.
#' @details The variance matrix of the estimates is found from the Hessian of the log-likelihood, and so may only be a rough guide.
#' @details Optimization is done by \code{\link{optim}}. It will work best if the columns in \code{xreg} are roughly scaled to zero mean and unit variance, but does attempt to estimate suitable scalings.
#' @details If \code{train.prop} is smaller than 1, the function will only treat the training part of the series as past data. When applying `\code{tsforecast}` or `\code{predict}`, the forecast will start after the end of the training part of the original series.
#' @returns A list of class `\code{tsarima}` with components: 
#' @return \item{coef}{a vector of AR, MA and regression coefficients, which can be extracted by the \code{\link{coef}} method.}
#' @return \item{const}{a value of the model's constant term. Return \code{NULL} if \code{include.const = FALSE}.}
#' @return \item{sigma2}{the maximum likelihood estimate of the white noise variance.}
#' @return \item{var.coef}{he estimated variance matrix of the coefficients \code{coef}, which can be extracted by the \code{vcov} method.}
#' @return \item{loglik}{the maximized log-likelihood (of the differenced data), or the approximation to it used.}
#' @return \item{aic, aicc, bic}{the AIC, AICc, and BIC values corresponding to the log-likelihood. Only valid for method = `ML` fits.}
#' @return \item{error}{a list of prediction error estimators, including \code{$ME} for mean error, \code{$RMSE} for root mean squared error, \code{$MAE} for mean absolute error, \code{$MPE} for mean percentage error, \code{$MAPE} for mean absolute percentage error, \code{$MASE} for mean absolute scaled error, \code{$MASE.S} for seasonal mean absolute scaled error, and \code{$ACF1} for lag 1 autocorrelation.}
#' @return \item{arma}{a vector of the ARIMA order: \eqn{(p, P, q, Q, \ell, d, D)}.}
#' @return \item{train.prop}{proportion of training data.}
#' @return \item{x}{data of the original series.}
#' @return \item{x.time}{list of time in which the series values were observed.} 
#' @return \item{x.timegap}{time gap between the series and forecasted values.}
#' @return \item{x.name}{name of the time series for which forecasts was requested.}
#' @return \item{x.time.used}{list of time in which the series values were used for model fitting. It will be the same as \code{x.time} if \code{train.prop = 1}.}
#' @return \item{x.used}{data of the original series which were used for model fitting. It will be the same as \code{x} if \code{train.prop = 1}.}
#' @return \item{fitted}{a vector of fitted series values. If \code{train.prop} is smaller than 1, it will have the same length as \code{x.used}.}
#' @return \item{residuals}{a vector of the series residuals. If \code{train.prop} is smaller than 1, it will have the same length as \code{x.used}.}
#' @return \item{exp.fitted}{a vector of fitted series values after inverted from log-transformation. Only available if \code{log = TRUE}. If \code{train.prop} is smaller than 1, it will have the same length as \code{x.used}.}
#' @return \item{exp.residuals}{a vector of the series residuals after inverted from log-transformation. Only available if \code{log = TRUE}. If \code{train.prop} is smaller than 1, it will have the same length as \code{x.used}.}
#' @return \item{log}{logical. Indicates whether series values are log-transformed for model fitting or not.}
#' @return \item{call}{the matched call.}
#' @return \item{series}{series name \code{x} in match call.}
#' @return \item{code}{the convergence value returned by \code{\link{optim}}.}
#' @return \item{nobs}{the number of `used` observations for the fitting, can also be extracted via \code{\link{nobs}} and is used by \code{\link{BIC}}.}
#' @return \item{model}{a list representing the Kalman filter used in the fitting. See \code{\link{KalmanLike}}.}
#' @return \item{mcleod.li.test}{resulting chi-square test statistics and the corresponding p-values of the McLeod-Li test for ARCH effect. Only available if \code{arch.test = TRUE}.}
#' @return \item{model.test}{a list of information regarding the prediction of the testing data including `\code{x.test}` (part of `\code{x}` used for testing), `\code{fitted.test}` (predicted values of the testing data), `\code{residuals.test}` (prediction error of the testing data), and `\code{error.test}` (prediction error measurements based on the testing data). Only available if \code{train.prop} is smaller than 1.}
#' @section Fitting methods: 
#' The exact likelihood is computed via a state-space representation of the ARIMA process, and the innovations and their variance found by a Kalman filter. The initialization of the differenced ARMA process uses stationarity and is based on Gardner et. al. (1980). For a differenced process the non-stationary components are given a diffuse prior (controlled by \code{kappa}). Observations which are still controlled by the diffuse prior (determined by having a Kalman gain of at least \code{1e4}) are excluded from the likelihood calculations. (This gives comparable results to \code{\link{arima0}} in the absence of missing values, when the observations excluded are precisely those dropped by the differencing.)
#' \cr\cr Missing values are allowed, and are handled exactly in method `\code{ML}`.
#' \cr\cr If \code{transform.pars} is \code{TRUE}, the optimisation is done using an alternative parametrization which is a variation on that suggested by Jones (1980) and ensures that the model is stationary. For an AR(p) model the parametrisation is via the inverse tanh of the partial autocorrelations: the same procedure is applied (separately) to the AR and seasonal AR terms. The MA terms are not constrained to be invertible during optimisation, but they will be converted to invertible form after optimisation if \code{transform.pars} is \code{TRUE}.
#' \cr\cr Conditional sum-of-squares is provided mainly for expositional purposes. This computes the sum of squares of the fitted innovations from observation \code{n.cond} on, (where \code{n.cond} is at least the maximum lag of an AR term), treating all earlier innovations to be zero. Argument \code{n.cond} can be used to allow comparability between different fits. The `part log-likelihood` is the first term, half the log of the estimated mean square. Missing values are allowed, but will cause many of the innovations to be missing.
#' \cr\cr When regressors are specified, they are orthogonalised prior to fitting unless any of the coefficients is fixed. It can be helpful to roughly scale the regressors to zero mean and unit variance.
#' @example inst/examples/tsarima/tsarima.R
#' @references Box, G. E. P., & Jenkins, G. M. (1970). Time series analysis: Forecasting and control. Holden-Day.
#' @references Hyndman, R. J., & Athanasopoulos, G. (2021). Forecasting: Principles and practice (3rd ed.). OTexts. \cr \url{https://otexts.com/fpp3/} 
#' @references Hyndman, R. J., Athanasopoulos, G., Bergmeir, C., Caceres, G., Chhay, L., O'Hara-Wild, M., Petropoulos, F., Razbash, S., Wang, E., & Yasmeen, F. (2025). \emph{forecast: Forecasting functions for time series and linear models}. R package version 8.24.0, \cr \url{https://pkg.robjhyndman.com/forecast/}.
#' @references Hyndman, R. J., & Khandakar, Y. (2008). Automatic time series forecasting: the forecast package for R. \emph{Journal of Statistical Software}, \strong{27}(3), 1-22. \doi{10.18637/jss.v027.i03}.
#' @references Brockwell, P. J., & Davis, R. A. (1996). Introduction to Time Series and Forecasting. Springer, New York. Sections 3.3 and 8.3.
#' @references Durbin, J., & Koopman, S. J. (2001). Time Series Analysis by State Space Methods. Oxford University Press.
#' @references Gardner, G, Harvey, A. C., & Phillips, G. D. A. (1980). Algorithm AS 154: An algorithm for exact maximum likelihood estimation of autoregressive-moving average models by means of Kalman filtering. \emph{Applied Statistics}, \strong{29}, 311-322. \doi{10.2307/2346910}.
#' @references Harvey, A. C. (1993). Time Series Models. 2nd Edition. Harvester Wheatsheaf. Sections 3.3 and 4.4.
#' @references Jones, R. H. (1980). Maximum likelihood fitting of ARMA models to time series with missing observations. \emph{Technometrics}, \strong{22}, 389-395. \doi{10.2307/1268324}.
#' @references Ripley, B. D. (2002). Time series in R 1.5.0. \emph{R News}, \strong{2}(2), 2-7. \url{https://www.r-project.org/doc/Rnews/Rnews_2002-2.pdf}
#' @seealso \link{arima}
#' @importFrom stats KalmanForecast
#' @importFrom stats makeARIMA
#' @importFrom stats KalmanLike
#' @importFrom stats KalmanRun
#' @importFrom stats KalmanSmooth
#' @importFrom stats ts
#' @importFrom stats frequency
#' @importFrom stats deltat
#' @importFrom stats coef
#' @importFrom stats na.fail
#' @importFrom stats na.omit
#' @importFrom stats na.exclude
#' @importFrom stats na.pass
#' @importFrom stats as.ts
#' @importFrom stats optim
#' @importFrom stats arima
#' @export
tsarima <- function(x, order = c(0L, 0L, 0L), seasonal = list(order = c(0L, 0L, 0L), period = NA), 
                    xreg = NULL, include.const = TRUE, log = FALSE, train.prop = 1, arch.test = FALSE, 
                    transform.pars = TRUE, fixed = NULL, init = NULL, method = c("CSS-ML", "ML", "CSS"), SSinit = c("Gardner1980", "Rossignol2011"),
                    optim.method = "BFGS", optim.control = list(), kappa = 1e+06)
{
    if (train.prop > 1 | train.prop < 0)
    {
        warning("Proportion of training data must be between 0 and 1. Automatically set to 1.")
        train.prop <- 1
    }
    method <- match.arg(method)
    norig <- length(x)
    xfreq <- frequency(x)
    xdate <- tstime(x)
    series <- deparse1(substitute(x))
    if (!is.null(xreg)) 
    {
        if (!is.numeric(xreg)) 
        {
            stop("xreg should be a numeric matrix or a numeric vector")
        }
        xreg <- as.matrix(xreg)
        if (is.null(colnames(xreg))) 
        {
            colnames(xreg) <- if (ncol(xreg) == 1) {"xreg"} else {paste("xreg", 1:ncol(xreg), sep = "")}
        }
        xnam <- colnames(xreg)
    }
    else
    {
        xnam <- NULL
    }
    if (!is.list(seasonal)) 
    {
        if (xfreq <= 1) 
        {
            seasonal <- list(order = c(0, 0, 0), period = NA)
            if (norig <= order[2L]) {stop("Not enough data to fit the model")}
        }
        else 
        {
            seasonal <- list(order = seasonal, period = xfreq)
            if (norig <= order[2L] + seasonal$order[2L] * seasonal$period) {stop("Not enough data to fit the model")}
        }
    }
    if (include.const) 
    {
        intercept <- driftx(x, order = order[2], order.D = seasonal$order[2], period = if(is.null(seasonal$period)) {xfreq} else {seasonal$period})
        xreg_model <- cbind(xreg, intercept)
        colnames(xreg_model) <- c(xnam, "intercept")
    }
    else
    {
        xreg_model <- xreg    
    }
    if (log) {xreg_model <- log(xreg_model)}
    if (include.const == TRUE && (order[2] + seasonal$order[2] > 0))
    {
        SSinit <- "Rossignol2011"
        optim.method = "L-BFGS-B"
    }
    xmodel <- if (log) {log(x)} else {x}
    trainlen <- round(train.prop * length(xmodel), 0)
    xdate_train <- xdate$time[1:trainlen]
    xmodel_train <- tsattrcopy(x = xmodel[1:trainlen], x.orig = x)
    xreg_modeltrain <- xreg_model[1:trainlen, , drop = FALSE]
    arima.fit <- stats::arima(x = xmodel_train, order = order, seasonal = seasonal, xreg = xreg_modeltrain, include.mean = FALSE, transform.pars = transform.pars, fixed = fixed, init = init, method = method, SSinit = SSinit, optim.method = optim.method, optim.control = optim.control, kappa = kappa)
    coef <- arima.fit$coef
    coef[(grepl("ma", names(coef)))] <- -coef[(grepl("ma", names(coef)))]
    const <- if (arima.fit$arma[1] + arima.fit$arma[3] > 0 && include.const == TRUE) {(1 - sum(coef[(grepl("ar", names(coef)))])) * coef["intercept"]} else {NULL}
    fitted <- xmodel_train - arima.fit$residuals
    attributes(arima.fit$residuals) <- attributes(xmodel_train)
    attributes(fitted) <- attributes(xmodel_train)
    error <- tsmodeleval(list(x = xmodel_train, fitted = fitted))
    value <- -2 * arima.fit$loglik
    aic <- arima.fit$aic
    aicc <- aic + 2 * sum(arima.fit$mask) * (arima.fit$nobs/(arima.fit$nobs - sum(arima.fit$mask) - 1) - 1)
    bic <- aic + sum(arima.fit$mask) * (log(arima.fit$nobs) - 2)
    x.name <- if (is.null(tsname(x))) {series} else {tsname(x)}
    out <- c(list(coef = coef, const = const, sigma2 = arima.fit$sigma2, var.coef = arima.fit$var, mask = arima.fit$mask, loglik = arima.fit$loglik,
                  aic = aic, aicc = aicc, bic = bic, error = error, arma = arima.fit$arma, train.prop = train.prop, 
                  x = x, x.time = xdate$time, x.timegap = xdate$frequency, x.name = x.name, x.time.used = xdate_train, x.used = xmodel_train, fitted = fitted, residuals = arima.fit$residuals, xreg = xreg, xreg.used = xreg_modeltrain[, xnam]),
             if (log) {list(exp.residuals = exp(arima.fit$residuals), exp.fitted = exp(fitted))},
             list(log = log, include.const = include.const, call = match.call(), series = series, code = arima.fit$code, n.cond = arima.fit$n.cond, nobs = arima.fit$nobs, model = arima.fit$model))
    if (arch.test) {out$mcleod.li.test <- tsmltest(object = out, lag.max = min(10 * log10(length(xmodel_train)), 8))}
    if (train.prop < 1)
    {
        testlen <- length(xmodel) - trainlen
        xdate_test <- xdate$time[(trainlen + 1):length(xmodel)]
        xmodel_test <- tsconvert(x = xmodel[(trainlen + 1):length(xmodel)], t = xdate_test, x.name = tsname(x))
        xreg_modeltest <- if (!is.null(xreg)) {xreg_model[(trainlen + 1):length(xmodel), , drop = FALSE]} else {NULL}
        testval <- predict.tsarima(out, n.ahead = testlen, newxreg = xreg_modeltest[, xnam], se.fit = FALSE, log = FALSE)$pred
        testres <- as.numeric(xmodel_test) - as.numeric(testval)
        attributes(testval) <- attributes(xmodel_test)
        attributes(testres) <- attributes(xmodel_test)
        error_test <- tsmodeleval(list(x = xmodel_test, fitted = testval))
        out <- c(out, list(model.test = list(x.test = xmodel_test, fitted.test = testval, residuals.test = testres, error.test = error_test, xreg.test = xreg_modeltest[, xnam])))
    }
    out <- structure(out, class = c("tsarima"))
    return(out)
}

##### Print ARIMA Models #####
#' @rdname tsarima
#' @param digits the number of significant digits.
#' @param se logical. If \code{TRUE}, standard error will be included in displaying the result. Default is \code{TRUE}.
#' @param signif.stars logical. If \code{TRUE}, a shorthand used to indicate the statistical significance of a result will be displayed next to the p-values with *** for p < 0.001, ** for p < 0.01, * for p < 0.05, and . for p < 0.1. Default is \code{TRUE}.
#' @param ... other printing or summary parameters.
#' @importFrom stats printCoefmat
#' @importFrom stats pt
#' @exportS3Method 
print.tsarima <- function(x, digits = max(3L, getOption("digits") - 3L), se = TRUE, signif.stars = TRUE, ...)
{
    cat("\nCall:", deparse(x$call, width.cutoff = 75L), "", sep = "\n")
    if (length(x$coef))
    {
        cat("Coefficients:\n")
        coef <- round(x$coef, digits = digits)
        if (se && NROW(x$var.coef)) {
            ses <- rep.int(0, length(coef))
            ses[x$mask] <- round(sqrt(diag(x$var.coef)), digits = digits)
            tval <- abs(coef / ses)
            pval <- round(pt(abs(coef / ses), df = x$nobs - length(coef), lower.tail = FALSE) * 2, digits = digits)
            coef <- cbind(t(matrix(coef, 1L, dimnames = list(NULL, names(coef)))), s.e. = ses, tval = round(tval, digits = digits), pval = pval)
            colnames(coef) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
            coef <- printCoefmat(coef, print.gap = 2)
        }
        else
        {
            print.default(coef, print.gap = 2)
        }
    }
    if (!is.null(x$const))
    {
        cat("\nConstant estimate: ", x$const)
    }
    cm <- x$call$method
    if (is.null(cm) || cm != "CSS")
        cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits),
            ":  log likelihood = ", format(round(x$loglik, 2L)),
            ",  aic = ", format(round(x$aic, 2L)), "\n\n", sep = "")
    else cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits),
             ":  part log likelihood = ", format(round(x$loglik, 2)),
             "\n\n", sep = "")
    invisible(x)
}

##### Summarise ARIMA Models #####
#' @rdname tsarima
#' @param object a \code{tsarima} object for summary.
#' @param ... other printing or summary parameters.
#' @exportS3Method 
summary.tsarima <- function(object, digits = max(3L, getOption("digits") - 3L), se = TRUE, signif.stars = TRUE, ...)
{
    cat("\nSeries:", if (is.null(tsname(object$x))) {object$series} else {tsname(object$x)}, "\n")
    sarima <- sum(object$arma[c(3L, 4L, 7L)])
    if (sarima > 0)
    {
        modelnam <- "SARIMA"
        modelorder <- paste0("(", paste(object$arma[3], object$arma[7], object$arma[4], sep = ", "), ")[", object$arma[5], "]")
    }
    else
    {
        modelnam <- "ARIMA"
        modelorder <- ""
    }
    modelorder <- paste0("(", paste(object$arma[1], object$arma[6], object$arma[2], sep = ", "), ")", modelorder)
    cat("Model: ", paste(modelnam, modelorder), "\n\n")
    if (length(object$coef)) {
        cat("Coefficients:\n")
        coef <- round(object$coef, digits = digits)
        if (se && NROW(object$var.coef)) {
            ses <- rep.int(0, length(coef))
            ses[object$mask] <- round(sqrt(diag(object$var.coef)), digits = digits)
            tval <- abs(coef / ses)
            pval <- round(pt(abs(coef / ses), df = object$nobs - length(coef), lower.tail = FALSE) * 2, digits = digits)
            coef <- cbind(t(matrix(coef, 1L, dimnames = list(NULL, names(coef)))), s.e. = ses, tval = round(tval, digits = digits), pval = pval)
            colnames(coef) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
            coef <- printCoefmat(coef, print.gap = 2)
        }
        else
        {
            print.default(coef, print.gap = 2)
        }
    }
    if (!is.null(object$const))
    {
        cat("\nConstant estimate: ", object$const)
    }
    cm <- object$call$method
    if (is.null(cm) || cm != "CSS")
        cat("\nsigma^2 estimated as ", format(object$sigma2, digits = digits),
            ":  log likelihood = ", format(round(object$loglik, 2L)),
            "\nAIC = ", format(round(object$aic, 2L)), "   AICc = ", format(round(object$aicc, 2L)), "   BIC = ", format(round(object$bic, 2L)), "\n\n", sep = "")
    else cat("\nsigma^2 estimated as ", format(object$sigma2, digits = digits),
             ":  part log likelihood = ", format(round(object$loglik, 2)),
             "\n\n", sep = "")
    if (!is.null(object$mcleod.li.test))
    {
        cat("Test of ARCH effects:\n")
        arch.tab <- t(as.data.frame(object$mcleod.li.test))
        rownames(arch.tab) <- c("chi^2 value", "Pr(> chi^2)")
        colnames(arch.tab) <- seq(1, ncol(arch.tab))
        print(arch.tab, digits = digits, print.gap = 3)
        cat("\n")
    }
    cat("Error measures:\n")
    print(tsmodeleval(object), digits = digits, print.gap = 2)
    cat("\n")
    invisible(object)
}

##### Predict ARIMA Models #####
#' @rdname predict
#' @param newxreg new values of the regressors. Only necessary if ARIMA model is built with independent variables.
#' @param se.fit logical. If \code{TRUE}, standard error of each prediction will be calculated and included. Default is \code{TRUE}.
#' @param log optional. A logical value indicating whether the forecasted values are log-transformed and should be inverted back to the original series scale. If the object is an \code{tsarima} model and this parameter is omitted, the value will be taken over by the settings of the model given in object. Default is \code{NULL} here.
#' @exportS3Method 
predict.tsarima <- function (object, n.ahead = 1L, newxreg = NULL, se.fit = TRUE, alpha = 0.05, log = NULL, ...)
{
    xr <- object$call$xreg
    xreg <- if (!is.null(xr)) {eval.parent(xr)} else {NULL}
    if (!is.null(newxreg)) 
    {
        ncxreg <- ncol(xreg)
        if (ncol(newxreg) != ncxreg) {stop("'xreg' and 'newxreg' have different numbers of columns")}
        if (!is.numeric(newxreg)) {stop("newxreg should be numeric.")}
        newxreg <- as.matrix(newxreg)
    }
    else
    {
        ncxreg <- 0
    }
    n <- length(object$x)
    arma <- object$arma
    coefs <- object$coef
    coefs[(grepl("ma", names(coefs)))] <- -coefs[(grepl("ma", names(coefs)))]
    userlog <- if (is.null(log)) {FALSE} else {log}
    narma <- sum(arma[1L:4L])
    if (length(coefs) > narma) 
    {
        if (object$include.const) 
        {
            newconst <- tail(driftx(rep(1, n + n.ahead), order = arma[6L], order.D = arma[7L], period = arma[5L]), n.ahead)
            newxreg <- cbind(newxreg, mean = newconst)
            ncxreg <- ncxreg + 1L
        }
        if (userlog == TRUE | (is.null(log) & object$log == TRUE))
        {
            newxreg <- log(newxreg)
        }
        xm <- if (narma == 0) {drop(as.matrix(newxreg) %*% coefs)} else {drop(as.matrix(newxreg) %*% coefs[-(1L:narma)])}
    }
    else 
    {
        xm <- 0
    }
    if (arma[2L] > 0L) {
        ma <- coefs[arma[1L] + 1L:arma[2L]]
        if (any(Mod(polyroot(c(1, ma))) < 1))
            warning("MA part of model is not invertible")
    }
    if (arma[4L] > 0L) {
        ma <- coefs[sum(arma[1L:3L]) + 1L:arma[4L]]
        if (any(Mod(polyroot(c(1, ma))) < 1))
            warning("seasonal MA part of model is not invertible")
    }
    z <- KalmanForecast(n.ahead, object$model)
    ldate <- fperiod(object$x.used, n.ahead = n.ahead)
    pred <- tsconvert(x = z[[1L]] + xm, t = ldate)
    if (se.fit)
    {
        se <- tsconvert(x = sqrt(z[[2L]] * object$sigma2), t = ldate)
        df <- object$nobs - length(coefs)
        pred.err <- se * abs(qt(1L - alpha / 2L, df = df))
        cil <- pred - pred.err
        ciu <- pred + pred.err
        out <- list(pred = pred, se = se, cil = cil, ciu = ciu)
    }
    else
    {
        out <- list(pred = pred)
    }
    if (userlog == TRUE | (is.null(log) & object$log == TRUE))
    {
        userlog = TRUE
        out$pred <- exp(out$pred)
        if (se.fit)
        {
            out$cil <- exp(out$cil)
            out$ciu <- exp(out$ciu)
            out$se <- (out$ciu - out$cil) / (2 * abs(qt(1L - alpha / 2L, df = df)))
        }
    }
    tgap <- if (n.ahead == 1) {tsfreq(object$x)} else {tstimegap(ldate)}
    return(structure(c(list(pred.time = ldate, x.timegap = tgap, x.name = tsname(object$x)), out, list(n.ahead = n.ahead, log = userlog, alpha = alpha)), class = "tspredict"))
}

##### Generate Exponential Smoothing #####
#' Exponential Smoothing Forecasts
#' @description The `\code{tsesm}` function forecasts future values of a univariate time series using exponential smoothing.
#' @name tsesm
#' @rdname tsesm
#' @param x a univariate time series or a `\code{tsesm}` object.
#' @param order a specification of the exponential smoothing method. The available options are "\code{simple}", "\code{holt}", "\code{holt-winters}".
#' @param damped logical. If \code{TRUE}, a damped trend is used. Default is \code{FALSE}.
#' @param initial method used for selecting initial state values. Available options are `\code{optimal}` and `\code{simple}`. If `\code{optimal}`, the initial values are optimised along with the smoothing parameters using `\code{ets}`. If `\code{simple}`, the initial values are set to values obtained using simple calculations on the first few observations.
#' @param type specify whether the time series is `\code{additive}` or `\code{multiplicative}`.
#' @param alpha value of smoothing parameter for the level. If \code{NULL}, it will be estimated.
#' @param beta value of smoothing parameter for the trend. If \code{NULL}, it will be estimated.
#' @param gamma value of smoothing parameter for the seasonal component If \code{NULL}, it will be estimated.
#' @param lambda parameter of the Box-Cox transformation. If \code{lambda = "auto"}, a transformation is automatically selected using `\code{BoxCox.lambda}`. The transformation is ignored if \code{NULL}. Otherwise, data transformed before model is estimated.
#' @param phi value of damping parameter if \code{damped = TRUE}. If \code{NULL}, it will be estimated.
#' @param biasadj use adjusted back-transformed mean for Box-Cox transformations. If transformed data is used to produce forecasts and fitted values, a regular back transformation will result in median forecasts. If \code{biasadj == TRUE}, an adjustment will be made to produce mean forecasts and fitted values.
#' @param exp.trend logical. If \code{TRUE}, an exponential trend is fitted. Otherwise, the trend is (locally) linear. Default is \code{FALSE}.
#' @param seasonal.period number of seasons within a seasonal cycle. If \code{NULL}, the value of \code{frequency(x)} will be taken. Default is \code{NULL}.
#' @param train.prop a numerical value specifying the proportion of training data in the series. The value must be between 0 and 1. Default is \code{1}.
#' @author Ka Yui Karl Wu
#' @details The `\code{tsesm}` function uses the same mechanism for exponential smoothing like `\code{forecast::ses}`, `\code{forecast::holt}`, and `\code{forecast::hw}`. Instead of using three different functions, all of them are integrated in \code{tsesm}, and user can choose the method by specifying the value of the \code{order} parameter.
#' @details The option `\code{simple}` is the lowest exponential smoothing order, or the first exponential smoothing order (that's why the parameter is called `\code{order}` and not `method`). It can be used to forecast series with only level \eqn{(\ell)} component. The corresponding forecasting formula is given by:
#' @details \deqn{\tilde{x}_{t+h|t} = \ell_t = \alpha x_{t-1} + (1-\alpha)\ell_{t-1},}
#' @details where \eqn{\alpha} is the smoothing parameter for the level component, and \eqn{h} is the number of periods ahead for forecasting.
#' @details With the `\code{holt}` option, the second exponential smoothing order can be chosen. It is suitable for the forecasting of series with level \eqn{(\ell)} and trend \eqn{(b)} components.
#' @details \deqn{\ell_t = \alpha x_{t} + (1-\alpha)(\ell_{t-1} + b_{t-1})}
#' @details \deqn{b_t = \beta(\ell_{t}-\ell_{t-1}) + (1-\beta)b_{t-1}}
#' @details \deqn{\tilde{x}_{t+h|t} = \ell_{t}+hb_{t}}
#' @details where \eqn{\alpha} and \eqn{\beta} are the smoothing parameters for the level and trend components, respectively, \eqn{h} is the number of periods ahead for forecasting.
#' @details The third exponential smoothing order can be specified by the `\code{holt-winters}` option. It can be used to forecast time series with level \eqn{(\ell)}, trend \eqn{(b)}, and seasonal \eqn{(s)} components.
#' @details \deqn{\ell_t = \alpha (x_{t} - s_{t-m}) + (1-\alpha)(\ell_{t-1} + b_{t-1})}
#' @details \deqn{b_t = \beta(\ell_{t}-\ell_{t-1}) + (1-\beta)b_{t-1}}
#' @details \deqn{s_t = \gamma(x_{t}-\ell_{t-1}-b_{t-1}) + (1-\gamma)s_{t-m}}
#' @details \deqn{\tilde{x}_{t+h|t} = \ell_{t}+hb_{t}+s_{t+h-m(k+1)}}
#' @details where \eqn{\alpha}, \eqn{\beta}, and \eqn{\gamma} are the smoothing parameters for the level, trend, and seasonal components, respectively, \eqn{h} is the number of periods ahead for forecasting, \eqn{m} is the number of periods within a seasonal cycle, and \eqn{k} is the integer part of \eqn{(h-1)/m}, which ensures that the estimates of the seasonal indices used for forecasting come from the final period of the series.
#' @details If \code{train.prop} is smaller than 1, the function will only treat the training part of the series as past data. When applying `\code{tsforecast}` or `\code{predict}`, the forecast will start after the end of the training part of the original series.
#' @returns A list of class "\code{tsesm}" with components: 
#' @return \item{coef}{a vector of smoothing parameters.}
#' @return \item{m}{number of seasons in the series, usually equivalent to the series' frequency obtained by \code{tsfreq}.}
#' @return \item{components}{values used for fitting exponential smoothing model.}
#' @return \item{states}{estimated values of all smoothing parameter for each observation.}
#' @return \item{initstate}{initial values of the smoothing parameters.}
#' @return \item{sigma2}{residual variance.}
#' @return \item{loglik}{the maximized log-likelihood, or the approximation to it used.}
#' @return \item{aic, aicc, bic}{the AIC, AICc, and BIC values corresponding to the log-likelihood. Only valid for method = "ML" fits.}
#' @return \item{x}{original series data or a `\code{tsesm}` object.}
#' @return \item{x.time}{list of time in which the series values were observed.} 
#' @return \item{x.timegap}{time gap between the series and forecasted values.}
#' @return \item{x.name}{name of the time series for which forecasts was requested.}
#' @return \item{fitted}{a vector of fitted series values.}
#' @return \item{residuals}{a vector of the series residuals.}
#' @return \item{damped}{logical value indicating whether damped trend was used or not.}
#' @return \item{initial}{method used for selecting initial state values.}
#' @return \item{type}{indicator of an additive or multiplicative time series.}
#' @return \item{exp.trend}{indicator of the use of exponential trend.}
#' @return \item{lambda}{parameter value of the Box-Cox transformation.}
#' @return \item{biasadj}{indicator of the use of adjusted back-transformed mean for Box-Cox transformations.}
#' @return \item{series}{series name \code{x} in match call.}
#' @return \item{call}{the matched call.}
#' @return \item{error}{a list of prediction error estimators, including \code{$ME} for mean error, \code{$RMSE} for root mean squared error, \code{$MAE} for mean absolute error, \code{$MPE} for mean percentage error, \code{$MAPE} for mean absolute percentage error, \code{$MASE} for mean absolute scaled error, \code{$MASE.S} for seasonal mean absolute scaled error, and \code{$ACF1} for lag 1 autocorrelation.}
#' @return \item{model.test}{a list of information regarding the prediction of the testing data including `\code{x.test}` (part of `\code{x}` used for testing), `\code{fitted.test}` (predicted values of the testing data), `\code{residuals.test}` (prediction error of the testing data), and `\code{error.test}` (prediction error measurements based on the testing data). Only available if \code{train.prop} is smaller than 1.}
#' @example inst/examples/tsesm/tsesm.R
#' @references Box, G. E. P., & Jenkins, G. M. (1970). Time series analysis: Forecasting and control. Holden-Day.
#' @references Hyndman, R. J., Athanasopoulos, G. (2021). Forecasting: Principles and practice (3rd ed.). OTexts. \cr \url{https://otexts.com/fpp3/} 
#' @references Hyndman, R. J., Athanasopoulos, G., Bergmeir, C., Caceres, G., Chhay, L., O'Hara-Wild, M., Petropoulos, F., Razbash, S., Wang, E., Yasmeen, F. (2025). \cr \emph{forecast: Forecasting functions for time series and linear models}. R package version 8.24.0, \cr \url{https://pkg.robjhyndman.com/forecast/}.
#' @references Hyndman, R. J., Khandakar, Y. (2008). Automatic time series forecasting: the forecast package for R. \emph{Journal of Statistical Software}, \strong{27}(3), 1-22.
#' @import forecast
#' @importFrom forecast ses
#' @importFrom forecast holt
#' @importFrom forecast hw
#' @importFrom stats start
#' @export
tsesm <- function(x, order = c("simple", "holt", "holt-winters"), damped = FALSE, initial = c("optimal", "simple"), type = c("additive", "multiplicative"),
                  alpha = NULL, beta = NULL, gamma = NULL, lambda = NULL, phi = NULL, biasadj = FALSE, exp.trend = FALSE, seasonal.period = NULL, train.prop = 1)
{
    initial <- match.arg(initial)
    order <- match.arg(order)
    type <- match.arg(type)
    series <- deparse1(substitute(x))
    xdate <- tstime(x)
    xmodel <- x
    if (!is.null(seasonal.period))
    {
        if(seasonal.period != frequency(xmodel)) 
        {
            xmodel <- ts(x, start = start(x), frequency = seasonal.period)
        }
    }
    trainlen <- round(train.prop * length(xmodel), 0)
    xdate_train <- xdate$time[1:trainlen]
    xmodel_train <- tsattrcopy(x = xmodel[1:trainlen], x.orig = xmodel)
    if (order == "simple")
    {
        arglist <- list(y = xmodel_train, initial = initial, alpha = alpha, lambda = lambda, biasadj = biasadj)
    }
    else if (order == "holt")
    {
        arglist <- list(y = xmodel_train, damped = damped, initial = initial, exponential = exp.trend, alpha = alpha, beta = beta, phi = phi, lambda = lambda, biasadj = biasadj)
    }
    else if (order == "holt-winters")
    {
        arglist <- list(y = xmodel_train, seasonal = type, damped = damped, initial = initial, exponential = exp.trend, alpha = alpha, beta = beta, gamma = gamma, phi = phi, lambda = lambda, biasadj = biasadj)
    }
    esm.fun <- if (order == "simple") {"ses"} else if (order == "holt") {"holt"} else if (order == "holt-winters") {"hw"}
    esm <- suppressWarnings(suppressMessages(do.call(esm.fun, arglist)))
    esm$model$states <- ts(esm$model$states, start = start(x) - deltat(x), frequency = frequency(x))
    parincl <- c("par", "m", "components", "states", "initstate", "sigma2")
    if (initial == "optimal") {parincl <- c(parincl, "loglik", "aic", "bic", "aicc")}
    xmodel_train <- tsattrcopy(xmodel_train, x.orig = x)
    xres <- tsattrcopy(esm$residuals, x.orig = x)
    xfit <- tsattrcopy(esm$fitted, x.orig = x)
    x.name <- if (is.null(tsname(x))) {deparse1(substitute(x))} else {tsname(x)}
    esm.out <- esm$model[parincl]
    names(esm.out)[names(esm.out) == "par"] <- "coef"
    out <- c(esm.out, list(x = x, x.time = xdate$time, x.timegap = xdate$frequency, x.name = x.name, 
                           train.prop = train.prop, x.used = xmodel_train, x.time.used = xdate_train, fitted = xfit, residuals = xres,
                           damped = damped, initial = initial, type = type, exp.trend = exp.trend, lambda = lambda, biasadj = biasadj, 
                           series = series, call = match.call(), error = tsmodeleval(list(x = xmodel_train, fitted = xfit))))
    if (train.prop < 1)
    {
        testlen <- length(xmodel) - trainlen
        xdate_test <- xdate$time[(trainlen + 1):length(xmodel)]
        xmodel_test <- tsconvert(x = xmodel[(trainlen + 1):length(xmodel)], t = xdate_test, x.name = tsname(x))
        testval <- predict.tsesm(out, n.ahead = testlen, se.fit = FALSE)$pred
        testres <- as.numeric(xmodel_test) - as.numeric(testval)
        attributes(testval) <- attributes(xmodel_test)
        attributes(testres) <- attributes(xmodel_test)
        error_test <- tsmodeleval(list(x = xmodel_test, fitted = testval))
        out <- c(out, list(model.test = list(x.test = xmodel_test, fitted.test = testval, residuals.test = testres, error.test = error_test)))
    }
    return(structure(out, class = c("tsesm")))
}

##### Print Exponential Smoothing Models #####
#' @rdname tsesm
#' @param ... other printing parameters
#' @exportS3Method 
print.tsesm <- function(x, ...)
{
    orgmatch <- match.call()
    arglist <- as.list(orgmatch)[-1]
    argnam <- names(arglist)
    digits <- if ("digits" %in% argnam) {arglist$digits} else {3L}
    scipenval <- getOption("scipen")
    options(scipen = 999)
    on.exit(options(scipen = scipenval))
    cat("\nCall:", deparse(x$call, width.cutoff = 75L), "", sep = "\n")
    if (length(x$coef))
    {
        cat("Smoothing Parameter:\n")
        smparnam <- c("alpha", "beta", "gamma", "phi", "lambda")
        smpar <- x$coef[names(x$coef) %in% smparnam]
        if (length(smpar))
        {
            for (j in names(smpar))
            {
                if (!is.na(smpar[j]))
                {
                    val <- round(smpar[j], digits = digits)
                    cat(j, rep(" ", 13 - nchar(j)), format(val, width = nchar(val), justify = "right"), "\n", sep = "")
                }
            }
        }
        if (!is.null(x$lambda))
        {
            cat("\nBox-Cox transformation parameter:\n")
            cat("lambda", rep(" ", 7), format(round(x$lambda, digits = digits), width = nchar(round(x$lambda, digits = digits)), justify = "right"), "\n", sep = "")
        }
        if (length(x$initstate))
        {
            cat("\nInitial Values: \n")
            hwpar <- x$initstate
            names(hwpar)[names(hwpar) == "l"] <- "level"
            names(hwpar)[names(hwpar) == "b"] <- "trend"
            names(hwpar) <- gsub("s", "season ", names(hwpar))
            maxvalwidth <- max(nchar(as.character(round(hwpar, digits = digits)))) + 2
            for (k in names(hwpar))
            {
                cat(k, rep(" ", 10 - nchar(k)), format(round(hwpar[k], digits = digits), width = maxvalwidth, justify = "right"), "\n", sep = "")
            }
        }
    }
    if (x$initial == "optimal")
    {
        cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits),
            ":  log likelihood = ", format(round(x$loglik, 2L)),
            "\nAIC = ", format(round(x$aic, 2L)), "   AICc = ", format(round(x$aicc, 2L)), "   BIC = ", format(round(x$bic, 2L)),
            "\n\n", sep = "")
    }
    else
    {
        cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits), "\n\n", sep = "")
    }
    invisible(x)
}

##### Summarise Exponential Smoothing Models #####
#' @rdname tsesm
#' @param object a `\code{tsesm}` object to summarise.
#' @param ... other printing parameters
#' @exportS3Method 
summary.tsesm <- function(object, ...)
{
    orgmatch <- match.call()
    arglist <- as.list(orgmatch)[-1]
    argnam <- names(arglist)
    digits <- if ("digits" %in% argnam) {arglist$digits} else {3L}
    scipenval <- getOption("scipen")
    options(scipen = 999)
    on.exit(options(scipen = scipenval))
    cat("\nCall:", deparse(object$call, width.cutoff = 75L), "", sep = "\n")
    if (length(object$coef))
    {
        cat("Smoothing Parameter:\n")
        smparnam <- c("alpha", "beta", "gamma", "phi", "lambda")
        smpar <- object$coef[names(object$coef) %in% smparnam]
        if (length(smpar))
        {
            for (j in names(smpar))
            {
                if (!is.na(smpar[j]))
                {
                    val <- round(smpar[j], digits = digits)
                    cat(j, rep(" ", 13 - nchar(j)), format(val, width = nchar(val), justify = "right"), "\n", sep = "")
                }
            }
        }
        if (!is.null(object$lambda))
        {
            cat("\nBox-Cox transformation parameter:\n")
            cat("lambda", rep(" ", 7), format(round(object$lambda, digits = digits), width = nchar(round(object$lambda, digits = digits)), justify = "right"), "\n", sep = "")
        }
        if (length(object$initstate))
        {
            cat("\nInitial Values: \n")
            hwpar <- object$initstate
            names(hwpar)[names(hwpar) == "l"] <- "level"
            names(hwpar)[names(hwpar) == "b"] <- "trend"
            names(hwpar) <- gsub("s", "season ", names(hwpar))
            maxvalwidth <- max(nchar(as.character(round(hwpar, digits = digits)))) + 2
            for (k in names(hwpar))
            {
                cat(k, rep(" ", 10 - nchar(k)), format(round(hwpar[k], digits = digits), width = maxvalwidth, justify = "right"), "\n", sep = "")
            }
        }
    }
    if (object$initial == "optimal")
    {
        cat("\nsigma^2 estimated as ", format(object$sigma2, digits = digits),
            ":  log likelihood = ", format(round(object$loglik, 2L)),
            "\nAIC = ", format(round(object$aic, 2L)), "   AICc = ", format(round(object$aicc, 2L)), "   BIC = ", format(round(object$bic, 2L)),
            "\n\n", sep = "")
    }
    else
    {
        cat("\nsigma^2 estimated as ", format(object$sigma2, digits = digits), "\n\n", sep = "")
    }
    cat("Error measures:\n")
    print(tsmodeleval(object), digits = digits, print.gap = 2)
    cat("\n")
    invisible(object)
}

##### Predict ESM Models #####
#' Predict Time Series Values
#' @description The function `\code{predict}` is generic and predicts past/future values of a time series.
#' @name predict
#' @rdname predict
#' @param object a time series or time series model for which prediction is required.
#' @param n.ahead number of forecasting periods. Default is \code{1}.
#' @param alpha significance level. (1 - \code{alpha}) indicates is the confidence level of the prediction interval. Default is \code{0.05}.
#' @param ... additional arguments affecting the forecasts produced.
#' @returns An object of class "\code{tspredict}".
#' @returns The function \code{print} is used to obtain and print the prediction results, including the predicted values, the corresponding standard errors, as well as the lower and upper limit of the prediction intervals.
#' @returns An object of class "\code{tspredict}" is a list usually containing the following elements: 
#' @return \item{x.time}{list of time in which the series values were observed.} 
#' @return \item{x.timegap}{time gap between the series and forecasted values.}
#' @return \item{x.name}{name of the time series for which forecasts was requested.}
#' @return \item{pred}{predicted past values and forecasted future values.}
#' @return \item{se}{standard errors of the forecasted values.}
#' @return \item{cil, ciu}{lower and upper limits of the prediction interval.}
#' @return \item{n.ahead}{number of forecasting periods.}
#' @return \item{log}{logical. Indicates whether series values are log-transformed for model fitting or not. (Only available for class "\code{tsarima}")}
#' @return \item{alpha}{significance level.}
#' @author Ka Yui Karl Wu
#' @references Box, G. E. P., & Jenkins, G. M. (1970). Time series analysis: Forecasting and control. Holden-Day.
#' @references Hyndman, R. J., & Athanasopoulos, G. (2021). Forecasting: Principles and practice (3rd ed.). OTexts. \cr \url{https://otexts.com/fpp3/} 
#' @example inst/examples/predict/predict.R
#' @importFrom stats qt
#' @importFrom forecast forecast.ets
#' @exportS3Method 
predict.tsesm <- function(object, n.ahead = 1L, se.fit = TRUE, alpha = 0.05, ...)
{
    names(object)[names(object) == "coef"] <- "par"
    if (object$train.prop < 1) {object$x <- object$x.used}
    f <- forecast.ets(object, h = n.ahead, level = (1 - alpha) * 100)
    ldate <- fperiod(object$x.used, n.ahead = n.ahead)
    df <- length(object$x) - length(object$coef)
    mu <- tsconvert(x = f$mean, t = ldate)
    if (se.fit)
    {
        se <- tsconvert(abs((f$mean - f$lower) / qnorm(1 - alpha / 2)), t = ldate)
        pred.err <- se * abs(qt(1L - alpha / 2L, df = df))
        cil <- tsconvert(x = mu - pred.err, t = ldate)
        ciu <- tsconvert(x = mu + pred.err, t = ldate)
    }    
    tgap <- if (n.ahead == 1) {tsfreq(object$x.used)} else {tstimegap(ldate)}
    out <- c(list(pred.time = ldate, x.timegap = tgap, x.name = tsname(object$x), pred = mu), if (se.fit) {list(se = se, cil = cil, ciu = ciu)}, list(n.ahead = n.ahead, alpha = alpha))
    return(structure(out, class = "tspredict"))
}

##### Print Prediction #####
#' @rdname predict
#' @param x a `\code{tspredict}` object.
#' @exportS3Method 
print.tspredict <- function(x, ...)
{
    args <- as.list(match.call())
    args$object <- x
    do.call(tsfprint, args[-1L])
}

##### Generate ACF #####
#' Auto- Covariance and -Correlation Function Estimation
#' @description The function `\code{tsacf}` computes (and by default plots) estimates of the autocorrelation function. Function \code{pacf} is the function used for the partial autocorrelations, and function \code{tsacov} computes the autocovariance function.
#' @name tsacf
#' @rdname tsacf
#' @param x,y univariate time series object(s) or a numeric vector(s) or matrix/matrices, or a `\code{tsacf}` object.
#' @param lag.max maximum lag at which to calculate the acf. Default is \code{8}. Will be automatically limited to one less than the number of observations in the series. If the series has less than 8 observations.
#' @param type character string giving the type of acf to be computed. Allowed values are "correlation" (the default), "covariance", "partial", "cross-correlation", or "cross-covariance". Will be partially matched.
#' @param show.plot logical. If \code{TRUE}, the acf/pacf/acov will be plotted. Default is \code{TRUE}.
#' @param na.action function to be called to handle missing values. Default is \code{na.omit}. 
#' @param demean logical. If \code{TRUE}, the covariances will be about the sample means. Default is \code{TRUE}.
#' @param alpha significance level. (1 - \code{alpha}) indicates is the confidence level of the prediction interval. Default is \code{0.05}.
#' @param x.name name of the series. If omitted here, the series name found by \code{tsname()} will be taken over here. If \code{tsname()} is \code{NULL}, the variable name will be used instead. Default is \code{NULL}. 
#' @param title character string indicating the plot title.
#' @author Ka Yui Karl Wu
#' @details For \code{type = "correlation"} and \code{"covariance"}, the estimates are based on the sample covariance of \eqn{x_t} and \eqn{x_{t-k}} (lag 0 autocorrelation is fixed at 1 by convention.). For \code{"cross-correlation"} and \code{"cross-covariance"}, the estimates are based on the sample covariance of \eqn{x_t} and \eqn{y_{t-k}}.
#' @details By default, no missing values are allowed. However, by default, \code{na.action = na.omit}, the covariances are computed only from complete cases. This means that the estimate computed may well not be a valid autocorrelation sequence, and may contain missing values. Missing values are not allowed when computing the PACF of a multivariate time series.
#' @details The partial correlation coefficient is estimated by fitting autoregressive models of successively higher orders up to \code{lag.max}.
#' @details The lag is returned and plotted in units of time, and not numbers of observations.
#' @details Different from \code{\link{acf}}, the lags here are not converted based on the seasonal cycle length. It simply reflects the time lag \eqn{k} between \eqn{X_t} and \eqn{X_{t-k}}. Furthermore, the ACF/PACF/ACOV/CCF/CCOV plots are created using \code{ggplot2}.
#' @details The generic functions \code{plot} and \code{print} have both methods for objects of class "\code{tsacf}".
#' @details For `\code{cross-correlation}` and `\code{cross-covariance}`, positive lags indicate that y is shifted forward in time relative to x, while negative lags indicate y is shifted backward.
#' @returns An object of class "\code{tsacf}", which is a list with the following elements:
#' @return \item{plot}{bar chart of the estimated ACF/PACF/ACOV. Only available if \code{save.plot} is \code{TRUE}.}
#' @return \item{acf, pacf, acov, ccf, ccov}{An array with the same dimensions as lag containing the estimated ACF/PACF/ACOV/CCF/CCOV.}
#' @return \item{clim}{upper limits of the estimated confidence intervals for each lag. Since confidence intervals are symmetrical around 0, the lower limits are simply the negative values of \code{clim}.}
#' @return \item{type}{type of correlation (same as the \code{type} argument).}
#' @return \item{n.used}{number of observations in the time series.}
#' @return \item{lag}{lags at which the acf is estimated.}
#' @return \item{x.name}{name of the time series for which forecasts was requested.}
#' @return \item{alpha}{significance level.}
#' @example inst/examples/tsacf/tsacf.R
#' @references Box, G. E. P., & Jenkins, G. M. (1970). Time series analysis: Forecasting and control. Holden-Day.
#' @references Hyndman, R. J., & Athanasopoulos, G. (2021). Forecasting: Principles and practice (3rd ed.). OTexts. \cr \url{https://otexts.com/fpp3/} 
#' @references Venables, W. N., & Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth Edition. Springer-Verlag.
#' @usage ## Autocorrelation Function (ACF)
#' tsacf(
#'   x, y = NULL, lag.max = 8, 
#'   type = c("correlation", "covariance", "partial", 
#'            "cross-correlation", "cross-covariance"), 
#'   show.plot = TRUE, na.action = na.omit, 
#'   demean = TRUE, alpha = 0.05, 
#'   x.name = NULL, title = NULL)
#' @importFrom stats acf
#' @export
tsacf <- function(x, y = NULL, lag.max = 8, type = c("correlation", "covariance", "partial", "cross-correlation", "cross-covariance"), show.plot = TRUE, na.action = na.omit, demean = TRUE, alpha = 0.05, x.name = NULL, title = NULL)
{
    type <- match.arg(type)
    x.name <- if (is.null(x.name) & is.null(tsname(x))) {deparse1(substitute(x))} else if (is.null(x.name) & !is.null(tsname(x))) {tsname(x)} else {x.name}
    lag.max <- min(length(x) - 1L, lag.max)
    ccftype <- c("cross-correlation", "cross-covariance")
    lag <- (if (type == "partial") {1L} else if (type %in% ccftype) {-lag.max} else {0L}):lag.max
    if (type %in% ccftype)
    {
        fun <- "ccf"
        ctype <- if (type == "cross-correlation") {"correlation"} else if (type == "cross-covariance") {"covariance"}
        acfdat <- list(x = x, y = y)
    }
    else
    {
        fun <- "acf"
        ctype <- type
        acfdat <- list(x = x)
    }
    x.acf <- do.call(fun, args = c(acfdat, list(plot = FALSE, type = ctype, lag.max = lag.max, na.action = na.action)))
    acf.res <- as.vector(x.acf$acf)
    if (type == "correlation" & !x.acf$lag[1L] == 0L) {acf.res <- c(1L, acf.res)}
    names(acf.res) <- lag
    if (type == "correlation")
    {
        rho <- acf.res[2L:(length(acf.res) - 1L)]
        bartapprox <- 1 + cumsum(2 * c(0, rho) ^ 2)
        x.clim <- qnorm(1 - alpha / 2) * sqrt(bartapprox / x.acf$n.used)
        x.clim <- c(NA, x.clim)
        names(x.clim) <- 0:lag.max
    }
    else if (type == "partial")
    {
        x.clim <- rep(qnorm(1 - alpha / 2) * sqrt(1 / x.acf$n.used), lag.max)
        names(x.clim) <- names(acf.res)
    }
    else if (type == "cross-correlation")
    {
        rhoxy <- acf.res[2L:(length(acf.res) - 1L)]
        xacf <- acf(x, plot = FALSE, lag.max = lag.max, na.action = na.action)$acf
        rhox <- xacf[2L:(length(xacf) - 1L)]
        yacf <- acf(y, plot = FALSE, lag.max = lag.max, na.action = na.action)$acf
        rhoy <- yacf[2L:(length(yacf) - 1L)]
        x.clim <- rep(qnorm(1 - alpha / 2) * sqrt(1 / x.acf$n.used), length(lag))
        names(x.clim) <- (-lag.max):lag.max
    }
    if (ctype != "covariance")
    {
        x.clim[abs(x.clim) > 1] <- 1L
    }
    else
    {
        x.clim <- NULL
    }
    acf_list <- if (type == "correlation") {list(acf = acf.res)} else if (type == "partial") {list(pacf = acf.res)} else if (type == "covariance") {list(acov = acf.res)} else if (type == "cross-covariance") {list(ccov = acf.res)} else if (type == "cross-correlation") {list(ccf = acf.res)}
    out <- c(acf_list, list(clim = x.clim, type = type, n.used = x.acf$n.used, lag = lag, x.name = x.name, alpha = alpha))
    if (show.plot) {plot.tsacf(out, title = title)}
    structure(out, class = c("tsacf"))
}

##### Generate PACF #####
#' @rdname tsacf
#' @param ... same parameters used in \code{tsacf}.
#' @usage ## Partial Autocorrelation Function (PACF)
#' tspacf(...)
#' @export
tspacf <- function(...)
{
    tsacf(type = "partial", ...)
}

##### Generate Autocovariance #####
#' @rdname tsacf
#' @param ... same parameters used in \code{tsacf}.
#' @usage ## Autocovariance (ACOV)
#' tsacov(...)
#' @export
tsacov <- function(...)
{
    tsacf(type = "covariance", ...)
}

##### Generate CCF #####
#' @rdname tsacf
#' @param ... same parameters used in \code{tsacf}.
#' @usage ## cross-correlation (CCF)
#' tsccf(...)
#' @export
tsccf <- function(...)
{
    tsacf(type = "cross-correlation", ...)
}

##### Generate Cross-covariance #####
#' @rdname tsacf
#' @param ... same parameters used in \code{tsacf}.
#' @usage ## cross-ovariance (CCOV)
#' tsccov(...)
#' @export
tsccov <- function(...)
{
    tsacf(type = "cross-covariance", ...)
}

##### Print tsacf Object #####
#' @rdname tsacf
#' @param digits number of decimal digits displayed in the results.
#' @param ... other printing or plotting parameters.
#' @exportS3Method 
print.tsacf <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    type <- x$type
    lag <- x$lag
    acf.res <- if (type == "correlation") {x$acf} else if (type == "partial") {x$pacf} else if (type == "covariance") {x$acov} else if (type == "cross-covariance") {x$ccov} else if (type == "cross-correlation") {x$ccf}
    x.clim <- x$clim
    ylabel <- switch(type, "correlation" = "ACF", "partial" = "PACF", "covariance" = "Autocovariance", "cross-covariance" = "Cross-covariance", "cross-correlation" = "CCF")
    cilevel <- paste0("CI(", (1 - x$alpha) * 100, "%)")
    out <- list(lag, acf.res)
    if (type == "correlation" | type == "partial" | type == "cross-correlation") {out <- c(out, list(-x.clim, x.clim))}
    out <- as.data.frame(out)
    colnames(out) <- c("Lags", ylabel, if (type == "correlation" | type == "partial" | type == "cross-correlation") {c(paste(cilevel, "Lower"), paste(cilevel, "Upper"))})
    print(out, digits = digits, row.names = FALSE, print.gap = 3)
}

##### Generate ACF/PACF/ACOV Plot with ggplot #####
#' @rdname tsacf
#' @param ... other printing or plotting parameters.
#' @exportS3Method 
plot.tsacf <- function(x, title = NULL, ...)
{
    type <- x$type
    lag <- x$lag
    acf.res <- if (type == "correlation") {x$acf} else if (type == "partial") {x$pacf} else if (type == "covariance") {x$acov} else if (type == "cross-covariance") {x$ccov} else if (type == "cross-correlation") {x$ccf}
    x.clim <- x$clim
    ylabel <- switch(type, "correlation" = "ACF", "partial" = "PACF", "covariance" = "Autocovariance", "cross-covariance" = "Cross-covariance", "cross-correlation" = "CCF")
    ymin <- min(acf.res)
    ymax <- max(acf.res)
    if(type == "partial")
    {
        lag <- c(0L, lag)
        acf.res <- c(0, acf.res)
        x.clim <- c(NA, x.clim)
        names(acf.res)[1L] <- "0"
        names(x.clim)[1L] <- "0"
    }
    if (type == "correlation" | type == "partial" | type == "cross-correlation")
    {
        if (sign(ymin) != sign(ymax))
        {
            ylim <- c(-1, 1)
            yclim <- list(-x.clim, x.clim)
        }
        else if (ymin < 0 & ymax <= 0)
        {
            ylim <- c(-1, 0)
            yclim <- list(-x.clim, rep(0, length(x.clim)))
        }
        else if (ymin >= 0 & ymax > 0)
        {
            ylim <- c(0, 1)
            yclim <- list(rep(0, length(x.clim)), x.clim)
        }
        ybreaks <- seq(ylim[1], ylim[2], 0.2)
        yscale <- ggplot2::scale_y_continuous(name = ylabel, limits = ylim, breaks = ybreaks, labels = scales::comma)
    }
    else if (type == "covariance" | type == "cross-covariance")
    {
        yextr <- max(abs(ymin), ymax)
        ybase <- if (max(yextr) > 1) {10 ^ floor(log10(max(yextr)))} else {round(yextr, 1)}
        ymaxlim <- if (ybase > 1) {floor(max(yextr) / ybase) * ybase} else {1}
        ylim <- if (sign(ymin) != sign(ymax)) {c(-ymaxlim, ymaxlim)} else if (ymin < 0  & ymax <= 0) {c(-ymaxlim, 0)} else if (ymin >= 0  & ymax > 0) {c(0, ymaxlim)}
        ybreaks <- seq(from = ylim[1], to = ylim[2], length.out = 11)
        yscale <- ggplot2::scale_y_continuous(name = ylabel, labels = scales::comma)
        x.clim <- NULL
    }
    clribbon <- if (type == "covariance" | type == "cross-covariance") {NULL} else {ggplot2::geom_ribbon(mapping = aes(x = lag, ymin = yclim[[1]], ymax = yclim[[2]]), colour = "royalblue", fill = "royalblue", alpha = 0.5, na.rm = TRUE)}
    zerovline <- if (type == "cross-correlation" | type == "cross-covariance") {ggplot2::geom_vline(aes(xintercept = 0), colour = "grey40", linewidth = 0.7, linetype = "dashed")} else {NULL}
    if (is.null(title)) {title <- paste(if (type == "correlation") {"ACF"} else if (type == "partial") {"PACF"} else if (type == "covariance") {"Autocovariance"} else if (type == "cross-covariance") {"Cross-covariance"} else if (type == "cross-correlation") {"CCF"}, "Plot of", x$x.name)}
    x.plot <- ggplot2::ggplot(mapping = aes(x = lag, y = acf.res)) +
        ggplot2::geom_bar(stat = "identity", fill = "darkgrey", width = 0.5) + clribbon + zerovline + ggplot2::xlab("Lags") + yscale +
        ggplot2::ggtitle(title) + ggplot2::theme(plot.title = element_text(hjust = 0.5)) + ggplot2::theme(axis.title.y = element_text(angle = 90, vjust = 2)) +
        ggplot2::theme(plot.margin = margin(t = 10, r = 10, b = 5, l = 10)) + ggplot2::geom_hline(aes(yintercept = 0), color = "grey40", linewidth = 0.7)
    print(x.plot)
}

##### Time Series Decomposition #####
#' Decompose a Time Series 
#' @description Decompose a time series into trend, cyclical, seasonal and irregular components. Deals with additive or multiplicative components.
#' @name tsdecomp
#' @rdname tsdecomp
#' @param x a time series for which decomposition is required or a `\code{tsdecomp}` object.
#' @param type type of the time series components. Available options are `\code{additive}` and `\code{multiplicative}`. Can be abbreviated. Default is `\code{additive}`.
#' @param trend.method estimating method of the trend. Available options are `\code{lm}` (linear) and `\code{loess}` (non-linear). Default is `\code{lm}`.
#' @param tcc.order moving average order for the estimation of the trend-cycle component.
#' @param x.name name of the series. If omitted here, the series name found by \code{tsname()} will be taken over here. If \code{tsname()} is \code{NULL}, the variable name will be used instead. Default is \code{NULL}.
#' @param show.plot logical. If \code{TRUE}, forecasting plot will be displayed directly. Default is \code{TRUE}.
#' @details The additive model used is:
#' \deqn{Y_t = T_t + C_t + S_t + I_t}    
#' The multiplicative model used is:
#' \deqn{Y_t = T_t \cdot C_t \cdot S_t \cdot I_t}    
#' The function first determines the trend-cycle component using a moving average, and removes it from the time series. Then, the seasonal figure is computed by averaging, for each time unit, over all periods. The seasonal figure is then centred. Finally, the error component is determined by removing trend and seasonal figure (recycled as needed) from the original time series.
#' \cr\cr This only works well if `\code{x}` covers an integer number of complete periods.
#' @returns An object of class "\code{tsdecomp}" with following components:
#' @return \item{x}{original series data}
#' @return \item{x.time}{list of time in which the series values were observed.} 
#' @return \item{x.timegap}{time gap between the series and forecasted values.}
#' @return \item{x.name}{name of the time series for which forecasts was requested.}
#' @return \item{trend}{value of the trend component for each observation.}
#' @return \item{cycle}{value of the cyclical component for each observation.} 
#' @return \item{trend.cycle}{trend-cyclical component value of each observation.}
#' @return \item{detrended}{value of each observation after removing the trend component.}
#' @return \item{seasonal}{value of the seasonal component for each observation.}
#' @return \item{seasonal.adjusted}{value of each observation after removing the seasonal component}
#' @return \item{random}{value of the irregular component for each observation.}
#' @return \item{seasonal.effect}{value expressing the estimated overall effect of each season in the time series.}
#' @return \item{type}{type of the time series decomposition.}
#' @author Ka Yui Karl Wu
#' @references Hyndman, R. J., & Athanasopoulos, G. (2021). Forecasting: Principles and practice (3rd ed.). OTexts. \cr \url{https://otexts.com/fpp3/} 
#' @references Kendall, M., & Stuart, A. (1983) The Advanced Theory of Statistics, Vol.3, Griffin. pp. 410-414.
#' @example inst/examples/tsdecomp/tsdecomp.R
#' @importFrom stats lm
#' @importFrom stats loess
#' @importFrom stats aggregate
#' @export
tsdecomp <- function(x, type = c("additive", "multiplicative"), trend.method = c("lm", "loess"), tcc.order = 3, x.name = NULL, show.plot = TRUE, ...)
{
    type <- match.arg(type)
    trend.method <- match.arg(trend.method)
    sc.exist <- if (tsfreq(x) == "day") {FALSE} else {ifelse(frequency(x) != 1, TRUE, FALSE)}
    trend.d <- data.frame(list(t = 1:length(x), x = x))
    if (trend.method == "lm")
    {
        trend.lm <- lm(x ~ t, data = trend.d)
    }
    else if (trend.method == "loess")
    {
        trend.lm <- loess(x ~ t, data = trend.d)
        
    }
    x.name <- if (is.null(x.name) & !is.null(attr(x, "series.name"))) {tsname(x)} else if (!is.null(x.name)) {x.name} else {NULL}
    xtime <- tstime(x)
    tc <- unname(trend.lm$fitted)
    tcc <- filter(x, filter = rep(1 / tcc.order, tcc.order), sides = 2L)
    if (tcc.order %% 2 == 0) {tcc <- filter(tcc, filter = rep(0.5, 2L), sides = 1L)}
    cc <- if (type == "additive") {tcc - tc} else {tcc / tc}
    sic <- if (type == "additive") {x - tcc} else {x / tcc}
    dtc <- if (type == "additive") {x - tc} else {x / tc}
    if (sc.exist)
    {
        t1 <- tstime(x)$time[1]
        wdlist <- list(Sunday = 0, Monday = 1, Tuesday = 2, Wednesday = 3, Thursday = 4, Friday = 5, Saturday = 6)
        ctype <- attr(x, "seasonal.cycle")
        callfunc <- if (ctype == "min") {"minute"} else if (ctype == "sec") {"second"} else if (ctype == "weekday") {"weekdays"} else {ctype}
        season.start <- do.call(callfunc, list(t1))
        if (ctype == "weekday") {season.start <- wdlist[season.start]}
        season.init <- list(month = 1, week = 1, quarter = 1, day = 1, weekdays = 0, hour = 0, min = 0, sec = 0)
        season.total <- list(month = 12, week = 52, quarter = 4, day = 365, weekdays = 6, hour = 23, min = 59, sec = 59)
        season.seq <- c(season.start:as.numeric(season.total[ctype]), if (season.start != as.numeric(season.init[ctype])) {as.numeric(season.init[ctype]):(season.start - 1)} else {c()})
        season.nam <- rep(season.seq, length.out = length(x))
        season.d <- data.frame(list(Season = season.nam, Effect = sic))
        sc.effect <- aggregate(Effect ~ Season, data = season.d, FUN = mean)
        sc.effect$Season <- season.seq
        sc.effect$Effect <- if (type == "additive") {sc.effect$Effect - mean(sc.effect$Effect)} else {sc.effect$Effect / mean(sc.effect$Effect)}
        sc <- sc.effect$Effect[match(season.d$Season, sc.effect$Season)]
        ic <- if (type == "additive") {sic - sc} else {sic / sc}
        sadj <- if (type == "additive") {x - sc} else {x / sc}
    }
    else
    {
        ic <- sic
    }
    x.name <- if (is.null(x.name) & is.null(tsname(x))) {deparse1(substitute(x))} else if (is.null(x.name) & !is.null(tsname(x))) {tsname(x)} else {x.name}
    numcomp <- c(list(trend = tc, cycle = cc, trend.cycle = tcc, detrended = dtc), if (sc.exist) {list(seasonal = sc, seasonal.adjusted = sadj)}, list(random = ic))
    for (j in names(numcomp)) {attributes(numcomp[[j]]) <- attributes(x)}
    out <- c(list(x = x, x.time = xtime$time, x.timegap = xtime$frequency, x.name = x.name), numcomp, if (sc.exist) {list(seasonal.effect = sc.effect)}, list(type = type))
    if (show.plot)
    {
        plot.tsdecomp(out, ...)
    }
    return(structure(out, class = "tsdecomp"))
}

##### Print Time Series Decomposition #####
#' @rdname tsdecomp
#' @param decomp.incl time series components that should be printed. Available options are `\code{all}` (default), `\code{tc}` (trend), `\code{tcc}` (trend-cycles), `\code{tcc}` (cycles), `\code{detrend}`, `\code{ic}` (irregular), `\code{scadj}` (seasonally adjusted), `\code{sc}` (seasonality), and `\code{sceffect}` (seasonal effects).
#' @exportS3Method 
print.tsdecomp <- function(x, decomp.incl = c("all", "tc", "tcc", "cc", "detrend", "ic", "scadj", "sc", "sceffect"), ...)
{
    if (any(c("all", "tc") %in% decomp.incl))
    {
        cat("Trend:\n")
        cat("======\n")
        print(x$trend)
    }
    if (any(c("all", "cc") %in% decomp.incl))
    {
        cat("\nCycles:\n")
        cat("=======\n")
        print(x$cycle)
    }
    if (any(c("all", "tcc") %in% decomp.incl))
    {
        cat("\nTrend-Cycles:\n")
        cat("=============\n")
        print(x$trend.cycle)
    }
    if (any(c("all", "sc") %in% decomp.incl))
    {
        if ("seasonal" %in% names(x))
        {
            cat("\nSeasonality:\n")
            cat("============\n")
            print(x$seasonal)
        }
    }
    if (any(c("all", "ic") %in% decomp.incl))
    {
        cat("\nIrregular:\n")
        cat("==========\n")
        print(x$random)
    }
    if (any(c("all", "detrend") %in% decomp.incl))
    {
        cat("\nDetrended:\n")
        cat("==========\n")
        print(x$detrended)
    }
    if (any(c("all", "scadj") %in% decomp.incl))
    {
        if ("seasonal.adjusted" %in% names(x)) 
        {
            cat("\nSeasonally Adjusted:\n")
            cat("====================\n")
            print(x$seasonal.adjusted)
        }
    }
    if (any(c("all", "sceffect") %in% decomp.incl))
    {
        if ("seasonal.effect" %in% names(x)) 
        {
            cat("\nSeasonal Effects:\n")
            cat("=================\n")
            print(x$seasonal.effect)
        }
    }
}

##### Generate Plots for Time Series Decomposition #####
#' @rdname tsdecomp
#' @param plot.incl time series components that should be plotted. Available options are `\code{all}` (default), `\code{tc}` (trend), `\code{tcc}` (trend-cycles), `\code{tcc}` (cycles), `\code{detrend}`, `\code{ic}` (irregular), `\code{scadj}` (seasonally adjusted), `\code{sc}` (seasonality), and `\code{sceffect}` (seasonal effects). Ignored if \code{show.plot = FALSE}.
#' @param ... parameter values that can affect the time series decomposition plots.
#' @details The function `\code{plot}` generates the following plots: 
#' \tabular{lcl}{\code{Trend} \tab \tab a time series line plot together with a trend line. \cr
#' \code{Trend-Cycles} \tab \tab a time series line plot together with the trend-cycle component (moving averages). \cr
#' \code{Cycles} \tab \tab a line plot of the trend-cycle component (moving averages). \cr
#' \code{Detrended} \tab \tab a line plot of the time series after removing the trend component. \cr
#' \code{Irregular} \tab \tab a line plot of the irregular component. \cr
#' \code{Seasonally Adjusted} \tab \tab a line plot of the time series after removing the seasonal component. \cr
#' \code{Seasonal} \tab \tab a line plot of the seasonal component. \cr
#' \code{Seasonal Effect} \tab \tab a line plot of the overall estimated effect of each season in the time series. \cr
#' }
#' @exportS3Method 
plot.tsdecomp <- function(x, plot.incl = c("all", "tc", "tcc", "cc", "detrend", "ic", "scadj", "sc", "sceffect"), ...)
{
    pwidth <- 0.7
    if (any(c("all", "tc") %in% plot.incl)) {tslineplot(x$x, pred = x$trend, title = "Trend", x.name = x$x.name, pred.name = "Trend", pred.lwidth = pwidth)}
    if (any(c("all", "tcc") %in% plot.incl)) {tslineplot(x$x, pred = x$trend.cycle, title = "Trend-Cycles", x.name = x$x.name, pred.name = "Trend & Cycles", pred.lwidth = pwidth)}
    if (any(c("all", "cc") %in% plot.incl)) {tslineplot(x$cycle, title = "Cycles", x.name = x$x.name, x.lwidth = pwidth, x.col = "steelblue4")}
    if (any(c("all", "detrend") %in% plot.incl)) {tslineplot(x$detrended, title = "Detrended", x.name = x$x.name, pred.name = "Series without trend", x.col = "steelblue4", x.lwidth = pwidth)}
    if (any(c("all", "ic") %in% plot.incl)) {tslineplot(x$random, title = "Irregular", x.name = x$x.name, pred.name = "Irregular Component", x.col = "steelblue4", x.lwidth = pwidth)}
    if ((any(c("all", "scadj") %in% plot.incl)) & ("seasonal.adjusted" %in% names(x))) {tslineplot(x$x, pred = x$seasonal.adjusted, title = "Seasonally Adjusted", x.name = x$x.name, pred.name = "Series without Seasonality", pred.lwidth = pwidth)}
    if ((any(c("all", "sc") %in% plot.incl)) & ("seasonal" %in% names(x))) {tslineplot(x$seasonal, title = "Seasonal", x.name = x$x.name, pred.name = "Seasonality", x.col = "steelblue4", x.lwidth = pwidth)}
    if ((any(c("all", "sceffect") %in% plot.incl)) & ("seasonal.effect" %in% names(x))) {tslineplot(x = x$seasonal.effect$Effect, t = x$seasonal.effect$Season, title = "Seasonal Effect", x.name = x$x.name, pred.name = "Seasonal Effects", t.name = "Seasons", x.lwidth = pwidth, x.col = "steelblue4", t.numbreak = 12, t.text.angle = 0)}
    if (all(plot.incl %in% c("all", "tc", "tcc", "cc", "detrend", "ic", "scadj", "sc", "sceffect")) == FALSE)
    {
        warning("Some of the `plot.incl` entries are inappropriate. Corresponding plots omitted.")
    }
}

##### Time Series Exploration #####
#' Explore a Time Series Numerically and Graphically
#' @description The function `\code{tsexplore}` generates statistics, various tests, and graphs regarding the location, deviation, distribution of a time series. 
#' @name tsexplore
#' @rdname tsexplore
#' @param x a time series to be explored or a `\code{tsexplore}` object.
#' @param show.plot logical. If \code{TRUE}, all exploration charts will be displayed directly. Default is \code{TRUE}.
#' @param x.name name of the series. If omitted here, the series name found by \code{tsname()} will be taken over here. If \code{tsname()} is \code{NULL}, the variable name will be used instead. Default is \code{NULL}.
#' @param mu test value specified under the null hypothesis of the t-test for the mean location. Default is \code{0}.
#' @param adf.lag number of AR lags included in the ADF test. Default is \code{0}.
#' @param lag.max maximum lag at which to calculate the acf. Default is \code{8}. Will be automatically limited to one less than the number of observations in the series. If the series has less than 8 observations.
#' @returns An object of class "\code{tsexplore}" with following components:
#' @return \item{x}{original series data}
#' @return \item{x.time}{list of time in which the series values were observed.} 
#' @return \item{x.timegap}{time gap between the series and forecasted values.}
#' @return \item{x.name}{name of the time series for which forecasts was requested.}
#' @return \item{stats}{a list of statistics and test results conducted on the time series}
#' @section Details of the `\code{stats}` component:
#' The following statistics and test results are stored in the component `\code{stats}` of the `\code{tsexplore}` object: 
#' \tabular{lcl}{\code{statistics} \tab \tab \code{n} (number of observations), \code{nvalid} (number of valid observations), \code{sum}, \code{mean}, \code{median}, \code{skewness}, \code{kurtosis}, \code{cv} (coefficient of variation) \cr
#' \code{variability} \tab \tab \code{variance}, \code{sd} (standard deviation), \code{range}, \code{iqr} (interquartile range) \cr
#' \code{quantiles} \tab \tab \code{minimum}, \code{q1} (1st quartile), \code{median}, \code{q3} (3rd quartile), \code{maximum} \cr
#' \code{autocorrelation} \tab \tab \code{acf} (autocorrelation function), \code{pacf} (partial autocorrelation function) - from lag 0 (ACF) or lag 1 (PACF) until \code{lag.max} \cr
#' \code{tests} \tab \tab \code{location} (t-test), \code{normality} (Shapiro-Wilk-test), \code{stationarity} (ADF-test), \code{independence} (Ljung-Box-test) - each test contains the test statistics (\code{statistics} and the p-value (\code{p.value})) \cr
#' }
#' @author Ka Yui Karl Wu
#' @references Hyndman, R. J., & Athanasopoulos, G. (2021). Forecasting: Principles and practice (3rd ed.). OTexts. \cr \url{https://otexts.com/fpp3/} 
#' @example inst/examples/tsexplore/tsexplore.R
#' @importFrom tseries adf.test
#' @importFrom stats quantile
#' @importFrom stats sd
#' @importFrom stats shapiro.test
#' @importFrom stats Box.test
#' @importFrom stats t.test
#' @importFrom stats median
#' @importFrom stats var
#' @export
tsexplore <- function(x, show.plot = TRUE, x.name = NULL, mu = 0, adf.lag = 0, lag.max = 8, ...)
{
    orgmatch <- match.call()
    arglist <- as.list(orgmatch)[-1]
    argnam <- names(arglist)
    x.name <- if (is.null(x.name) & is.null(tsname(x))) {deparse1(substitute(x))} else if (is.null(x.name) & !is.null(tsname(x))) {tsname(x)} else {x.name}
    xtime <- tstime(x)
    x.ttest <- suppressWarnings(suppressMessages(t.test(x, mu = mu)))
    xtemp <- if (length(x) > 5000) {tail(x, 5000)} else {x}
    x.sptest <- suppressWarnings(suppressMessages(shapiro.test(xtemp)))
    x.adftest <- suppressWarnings(suppressMessages(tseries::adf.test(na.omit(x), alternative = "stationary", k = adf.lag)))
    x.lbtest <- suppressWarnings(suppressMessages(Box.test(x, type = "Ljung-Box")))
    loctest <- list(statistic = x.ttest$statistic, p.value = x.ttest$p.value, null.value = mu)
    normtest <- list(statistic = x.sptest$statistic, p.value = x.sptest$p.value)
    stattest <- list(statistic = x.adftest$statistic, p.value = x.adftest$p.value)
    indtest <- list(statistic = x.lbtest$statistic, p.value = x.lbtest$p.value)
    x.test <- list(location = loctest, normality = normtest, stationarity = stattest, independence = indtest)
    tsmean <- mean(x, na.rm = TRUE)
    tsdev <- x - tsmean
    nreal <- length(na.omit(x))
    tsskew <- (sum(tsdev ^ 3L, na.rm = TRUE) / nreal) / ((sum(tsdev ^ 2L, na.rm = TRUE) / nreal) ^ 1.5)
    tskurt <- (sum(tsdev ^ 4L, na.rm = TRUE) / nreal) / ((sum(tsdev ^ 2L, na.rm = TRUE) / nreal) ^ 2L) - 3L
    x.stat <- list(n = length(x), nvalid = nreal, sum = sum(x, na.rm = TRUE), mean = tsmean, median = median(x, na.rm = TRUE), skewness = tsskew, kurtosis = tskurt, cv = sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
    x.var <- list(variance = var(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), range = diff(range(x, na.rm = TRUE)), iqr = IQR(x, na.rm = TRUE))
    x.quantiles <- list(minimum = min(x, na.rm = TRUE), q1 = quantile(x, prob = 0.25, na.rm = TRUE), median = median(x, na.rm = TRUE), q3 = quantile(x, prob = 0.75, na.rm = TRUE), maximum = max(x, na.rm = TRUE))
    x.acf <- list(acf = tsacf(x = x, lag.max = lag.max, show.plot = FALSE)$acf, pacf = tspacf(x = x, lag.max = lag.max, show.plot = FALSE)$pacf)
    xcomp.stat <- list(statistics = x.stat, variability = x.var, quantiles = x.quantiles, autocorrelation = x.acf, tests = x.test)
    out <- c(list(x = x, x.time = xtime$time, x.timegap = xtime$frequency, x.name = x.name), list(stats = xcomp.stat))
    if (show.plot)
    {
        plot_arg <- c("histbin", "lwidth", "pwidth", "x.col", "extra.col", if ("plot.incl" %in% argnam) {"plot.incl"})
        user_plot_arglist <- arglist[plot_arg[plot_arg %in% argnam]]
        plot_arglist <- c(list(x = out), user_plot_arglist)
        do.call(plot.tsexplore, args = plot_arglist)
    }
    return(structure(out, class = c("tsexplore")))
}

##### Print Time Series Exploration #####
#' @rdname tsexplore
#' @param trend a character string indicating whether and how the trend line should be fitted in the time series line plot. Available options are `\code{linear}`, `\code{smooth}`, `\code{none}`. Default is `\code{linear}`.
#' @param histbin a numeric value to specify the number of bins in the histogram. Can be omitted. Default is \code{15}.
#' @param lwidth line width of the series line plot. Default is \code{0.7}.
#' @param pwidth size of the markers in the QQ plot. Default is \code{0.7}.
#' @param x.col line colour of the time series line plot. Default is `\code{darkgrey}`.
#' @param extra.col colour of extra information in the plots. Default is `\code{red}`.
#' @param plot.incl time series components that should be plotted. Available options are `\code{all}` (default), `\code{line}` (line plot), `\code{hist}` (histogram), `\code{box}` (boxplot), `\code{qq}` (QQ plot), `\code{acf}` (ACF plot), and `\code{pacf}` (PACF plot). Ignored if \code{show.plot = FALSE}.
#' @param ... parameter values that can affect the plots created for time series exploration.
#' @exportS3Method 
plot.tsexplore <- function(x, trend = c("linear", "smooth", "none"), histbin = 15, 
                           lwidth = 0.7, pwidth = 0.7, x.col = "darkgrey", extra.col = "red", 
                           plot.incl = c("all", "line", "hist", "box", "qq", "acf", "pacf"), ...)
{
    trend <- match.arg(trend)
    if (any(c("all", "line") %in% plot.incl)) {tslineplot(x = x$x, t = x$x.time, x.name = x$x.name, trend = trend, x.lwidth = lwidth, pred.lwidth = lwidth, x.col = x.col, pred.col = extra.col)}
    if (any(c("all", "hist") %in% plot.incl)) {tshistogram(x = x$x, density = TRUE, x.name = x$x.name, density.lwidth = lwidth, x.col = x.col, density.col = extra.col, bins = histbin)}
    if (any(c("all", "box") %in% plot.incl)) {tsboxplot(x = x$x, x.name = x$x.name, x.col = x.col, mean.col = extra.col)}
    if (any(c("all", "qq") %in% plot.incl)) {tsqqplot(x = x$x, x.name = x$x.name, qq.lwidth = lwidth, qq.pwidth = pwidth, qq.col = x.col, qqline.col = extra.col)}
    lag.max <- length(x$stats$autocorrelation$pacf)
    if (any(c("all", "acf") %in% plot.incl)) {tsacf(x = x$x, x.name = x$x.name, lag.max = lag.max)}
    if (any(c("all", "pacf") %in% plot.incl)) {tspacf(x = x$x, x.name = x$x.name, lag.max = lag.max)}
    if (all(plot.incl %in% c("all", "line", "hist", "box", "qq", "acf", "pacf")) == FALSE)
    {
        warning("Some of the `plot.incl` entries are inappropriate. Corresponding plots omitted.")
    }
}

##### Print Time Series Exploration #####
#' @rdname tsexplore
#' @param digits the number of significant digits.
#' @param stats.incl time series statistics that should be printed. Available options are `\code{all}` (default), `\code{stats}` (statistics), `\code{var}` (variability), `\code{qtls}` (quantiles), `\code{acf}` (autocorrelation), and `\code{tests}` (tests).
#' @exportS3Method 
print.tsexplore <- function(x, digits = max(3L, getOption("digits") - 3L), stats.incl = c("all", "stats", "var", "qtls", "acf", "tests"), ...)
{
    scipenval <- getOption("scipen")
    options(scipen = 999)
    on.exit(options(scipen = scipenval))
    stats <- x$stats
    cat("\nSeries: ", x$x.name, "\n", sep = "")
    cat("\nNumber of Observations:       ", stats$statistics$n, "\n", sep = " ")
    cat("Number of Valid Observations: ", stats$statistics$nvalid, "\n", sep = " ")
    if (!is.null(stats))
    {
        if (any(c("all", "stats") %in% stats.incl)) 
        {
            cat("\nStatistics:\n", rep("=", 11L),"\n", sep = "")
            stat <- as.data.frame(stats$statistics[names(stats$statistics) != c("n", "nvalid")])
            colnames(stat) <- c("Sum", "Mean", "Median", "Skewness", "Kurtosis", "Coef. of Variation")
            print(stat, print.gap = 4L, digits = digits, row.names = FALSE)
        }
        if (any(c("all", "var") %in% stats.incl))
        {
            cat("\nVariability:\n", rep("=", 12L),"\n", sep = "")
            disp <- data.frame(stats$variability)
            colnames(disp) <- c("Variance", "Std. Deviation", "Range", "Interquartile Range")
            print(disp, row.names = FALSE, print.gap = 4L, digits = digits)
        }
        if (any(c("all", "qtls") %in% stats.incl))
        {
            cat("\nQuantiles:\n", rep("=", 10L),"\n", sep = "")
            qrtl <- data.frame(stats$quantiles)
            colnames(qrtl) <- c("Minimum", "1st Quartile", "Median", "3rd Quartile", "Maximum")
            print(qrtl, row.names = FALSE, print.gap = 4L, digits = digits)
        }
        if (any(c("all", "acf") %in% stats.incl))
        {
            cat("\nAutocorrelation:\n", rep("=", 16L),"\n", sep = "")
            acf <- data.frame(stats$autocorrelation$acf, c(NA, stats$autocorrelation$pacf))
            colnames(acf) <- c("ACF", "PACF")
            print(round(t(acf), digits), print.gap = 3L)
        }
        if (any(c("all", "tests") %in% stats.incl))
        {
            cat("\nTests:\n", rep("=", 6L),"\n", sep = "")
            ltest <- stats$tests$location
            ntest <- stats$tests$normality
            stest <- stats$tests$stationarity
            itest <- stats$tests$independence
            test <- data.frame(type = c(paste("t-Test ", "(mu=", stats$tests$location$null.value, ")", sep = ""), "Shapiro-Wilk", "Augmented Dickey-Fuller", "Ljung-Box"), statistic = c(ltest$statistic, ntest$statistic, stest$statistic, itest$statistic), p.value = round(c(ltest$p.value, ntest$p.value, stest$p.value, itest$p.value), digits))
            colnames(test) <- c("Test", "Statistic", "p-Value")
            rownames(test) <- c("Location", "Normality", "Stationarity", "Independence")
            print(test, print.gap = 4L, digits = digits, right = FALSE)
        }
        cat("\n")
    }
    if (all(stats.incl %in% c("all", "stats", "var", "qtls", "acf", "tests")) == FALSE)
    {
        warning("Some of the `stats.incl` entries are inappropriate. Corresponding outputs omitted.")
    }
}

##### Box Plot #####
#' Box Plots
#' @description Produce box-and-whisker plot of a given univariate time series.
#' @name tsboxplot
#' @rdname tsboxplot
#' @param x a univariate time series object or a numeric vector or matrix.
#' @param title title of the box plot
#' @param x.name name of the series. If omitted here, the series name found by \code{tsname()} will be taken over here. If \code{tsname()} is \code{NULL}, the variable name will be used instead. Default is \code{NULL}.
#' @param x.col line colour of the box plot.
#' @param mean.col colour of the dot indicating the series mean.
#' @author Ka Yui Karl Wu
#' @details Since \code{x} is a univariate time series, no parallel box plots can be plotted here.
#' @details Missing values are ignored when forming boxplots.
#' @example inst/examples/tsboxplot/tsboxplot.R
#' @return A boxplot of \code{x} will be displayed with no further values or objects returned.
#' @references Chambers, J. M., Cleveland, W. S., Kleiner, B., & Tukey, P. A. (1983). Graphical Methods for Data Analysis. Wadsworth & Brooks/Cole.
#' @import scales
#' @importFrom scales comma
#' @export
tsboxplot <- function(x, title = NULL, x.name = NULL, x.col = "darkgrey", mean.col = "steelblue4")
{
    x.name <- if (is.null(x.name) & is.null(tsname(x))) {deparse1(substitute(x))} else if (is.null(x.name) & !is.null(tsname(x))) {tsname(x)} else {x.name}
    if (is.null(title)) {title <- paste("Box Plot of", x.name)}
    xplot <- ggplot2::ggplot(mapping = aes(y = x)) + ggplot2::geom_boxplot(outlier.color = "red", colour = x.col, na.rm = TRUE) +
        ggplot2::geom_point(mapping = aes(y = mean(x), x = 0L), color = mean.col) +
        ggplot2::scale_y_continuous(name = x.name, labels = scales::comma) +
        ggplot2::scale_x_continuous(name = "", breaks = NULL) +
        ggplot2::ggtitle(title) + ggplot2::theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) +
        ggplot2::theme(axis.title.x = element_text(vjust = -1L)) +
        ggplot2::theme(axis.title.y = element_text(angle = 90, vjust = 3L))
    suppressWarnings(print(xplot))
}

##### Time Series Line Plot #####
#' Time Series Line Plots
#' @description Produce line plot of a given univariate time series.
#' @name tslineplot
#' @rdname tslineplot
#' @param x a univariate time series object or a numeric vector or matrix.
#' @param t a vector or list of time periods in which the series values were observed. The length of \code{t} must be identical to the length of \code{x}.
#' @param pred a vector of predicted or forecasted values of a univariate time series. This parameter can be omitted. Default is \code{NULL}.
#' @param pred.t a vector of time periods in which the predicted or forecasted values of a univariate time series are estimated. This parameter can be omitted. Default is \code{NULL}. 
#' @param cil a vector of the prediction intervals' lower limits. Only necessary if \code{pred} and \code{pred.t} are provided. This parameter can be omitted. Default is \code{NULL}.
#' @param ciu a vector of the prediction intervals' upper limits. Only necessary if \code{pred} and \code{pred.t} are provided. This parameter can be omitted. Default is \code{NULL}.
#' @param ci.t a vector of the time periods in which the prediction intervals are estimated.This parameter can be omitted. Default is \code{NULL}.
#' @param trend indicate whether a trend line should be included in the time series plot. Available options are `\code{none}`, `\code{linear}`, and `\code{smooth}`. If `\code{linear}`, a straight trend line estimated by linear regression model will be included. If `\code{smooth}`, the trend line will be estimated by LOESS regression model. Default is \code{NULL}, indicating that no trend line should be displayed.
#' @param title title of the time series line plot.
#' @param x.name name of the series. If omitted here, the series name found by \code{tsname()} will be taken over here. If \code{tsname()} is \code{NULL}, the variable name will be used instead. Default is \code{NULL}.
#' @param pred.name name of the series' predicted/forecasted values. Only necessary if \code{pred} and \code{pred.t} are provided. Default is `\code{Predicted}`.
#' @param x.lwidth line width of the series line plot. Default is \code{0.7}.
#' @param pred.lwidth line width of the line plot for the predicted/forecasted values. Default is \code{0.7}.
#' @param x.col line colour of the time series line plot. Default is `\code{darkgrey}`.
#' @param pred.col line colour of the line plot for the predicted/forecasted values. Default is `\code{steelblue4}`.
#' @param ci.col area colour of the prediction intervals. Default is `\code{royalblue}`.
#' @param t.name name of the x-axis (time axis). Default is `\code{Date/Time}`.
#' @param t.text.angle angle of the tick labels on the x-axis (time axis). Default is \code{90} (vertical).
#' @param t.numbreak number of tick labels on the x-axis (time axis). Default is \code{10}.
#' @param ylim value limit of the y-axis. The values should be specified in \code{c(lower_limit, upper_limit)}, where \code{lower_limit} and \code{upper_limit} are the values of the smallest and largest number of the y-axis, respectively. Default is \code{NULL}.
#' @author Ka Yui Karl Wu
#' @details If \code{x} is a \code{ts} object, parameter \code{t} can be omitted. You can convert a vector, a matrix or data frame column using \code{tsconvert} to a \code{ts} object.
#' @example inst/examples/tslineplot/tslineplot.R
#' @return A line plot of \code{x} will be displayed with no further values or objects returned.
#' @importFrom utils install.packages
#' @import ggplot2
#' @export
tslineplot <- function(x, t = NULL, pred = NULL, pred.t = NULL, cil = NULL, ciu = NULL, ci.t = NULL, trend = c("none", "linear", "smooth"),
                       title = NULL, x.name = NULL, pred.name = "Predicted", x.lwidth = 0.7, pred.lwidth = 0.7,
                       x.col = "darkgrey", pred.col = "steelblue4", ci.col = "royalblue",
                       t.name = "Date/Time", t.text.angle = 90, t.numbreak = 10, ylim = NULL)
{
    x.name <- if (is.null(x.name) & is.null(tsname(x))) {deparse1(substitute(x))} else if (is.null(x.name) & !is.null(tsname(x))) {tsname(x)} else {x.name}
    trend <- match.arg(trend)
    if (is.null(t))
    {
        xtime <- tstime(x)
        tgap <- xtime$frequency
        t <- tstimeformat(xtime$time, timegap = tgap)
    }
    else
    {
        if (length(x) != length(t))
        {
            stop("The time variable must have the same length as and the series!")
        }
        tgap <- if (is.null(tsfreq(x))) {tstimegap(t)} else {tsfreq(x)}
        t <- tstimeformat(t, timegap = tgap)
    }
    if (is.null(pred))
    {
        orig.s <- ggplot2::geom_line(mapping = aes(x = t, y = x), color = x.col, linewidth = x.lwidth, na.rm = TRUE)
        add.s <- NULL
        add.legend <- NULL
        add.legend.pos <- NULL
    }
    else
    {
        pred.t <- if (is.null(pred.t)) {head(t, length(pred))} else {tgap <- if (is.null(tsfreq(x))) {tstimegap(t)} else {tsfreq(x)}; tstimeformat(pred.t, timegap = tgap)}
        orig.s <- ggplot2::geom_line(mapping = aes(x = t, y = x, color = "observed"), na.rm = TRUE, linewidth = x.lwidth)
        add.s <- ggplot2::geom_line(aes(x = pred.t, y = pred, colour = "predicted"), na.rm = TRUE, linewidth = pred.lwidth)
        add.legend <- ggplot2::scale_color_manual(values = c(observed = x.col, predicted = pred.col), labels = c(observed = x.name, predicted = pred.name))
        add.legend.pos <- ggplot2::theme(legend.position = "bottom")
    }
    if (!is.null(cil) & !is.null(ciu) & !is.null(ci.t))
    {
        citgap <- tstimegap(ci.t)
        ci.t <- tstimeformat(ci.t, timegap = citgap)
        add.ci <- ggplot2::geom_ribbon(alpha = 0.5, mapping = aes(x = ci.t, ymin = as.numeric(cil), ymax = as.numeric(ciu)), colour = ci.col, fill = ci.col, na.rm = TRUE)
    }
    else
    {
        add.ci <- NULL
    }
    add.t <- if (trend == "none") {NULL} else {ggplot2::geom_smooth(mapping = aes(x = t, y = x), na.rm = TRUE, formula = y ~ x, method = if (trend == "linear") {"lm"} else {"loess"}, se = FALSE, linewidth = 0.7, color = pred.col)}
    xbreak <- c("year" = 5, "quarter" = 4, "month" = 12, "week" = 52, "weekdays" = 7, "day" = 30.25, "hour" = 4, "min" = 5, "sec" = 5)
    xtick <- if (length(t) < length(pred.t)) {pred.t} else {t}
    xlen <- max(length(t), length(pred.t))
    xseq <- optimTimeBreaks(len = xlen, frequency = xbreak[tgap], numbreaks = t.numbreak)
    if (is.null(title)) {title <- paste("Series Line Plot of", x.name)}
    if (!is.null(ylim))
    {
        draw0 <- if (sign(ylim[1]) != sign(ylim[2])) {TRUE} else {FALSE}
    }
    else
    {
        ymax <- max(x, pred, ciu, cil, na.rm = TRUE)
        ymin <- min(x, pred, ciu, cil, na.rm = TRUE)
        draw0 <- if (sign(ymax) != sign(ymin)) {TRUE} else {FALSE}
    }
    zerohline <- if (draw0) {ggplot2::geom_hline(aes(yintercept = 0), color = "grey50", linetype = "dashed", linewidth = 0.7)} else {NULL}
    xplot <- ggplot2::ggplot(mapping = aes(group = 1L)) +
        ggplot2::ggtitle(title) + orig.s + add.s + add.ci + add.legend + add.t + ggplot2::labs(color = NULL) +
        ggplot2::scale_y_continuous(name = x.name, labels = scales::comma, limits = ylim) +
        ggplot2::scale_x_discrete(name = t.name, breaks = xtick[xseq]) + add.legend.pos +
        ggplot2::theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) +
        ggplot2::theme(axis.title.x = element_text(vjust = -1L), axis.text.x = element_text(angle = t.text.angle)) +
        ggplot2::theme(axis.title.y = element_text(angle = 90, vjust = 3L)) + zerohline
    suppressWarnings(print(xplot))
}

##### Histogram #####
#' Histograms
#' @description Produce a histogram of a given univariate time series.
#' @name tshistogram
#' @rdname tshistogram
#' @param x a univariate time series object or a numeric vector or matrix.
#' @param bins a numeric value to specify the number of bins in the histogram. Can be omitted. Default is \code{NULL}.
#' @param density logical. Indicate whether the density curve of the normal distribution should be included. Default is \code{FALSE}.
#' @param density.lwidth line width of the density curve in the output plot. Will be ignored if \code{density = FALSE}. Default is 0.7.
#' @param title title of the histogram. Default is \code{NULL}.
#' @param x.name name of the series. If omitted here, the series name found by \code{tsname()} will be taken over here. If \code{tsname()} is \code{NULL}, the variable name will be used instead. Default is \code{NULL}.
#' @param x.col colour of the histogram bars. Default is `\code{darkgrey}`.
#' @param density.col colour of the density curve. Will be ignored if \code{density = FALSE}. Default is `\code{steelblue4}`.
#' @author Ka Yui Karl Wu
#' @example inst/examples/tshistogram/tshistogram.R
#' @return A histogram of \code{x} will be displayed with no further values or objects returned.
#' @references Venables, W. N., & Ripley. B. D. (2002) Modern Applied Statistics with S. Springer.
#' @importFrom stats dnorm
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @importFrom stats IQR
#' @export
tshistogram <- function(x, bins = NULL, density = FALSE, density.lwidth = 0.7, title = NULL, x.name = NULL, x.col = "darkgrey", density.col = "steelblue4")
{
    x.name <- if (is.null(x.name) & is.null(tsname(x))) {deparse1(substitute(x))} else if (is.null(x.name) & !is.null(tsname(x))) {tsname(x)} else {x.name}
    x <- as.data.frame(x)
    colnames(x) <- x.name
    add.d <- if (density) {ggplot2::stat_function(fun = dnorm, args = list(mean = mean(x[[x.name]], na.rm = TRUE), sd = sd(x[[x.name]], na.rm = TRUE)), linewidth = density.lwidth, color = density.col)} else {NULL}
    if (is.null(title)) {title <- paste("Histogram of", x.name)}
    if (is.null(bins)) {bins <- round((diff(range(x[[x.name]]))) / (2L * IQR(x[[x.name]], na.rm = TRUE) / (length(x[[x.name]]) ^ (1L / 3L))), 0L)}
    xplot <- ggplot2::ggplot(data = x, mapping = aes(x = get(x.name))) + ggplot2::geom_histogram(aes(y = after_stat(density)), na.rm = TRUE, colour = "white", fill = x.col, bins = bins) +
        add.d + ggplot2::scale_x_continuous(name = x.name, labels = scales::comma) + ggplot2::scale_y_continuous(name = "Density", labels = scales::comma) +
        ggplot2::ggtitle(title) + ggplot2::theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) +
        ggplot2::theme(axis.title.x = element_text(vjust = -1L)) +
        ggplot2::theme(axis.title.y = element_text(angle = 90, vjust = 3L))
    suppressWarnings(print(xplot))
}

##### QQ-Plot #####
#' Quantile-Quantile Plots
#' @description `\code{tsqqplot}` is a function to produce a normal QQ plot of the values in \code{y}.
#' @name tsqqplot
#' @rdname tsqqplot
#' @param x a univariate time series object or a numeric vector or matrix.
#' @param title title of the QQ plot. Default is \code{NULL}.
#' @param x.name name of the series. If omitted here, the series name found by \code{tsname()} will be taken over here. If \code{tsname()} is \code{NULL}, the variable name will be used instead. Default is \code{NULL}.
#' @param qq.pwidth size of the markers in the QQ plot. Default is \code{0.7}.
#' @param qq.lwidth line width of the theoretical normal line in the QQ plot. Default is \code{0.7}.
#' @param qq.col colour of the data points in the QQ plot. Default is `\code{black}`.
#' @param qqline.col colour of the theoretical normal line in the QQ plot. Default is `\code{red}`.
#' @author Ka Yui Karl Wu
#' @example inst/examples/tsqqplot/tsqqplot.R
#' @return A QQ plot of \code{x} will be displayed with no further values or objects returned.
#' @references Switzer, P. (1976). Confidence procedures for two-sample problems. \emph{Biometrika}, \strong{63}(1), 13-25. \doi{10.1093/biomet/63.1.13}.
#' @export
tsqqplot <- function(x, title = NULL, x.name = NULL, qq.pwidth = 0.7, qq.lwidth = 0.7, qq.col = "black", qqline.col = "red")
{
    x.name <- if (is.null(x.name) & is.null(tsname(x))) {deparse1(substitute(x))} else if (is.null(x.name) & !is.null(tsname(x))) {tsname(x)} else {x.name}
    x <- as.data.frame(as.numeric(x))
    colnames(x) <- x.name
    if (is.null(title)) {title <- paste("QQ-Plot of", x.name)}
    xplot <- ggplot2::ggplot(data = x, aes(sample = get(x.name))) + ggplot2::stat_qq_line(colour = qqline.col, linewidth = qq.lwidth, na.rm = TRUE) + stat_qq(size = qq.pwidth, colour = qq.col, na.rm = TRUE) +
        ggplot2::scale_x_continuous(name = "Theoretical", labels = scales::comma) +
        ggplot2::scale_y_continuous(name = "Empirical", labels = scales::comma) +
        ggplot2::ggtitle(title) + ggplot2::theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) +
        ggplot2::theme(axis.title.x = element_text(vjust = -1L)) +
        ggplot2::theme(axis.title.y = element_text(angle = 90, vjust = 3L))
    suppressWarnings(print(xplot))
}

##### Scatter Plot #####
#' Scatter Plot
#' @description Produce a scatter plot of two given univariate time series.
#' @name tsscatterplot
#' @rdname tsscatterplot
#' @param x,y two univariate time series object or a numeric vector or matrix.
#' @param reg optional. A logical value indicating whether a trend line estimated by regression should be included in the scatter plot. Default is \code{FALSE}.
#' @param title title of the histogram. Default is \code{NULL}.
#' @param x.name name of the series `\code{x}`. If omitted here, the series name found by \code{tsname(x)} will be taken over here. If \code{tsname()} is \code{NULL}, the variable name will be used instead. Default is \code{NULL}.
#' @param y.name name of the series `\code{y}`. If omitted here, the series name found by \code{tsname(y)} will be taken over here. If \code{tsname()} is \code{NULL}, the variable name will be used instead. Default is \code{NULL}.
#' @param pwidth size of the markers in the scatter plot. Default is \code{1}.
#' @param pcol colour of the data points in the scatter plot. Default is `\code{steelblue4}`.
#' @param regwidth width of the trend line in the scatter plot. Default is \code{0.7}.
#' @param regcol colour of the trend line in the scatter plot. Default is `\code{red}`.
#' @author Ka Yui Karl Wu
#' @example inst/examples/tsscatterplot/tsscatterplot.R
#' @return A scatter plot of \code{x} (values on the x-axis) and \code{y} (values on the y-axis) will be displayed with no further values or objects returned.
#' @export
tsscatterplot <- function(x, y, reg = FALSE, title = NULL, x.name = NULL, y.name = NULL, pwidth = 1, pcol = "steelblue4", regwidth = 0.7, regcol = "red")
{
    x.name <- if (is.null(x.name) & is.null(tsname(x))) {deparse1(substitute(x))} else if (is.null(x.name) & !is.null(tsname(x))) {tsname(x)} else {x.name}
    y.name <- if (is.null(y.name) & is.null(tsname(y))) {deparse1(substitute(y))} else if (is.null(y.name) & !is.null(tsname(y))) {tsname(y)} else {y.name}
    xy <- as.data.frame(list(x = x, y = y))
    colnames(xy) <- c(x.name, y.name)
    if (is.null(title)) {title <- paste("Scatter Plot of", x.name, "and", y.name)}
    rline <- if (reg) {ggplot2::geom_smooth(mapping = aes(x = x, y = y), na.rm = TRUE, formula = y ~ x, method = "loess", se = FALSE, color = regcol, linewidth = regwidth)} else {NULL}
    xplot <- ggplot2::ggplot(data = xy, aes(x = x, y = y)) + ggplot2::geom_point(color = pcol, size = pwidth) + rline +
        ggplot2::scale_x_continuous(name = x.name, labels = scales::comma) +
        ggplot2::scale_y_continuous(name = y.name, labels = scales::comma) +
        ggplot2::ggtitle(title) + 
        ggplot2::theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) +
        ggplot2::theme(axis.title.x = element_text(vjust = -1L)) +
        ggplot2::theme(axis.title.y = element_text(angle = 90, vjust = 3L))
    suppressWarnings(print(xplot))
}

##### Test for ARCH Effect #####
#' McLeod-Li Test for ARCH Effect
#' @description The function `\code{tsmltest}` applies the McLeod-Li test to examine whether ARCH effect exists in the squared residuals of an ARIMA model.
#' @name tsmltest
#' @rdname tsmltest
#' @param object a univariate time series object or a numeric vector or matrix.
#' @param lag.max maximum lag at which to examine the ARCH effect. Default is \code{NULL}. Will be automatically calculated using the formula \eqn{10\cdot \log_{10}(n)}, where \eqn{n} is the series length, if the parameter is omitted.
#' @author Ka Yui Karl Wu
#' @example inst/examples/tsmltest/tsmltest.R
#' @references McLeod, A. I., & Li, W. K. (1983). Diagnostic Checking ARMA Time Series Models Using Squared-Residual Autocorrelations. \emph{Journal of Time Series Analysis}, \strong{4}(4), 269-273. \cr \doi{10.1111/j.1467-9892.1983.tb00373.x}.
#' @return A list with two elements:
#' @return \item{\code{test.value}}{chi-square test statistics of the McLeod-Li test for the lags 1 to \code{lag.max}.}
#' @return \item{\code{p.value}}{corresponding p-values of the McLeod-Li test for the lags 1 to \code{lag.max}.}
#' @importFrom stats pchisq
#' @export
tsmltest <- function(object, lag.max = NULL)
{
    if (is.null(lag.max)) {lag.max <- 10 * log10(length(object$x.used))}
    teststat = rep(NA, lag.max)
    pval = rep(NA, lag.max)
    x <- object$residuals
    for (j in 1:lag.max) 
    {
        cor <- tsacf(x ^ 2, lag.max = j, show.plot = FALSE)
        n <- sum(!is.na(x))
        obs <- cor$acf[1:j]
        teststat[j] <- n * (n + 2) * sum(1/seq.int(n - 1, n - j) * obs ^ 2)
        pval[j] <- 1 - pchisq(teststat[j], j)
    }
    out <- list(test.value = teststat, p.value = pval)
    return(out)
}

##### Convert Data to Time Series #####
#' Convert One-Dimensional Data to Time Series
#' @description The function `\code{tsconvert}` is used to convert any one-dimensional vector/list into a time-series objects.
#' @name tsconvert
#' @rdname tsconvert
#' @param x a univariate time series object or a numeric vector or matrix.
#' @param t a vector or list of time in which the series values were observed.
#' @param frequency the number of observations per unit of time. If omitted, R will identify the frequency based on the time vector specified in \code{t}. Default is \code{NULL}.
#' @param format time format specified for the variable provided in \code{t}. If omitted, R will identify this automatically. Default is \code{NULL}.
#' @param x.name name of the series. If omitted here, the series name found by \code{tsname()} will be taken over here. If \code{tsname()} is \code{NULL}, the variable name will be used instead. Default is \code{NULL}.
#' @author Ka Yui Karl Wu
#' @details The function \code{tsconvert} is used to convert vectors/lists into time-series objects. These are vectors or matrices which inherit from class "ts" (and have additional attributes) and represent data sampled at equispaced points in time. Time series must have at least one observation, and although they need not be numeric there is very limited support for non-numeric series.
#' @details Argument \code{frequency} indicates the sampling frequency of the time series, with the default value 1 indicating one sample in each unit time interval. For example, one could use a value of 7 for frequency when the data are sampled daily, and the natural time period is a week, or 12 when the data are sampled monthly and the natural time period is a year. Values of 4 and 12 are assumed in (e.g.) print methods to imply a quarterly and monthly series respectively. Note that \code{frequency} does not need to be a whole number: for example, \code{frequency = 0.2} would imply sampling once every five time units. 
#' @example inst/examples/tsconvert/tsconvert.R
#' @return \code{x} is returned as a `\code{ts}` object with the attributes `\code{tsp}` (start, end, and frequency of \code{x}), `\code{series.name}` (series name, \emph{optional}), and `\code{seasonal.cycle}` (time gap between series observations).
#' @seealso \link{attributes}
#' @importFrom lubridate decimal_date
#' @export
tsconvert <- function(x, t, frequency = NULL, format = NULL, x.name = NULL)
{
    tsymbol <- c(":", "-", "/")
    tformat <- lapply(X = tsymbol, FUN = grepl, x = t)
    if (sum(tformat[[2]], tformat[[3]]) == 0 & sum(tformat[[1]]) > 0)
    {
        if (is.null(format)) format <- "%H:%M:%OS"
        timeformat <- "DateTime"
    }
    else if (sum(unlist(tformat)) == 0)
    {
        t <- as.character(t)
        if (is.null(format)) format <- "%Y"
        timeformat <- "OnlyYear"
    }
    else
    {
        if (is.null(format)) format <- c("%d/%m/%Y", "%Y-%m-%d", "%Y/%m/%d")
        timeformat <- "Normal"
    }
    t <- as.POSIXct(t, format = format)
    tfreq <- tstimegap(t)
    if (is.null(frequency))
    {
        timefreq <- c("year" = 1, "quarter" = 4, "month" = 12, "week" = 52.1429, "weekdays" = 7, "day" = 365.2422, "hour" = 24, "min" = 1 / 60, "sec" = 1 / 3600)
        frequency <- timefreq[tfreq]
    }
    if (timeformat == "Normal")
    {
        tsstart <- lubridate::decimal_date(t[1])
    }
    else if (timeformat == "OnlyYear")
    {
        tsstart <- lubridate::year(t[1])
    }
    else if (timeformat == "DateTime")
    {
        tsstart <- as.numeric(t[1])
    }
    x <- ts(x, start = tsstart, frequency = frequency)
    if (!is.null(x.name)) {x <- tsname(x, x.name)} else if (!is.null(attr(x, "series.name"))) {x <- tsname(x, attr(x, "series.name"))}
    attr(x, "seasonal.cycle") <- tfreq
    return(x)
}

##### Generate Differencing #####
#' Difference a Time Series
#' @description The function `\code{tsdiff}` generates the differenced series of a non-stationary time series.
#' @name tsdiff
#' @rdname tsdiff
#' @param x a univariate time series object or a numeric vector or matrix.
#' @param lag number of lags for non-seasonal differencing. Default is \code{1}.
#' @param order order of non-seasonal differencing. Default is \code{1}.
#' @param lag.D number of lags for seasonal differencing. Default is \code{0}.
#' @param order.D order of seasonal differencing. Default is \code{0}.
#' @author Ka Yui Karl Wu
#' @details The parameters \code{lag} and \code{lag.D} are only necessary if the lag difference for differencing is not 1 or \eqn{\ell}, the seasonal cycle length, respectively. If \code{order.D > 0} but \code{lag.D} is omitted, R will use the frequency of the series for this parameter.
#' @example inst/examples/tsdiff/tsdiff.R
#' @seealso \link{diff}
#' @return The differences of \code{x} will be returned with all the attributes being carried over. Unlike the function \code{\link{diff}}, the output series by \code{tsdiff} has the same length as the original series by adding \code{NA} to those observations at the beginning of the series for which no differencing can be carried out.
#' @export
tsdiff <- function(x, lag = 1L, order = 1L, lag.D = 0L, order.D = 0L)
{
    x.d <- x
    missnum <- 0L
    if (missing(order) & missing(lag) & (!missing(order.D) | !missing(lag.D))) {lag <- order <- 0L}
    if (order.D > 0L & lag.D == 0L) {lag.D <- frequency(x)}
    if (order.D == 0L & lag.D > 0L) {order.D <- 1L}
    if (order > 0L & lag == 0L) {lag <- 1L}
    if (order == 0L & lag > 0L) {order <- 1L}
    orders <- list(d = order, D = order.D)
    lags <- list(d = lag, D = lag.D)
    for (j in 1L:length(orders))
    {
        diff.o <- if (is.na(orders[j])) {0L} else {as.integer(orders[j])}
        lagnum <- 0
        if (diff.o > 0L)
        {
            lagnum <- if (is.na(lags[j])) {0L} else {as.integer(lags[j])}
            x.d <- diff(x.d, lag = lagnum, differences = diff.o)
        }
        missnum <- missnum + diff.o * lagnum
    }
    x.d <- c(rep(NA, missnum), x.d)
    attributes(x.d) <- attributes(x)
    return(x.d)
}

##### Generate Time Lags #####
#' Lag a Time Series
#' @description The function `\code{tslag}` back- or foreshifts a time series and generates the lagged version of it.
#' @name tslag
#' @rdname tslag
#' @param x a univariate time series object.
#' @param lag number of lags (in units of observations). Default is \code{1}.
#' @author Ka Yui Karl Wu
#' @return The same time series object as \code{x} after being back- or foreshifted for the number of periods specified in \code{lag}.
#' @details Note the sign of \code{lag}: a series lagged by a positive \code{lag} starts earlier.
#' @seealso \link{lag}
#' @example inst/examples/tslag/tslag.R
#' @export
tslag <- function(x, lag = 1L)
{
    if (!is.ts(x)) {stop("x must be a time series.")}
    xlen <- length(x)
    cutlen <- xlen - abs(lag)
    xlag <- if (lag > 0) {c(rep(NA, abs(lag)), head(x, cutlen))} else if (lag < 0) {xlag <- c(tail(x, cutlen), rep(NA, abs(lag)))} else {x}
    attributes(xlag) <- attributes(x)
    return(xlag)
}

##### Get the Frequency of a Time Series #####
#' Extract Information of a Time Series
#' @description The function `\code{tsfreq}` extract the frequency of a time series, while `\code{tstime}` extract the time periods in which the series data were observed, and `\code{tstimegap}` returns the time gap between the observations of a time series.
#' @name ts-functions
#' @rdname ts-functions
#' @param x a univariate time series object or a numeric vector or matrix.
#' @author Ka Yui Karl Wu
#' @example inst/examples/tsfreq/tsfreq.R
#' @return For `\code{tstime}`, a list will be returned to the user with two elements: \code{time} (observation time) and \code{frequnecy} (observation frequency).
#' @return For `\code{tsfreq}`, R extracts the attribute `\code{seasonal.cycle}` from the time series object \code{x}.
#' @return For `\code{tstimegap}`, R calculates the time gap between the time periods stored in the vector \code{t}.
#' @return So, if \code{x} and \code{t} are consistent and refer to the data and time of the same time series, the results of `\code{tsfreq}` and `\code{tstimegap}` as well as the \code{frequnecy} element of `\code{tstime}` must be identical.
#' @usage ## Extract frequency of a time series
#' tsfreq(x)
#' @export
tsfreq <- function(x)
{
    return(attr(x, "seasonal.cycle"))
}

##### Extract Time Series Date #####
#' @rdname ts-functions
#' @usage ## Extract observation time periods of a time series
#' tstime(x)
#' @example inst/examples/tstime/tstime.R
#' @import lubridate
#' @importFrom lubridate date_decimal
#' @importFrom lubridate days
#' @importFrom stats time
#' @importFrom stats is.ts
#' @export
tstime <- function(x)
{
    n <- length(x)
    if (is.ts(x))
    {
        xtime <- time(x)
        xfreq <- tsfreq(x)
        if (all(log10(xtime) < 4))
        {
            if (xfreq == "week")
            {
                sdays <- 7 * (0:(n - 1))
                sstart <- lubridate::date_decimal(xtime[1] + 0.00000001)
                xdate <- sstart + lubridate::days(sdays)
            }
            else
            {
                xtimeadj <- xtime + 0.0015
                xyear <- floor(xtimeadj)
                multifreq <- if (xfreq == "day") {ifelse(xyear %% 400 == 0L | xyear %% 4L == 0L, 366, 365)} else {rep(frequency(x), n)}
                xrest <- floor((xtimeadj - xyear) * multifreq) + 1L
                if (xfreq == "quarter") xrest <- (xrest - 1L) * 3L + 1L
                xdate <- paste0(xyear, "-", if (xfreq == "year") {rep("01", n)} else {sprintf("%02g", xrest)}, if ((xfreq != "day")) {rep("-01", n)})
            }
            tz <- "UTC"
            xformat <- paste0("%Y", if (xfreq == "day") {"-%j"} else {"-%m-%d"})
        }
        else
        {
            xdate <- as.integer(xtime)
            tz <- ""
            xformat = NULL
        }
        xposixt <- as.POSIXct(xdate, tz = tz, format = xformat)
        xfreq <- if (!is.null(tsfreq(x))) {tsfreq(x)} else {suppressMessages(tstimegap(xposixt))}
        out <- list(time = xposixt, frequency = xfreq)
    }
    else
    {
        warning("Time series is required to extract the observation dates of the data!")
        out <- NULL
    }
    return(out)
}

##### Extract the frequency of a date vector #####
#' @rdname ts-functions
#' @param t a vector or list of time in which the series values were observed.
#' @usage ## Extract time gaps between a time series' observations
#' tstimegap(t)
#' @example inst/examples/tstimegap/tstimegap.R
#' @importFrom lubridate year
#' @importFrom lubridate quarter
#' @importFrom lubridate month
#' @importFrom lubridate week
#' @importFrom lubridate day
#' @importFrom lubridate hour
#' @importFrom lubridate minute
#' @importFrom lubridate second
#' @importFrom lubridate is.POSIXt
#' @importFrom stats is.ts
#' @export
tstimegap <- function(t)
{
    if (!lubridate::is.POSIXt(t)) {t <- as.POSIXct(t, tryFormats = c("%Y", "%d/%m/%Y", "%H:%M:%OS", "%Y-%m-%d %H:%M:%OS", "%Y/%m/%d %H:%M:%OS", "%Y-%m-%d %H:%M", "%Y/%m/%d %H:%M", "%Y-%m-%d", "%Y/%m/%d"))}
    tcomp <- list(year = lubridate::year(t), quarter = lubridate::quarter(t), month = lubridate::month(t), week = lubridate::week(t), day = lubridate::day(t), weekdays = weekdays(t), hour = lubridate::hour(t), min = lubridate::minute(t), sec = lubridate::second(t))
    weekdays_list <- c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")
    tcomp$weekdays <- match(tcomp$weekdays, weekdays_list) - 1
    tdiff <- lapply(tcomp, FUN = tsdiff)
    for (j in names(tdiff))
    {
        if (!any(tdiff[[j]] == 0L, na.rm = TRUE))
        {
            seq.by <- j
            break
        }
        
    }
    return(seq.by)
}

##### Name a Time Series #####
#' @description The function `\code{tsname}` can be used to extract or specify the name of a time series.
#' @rdname ts-functions
#' @usage ## Get or set the name of a time series
#' tsname(x, x.name = NULL)
#' @param x a univariate time series object or a numeric vector or matrix.
#' @param x.name a new name for \code{x}. If the parameter is omitted, the current name of the time series will be returned to the user.
#' @author Ka Yui Karl Wu
#' @details To set a new name for a time series, the function must be assigned to an object. Otherwise, the new name will not be taken over.
#' @example inst/examples/tsname/tsname.R
#' @return If \code{x.name} is \code{NULL}, the attribute \code{series.name} of \code{x} will be returned. Otherwise, the series will be returned with a new value for the attribute \code{series.name} specified by \code{x.name}.
#' @export
tsname <- function(x, x.name = NULL)
{
    if (is.null(x.name))
    {
        return(attr(x, "series.name"))
    }
    else
    {
        attr(x, "series.name") <- x.name
        return(x)
    }
}

##### Copy Attributes of a Time Series #####
#' @rdname ts-functions
#' @param x.orig a univariate time series object whose attributes will be transferred to \code{x}.
#' @usage ## Copy Attributes from a Time Series to Another
#' tsattrcopy(x, x.orig)
#' @example inst/examples/tsattrcopy/tsattrcopy.R
#' @return For \code{tsattrcopy}, the function does the same as \code{\link{attributes}}. However, \code{\link{attributes}} only works if both \code{x} and \code{x.orig} share the same length, whereas \code{tsattrcopy} does not require this property and returns \code{x} with all the attributes originated from the series \code{x.orig}.
#' @seealso \code{\link{time}}, \code{\link{frequency}}, \code{\link{attributes}}
#' @export
tsattrcopy <- function(x, x.orig)
{
    if (!is.ts(x.orig)) {stop("'x.orig' must be a time series object")}
    if (is.matrix(x))
    {
        if (ncol(x) > 1) {print("More than 1 column found in 'x'. Only the first column will be considered.")}
        x <- as.vector(x[, 1])
    }
    if (length(x) == length(x.orig))
    {
        attributes(x) <- attributes(x.orig)
    }
    else
    {
        x <- ts(data = x, start = start(x.orig), frequency = frequency(x.orig))
        attr(x, "series.name") <- attr(x.orig, "series.name")
        attr(x, "seasonal.cycle") <- attr(x.orig, "seasonal.cycle")
    }
    return(x)
}

##### Generate Simple Moving Average #####
#' Generate Moving Averages of a Time Series
#' @description The function `\code{tsmovav}` calculates the moving averages of a time series.
#' @name tsmovav
#' @rdname tsmovav
#' @param x a univariate time series object or a numeric vector or matrix, or a `tsmovav` object.
#' @param order moving average order. Default is \code{3}.
#' @param type type of moving average to be calculated. Available options are "\code{backward}" and "\code{center}". While \code{backward} assigns the moving averages to the next period after the averaging window, which is more useful for forecasting purpose, \code{center} assigns the moving averages to the middle period of the averaging window, which is more suitable for time series smoothing. Default is `\code{backward}`.
#' @param n.ahead number of forecasting periods. Only useful if "\code{type = backward}". Default is \code{0}.
#' @param x.name a new name for \code{x}. If the parameter is omitted, the current name of the time series will be returned to the user.
#' @param show.plot logical. If \code{TRUE}, the smoothing/forecasting plot will be displayed directly. Default is \code{TRUE}.
#' @author Ka Yui Karl Wu
#' @details Centred moving averages are better suited for smoothing a time series than for forecasting. By definition, each moving average is aligned with the midpoint of its averaging window. When the number of periods in the averaging window (i.e., the moving average order) is even, the averages calculated for the two central positions must be combined. Specifically, the mean of these two middle moving averages is assigned to the central period that lies closest to the true midpoint of the series, forming its final centred moving average.
#' @details Mathematically, for odd number order \eqn{r}:
#' @details \deqn{\tilde{y}_t = \dfrac{y_{t-\frac{r-1}{2}}+\ldots+y_{t+\frac{r+1}{2}}}{r}}
#' @details For even number order \eqn{r}:
#' @details \deqn{\tilde{y}_t = \dfrac{0.5y_{t-\frac{r}{2}}+y_{t-\frac{r}{2}+1}+\ldots+y_{t+\frac{r}{2}-1}+0.5y_{t+\frac{r}{2}}}{2r}}
#' @details Backward moving average is a forecasting method that assigns each computed average to the period immediately following the observation window. This approach works the same regardless of whether the moving average order is odd or even.
#' @details \deqn{\tilde{y}_t = \dfrac{y_{t-1}+\ldots+x_{y-r}}{r}}
#' @return \item{x}{original series data}
#' @return \item{x.time}{list of time in which the series values were observed.} 
#' @return \item{x.timegap}{time gap between the series and forecasted values.}
#' @return \item{x.name}{name of the time series for which forecasts was requested.}
#' @return \item{pred}{predicted past values and forecasted future values.}
#' @return \item{pred.time}{list of time in which the predictions/forecasts were estimated.} 
#' @return \item{pred.name}{name of the series containing the predicted/forecasted values.}
#' @return \item{se}{standard errors of the forecasted values.}
#' @return \item{cil, ciu}{lower and upper limits of the prediction interval.}
#' @return \item{n.ahead}{number of forecasting periods.}
#' @return \item{forecast.incl}{indication of the series part that should be predicted or forecasted.}
#' @return \item{log}{logical. Indicates whether series values are log-transformed for model fitting or not.}
#' @return \item{alpha}{significance level.}
#' @return \item{order}{moving average order.}
#' @return \item{type}{type of moving average.}
#' @references Hyndman, R. J., & Athanasopoulos, G. (2021). Forecasting: Principles and practice (3rd ed.). OTexts. \cr \url{https://otexts.com/fpp3/} 
#' @example inst/examples/tsmovav/tsmovav.R
#' @importFrom stats filter
#' @importFrom utils head
#' @importFrom utils tail
#' @export
tsmovav <- function(x, order = 3, type = c("backward", "center"), n.ahead = 0, x.name = NULL, show.plot = TRUE)
{
    movavtype <- tolower(match.arg(type))
    x.name <- if (is.null(x.name) & is.null(tsname(x))) {deparse1(substitute(x))} else if (is.null(x.name) & !is.null(tsname(x))) {tsname(x)} else {x.name}
    if (movavtype == "backward")
    {
        x.arima <- tsarima(x, order = c(order, 0, 0), include.const = FALSE, method = "CSS", fixed = rep(1 / order, order))
        plot.style <- if (n.ahead > 0) {"all"} else {"predict"}
        n.ahead <- if (n.ahead <= 0) {1} else {n.ahead}
        x.ma.fc <- get_forecast(object = x.arima, n.ahead = n.ahead, show.plot = FALSE, x.name = x.name, pred.name = "Moving Averages")
        x.ma.fc$pred[1:order] <- rep(NA, order)
    }
    else if (movavtype == "center")
    {
        x.ma.1 <- filter(x, filter = rep(1 / order, order), sides = 2)
        if (order %% 2 == 0) {
            x.ma.2 <- c(NA, head(x.ma.1, length(x) - 1))
            x.ma <- (x.ma.1 + x.ma.2) / 2
        }
        else
        {
            x.ma <- x.ma.1
        }
        if (n.ahead > 0) {warning("Warning: No forecasting can be conducted for centred moving averages.")}
        x.ma.cil <- x.ma.ciu <- x.ma.se <- NULL
        xtime <- tstime(x)
        attributes(x.ma) <- attributes(x)
        x.ma.fc <- list(x = x, x.time = xtime$time, x.timegap = xtime$frequency, x.name = x.name, pred = x.ma, pred.time = xtime$time, pred.name = "Moving Averages")
    }
    out <- c(x.ma.fc, list(order = order, type = movavtype))
    if (show.plot)
    {
        plot.tsmovav(out)
    }
    return(structure(out, class = "tsmovav"))
}

##### Print Simple Moving Average #####
#' @rdname tsmovav
#' @param digits the number of significant digits.
#' @param ... other printing or plotting parameters.
#' @exportS3Method
print.tsmovav <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    args <- as.list(match.call())
    args$object <- x
    do.call(tsfprint, args[-1L])
}

##### Plot Simple Moving Average #####
#' @rdname tsmovav
#' @param title title of the moving average plot. Default is \code{NULL}.
#' @exportS3Method 
plot.tsmovav <- function(x, title = NULL, ...)
{
    if (is.null(title))
    {
        title <- if (x$type == "backward") {"Backward Moving Average"} else if (x$type == "center") {"Centred Moving Average"}
    }
    tslineplot(x = x$x, t = x$x.time, pred = x$pred, pred.t = x$pred.time, 
               cil = x$cil, ciu = x$ciu, ci.t = x$pred.time, 
               title = title, pred.col = "red", x.name = x$x.name, pred.name = x$pred.name, ...)
}

##### Model Evaluation #####
#' Goodness of Fit of a Time Series Model
#' @description The function `\code{tsmodeleval}` can be used to evaluate the goodness of fit of a time series model.
#' @name tsmodeleval
#' @rdname tsmodeleval
#' @param object a time series model of class `\code{tsarima}`, `\code{tsesm}`, or `\code{tsmovav}`. It can also be a list with at least the two elements: `\code{x}` and `\code{fitted}`.
#' @author Ka Yui Karl Wu
#' @returns A list with the following model evaluation criteria:
#' @return \item{ME}{mean error}
#' @return \item{RMSE}{Root mean square error}
#' @return \item{MAE}{mean absolute error}
#' @return \item{MPE}{mean percentage error}
#' @return \item{MAPE}{mean absolute percentage error}
#' @return \item{MASE}{mean absolute scaled error}
#' @return \item{MASE.S}{seasonal mean absolute scaled error}
#' @return \item{ACF1}{lag 1 autocorrelation}
#' @references Hyndman, R. J., Athanasopoulos, G. (2021). Forecasting: Principles and practice (3rd ed.). OTexts. \cr \url{https://otexts.com/fpp3/} 
#' @example inst/examples/tsmodeleval/tsmodeleval.R
#' @export
tsmodeleval <- function(object)
{
    if (inherits(object, "tsarima") | inherits(object, "tsesm"))
    {
        mfit <- object$error
        rownames(mfit) <- paste0("Training set (", object$train.prop * 100, "%)")
        if ("model.test" %in% names(object))
        {
            testfit <- object$model.test$error.test
            rownames(testfit) <- paste0("Testing set  (", (1 - object$train.prop) * 100, "%)")
            mfit <- rbind(mfit, testfit)
        }
    }
    else
    {
        if ("x.used" %in% names(object))
        {
            x <- object$x.used
        }
        else
        {
            logts <- "log" %in% names(object)
            x <- if (logts) {log(object$x)} else {object$x}
        }
        if ("fitted" %in% names(object))
        {
            fitted <- object$fitted
        }
        else if ("pred" %in% names(object))
        {
            fitted <- object$pred
        }
        else
        {
            stop("Model evaluation is not possible. Fitted values missing.")
        }
        xfreq <- max(frequency(x), frequency(fitted))
        error <- x - fitted
        me <- mean(error, na.rm = TRUE)
        mae <- mean(abs(error), na.rm = TRUE)
        rmse <- sqrt(mean(error ^ 2L, na.rm = TRUE))
        mpe <- mean(error / x, na.rm = TRUE) * 100
        mape <- mean(abs(error) / x, na.rm = TRUE) * 100
        mase <- mean(abs(error), na.rm = TRUE) / mean(abs(diff(x)), na.rm = TRUE)
        irseasonal <- c("day", "week", "min", "sec")
        if (xfreq > 1 & !attr(x, "seasonal.cycle") %in% irseasonal) {mases <- mean(abs(error)) / mean(abs(diff(x, xfreq)))}
        acf1 <- tsacf(error, lag.max = 1, show.plot = FALSE)
        mfit <- as.data.frame(c(if ("aic" %in% names(object)) {list(AIC = object$aic)}, list(ME = me, RMSE = rmse, MAE = mae, MPE = mpe, MAPE = mape, MASE = mase), if (xfreq > 1 & !attr(x, "seasonal.cycle") %in% irseasonal) {list(MASE.S = mases)}, list(ACF1 = acf1$acf[2])))
        rownames(mfit) <- "Error"
    }
    return(mfit)
}

##### Extract Time Series Date #####
tstimeformat <- function(x, timegap)
{
    if (timegap == "year")
    {
        format1 <- "%Y"
        format2 <- ""
    }
    else if (timegap == "quarter")
    {
        format1 <- "%Y"
        format2 <- paste0("Q", ceiling(lubridate::month(x) / 3))
    }
    else if (timegap == "month")
    {
        format1 <- "%Y/%m"
        format2 <- ""
    }
    else if (timegap == "week")
    {
        format1 <- "%Y-%m-%d %U"
        format2 <- ""
    }
    else if (timegap == "weekdays")
    {
        format1 <- "%Y-%m-%d %A"
        format2 <- ""
    }
    else if (timegap == "day")
    {
        format1 <- "%Y-%m-%d"
        format2 <- ""
    }
    else if (timegap == "hour")
    {
        format1 <- "%Y-%m-%d %H:%M:%OS"
        format2 <- ""
    }
    else if (timegap == "min" | timegap == "sec")
    {
        format1 <- "%H:%M:%OS"
        format2 <- ""
    }
    xtime <- paste0(format(x, format = format1), format2)
    return(xtime)
}

##### Calculate Drift Term #####
driftx <- function(x, order = 0, order.D = 0, period = NULL)
{
    if (is.null(period)) {period <- frequency(x)}
    n <- length(x)
    const <- rep(1, n)
    if (order.D > 0)
    {
        ncycle <- ceiling(n / period)
        sconst <- rep(1, ncycle)
        for (i in 1:order.D)
        {
            sconst <- cumsum(sconst)
        }
        const <- rep(sconst, each = period, length.out = n)
    }
    if (order > 0)
    {
        for (i in 1:order)
        {
            const <- cumsum(const)
        }
    }
    return(const)
}

##### Airport Dataset #####
#' Airport Travellers Time Series.
#'
#' A dataset containing the total number of monthly travellers recorded
#' in an unnamed airport between January 2013 and December 2019.
#'
#' @format A data frame with 84 rows and 2 variables:
#' \describe{
#'   \item{Date}{Months in which the traveller numbers are recorded.}
#'   \item{Travellers}{Total travellers in the airport of that month.}
#'   \item{AvgTemp}{Monthly average temperature (in °C) of the city where the airport is located.}
#'   \item{AvgRain}{Monthly average rain fall (in mm) in the city where the airport is located.}
#' }
#' @examples
#' data(airport)
#' head(airport)
"airport"