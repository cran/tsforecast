tsforecast(tsarima(airport$Travellers, order = c(1, 1, 0), 
                   seasonal = c(0, 1, 1), log = TRUE, include.const = TRUE), 
           n.ahead = 6, forecast.incl = "all", nobs.incl = 12, 
           title = "Forecast of Travellers by ARIMA")

tsforecast(tsesm(airport$Travellers, order = "holt-winters"), 
           n.ahead = 6, title = "Forecast of Travellers by Holt-Winters")
