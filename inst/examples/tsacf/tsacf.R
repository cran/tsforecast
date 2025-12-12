# Autocorrelation (ACF)
airport.acf <- tsacf(airport$Travellers, lag.max = 24, show.plot = FALSE)
print(airport.acf, digits = 4)
plot(airport.acf)

# Partial Autocorrelation (PACF)
tspacf(airport$Travellers, lag.max = 24)

# Autocovariance (ACOV)
tsacov(airport$Travellers, lag.max = 24, show.plot = FALSE)

# Cross-Correlation (CCF)
tsccf(airport$AvgRain, airport$Travellers, lag.max = 24)

# Cross-Covariance (CCOV)
tsccov(airport$AvgRain, airport$Travellers, lag.max = 24)
