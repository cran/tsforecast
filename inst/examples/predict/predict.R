predict(tsarima(airport$Travellers, order = c(1, 1, 0), 
                seasonal = c(0, 1, 1), log = TRUE, include.const = TRUE),
        n.ahead = 6, alpha = 0.05)

predict(tsesm(airport$Travellers, order = "holt-winters"), 
        n.ahead = 6, alpha = 0.05)
