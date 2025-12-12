travellers_d <- tsdiff(x = airport$Travellers, order = 1, order.D = 1)
tsexplore(travellers_d, lag.max = 60)
