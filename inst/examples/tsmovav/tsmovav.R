## Backward Moving Average
tsmovav(airport$Travellers, order = 12, type = "backward", n.ahead = 6, show.plot = TRUE)

## Centered Moving Average
tsmovav(airport$Travellers, order = 12, type = "center")
