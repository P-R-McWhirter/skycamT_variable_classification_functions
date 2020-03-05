classvalid <- function(data, cut = c(0.8, 5)){
  
  data <- data[which(data$Predicted.Probability >= cut[1] & data$Anomaly.Score <= cut[2]),]
  
  data
  
}