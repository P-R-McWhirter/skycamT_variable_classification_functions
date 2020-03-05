colimputetest <- function(data, maxiter = 10, ntree = 100, seed = 1){
  
  set.seed(seed)
  
  samp <- sample(1:nrow(data), floor(nrow(data)/4))
  
  predata <- data$B.R.Colour[samp]
  
  data$B.R.Colour[samp] <- NA
  
  data <- colimpute(data)
  
  postdata <- data$B.R.Colour[samp]
  
  out <- list(predata = predata, postdata = postdata, samp = samp)
  
}