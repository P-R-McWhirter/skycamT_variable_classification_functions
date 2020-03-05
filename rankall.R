rankall <- function(outcome, num = "best") {
  
  data <- read.csv("outcome-of-care-measures.csv", colClasses = "character")
  
  valid <- c("heart attack", "heart failure", "pneumonia")
  
  if(!(outcome %in% valid)) {
    err = "invalid outcome"
    stop(err, call. = TRUE, domain = NULL)
    geterrmessage()
  }
  
  if (outcome == "heart attack"){
    pos = 11
  }
  else if (outcome == "heart failure"){
    pos = 17
  }
  else {
    pos = 23
  }
  
  data[, pos] <- as.numeric(data[, pos])
  
  bad <- is.na(data[,pos])
  
  data <- data[!bad,]
  
  range(data[, pos], na.rm = TRUE)
  
  if (num == "worst"){
    numb = length(data[, 2])
  }
  else if (num == "best"){
    nums = 1
  }
  
  if (numb > length(data[, 2])){
    return(NA)
  }
  
  out <- split(data , f = data$State)
  
  resName = c()
  
  resState = c()
  
  for (i in 1:length(out)){
    
    
  
    out[[i]] <- out[[i]][order(out[[i]]$State, out[[i]][,pos], out[[i]][,2]),]
  
    bes <- which.min(out[[i]][, pos])
  
    if (num == "worst"){
      numb = length(out[[i]][, 2])
      resName <- c(resName, out[[i]][numb, 2])
      resState <- c(resState, out[[i]][1, 7])
    }
    else if (num == "best"){
      resName <- c(resName, out[[i]][1, 2])
      resState <- c(resState, out[[i]][1, 7])
    }
    else{
      resName <- c(resName, out[[i]][num, 2])
      resState <- c(resState, out[[i]][1, 7])
    }
  }
  
  df <- data.frame(resName, resState)
  colnames(df) <- c("hospital", "state")
  rownames(df) <- resState
  
  df
  
}