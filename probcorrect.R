probcorrect <- function(data, type = "class", levs = "", seedno = 20, pop = 100, pairups = 20, nogen = 50, crossover = 0.65, mutation = 0.3, fdif = 0.6, dfrac = 0.7){
  
  # This function uses a genetic algorithm to optimise A and B using Brier Score for the
  # correction of the probabilities from a classifier.
  
  if (type == "class"){
    
    trained <- paste(data$Train.Class, ".prob", sep="")
    
    probs <- data[,17:ncol(data)]
    
  }
  else if (type == "crossval"){
    
    trlevs <- paste(levs, ".prob", sep="")
    
    colnames(data) <- c("Type", trlevs)
    
    trained <- trlevs[as.numeric(data[,1])]
    
    probs <- data[,2:ncol(data)]
    
  }
  else{
    
    stop(paste("'", type, "' is not a valid input type."))
    
  }
  
  trainclass <- matrix(0, nrow = nrow(probs), ncol = ncol(probs))
  
  for (z in 1:nrow(probs)){
    
    trainclass[z,which(colnames(probs) == trained[z])] <- 1
    
  }
  
  # Now that we have the probability vectors, we need to determine the coefficients for 
  # a scalar named r which is used to adjust the probabilities.
  
  oldbrierscore <- (1/nrow(probs)) * sum(sum((trainclass - probs)^2.0))
  
  genres <- briergen(probs = probs, trainclass = trainclass, seedno = seedno, pop = pop, pairups = pairups, nogen = nogen, crossover = crossover, mutation = mutation, fdif = fdif, dfrac = dfrac)
  
  A <- genres$A
  
  B <- genres$B
  
  problargest <- as.numeric(as.vector(apply(probs, 1, max)))
  
  prob2ndlargest <- as.numeric(as.vector(apply(probs, 1, function(x) sort(x)[dim(probs)[2]-1])))
  
  probmarg <- problargest - prob2ndlargest
  
  r <- 1 / (1 + exp(A*probmarg + B))
  
  corprobs <- probs*(1-r)
  
  maxprobpos <- as.numeric(as.vector(match(colnames(probs)[max.col(probs,ties.method="first")], names(probs))))
  
  for (i in 1:nrow(probs)){
    
    corprobs[i,maxprobpos[i]] <- probs[i,maxprobpos[i]] + r[i]*(1-probs[i,maxprobpos[i]])
    
  }
  
  brierscore <- (1/nrow(corprobs)) * sum(sum((trainclass - corprobs)^2.0))
  
  print(paste(A, B))
  
  print(paste(oldbrierscore, brierscore))
  
  out <- list(corprobs = corprobs, A = A, B = B)
  
}



probcorother <- function(data, type = "class", levs = "", A, B, cutoffs = NULL){
  
  # This function uses a genetic algorithm to optimise A and B using Brier Score for the
  # correction of the probabilities from a classifier.
  
  corcutoffs <- NULL
  
  if (type == "class"){
    
    trained <- paste(data$Train.Class, ".prob", sep="")
    
    probs <- data[,17:ncol(data)]
    
  }
  else if (type == "vardet"){
    
    trained <- paste(data$Train.Class, ".prob", sep="")
    
    probs <- data[,12:ncol(data)]
    
  }
  else if (type == "crossval"){
    
    trlevs <- paste(levs, ".prob", sep="")
    
    colnames(data) <- c("Type", trlevs)
    
    trained <- trlevs[as.numeric(data[,1])]
    
    probs <- data[,2:ncol(data)]
    
  }
  else{
    
    stop(paste("'", type, "' is not a valid input type.", sep=""))
    
  }
  
  trainclass <- matrix(0, nrow = nrow(probs), ncol = ncol(probs))
  
  for (z in 1:nrow(probs)){
    
    trainclass[z,which(colnames(probs) == trained[z])] <- 1
    
  }
  
  # Now that we have the probability vectors, we need to determine the coefficients for 
  # a scalar named r which is used to adjust the probabilities.
  
  oldbrierscore <- (1/nrow(probs)) * sum(sum((trainclass - probs)^2.0))
  
  problargest <- as.numeric(as.vector(apply(probs, 1, max)))
  
  prob2ndlargest <- as.numeric(as.vector(apply(probs, 1, function(x) sort(x)[dim(probs)[2]-1])))
  
  probmarg <- problargest - prob2ndlargest
  
  r <- 1 / (1 + exp(A*probmarg + B))
  
  corprobs <- probs*(1-r)
  
  maxprobpos <- as.numeric(as.vector(match(colnames(probs)[max.col(probs,ties.method="first")], names(probs))))
  
  if (!is.null(cutoffs)){
    
    cutoffs <- matrix(cutoffs, nrow = nrow(probs), ncol = length(cutoffs), byrow = TRUE)
    
    corcutoffs <- cutoffs*(1-r)
    
  }
  
  for (i in 1:nrow(probs)){
    
    corprobs[i,maxprobpos[i]] <- probs[i,maxprobpos[i]] + r[i]*(1-probs[i,maxprobpos[i]])
    
    if (!is.null(cutoffs)){
      
      corcutoffs[i,maxprobpos[i]] <- cutoffs[i,maxprobpos[i]] + r[i]*(1-cutoffs[i,maxprobpos[i]])
      
    }
    
  }
  
  brierscore <- (1/nrow(corprobs)) * sum(sum((trainclass - corprobs)^2.0))
  
  print(paste(A, B))
  
  print(paste(oldbrierscore, brierscore))
  
  if (type == "vardet"){
    
    corprobs <- cbind(data[,1:11], corprobs)
    
    corprobs$Predicted.Probability <- as.numeric(as.vector(apply(corprobs[,12:ncol(data)], 1, max)))
    
  }
  else{
    
    corprobs <- cbind(data[,1:16], corprobs)
    
    corprobs$Predicted.Probability <- as.numeric(as.vector(apply(corprobs[,17:ncol(data)], 1, max)))
  
  }
  
  out <- list(corprobs = corprobs, corcutoffs = corcutoffs)
  
  out
  
}