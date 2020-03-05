randseedtest <- function(name, seed, redo = 10){
  
  n <- length(seed)
  
  pers <- rep(0, n)
  
  for (i in 1:n){
  
    featvec <- featgen(AAVSO_trainset[AAVSO_trainset$Name == name,], radii = 0.042, quiet = F, seed = seed[i], lctrenddata = lctrenddata, sigk = 4, nonper = F, redo = redo)
  
    aavper <- featvec$AAVSO.Period
    
    pers[i] <- featvec$Period
    
  }
  
  print(as.numeric(as.vector(aavper)))
  
  print(mean(as.numeric(as.vector(pers))))
  
  print(sd(as.numeric(as.vector(pers))))
  
  pers
  
}