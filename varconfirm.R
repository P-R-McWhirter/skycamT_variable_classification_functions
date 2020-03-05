varconfirm <- function(){
  
  varpresent <- rep(0, nrow(donevarset_featgen_classified_set_selcut_good))
  
  varpos <- rep(0, nrow(donevarset_featgen_classified_set_selcut_good))
  
  for (i in 1:nrow(donevarset_featgen_classified_set_selcut_good)){
    
    if (i %% 1000 == 0){
      
      print(i)
      
    }
    
    geodist <- sqrt((donevarset_featgen_classified_set_selcut_good$Declination[i] - gaiavariables$DEC)^2.0 + (donevarset_featgen_classified_set_selcut_good$Right.Ascension[i] - gaiavariables$RA)^2.0)
    
    mindist <- min(geodist)*60*60
    
    minres <- which.min(geodist)
    
    if (mindist <= 148){
      
      varpresent[i] <- 1
      
      varpos[i] <- minres
      
    }
    
  }
  
  varpresent
  
  varpos
  
}

varapply <- function(data){
  
  data$GaiaClass <- as.character(as.vector(data$GaiaClass))
  
  for (i in 1:nrow(data)){
    
    geodist <- sqrt((data$Declination[i] - gaiavariables$DEC)^2.0 + (data$Right.Ascension[i] - gaiavariables$RA)^2.0)
    
    mindist <- min(geodist)*60*60
    
    minres <- which.min(geodist)
    
    if (mindist <= 300){
    
      data$GaiaClass[i] <- as.character(gaiavariables$Class[minres])
    
    }
    else{
      
      data$GaiaClass[i] <- NA
      
    }
    
  }
  
  data
  
}

