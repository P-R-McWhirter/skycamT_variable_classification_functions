photcompute <- function(fieldwise = F){
  
  data <- list()
  
  for (i in 1:120){
    
    data[[i]] <- read.table(paste("F:/Documents/BLAP/catalogue/h_e_20180801_45_", i, "_1_1.cat", sep=""), quote="\"", stringsAsFactors=FALSE)
    
  }
  
  if (fieldwise == F){
    
    out <- list()
    
    for (i in 1:559){
      
      print(i)
      
      out[[i]] <- data[[1]][i,,drop=F]
      
      for (j in 2:120){
        
        out[[i]] <- rbind(out[[i]], data[[j]][i,])
        
      }
      
    }
    
    out
    
  }
  else{
    
    data
    
  }
  
}





photcompute2 <- function(fieldwise = F, band = "g"){
  
  data <- list()
  
  if (band == "g"){
    
    z <- 54
    
    j <- 1
  
    for (i in 1:45){
    
      data[[j]] <- read.table(paste("F:/Documents/BLAP/catalogue2/h_e_20180802_", z, "_1_1_1.cat", sep=""), quote="\"", stringsAsFactors=FALSE)
    
      data[[j+1]] <- read.table(paste("F:/Documents/BLAP/catalogue2/h_e_20180802_", z, "_2_1_1.cat", sep=""), quote="\"", stringsAsFactors=FALSE)
    
      j <- j + 2
    
      z <- z + 2
    
    }
    
  }
  else{
    
    z <- 55
    
    for (i in 1:45){
      
      data[[i]] <- read.table(paste("F:/Documents/BLAP/catalogue2/h_e_20180802_", z, "_1_1_1.cat", sep=""), quote="\"", stringsAsFactors=FALSE)
      
      z <- z + 2
      
    }
    
  }
  
  if (fieldwise == F){
    
    out <- list()
    
    if (band == "g"){
    
      for (i in 1:537){
      
        print(i)
      
        out[[i]] <- data[[1]][i,,drop=F]
      
        for (j in 2:90){
        
          out[[i]] <- rbind(out[[i]], data[[j]][i,])
        
        }
      
      }
      
    }
    else{
      
      for (i in 1:932){
        
        print(i)
        
        out[[i]] <- data[[1]][i,,drop=F]
        
        for (j in 2:45){
          
          out[[i]] <- rbind(out[[i]], data[[j]][i,])
          
        }
        
      }
      
    }
    
    out
    
  }
  else{
    
    data
    
  }
  
}





photcompute3 <- function(fieldwise = F, band = "i"){
  
  data <- list()
  
  if (band == "i"){
    
    j <- 1
    
    for (i in seq(33, 62, by = 2)){
      
      data[[j]] <- read.table(paste("F:/Documents/BLAP/catalogue3/h_e_20180807_", i, "_1_1_1.cat", sep=""), quote="\"", stringsAsFactors=FALSE)
      
      j <- j + 1
      
    }
    
  }
  else{
    
    j <- 1
    
    for (i in seq(33, 62, by = 2)+1){
      
      data[[j]] <- read.table(paste("F:/Documents/BLAP/catalogue3/h_e_20180807_", i, "_1_1_1.cat", sep=""), quote="\"", stringsAsFactors=FALSE)
      
      j <- j + 1
      
    }
    
  }
  
  if (fieldwise == F){
    
    out <- list()
    
    if (band == "i"){
      
      for (i in 1:1457){
        
        if (i %% 100 == 0){
          
          print(i)
          
        }
        
        out[[i]] <- data[[1]][i,,drop=F]
        
        for (j in 2:15){
          
          out[[i]] <- rbind(out[[i]], data[[j]][i,])
          
        }
        
      }
      
    }
    else{
      
      for (i in 1:1313){
        
        if (i %% 100 == 0){
          
          print(i)
          
        }
        
        out[[i]] <- data[[1]][i,,drop=F]
        
        for (j in 2:15){
          
          out[[i]] <- rbind(out[[i]], data[[j]][i,])
          
        }
        
      }
      
    }
    
    out
    
  }
  else{
    
    data
    
  }
  
}






photcompute4 <- function(fieldwise = F, band = "g"){
  
  data <- list()
  
  if (band == "g"){
    
    j <- 1
    
    for (i in seq(18, 108, by = 2)){
      
      data[[j]] <- read.table(paste("F:/Documents/BLAP/catalogue4/h_e_20180826_", i, "_1_1_1.cat", sep=""), quote="\"", stringsAsFactors=FALSE)
      
      data[[j+1]] <- read.table(paste("F:/Documents/BLAP/catalogue4/h_e_20180826_", i, "_2_1_1.cat", sep=""), quote="\"", stringsAsFactors=FALSE)
      
      data[[j+2]] <- read.table(paste("F:/Documents/BLAP/catalogue4/h_e_20180826_", i, "_3_1_1.cat", sep=""), quote="\"", stringsAsFactors=FALSE)
      
      j <- j + 3
      
    }
    
  }
  else{
    
    j <- 1
    
    for (i in seq(18, 108, by = 2)+1){
      
      data[[j]] <- read.table(paste("F:/Documents/BLAP/catalogue4/h_e_20180826_", i, "_1_1_1.cat", sep=""), quote="\"", stringsAsFactors=FALSE)
      
      j <- j + 1
      
    }
    
  }
  
  if (fieldwise == F){
    
    out <- list()
    
    if (band == "g"){
      
      for (i in 1:2987){
        
        if (i %% 100 == 0){
          
          print(i)
          
        }
        
        out[[i]] <- data[[1]][i,,drop=F]
        
        for (j in 2:138){
          
          out[[i]] <- rbind(out[[i]], data[[j]][i,])
          
        }
        
      }
      
    }
    else{
      
      for (i in 1:4344){
        
        if (i %% 100 == 0){
          
          print(i)
          
        }
        
        out[[i]] <- data[[1]][i,,drop=F]
        
        for (j in 2:46){
          
          out[[i]] <- rbind(out[[i]], data[[j]][i,])
          
        }
        
      }
      
    }
    
    out
    
  }
  else{
    
    data
    
  }
  
}







photcompute5 <- function(fieldwise = F, band = "g"){
  
  data <- list()
  
  if (band == "g"){
    
    j <- 1
    
    for (i in seq(18, 108, by = 2)){
      
      data[[j]] <- read.table(paste("F:/Documents/BLAP/catalogue5/h_e_20180826_", i, "_1_1_1.cat", sep=""), quote="\"", stringsAsFactors=FALSE)
      
      data[[j+1]] <- read.table(paste("F:/Documents/BLAP/catalogue5/h_e_20180826_", i, "_2_1_1.cat", sep=""), quote="\"", stringsAsFactors=FALSE)
      
      data[[j+2]] <- read.table(paste("F:/Documents/BLAP/catalogue5/h_e_20180826_", i, "_3_1_1.cat", sep=""), quote="\"", stringsAsFactors=FALSE)
      
      j <- j + 3
      
    }
    
  }
  else{
    
    j <- 1
    
    for (i in seq(18, 108, by = 2)+1){
      
      data[[j]] <- read.table(paste("F:/Documents/BLAP/catalogue5/h_e_20180826_", i, "_1_1_1.cat", sep=""), quote="\"", stringsAsFactors=FALSE)
      
      j <- j + 1
      
    }
    
  }
  
  if (fieldwise == F){
    
    out <- list()
    
    if (band == "g"){
      
      for (i in 1:3037){
        
        if (i %% 100 == 0){
          
          print(i)
          
        }
        
        out[[i]] <- data[[1]][i,,drop=F]
        
        for (j in 2:138){
          
          out[[i]] <- rbind(out[[i]], data[[j]][i,])
          
        }
        
      }
      
    }
    else{
      
      for (i in 1:4991){
        
        if (i %% 100 == 0){
          
          print(i)
          
        }
        
        out[[i]] <- data[[1]][i,,drop=F]
        
        for (j in 2:46){
          
          out[[i]] <- rbind(out[[i]], data[[j]][i,])
          
        }
        
      }
      
    }
    
    out
    
  }
  else{
    
    data
    
  }
  
}