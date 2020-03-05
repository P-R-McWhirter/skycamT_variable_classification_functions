lcloader <- function(folder, alert = 500) {
  
  files <- list.files(path = folder, pattern = "*.csv")
  
  first = as.matrix(read.csv(paste(folder, "/", files[1], sep = ""), header = FALSE))
  
  x <- matrix(0, nrow = length(files), ncol = dim(first)[1]+1)
  
  for (i in 1:length(files)) {
    
    if (i %% alert == 0){
      
      print(paste("Reading file: ", i, ".", sep = ""))
      
    }
    
    vec = as.matrix(read.csv(paste(folder, "/", files[i], sep = ""), header = FALSE))
    
    x[i,] <- cbind(gsub(".csv", '', files[i]), t(vec))
    
  }
  
  x <- as.data.frame(x)
  
  x
  
}