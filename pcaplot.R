pcaplot <- function(x, pca, num, spur = F, pero = 1.0) {
  
  x <- x[1,]
  
  spline <- lcsplinefit(x, spur = spur, pero = pero)
  
  name <- spline$Name
  
  type <- spline$Type
  
  period <- spline$Period
  
  if (spur == T){
    
    period <- pero
    
  }
  
  spline <- spline[8:106]
  
  splinepca <- as.matrix(spline) %*% as.matrix((pca$rotation[,1:num]))
  
  print(splinepca)
  
  spline2 <- as.matrix(splinepca) %*% t(as.matrix(pca$rotation[,1:num]))
  
  plot(seq(-0.45, 0.5, length.out = 99), spline2, type = "l", main = paste(name, " (", type, ") folded light curve with ", num, " PCAs", sep = ""),
       ylab = "Normalised Magnitude", xlab = paste("Phase at a period of ", period, " days.", sep = ""), ylim = c(1, 0))
  
  lines(seq(-0.45, 0.5, length.out = 99), spline, col = "red")
  
}