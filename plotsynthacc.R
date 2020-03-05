plotsynthacc <- function(data, grape = T, legpos = "bottomright", title = ""){
  
  ep <- seq(0, 1, by = 0.01)
  
  hit <- rep(0, length(ep))
  
  mult <- rep(0, length(ep))
  
  smult <- rep(0, length(ep))
  
  oneal <- rep(0, length(ep))
  
  halfal <- rep(0, length(ep))
  
  unknown <- rep(0, length(ep))
  
  for (i in 1:length(ep)){
    
    test <- synthpermodes(data, ep = ep[i], grape = grape, quiet = T)
    
    hit[i] <- sum(test$Hit)/nrow(test)
    
    mult[i] <- sum(test$Multiple)/nrow(test)
    
    smult[i] <- sum(test$Sub.Multiple)/nrow(test)
    
    oneal[i] <- sum(test$Alias.1)/nrow(test)
    
    halfal[i] <- sum(test$Alias.0.5)/nrow(test)
    
    unknown[i] <- sum(test$Unknown)/nrow(test)
    
  }
  
  dots <- which((1:length(ep)) %% 3 == 0)
  
  plot(ep, hit, type = "l", ylim = c(0, 1), xlim = c(0, 1), lwd = 2, col = "red", main = paste(title, " - Performance vs tolerance", sep=""), ylab = "Performance (Hit Rate)", xlab = "Tolerance", pch = 19, cex.lab=2, cex.axis=2, cex.main=1.7, cex.sub=2)

  lines(ep[dots], hit[dots], type = "p", pch = 19, col = "red")
  
  lines(ep, mult, type = "l", lwd = 2, col = "blue")
  
  lines(ep[dots], mult[dots], type = "p", pch = 18, col = "blue")
  
  lines(ep, smult, type = "l", lwd = 2, col = "dark green")
  
  lines(ep[dots], smult[dots], type = "p", pch = 17, col = "dark green")
  
  lines(ep, oneal+halfal, type = "l", lwd = 2, col = "orange")
  
  lines(ep[dots], (oneal+halfal)[dots], type = "p", pch = 15, col = "orange")
  
  abline(a=0, b=1, lty=2, lwd=3)
  
  pchset <- c(19, 18, 17, 15)
  
  legend(legpos, inset = .05, cex = 2, title = "Modes", c("Hit", "Multiple", "Submultiple", "Alias"), lty = c(1,1), lwd = c(2,2), col = c("red", "blue", "dark green", "orange"), pch = pchset, bg = "grey96")
  
}






plotsynthhist <- function(data1, data2, breaks = 50, xlim = c(0, 1500), double = FALSE, legpos = "topright", maintitle = ""){
  
  if (double == TRUE){
    
    data1$Vuong.Period <- 2*data1$Vuong.Period
    
    data2$Vuong.Period <- 2*data2$Vuong.Period
    
  }
  
  p1 <- hist(data1$Vuong.Period, breaks = max(data1$Vuong.Period) / breaks, plot = F)
  
  p1$counts[which(p1$counts > 0)] <- log10(p1$counts[which(p1$counts > 0)])
  
  p2 <- hist(data2$Vuong.Period, breaks = max(data2$Vuong.Period) / breaks, plot = F)
  
  p2$counts[which(p2$counts > 0)] <- log10(p2$counts[which(p2$counts > 0)])
  
  p3 <- hist(data1$Input.Period, breaks = max(data1$Input.Period) / breaks, plot = F)
  
  p3$counts <- log10(p3$counts)
  
  plot( p1, col=rgb(1,0,0,1/4), xlim = xlim, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, main = paste("Histogram of ", maintitle, sep = ""), xlab = "Estimated Period (days)", ylab = "Log Frequency")
  
  plot( p2, col=rgb(0,0,1,1/4), add=T)
  
  plot( p3, col=rgb(1,1,1,0), add=T)
  
  legend(legpos, inset = .05, cex = 2, title = "Method", c("GRAPE", "BGLS Periodogram"), lty = c(1,1), lwd = c(4,4), col = c(rgb(1,0,0,1/2), rgb(0,0,1,1/2)), pch = "", bg = "grey96")
  
  
}