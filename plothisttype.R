plothisttype <- function(data, type, step = 0.1, sub = NULL, ylim = c(0, 350), xlim = c(0, 2), legpos = "topright", main = "", xlab = ""){
  
  breaks <- seq(min(as.numeric(data))-10, max(as.numeric(data))+10, by = step)
  
  type <- as.factor(type)
  
  len <- length(levels(type))
  
  cols <- rainbow(len, alpha = 0.25)
  
  legcols <- rainbow(len, alpha = 0.5)
  
  hists <- list()
  
  if (is.null(sub)){
  
    for (i in 1:len){
    
      hists[[i]] <- hist(data[which(type == levels(type)[i])], breaks = breaks, plot = F)
    
      tot <- length(data[which(type == levels(type)[i])])
    
      hists[[i]]$counts <- hists[[i]]$counts / tot
    
    }
    
  }
  else{
    
    for (i in sub){
      
      hists[[i]] <- hist(data[which(type == levels(type)[i])], breaks = breaks, plot = F)
      
      tot <- length(data[which(type == levels(type)[i])])
      
      hists[[i]]$counts <- hists[[i]]$counts / tot
      
    }
    
  }
  
  
  
  if (is.null(sub)){
    
    plot( hists[[1]], col=cols[1], ylim = ylim, xlim = xlim, cex.lab=1.2, 
          cex.axis=1.2, cex.main=1.2, cex.sub=1.2, main = main,
          xlab = xlab, ylab = "Density")
  
    for (i in 2:len){
    
      plot(hists[[i]], col = cols[i], add=T)
    
    }
    
    legend(legpos, inset = .05, cex = 1.2, title = "Type", levels(type),
           lty = c(1,1), lwd = c(4,4), col = legcols, pch = "", bg = "grey96")
    
  }
  else{
    
    plot( hists[[sub[1]]], col=cols[sub[1]], ylim = ylim, xlim = xlim, cex.lab=1.2, 
          cex.axis=1.2, cex.main=1.2, cex.sub=1.2, main = main,
          xlab = xlab, ylab = "Density")
    
    for (i in sub[-1]){
      
      plot(hists[[i]], col = cols[i], add=T)
      
    }
    
    legend(legpos, inset = .05, cex = 1.2, title = "Type", levels(type)[sub],
           lty = c(1,1), lwd = c(4,4), col = legcols[sub], pch = "", bg = "grey96")
    
  }
  
}






plotpcaset <- function(data, p, seed = 1, lev = 10, ylim = c(1, -0.4), legpos = "bottomright"){
  
  polyset <- try(lcgbinpolyfit(data, spur = T, pero = p, delta = 1, seedno = seed), TRUE)[38:136]
  
  set <- list()
  
  print((as.matrix(polyset)) %*% as.matrix(polyfit2.pca$rotation[,1:10]))
  
  for (i in 1:lev){
    
    hold <- (as.matrix(polyset)) %*% as.matrix(polyfit2.pca$rotation[,1:i])
    
    set[[i]] <- as.vector(as.matrix(polyfit2.pca$rotation[,1:i]) %*% t(as.matrix(hold)))
    
  }
  
  plot(seq(-0.49, 0.49, by = 0.01), set[[1]], type = "l", lwd = 2, col = rainbow(lev, alpha = 0.5)[1], ylim = ylim, main = paste(data$Name, " Principal Component Reconstruction", sep=""), ylab = "Mean-Zeroed Apparent Magnitude", xlab = paste("Phase at a period of ", p, " days", sep = ""))
  
  for (i in 2:lev){
    
    lines(seq(-0.49, 0.49, by = 0.01), set[[i]], lwd = 2, col = rainbow(lev, alpha = 0.5)[i])
    
  }
  
  hold <- (as.matrix(polyset)) %*% as.matrix(polyfit2.pca$rotation)
  
  set <- as.vector(as.matrix(polyfit2.pca$rotation) %*% t(as.matrix(hold)))
  
  lines(seq(-0.49, 0.49, by = 0.01), set, lwd = 2, col = "black")
  
  legend(legpos, inset = .05, cex = 1.2, title = "PCA", c("PCA1", paste("PCA1-", 2:10, sep = ""), "Original"),
         lty = c(1,1), lwd = c(4,4), col = c(rainbow(10, alpha = 0.5), "black"), pch = "", bg = "grey96")
  
}