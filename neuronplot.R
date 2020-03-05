neuronplot <- function(x) {
  
  x <- as.matrix(x)
  
  parnum <- dim(x)[1]
  
  sqparn <- sqrt(parnum)
  
  neunum <- dim(x)[2]
  
  for (i in 1:neunum){
    
    neuron <- as.matrix(x[,i])
    
    dim(neuron) <- c(sqparn, sqparn)
    
    neuron <- t(neuron)
    
    neuron <- (neuron / (max(abs(neuron)) * 2)) + 0.5
    
    par(mfrow=c(1, 1))
    
    png(paste("F:/Documents/PhD/Neurons/Neuron-", i, ".png", sep=""))
    
    par(mar = rep(0, 4))
    
    image(neuron, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    
    dev.off()
    
  }
  
  dim(x) <- c()
  
}