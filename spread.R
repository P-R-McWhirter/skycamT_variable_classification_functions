spread <- function(yy, y, n, x, m) {
  
  nfac=c(1,1,2,6,24,120,720,5040,40320,362880)
  
  ix <- x
  
  if (x == as.single(x)){
    
    yy[ix] <- yy[ix] + y
    
  }
  else
  {
    
    ilo = min(max(round(x - 0.5*m + 1), 1), n-m+1)
    
    ihi <- iho + m - 1
    
    nden <- nfac(m)
    
    fac <- x - ilo
    
    j <- (ilo+1):ihi
    
    fac <- fac * cumprod(x-j)
    
    fac <- fac(length(fac))
    
    yy(ihi) <- yy(ihi) + y*fac/(nden*(x-ihi))
    
    for (j in seq(from = (ihi - 1), to = ilo, by = -1)){
      
      if ((nden/(j+1-ilo)) < 0){
        
        nden <- ceiling(nden/(j+1-ilo))*(j - ihi)
        
      }
      else
      {
        
        nden <- floor(nden/(j+1-ilo))*(j - ihi)
        
      }
      
      yy(j) <- yy(j) + y*fac/(nden*(x-j))
      
    }
    
  }
  
  yy
  
}