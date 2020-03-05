shapiro <- function(x){
  
    DNAME <- deparse(substitute(x))
    stopifnot(is.numeric(x))
    x <- sort(x[complete.cases(x)])
    n <- length(x)
    rng <- x[n] - x[1L]
    if (rng == 0) 
      stop("all 'x' values are identical")
    if (rng < 1e-10) 
      x <- x/rng
    res <- .Call(C_SWilk, x)
    RVAL <- list(statistic = c(W = res[1]), p.value = res[2], 
                 method = "Shapiro-Wilk normality test", data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  
}