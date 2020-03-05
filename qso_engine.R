qso_engine <- function(time, data, error, ltau = 3, lvar = -1.7, sys_err = 0, return_model = FALSE) {

  # Calculates the fit quality of a damped random walk to a qso lightcurve.
  # Written by N. Butler (nat@astro.berkeley.edu), Feb. 2010.
  # Version 1.0, Translated to R, P.R. McWhirter May 2016.
  
  # The formalism is from Rybicki & Press (1994; arXiv:comp-gas/9405004)
  
  # Data are modelled with a covariance function
  # Lij = 0.5*var*tau*exp(-|time_i-time_j|/tau) .
  
  # Input:
  # time - measurement times, typically days
  # data - measured magnitudes
  # error - uncertainty in measured magnitudes
  
  # Output (dictionary):
    
  # chi2/nu - classical variability measure
  # chi2_qso/nu - for goodness of fit given fixed parameters
  # chi2_qso/nu_extra - for parameter fitting, add to chi2/nu
  # chi^2/nu_NULL - expected chi2/nu for non-qso variable
  
  # signif_qso - significance chi^2/nu<chi^2/nu_NULL (rule out false alarm)
  # signif_not_qso - significance chi^2/nu>1 (rule out qso)
  # signif_vary - significance that source is variable
  # class - resulting source type (ambiguous, not_qso, qso)
  
  #model - time series prediction for each datum given all others (iff return_model==True)
  #dmodel - model uncertainty, including uncertainty in data
  
  #Notes:
  #T = L^(-1)
  #Data variance is D
  #Full covariance C^(-1) = (L+D)^(-1) = T [T+D^(-1)]^(-1) D^(-1)
  #Code takes advantage of the tridiagonality of T and T+D^(-1).
  
  out_dict <- list(chi2_qso = 999, chi2_qso_extra = 0, signif_qso = 0, signif_not_qso = 0,
                   signif_vary = 0, chi2_qso_null = 0, chi2divnu = 0, nu = 0, class = "ambiguous")
  
  lvar0 <- log10(0.5) + lvar + ltau
  
  ln <- length(data)
  dt <- abs(time[2:ln] - time[1:(ln-1)])
  
  # first make sure all dt>0
  g <- which(dt > 0, arr.ind = TRUE)
  lg <- length(g)
  # must have atleast 2 data points
  if (lg <= 0){
    return(out_dict)
  }
  
  if (return_model == TRUE){
    model <- data
    dmodel <- -1.0 * error
  }
  
  if (lg < ln){
    dt <- dt[g]
    gg <- rep(0, (lg+1))
    gg[1:length(gg)] <- g+1
    dat <- data[gg]
    wt <- 1.0 / (sys_err^2 + error[gg]^2)
    ln <- lg+1
  }
  else{
    dat <- 1.0 * data
    wt <- 1.0 / (sys_err^2 + error^2)
  }
  
  out_dict$nu <- ln - 1
  varx <- var(dat)
  dat0 <- sum(dat * wt) / sum(wt)
  out_dict$chi2divnu <- sum((dat - dat0)^2 * wt) / out_dict$nu

  # define tridiagonal matrix T = L^(-1)
  # sparse matrix form: ab[u + i - j, j] == a[i,j]  i <= j, (here u=1)
  T = matrix(0, nrow = 2, ncol = ln)
  arg <- dt * exp(-1.0 * log(10) * ltau)
  ri <- exp(-1.0 * arg)
  ei <- ri / (1.0 + ri) / (1.0 - ri)
  T[0,] <- -1.0 * ei
  print(dim(T))
  print(length(ri))
  T[1,(ln:1)] <- 1.0 + ri * ei
  T[1,] <- T[1,] + ri * ei
  T[1,ln] <- T[1,ln] + 1
  T0 <- median(T[1,])
  T <- T / T0
  
  # equation for chi2_qso is [ (dat-x0)T Tp^(-1) D^(-1) (dat-x0) ] , where Tp = T+D^(-1) and D^(-1) = wt
  fac <- exp(log(10) * lvar0) / T0
  Tp <- 1.0 * T
  Tp[1,] <- Tp[1,] + wt * fac
  # solve Tp*z=y for z (y=wt*dat)
  
  out_dict
}