el2.test.wts <- function (u,v,wu,wv,mu0,nu0,indicmat,mean,lamOld=0) {
#If mean is not a scalar then stop
    if (length(mean) != 1) 
        stop("mean must be a scalar")

#Calculate scalars to be used in calculations
    sumwu <- sum(wu)
    sumwv <- sum(wv)
    nu <- length(u)
    nv <- length(v)
  
#Calculate matrix and vectors to be used in calculations
    indic4mu <- nu0 %*% t(indicmat) 
    indic4nu <- mu0 %*% indicmat  
#Calculate the delta for lam search
    du <- 0.02 * sumwu/abs(sum(indic4mu))
    dv <- 0.02 * sumwv/abs(sum(indic4nu))
    dd <- min(du, dv)

#Define lamfun, where lamfun(true-lam) = 0
    lamfun <- function(lam, wu, wv, sumwu, sumwv, indic4mu,
      indic4nu, indicmat) {
    mu <- wu/(sumwu+lam*indic4mu)
    nu <- wv/(sumwv+lam*indic4nu)
    return(mu %*% indicmat %*% t(nu))
      }

    flamOld <- lamfun(lamOld, wu, wv, sumwu, sumwv, indic4mu, indic4nu,indicmat)
	
	Diff <- lamfun(lamOld+dd, wu, wv, sumwu, sumwv, indic4mu, indic4nu,indicmat) - 
	            lamfun(lamOld-dd, wu, wv, sumwu, sumwv, indic4mu, indic4nu,indicmat)    
	while( (abs(Diff) < 1e-04) & (dd < 1e03) ) {
	 dd <- dd*2
	 Diff <- lamfun(lamOld+dd, wu, wv, sumwu, sumwv, indic4mu, indic4nu,indicmat) - 
	            lamfun(lamOld-dd, wu, wv, sumwu, sumwv, indic4mu, indic4nu,indicmat)
	}

	
    if ( abs(flamOld) < 1e-08 )
      lam <- lamOld
	  else 
      { 
	if ( (flamOld > 0) & (Diff > 0) ) {
            lo <- lamOld-dd
            up <- lamOld
            while ( (lamfun(lo, wu, wv, sumwu, sumwv, indic4mu,
                indic4nu, indicmat) > 0) & (abs(lo) < 1e05) ) {lo <- lo - dd} 
        }
     if( (flamOld < 0) & (Diff < 0) ) {
            lo <- lamOld-dd
            up <- lamOld
            while ( (lamfun(lo, wu, wv, sumwu, sumwv, indic4mu,
                indic4nu, indicmat) < 0) & (abs(lo) < 1e05) ) {lo <- lo - dd} 	 
        }		
		
     if( (flamOld < 0) & (Diff > 0) ) {
            lo <- lamOld
            up <- lamOld+dd
            while ( (lamfun(up, wu, wv, sumwu, sumwv, indic4mu,
                indic4nu, indicmat) < 0) & (abs(up) < 1e05) ) {up <- up + dd}
        }
			
     if( (flamOld > 0) & (Diff < 0)  )  {
            lo <- lamOld
            up <- lamOld+dd
            while ( (lamfun(up, wu, wv, sumwu, sumwv, indic4mu,
                indic4nu, indicmat) > 0) & (abs(up) < 1e05) ) {up <- up + dd}
          }
	  #Find lam using uniroot.  Changed 12/2021 by M. Zhou
      lam <- uniroot(lamfun, lower = lo, upper = up,
              tol = 1e-09, wu=wu, wv=wv, sumwu=sumwu, sumwv=sumwv,
              indic4mu=indic4mu, indic4nu=indic4nu, indicmat=indicmat)$root   
		 }

#Calculate updated mu1,nu1 using the lagrangian lam
    mu1 <- wu/(sumwu + lam * nu0 %*% t(indicmat))
    nu1 <- wv/(sumwv + lam * mu0 %*% indicmat)

#list the original data & weights plus the p_i, lam0, and mean
    list(u=u, wu=wu, jumpu=mu1, v=v, wv=wv, jumpv=nu1, lam=lam, mean)
}