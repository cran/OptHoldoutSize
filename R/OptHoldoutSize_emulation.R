################################################################################
## R script for functions in OptHoldoutSize package                           ##
################################################################################
##
## Sam Emerson, James Liley
## October 2021
##


######################################################################
## Functions for Bayesian emulation of the learning curve           ##
######################################################################


##' Covariance function for Gaussian process
##'
##'
##' @export
##' @name cov_fn
##' @description Radial kernel covariance function for Gaussian process.
##'
##' Used for a Gaussian process `GP(m,k(.,.))`, an instance X of which has covariance k(n,n') between X(n) and X(n').
##'
##' Covariance function is parametrised by `var_u` (general variance) and `k_width` (kernel width)
##'
##' @keywords emulation
##' @param n Argument 1: kernel is a function of `ndash-n`
##' @param ndash Argument 2: kernel is a function of `ndash-n`
##' @param var_u Global variance
##' @param k_width Kernel width
##' @return Covariance value
##' @examples
##'
##' # We will sample from Gaussian processes
##' #  GP(0,k(.,.)=cov_fn(.,.;var_u,theta))
##' # at these values of n
##' nvals=seq(1,300,length=100)
##'
##' # We will consider two theta values
##' kw1=10; kw2=30
##'
##' # We will consider two var_u values
##' var1=1; var2=10
##'
##' # Covariance matrices
##' cov11=outer(nvals,nvals,function(n,ndash) cov_fn(n,ndash,var_u=var1,
##'   k_width=kw1))
##' cov12=outer(nvals,nvals,function(n,ndash) cov_fn(n,ndash,var_u=var1,
##'   k_width=kw2))
##' cov21=outer(nvals,nvals,function(n,ndash) cov_fn(n,ndash,var_u=var2,
##'   k_width=kw1))
##' cov22=outer(nvals,nvals,function(n,ndash) cov_fn(n,ndash,var_u=var2,
##'   k_width=kw2))
##'
##' # Dampen slightly to ensure positive definiteness
##' damp=1e-5
##' cov11=(1-damp)*(1-diag(length(nvals)))*cov11 + diag(length(nvals))*cov11
##' cov12=(1-damp)*(1-diag(length(nvals)))*cov12 + diag(length(nvals))*cov12
##' cov21=(1-damp)*(1-diag(length(nvals)))*cov21 + diag(length(nvals))*cov21
##' cov22=(1-damp)*(1-diag(length(nvals)))*cov22 + diag(length(nvals))*cov22
##'
##' # Sample
##' set.seed(35243)
##' y11=rmnorm(1,mean=0,varcov=cov11)
##' y12=rmnorm(1,mean=0,varcov=cov12)
##' y21=rmnorm(1,mean=0,varcov=cov21)
##' y22=rmnorm(1,mean=0,varcov=cov22)
##'
##' # Plot
##' rr=max(abs(c(y11,y12,y21,y22)))
##' plot(0,xlim=range(nvals),ylim=c(-rr,rr+10),xlab="n",
##'   ylab=expression("GP(0,cov_fn(.,.;var_u,theta))"))
##' lines(nvals,y11,lty=1,col="black")
##' lines(nvals,y12,lty=2,col="black")
##' lines(nvals,y21,lty=1,col="red")
##' lines(nvals,y22,lty=2,col="red")
##' legend("topright",c("k_width=10, var_u=1", "k_width=30, var_u=1",
##'   "k_width=10, var_u=10","k_width=30, var_u=10"),
##'   lty=c(1,2,1,2),col=c("black","black","red","red"))
##'
cov_fn = function(n,ndash,var_u,k_width){
  return(var_u*exp(-(((n-ndash)/k_width)^2)))
}



##' Updating function for mean.
##'
##' @export
##' @name mu_fn
##' @description Posterior mean for emulator given points `n`.
##' @param n Set of training set sizes to evaluate
##' @param nset Training set sizes for which k2() has been evaluated
##' @param k2 Estimated k2() values at training set sizes `nset`
##' @param var_k2 Variance of error in k2() estimate at each training set size.
##' @param N Total number of samples on which the model will be fitted/used
##' @param k1 Mean cost per sample with no predictive score in place
##' @param var_u Marginal variance for Gaussian process kernel. Defaults to 1e7
##' @param k_width Kernel width for Gaussian process kernel. Defaults to 5000
##' @param k2form Functional form governing expected cost per sample given sample size. Should take two parameters: n (sample size) and theta (parameters). Defaults to function `powerlaw`.
##' @param theta Current estimates of parameter values for k2form. Defaults to the MLE power-law solution corresponding to n,k2, and var_k2.
##' @return Vector Mu of same length of n where Mu_i=mean(posterior(cost(n_i)))
##' @examples
##'
##' # Suppose we have population size and cost-per-sample without a risk score as follows
##' N=100000
##' k1=0.4
##'
##' # Kernel width and variance for GP
##' k_width=5000
##' var_u=8000000
##'
##' # Suppose we begin with k2() estimates at n-values
##' nset=c(10000,20000,30000)
##'
##' # with cost-per-individual estimates
##' # (note that since empirical k2(n) is non-monotonic, it cannot be perfectly
##' #  approximated with a power-law function)
##' k2=c(0.35,0.26,0.28)
##'
##' # and associated error on those estimates
##' var_k2=c(0.02^2,0.01^2,0.03^2)
##'
##' # We estimate theta from these three points
##' theta=powersolve(nset,k2,y_var=var_k2)$par
##'
##' # We will estimate the posterior at these values of n
##' n=seq(1000,50000,length=1000)
##'
##' # Mean and variance
##' p_mu=mu_fn(n,nset=nset,k2=k2,var_k2 = var_k2, N=N,k1=k1,theta=theta,
##'            k_width=k_width,var_u=var_u)
##' p_var=psi_fn(n,nset=nset,N=N,var_k2 = var_k2,k_width=k_width,var_u=var_u)
##'
##' # Plot
##' plot(0,xlim=range(n),ylim=c(20000,60000),type="n",
##'      xlab="Training/holdout set size",
##'      ylab="Total cost (= num. cases)")
##' lines(n,p_mu,col="blue")
##' lines(n,p_mu - 3*sqrt(p_var),col="red")
##' lines(n,p_mu + 3*sqrt(p_var),col="red")
##' points(nset,k1*nset + k2*(N-nset),pch=16,col="purple")
##' lines(n,k1*n + powerlaw(n,theta)*(N-n),lty=2)
##' segments(nset,k1*nset + (k2 - 3*sqrt(var_k2))*(N-nset),
##'          nset,k1*nset + (k2 + 3*sqrt(var_k2))*(N-nset))
##' legend("topright",
##'        c(expression(mu(n)),
##'          expression(mu(n) %+-% 3*sqrt(psi(n))),
##'          "prior(n)",
##'          "d",
##'          "3SD(d|n)"),
##'        lty=c(1,1,2,NA,NA),lwd=c(1,1,1,NA,NA),pch=c(NA,NA,NA,16,124),
##'        pt.cex=c(NA,NA,NA,1,1),
##'        col=c("blue","red","black","purple","black"),bg="white")
mu_fn= function(n,nset,k2,var_k2,N,k1,
  var_u=1e7,
  k_width=5000,
  k2form=powerlaw,
  theta=powersolve(nset,k2,y_var=var_k2)$par){

  cost_fn=function(n,theta) k1*n + k2form(n,theta)*(N-n)

  mean_n <- cost_fn(n,theta) # Expected cost at values n according to current parameters theta
  cov_n_nset <- outer(n,nset,cov_fn,var_u=var_u,k_width=k_width) # GP covariance between nstar values and current values n
  cov_nset_nset <- outer(nset,nset,cov_fn,var_u=var_u,k_width=k_width) # Internal covariance matrix for existing values n
  var_k2_mat <- diag(var_k2*(N-nset)^2) # variance on the scale of cost for whole population

  xd=k1*nset + k2*(N-nset) # outputs on scale of whole-dataset cost

  mean_nset <- cost_fn(nset,theta)
  #browser()
  return(as.numeric(mean_n + cov_n_nset%*%solve(cov_nset_nset + var_k2_mat)%*%(xd - mean_nset)))
}


##' Updating function for variance.
##'
##' @export
##' @name psi_fn
##' @description Posterior variance for emulator given points `n`.
##' @param n Set of training set sizes to evaluate at
##' @param nset Training set sizes for which k2() has been evaluated
##' @param var_k2 Variance of error in k2() estimate at each training set size.
##' @param N Total number of samples on which the model will be fitted/used. Only used to rescale var_k2
##' @param var_u Marginal variance for Gaussian process kernel. Defaults to 1e7
##' @param k_width Kernel width for Gaussian process kernel. Defaults to 5000
##' @return Vector Psi of same length of n where Psi_i=var(posterior(cost(n_i)))
##' @examples
##'
##' # See examples for `mu_fn`
##'
psi_fn =function(n,nset,var_k2,N,
  var_u=1e7,
  k_width=5000){

  cov_n_n <- cov_fn(n,n,var_u=var_u,k_width=k_width)

  cov_n_nset <- outer(n,nset,cov_fn,var_u=var_u,k_width=k_width)
  cov_nset_nset <- outer(nset,nset,cov_fn,var_u=var_u,k_width=k_width)
  var_k2_mat <- diag(var_k2*(N-nset)^2)
  xsol=solve(cov_nset_nset + var_k2_mat)
  return(cov_n_n - apply(cov_n_nset,1,function(x) t(x)%*%xsol%*%x))
}


##' Expected improvement
##'
##' @export
##' @name exp_imp_fn
##' @description Expected improvement
##'
##' Essentially chooses the `next point` to add to `n`, called `n*`, in order to minimise the expectation of `cost(n*)`.
##'
##' @param n Set of training set sizes to evaluate
##' @param nset Training set sizes for which a cost has been evaluated
##' @param k2 Estimates of k2() at training set sizes `nset`
##' @param var_k2 Variance of error in k2() estimates at each training set size.
##' @param N Total number of samples on which the model will be fitted/used
##' @param k1 Mean vost per sample with no predictive score in place
##' @param var_u Marginal variance for Gaussian process kernel. Defaults to 1e7
##' @param k_width Kernel width for Gaussian process kernel. Defaults to 5000
##' @param k2form Functional form governing expected cost per sample given sample size. Should take two parameters: n (sample size) and theta (parameters). Defaults to function `powerlaw`.
##' @param theta Current estimates of parameter values for k2form. Defaults to the MLE power-law solution corresponding to n,k2, and var_k2.
##' @return Value of expected improvement at values n
##' @examples
##'
##' # Set seed.
##' set.seed(24015)
##'
##' # Kernel width and Gaussian process variance
##' kw0=5000
##' vu0=1e7
##'
##' # Include legend on plots or not; inclusion can obscure plot elements on small figures
##' inc_legend=FALSE
##'
##' # Suppose we have population size and cost-per-sample without a risk score as follows
##' N=100000
##' k1=0.4
##'
##' # Suppose that true values of a,b,c are given by
##' theta_true=c(10000,1.2,0.2)
##' theta_lower=c(1,0.5,0.1) # lower bounds for estimating theta
##' theta_upper=c(20000,2,0.5) # upper bounds for estimating theta
##'
##'
##'
##' # We start with five random holdout set sizes (nset0),
##' #  with corresponding cost-per-individual estimates k2_0 derived
##' #  with various errors var_k2_0
##' nstart=4
##' vwmin=0.001; vwmax=0.005
##' nset0=round(runif(nstart,1000,N/2))
##' var_k2_0=runif(nstart,vwmin,vwmax)
##' k2_0=rnorm(nstart,mean=powerlaw(nset0,theta_true),sd=sqrt(var_k2_0))
##'
##' # We estimate theta from these three points
##' theta0=powersolve(nset0,k2_0,y_var=var_k2_0,lower=theta_lower,upper=theta_upper,
##'   init=theta_true)$par
##'
##' # We will estimate the posterior at these values of n
##' n=seq(1000,N,length=1000)
##'
##' # Mean and variance
##' p_mu=mu_fn(n,nset=nset0,k2=k2_0,var_k2 = var_k2_0, N=N,k1=k1,theta=theta0,k_width=kw0,
##'   var_u=vu0)
##' p_var=psi_fn(n,nset=nset0,N=N,var_k2 = var_k2_0,k_width=kw0,var_u=vu0)
##'
##' # Plot
##' yrange=c(-30000,100000)
##' plot(0,xlim=range(n),ylim=yrange,type="n",
##'   xlab="Training/holdout set size",
##'   ylab="Total cost (= num. cases)")
##' lines(n,p_mu,col="blue")
##' lines(n,p_mu - 3*sqrt(p_var),col="red")
##' lines(n,p_mu + 3*sqrt(p_var),col="red")
##' points(nset0,k1*nset0 + k2_0*(N-nset0),pch=16,col="purple")
##' lines(n,k1*n + powerlaw(n,theta0)*(N-n),lty=2)
##' lines(n,k1*n + powerlaw(n,theta_true)*(N-n),lty=3,lwd=3)
##' if (inc_legend) {
##'   legend("topright",
##'     c(expression(mu(n)),
##'       expression(mu(n) %+-% 3*sqrt(psi(n))),
##'       "prior(n)",
##'       "True",
##'       "d"),
##'     lty=c(1,1,2,3,NA),lwd=c(1,1,1,3,NA),pch=c(NA,NA,NA,NA,16),pt.cex=c(NA,NA,NA,NA,1),
##'     col=c("blue","red","black","purple"),bg="white")
##' }
##'
##' ## Add line corresponding to recommended new point
##' exp_imp_em <- exp_imp_fn(n,nset=nset0,k2=k2_0,var_k2 = var_k2_0, N=N,k1=k1,
##'   theta=theta0,k_width=kw0,var_u=vu0)
##' abline(v=n[which.max(exp_imp_em)])
##'
exp_imp_fn = function(n,nset,k2,var_k2,N,k1,
  var_u=1e7,
  k_width=5000,
  k2form=powerlaw,
  theta=powersolve(nset,k2,y_var=var_k2)$par){
  mu_val=mu_fn(n=n,nset=nset,k2=k2,var_k2=var_k2,N=N,k1=k1,var_u=var_u,k_width=k_width,k2form=k2form,theta=theta)
  d_min <- min(k1*nset + k2*(N-nset))
  psi_val <- psi_fn(n=n,nset=nset,var_k2=var_k2,N=N,var_u=var_u,k_width=k_width)
  Z <- (-mu_val+d_min)/sqrt(psi_val)
  EI <- ((-mu_val+d_min)*pnorm(Z)) + (sqrt(psi_val)*dnorm(Z))
  return(EI)
}









##' Measure of error for emulation-based OHS emulation
##'
##' @export
##' @name error_ohs_emulation
##' @description Measure of error for semiparametric (emulation) based estimation of optimal holdout set sizes.
##'
##' Returns a set of values of n for which a `1-alpha` credible interval for cost at includes a lower value than the cost at the estimated optimal holdout size.
##'
##' This is not a confidence interval, credible interval or credible set for the OHS, and is prone to misinterpretation.
##'
##' @keywords estimation emulation
##' @param nset Training set sizes for which k2() has been evaluated
##' @param k2 Estimated k2() at training set sizes `nset`
##' @param var_k2 Variance of error in k2() estimate at each training set size.
##' @param N Total number of samples on which the model will be fitted/used
##' @param k1 Mean cost per sample with no predictive score in place
##' @param alpha Use 1-alpha credible interval. Defaults to 0.1.
##' @param var_u Marginal variance for Gaussian process kernel. Defaults to 1e7
##' @param k_width Kernel width for Gaussian process kernel. Defaults to 5000
##' @param k2form Functional form governing expected cost per sample given sample size. Should take two parameters: n (sample size) and theta (parameters). Defaults to function `powerlaw`.
##' @param theta Current estimates of parameter values for k2form. Defaults to the MLE power-law solution corresponding to n,k2, and var_k2.
##' @param npoll Check npoll equally spaced values between 1 and N for minimum. If NULL, check all values (this can be slow). Defaults to 1000
##' @return Vector of values `n` for which 1-alpha credible interval for cost `l(n)` at n contains mean posterior cost at estimated optimal holdout size.
##' @examples
##'
##'  # Set seed
##' set.seed(57365)
##'
##' # Parameters
##' N=100000;
##' k1=0.3
##' A=8000; B=1.5; C=0.15; theta=c(A,B,C)
##'
##' # True mean function
##' k2_true=function(n) powerlaw(n,theta)
##'
##' # True OHS
##' nx=1:N
##' ohs_true=nx[which.min(k1*nx + k2_true(nx)*(N-nx))]
##'
##' # Values of n for which cost has been estimated
##' np=50 # this many points
##' nset=round(runif(np,1,N))
##' var_k2=runif(np,0.001,0.0015)
##' k2=rnorm(np,mean=k2_true(nset),sd=sqrt(var_k2))
##'
##' # Compute OHS
##' res1=optimal_holdout_size_emulation(nset,k2,var_k2,N,k1)
##'
##' # Error estimates
##' ex=error_ohs_emulation(nset,k2,var_k2,N,k1)
##'
##' # Plot
##' plot(res1)
##'
##' # Add error
##' abline(v=ohs_true)
##' abline(v=ex,col=rgb(1,0,0,alpha=0.2))
##'
##' # Show justification for error
##' n=seq(1,N,length=1000)
##' mu=mu_fn(n,nset,k2,var_k2,N,k1); psi=pmax(0,psi_fn(n, nset, var_k2, N)); Z=-qnorm(0.1/2)
##' lines(n,mu - Z*sqrt(psi),lty=2,lwd=2)
##' legend("topright",
##'     c("Err. region",expression(paste(mu(n)- "z"[alpha/2]*sqrt(psi(n))))),
##'     pch=c(16,NA),lty=c(NA,2),lwd=c(NA,2),col=c("pink","black"),bty="n")
error_ohs_emulation=function(nset,k2,var_k2,N,k1,
                             alpha=0.1,
                             var_u=1e7,
                             k_width=5000,
                             k2form=powerlaw,
                             theta=powersolve(nset,k2,y_var=var_k2)$par,
                             npoll=1000){

  # Candidate values n
  if (!is.null(npoll)) n=seq(1,N,length=npoll) else n=1:N

  # mu and psi
  xmu=mu_fn(n=n,nset=nset,k2=k2,var_k2=var_k2,N=N,k1=k1,var_u=var_u,
            k_width=k_width,k2form=k2form,theta=theta)
  xpsi=pmax(0,psi_fn(n=n, nset=nset, var_k2=var_k2, N=N, var_u=var_u,
                     k_width=k_width))

  # Compute minimum
  w=which.min(xmu)
  ohs=n[w]
  est_min=xmu[w]

  z=-qnorm(alpha/2) # Need to be this many standard deviations below the mean

  # Values of n for which 1-alpha credible interval for l(n) includes est_min
  n_cont=n[which(xmu - z*sqrt(xpsi) <= est_min)]

  return(n_cont)
}
