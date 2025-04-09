################################################################################
## R script for functions in OptHoldoutSize package                           ##
################################################################################
##
## Sami Haidar-Wehbe, Sam Emerson, James Liley
## October 2021
##



################################################################################
## Optimal holdout size, confidence interval, and gradient                    ##
################################################################################

##' Estimate optimal holdout size under parametric assumptions
##'
##'
##' @export
##' @name optimal_holdout_size
##' @description Compute optimal holdout size for updating a predictive score given appropriate parameters of cost function
##'
##' Evaluates empirical minimisation of cost function ``l(n;k1,N,theta) = k1 n + k2form(n;theta) (N-n)``.
##'
##' The function will return `Inf` if no minimum exists. It does not check if the minimum is unique, but this can be guaranteed using the assumptions for theorem 1 in the manuscript.
##'
##' This calls the function `optimize` from package `stats`.
##' @keywords estimation
##' @param N Total number of samples on which the predictive score will be used/fitted. Can be a vector.
##' @param k1 Cost value in the absence of a predictive score. Can be a vector.
##' @param theta Parameters for function k2form(n) governing expected cost to an individual sample given a predictive score fitted to n samples. Can be a matrix of dimension n x n_par, where n_par is the number of parameters of k2.
##' @param k2form Function governing expected cost to an individual sample given a predictive score fitted to n samples. Must take two arguments: n (number of samples) and theta (parameters). Defaults to a power-law form ``powerlaw(n,c(a,b,c))=a n^(-b) + c``.
##' @param round_result Set to TRUE to solve over integral sizes
##' @param ... Passed to function `optimize`
##' @return List/data frame of dimension (number of evaluations) x (4 + n_par) containing input data and results. Columns size and cost are optimal holdout size and cost at this size respectively. Parameters N, k1, theta.1, theta.2,...,theta.n_par are input data.
##' @examples
##'
##' # Evaluate optimal holdout set size for a range of values of k1 and two values of
##' #   N, some of which lead to infinite values
##' N1=10000; N2=12000
##' k1=seq(0.1,0.5,length=20)
##' A=3; B=1.5; C=0.15; theta=c(A,B,C)
##'
##' res1=optimal_holdout_size(N1,k1,theta)
##' res2=optimal_holdout_size(N2,k1,theta)
##'
##' oldpar=par(mfrow=c(1,2))
##' plot(0,type="n",ylim=c(0,500),xlim=range(res1$k1),xlab=expression("k"[1]),
##'   ylab="Optimal holdout set size")
##'   lines(res1$k1,res1$size,col="black")
##'   lines(res2$k1,res2$size,col="red")
##'   legend("topright",as.character(c(N1,N2)),title="N:",col=c("black","red"),lty=1)
##' plot(0,type="n",ylim=c(1500,1600),xlim=range(res1$k1),xlab=expression("k"[1]),
##'   ylab="Minimum cost")
##'   lines(res1$k1,res1$cost,col="black")
##'   lines(res2$k1,res2$cost,col="red")
##'   legend("topleft",as.character(c(N1,N2)),title="N:",col=c("black","red"),lty=1)
##'
##' par(oldpar)
optimal_holdout_size=function(
  N,
  k1,
  theta,
  k2form = powerlaw,
  round_result=FALSE,
  ...
) {

  ## Error handlers
  if (!is.numeric(c(N,k1))) stop("Parameters N and k1 must be numeric")
  if (!(is.numeric(theta)|is.matrix(theta))) stop("Parameter theta must be a vector or matrix")
  if (!is.function(k2form)) stop("Parameter k2form must be a function taking two arguments: n and theta")
  if (length(as.list(args(k2form)))!=3) stop("Parameter k2form must be a function taking two arguments: n and theta")
  if ((length(N)>1 | length(k1)>1)|!is.null(dim(theta))) {
    n=max(length(N),length(k1),dim(theta)[1])
    if (is.null(dim(theta))) theta=t(matrix(theta,length(theta),n))
    if (!all(c(length(k1),length(N),dim(theta)[1]) %in% c(1,n))) stop("If vectors, N and k1 must have the same length, or have length one. Parameter theta must either be a vector or a matrix of dimension length(n) x n_par, where n_par is the number of parameters of k2.")
    if (length(N)==1) N=rep(N,n)
    if (length(k1)==1) k1=rep(k1,n)
  } else theta=matrix(theta,1,length(theta))

  par_mat=cbind(N,k1,theta)
  np=dim(par_mat)[2]

    ## Solve
  out=suppressWarnings(apply(par_mat,1,function(x) {
    unlist(optimize(function(n) x[2]*n + k2form(n,x[3:np])*(x[1]-n),c(1,N),...))
  }))

  ## Post-processing
  out=cbind(t(out),par_mat)
  out=as.data.frame(out)
  colnames(out)=c("size","cost","N","k1",paste0("theta.",1:dim(theta)[2]))

  ## Replace errors with infinity
  no_min=which((out$size> (N-1))|(out$size < 2))
  out$size[no_min]=Inf
  out$cost[no_min]=Inf

  # Round result if required. Minimum will be at either floor(size) or ceiling(size)
  if (round_result) {
    size1=floor(out$size)
    size2=ceiling(out$size)
    k2_1=apply(cbind(size1,theta),1,function(x) k2form(x[1],x[2:(1+dim(theta)[2])]))
    k2_2=apply(cbind(size2,theta),1,function(x) k2form(x[1],x[2:(1+dim(theta)[2])]))
    cost1=out$k1*size1 + k2_1*(out$N-size1)
    cost2=out$k1*size2 + k2_1*(out$N-size2)

    w=which(cost2<cost1)
    out_size=size1; out_cost=cost1
    out_size[w]=size2[w]; out_cost[w]=cost2[w]

    out$size=out_size
    out$cost=out_cost
  }

  class(out)=c("optholdoutsize",class(out))
  return(out)
}

##' Plot estimated cost function
##'
##' @export
##' @name plot.optholdoutsize
##' @keywords estimation
##' @description Plot estimated cost function, when parametric method is used for estimation.
##'
##' Draws cost function as a line and indicates minimum. Assumes a power-law form of k2 unless parameter k2 is set otherwise.
##'
##' @param x Object of type `optholdoutsize`
##' @param ... Other arguments passed to `plot()` and `lines()`
##' @param k2form Function governing expected cost to an individual sample given a predictive score fitted to n samples. Must take two arguments: n (number of samples) and theta (parameters). Defaults to a power-law form ``powerlaw(n,c(a,b,c))=a n^(-b) + c``.
##' @return No return value; draws a plot only.
##' @examples
##'
##' # Simple example
##'
##' N=100000;
##' k1=0.3
##' A=8000; B=1.5; C=0.15; theta=c(A,B,C)
##'
##' res1=optimal_holdout_size(N,k1,theta)
##'
##' plot(res1)
##'
plot.optholdoutsize=function(x,...,k2form=powerlaw) {

  N=x$N; k1=x$k1; size=x$size; cost=x$cost;
  if (length(N)==1) {
    theta=unlist(x[5:length(x)])
    xx=seq(1,N,length.out=1000)
    yy= k1*xx + k2form(xx,theta)*(N-xx)
    costN=yy[length(xx)]

    plot(xx,yy,type="l",ylim=c(2*cost-costN,costN + cost),xlab="Holdout size",ylab="Total cost",...)
    abline(v=size,lty=2); abline(h=cost,lty=2)
    points(size,cost,pch=16,col="black")
    legend("bottomright",c("Est. cost","OHS"),lty=c(1,NA),pch=c(NA,16),col=c("black","black"),bty="n")
  } else {
    Nmax=max(N)
    theta=as.data.frame(x[5:length(x)])

    xx=seq(1,Nmax,length.out=1000)
    yy=matrix(0,length(N),length(xx))
    for (i in 1:length(N)) yy[i,]=k1[i]*xx + k2form(xx,as.numeric(theta[i,]))*(N[i]-xx)

    costN=yy[,length(xx)]

    plot(0,type="n",xlim=c(1,Nmax),ylim=c(2*min(cost)-mean(costN),max(costN) + min(cost)),
      xlab="Holdout size",ylab="Total cost",...)
    for (i in 1:length(N)) lines(xx,yy[i,],...)
    points(size,cost,pch=16)
    legend("bottomright",c("Est. cost","OHS"),lty=c(1,NA),pch=c(NA,16),col=c("black","black"),bty="n")
  }
}





##' Estimate optimal holdout size under semi-parametric assumptions
##'
##'
##' @export
##' @name optimal_holdout_size_emulation
##' @keywords estimation emulation
##' @description Compute optimal holdout size for updating a predictive score given a set of training set sizes and estimates of mean cost per sample at those training set sizes.
##'
##' This is essentially a wrapper for function `mu_fn()`.
##'
##'
##' @param nset Training set sizes for which a cost has been evaluated
##' @param k2 Estimated values of k2() at training set sizes `nset`
##' @param var_k2 Variance of error in k2 estimate at each training set size.
##' @param N Total number of samples on which the model will be fitted/used
##' @param k1 Mean cost per sample with no predictive score in place
##' @param var_u Marginal variance for Gaussian process kernel. Defaults to 1e7
##' @param k_width Kernel width for Gaussian process kernel. Defaults to 5000
##' @param k2form Functional form governing expected cost per sample given sample size. Should take two parameters: n (sample size) and theta (parameters). Defaults to function `powerlaw`.
##' @param theta Current estimates of parameter values for k2form. Defaults to the MLE power-law solution corresponding to n,k2, and var_k2.
##' @param npoll Check npoll equally spaced values between 1 and N for minimum. If NULL, check all values (this can be slow). Defaults to 1000
##' @param ... Passed to function `optimise()`
##' @return Object of class 'optholdoutsize_emul' with elements "cost" (minimum cost),"size" (OHS),"nset","k2","var_k2","N","k1","var_u","k_width","theta" (parameters)
##' @examples
##'
##' # See examples for mu_fn()
optimal_holdout_size_emulation= function(nset,k2,var_k2,N,k1,
  var_u=1e7,
  k_width=5000,
  k2form=powerlaw,
  theta=powersolve_general(nset,k2,y_var=var_k2)$par,
  npoll=1000,
  ...){

  if (!is.null(npoll)) n=seq(1,N,length=npoll) else n=1:N
  xv=mu_fn(n=n,nset=nset,k2=k2,var_k2=var_k2,N=N,k1=k1,
           var_u=var_u,k_width=k_width,k2form=k2form,
           theta=theta)
  w=which.min(xv)

  out=list(xv[w],n[w],nset,k2,var_k2,N,k1,var_u,k_width,theta)
  ## The 'optimise' below doesn't tend to work bimodally
  #minf=function(n) mu_fn(n,nset,k2,var_k2,N,k1,var_u,k_width,mean_fn,theta)
  #out=optimise(minf,c(1,N))

  names(out)=c("cost","size","nset","k2","var_k2","N","k1","var_u","k_width","theta")
  class(out)=c("optholdoutsize_emul",class(out))
  return(out)
}

##' Plot estimated cost function using emulation (semiparametric)
##'
##' @export
##' @name plot.optholdoutsize_emul
##' @keywords estimation
##' @description Plot estimated cost function, when semiparametric (emulation) method is used for estimation.
##'
##' Draws posterior mean of cost function as a line and indicates minimum. Also draws mean +/- 3 SE.
##'
##' Assumes a power-law form of k2 unless parameter k2 is set otherwise.
##'
##' @param x Object of type `optholdoutsize_emul`
##' @param ... Other arguments passed to `plot()`
##' @param k2form Function governing expected cost to an individual sample given a predictive score fitted to n samples. Must take two arguments: n (number of samples) and theta (parameters). Defaults to a power-law form ``powerlaw(n,c(a,b,c))=a n^(-b) + c``.
##' @return No return value; draws a plot only.
##' @examples
##'
##' # Simple example
##'
##' # Parameters
##' N=100000;
##' k1=0.3
##' A=8000; B=1.5; C=0.15; theta=c(A,B,C)
##'
##' # True mean function
##' k2_true=function(n) powerlaw(n,theta)
##'
##' # Values of n for which cost has been estimated
##' np=50 # this many points
##' nset=round(runif(np,1,N))
##' var_k2=runif(np,0.001,0.002)
##' k2=rnorm(np,mean=k2_true(nset),sd=sqrt(var_k2))
##'
##' # Compute OHS
##' res1=optimal_holdout_size_emulation(nset,k2,var_k2,N,k1)
##'
##' # Plot
##' plot(res1)
plot.optholdoutsize_emul=function(x,...,k2form=powerlaw) {

  # Plot at these values of n
  nn=seq(1,x$N,length.out=1000)

  # Mean function
  xv=mu_fn(n=nn,nset=x$nset,k2=x$k2,var_k2=x$var_k2,
           N=x$N,k1=x$k1,var_u=x$var_u,k_width=x$k_width,
           k2form=k2form,x$theta)

  # Variance function
  psiv=pmax(0,psi_fn(n=nn, nset=x$nset, var_k2=x$var_k2,
                     N=x$N, var_u=x$var_u, k_width=x$k_width))

  # 3 SD bounds
  xlo=xv - 3*sqrt(psiv)
  xhi=xv + 3*sqrt(psiv)

  # m(n)
  xm=x$k1*nn + k2form(nn,x$theta)*(x$N-nn)

  xq=quantile(c(xlo,xhi),c(0.1,0.9)); xmid=median(xv)
  yr=c(xmid - 2*(xmid-xq[1]),xmid + 2*(xq[2]-xmid))

  plot(0,type="n",xlab="Holdout size",ylab="Total cost",xlim=c(1,x$N),ylim=yr,...)
  lines(nn,xm,lty=2,col="black")
  lines(nn,xv,col="blue")
  lines(nn,xlo,col="red")
  lines(nn,xhi,col="red")
  points(x$nset,x$k1*x$nset + x$k2*(x$N - x$nset),col="purple",pch=16)
  abline(h=x$cost,lty=2); abline(v=x$size,lty=2)
  legend("bottomright",
    c(expression(mu(n)),
      expression(mu(n) %+-% 3*sqrt(psi(n))),
      "m(n)",
      "d"),
    lty=c(1,1,2,NA),lwd=c(1,1,1,NA),pch=c(NA,NA,NA,16),pt.cex=c(NA,NA,NA,1),
    col=c("blue","red","black","purple"),bg="white",bty="n")
}




##' Confidence interval for optimal holdout size, when estimated using parametric method
##'
##' @export
##' @name ci_ohs
##' @description Compute confidence interval for optimal holdout size given either a standard error covariance matrix or a set of m estimates of parameters.
##'
##' This can be done either asymptotically, using a method analogous to the Fisher information matrix, or empirically (using bootstrap resampling)
##'
##' If sigma (covariance matrix) is specified and method='bootstrap', a confidence interval is generated assuming a Gaussian distribution of (N,k1,theta). To estimate a confidence interval assuming a non-Gaussian distribution, simulate values under the requisite distribution and use then as parameters N,k1, theta, with sigma set to NULL.
##'
##' @keywords estimation
##' @param N Vector of estimates of total number of samples on which the predictive score will be used/fitted, or single estimate
##' @param k1 Vector of estimates of cost value in the absence of a predictive score, or single number
##' @param theta Matrix of estimates of parameters for function k2form(n) governing expected cost to an individual sample given a predictive score fitted to n samples. Can be a matrix of dimension n x n_par, where n_par is the number of parameters of k2.
##' @param alpha Construct 1-alpha confidence interval. Defaults to 0.05
##' @param k2form Function governing expected cost to an individual sample given a predictive score fitted to n samples. Must take two arguments: n (number of samples) and theta (parameters). Defaults to a power-law form ``k2(n,c(a,b,c))=a n^(-b) + c``.
##' @param grad_nstar Function giving partial derivatives of optimal holdout set, taking three arguments: N, k1, and theta. Only used for asymptotic confidence intervals. F NULL, estimated empirically
##' @param sigma Standard error covariance matrix for (N,k1,theta), in that order. If NULL, will derive as sample covariance matrix of parameters. Must be of the correct size and positive definite.
##' @param n_boot Number of bootstrap resamples for empirical estimate.
##' @param seed Random seed for bootstrap resamples. Defaults to NULL.
##' @param mode One of 'asymptotic' or 'empirical'. Defaults to 'empirical'
##' @param ... Passed to function `optimize`
##' @return A vector of length two containing lower and upper limits of confidence interval.
##' @examples
##'
##' ## Set seed
##' set.seed(493825)
##'
##' ## We will assume that our observations of N, k1, and theta=(a,b,c) are
##' ##  distributed with mean mu_par and variance sigma_par
##' mu_par=c(N=10000,k1=0.35,A=3,B=1.5,C=0.1)
##' sigma_par=cbind(
##'   c(100^2,       1,      0,       0,       0),
##'   c(    1,  0.07^2,      0,       0,       0),
##'   c(    0,       0,  0.5^2,    0.05,  -0.001),
##'   c(    0,       0,   0.05,   0.4^2,  -0.002),
##'   c(    0,       0, -0.001,  -0.002,  0.02^2)
##' )
##'
##' # Firstly, we make 500 observations
##' par_obs=rmnorm(500,mean=mu_par,varcov=sigma_par)
##'
##' # Optimal holdout size and asymptotic and empirical confidence intervals
##' ohs=optimal_holdout_size(N=mean(par_obs[,1]),k1=mean(par_obs[,2]),
##'   theta=colMeans(par_obs[,3:5]))$size
##' ci_a=ci_ohs(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],alpha=0.05,
##'   seed=12345,mode="asymptotic")
##' ci_e=ci_ohs(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],alpha=0.05,
##'   seed=12345,mode="empirical")
##'
##'
##' # Assess cover at various m
##' m_values=c(20,30,50,100,150,200,300,500,750,1000,1500)
##' ntrial=5000
##' alpha_trial=0.1 # use 90% confidence intervals
##' nstar_true=optimal_holdout_size(N=mu_par[1],k1=mu_par[2],
##'   theta=mu_par[3:5])$size
##'
##' ## The matrices indicating cover take are included in this package but take
##' ##  around 30 minutes to generate. They are generated using the code below
##' ##  (including random seeds).
##' data(ci_cover_a_yn)
##' data(ci_cover_e_yn)
##'
##' if (!exists("ci_cover_a_yn")) {
##'   ci_cover_a_yn=matrix(NA,length(m_values),ntrial) # Entry [i,j] is 1 if ith
##'   ##  asymptotic CI for jth value of m covers true nstar
##'   ci_cover_e_yn=matrix(NA,length(m_values),ntrial) # Entry [i,j] is 1 if ith
##'   ##  empirical CI for jth value of m covers true nstar
##'
##'   for (i in 1:length(m_values)) {
##'     m=m_values[i]
##'     for (j in 1:ntrial) {
##'       # Set seed
##'       set.seed(j*ntrial + i + 12345)
##'
##'       # Make m observations
##'       par_obs=rmnorm(m,mean=mu_par,varcov=sigma_par)
##'       ci_a=ci_ohs(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],
##'         alpha=alpha_trial,mode="asymptotic")
##'       ci_e=ci_ohs(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],
##'         alpha=alpha_trial,mode="empirical",n_boot=500)
##'
##'       if (nstar_true>ci_a[1] & nstar_true<ci_a[2]) ci_cover_a_yn[i,j]=1 else
##'          ci_cover_a_yn[i,j]=0
##'       if (nstar_true>ci_e[1] & nstar_true<ci_e[2]) ci_cover_e_yn[i,j]=1 else
##'          ci_cover_e_yn[i,j]=0
##'     }
##'     print(paste0("Completed for m = ",m))
##'   }
##'
##' save(ci_cover_a_yn,file="data/ci_cover_a_yn.RData")
##' save(ci_cover_e_yn,file="data/ci_cover_e_yn.RData")
##'
##' }
##'
##' # Cover at each m value and standard error
##' cover_a=rowMeans(ci_cover_a_yn)
##' cover_e=rowMeans(ci_cover_e_yn)
##' zse_a=2*sqrt(cover_a*(1-cover_a)/ntrial)
##' zse_e=2*sqrt(cover_e*(1-cover_e)/ntrial)
##'
##'
##' # Draw plot. Convergence to 1-alpha cover is evident. Cover is not far from
##' #   alpha even at small m.
##'
##' plot(0,type="n",xlim=range(m_values),ylim=c(0.7,1),xlab=expression("m"),
##'   ylab="Cover")
##'
##' # Asymptotic cover and 2*SE pointwise envelope
##' polygon(c(m_values,rev(m_values)),c(cover_a+zse_a,rev(cover_a-zse_a)),
##'   col=rgb(0,0,0,alpha=0.3),border=NA)
##' lines(m_values,cover_a,col="black")
##'
##' # Empirical cover and 2*SE pointwiseenvelope
##' polygon(c(m_values,rev(m_values)),c(cover_e+zse_e,rev(cover_e-zse_e)),
##'   col=rgb(0,0,1,alpha=0.3),border=NA)
##' lines(m_values,cover_e,col="blue")
##'
##' abline(h=1-alpha_trial,col="red")
##' legend("bottomright",c("Asym.","Emp.",expression(paste("1-",alpha))),lty=1,
##'   col=c("black","blue","red"))
##'
ci_ohs=function(
  N,
  k1,
  theta,
  alpha=0.05,
  k2form = powerlaw,
  grad_nstar=NULL,
  sigma=NULL,
  n_boot=1000,
  seed=NULL,
  mode="empirical",
  ...
) {

  ## Error handlers
  if (length(N)<5 & is.null(sigma)) stop("At least five samples necessary for estimation if parameter sigma not specified")
  if (!(is.numeric(N) & is.numeric(k1) & is.numeric(theta))) stop("Parameters N, k1 and theta must be numeric")
  if (!is.function(k2form)) stop("Parameter k2form must be a function taking two arguments: n and theta")
  if (length(as.list(args(k2form)))!=3) stop("Parameter k2 must be a function taking two arguments: n and theta")
  if (!is.null(grad_nstar)) {
    if (!is.function(grad_nstar)) stop("Parameter grad_nstar must be a function taking three arguments: N, k1, and theta")
    if (length(as.list(args(grad_nstar)))!=4) stop("Parameter grad_nstar must be a function taking three arguments: N, k1, and theta")
  }
  if (!(mode %in% c("asymptotic","empirical"))) stop("Parameter mode must be either 'asymptotic' or 'empirical'")
  if (mode=="empirical") if (n_boot*alpha < 10) stop("Parameter nboot too low for this level of alpha")
  if (!is.numeric(alpha)) stop("Parameter alpha must be numeric")
  if (alpha <= 0 | alpha>=1) stop("Parameter alpha must be >0 and <1")

  if (mode=="asymptotic") {
    if (is.null(sigma)) {
      par_mat=cbind(N,k1,theta)
      mu=colMeans(par_mat)
      sigma_hat=var(par_mat)
    } else {
      mu=c(N,k1,theta)
      sigma_hat=sigma
    }

    ## Determine function for grad_nstar if not specified
    if (is.null(grad_nstar)) {
      grad_nstar=function(N,k1,theta) {
        dx=1e-5 # estimate gradient by using this difference
        par2=outer(rep(1,2+length(theta)),c(N,k1,theta))
        par2=par2 + dx*diag(dim(par2)[1]) # parameters, shifted by dx one-at-a-time
        ohs2=optimal_holdout_size(N=par2[,1],k1=par2[,2],theta=par2[,3:dim(par2)[2]],k2form=k2form)$size
        ohs1=optimal_holdout_size(N,k1,theta,k2form=k2form)$size
        return(t((ohs2-ohs1)/dx))
      }
    }

    # Parameters for asymptotic confidence interval
    nstar=optimal_holdout_size(mu[1],mu[2],mu[3:length(mu)],k2form=k2form)$size
    beta_est=grad_nstar(mu[1],mu[2],mu[3:length(mu)])
    z_a=-qnorm(alpha/2)
    m=length(N)
    # Variation from estimated value in asymptotic confidence interval
    di=z_a*sqrt(((beta_est) %*% sigma_hat %*% t(beta_est))/m)

    cx=nstar+c(-di,di)
    names(cx)=c("lower","upper")
    return(cx)
  }

  if (mode=="empirical") {

    if (is.null(sigma)) {
      par_mat=cbind(N,k1,theta)
      mu=colMeans(par_mat)
      sigma_hat=var(par_mat)
      m=dim(par_mat)[1]
    } else {
      mu=c(N,k1,theta)
      sigma_hat=sigma
      m=n_boot

      # Sample, allowing that some parameters may have zero variance
      w=which(sigma_hat[cbind(1:length(mu),1:length(mu))]<1e-20)
      wc=setdiff(1:length(mu),w)
      if (length(w)>0) {
        psub=outer(rep(1,n_boot),mu[w])
        psubc=rmnorm(n_boot,mean=mu[wc],varcov=sigma_hat[wc,wc])
        par_mat=outer(rep(1,n_boot),rep(0,length(mu)))
        par_mat[,w]=psub; par_mat[,wc]=psubc
      } else par_mat=rmnorm(n_boot,mean=mu,varcov=sigma_hat)
    }

    # Random seed
    if (!is.null(seed)) set.seed(seed)

    # Populate with mean parameters from bootstrap resamples
    ci_mat=matrix(0,n_boot,dim(par_mat)[2])
    for (i in 1:n_boot) {
      sboot=sample(m,m,replace=TRUE)
      ci_mat[i,]=colMeans(par_mat[sboot,])
    }

    # Compute OHSs for bootstrap resamples
    ohs_boot=optimal_holdout_size(ci_mat[,1],ci_mat[,2],ci_mat[,3:dim(ci_mat)[2]],k2form=k2form)

    # Estimate confidence interval as quantile
    cx=quantile(ohs_boot$size,c(alpha/2,1-(alpha/2)))
    names(cx)=c("lower","upper")
    return(cx)
  }
}








##' Confidence interval for minimum total cost, when estimated using parametric method
##'
##' @export
##' @name ci_mincost
##' @description Compute confidence interval for cost at optimal holdout size given either a standard error covariance matrix or a set of m estimates of parameters.
##'
##' This can be done either asymptotically, using a method analogous to the Fisher information matrix, or empirically (using bootstrap resampling)
##'
##' If sigma (covariance matrix) is specified and method='bootstrap', a confidence interval is generated assuming a Gaussian distribution of (N,k1,theta). To estimate a confidence interval assuming a non-Gaussian distribution, simulate values under the requisite distribution and use then as parameters N,k1, theta, with sigma set to NULL.
##'
##' @keywords estimation
##' @param N Vector of estimates of total number of samples on which the predictive score will be used/fitted, or single estimate
##' @param k1 Vector of estimates of cost value in the absence of a predictive score, or single number
##' @param theta Matrix of estimates of parameters for function k2form(n) governing expected cost to an individual sample given a predictive score fitted to n samples. Can be a matrix of dimension n x n_par, where n_par is the number of parameters of k2.
##' @param alpha Construct 1-alpha confidence interval. Defaults to 0.05
##' @param k2form Function governing expected cost to an individual sample given a predictive score fitted to n samples. Must take two arguments: n (number of samples) and theta (parameters). Defaults to a power-law form ``k2(n,c(a,b,c))=a n^(-b) + c``.
##' @param grad_mincost Function giving partial derivatives of minimum cost, taking three arguments: N, k1, and theta. Only used for asymptotic confidence intervals. F NULL, estimated empirically
##' @param sigma Standard error covariance matrix for (N,k1,theta), in that order. If NULL, will derive as sample covariance matrix of parameters. Must be of the correct size and positive definite.
##' @param n_boot Number of bootstrap resamples for empirical estimate.
##' @param seed Random seed for bootstrap resamples. Defaults to NULL.
##' @param mode One of 'asymptotic' or 'empirical'. Defaults to 'empirical'
##' @param ... Passed to function `optimize`
##' @return A vector of length two containing lower and upper limits of confidence interval.
##' @examples
##'
##' ## Set seed
##' set.seed(574635)
##'
##' ## We will assume that our observations of N, k1, and theta=(a,b,c) are
##' ##  distributed with mean mu_par and variance sigma_par
##' mu_par=c(N=10000,k1=0.35,A=3,B=1.5,C=0.1)
##' sigma_par=cbind(
##'   c(100^2,       1,      0,       0,       0),
##'   c(    1,  0.07^2,      0,       0,       0),
##'   c(    0,       0,  0.5^2,    0.05,  -0.001),
##'   c(    0,       0,   0.05,   0.4^2,  -0.002),
##'   c(    0,       0, -0.001,  -0.002,  0.02^2)
##' )
##'
##' # Firstly, we make 500 observations
##' par_obs=rmnorm(500,mean=mu_par,varcov=sigma_par)
##'
##' # Minimum cost and asymptotic and empirical confidence intervals
##' mincost=optimal_holdout_size(N=mean(par_obs[,1]),k1=mean(par_obs[,2]),
##'   theta=colMeans(par_obs[,3:5]))$cost
##' ci_a=ci_mincost(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],alpha=0.05,
##'   seed=12345,mode="asymptotic")
##' ci_e=ci_mincost(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],alpha=0.05,
##'   seed=12345,mode="empirical")
##'
##'
##' # Assess cover at various m
##' m_values=c(20,30,50,100,150,200,300,500,750,1000,1500)
##' ntrial=5000
##' alpha_trial=0.1 # use 90% confidence intervals
##' mincost_true=optimal_holdout_size(N=mu_par[1],k1=mu_par[2],
##'   theta=mu_par[3:5])$cost
##'
##' ## The matrices indicating cover take are included in this package but take
##' ##  around 30 minutes to generate. They are generated using the code below
##' ##  (including random seeds).
##' data(ci_cover_cost_a_yn)
##' data(ci_cover_cost_e_yn)
##'
##' if (!exists("ci_cover_cost_a_yn")) {
##'   ci_cover_cost_a_yn=matrix(NA,length(m_values),ntrial) # Entry [i,j] is 1
##'   #  if ith asymptotic CI for jth value of m covers true mincost
##'   ci_cover_cost_e_yn=matrix(NA,length(m_values),ntrial) # Entry [i,j] is 1
##'   #  if ith empirical CI for jth value of m covers true mincost
##'
##'   for (i in 1:length(m_values)) {
##'     m=m_values[i]
##'     for (j in 1:ntrial) {
##'       # Set seed
##'       set.seed(j*ntrial + i + 12345)
##'
##'       # Make m observations
##'       par_obs=rmnorm(m,mean=mu_par,varcov=sigma_par)
##'       ci_a=ci_mincost(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],
##'         alpha=alpha_trial,mode="asymptotic")
##'       ci_e=ci_mincost(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],
##'         alpha=alpha_trial,mode="empirical",n_boot=500)
##'
##'       if (mincost_true>ci_a[1] & mincost_true<ci_a[2])
##'         ci_cover_cost_a_yn[i,j]=1 else ci_cover_cost_a_yn[i,j]=0
##'       if (mincost_true>ci_e[1] & mincost_true<ci_e[2])
##'         ci_cover_cost_e_yn[i,j]=1 else ci_cover_cost_e_yn[i,j]=0
##'     }
##'     print(paste0("Completed for m = ",m))
##'   }
##'
##' save(ci_cover_cost_a_yn,file="data/ci_cover_cost_a_yn.RData")
##' save(ci_cover_cost_e_yn,file="data/ci_cover_cost_e_yn.RData")
##'
##' }
##'
##' # Cover at each m value and standard error
##' cover_a=rowMeans(ci_cover_cost_a_yn)
##' cover_e=rowMeans(ci_cover_cost_e_yn)
##' zse_a=2*sqrt(cover_a*(1-cover_a)/ntrial)
##' zse_e=2*sqrt(cover_e*(1-cover_e)/ntrial)
##'
##'
##' # Draw plot. Convergence to 1-alpha cover is evident. Cover is not far from
##' #   alpha even at small m.
##'
##' plot(0,type="n",xlim=range(m_values),ylim=c(0.7,1),xlab=expression("m"),
##'   ylab="Cover")
##'
##' # Asymptotic cover and 2*SE pointwise envelope
##' polygon(c(m_values,rev(m_values)),c(cover_a+zse_a,rev(cover_a-zse_a)),
##'   col=rgb(0,0,0,alpha=0.3),border=NA)
##' lines(m_values,cover_a,col="black")
##'
##' # Empirical cover and 2*SE pointwiseenvelope
##' polygon(c(m_values,rev(m_values)),c(cover_e+zse_e,rev(cover_e-zse_e)),
##'   col=rgb(0,0,1,alpha=0.3),border=NA)
##' lines(m_values,cover_e,col="blue")
##'
##' abline(h=1-alpha_trial,col="red")
##' legend("bottomright",c("Asym.","Emp.",expression(paste("1-",alpha))),lty=1,
##'   col=c("black","blue","red"))
##'
ci_mincost=function(
  N,
  k1,
  theta,
  alpha=0.05,
  k2form = powerlaw,
  grad_mincost=NULL,
  sigma=NULL,
  n_boot=1000,
  seed=NULL,
  mode="empirical",
  ...
) {

  ## Error handlers
  if (length(N)<5 & is.null(sigma)) stop("At least five samples necessary for estimation if parameter sigma not specified")
  if (!(is.numeric(N) & is.numeric(k1) & is.numeric(theta))) stop("Parameters N, k1 and theta must be numeric")
  if (!is.function(k2form)) stop("Parameter k2form must be a function taking two arguments: n and theta")
  if (length(as.list(args(k2form)))!=3) stop("Parameter k2 must be a function taking two arguments: n and theta")
  if (!is.null(grad_mincost)) {
    if (!is.function(grad_mincost)) stop("Parameter grad_mincost must be a function taking three arguments: N, k1, and theta")
    if (length(as.list(args(grad_mincost)))!=4) stop("Parameter grad_mincost must be a function taking three arguments: N, k1, and theta")
  }
  if (!(mode %in% c("asymptotic","empirical"))) stop("Parameter mode must be either 'asymptotic' or 'empirical'")
  if (mode=="empirical") if (n_boot*alpha < 10) stop("Parameter nboot too low for this level of alpha")
  if (!is.numeric(alpha)) stop("Parameter alpha must be numeric")
  if (alpha <= 0 | alpha>=1) stop("Parameter alpha must be >0 and <1")

  if (mode=="asymptotic") {
    if (is.null(sigma)) {
      par_mat=cbind(N,k1,theta)
      mu=colMeans(par_mat)
      sigma_hat=var(par_mat)
    } else {
      mu=c(N,k1,theta)
      sigma_hat=sigma
    }

    ## Determine function for grad_mincost if not specified
    if (is.null(grad_mincost)) {
      grad_mincost=function(N,k1,theta) {
        dx=1e-5 # estimate gradient by using this difference
        par2=outer(rep(1,2+length(theta)),c(N,k1,theta))
        par2=par2 + dx*diag(dim(par2)[1]) # parameters, shifted by dx one-at-a-time
        cost2=optimal_holdout_size(N=par2[,1],k1=par2[,2],theta=par2[,3:dim(par2)[2]],k2form=k2form)$cost
        cost1=optimal_holdout_size(N,k1,theta,k2form=k2form)$cost
        return(t((cost2-cost1)/dx))
      }
    }

    # Parameters for asymptotic confidence interval
    xcost=optimal_holdout_size(mu[1],mu[2],mu[3:length(mu)],k2form=k2form)$cost
    gamma_est=grad_mincost(mu[1],mu[2],mu[3:length(mu)])
    z_a=-qnorm(alpha/2)
    m=length(N)
    # Variation from estimated value in asymptotic confidence interval
    di=z_a*sqrt(((gamma_est) %*% sigma_hat %*% t(gamma_est))/m)

    cx=xcost+c(-di,di)
    names(cx)=c("lower","upper")
    return(cx)
  }

  if (mode=="empirical") {

    if (is.null(sigma)) {
      par_mat=cbind(N,k1,theta)
      mu=colMeans(par_mat)
      sigma_hat=var(par_mat)
      m=dim(par_mat)[1]
    } else {
      mu=c(N,k1,theta)
      sigma_hat=sigma
      m=n_boot

      # Sample, allowing that some parameters may have zero variance
      w=which(sigma_hat[cbind(1:length(mu),1:length(mu))]<1e-20)
      wc=setdiff(1:length(mu),w)
      if (length(w)>0) {
        psub=outer(rep(1,n_boot),mu[w])
        psubc=rmnorm(n_boot,mean=mu[wc],varcov=sigma_hat[wc,wc])
        par_mat=outer(rep(1,n_boot),rep(0,length(mu)))
        par_mat[,w]=psub; par_mat[,wc]=psubc
      } else par_mat=rmnorm(n_boot,mean=mu,varcov=sigma_hat)
    }

    # Random seed
    if (!is.null(seed)) set.seed(seed)

    # Populate with mean parameters from bootstrap resamples
    ci_mat=matrix(0,n_boot,dim(par_mat)[2])
    for (i in 1:n_boot) {
      sboot=sample(m,m,replace=TRUE)
      ci_mat[i,]=colMeans(par_mat[sboot,])
    }

    # Compute min costs for bootstrap resamples
    cost_boot=optimal_holdout_size(ci_mat[,1],ci_mat[,2],ci_mat[,3:dim(ci_mat)[2]],k2form=k2form)

    # Estimate confidence interval as quantile
    cx=quantile(cost_boot$cost,c(alpha/2,1-(alpha/2)))
    names(cx)=c("lower","upper")
    return(cx)
  }
}










##' Power law function
##'
##' @export
##' @name powerlaw
##' @description Power law function for modelling learning curve (taken to mean change in expected loss per sample with training set size)
##'
##' Recommended in [review of learning curve forms](https://arxiv.org/abs/2103.10948)
##'
##' If `theta=c(a,b,c)` then models as `a n^(-b) + c`. Note `b` is negated.
##'
##' Note that `powerlaw(n,c(a,b,c))` has limit `c` as `n` tends to infinity, if `a,b > 0`
##'
##' @param n Set of training set sizes to evaluate
##' @param theta Parameter of values
##' @return Vector of values of same length as `n`
##' @examples
##'
##' ncheck=seq(1000,10000)
##' plot(ncheck, powerlaw(ncheck, c(5e3,1.2,0.3)),type="l",xlab="n",ylab="powerlaw(n)")
##'
powerlaw=function(n,theta)  (theta[1] *n^(-theta[2]) + theta[3])





##' Gradient of optimal holdout size (power law)
##'
##'
##' @export
##' @name grad_nstar_powerlaw
##' @description Compute gradient of optimal holdout size assuming a power-law form of k2
##'
##' Assumes cost function is ``l(n;k1,N,theta) = k1 n + k2(n;theta) (N-n)`` with ``k2(n;theta)=k2(n;a,b,c)= a n^(-b) + c``
##' @keywords estimation
##' @param N Total number of samples on which the predictive score will be used/fitted. Can be a vector.
##' @param k1 Cost value in the absence of a predictive score. Can be a vector.
##' @param theta Parameters for function k2(n) governing expected cost to an individual sample given a predictive score fitted to n samples. Can be a matrix of dimension n x n_par, where n_par is the number of parameters of k2.
##' @return List/data frame of dimension (number of evaluations) x 5 containing partial derivatives of nstar (optimal holdout size) with respect to N, k1, a, b, c respectively.
##' @examples
##'
##' # Evaluate optimal holdout set size for a range of values of k1, and compute
##' #  derivative
##' N=10000;
##' k1=seq(0.1,0.5,length=20)
##' A=3; B=1.5; C=0.15; theta=c(A,B,C)
##'
##' nstar=optimal_holdout_size(N,k1,theta)
##' grad_nstar=grad_nstar_powerlaw(N,k1,theta)
##'
##' plot(0,type="n",ylim=c(-2000,500),xlim=range(k1),xlab=expression("k"[1]),
##'   ylab="Optimal holdout set size")
##' lines(nstar$k1,nstar$size,col="black")
##' lines(nstar$k1,grad_nstar[,2],col="red")
##' legend("bottomright",c(expression("n"["*"]),
##'     expression(paste(partialdiff[k1],"n"["*"]))),
##'     col=c("black","red"),lty=1)
##'
grad_nstar_powerlaw=function(
  N,
  k1,
  theta
) {

  ## Error handlers
  if (!is.numeric(c(N,k1))) stop("Parameters N and k1 must be numeric")
  if (!(is.numeric(theta)|is.matrix(theta))) stop("Parameter theta must be a vector or matrix")
  if ((length(N)>1 | length(k1)>1)|!is.null(dim(theta))) {
    n=max(length(N),length(k1),dim(theta)[1])
    if (is.null(dim(theta))) theta=t(matrix(theta,length(theta),n))
    if (!all(c(length(k1),length(N),dim(theta)[1]) %in% c(1,n))) stop("If vectors, N and k1 must have the same length, or have length one. Parameter theta must either be a vector or a matrix of dimension length(n) x n_par, where n_par is the number of parameters of k2.")
    if (length(N)==1) N=rep(N,n)
    if (length(k1)==1) k1=rep(k1,n)
  } else theta=matrix(theta,1,length(theta))
  if (dim(theta)[2]!=3) stop("Parameter theta must have length 3 or dimension n x 3")


  ns=optimal_holdout_size(N,k1,theta)$size

  A=theta[,1]; B=theta[,2]; C=theta[,3]


  # partial derivatives of ns*
  dnstarda = (B*N*ns - (B-1)*(ns^2))/(B*(B+1)*A*N - B*(B-1)*A*ns)
  dnstardb = -(ns*N*( -1 + B*log(ns))-(ns^2)*(- 1 + (B-1)*log(ns)))/(B*(B+1)*N - B*(B-1)*ns)
  dnstardc = 1/(B*(B+1)*A*N*ns^(-B-2) - B*(B-1)*A*ns^(-B-1))
  dnstardk1 = -1/(B*(B+1)*B*N*ns^(-B-2) - B*(B-1)*B*ns^(-B-1))
  dnstardN = B*ns/(B*(B+1)*N - B*(B-1)*ns)

  return(cbind(
    dnstardN,
    dnstardk1,
    dnstarda,
    dnstardb,
    dnstardc
  ))
}






##' Gradient of minimum cost (power law)
##'
##'
##' @export
##' @name grad_mincost_powerlaw
##' @description Compute gradient of minimum cost assuming a power-law form of k2
##'
##' Assumes cost function is ``l(n;k1,N,theta) = k1 n + k2(n;theta) (N-n)`` with ``k2(n;theta)=k2(n;a,b,c)= a n^(-b) + c``
##' @keywords estimation
##' @param N Total number of samples on which the predictive score will be used/fitted. Can be a vector.
##' @param k1 Cost value in the absence of a predictive score. Can be a vector.
##' @param theta Parameters for function k2(n) governing expected cost to an individual sample given a predictive score fitted to n samples. Can be a matrix of dimension n x n_par, where n_par is the number of parameters of k2.
##' @return List/data frame of dimension (number of evaluations) x 5 containing partial derivatives of nstar (optimal holdout size) with respect to N, k1, a, b, c respectively.
##' @examples
##' # Evaluate minimum for a range of values of k1, and compute derivative
##' N=10000;
##' k1=seq(0.1,0.5,length=20)
##' A=3; B=1.5; C=0.15; theta=c(A,B,C)
##'
##' mincost=optimal_holdout_size(N,k1,theta)
##' grad_mincost=grad_mincost_powerlaw(N,k1,theta)
##'
##' plot(0,type="n",ylim=c(0,1560),xlim=range(k1),xlab=expression("k"[1]),
##'   ylab="Optimal holdout set size")
##' lines(mincost$k1,mincost$cost,col="black")
##' lines(mincost$k1,grad_mincost[,2],col="red")
##' legend(0.2,800,c(expression(paste("l(n"["*"],")")),
##'                        expression(paste(partialdiff[k1],"l(n"["*"],")"))),
##'     col=c("black","red"),lty=1,bty="n")
##'
grad_mincost_powerlaw=function(
  N,
  k1,
  theta
) {

  ## Error handlers
  if (!is.numeric(c(N,k1))) stop("Parameters N and k1 must be numeric")
  if (!(is.numeric(theta)|is.matrix(theta))) stop("Parameter theta must be a vector or matrix")
  if ((length(N)>1 | length(k1)>1)|!is.null(dim(theta))) {
    n=max(length(N),length(k1),dim(theta)[1])
    if (is.null(dim(theta))) theta=t(matrix(theta,length(theta),n))
    if (!all(c(length(k1),length(N),dim(theta)[1]) %in% c(1,n))) stop("If vectors, N and k1 must have the same length, or have length one. Parameter theta must either be a vector or a matrix of dimension length(n) x n_par, where n_par is the number of parameters of k2.")
    if (length(N)==1) N=rep(N,n)
    if (length(k1)==1) k1=rep(k1,n)
  } else theta=matrix(theta,1,length(theta))
  if (dim(theta)[2]!=3) stop("Parameter theta must have length 3 or dimension n x 3")


  ns=optimal_holdout_size(N,k1,theta)$size

  A=theta[,1]; B=theta[,2]; C=theta[,3]


  # partial derivatives of ns*
  dmincostda = (N-ns)*(ns^(-B))
  dmincostdb = -log(ns)*(N-ns)*A*(ns^(-B))
  dmincostdc = N-ns
  dmincostdk1 = ns
  dmincostdN = A*(ns^(-B)) + C

  return(cbind(
    dmincostdN,
    dmincostdk1,
    dmincostda,
    dmincostdb,
    dmincostdc
  ))
}






##' Fit power law curve
##'
##' @export
##' @name powersolve
##' @description  Find least-squares solution: MLE of (a,b,c) under model
##'  ``y_i = a x_i^-b + c + e_i``;
##'  ``e_i~N(0,y_var_i)``
##'
##' @keywords estimation aspre
##' @param x X values
##' @param y Y values
##' @param init Initial values of (a,b,c) to start. Default c(20000,2,0.1)
##' @param y_var Optional parameter giving sampling variance of each y value. Defaults to 1.
##' @param estimate_s Parameter specifying whether to also estimate s (as above). Defaults to FALSE (no).
##' @param ... further parameters passed to optim. We suggest specifying lower and upper bounds for (a,b,c); e.g. lower=c(1,0,0),upper=c(10000,3,1)
##' @return List (output from `optim`) containing MLE values of (a,b,c)
##' @examples
##'
##' # Retrieval of original values
##' A_true=2000; B_true=1.5; C_true=0.3; sigma=0.002
##'
##' X=1000*abs(rnorm(10000,mean=4))
##' Y=A_true*(X^(-B_true)) + C_true + rnorm(length(X),sd=sigma)
##'
##' c(A_true,B_true,C_true)
##' powersolve(X[1:10],Y[1:10])$par
##' powersolve(X[1:100],Y[1:100])$par
##' powersolve(X[1:1000],Y[1:1000])$par
##' powersolve(X[1:10000],Y[1:10000])$par
powersolve=function(x,y,init=c(20000,2,0.1),y_var=rep(1,length(y)),estimate_s=FALSE,...) {
  if (!estimate_s) {
    fabc=function(abc) {
      out=sum( ((y- (abc[1]*(x^(-abc[2])) + abc[3]))^2)/(y_var))
      if (is.finite(out)) return(out) else return(1e10)
    }
    out=suppressWarnings(optim(par=init,fn=fabc,control=list(parscale=init),...))
  } else {
    fabcs=function(abcs) {
      out=-(sum( -((y- (abcs[1]*(x^(-abcs[2])) + abcs[3]))^2 / (2*y_var*(abcs[4]^2))) - log(sqrt(2*3.1415*y_var)*abcs[4])))
      if (is.finite(out)) return(out) else return(1e10)
    }
    out=suppressWarnings(optim(par=c(init,0.05),fn=fabcs,control=list(parscale=c(init,0.5)),...))
  }
  return(out)
}





##' General solver for power law curve
##'
##' @export
##' @name powersolve_general
##' @description  Find least-squares solution: MLE of (a,b,c) under model
##'  ``y_i = a x_i^-b + c + e_i``;
##'  ``e_i~N(0,y_var_i)``
##'
##' Try a range of starting values and refine estimate.
##'
##' Slower than a single call to ``powersolve()``
##'
##' @keywords estimation aspre
##' @param x X values
##' @param y Y values
##' @param y_var Optional parameter giving sampling variance of each y value. Defaults to 1.
##' @param ... further parameters passed to optim. We suggest specifying lower and upper bounds for (a,b,c); e.g. lower=c(1,0,0),upper=c(10000,3,1)
##' @return List (output from `optim`) containing MLE values of (a,b,c)
##' @examples
##'
##' # Retrieval of original values
##' A_true=2000; B_true=1.5; C_true=0.3; sigma=0.002
##'
##' X=1000*abs(rnorm(10000,mean=4))
##' Y=A_true*(X^(-B_true)) + C_true + rnorm(length(X),sd=sigma)
##'
##' c(A_true,B_true,C_true)
##' powersolve_general(X[1:10],Y[1:10])$par
##' powersolve_general(X[1:100],Y[1:100])$par
##' powersolve_general(X[1:1000],Y[1:1000])$par
powersolve_general=function(x,y,y_var=rep(1,length(x)),...) {
  # Values of a,b,c to trial
  atry=c(0.5,1,5,10,100,1000); btry=c(0.1,0.5,1,1.5,2); ctry=c(0,0.01,0.5)

  # Current estimate
  cur=c(1,1,1); Z=powersolve(x,y,y_var,init=cur,...); best=Z$par; top=Z$value
  for (aa in atry) for (bb in btry) for (cc in ctry) {
    Z=powersolve(x,y,y_var,init=cur,...);
    if (Z$value<top) {
      top=Z$value
      best=Z$par
    }
  }

  for (i in 1:10) best=powersolve(x,y,y_var,init=best,...)$par
  return(powersolve(x,y,y_var,init=best,...))
}





##' Standard error matrix for learning curve parameters (power law)
##'
##'
##' @export
##' @name powersolve_se
##' @description Find approximate standard error matrix for ``(a,b,c)`` under power law model for learning curve.
##'
##' Assumes that
##'
##'   ``y_i= a x_i^-b + c + e, e~N(0,s^2 y_var_i^2)``
##'
##' Standard error can be computed either asymptotically using Fisher information (`method='fisher'`) or boostrapped (`method='bootstrap'`)
##'
##' These estimate different quantities: the asymptotic method estimates
##'
##' ``Var[MLE(a,b,c)|X,y_var]``
##'
##' and the boostrap method estimates
##'
##' ``Var[MLE(a,b,c)]``.
##'
##' @keywords estimation aspre
##' @param x X values (typically training set sizes)
##' @param y Y values (typically observed cost per individual/sample)
##' @param method One of 'fisher' (for asymptotic variance via Fisher Information) or 'bootstrap' (for Bootstrap)
##' @param init Initial values of (a,b,c) to start when computing MLE. Default c(20000,2,0.1)
##' @param y_var Optional parameter giving sampling variance of each y value. Defaults to 1.
##' @param n_boot Number of bootstrap resamples. Only used if method='bootstrap'. Defaults to 1000
##' @param seed Random seed for bootstrap resamples. Defaults to NULL.
##' @param ... further parameters passed to optim. We suggest specifying lower and upper bounds; since optim is called on (a*1000^-b,b,c), bounds should be relative to this; for instance, lower=c(0,0,0),upper=c(100,3,1)
##' @return Standard error matrix; approximate covariance matrix of MLE(a,b,c)
##' @examples
##'
##' A_true=10; B_true=1.5; C_true=0.3; sigma=0.1
##'
##' set.seed(31525)
##'
##' X=1+3*rchisq(10000,df=5)
##' Y=A_true*(X^(-B_true)) + C_true + rnorm(length(X),sd=sigma)
##'
##' # 'Observations' - 100 samples
##' obs=sample(length(X),100,rep=FALSE)
##' Xobs=X[obs]; Yobs=Y[obs]
##'
##' # True covariance matrix of MLE of a,b,c on these x values
##' ntest=100
##' abc_mat_xfix=matrix(0,ntest,3)
##' abc_mat_xvar=matrix(0,ntest,3)
##' E1=A_true*(Xobs^(-B_true)) + C_true
##' for (i in 1:ntest) {
##'   Y1=E1 + rnorm(length(Xobs),sd=sigma)
##'   abc_mat_xfix[i,]=powersolve(Xobs,Y1)$par # Estimate (a,b,c) with same X
##'
##'   X2=1+3*rchisq(length(Xobs),df=5)
##'   Y2=A_true*(X2^(-B_true)) + C_true + rnorm(length(Xobs),sd=sigma)
##'   abc_mat_xvar[i,]=powersolve(X2,Y2)$par # Estimate (a,b,c) with variable X
##' }
##'
##' Ve1=var(abc_mat_xfix) # empirical variance of MLE(a,b,c)|X
##' Vf=powersolve_se(Xobs,Yobs,method='fisher') # estimated SE matrix, asymptotic
##'
##' Ve2=var(abc_mat_xvar) # empirical variance of MLE(a,b,c)
##' Vb=powersolve_se(Xobs,Yobs,method='bootstrap',n_boot=200) # estimated SE matrix, bootstrap
##'
##' cat("Empirical variance of MLE(a,b,c)|X\n")
##' print(Ve1)
##' cat("\n")
##' cat("Asymptotic variance of MLE(a,b,c)|X\n")
##' print(Vf)
##' cat("\n\n")
##' cat("Empirical variance of MLE(a,b,c)\n")
##' print(Ve2)
##' cat("\n")
##' cat("Bootstrap-estimated variance of MLE(a,b,c)\n")
##' print(Vb)
##' cat("\n\n")
##'
powersolve_se=function(x,y,method='fisher',init=c(20000,2,0.1),y_var=rep(1,length(y)),n_boot=1000,seed=NULL,...) {
  if (method=="fisher") {
    ## Fisher Information Matrix - straightforward to compute analytically
    FI_mat=function(x,v,a,b,c,s)
      -(1/(s^2 * v))*cbind(
        c(-x^(-2*b),a*(x^(-2*b))*log(x),-x^(-b),0),
        c(a*(x^(-2*b))*log(x),-(a^2)*(x^(-2*b))*(log(x)^2),a*(x^(-b))*log(x),0),
        c(-x^(-b),a*(x^(-b))*log(x),-1,0),
        c(0,0,0,-2*v))

    ## Need to estimate s as well; call powersolve
    abcs=powersolve(x,y,y_var=y_var,init=init,estimate_s=TRUE,...)$par

    fmat=matrix(0,4,4); for (i in 1:length(x)) fmat=fmat + FI_mat(x[i],y_var[i],abcs[1],abcs[2],abcs[3],abcs[4])

    # Check if fmat is invertible
    if ("matrix" %in% class(try(solve(fmat),silent=T))) {
      vmat=solve(fmat)[1:3,1:3]
      return(vmat)
    } else {
      return(matrix(NA,3,3))
    }
  }
  if (method=="bootstrap") {

    # set seed
    if (!is.null(seed)) set.seed(seed)

    # Compute bootstrap resamples
    xboot=matrix(0,n_boot,3)
    for (i in 1:n_boot) {
      sub=sample(length(x),length(x),replace=TRUE)
      x1=init
      for (j in 1:5) x1=powersolve(x[sub],y[sub],y_var=y_var[sub],init=x1,estimate_s=FALSE,...)$par
      xboot[i,]=x1
    }
    vmat=var(xboot)
    return(vmat)
  }
}






##' Finds best value of n to sample next
##'
##' @export
##' @name next_n
##' @description Recommends a value of `n` at which to next evaluate individual cost in order to most accurately estimate optimal holdout size. Currently only for use with a power-law parametrisation of k2.
##'
##' Approximately finds a set of n points which, given estimates of cost, minimise width of 95% confidence interval around OHS. Uses a greedy algorithm, so various parameters can be learned along the way.
##'
##' Given existing training set size/k2 estimates `nset` and `k2`, with `var_k2[i]=variance(k2[i])`, finds, for each candidate point `n[i]`, the median width of the 90% confidence interval for OHS if
##'
##' ``nset <- c(nset,n[i])``
##' ``var_k2 <- c(var_k2,mean(var_k2))``
##' ``k2 <- c(k2,rnorm(powerlaw(n[i],theta),variance=mean(var_k2)))``
##'
##' @param n Set of training set sizes to evaluate
##' @param nset Training set sizes for which a loss has been evaluated
##' @param k2 Estimated k2() at training set sizes `nset`
##' @param var_k2 Variance of error in k2() estimate at each training set size.
##' @param mode Mode for calculating OHS CI (passed to `ci_ohs`): 'asymptotic' or 'empirical'
##' @param N Total number of samples on which the model will be fitted/used
##' @param k1 Mean loss per sample with no predictive score in place
##' @param nmed number of times to re-evaluate d and confidence interval width.
##' @param ... Passed to `powersolve` and `powersolve_se`
##' @return Vector `out` of same length as `n`, where `out[i]` is the expected width of the 95% confidence interval for OHS should `n` be added to `nset`.
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
##' nstart=10
##' vwmin=0.001; vwmax=0.005
##' nset0=round(runif(nstart,1000,N/2))
##' var_k2_0=runif(nstart,vwmin,vwmax)
##' k2_0=rnorm(nstart,mean=powerlaw(nset0,theta_true),sd=sqrt(var_k2_0))
##'
##' # We estimate theta from these three points
##' theta0=powersolve(nset0,k2_0,y_var=var_k2_0,lower=theta_lower,upper=theta_upper,init=theta_true)$par
##'
##' # We will estimate the posterior at these values of n
##' n=seq(1000,N,length=1000)
##'
##' # Mean and variance
##' p_mu=mu_fn(n,nset=nset0,k2=k2_0,var_k2 = var_k2_0, N=N,k1=k1,theta=theta0,k_width=kw0,var_u=vu0)
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
##' ## Add line corresponding to recommended new point. This is slow.
##' nn=seq(1000,N,length=20)
##' exp_imp <- next_n(nn,nset=nset0,k2=k2_0,var_k2 = var_k2_0, N=N,k1=k1,nmed=10,
##'                      lower=theta_lower,upper=theta_upper)
##' abline(v=nn[which.min(exp_imp)])
##'
##'
next_n=function(n,nset,k2,N,k1,nmed=100,var_k2=rep(1,length(nset)),mode="asymptotic",...) {
  out=rep(0,length(n))
  theta=powersolve_general(nset,k2,y_var=var_k2)$par
  for (i in 1:length(n)) {
    out_i=rep(0,nmed)
    for (j in 1:nmed) {
      # Data appended with new candidate point
      nsetx=c(nset,n[i]);
      var_k2x=c(var_k2,mean(var_k2));
      k2x=c(k2,rnorm(1,mean=powerlaw(n[i],theta),sd=sqrt(mean(var_k2))))

      # Parameter estimates with new candidate point
      thetax=powersolve(nsetx,k2x,y_var=var_k2x,init=theta,...)$par
      covx=powersolve_se(nsetx,k2x,y_var=var_k2x,method="fisher",init=theta,...)
      cov_all=matrix(0,5,5); cov_all[3:5,3:5]=covx

      if (!is.na(cov_all[3,3])) {
        cx=ci_ohs(N,k1,thetax,sigma=cov_all,mode=mode,grad_nstar=grad_nstar_powerlaw)
        out_i[j]=max(cx)-min(cx)
      } else out[j]=Inf
    }
    out[i]=median(out_i,na.rm=TRUE)
  }
  return(out)
}

