% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OptHoldoutSize_estimation.R
\name{optimal_holdout_size_emulation}
\alias{optimal_holdout_size_emulation}
\title{Estimate optimal holdout size under semi-parametric assumptions}
\usage{
optimal_holdout_size_emulation(
  nset,
  k2,
  var_k2,
  N,
  k1,
  var_u = 1e+07,
  k_width = 5000,
  k2form = powerlaw,
  theta = powersolve_general(nset, k2, y_var = var_k2)$par,
  npoll = 1000,
  ...
)
}
\arguments{
\item{nset}{Training set sizes for which a cost has been evaluated}

\item{k2}{Estimated values of k2() at training set sizes \code{nset}}

\item{var_k2}{Variance of error in k2 estimate at each training set size.}

\item{N}{Total number of samples on which the model will be fitted/used}

\item{k1}{Mean cost per sample with no predictive score in place}

\item{var_u}{Marginal variance for Gaussian process kernel. Defaults to 1e7}

\item{k_width}{Kernel width for Gaussian process kernel. Defaults to 5000}

\item{k2form}{Functional form governing expected cost per sample given sample size. Should take two parameters: n (sample size) and theta (parameters). Defaults to function \code{powerlaw}.}

\item{theta}{Current estimates of parameter values for k2form. Defaults to the MLE power-law solution corresponding to n,k2, and var_k2.}

\item{npoll}{Check npoll equally spaced values between 1 and N for minimum. If NULL, check all values (this can be slow). Defaults to 1000}

\item{...}{Passed to function \code{optimise()}}
}
\value{
Object of class 'optholdoutsize_emul' with elements "cost" (minimum cost),"size" (OHS),"nset","k2","var_k2","N","k1","var_u","k_width","theta" (parameters)
}
\description{
Compute optimal holdout size for updating a predictive score given a set of training set sizes and estimates of mean cost per sample at those training set sizes.

This is essentially a wrapper for function \code{mu_fn()}.
}
\examples{

# See examples for mu_fn()
}
\keyword{emulation}
\keyword{estimation}
