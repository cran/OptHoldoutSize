% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OptHoldoutSize_aspre.R
\name{aspre}
\alias{aspre}
\title{Computes ASPRE score}
\usage{
aspre(X)
}
\arguments{
\item{X}{matrix, assumed to be output of sim_random_aspre with parameter params=params_aspre and transformed using add_aspre_interactions}
}
\value{
vector of scores.
}
\description{
Computes ASPRE model prediction on a matrix \code{X} of covariates

Full ASPRE model from https://www.nejm.org/doi/suppl/10.1056/NEJMoa1704559/suppl_file/nejmoa1704559_appendix.pdf

Model is to predict gestational age at PE; that is, a higher score indicates a lower PE risk, so coefficients are negated for model to predict PE risk.
}
\examples{

# Load ASPRE related data
data(params_aspre)

X=sim_random_aspre(1000,params_aspre)
Xnew=add_aspre_interactions(X)

aspre_score=aspre(Xnew)

plot(density(aspre_score))

}
\keyword{aspre}
