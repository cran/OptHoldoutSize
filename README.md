## OptHoldoutSize: an R package for estimating the optimal holdout set size for a predictive risk score to be deployed in a population.

This R package implements procedures for estimating an 'optimal holdout size' for a predictive score in order for it to be safely updated. Procedures are detailed in the manuscript 'Optimal sizing of a holdout set for safe predictive model updating' by Sami Haidar-Wehbe, Samuel R. Emerson, Louis J. M. Aslett, and James Liley.

When a predictive risk score for binary outcome $Y$ given covariates $X$ is deployed in a population, it may be used to guide interventions so as to avoid $Y$. This makes it difficult to update the predictive score safely, since $X$ can influence incidence of $Y$ in two ways: through the system being modelled, or through the predictive score itself. 

A simple way to safely update a predictive is to with-hold calculation of the risk score for a proportion of the population maintained as a 'holdout' set. The predictive score can then be updated using data $X$, $Y$ from this holdout set. A question naturally arises over how large this hold-out set should be: too small, and a new predictive score cannot be trained sufficiently accurately; too large, and too many members of the population miss out on potential benefits of the risk score.

To download and install this package, use

```
install.packages("OptHoldoutSize")
library(OptHoldoutSize)
```

For examples demonstrating use of this package, see vignettes ```simulated_example``` and ```ASPRE_example```. For a comparison of the two major algorithms implemented in this package, see vignette ```comparison_of_algorithms```.
