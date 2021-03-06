% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_fitting_functions.R
\name{regress.back.fitness.vs.effect}
\alias{regress.back.fitness.vs.effect}
\title{Linear regression of background fitness against effects}
\usage{
regress.back.fitness.vs.effect(fitness.of.backs, effects.matrix, n.muts)
}
\arguments{
\item{fitness.of.backs}{Matrix with background fitness for \code{effects.matrix}}

\item{effects.matrix}{Effect of mutation in that genotype}

\item{n.muts}{Number of mutations}
}
\value{
List: \cr
[[1]] \code{p.vals} P-value for linear regression for each mutation \cr
[[2]] \code{lm.intercepts} Intercepts for each mutation
[[3]] \code{lm.slopes} Slopes for each mutation
[[4]] \code{P} Sum of the logs of the p-values.
[[5]] \code{fitness.of.backs} Matrix with fitness of backgrounds when each mutation (columns) is added to each genotype (rows).
[[6]] \code{effects.matrix} Matrix with fitness effect when given mutation (column) is added to given create genotype (row).
If there is insufficient data to do regression, a warning is returned.
}
\description{
Linear regression of background fitness against effects
}
\details{
For each mutation, function does simple linear regression using \code{lm()}.
The product of logged p-values (\code{P}) is the summary statistic used in model selection.
When the model that generated the data and the model analyzing the data match, the expected slope
of the regression line is zero.
}

