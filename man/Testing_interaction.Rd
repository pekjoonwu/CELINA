% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CELINA.R
\name{Testing_interaction}
\alias{Testing_interaction}
\title{Test heterogeneity term using null model for one pair. Used by Testing_interaction_all}
\usage{
Testing_interaction(Y, X, covariates, kernelmat_approx_U)
}
\arguments{
\item{Y}{Gene expression vector}

\item{X}{celltype proportion for testing}

\item{covariates}{covariates for the null model}

\item{kernelmat_approx_U}{The calculated approximated kernel matrix for spatial relationship between locations.}
}
\description{
Test heterogeneity term using null model for one pair. Used by Testing_interaction_all
}
