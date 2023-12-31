% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitSPT.R
\name{fitSPT}
\alias{fitSPT}
\title{Nonparametric MLE for semi-parametric transformation regression model (SPT)}
\usage{
fitSPT(data, time = "time", status = "status", formula = ~1, Gfun = "PH")
}
\arguments{
\item{data}{a data frame containing variables names in the \code{formula}.}

\item{time}{a character string specifying the variable name in the \code{data} for the event time.}

\item{status}{a character string specifying the variable name in the \code{data} for the censoring indicator: \code{1} indicates the time is observed, and \code{0} indicates censored.}

\item{formula}{an object of class \code{\link[stats]{formula}}: a symbolic description of the model to be fitted.}

\item{Gfun}{a character string specifying the link function, \code{Gfun = "PH"} (default) or \code{Gfun = "PO"} in the SPT model.}
}
\value{
a list of the following components:
\describe{
 \item{\code{beta}}{a data frame containing the estimate, model-based standard eror (SE), and robust SE for regresson coefficients.}
 \item{\code{dLambda}}{a data frame containing the estimate, model-based standard eror (SE), and robust SE for the jump size of the baseline function at observed time points.}
 \item{\code{Lambda}}{a data frame containing the estimate, model-based standard eror (SE), and robust SE for the baseline function at observed time points.}
 \item{\code{varcov}}{a list containing \code{model}, the model-based variance-covariance matrix, and \code{robust}, the robust variance-covariance matrix for all the parameters.}
 \item{\code{Psi.theta}}{a data frame containing the estimated asymptotic expansion of the NPMLE estimator.}
 \item{\code{call}}{a list containing the specified values of the arguments.}
 \item{\code{convergence}}{a logical value indicating whether the maximization of the log-likelihood converges or not.}
 \item{\code{niter}}{an interger which is the number of iterations for the maximizaiton of the log-likelihood.}
}
}
\description{
Nonparametric MLE for semi-parametric transformation regression model (SPT)
}
\examples{
\code{data(BMT, package = "SemiCompRisks")}
\code{data = BMT \%>\%
  mutate(g = factor(g, levels = c(2, 3, 1),
                    labels = c("AML-low", "AML-high", "ALL")))}
\code{fitSPT(data, time = "T1", status = "delta1",
              formula = ~ g, Gfun = "PH")$beta}
}
