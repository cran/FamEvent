\name{summary.penmodel}
\alias{summary.penmodel}

\title{
Summary of a penetrance model fit
}
\description{
Summary method for \code{penmodel}. Provides summary of a penetrance model fit including model parameter and penetrance estimates.
}
\usage{
\method{summary}{penmodel}(object, \ldots)
}
\arguments{
  \item{object}{
An object class of \code{penmodel}, a model fit by the functions \code{penmodel} or \code{penmodelEM}.
}
 \item{...}{not used }

}
\value{
It displays the following summary values:
\item{parameters}{Parameter estimates of baseline parameters (\eqn{\lambda, \rho}) and regression coefficients for gender and mutation status (\eqn{\beta_s, \beta_g}) including their standard errors and also robust standard errors.}
\item{pen70}{ Penetrance estimates by age 70 specific to gender and mutation-status subgroups including their standard errors and 95\% confidence intervals.}
}

\author{
Yun-Hee Choi	
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{penmodel}, \link{plot.penmodel}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ summary }
