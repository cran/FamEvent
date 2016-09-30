\name{plot.penmodel}
\alias{plot.penmodel}
\title{
Plot estimated penetrance curves from a model fit }
\description{
Plot method for \code{penmodel}. Plots penetrance curves estimated from the model fit and overlays non-parametric penetrance curves estimated from the data without proabands.
}
\usage{
\method{plot}{penmodel}(x, \ldots)
}
\arguments{
  \item{x}{
An object class of \code{penmodel}, a fitted model by \code{penmodel} function.
}
 \item{...}{additional arguments}

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
\code{\link{penmodel}, \link{summary.penmodel}, \link{simfam}}	
}

\examples{

# Simulated family data  

fam <- simfam(N.fam=300, design="pop+", base.dist="Weibull", variation="none",
       base.parms=c(0.01,3), vbeta=c(-1.13, 2.35), allelefreq=0.02, agemin=20)
 
# Fit family data

fit <- penmodel(parms=c(0.01, 3), vbeta=c(-1.13, 2.35), data=fam, 
       design="pop+", base.dist="Weibull")
       
# Plot penetrance function curves

plot(fit)
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Penetrance }
