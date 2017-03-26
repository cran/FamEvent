\name{plot.penmodel}
\alias{plot.penmodel}
\title{
Plot method for \code{penmodel}
}
\description{
Plots penetrance curves estimated from the fitted penetrance model and overlays non-parametric penetrance curves estimated from the data without proabands.
}
\usage{
\method{plot}{penmodel}(x, print=TRUE, mark.time=FALSE, conf.int=FALSE, \ldots)
}
\arguments{
  \item{x}{
An object class of \code{penmodel}, a fitted model by \code{\link{penmodel}} or \code{\link{penmodelEM}} functions.
}
 \item{print}{Logical; if \code{TRUE}, displays parameter estimates and penetrance estimates by age 70.}
 \item{mark.time}{Logical; if \code{TRUE}, curves are marked at each censoring time, otherwise, no labeling is done.}
 \item{conf.int}{Logical; if \code{TRUE}, displays 95\% confidence intervals for both parametric and non-parametric penetrance estimates for each subgroup and returns their \code{lower} and \code{upper} limits.}
 \item{...}{Other parameters to be passed through to plotting functions.}

}
\details{
The 95\% confidence intervals for the parametric penetrance curves were
obtained based on 1000 simulations of the parameters,
assuming a multivariate normal distribution for the estimated
parameters with their variance covariance matrix. See \code{\link{penci}} for more details.
}
\value{
It displays the following summary values:
\item{coefficients}{Parameter estimates of transformed baseline parameters (\eqn{\lambda, \rho}) and regression coefficients for gender and mutation status (\eqn{\beta_s, \beta_g}).}
\item{pen70}{ Penetrance estimates by age 70 specific to gender and mutation-status subgroups.}
\item{x.age}{Vetor of ages of onsest ranging from \code{agemin} to 80 years}
\item{pen}{Penetrance estimates at each age in \code{x.age}, specific to gender and mutation-status subgroups.}
\item{lower}{Lower limits of 95\% confidence interval estimates for penetrance at each age in \code{x.age}, specific to gender and mutation status subgroups.}
\item{upper}{Upper limits of 95\% confidence interval estimates for penetrance at each age in \code{x.age}, specific to gender and mutation status subgroups.}
}

\author{
Yun-Hee Choi	
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{penmodel}, \link{print.penmodel}, \link{penmodelEM}, \link{summary.penmodel}, \link{print.summary.penmodel}}, 

\code{\link{simfam}, \link{penci}, \link{penf}}	
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
\keyword{Penetrance plot}
