\name{penci}
\alias{penci}
\title{Estimate the confidence intervals for the penetrances}
\description{
  Estimates the simulation-based confidence intervals for the penetrances and the standard errors of the penetrance estimates at a given age, specific to gender and mutation status subgroups.
}  
\usage{
  penci(est, cov, age=70, base.dist="Weibull", frailty.dist=NULL, agemin, n=1000)
}
\arguments{
  \item{est}{
    Parameter estimates of transformed baseline parameters (\eqn{\lambda, \rho}) and 
    regression coefficients for gender and mutation status (\eqn{\beta_s, \beta_g})
    for the assumed penetrance model.
  }
  \item{cov}{
    Variance-covariance matrix of the parameter estimates.
  }
  \item{age}{
    Specified age at which the penetrance is computed. 
  }
  \item{base.dist}{
    Choice of baseline hazard distribution for the penetrance function. Possible choices are \code{"Weibull"}, \code{"loglogistic"}, \code{"Gompertz"}, \code{"lognormal",} or \code{"gamma"}. Default is \code{"Weibull"}.
  }
  \item{frailty.dist}{
    Choice of frailty distribution, either \code{"gamma"} or \code{"lognormal"}, if the penetrance function is based on a shared frailty model. Otherwise, \code{frailty.dist = NULL} is set as default.
  }
  \item{agemin}{Minimum age that the pentrance function starts.}
  \item{n}{Number of Monte-Carlo simulations for calculating standard errors and 95\% confidence intervals for the penetrance estimate at a given age. Default value is \code{n = 1000}.
  }
    
}
\details{
  
  Calculations of standard errors of the penetrance estimates and 95\% confidence intervals (CIs) for the penetrance at a given age are based on Monte-Carlo simulations of the estimated penetrance model. 
  
  A multivariate normal distribution is assumed for the parameter estimates, 
  and \eqn{n} sets of the parameters are generated from the multivariate normal distribution
  with the parameter estimates and their variance and covariance matrix.
  For each simulated set, a penetrance estimate is calculated at a given age by substituting the simulated parameters into the penetrance function.
  
  The standard error of the penetrance estimate at a given age is calculated by the standard deviation of penetrance estimates obtained from \eqn{n} simulations.
  
  The 95\% CI for the penetrance at a given age is calculated using the 2.5th and 97.5th percentiles of the penetrance estimates obtained from \eqn{n} simulations.
  
}
\value{
  Returns an object including the following values:
    
  \item{Estimate}{Penetrance estimates (\%) at the specified age, specific to gender and mutation status subgroups, based on the assumed penetrance model.}    
  \item{SE}{Simulation-based standard errors of the penetrance estimates, specific to gender and mutation status subgroups.}
  \item{lower}{Simulation-based 2.5th percentile of the penetrance estimates, specific to gender and mutation status subgroups.}
  \item{upper}{Simulation-based 97.5th percentile of the penetrance estimates, specific to gender and mutation status subgroups.}
}
\author{
  Yun-Hee Choi
}
\seealso{
  \code{\link{penmodelEM}, \link{penmodel}, \link{penplot}, \link{penf}, \link{print.penmodel},\link{summary.penmodel},}
  
  \code{\link{print.summary.penmodel}, \link{plot.penmodel}} 
}
\examples{
  
  # Family data simulated from population-based design using a Weibull baseline hazard 
  
  fam <- simfam(N.fam=300, design="pop+", variation="none", base.dist="Weibull", 
                base.parms=c(0.01,3), vbeta=c(-1.13, 2.35), agemin=20, allelefreq=0.02)
  
  # Penetrance model fit for simulated family data
  
  fit <- penmodel(parms=c(0.01, 3), vbeta=c(-1.13, 2.35), data=fam, 
                  design="pop+", base.dist="Weibull")
  
  # 95\% confidence intervals for the penetrance at age 50 based on 1000 simulations
  
  penci(fit$coefficients, fit$varcov, age=50, base.dist="Weibull", agemin=20, n=1000)
  

}
\keyword{Penetrance confidence interval}
