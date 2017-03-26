\name{penf}
\alias{penf}
\title{Penetrance function}
\description{
  Calculates the penetrance for a given age, gender and mutation status
  based on the assumed penetrance model.
}  
\usage{
  penf(est, age, sex, mut, base.dist, frailty.dist=NULL, agemin)
}
\arguments{
  \item{est}{
    Parameter estimates for the assumed penetrance model including the transformed baseline parameters 
    (\eqn{\lambda, \rho}),  regression coefficients for gender and mutation status (\eqn{\beta_s, \beta_g})
    and a frailty parameter (\eqn{\kappa}) if \code{frailty.dist} is specified.
  }
  \item{age}{
    Vector of ages or a single value of age at which the penetrance function is evaluated. 
  }
  \item{sex}{1 for male, 0 for female.}
  \item{mut}{1 for mutation carrier, 0 for mutation noncarrier.}
  \item{base.dist}{
    Choice of baseline hazard distribution for the penetrance function. Possible choices are \code{"Weibull"}, \code{"loglogistic"}, \code{"Gompertz"}, \code{"lognormal",} or \code{"gamma"}. Default is \code{"Weibull"}.
  }
  \item{frailty.dist}{
    Choice of frailty distribution, either \code{"gamma"} or \code{"lognormal"}, if the penetrance function is based on a shared frailty model. Otherwise, \code{frailty.dist = NULL} is set as default.
  }
  \item{agemin}{Minimum age that the pentrance function starts.}

}
\details{
The penetrance function is defined as the probability of developing disease by time \eqn{t} given gender (\eqn{x_s}) and mutation status (\eqn{x_g}),
\deqn{ P(T < t|x_s, x_g) = 1 - S(t; x_s, x_g) , }
where \eqn{S(t; x_s, x_g)}  is the survival distribution based on a proportional hazards model with a specified baseline hazard distribution or based on a shared frailty mode with specified frailty and baseline hazard distributions.

Proportional hazards model:

\deqn{ h(t|x_s, x_g) = h_0(t) \exp(\beta_s x_s+\beta_g x_g) }
where \eqn{h_0(t)} is the baseline hazards function specified by \code{base.dist}, which depends on the shape and scale parameters, \eqn{\lambda} and \eqn{\rho}; \eqn{x_s} indicates male (1) and female (0) and \eqn{x_g} indicates carrier (1) or non-carrier (0) of a gene of interest (major gene).

Shared frailty model:

\deqn{ h(t|z, x_s, x_g) = z h_0(t) \exp(\beta_s x_s+\beta_g x_g) }
where \eqn{z} is the shared frailty whose distribution is specified by \code{frailty.dist} with associated parameter \eqn{\kappa}.
The marginal survival function is obtained by integrating the conditional survival distribution, \eqn{S(s; z, x_s, x_g) }, over frailty distribution, i.e.,
\deqn{S(t; x_s, x_g) = \int_0^\infty S(s; z, x_s, x_g) dG(s) }
where \eqn{G(s)} is the cumulative distribution function for the frailty.


}

\value{
  Returns a penetrance value evaluated at given age, gender and mutation status based on the specified penetrance model.
}
\author{
  Yun-Hee Choi
}
\seealso{
  \code{\link{penmodelEM}, \link{penmodel}, \link{penplot}, \link{print.penmodel},\link{summary.penmodel},}
  
  \code{\link{print.summary.penmodel}, \link{plot.penmodel}, \link{penci}} 
}
\examples{
  
  # Family data simulated from population-based design using a Weibull baseline hazard 
  
  fam <- simfam(N.fam=300, design="pop+", variation="none", base.dist="Weibull", 
                base.parms=c(0.01,3), vbeta=c(-1.13, 2.35), agemin=20, allelefreq=0.02)
  
  # Penetrance model fit for simulated family data
  
  fit <- penmodel(parms=c(0.01, 3), vbeta=c(-1.13, 2.35), data=fam, 
                  design="pop+", base.dist="Weibull")
  
  # Computing the penetrance at age 50 for male (sex=1) carriers (mut=1) from the assumed 
  # penetrance model based on Weibull baseline hazard.
  
  penf(fit$coefficients, age=50, sex=1, mut=1, base.dist="Weibull", agemin=20)
  
  
}
\keyword{Penetrance function}
