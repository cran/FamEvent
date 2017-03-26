\name{penmodel}
\alias{penmodel}
\title{Estimate the penetrance model and penetrance curves}
\description{
Fits a penetrance model for family data based on a prospective likelihood with ascertainment correction and provides parameter estimates as well as the gender- and mutation-specific penetrance estimates. 
}
\usage{
penmodel(parms, vbeta, data, design="pop", base.dist="Weibull", robust=FALSE)

}
\arguments{
  \item{parms}{
Vector of initial values for baseline parameters.
\code{parms=c(lambda, rho)}, where lambda and rho are the initial values for the scale and shape parameters, respectively. For the \code{"lognormal"} baseline distribution, \code{rho} > 0; for the other baselinse distributions, \code{lambda} > 0 and \code{rho} > 0.
}
\item{vbeta}{
 Vector of initial values for regression coefficients for gender and majorgene;
 \code{vbeta=c(beta.s, beta.g)}.   
}
\item{data}{
  Data frame generated from \code{\link{simfam}} or data frame 
  containing specific variables: \code{famID}, \code{indID}, \code{generation}, \code{gender}, \code{currentage}, \code{mgene}, \code{time}, \code{status} and \code{weight} with \code{attr(data,"agemin")} specified.
}
\item{design}{
Study design of the family data.  Possible choices are: \code{"pop", "pop+", "cli"}, \code{"cli+"} or \code{"twostage"}, where \code{"pop"} is for the population-based design with affected probands whose mutation status can be either carrier or non-carrier, \code{"pop+"} is similar to \code{"pop"} but with mutation carrier probands, \code{"cli"} is for the clinic-based design that includes affected probands with at least one parent and one sib affected, \code{"cli+"} is similar to \code{"cli"} but with mutation carrier probands, and \code{"twostage"} is for the two-stage design with oversampling of high risks families. Default is \code{"pop"}.
	}
 \item{base.dist}{
Choice of baseline hazard distribution to fit. Possible choices are: \code{"Weibull"}, \code{"loglogistic"}, \code{"Gompertz"}, \code{"lognormal",} or \code{"gamma"}. Default is \code{"Weibull"}.
}
 \item{robust}{
 	Logical; if TRUE, use robust `sandwich' standard errors and variance covariance matrix, otherwise use conventional standard errors and variance covariance matrix.
 }

}
\details{

The penetrance model is fitted to family data with a specified baseline hazard distribution,
\deqn{ h(t|x_s, x_g) = h_0(t) \exp(\beta_s x_s+\beta_g x_g) }
where \eqn{h_0(t)} is the baseline hazards function specified by \code{base.dist}, which depends on the shape and scale parameters, \eqn{\lambda} and \eqn{\rho}; \eqn{x_s} indicates male (1) and female (0) and \eqn{x_g} indicates carrier (1) or non-carrier (0) of a gene of interest (major gene).

For family data arising from population- or clinic-based study designs (\code{design="pop", "pop+"}, \code{"cli"}, or \code{"cli+"}), the parameters of the penetrance model are estimated from the ascertainment-corrected prospective likelihood approach (Choi, Kopciuk and Briollais, 2008).

For family data arising from a two-stage study design (\code{design="twostage"}), model parameters are estimated based on the composite likelihood approach (Choi and Briollais, 2011)

Transformed baseline parameters (\eqn{\lambda, \rho}) were used for estimation; log tranformation was applied to both scale and shape parameters for \code{"Weibull"}, \code{"loglogistic"}, \code{"Gompertz"} and \code{"gamma"} baseline distributions. For \code{"lognormal"} baseline distribution, the log transformation was applied only to shape parameter \eqn{\rho}, not to \eqn{\lambda} which represents the location parameter in log-normal distribution.

Calculations of standard errors and 95\% confidence intervals for penetrance estimates by age 70 were based on the penetrances obtained from 1000 Monte-Carlo simulations of the estimated penetrance model; for more details, see \code{\link{penci}}.

}
\value{
Returns an object of class \code{'penmodel'}, including the following elements:

\item{coefficients}{
Parameter estimates of transformed baseline parameters (\eqn{\lambda, \rho}) and regression coefficients for gender and mutation status (\eqn{\beta_s, \beta_g}).
}
\item{varcov}{
Variance covariance matrix of parameter estimates. If \code{robust=TRUE}, robust `sandwich' variance covariance matrix is returned.
}
\item{se}{
Standard errors of parameter estimates. If \code{robust=TRUE}, robust 'sandwich' standard errors are returned.
}
\item{pen70.est}{
 Penetrance estimates by age 70 specific to gender and mutation-status subgroups.
}
\item{pen70.se}{
 Standard errors of penetrance estimates by age 70  specific to gender and mutation-status subgroups. 
}
\item{pen70.ci}{
 95\% confidence interval for penetrance estimates by age 70 specific to gender and mutation-status subgroups. 
}

\item{ageonset}{
	Vector of ages of onset ranging from \code{agemin} to 90 years.
}
\item{pen.maleCarr}{Vector of penetrance estimates for male carriers from \code{agemin} to 90 years.}
\item{pen.femaleCarr}{Vector of penetrance estimates for female carriers from \code{agemin} to 90 years.}
\item{pen.maleNonCarr}{Vector of penetrance estimates for male non-carriers from \code{agemin} to 90 years.}
\item{pen.femaleNonCarr}{Vector of penetrance estimates for female non-carriers from \code{agemin} to 90 years.}
\item{logLik}{Loglikelihood value for the fitted penetrance model.}
}
\references{
Choi, Y.-H., Kopciuk, K. and Briollais, L. (2008) Estimating Disease Risk Associated Mutated Genes in Family-Based Designs, Human Heredity 66, 238-251

Choi, Y.-H. and Briollais (2011) An EM Composite Likelihood Approach for Multistage Sampling of Family Data with Missing Genetic Covariates, Statistica Sinica 21, 231-253
}
\author{
Yun-Hee Choi
}
\seealso{
\code{\link{penmodelEM}, \link{simfam}, \link{penplot}, \link{print.penmodel},\link{summary.penmodel},}

\code{\link{print.summary.penmodel}, \link{plot.penmodel}, \link{penci}, \link{penf}} 
}
\examples{

# Family data simulated from population-based design using a Weibull baseline hazard 

fam <- simfam(N.fam=300, design="pop+", variation="none", base.dist="Weibull", 
       base.parms=c(0.01,3), vbeta=c(-1.13, 2.35), agemin=20, allelefreq=0.02)
 
# Penetrance model fit for simulated family data

fit <- penmodel(parms=c(0.01, 3), vbeta=c(-1.13, 2.35), data=fam, 
       design="pop+", base.dist="Weibull")

# Summary of the model parameter and penetrance estimates from model fit

summary(fit)

# Generate the lifetime penetrance curves from model fit for specific gender and 
# mutation status groups along with their non-parametric penetrance curves 
     
plot(fit)

}
\keyword{Penetrance model}
