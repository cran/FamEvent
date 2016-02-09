\name{penmodel}
\alias{penmodel}
\title{Estimate the penetrance model and penetrance curves}
\description{
This function fits the penetrance model to the family data based on a prospective likelihood with ascertainment correction and provides parameter estimates as well as the gender- and mutation-specific penetrance estimates. 
}
\usage{
penmodel(parms, vbeta, data, design="pop", base.dist="Weibull")

}
\arguments{
  \item{parms}{
Vector of initial values for baseline parameters.
\code{parms=c(lambda, rho)}, where lambda and rho are the initial values for the scale and shape parameters, respectively; lambda > 0, rho > 0.
}
\item{vbeta}{
 Vector of initial values for regression coefficients for gender and majorgene.
 \code{vbeta=c(beta.s, beta.g)}.   
}
\item{data}{
Family data that contains specific variables famID, indID, generation, gender, currentage, mgene, time, status, weight; data structure should follow the format of the data generated from \code{\link{simfam}}.
}
\item{design}{
Study design of the family data.  Possible choices are: \code{"pop", "pop+", "cli"}, \code{"cli+"} or \code{"twostage"}, where \code{"pop"} is for the population-based design with affected probands whose mutation status can be either carrier or non-carrier, \code{"pop+"} is similar to \code{"pop"} but with mutation carrier probands, \code{"cli"} is for the clinic-based design that includes affected probands with at least one parent and one sib affected, \code{"cli+"} is similar to \code{"cli"} but with mutation carrier probands, and \code{"twostage"} is for the two-stage design with oversampling of high risks families. Default is \code{"pop"}.
	}
 \item{base.dist}{
Choice of baseline hazard distribution to fit. Possible choices are: \code{"Weibull"}, \code{"loglogistic"}, \code{"Gompertz"}, \code{"lognormal",} or \code{"gamma"}. Default is \code{"Weibull"}.
}

}
\details{

The penetrance model is fitted to family data with a specified baseline hazard distribution,
\deqn{ h(t|x_s, x_g) = h_0(t) \exp(\beta_s x_s+\beta_g x_g) }
where \eqn{h_0(t)} is the baseline hazards function specified by \code{base.dist}, which depends on the shape and scale parameters, \eqn{\lambda} and \eqn{\rho}; \eqn{x_s} indicates male (1) and female (0) and \eqn{x_g} indicates carrier (1) or non-carrier (0) of a gene of interest (major gene).

For family data arising from population- or clinic-based study designs (\code{design="pop", "pop+"}, \code{"cli"}, or \code{"cli+"}), the parameters of the penetrance model are estimated from the ascertainment-corrected prospective likelihood approach (Choi, Kopciuk and Briollais, 2008).

For family data arising from a two-stage study design (\code{design="twostage"}), model parameters are estimated based on the composite likelihood approach (Choi and Briollais, 2011)

}
\value{
An object of class \code{penmodel}, a list including elements

\item{parms.est}{
Parameter estimates of baseline parameters (\eqn{\lambda, \rho}) and regression coefficients for gender and mutation status (\eqn{\beta_s, \beta_g}) including their standard errors and also robust standard errors.
}
\item{parms.cov}{
Covariance matrix of parameter estimates.
}
\item{parms.se}{
Standard errors of parameter estimates.
}
\item{parms.rcov}{
Robust (sandwich) covariance matrix of parameter estimates.
}
\item{parms.rse}{
Robust standard errors of parameter estimates.
}

\item{pen70.est}{
 Penetrance estimates by age 70 specific to gender and mutation-status subgroups.
}
\item{pen70.se}{
 Standard errors of penetrance estimates by age 70  specific to gender and mutation-status subgroups. 
}
\item{pen70.ci}{
 95\% confidence interval for penetrance estimates by age 70 specific to gender and mutation-status groups. 
}

\item{ageonset}{
	Vector of ages of onset ranging from \code{agemin} to 80 years.
}
\item{pen.maleCarr}{Vector of penetrance estimates for male carriers from \code{agemin} to 80 years.}
\item{pen.femaleCarr}{Vector of penetrance estimates for female carriers from \code{agemin} to 80 years.}
\item{pen.maleNonCarr}{Vector of penetrance estimates for male non-carriers from \code{agemin} to 80 years.}
\item{pen.femaleNonCarr}{Vector of penetrance estimates for female non-carriers from \code{agemin} to 80 years.}
}
\references{
Choi, Y.-H., Kopciuk, K. and Briollais, L. (2008) Estimating Disease Risk Associated Mutated Genes in Family-Based Designs, Human Heredity 66, 238-251

Choi, Y.-H. and Briollais (2011) An EM Composite Likelihood Approach for Multistage Sampling of Family Data with Missing Genetic Covariates, Statistica Sinica 21, 231-253
}
\author{
Yun-Hee Choi
}
\seealso{
\code{\link{simfam}, \link{penplot}, \link{summary.penmodel}, \link{plot.penmodel}} 
}
\examples{

# Family data simulated from population-based design using a Weibull baseline hazard 
# and inducing residual familial correlation by a shared second gene

fam <- simfam(N.fam=300, design="pop+", variation="secondgene", 
       base.dist="Weibull", allelefreq=c(0.02, 0.2), 
       base.parms=c(0.01,3), vbeta=c(-1.13, 2.35, 0.5), agemin=20)
 
# Penetrance model fit for simulated family data

fit <- penmodel(parms=c(0.01, 3), vbeta=c(-1.13, 2.35), data=fam, 
       design="pop+", base.dist="Weibull")

# Summary of the model parameter and penetrance estimates from model fit

summary(fit)

# Generate the lifetime penetrance curves from model fit for specific gender and 
# mutation status groups along with their non-parametric penetrance curves 
     
plot(fit)

}
\keyword{Penetrance}
