\name{penmodelEM}
\alias{penmodelEM}
\title{EM algorithm for estimating the penetrance model with missing genotypes}
\description{
Fits a penetrance model for family data with missing genotypes via the EM algorithm and provides model parameter estimates and corresponding gender- and genotype-specific penetrance estimates. 
}
\usage{
penmodelEM(parms, vbeta, data, design="pop", base.dist="Weibull", 
           robust=FALSE, method="data", mode="dominant", q=0.02)

}
\arguments{
  \item{parms}{
Vector of initial values for baseline parameters.
\code{parms=c(lambda, rho)}, where lambda and rho are the initial values for the scale and shape parameters, respectively. For the \code{"lognormal"} baseline distribution, rho > 0; for the other baselinse distributions, lambda > 0 and rho > 0.
}
\item{vbeta}{
 Vector of initial values for the regression coefficients for gender and majorgene,
 \code{vbeta=c(beta.s, beta.g)}.   
}
\item{data}{
  Data frame generated from \code{\link{simfam}} or data frame 
  containing specific variables: \code{famID}, \code{indID}, \code{generation}, \code{gender}, \code{currentage}, \code{mgene}, \code{time}, \code{status} and \code{weight} with \code{attr(data,"agemin")} specified.
}
\item{design}{
Study design of the family data. Possible choices are: \code{"pop"}, \code{"pop+"}, \code{"cli"}, \code{"cli+"} or \code{"twostage"}, where \code{"pop"} is for the population-based design with affected probands whose mutation status can be either carrier or non-carrier, \code{"pop+"} is similar to \code{"pop"} but with mutation carrier probands, \code{"cli"} is for the clinic-based design that includes affected probands with at least one parent and one sib affected, \code{"cli+"} is similar to \code{"cli"} but with mutation carrier probands, and \code{"twostage"} is for the two-stage design with oversampling of high risks families. Default is \code{"pop"}.
	}
 \item{base.dist}{
Choice of baseline hazard distribution to fit. Possible choices are: \code{"Weibull"}, \code{"loglogistic", "Gompertz", "lognormal",} or \code{"gamma"}. Default is \code{"Weibull"}.
}
 \item{robust}{
 	Logical; if TRUE, use robust 'sandwich' standard errors and variance matrix, otherwise use conventional standard errors and variance matrix.
 }

 \item{method}{
Choice of methods for calculating the carrier probabilities for individuals with missing mutation status. Possible choices are \code{"data"} for empirical calculation of the carrier probabilities based on the observed carriers' statuses in the entire sample, specific to generation and proband's mutation status or \code{"mendelian"} for calculating carrier probabilities based on Mendelian transmission probabilies with the given allele frequency and mutation statuses observed in the family. Default is \code{"data"}.
 
 If \code{method="mendelian"}, specify both \code{mode} of the inheritance and the allele frequency \code{q}. 
 
}
 \item{mode}{
Choice of modes of inheritance for calculating carrier probabilies for individuals with missing mutation status. Possible choices are \code{"dominant"} for dominant model or \code{"recessive"} for recessive model. Default is \code{"dominant"}.
}

 \item{q}{
Frequency of the disease causing allele used for calculating carrier pobabilities. The value should be between 0 and 1. If \code{NULL}, the estimated allele frequency from data will be used. Default value  is \code{0.02}.
}



}
\details{

The expectation and maximization (EM) algorithm is applied for making inference about the missing genotypes. In the expectation step, for individuals with unknown carrier status, we first compute their carrier probabilities given their family's observed phenotype and genotype information based on current estimates of parameters \eqn{\theta}
\deqn{ w_{fi} = P(X_{fi}=1|Y_{fi}, X^o_f) ,}
where  \eqn{X_{fi}} represents the mutation carrier status and \eqn{Y_{fi}} represents the phenotype (\eqn{t_{fi}, \delta_{fi}}) in terms of age at onset \eqn{t_{fi}} and disease status \eqn{\delta_{fi}} for individual \eqn{i} in family \eqn{f} and \eqn{X^o_f} represents the observed genotypes in family \eqn{f}.

 
Then, we obtain the conditional expectation of the log-likelihood function of the complete data given the observed data as a weighted log-likelihood, which has the form
\deqn{ E_{\theta} [\ell (\theta) | Y, X^o)] =  \sum_f^n \sum_i^{n_f} \ell_{fi}(\theta | X_{fi}=1) w_{fi} +  \ell_{fi}(\theta | X_{fi}=0) (1-w_{fi}), }

In the maximization step, the updated parameter estimates are obtained by maximizing the weighted log likelihood computed in the E-step.

These expectation and maximization steps iterate until convergence to obtain the maximum likelihood estimates.

See more details in Choi and Briollais (2011) or Choi et al. (2014).

Transformed baseline parameters (\eqn{\lambda, \rho}) were used for estimation; see \code{\link{penmodel}} for details.
}
\value{
Returns an object of class \code{'penmodel'}, including the following elements:

\item{coefficients}{
Parameter estimates of transformed baseline parameters (\eqn{\lambda, \rho}) and regression coefficients for gender and mutation status (\eqn{\beta_s, \beta_g}) including their standard errors and also their robust standard errors.
}
\item{varcov}{
Variance covariance matrix of parameter estimates. If \code{robust=TRUE}, robust `sandwich' variance covariance matrix is returned.
}
\item{se}{
Standard errors of parameter estimates. If \code{robust=TRUE}, robust `sandwich' standard errors are returned.
}
\item{pen70.est}{
 Penetrance estimates by age 70 specific to gender and mutation-status subgroups.
}
\item{pen70.se}{
 Standard errors of penetrance estimates by age 70 specific to gender and mutation-status subgroups. 
}
\item{pen70.ci}{
 95\% confidence interval estimates of penetrance by age 70 specific to gender and mutation-status subgroups. 
}

\item{ageonset}{
	Vector of ages of onset ranging from \code{agemin} to 90 years.
}
\item{pen.maleCarr}{Vector of penetrance estimates for male carriers from \code{agemin} to 90 years.}
\item{pen.femaleCarr}{Vector of penetrance estimates for female carriers from \code{agemin} to 90 years.}
\item{pen.maleNoncarr}{Vector of penetrance estimates for male non-carriers from \code{agemin} to 90 years.}
\item{pen.femaleNoncarr}{Vector of penetrance estimates for female non-carriers from \code{agemin} to 90 years.}
\item{logLik}{Loglikelihood value for the fitted penetrance model.}
}
\references{
Choi, Y.-H. and Briollais, L. (2011) An EM composite likelihood approach for multistage sampling of family data with missing genetic covariates, Statistica Sinica 21, 231-253.

Choi, Y.-H., Briollais, L., Green, J., Parfrey, P., and Kopciuk, K. (2014) Estimating successive cancer risks in Lynch Syndrome families using a progressive three-state model, Statistics in Medicine 33, 618-638.
}
\author{
Yun-Hee Choi
}
\seealso{
\code{\link{simfam}}, \code{\link{penmodel}},  \code{\link{print.penmodel}}, \code{\link{summary.penmodel}},  \code{\link{print.summary.penmodel}}, 

\code{\link{plot.penmodel}},\code{\link{carrierprob}} 
}

\examples{

# Family data simulated with 30\% of members missing their genetic information.

fam <- simfam(N.fam=100, design="pop+", base.dist="Weibull", base.parms=c(0.01,3),
       vbeta=c(-1.13, 2.35), agemin=20, allelefreq=0.02, mrate=0.3)
 
# EM algorithm for fitting family data with missing genotypes

fit <- penmodelEM(parms=c(0.01, 3), vbeta=c(-1.13, 2.35), data=fam, design="pop+",
       base.dist="Weibull", method="mendelian", mode="dominant", q=NULL)

# Summary of the model parameter and penetrance estimates from model fit 
# by penmodelEM 

summary(fit)

# Generate the lifetime penetrance curves from model fit for gender and 
# mutation status groups along with their non-parametric penetrance curves 
# based on observed data
 
plot(fit)

}
\keyword{Penetrance model}
