\name{carrierprobpheno}
\alias{carrierprobpheno}
\title{
Mutation carrier probabilities for individuals with missing gentoypes conditional on observed phenotypes and genotypes 
}
\description{
This function computes model-based carrier probabilities for individuals with missing genotypes based on observed pheotypes and mutation status of family members using the fitted penetrance model.
}
\usage{
carrierprobpheno(fit, method="data", mode="dominant", q=0.02)
}
\arguments{
  \item{fit}{
An object of class \code{penmodel}, a fitted model by \code{penmodelEM} function for inferring missing mutation status in the family.
}
 \item{method}{
Choice of methods to calculate the carrier probability based on observed mutation status only.  Possible choices are \code{"data"} for empirical calculation of the carrier probabilities based on the observed carrier status in the entire sample, specific to generation and proband's mutation status or \code{"mendelian"} using Mendelian transmission probabilies based on observed mutation status in the family with a fixe allele frequency. Default is \code{"data"}.
 
 If \code{method="mendelian"}, specify both the \code{mode} of inheritance and the allele frequency \code{q}. 
 
}
 \item{mode}{
Choice of modes of inheritance. Possible choices are \code{"dominant"} for dominant model or \code{"recessive"} for recessive model. Default is \code{"dominant"}.
}

 \item{q}{
The frequency of the disease causing allele. The value should be between 0 and 1. If \code{NULL}, the estimated allele frequency from data will be used. Default value  is \code{0.02}.
}
}
\details{

The carrier probability for individual \eqn{i} conditional on the observed phenotype and carrier status of family members is calculated by 
\deqn{ P(X_i = 1 | Y_i, X^o ) = \frac{P(Y_i | X_i=1) P(X_i=1|X^o)}{P(Y_i | X_i=1) P(X_i=1|X^o) + P(Y_i| X_i=0) P(X_i=0|X^o)},  }
where \eqn{X_i} indicates the carrier status of individual \eqn{i} and \eqn{X^o} represents the observed carrier statuses in his or her family members; \eqn{Y_i} represents the observed phenotype \eqn{(t_i, \delta_i)}  of individual \eqn{i} in terms of age at onset \eqn{t_i} and disease status indicator \eqn{\delta_i} with 1 used for affected individuals and 0 for unaffected individuals. 

}
\value{
Returns a data frame with a vector of carrier probabilities called \code{carrp.pheno} added after the last column of the family data.
}
\author{
Yun-Hee Choi
}
\seealso{
\code{\link{simfam}, \link{penmodelEM}, \link{plot.simfam}}
}
\examples{

# Simulated family data with 30\% of members missing their genetic information.

fam <- simfam(N.fam=100, design="pop+", variation="secondgene", 
       base.dist="Weibull", allelefreq=c(0.02, 0.2), mrate=0.3,
       base.parms=c(0.01,3), vbeta=c(-1.13, 2.35, 0.5), agemin=20)
 
# EM algorithm for fitting family data with missing genotypes assuming a Weibull
# baseline hazard and dominant mode of Mendelian inheritance for a major gene.

fit <- penmodelEM(parms=c(0.01, 3), vbeta=c(-1.13, 2.35), data=fam, design="pop+",
       base.dist="Weibull", method="mendelian", mode="dominant")

# Carrier probability obtained by conditioning on the observed phenotypes, 
# assuming a dominant Mendelian mode of inheritance

fam.added <- carrierprobpheno(fit, method="mendelian", mode="dominant", q=0.02)

# pedigree plot for family 1 displaying carrier probabilities

plot.simfam(fam.added, famid=1)

}
\keyword{Carrier probability}
