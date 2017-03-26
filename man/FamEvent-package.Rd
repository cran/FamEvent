\name{FamEvent-package}
\alias{FamEvent-package}
\docType{package}
\title{
Family age-at-onset data simulation and penetrance estimation
}
\description{

Family-based studies are used to characterize the disease risk associated
with being a carrier of a major gene. When the disease risk can vary with 
age of onset, penetrance or disease risk functions need to provide age-dependent
estimates of this disease risk over lifetime.
This FamEvent package can generate age-at-onset data in the context of familial studies,
with correction for ascertainment (selection) bias arising from a specified
study design based on proband's mutation and disease statuses.
Possible study designs are: \code{"pop"} for population-based design where
families are ascertained through affected probands, \code{"pop+"} are similar to
\code{"pop"} but probands are also known mutation carriers, \code{"cli"} for
clinic-based design that includes affected probands with at least one
parent and one sib affected, \code{"cli+"} are similar to \code{"cli"} but probands are
also known mutation carrriers. And \code{"twostage"} for two-stage design that randomly samples families from
the population in the first stage and oversamples high risk families that includes at least two
affected members in the family at the second stage.

Ages at disease onset are generated specific to family members' gender and
mutation status according to the specified study design with residual
familial correlations induced by either a shared frailty or a second gene.
For estimating age at onset risks with family data, an ascertainment
corrected prospective likelihood approach is used to account for the
population or clinic-based study designs while a composite likelihood
approach is used for the two-stage sampling design. The Expectation and
Maximization (EM) algorithm has been implemented for inferring missing
genotypes conditional on observed genotypes and phenotypes in the families. For family members who have missing
genotypes, their carrier probabilities are obtained either from the fitted
model or from Mendelian transmission probabilities. This package also
provides functions to plot the age-dependent penetrance curves estimated
parametrically from the fitted model or non-parametrically from the data,
pedigree plots of simulated families and penetrance function curves for
carriers and non-carriers of a major and second gene based on specified
parameter values.
%\packageDescription{FamilyData}
}
\author{
Yun-Hee Choi, Karen Kopciuk, Laurent Briollais, Wenqing He

Maintainer: Yun-Hee Choi < yun-hee.choi@schulich.uwo.ca >
}
\references{
Choi, Y.-H., Kopciuk, K. and Briollais, L. (2008) Estimating Disease Risk Associated Mutated Genes in Family-Based Designs, Human Heredity 66, 238-251

Choi, Y.-H. and Briollais (2011) An EM Composite Likelihood Approach for Multistage Sampling of Family Data with Missing Genetic Covariates, Statistica Sinica 21, 231-253
}
\keyword{ package }
\seealso{
\code{\link{simfam}, \link{summary.simfam}, \link{plot.simfam}, \link{penplot},\link{carrierprob},}
	
\code{\link{penmodel},  \link{penmodelEM}, \link{print.penmodel}, \link{summary.penmodel},}

\code{\link{print.summary.penmodel}, \link{plot.penmodel}} 


}
\examples{

#  Simulate family data 

fam <- simfam(N.fam=100, design="pop+", variation="none", base.dist="Weibull", 
       base.parms=c(0.01, 3), vbeta=c(-1.13, 2.35), allelefreq=0.02)

# summary of simulated family data

summary(fam) 

# Pedigree plots for family 1 and 2

plot(fam, famid=c(1,2))

# penetrance function plots given model parameter values for Weibull baseline 

penplot(base.parms=c(0.01,3), vbeta=c(-1.3, 2.35), base.dist="Weibull", 
        variation="none", agemin=20)

# model fit of family data 

fit <- penmodel(parms=c(0.01, 3), vbeta=c(-1.13, 2.35), data=fam, 
       design="pop+", base.dist="Weibull")

# summary of estimated model parameters and penetrance estimates 

summary(fit)

# penetrance curves useful for model checking 

plot(fit) 

}
