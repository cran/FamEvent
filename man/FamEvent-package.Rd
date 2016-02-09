\name{FamEvent-package}
\alias{FamEvent-package}
\docType{package}
\title{
Family age-at-disease onset data simulation and penetrance function estimation
}
\description{
This package provides functions to generate time-to-event data in the context of familial studies and to fit familial time-to-event data by correcting for ascertainment bias arising from a specified study design based on probands' mutation and disease statuses. The options on study designs are: \code{"pop"} for population-based design where families are ascertained through affected probands, \code{"pop+"} for population-based design with affected and mutation carrier probands, \code{"cli"} for clinic-based design that includes affected probands with at least one parent and one sib affected, \code{"cli+"} for clinic-based design that includes affected and mutation carrier probands with at least one parent and one sib affected, and \code{"twostage"} for two-stage design that randomly samples families from the population in the first stage and oversamples high risk families that include at least two affected members in the family in the second stage. 

For generating family data, ages at onset are generated specific to family members' gender and mutation status according to a specified study design and residual familial correlations are induced by either a shared frailty or a second gene. For fitting family data, an ascertainment corrected prospective likelihood approach is used to account for the study design and a composite likelihood approach is used for the two-stage sampling design. The Expectation and Maximization (EM) algorithm has been implemented for inferring missing genotypes conditional on observed genotypes and phenotypes in the families. For those who have missing genotypes, their carrier probabilities are obtained either from the fitted model or from Mendelian trasmission probabilities. 
This package also provides functions to plot the age-dependent penetrance curves estimated parametrically from the fitted model or non-parametrically from the data. 

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
\code{\link{simfam}, \link{summary.simfam}, \link{plot.simfam}, \link{penplot},}
	
\code{\link{penmodel}, \link{summary.penmodel}, \link{plot.penmodel}} 
}
\examples{

#  Simulate family data 

fam <- simfam(N.fam=100, design="pop+", variation="frailty", base.dist="Weibull", 
       frailty.dist="gamma", depend=1, base.parms=c(0.01, 3), vbeta=c(-1.13, 2.35), 
       allelefreq=0.02)

# summary of simulated family data

summary(fam) 

# Pedigree plots for family 1 and 2

plot(fam, famid=c(1,2))

# penetrance function plots given model parameter values for specified 
# baseline and frailty distributions

penplot(base.parms=c(0.01,3), vbeta=c(-1.3, 2.35), base.dist="Weibull", 
        frailty.dist="gamma", variation="frailty", depend=1, agemin=20)

# model fit of family data 

fit <- penmodel(parms=c(0.01, 3), vbeta=c(-1.13, 2.35), data=fam, 
       design="pop+", base.dist="Weibull")

# summary of estimated model parameters and penetrance estimates 

summary(fit)

# penetrance curves useful for model checking 

plot(fit) 

}
