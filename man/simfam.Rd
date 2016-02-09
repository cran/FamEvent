\name{simfam}
\alias{simfam}
\title{
Generate familial time-to-event data
}
\description{
This function generates familial time-to-event data for specified study design, genetic model and source of residual familial correlation; the generated data frame also contains family structure (individual's id, father id, mother id, relationship to proband, generation), gender, current age, genotypes of major or second genes.
}
\usage{
simfam(N.fam, design="pop", variation="none", depend=1, 
       base.dist="Weibull", frailty.dist="gamma", base.parms, vbeta, 
       allelefreq=c(0.02, 0.2), dominant.m=TRUE, dominant.s=TRUE,
       mrate=0, hr=0, age1=c(65,2.5), age2=c(45,2.5), agemin=20)
% \method{summary}{simulation}(object, \dots)
}

\arguments{
  \item{N.fam}{
Number of families to generate.
}
  \item{design}{
The family based study design used in the simulations. Possible choices are: \code{"pop", "pop+", "cli", "cli+"} or \code{"twostage"}, where \code{"pop"} is for the population-based design that families are ascertained by affected probands, \code{"pop+"} is similar to \code{"pop"} but with mutation carrier probands, \code{"cli"} is for the clinic-based design that includes affected probands with at least one parent and one sib affected, \code{"cli+"} is similar to \code{"cli"} but with mutation carrier probands and \code{"twostage"} for two-stage design that randomly samples families from the population in the first stage and oversamples high risk families in the second stage that include at least two affected members in the family. Default is \code{"pop"}.
}
  \item{variation}{
Source of residual familial correlation.  Possible choices are: \code{"frailty"} for frailty shared within families, \code{"secondgene"} for second gene variation, or \code{"none"} for no residual familial correlation. Default is \code{"none"}.
}

  \item{depend}{
Variance of the frailty distribution. Dependence within families increases with depend value.  Default value is 1. 
}

 \item{base.dist}{
Choice of baseline hazard distribution.  Possible choices are: \code{"Weibull"}, \code{"loglogistic"}, \code{"Gompertz"}, \code{"lognormal",} or \code{"gamma"}. Default is \code{"Weibull"}.
}

 \item{frailty.dist}{
Choice of frailty distribution.  Possible choices are: \code{"gamma"} for gamma distribution or \code{"lognormal"} for log normal distribution. Default is \code{"gamma"}.
}
  \item{base.parms}{
Vector of parameter values for baseline hazard function. 

\code{base.parms=c(lambda, rho)}, where lambda and rho are the shape and scale parameters, respectively.
}
  \item{vbeta}{
 Vector of parameter values for gender, majorgene, and secondgene. 
 }
  \item{allelefreq}{
Vector of population allele frequencies of major and second disease gene alleles. Frequencies must be between 0 and 1. Default frequencies are 0.02 for major gene allele and 0.2 for second gene allele, \code{allelefreq=c(0.02, 0.2)}
}
  \item{dominant.m}{logical; if TRUE, the genetic model of major gene is dominant, otherwise recessive.
}
  \item{dominant.s}{logical; if TRUE, the genetic model of second gene is dominant, otherwise recessive.
}
  \item{mrate}{
Proportion of missing genotypes, value between 0 and 1. Default value is 0.
}
  \item{hr}{
Proportion of high risk families, which include at least two affected members, to be sampled from the two stage sampling. This value should be specified when \code{design="twostage"} is used. Default value is 0. Value should lie between 0 and 1.
}
  \item{age1}{
Vector of mean and standard deviation for the current age of generation 1  or grandparents. Default values are mean of 65 years and standard deviation of 2.5 years, \code{age1=c(65,2.5)}.
}
  \item{age2}{
Vector of mean and standard deviation for the current age of generation 2 or proband generation. Default values are mean of 45 years and standard deviation of 2.5 years, \code{age2=c(45,2.5)}.
}
\item{agemin}{
Minimum age of disease onset. Default is 20 years of age.
}
}
\details{
The \code{design} argument defines the type of family based design to be simulated. Two variants of the population-based and clinic-based design can be chosen: \code{"pop"} when proband is affected, \code{"pop+"} when proband is affected mutation carrier, \code{"cli"} when proband is affected and at least one parent and one sibling are affected, \code{"cli+"} when proband is affected mutation-carrier and at least one parent and one sibling are affected. The two-stage design, \code{"twostage"}, is used to oversample high risk families, where the proportion of high risks families to include in the sample is specified by \code{hr}. High risk families often include multiple (at least two) affected members in the family.

Age at onset is generated from the penetrance model where residual familial correlation is induced by either a latent random variable called "frailty"" or a second gene shared by family members.

The penetrance model with a shared frailty model has the form
\deqn{ h(t|Z) = h_0(t-t_0) Z \exp(\beta_s x_s + \beta_{g1} x_{g1}) }
where \eqn{Z} represents a frailty shared within families and follows either a gamma or log-normal distribution; \eqn{t_0} is a minimum age of disease onset; \eqn{x_s} indicates males (1) and females (0) and \eqn{x_{g1}} indicates carriers (1) and non-carriers (0) of major gene mutation.

The penetrance model with a second gene variation has the form
\deqn{ h(t|Z) = h_0(t-t_0) \exp(\beta_s x_s + \beta_{g1} x_{g1} + \beta_{g2} x_{g2}) }
where \eqn{x_{g2}} indicates carriers (1) and non-carriers (0) of a second gene mutation.


The current ages for each generation are simulated assuming normal distributions. However, the probands' ages are generated using a left truncated normal distribution as their ages cannot be less than the minimum age of onset. The mean age difference between each generation and their parents is specified as at least 20 years apart. 

}
\value{
The function returns a data frame which contains:
\item{famID}{Family identification number (id).}
\item{indID}{Individual id.}
\item{gender}{Gender indicator: 1 for males, 0 for females.}
\item{motherID}{Mother id number.}
\item{fatherID}{Father id number.}
\item{proband}{Proband indicator: 1 if the individual is the proband, 0 otherwise.}
\item{generation}{Individuals generation:1=parents of probands,2=probands and siblings,3=children of probands and siblings.}
\item{majorgene}{Genotype of major gene: 1=AA, 2=Aa, 3=aa where A is disease gene.}
\item{secondgene}{Genotype of second gene: 1=BB, 2=Bb, 3=bb where B is disease gene.}
\item{ageonset}{Age at disease onset.}
\item{currentage}{Current age.}
\item{time}{Minimum time between current age and age at onset. }
\item{status}{Disease status: 1 for affected and 0 for unaffected (censored).}
\item{mgene}{Carrier status of major gene which can possibly be missing: 1 for carrier, 2 for non-carrier, NA for missing carrier status}
\item{relation}{Family members' relationship with the proband is as follows}
  \tabular{ll}{
    1  \tab Proband (self) \cr
    2  \tab Brother or sister \cr
    3  \tab Son or daughter \cr
    4  \tab Parent \cr
    5  \tab Nephew or niece \cr
    6  \tab Husband \cr
    7  \tab Brother or sister in law
  }

\item{fsize}{Family size including parents, siblings and children of the proband and the siblings.}
\item{naff}{Number of affected members in family.}
\item{weight}{Sampling weights.}


}
\references{
Choi, Y.-H., Kopciuk, K. and Briollais, L. (2008) Estimating Disease Risk Associated Mutated Genes in Family-Based Designs, Human Heredity 66, 238-251

Choi, Y.-H. and Briollais (2011) An EM Composite Likelihood Approach for Multistage Sampling of Family Data with Missing Genetic Covariates, Statistica Sinica 21, 231-253
}
\author{
Yun-Hee Choi, Wenqing He
}
\seealso{
\code{\link{summary.simfam}, \link{plot.simfam}, \link{penplot}}
}
\examples{

## Example 1: simulate family data from population-based design using
#  a Weibull distribution for the baseline hazard and inducing 
#  residual familial correlation through a shared gamma frailty.

fam <- simfam(N.fam=100, design="pop+", variation="frailty", 
       base.dist="Weibull", frailty.dist="gamma", depend=1, 
       allelefreq=0.02, base.parms=c(0.01,3), vbeta=c(-1.13, 2.35))

head(fam) 

#   famID indID gender motherID fatherID proband generation majorgene secondgene
# 1     1     1      1        0        0       0          1         2          0
# 2     1     2      0        0        0       0          1         3          0
# 3     1     3      0        2        1       1          2         2          0
# 4     1     4      1        0        0       0          0         3          0
# 5     1     7      0        3        4       0          3         2          0
# 6     1     8      1        3        4       0          3         3          0
#   ageonset currentage time status mgene relation fsize naff weight
# 1       70         68   68      0     1        4    11    1      1
# 2      110         68   68      0     0        4    11    1      1
# 3       36         40   36      1     1        1    11    1      1
# 4      212         50   50      0     0        6    11    1      1
# 5       79         19   19      0     1        3    11    1      1
# 6      169         16   16      0     0        3    11    1      1

summary(fam)

plot(fam, famid=c(1:2)) # pedigree plots for families with IDs=1 and 2

## Example 2: simulate family data from two stage design to include 
#  30\% of high risk families in the sample. 

fam <- simfam(N.fam=100, design="twostage", variation="frailty", 
       base.dist="Weibull", frailty.dist="gamma", depend=1, hr=0.3,
       base.parms=c(0.01,3), vbeta=c(-1.13, 2.35), allelefreq=0.02)

summary(fam)

}
\keyword{Simulation}