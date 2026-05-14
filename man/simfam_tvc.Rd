\name{simfam_tvc}
\alias{simfam_tvc}
\title{
Generate familial time-to-event data with a time-varying covariate
}
\description{
Generates familial time-to-event data with a time-varying covariate for specified study design, genetic model and source of residual familial correlation; the generated data frame also contains family structure (individual's id, father id, mother id, relationship to proband, generation), gender, current age, genotypes of major or second genes.
}
\usage{
simfam_tvc(N.fam, design = "pop", variation = "none", interaction = FALSE, 
       add.x = FALSE, x.dist = NULL, x.parms = NULL, depend = NULL, 
       add.tvc = FALSE, tvc.type = "PE", tvc.range = NULL, tvc.parms = 1,
       base.dist = "Weibull", frailty.dist = NULL, base.parms = c(0.016, 3), 
       vbeta = c(1, 1), 
       allelefreq = 0.02, dominant.m = TRUE, dominant.s = TRUE,
       mrate = 0, hr = 0, probandage = c(45, 2), agemin = 20, agemax = 100)
% \method{summary}{simulation}(object, \dots)
}

\arguments{
  \item{N.fam}{
Number of families to generate.
}

\item{design}{
Family-based study design used in the simulations. Possible choices are \code{"pop"}, \code{"pop+"}, \code{"cli"}, \code{"cli+"}, \code{"twostage"}, and \code{"noasc"}.  

\code{"pop"} is a population-based design with affected probands; \code{"pop+"} is similar to \code{"pop"} but with mutation-carrying probands; \code{"cli"} is a clinic-based design with affected probands and at least one affected parent and one affected sibling; \code{"cli+"} is similar to \code{"cli"} but with mutation-carrying probands; \code{"twostage"} is a two-stage design that oversamples high-risk families (with at least two affected members); and \code{"noasc"} corresponds to simple random sampling without ascertainment correction. Default is \code{"pop"}.
}

  \item{variation}{
Source of residual familial correlation.  Possible choices are \code{"frailty"} for frailty shared within families, \code{"secondgene"} for a second gene variation, \code{"kinship"} for correlated frailties based on a kinship matrix, and \code{"none"} for no residual familial correlation. Default is \code{"none"}.
}
 \item{interaction}{
Logical; if \code{TRUE}, allows the interaction between gender and mutation status. Default is \code{FALSE}.
}
 \item{add.x}{
Logical; if \code{TRUE}, generates a covariate in addition to gender and majorgene. Default is \code{FALSE}. 
 }
\item{x.dist}{
Distribution of the covairate. Possible choices to generate the covariate are \code{"normal"} from normal distribution and \code{"binomial"} from binominal distribution. Default is \code{NULL}.
 }
\item{x.parms}{
Parameter values for the specified distribution of the covariate. \code{x.parms = c(mean, sd)} should be specified for \code{x.dist = "normal"} and \code{x.parms = c(size, probs)} should be specified for \code{"binomial"}; when \code{size = 1}, it generates binary covariate. Default is \code{NULL}.
}
 \item{depend}{
Inverse of variance of the frailty distribution. Dependence within families decreases with depend value.  Default is \code{NULL}. Value should be specified as a positive real number when \code{variation = "frailty"} or \code{variation = "kinship"}. 
}
\item{add.tvc}{
Logical; if \code{TRUE}, generates a time varying covariate. Default is \code{FALSE}. 
}
\item{tvc.type}{
Choice of time-varying covariate model. Possible choices are \code{"PE"} and \code{"CO"}.
Default is \code{"PE"}.
}
\item{tvc.range}{
Range of ages at which the time-varying covariate occurs. 
Default is \code{NULL}.
}
\item{tvc.parms}{
Vector of parameter values used for the time-varying covariate model. 
Default value is 1.
}
\item{base.dist}{
Choice of baseline hazard distribution. Possible choices are \code{"Weibull"}, \code{"loglogistic"}, \code{"Gompertz"}, \code{"lognormal"} \code{"gamma"}, and \code{"logBurr"}. Default is \code{"Weibull"}.
}

 \item{frailty.dist}{
Choice of frailty distribution.  Possible choices are \code{"gamma"} and \code{"lognormal"} when \code{variation = "frailty"}. Default is \code{NULL}.
}
  \item{base.parms}{
Vector of parameter values for the specified baseline hazard function. \code{base.parms = c(lambda, rho)} should be specified for \code{base.dist = "Weibull"}, \code{"loglogistic"}, \code{"Gompertz"}, \code{"gamma"}, and \code{"lognormal"}. For \code{base.dist = "logBurr"}, three parameters should be specified \code{base.parms = c(lambda, rho, eta)}. 
Default value is \code{base.parms = c(0.016, 3)} for \code{base.dist = "Weibull"}.
}
  \item{vbeta}{
Vector of regression coefficients for gender, majorgene, interaction between gender and majorgene (if \code{interaction = TRUE}), secondgene (if \code{variation = "secondgene"}), and additional covariate (if \code{add.x = TRUE}).}

 \item{allelefreq}{
Population allele frequencies of major disease gene. Value should be between 0 and 1. 
Vector of population allele frequencies of major and second disease genes should be provided when \code{variation = "secondgene"}. Default value is \code{allelefreq = 0.02}.}
  \item{dominant.m}{Logical; if \code{TRUE}, the genetic model of major gene is dominant, otherwise recessive.}
  \item{dominant.s}{Logical; if \code{TRUE}, the genetic model of second gene is dominant, otherwise recessive.}
  \item{mrate}{
Proportion of missing genotypes, value between 0 and 1. Default value is 0.}
  \item{hr}{
Proportion of high risk families, which include at least two affected members, to be sampled from the two stage sampling. This value should be specified when \code{design = "twostage"}. Default value is 0. Value should lie between 0 and 1.
}
\item{probandage}{
Vector of mean and standard deviation for the proband age. Default values are mean of 45 years and standard deviation of 2 years, \code{probandage = c(45, 2)}.
}
\item{agemin}{
Minimum age of disease onset or minimum age. Default is 20 years of age.
}
\item{agemax}{
Maximum age of disease onset or maximum age. Default is 100 years of age.
}

}
\details{

\bold{Time-varying covariate}

When \code{add.tvc = TRUE}, the time at which the time-varying covariate (TVC) occurs, \code{tvc.age}, is generated from a uniform distribution with the range specified by \code{tvc.range}. A vector of minimum and maximum ages for the TVC should be specified in \code{tve.range}. When \code{tvc.range = NULL}, \code{agemin} and \code{agemax} are used as the range. In addition, \code{tvc.type} should be either \code{"PE"} or \code{"CO"} and the parameter values for the specified TVC type should be provided in \code{tvc.parms}.

\code{tvc.type = "PE"} represents a permanent exposure model for TVC which assumes that the effect of the TVC stays constant after \code{tvc.age}. The \code{tvc.parms} for the PE model should be specified as a single value, which represents log hazard ratio.

\code{tvc.type = "CO"} represents the Cox and Oaks model for TVC which assumes that the effect of the TVC decays exponentially over time in the form
\eqn{\beta \exp(-(t - t^\ast)\eta ) + \eta_0,}{beta*exp( -(t - t*)*eta ) + eta0,}
where \eqn{t^\ast}{t*} is the time at which the TVC occurs.
The \code{tvc.parms} for the CO model should be specified by a vector of three parameters consisting of \code{c(beta, eta, eta0)}.

\bold{Family-based study design}

The \code{design} argument specifies the type of family based design to be simulated. Two variants of the population-based and clinic-based design can be chosen: \code{"pop"}, where the proband is affected; \code{"pop+"}, where the proband is affected and a mutation carrier; \code{"cli"}, where the proband is affected and at least one parent and one sibling are also affected; \code{"cli+"}, where the proband is affected, is a mutation carrier, and has at least one affected parent and one affected sibling. 
The two-stage design, \code{"twostage"}, is used to oversample high risk families, where the proportion of high risks families included in the sample is specified by \code{hr}. High risk families often include multiple (at least two) affected members. Finally, \code{design = "noasc"} specifies a design with no ascertainment criteria applied.


\bold{Penetrance model}

The ages at onset are generated from the following penetrance models depending on the choice of \code{variation = "none", "frailty", "secondgene", "kinship".}. When \code{variation = "none"}, the ages at onset are independently generated from the proportional hazard model conditional on the gender and carrier status of major gene mutation, \eqn{X = (x_s, x_g)}{X = c(xs, xg)}.  

The ages at onset correlated within families are generated from the shared frailty model (\code{variation = "frailty"}), the correlated shared frailty model with kinship matrix (\code{variation = "kinship"}), or the two-gene model (\code{variation = "secondene"}), where the residual familial correlation is induced by a frailty or a second gene, respectively, shared within the family.

\emph{The proportional hazard model}

\deqn{ h(t|X) = h_0(t - t_0) \exp(\beta_s x_s+\beta_g x_g),}{ h(g|X) = h0(t - t0) * exp(beta.s * xs + beta.g * xg),}

where \eqn{h_0(t)}{h0(t)} is the baseline hazard function, \eqn{t_0}{t0} is a minimum age of disease onset, \eqn{x_s}{xs} and \eqn{x_g}{xg} indicate male (1) or female (0) and carrier (1) or non-carrier (0) of a main gene of interest, respectively.


\emph{The shared frailty model}

\deqn{ h(t|X,Z) = h_0(t - t_0) Z \exp(\beta_s x_s+\beta_g x_g),}{ h(g|X,Z) = h0(t - t0) * Z * exp(\betas * xs + \betag * xg),}

where \eqn{h_0(t)}{h0(t)} is the baseline hazard function, \eqn{t_0}{t0} is a minimum age of disease onset, \eqn{Z} represents a frailty shared within families and follows either a gamma or log-normal distribution, \eqn{x_s}{xs} and \eqn{x_g}{xg} indicate male (1) or female (0) and carrier (1) or non-carrier (0) of a main gene of interest, respectively.

\emph{The correlated shared frailty model with kinship matrix}

\deqn{ h(t | X, Z) = h_0(t - t_0) Z \exp(\beta_s x_s+\beta_g x_g),}{ h(t | X, Z) = h0(t - t0) * Z * exp(\betas * xs + \betag * xg),}

where \eqn{h_0(t)}{h0(t)} is the baseline hazard function, \eqn{t_0}{t0} is a minimum age of disease onset, \eqn{Z}{Z} represents a vector of frailties following a multivariate log-normal distribution with mean 0 and variance \eqn{2*K*depend}{2*K*depend}, where \eqn{K} represents the kinship matrix, \eqn{x_s}{xs} and \eqn{x_g}{xg} indicate male (1) or female (0) and carrier (1) or non-carrier (0) of a main gene of interest, respectively.

\emph{The two-gene model}

\deqn{ h(t|X,Z) = h_0(t - t_0) Z \exp(\beta_s x_s+\beta_1 x_1 + \beta_2 x_2),}{ h(g|X,Z) = h0(t - t0) * Z * exp(betas * xs + beta1 * x1 + beta2 * x2),}

where \eqn{x_1, x_2}{x1,x2} indicate carriers (1) and non-carriers (0) of a major gene and of second gene mutation, respectively.

The current ages for each generation are simulated assuming normal distributions. However, the probands' ages are generated using a left truncated normal distribution as their ages cannot be less than the minimum age of onset. The average age difference between each generation and their parents is specified as 20 years apart. 

}
\value{
Returns an object of class \code{'simfam'}, a data frame which contains:
\item{famID}{Family identification (ID) numbers.}
\item{indID}{Individual ID numbers.}
\item{gender}{Gender indicators: 1 for males, 0 for females.}
\item{motherID}{Mother ID numbers.}
\item{fatherID}{Father ID numbers.}
\item{proband}{Proband indicators: 1 if the individual is the proband, 0 otherwise.}
\item{generation}{Individuals generation: 1 = parents of probands, 2 = probands and siblings, 3 = children of probands and siblings.}
\item{majorgene}{Genotypes of major gene: 1 = AA, 2 = Aa, 3 = aa where A is disease gene.}
\item{secondgene}{Genotypes of second gene: 1 = BB, 2 = Bb, 3 = bb where B is disease gene.}
\item{ageonset}{Ages at disease onset in years.}
\item{currentage}{Current ages in years.}
\item{time}{Ages at disease onset for the affected or ages of last follow-up for the unaffected. }
\item{status}{Disease statuses: 1 for affected, 0 for unaffected (censored).}
\item{mgene}{Major gene mutation indicators: 1 for mutated gene carriers, 0 for mutated gene noncarriers, or \code{NA} if missing.}
\item{newx}{Additional covariate when \code{add.x = TRUE}.}
\item{tvc.age}{Age at which the time-varying covariate occurs when \code{add.tvc = TRUE}.}
\item{tvc.status}{TVC status: 1 if \code{tvc.age < time}, 0 otherwise.}
\item{relation}{Family members' relationship with the proband:
  \tabular{ll}{
    1  \tab Proband (self) \cr
    2  \tab Brother or sister \cr
    3  \tab Son or daughter \cr
    4  \tab Parent \cr
    5  \tab Nephew or niece \cr
    6  \tab Spouse \cr
    7  \tab Brother or sister in law
  }
  }
\item{fsize}{Family size including parents, siblings and children of the proband and the siblings.}
\item{naff}{Number of affected members in family.}
\item{weight}{Sampling weights.}


}
\references{
Choi, Y.-H., Briollais, L., He, W. and Kopciuk, K. (2021) FamEvent: An R Package for Generating and Modeling Time-to-Event Data in Family Designs, Journal of Statistical Software 97 (7), 1-30. doi:10.18637/jss.v097.i07

Choi, Y.-H., Kopciuk, K. and Briollais, L. (2008) Estimating Disease Risk Associated Mutated Genes in Family-Based Designs, Human Heredity 66, 238-251.

Choi, Y.-H. and Briollais (2011) An EM Composite Likelihood Approach for Multistage Sampling of Family Data with Missing Genetic Covariates, Statistica Sinica 21, 231-253.
}
\author{
Yun-Hee Choi
}
\seealso{
\code{\link{summary.simfam_tvc}, \link{plot.simfam_tvc}}
}
\examples{

## Example: simulate family data with TVC based on CO model.

set.seed(4321)
fam <- simfam_tvc(N.fam = 10, design = "pop", variation = "frailty", 
       base.dist = "Weibull", frailty.dist = "gamma", depend = 1, 
       add.tvc = TRUE, tvc.type = "CO", tvc.range = c(30,60), 
       tvc.parms = c(1, 0.1, 0), allelefreq = 0.02, 
       base.parms = c(0.01, 3), vbeta = c(-1.13, 2.35))

\dontrun{
> head(fam) 
  famID indID gender motherID fatherID proband generation majorgene secondgene ageonset
1     1     1      1        0        0       0          1         2          0 61.80566
2     1     2      0        0        0       0          1         3          0 61.56996
3     1     3      0        2        1       1          2         2          0 39.42050
4     1     4      1        0        0       0          0         3          0 90.17320
5     1    13      0        3        4       0          3         3          0 51.49538
6     1    14      0        3        4       0          3         3          0 75.97238
  currentage     time status mgene  tvc.age tvc.status relation fsize naff weight
1   68.26812 61.80566      1     1 59.16387          1        4    29    3      1
2   68.60174 61.56996      1     0 39.45786          1        4    29    3      1
3   47.05410 39.42050      1     1 35.01941          1        1    29    3      1
4   44.86501 44.86501      0     0 58.67013          0        6    29    3      1
5   22.73075 22.73075      0     0 30.19254          0        3    29    3      1
6   22.71399 22.71399      0     0 40.66258          0        3    29    3      1

> summary(fam) 
Study design:           pop: population-based study with affected probands 
Baseline distribution:  Weibull 
Frailty distribution:   gamma 
Number of families:                     10 
Average number of affected per family:  3.1 
Average number of carriers per family:  3.4 
Average family size:                    16.3 
Average age of onset for affected:      48.19 
Average number of TVC event per family: 4 
Sampling weights used:                  1
}


}
\keyword{Simulation}