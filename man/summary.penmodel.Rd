\name{summary.penmodel}
\alias{summary.penmodel}

\title{
Summary method for class \code{penmodel}
}
\description{
Provides a summary of a fitted penetrance model.
}
\usage{
\method{summary}{penmodel}(object, correlation=FALSE, \ldots)
}

\arguments{
  \item{object}{
An object class of \code{'penmodel'}, a fitted model by \code{\link{penmodel}} or \code{\link{penmodelEM}} functions.
}
 \item{correlation}{Logical; if TRUE, returns the correlation matrix of the estimated parameters.}
 \item{...}{Further arguments passed to or from other methods.}
}

\value{
Returns the object of class \code{'summary.penmodel'}, including the following summary values:
\item{coefficients}{4 x 4 matrix with columns for parameter estimates of transformed baseline parameters (\eqn{\lambda, \rho}) and regression coefficients for gender and mutation status (\eqn{\beta_s, \beta_g}), their standard errors (or robust standard errors if \code{robust=TRUE} was selected when fitting the penetrance model), \eqn{t}-statistics and corresponding two-sided \eqn{p}-values.}
\item{varcov}{4 x 4 variance covariance matrix of the parameter estimates.}
\item{correlation}{ Correlation matrix corresponding to the specified variance covariance matrix, if \code{correlation=TRUE} is specified.}
\item{pen70}{ Penetrance estimates by age 70 specific to gender and mutation-status subgroups including their standard errors and 95\% confidence intervals.}
}

\examples{
# Simulated family data 

fam <- simfam(N.fam=300, design="pop+", variation="none", base.dist="Weibull", 
       base.parms=c(0.01,3), vbeta=c(-1.13, 2.35), agemin=20, allelefreq=0.02)
 
# Penetrance model fit for the simulated family data

fit <- penmodel(parms=c(0.01, 3), vbeta=c(-1.13, 2.35), data=fam, 
       design="pop+", base.dist="Weibull")

# Summary of the model parameter and penetrance estimates from model fit

summary(fit)

\dontrun{
Coefficients: 
            Estimate Std. Error t value Pr(>|t|)  
log(lambda)   -4.637    0.07989 -58.048  0.01097 *
log(rho)       1.094    0.03910  27.980  0.02274 *
beta.sex      -1.315    0.16298  -8.066  0.07852 .
beta.gene      2.545    0.21423  11.881  0.05346 .
Signif. codes:   0 '***'  0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

Penetrance (\%) by age 70: 
         Male Carrier Female Carrier Male Noncarrier Female Noncarrier
Estimate        32.47          76.82            3.03             10.83
SE               3.85           3.22            0.74              2.17

95\% Confidence intervals on the penetrances: 
           Male Carrier Female Carrier Male Noncarrier Female Noncarrier
lowerlimit        25.59          70.43            1.94              7.43
upperlimit        40.70          83.08            4.85             15.87
}

}
\author{
Yun-Hee Choi	
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{penmodel}, \link{penmodelEM}, \link{print.penmodel}, \link{print.summary.penmodel} \link{plot.penmodel}}

Function \code{\link{coef}} will extract the matrix of coefficients with standard errors, \eqn{t}-statistics and \eqn{p}-values.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ summary }
