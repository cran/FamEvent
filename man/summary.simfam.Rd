\name{summary.simfam}
\alias{summary.simfam}

\title{
Summary of simulated data
}
\description{
Summary method for \code{simfam}. Provides summaries of simulated data.
}
\usage{
\method{summary}{simfam}(object, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object class of \code{simfam} generated from \code{simfam} function
}
 \item{...}{
 	Additional arguments
 }

}
\value{
It displays the following summary values:
\item{design}{Study Design.}
\item{numfam}{Number of families simulated.}
\item{avgnumaffec}{Average number of affected individuals per family.}
\item{avgnumcarr}{Average number of mutation carriers per family.}
\item{avgfamsize}{Average family size.}
\item{aveageonset}{Average age of onset for affected individuals.}
}

\author{
Yun-Hee Choi	
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{simfam}}
}
\examples{
fam <- simfam(N.fam=100, design="pop", variation="frailty", depend=1,
       frailty.dist="gamma", base.dist="Weibull", base.parms=c(0.01, 3), 
       vbeta=c(-1.13, 2.35))

summary(fam)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ summary }
