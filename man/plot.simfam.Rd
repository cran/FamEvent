\name{plot.simfam}
\alias{plot.simfam}
\title{Plot pedigrees}
\description{
Plot method for \code{simfam}. Provides pedigree plots for specified families generated from \code{simfam} function with option to save plots into a pdf file.
}
\usage{
\method{plot}{simfam}(x, famid, pdf=FALSE, file=NULL, \ldots)
}
\arguments{
 \item{x}{An object created by the function \code{simfam}.}
 \item{famid}{a list of family IDs to plot.}
 \item{pdf}{logical; if TRUE, pedigree plots are saved in a pdf file. If FALSE, plot pedigrees on current plotting device. Default=FALSE.}
 \item{file}{a pdf file name to save the pedigree plots; Default=NULL, if NULL, the pedigrees are saved with file name "pedigreeplot.pdf".}
 \item{...}{additional arguments passed on to the plot function.}
 
}
\value{
Pedigrees of specified families created by the function \code{simfam} are plotted along with the affection and carrier mutation statuses of family members. Probands from each pedigree are indicated using red color.

When it is used with \code{\link{carrierprobgeno}} or \code{\link{carrierprobpheno}} functions, the plot function displays the carrier probabilities for those with missing carrier status.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{simfam}, \link{summary.simfam}}	
}

\examples{

# Simulated family data

fam <- simfam(N.fam=300, design="pop+", variation="secondgene", 
       base.dist="Weibull", allelefreq=c(0.02, 0.2), 
       base.parms=c(0.01,3), vbeta=c(-1.13, 2.35, 0.5), agemin=20)
 
# Pedigree plots for first five simulated families

plot.simfam(fam, famid=c(1:5), pdf=TRUE, file="pedigrees.pdf")
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Pedigree plot }
