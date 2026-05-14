# joint probability of genotype and phenotype
cprob <- function(theta, X, y, delta, p, base.dist, cuts, nbase)
{
  
  nX <- dim(X)[2]
  vbeta <- theta[(nbase+1):(nbase+nX)]  
  if(base.dist=="lognormal") etheta <-  c(theta[1], exp(theta[2]))
  else etheta <- exp(theta[1:nbase])
  
  xbeta <- c(X%*%vbeta)
  
  haz <- hazards(base.dist, y, etheta, cuts)*exp(xbeta)
  Haz <- cumhaz(base.dist, y, etheta, cuts)*exp(xbeta)
  out <- ifelse(delta==1, haz*exp(-Haz)*p, exp(-Haz)*p)
  #return((haz^delta)*exp(-Haz)*p)
  return(out)
}
