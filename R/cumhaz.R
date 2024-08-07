cumhaz  <-  function(dist="Weibull", t, parms, cuts=NULL){

if(dist=="Weibull")	 	chaz <- (parms[1]*t)^parms[2] 
else if(dist=="loglogistic")	chaz <- log(1+(parms[1]*t)^parms[2])
else if(dist=="Gompertz") chaz <- parms[1]*(exp(parms[2]*t)-1)/parms[2]
else if(dist=="lognormal") chaz <- -log(1-pnorm((log(t)-parms[1])/parms[2]))
else if(dist=="gamma") chaz <- -log(1-pgamma(t,shape=parms[2], scale=1/parms[1] ))
else if(dist=="logBurr") chaz <- parms[3]*log(1+(parms[1]*t)^parms[2]/parms[3])
else if(dist=="piecewise") chaz <- Hpch(t, cuts=cuts, levels=parms)
else stop("Unrecognized baseline distribution")
return(chaz)
}