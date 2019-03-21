familyDesign <-
function(n=1000, affectnum=0,  m.carrier= 0, dominant.m = TRUE, dominant.s = TRUE, 
         interaction = TRUE, depend=1, 
		      base.dist="Weibull", frailty.dist="gamma", 
          vbeta= c(-1.126, 2.55, 1.6), parms=c(0.016, 3), variation="none", 
          allelefreq=c(0.02, 0.2), mrate=0.1, probandage=c(45,1.5), 
          agemin=20, agemax=100) 
{
data<-numeric()
cumind<-0
i<- 1
j<- 0
while (i <= n) {
  j <- j + 1
  dat<-try(familyStructure(i,cumind=cumind, m.carrier=m.carrier, 
  		base.dist=base.dist, frailty.dist=frailty.dist, interaction=interaction,
	    depend=depend, parms=parms, vbeta=vbeta, dominant.m=dominant.m, dominant.s=dominant.s,
	    variation=variation, allelefreq=allelefreq, mrate=mrate, probandage=probandage,
	    agemin=agemin, agemax=agemax))
	    
  if(is.null(attr(dat, "class"))){
	   # At least one parent in first gen and two sibs in the second gen should be affected
     #[,7]=generation, [,6]=proband, [,13]=status
	
	if(affectnum==3) until <- ifelse(sum(dat[dat[,7]==1,13]) > 0 & sum(dat[dat[,7]==2,13]) > 1, TRUE, FALSE)
	else if(affectnum==2) until <- ifelse(sum(dat[, 13]) >= affectnum, TRUE, FALSE)
#	else if(affectnum==1) until <- ifelse(dat[dat[,6]==1,13] == 1, TRUE, FALSE)
	else if(variation=="frailty") until <- ifelse( dat[dat[,6]==1,13] == 1, TRUE, FALSE) 
	else until <- TRUE
	
			if(!is.null(dim(dat))){
	      if(nrow(dat)>0 ){
	    	  if(until){
		   	    data<-rbind(data, dat)
        		cumind<-cumind+nrow(dat)
        		i<-i+1
	   		  }    
		    }
	    }
} # close "is.null(attr(dat, "class"))"
} # close while

data
}
