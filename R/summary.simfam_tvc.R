summary.simfam_tvc <- function(object, digits = max(3, getOption("digits") - 3), ...){

  savedig <- options(digits = digits)
  on.exit(options(savedig))
  design <- attr(object, "design")
  variation <- attr(object, "variation") 
  frailty.dist <- attr(object, "frailty.dist") 
  base.dist <- attr(object, "base.dist") 

if(!is.null(design)){  
design.name <- c("pop: population-based study with affected probands",
"pop+: population-based study with affected and mutation carrier probands",
"cli: clinic-based study with affected probands",
"cli+: clinic-based study with affected and mutation carrier probands",
"twostage: two-stage design",
"noasc: no ascertainment")[match(design, c("pop","pop+","cli","cli+","twostage", "noasc"))]

 cat("Study design:          ", design.name,"\n")
 }
if(!is.null(base.dist)) cat("Baseline distribution: ", base.dist,"\n")
if(!is.null(variation)) {
    if(length(variation)==1){
      if(variation == "frailty") 
        cat("Frailty distribution:  ", attr(object,"frailty.dist"),"\n")
      else if(variation=="kinship")
        cat("Frailty distribution:                   Lognormal with Kinship matrix \n")
      else if(variation=="IBD")
        cat("Frailty distribution:                   Lognormal with IBD matrix \n")
      else if(variation=="none")
        cat("Frailty distribution:                   No frailties \n")
    }
    else if(length(variation)==2)
      cat("Frailty distribution:                   Lognormal with Kinship and IBD matrices \n")
}
  
if(is.null(object$weight)) object$weight <- 1

	k <- ifelse(is.null(design), 1, ifelse(design=="twostage",2,1))
	for(i in 1:k){
		  if(is.null(design)) data <- object
  		else if(design=="twostage"){
	 		  if(i==1){ 
	 			data <- object[object$weight==1,]
				cat("\n ## High risk families", "\n")
	 		  }
	 		  else{
	 			data <- object[object$weight!=1,]
				cat("\n ## Low risk families", "\n")
	 		  }
  		}
		  else data <- object	 

numfam <- length(data$famID[data$proband==1])
avgnumaffec <- sum(data$status, na.rm=TRUE)/numfam
avgnumcarr <- sum(data$mgene,na.rm=TRUE)/numfam
if(!is.null(data$fsize)) avgfamsize <- mean(data$fsize[data$proband==1])
else avgfamsize <- mean(aggregate(data$indID, by=list(data$famID), length)[,2])

avgageonset <- mean(data$time[data$status==1], na.rm=TRUE)
wt <- unique(data$weight)

if(!is.null(data$tvc.status)) {
  avgnumtvc <- sum(data$tvc.status, na.rm=TRUE)/numfam
  ans <- list(num.fam=numfam, avg.num.affected=avgnumaffec, avg.num.carriers=avgnumcarr, 
              avg.family.size=avgfamsize, avgageonset=avgageonset,
              avg.num.tvc = avgnumtvc)
}
else{
ans <- list(num.fam=numfam, avg.num.affected=avgnumaffec, avg.num.carriers=avgnumcarr, 
            avg.family.size=avgfamsize, avgageonset=avgageonset)
}

if(i==1) ans.high <- ans

if(is.null(data$tvc.status)){
  nn<-c(
    "Number of families:                    ", 
    "Average number of affected per family: ", 
    "Average number of carriers per family: ", 
    "Average family size:                   ",
    "Average age of onset for affected:     ")
  for(j in 1:5)  cat(nn[j], ans[[j]], "\n")
}
else{
  nn<-c(
    "Number of families:                     ", 
    "Average number of affected per family:  ", 
    "Average number of carriers per family:  ", 
    "Average family size:                    ",
    "Average age of onset for affected:      ",
    "Average number of TVC events per family:")
  for(j in 1:6)  cat(nn[j], ans[[j]], "\n")
}

#cat("Sampling weights used:                 ", unique(data$weight))
#cat("\n")

if(i==2) ans<-cbind(HighRisk=ans.high, LowRisk=ans)

}

attr(ans,"design") <- design
attr(ans,"base.dist") <- base.dist
if(is.null(variation)) attr(ans, "frailty.dist") <- NULL
else if(length(variation)==1 & variation == "none") attr(ans,"frailty.dist") <- NULL
else  attr(ans,"frailty.dist") <- frailty.dist
attr(ans, "variation") <- variation
class(ans) <- "summary.simfam_tvc"
invisible(ans)

}
