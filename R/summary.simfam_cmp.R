summary.simfam_cmp <- function(object, digits = max(3, getOption("digits") - 3), ...){

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
"twostage: two-stage design")[match(design, c("pop","pop+","cli","cli+","twostage"))]

 cat("Study design:                          ", design,"\n")
}
  
if(length(base.dist)==2){ #competing events
  cat("Baseline distribution for event 1:     ", base.dist[1],"\n")
  cat("Baseline distribution for event 2:     ", base.dist[2],"\n")
  if(!is.null(variation)) {
	if(variation == "frailty") 
  cat("Frailty distribution:                  ", attr(object,"frailty.dist"),"\n")
	else if(variation=="secondgene") 
  cat("Residucal familial correlation indueced by a second gene variation \n")
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
avgnumdf1 <- sum(data$status==1, na.rm=TRUE)/numfam
avgnumdf2 <- sum(data$status==2, na.rm=TRUE)/numfam
avgnumcarr <- sum(data$mgene,na.rm=TRUE)/numfam
if(!is.null(data$fsize)) avgfamsize <- mean(data$fsize[data$proband==1])
else avgfamsize <- mean(aggregate(data$indID, by=list(data$famID), length)[,2])
avgageonset1 <- mean(data$time[data$status==1], na.rm=TRUE)
avgageonset2 <- mean(data$time[data$status==2], na.rm=TRUE)
wt <- unique(data$weight)

ans <- list(num.fam=numfam, avg.num.affected1=avgnumdf1, avg.num.affected2=avgnumdf2,
            avg.num.carriers=avgnumcarr, avg.family.size=avgfamsize, 
            avg.ageonset1=avgageonset1, avg.ageonset2=avgageonset2)

if(i==1) ans.high <- ans
nn<-c(
"Number of families:                    ", 
"Average number of event 1 per family:  ", 
"Average number of event 2 per family:  ", 
"Average number of carriers per family: ", 
"Average family size:                   ",
"Average age of onset for event 1:      ",
"Average age of onset for event 2:      ")

for(j in 1:7)  cat(nn[j], ans[[j]], "\n")

#cat("Sampling weights used:                 ", unique(data$weight))
#cat("\n")

if(i==2) ans<-cbind(HighRisk=ans.high, LowRisk=ans)

	}
}
  else{ # single event (competing = FALSE)
    
    cat("Baseline distribution:                 ", base.dist,"\n")
    if(!is.null(variation)) {
      if(variation == "frailty") 
        cat("Frailty distribution:                  ", attr(object,"frailty.dist"),"\n")
      else if(variation=="secondgene") 
        cat("Residucal familial correlation indueced by a second gene variation \n")
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
      avgnumaffec <- sum(data$status==1, na.rm=TRUE)/numfam
      avgnumcarr <- sum(data$mgene,na.rm=TRUE)/numfam
      if(!is.null(data$fsize)) avgfamsize <- mean(data$fsize[data$proband==1])
      else avgfamsize <- mean(aggregate(data$indID, by=list(data$famID), length)[,2])
      avgageonset <- mean(data$time[data$status==1], na.rm=TRUE)
      wt <- unique(data$weight)
      
      ans <- list(num.fam=numfam, avg.num.affected=avgnumaffec, 
                  avg.num.carriers=avgnumcarr, avg.family.size=avgfamsize, 
                  avg.ageonset=avgageonset)
      
      if(i==1) ans.high <- ans
      nn<-c(
        "Number of families:                    ", 
        "Average number of event 1 per family:  ", 
        "Average number of event 2 per family:  ", 
        "Average number of carriers per family: ", 
        "Average family size:                   ",
        "Average age of onset for event 1:      ",
        "Average age of onset for event 2:      ")
      
      for(j in 1:5)  cat(nn[j], ans[[j]], "\n")
      
      cat("Sampling weights used:                 ", unique(data$weight))
      cat("\n")
      
      if(i==2) ans<-cbind(HighRisk=ans.high, LowRisk=ans)
      
    }
    
  } # ends single event

attr(ans,"design") <- design
attr(ans,"base.dist") <- base.dist
if(is.null(variation)) attr(ans, "frailty.dist") <- NULL
else if(variation=="frailty") attr(ans,"frailty.dist") <- frailty.dist
else  attr(ans,"frailty.dist") <- NA
attr(ans, "variation") <- variation
class(ans) <- "summary.simfam_cmp"
invisible(ans)

}
