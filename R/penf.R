penf <- function(est, age, sex, mut, base.dist, agemin, frailty.dist=NULL){
  H <- cumhaz(base.dist, age-agemin, exp(est[1:2]))*exp(est[3]*sex+est[4]*mut) 
  if(is.null(frailty.dist)) return(1-exp(-H))
  else return(1-laplace(dist=frailty.dist, g=H, k=est[5]))
  
}