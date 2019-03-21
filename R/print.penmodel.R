print.penmodel <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  
  cat("Call: \n")
  cat("Penetrance model for",attr(x,"design"),"design using",
      attr(x,"base.dist"), "baseline distribution \n")
  cat("Minimum age at onset used: ", attr(x, "agemin"))
  cat("\n")
  cat("\nEstimates: \n")
  print(x$estimates)
  if(attr(x, "robust")) cat("Robust standard errors was obtained. \n")
  invisible(x)

}
