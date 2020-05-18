
#does the MLE estimation. inivs are a 3K-1 x NI initial values for the algorithm. nll is the
# negative log likelihood function, parameterized acording to desire.
MLestimation <- function(inivs,nll,data,method="Nelder-Mead",maxit=500,lower=-Inf,upper=Inf){
  maxll <- -Inf
  res <- NULL
  if (method %in% c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent")){
    NI <- dim(inivs)[2]
    for (i in 1:NI){
      s <- inivs[,i]
      temp <- optim(par=s,fn=nll,data=data,method=method,lower=lower,upper=upper,control=list(maxit=maxit))
      if (temp$value > maxll){
        res <- temp
        maxll <- res$value
      }
    }
  } else if (method == "cutom"){
    #will try to use a first round of BFGS with a single component, then feed the results with some extra random
    #parameters into Nelder-Mead with multiple components to see if it's better
  }
  return(res)
}

