# density function of a mix of beta distributions ###############
dbetamix <- function(x,p,a,b){
  params <- cbind(p,a,b)
  blop <- function(p,x){
    p[1]*dbeta(x,p[2],p[3])
  }
  temp <- apply(params,1,blop,x)
  y <- apply(temp,1,sum)
  return(y)
}



# sampling of a mix of beta distributions V1 ################
rbetamix <- function(r,p,a,b){
  n <- length(p)
  y <- c()
  idx <- sample(1:n,r,replace=TRUE,prob=p)
  for (i in 1:r){
    y <- c(y,rbeta(1,a[idx[i]],b[idx[i]]))
  }
  return(y)
}

rbetamix_old <- function(r,p,a,b){
  y <- c()
  for (i in 1:length(p)) {
    n <- r*p[i]
    temp <- rbeta(n,a[i],b[i])
    y <- c(y,temp)
  }
  if (length(y) < r){
    n <- r-length(y)
    temp <- rbeta(n,a[1],b[1])
    y <- c(y,temp)
  }
  return(y)
}

# mean of a mix of beta distributions ################
mbetamix <- function(p,a,b){
  sum(p*a/(a+b))
}

# negative log likelihood of mix of bedat distributions ############
# included only for completness
nllbetamix <- function(data,p,a,b){
  y <- dbetamix(data,p,a,b)
  y <- log(y)
  y <- -sum(y)
  return(y)
}