source("betamix.r")

# density on the 2 simplex where the first coord is a mix of beta
dsimplex2 <- function(x,pi,alpha,beta){
  y <- dbetamix(x[,1],pi,alpha,beta)
  return(y)
}

# sampling on the 2 simplex where the first coord is a mix of beta
rsimplex2 <- function(R,pi,alpha,beta){
  x1 <- rbetamix(R,pi,alpha,beta)
  x2 <- 1-x1
  x <- cbind(x1,x2)
  return(x)
}

# mean of the distribution on the 2 simplex where the first coord is a mix of beta
msimplex2 <- function(pi,alpha,beta){
  m1 <- mbetamix(pi,alpha,beta)
  m2 <- 1-m1
  m <- c(m1,m2)
  return(m)
}

#negative log likelihood of the above distribution
nllsimplex2 <- function(data,pi,alpha,beta){
  y <- dsimplex2(data,pi,alpha,beta)
  y <- log(y)
  y <- -sum(y)
  return(y)
}