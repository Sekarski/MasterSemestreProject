#get tilted density on the simplex for arbitray density and D
dtilted <- function(w,m,f,...){
  mw <- sweep(x=w, MARGIN=2, STATS=m, FUN="*")
  mTw <- apply(mw,1,sum)
  w_star <- sweep(x=mw, MARGIN=1, STATS=mTw, FUN="/")
  D <- length(m)
  dens <- prod(m) * f(w_star,...)/(D*(mTw)^(D+1))
  return(dens)
}


#get w's from wstar's samples
rtilted <- function(wstar,m){
  d <- sweep(x=wstar, MARGIN=2, STATS=m, FUN="/")
  s <- apply(d, 1, sum)
  w_proposed <- sweep(x=d, MARGIN=1, STATS=s, FUN="/")
  d <- sweep(x=w_proposed, MARGIN=2, STATS=m, FUN="*")
  s <- apply(d, 1, sum)
  #aceptance-rejection step
  u <- runif(nrow(w_proposed))
  i <- ( u<= min(m)/s)
  w_accepted <- w_proposed[i,]
  return(w_accepted)
}


# negative log likelihood  of a tilted density
nlltilted <- function(data,m,f,...){
  y <- dtilted(data,m,f,...)
  y <- log(y)
  y <- -sum(y)
  return(y)
}