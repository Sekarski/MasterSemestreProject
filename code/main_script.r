source("betamix.r") #loads the mix of beta distribution utilities
source("tilted.r") #loads the tilting of density utilites
source("simplex2.r") #loads the distribution on s2 where the first component in betamix
source("algorithm.r") #loads the ML esimating algorithm

####################################
# Auxiliary functions
####################################
## param transformation ##
parameterize <- function(params){
  if (length(params)==2){
    p <- c(1)
    a <- c(exp(params[1]))
    b <- c(exp(params[2]))
    return(c(p,a,b))
  }
  n_betas <- (length(params)+1)/3
  e <- params[1:(n_betas-1)]
  a <- exp(params[n_betas:(2*n_betas-1)])
  b <- exp(params[-(1:(2*n_betas-1))])
  p <- exp(e)
  temp <- sum(p)
  p <- c(1,p)
  p <- p/(1+temp)
  #  return(list(a=a,b=b,p=p))
  return(c(p,a,b))
}

#### wrapper function for nll of simplex2, using a different parametrization #######
ll_simplex2 <- function(params=c(1,1,1,1,1),data){
  new_params <- parameterize(params)
  len <- length(new_params)/3
  p <- new_params[1:len]
  a <- new_params[(len+1):(2*len)]
  b <- new_params[(2*len+1):(3*len)]
  nllsimplex2(data,p,a,b)
}

#### wrapper function for nll of tilted simplex2, using a different parametrization #########
ll_tilted <- function(params=c(1,1,1,1,1),data){
  new_params <- parameterize(params)
  len <- length(new_params)/3
  p <- new_params[1:len]
  a <- new_params[(len+1):(2*len)]
  b <- new_params[(2*len+1):(3*len)]
  m <- msimplex2(p,a,b)
  res <- nlltilted(data,m,dsimplex2,p,a,b)
  # if(!is.finite(res)){
  #   print(params)
  #   print(new_params)
  #   stop(!is.finite(res))
  # }
  
  return(res)
}


#### integrated square error for estimation of tilted mix of beta distribution on the 2-simplex
iserr <- function(p_t,a_t,b_t,p_e,a_e,b_e){
  dtrue <- function(x){
    m <- msimplex2(p_t,a_t,b_t)
    dtilted(x,m,dsimplex2,p_t,a_t,b_t)
  }
  destimated <- function(x){
    m <- msimplex2(p_e,a_e,b_e)
    dtilted(x,m,dsimplex2,p_e,a_e,b_e)
  }
  integrand <- function(x1){
    x2 <- 1-x1
    x <- cbind(x1,x2)
    (destimated(x)-dtrue(x))^2
  }
  return(integrate(integrand,0,1))
}

# Make a string for the mixture, useful for plot titles
make_beta_string <- function(p,a,b){
  if (length(p) == 1){
    return(paste("Beta(",a,",",b,")",sep=""))
  }
  else{
    beta_string <- c()
    for (i in 1:length(p)){
      beta_string <- c(beta_string,paste(p[i],"Beta(",a[i],",",b[i],")",sep=""))
    }
    s <- paste(beta_string,sep="",collapse="+")
    return(s)
  }
}

# makes a D dimentional lattice on the D-simplex
make_latice <- function(D,eps,spacing){
  if (D == 2){
    x1 <- seq(from=eps, to=1-eps, by=spacing)
    x2 <- 1 - x1
    return(cbind(x1,x2))
  }
  else{
    print("high dimensions not yet implemented")
  }
}

#constructs the relative path from the code folder to where the images will be saved
image_path_string <- function(p,a,b,tilted,K,plot_type,n,R){
  string <- "../img/"
  for (i in 1:length(p)){
    string <- paste(string,"p",gsub("\\.","",p[i]),"_a",gsub("\\.","",a[i]),"_b",gsub("\\.","",b[i]),"_",sep="")
  }
  string <- sub("_$","/",string)
  if (tilted){
    string <- paste(string,"tilted/",sep="")
  }
  else{
    string <- paste(string,"untilted/",sep="")
  }
  string <- paste(string,"K",K,"/",plot_type,sep="")
  return(string)
}

##########################################
# generates first set of plots for report
##########################################

#beta mix parameters
alpha <- cbind(c(2),c(0.9),c(0.9),c(2)) #2 1 #0.7 3
beta <- cbind(c(5),c(0.5),c(5),c(2)) #5 2 #2 4
pi <- cbind(c(1),c(1),c(1),c(1)) #0.2 0.8 #0.3 0.7


x <- make_latice(2,10^(-6),0.01)


par(mfrow=c(2,2))

for (i in 1:length(pi[1,])){
  p <- pi[,i]; a <- alpha[,i]; b <- beta[,i]
  
  title1 <- paste("Original and tilted density of ",make_beta_string(p,a,b),sep="")
  
  m <- msimplex2(p,a,b)
  #y coordinates of the original pdf
  y1_star <- dsimplex2(x,p,a,b)
  #y coordinates of the tilted pdf
  y1 <- dtilted(x,m,dsimplex2,p,a,b)
  
  plot(x[,1],y1,col="blue",type="l",lwd=3, xlim=c(0,1),ylim=c(0,2),ylab="",main=title1)
  points(x[,1],y1_star,col="red",type="l",lwd=3,lty=2)
  legend("topright", legend=c("v","v*"), col=c("blue","red"), lty=1:2, lwd=3)
  m1_e <- mean(x[,1]*y1)
  m1_star_e <- mean(x[,1]*y1_star)
  abline(v=m1_e,col="blue")
  abline(v=m1_star_e,col="red",lty=2)
}



############################
# generates 2nd set of plots
############################

#Helper function to plot the histograms
make_histo <- function(R,p,a,b){
  title2 <- paste("Tilt of ",make_beta_string(p,a,b))
  x <- make_latice(2,10^(-6),0.01)
  w_star <- rsimplex2(R,p,a,b)
  m <- msimplex2(p,a,b)
  w <- rtilted(w_star,m)
  hist(w[,1],breaks=50,prob=T,main=title2,xlab="x1")
  y1 <- dtilted(x,m,dsimplex2,p,a,b)
  lines(x[,1],y1,col="blue")
}


R <- 10^4
par(mfrow=c(2,2))

pi <- c(1); alpha <- c(2); beta <- c(5)
make_histo(R,pi,alpha,beta)
pi <- c(1); alpha <- c(0.9); beta <- c(0.5)
make_histo(R,pi,alpha,beta)
pi <- c(0.2,0.8); alpha <- c(2,1); beta <- c(5,2)
make_histo(R,pi,alpha,beta)
pi <- c(0.2,0.3,0.5); alpha <- c(2,1,0.5); beta <- c(5,2,3)
make_histo(R,pi,alpha,beta)


##############################################
# Parameter estimation
##############################################

# function that plots all the density estimates on the same graph, along with the original density
# _t is appended to anything meaning "true" and _e meaning "estimate" from her on out
plot_dens_estimates <- function(p_t,a_t,b_t,estimates,main,...){
  x <- make_latice(2,10^(-6),0.01)
  m_t <- msimplex2(p_t,a_t,b_t)
  y1_t <- dtilted(x,m_t,dsimplex2,p_t,a_t,b_t)
  plot(x[,1],y1_t,col="blue",lwd=3, type="l",xlim=c(0,1),ylim=c(0,2),ylab="", main=main,...)
  legend("topright", legend=c("True","Estimates"), col=c("blue","red"), lty=1, lwd=c(3,1))
  
  for (i in 1:length(estimates[,1])){
    nbetas <- length(estimates[1,])/3
    p_e <- c(estimates[i,1:nbetas])
    a_e <- c(estimates[i,(nbetas+1):(2*nbetas)])
    b_e <- c(estimates[i,(2*nbetas+1):(3*nbetas)])
    m_e <- msimplex2(p_e,a_e,b_e)
    y1_e <- dtilted(x,m_e,dsimplex2,p_e,a_e,b_e)
    lines(x[,1],y1_e,col="red")
  }
  lines(x[,1],y1_t,col="blue",lwd=3)
}

#Model params
alpha <- c(2,5)
beta <- c(3,2)
pi <- c(0.5,0.5)
#Number of observations
n <- 100 #50 100 200 500
#number of replicas
R <- 100 #200
#number of different initial values per replica
NI <- 10
#tilte for figures
smodel <- paste("Tilt of ",make_beta_string(pi,alpha,beta),sep="")
ftitle <- paste(smodel,", observations: ",n,", replicas: ",R)
#number of betas in the mix
K=2

mle_parized <- c()
n_unconverged <- 0 #will count the estimates that terminate before converging
nmd <- 0 #will count the number of Nelder-Mead degenerate cases
for (i in 1:R){
  w_star <- rsimplex2(n,pi,alpha,beta)
  m <- msimplex2(pi,alpha,beta)
  w <- rtilted(w_star,m)

  start <- runif((K*3-1)*NI,-0.5,0.5) #error: function cannot be evaluated at initial parameters
  dim(start) <- c((K*3-1),NI)
  
  #if I try using L-BFGS-B and bounds, I get problems...
  res <- MLestimation(inivs=start,nll=ll_tilted,data=w,method="Nelder-Mead",maxit=500) #lower=c(-Inf,-20,-20,-20,-20),upper=c(Inf,20,20,20,20)
  # res <- MLestimation(inivs=start,nll=ll_tilted,data=w,method="L-BFGS-B",maxit=500,lower=c(-10,-5,-5,-5,-5),upper=c(10,5,5,5,5))
  
  
  if (res$convergence == 1){n_unconverged <- n_unconverged + 1} #count the number of fits that reached maxit
  if (res$convergence == 10){nmd <- nmd + 1} #counts the number of fits where NM simplex id degenerate
  
  mle_parized <- rbind(mle_parized,res$par) #retrieve parameters in original scale
}
mle <- t(apply(mle_parized,1,parameterize)) #had to filp the result for conveinence


##################################
# Estimates Plots
##################################
par(mfrow=c(1,1))

#plot to RStudio output
plot_dens_estimates(pi,alpha,beta,mle,ftitle)
#plot to file
file <- image_path_string(pi,alpha,beta,T,K,"densities",n,R) #directory path. Works in Windos, have not tried it in MacOS or Linux
dir.create(file.path(".",file),recursive = T) #if the directory structure doesn't existe yet, creats it
pdf(paste(file,"/n",n,"_R",R,".pdf",sep="")) #opens an output flow to a pdf file
plot_dens_estimates(pi,alpha,beta,mle,ftitle) #write the plot to the pdf
dev.off() #closes and saves the pdf

#helper function for boxplot variable names
make_names <- function(K){
  p <- c()
  a <- c()
  b <- c()
  for (k in 1:K){
    p <- c(p,paste("pi",k,sep=""))
    a <- c(a,paste("alpha",k,sep=""))
    b <- c(b,paste("beta",k,sep=""))
  }
  return(c(p,a,b))
}

#plot to RStudio output
boxplot(log(mle), names=make_names(K),main=ftitle)
#plot to file
file <- image_path_string(pi,alpha,beta,T,K,"bxplots",n,R)
dir.create(file.path(".",file),recursive = T)
pdf(paste(file,"/n",n,"_R",R,".pdf",sep=""))
boxplot(log(mle), names=make_names(K),main=ftitle)
dev.off()


#plot to RStudio output
pairs(log(mle),main=ftitle)
#plot to file
file <- image_path_string(pi,alpha,beta,T,K,"pairs",n,R)
dir.create(file.path(".",file),recursive = T)
pdf(paste(file,"/n",n,"_R",R,".pdf",sep=""))
pairs(log(mle),main=ftitle)
dev.off()




