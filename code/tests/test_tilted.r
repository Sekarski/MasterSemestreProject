source("../tilted.r", chdir = TRUE)
source("../simplex2.r", chdir = TRUE)
library(testthat)

test_that("tilted distribution has mean 1/D",{
  p <- c(0.5,0.5)
  a <- c(2,0.8)
  b <- c(5,3)
  m <- msimplex2(p,a,b)
  eps <- 10e-6
  x1 <- seq(from=eps, to=1-eps, by=0.01)
  x2 <- 1-x1
  x <- cbind(x1,x2)
  y <- dtilted(x,m,dsimplex2,p,a,b)
  m_e <- apply(x*y,2,mean)
  expect_equivalent(m_e,c(0.5,0.5),tolerance=0.1)
  
  
  
})