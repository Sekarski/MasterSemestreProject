source("../simplex2.r", chdir = TRUE)
library(testthat)

test_that("1st component is betamix",{
  x1 <- seq(from=0, to=1, by=0.01)
  x2 <- 1-x1
  x <- cbind(x1,x2)
  a <- 2
  b <- 5
  expect_equal(dsimplex2(x,1,c(a),c(b)),dbetamix(x[,1],1,a,b),tolerance=10e-4)
  a <- 0.6
  b <- 4
  expect_equal(dsimplex2(x,1,c(a),c(b)),dbetamix(x[,1],1,a,b),tolerance=10e-4)
  a <- c(3,0.7)
  b <- c(3,0.2)
  p <- c(0.7,0.3)
  expect_equal(dsimplex2(x,p,a,b),dbetamix(x[,1],p,a,b),tolerance=10e-4)
})

test_that("msimplex2 is acurate",{
  p <- c(0.5,0.5)
  a <- c(2,0.8)
  b <- c(5,3)
  expect_equal(msimplex2(p,a,b),c(0.248,0.752), tolerance=0.001)
})

test_that("rsimplex2 is distributed correctly",{
  x1 <- seq(from=0, to=1, by=0.01)
  x2 <- 1-x1
  x <- cbind(x1,x2)
  R <- 10^5
  p <- c(0.5,0.5)
  a <- c(2,0.8)
  b <- c(5,3)
  y <- dsimplex2(x,p,a,b)
  w <- rsimplex2(R,p,a,b)
  hist(w[,1],breaks=100,prob=T)
  lines(x[,1],y)
  p <- c(0.9,0.1)
  a <- c(2,0.8)
  b <- c(5,0.2)
  y <- dsimplex2(x,p,a,b)
  w <- rsimplex2(R,p,a,b)
  hist(w[,1],breaks=100,prob=T)
  lines(x[,1],y)
  skip("numerical test not implemented yet. Check histograms")
})