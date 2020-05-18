source("../betamix.r", chdir = TRUE)
library(testthat)


test_that("mix of 1 is same as beta distribution",{
  x <- seq(from=0, to=1, by=0.01)
  a <- 2
  b <- 5
  expect_equal(dbetamix(x,1,c(a),c(b)),dbeta(x,a,b),tolerance=10e-4)
  a <- 0.6
  b <- 4
  expect_equal(dbetamix(x,1,c(a),c(b)),dbeta(x,a,b),tolerance=10e-4)
  a <- 3
  b <- 3
  expect_equal(dbetamix(x,1,c(a),c(b)),dbeta(x,a,b),tolerance=10e-4)
  a <- 0.7
  b <- 0.2
  expect_equal(dbetamix(x,1,c(a),c(b)),dbeta(x,a,b),tolerance=10e-4)
})

test_that("dbetamix integrats to 1",{
  p <- c(1)
  a <- c(2)
  b <- c(5)
  expect_equal(integrate(dbetamix,0,1,p,a,b)$value,1,tolerance=0.01)
  p <- c(1)
  a <- c(0.5)
  b <- c(3)
  expect_equal(integrate(dbetamix,0,1,p,a,b)$value,1,tolerance=0.01)
  p <- c(1)
  a <- c(0.3)
  b <- c(0.7)
  expect_equal(integrate(dbetamix,0,1,p,a,b)$value,1,tolerance=0.01)
  p <- c(0.3,0.3,0.4)
  a <- c(3,10,0.4)
  b <- c(2,0.9,0.9)
  expect_equal(integrate(dbetamix,0,1,p,a,b)$value,1,tolerance=0.01)
})

test_that("mbetamix is acurate",{
  p <- c(0.5,0.5)
  a <- c(2,0.8)
  b <- c(5,3)
  expect_equal(mbetamix(p,a,b),0.248, tolerance=0.001)
})

test_that("rbetamix is distributed correctly",{
  x <- seq(from=0, to=1, by=0.01)
  R <- 10^5
  p <- c(0.5,0.5)
  a <- c(2,0.8)
  b <- c(5,3)
  y <- dbetamix(x,p,a,b)
  w <- rbetamix(R,p,a,b)
  hist(w,breaks=100,prob=T)
  lines(x,y)
  p <- c(0.9,0.1)
  a <- c(2,0.8)
  b <- c(5,0.2)
  y <- dbetamix(x,p,a,b)
  w <- rbetamix(R,p,a,b)
  hist(w,breaks=100,prob=T)
  lines(x,y)
  skip("numerical test not implemented yet. Check histograms")
})