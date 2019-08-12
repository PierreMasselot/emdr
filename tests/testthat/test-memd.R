test_that("Reconstruction is correct for non noise-assisted EMD", {
  # Classical EMD
  x <- rnorm(100, 10, 15)
  res <- memd(x)
  expect_equal(apply(res, 1, sum), x)
  
  # Multivariate EMD
  x <- matrix(rnorm(300), nrow = 100, ncol = 3)
  res <- memd(x)
  expect_equivalent(apply(res, c(1,3), sum), x)
})

test_that("keep.noise works", {
  # Univariate EMD
  x <- rnorm(100, 10, 15)
  
  res <- memd(x, l = 2)
  expect_length(dim(res), 2)
  
  res <- memd(x, l = 2, keep.noise = T)
  expect_equal(dim(res)[3], 3)
  
  res <- memd(x, l = 4, keep.noise = T)
  expect_equal(dim(res)[3], 5)
  
  # Multivariate EMD
  x <- matrix(rnorm(300), nrow = 100, ncol = 3)
  res <- memd(x, l = 2)
  expect_equal(dim(res)[3], 3)
  
  res <- memd(x, l = 2, keep.noise = T)
  expect_equal(dim(res)[3], 5)
})

test_that("memd.stop has an effect", {
  # Univariate
  x <- rnorm(100, 10, 15)
  res <- memd(x)
  nextr <- sum(find.extrema(res[,ncol(res)])$nextrema)
  expect_lt(nextr, 2)
  
  # Multivariate
  x <- matrix(rnorm(300), nrow = 100, ncol = 3)
  res <- memd(x)
  nextr <- apply(res[,dim(res)[2],], 2, function(x){
    sum(find.extrema(x)$nextrema)
  })
  expect_equal(sum(nextr >= 3), 0)
})