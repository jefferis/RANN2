context('WANN')

test_that("constructing WANN and getting points", {
  data(kcpoints)
  w1=WANN(kcpoints[[1]])
  expect_equivalent(w1$getPoints(), kcpoints[[1]])
  w1_notree=WANN(kcpoints[[1]], FALSE)
})

test_that("Basic queries", {

  k <- 1
  e <- 0
  r <- NA

  p1=kcpoints[[1]]
  w1=WANN(p1)
  w1sq=w1$querySelf(k = k, eps = e, radius = r)
  expect_equivalent(w1sq$nn.dists, matrix(0, nrow(kcpoints[[1]])))
  expect_equal(w1$query(p1, k = k, eps = e, radius = r), w1sq)
  p2=kcpoints[[2]]
  expect_equal(w1$query(p2, k = k, eps = e, radius = r), nn2(p1, p2, k = k))
})

test_that("Basic fixed radius queries", {

  r <- 0.5
  k <- 3
  e <- 0

  p1=kcpoints[[1]]
  w1=WANN(p1)
  w1sq=w1$querySelf(k = k, eps = e, radius = r)
  expect_equal(w1sq$nn.dists, nn2(p1, k = k, eps = e, radius = r)$nn.dists)

})

test_that("Queries using WANN objects", {

  k <- 3
  e <- 0
  r <- NA

  p1=kcpoints[[1]]
  w1=WANN(p1)
  p2=kcpoints[[2]]
  w2=WANN(p2)
  expect_equal(w1$queryWANN(w2$.CppObject, k = k, eps = e, radius = r), nn2(p1, p2, k = k, radius = r))
})

test_that("Build and delete trees explicitly", {

  k <- 3
  e <- 0
  r <- NA

  p1=kcpoints[[1]]
  w1=WANN(p1)
  p2=kcpoints[[2]]
  w2=WANN(p2)
  expect_is(w1$querySelf(k = k, eps = e, radius = r), 'list')
  expect_is(w1$query(p2, k = k, eps = e, radius = r), 'list')
  expect_is(w1$queryWANN(w2$.CppObject, k = k, eps = e, radius = r), 'list')
})
