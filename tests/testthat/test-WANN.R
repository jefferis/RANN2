context('WANN')

test_that("constructing WANN and getting points", {
  data(kcpoints)
  w1=WANN(kcpoints[[1]])
  expect_equivalent(w1$getPoints(), kcpoints[[1]])
  w1_notree=WANN(kcpoints[[1]], FALSE)
})

test_that("Basic queries", {
  p1=kcpoints[[1]]
  w1=WANN(p1)
  w1sq=w1$querySelf(1,0)
  expect_equivalent(w1sq$nn.dists, matrix(0, nrow(kcpoints[[1]])))
  expect_equal(w1$query(p1, 1, 0), w1sq)
  p2=kcpoints[[2]]
  expect_equal(w1$query(p2, 1, 0), nn2(p1, p2, k=1))
})

test_that("Queries using WANN objects", {
  p1=kcpoints[[1]]
  w1=WANN(p1)
  p2=kcpoints[[2]]
  w2=WANN(p2)
  expect_equal(w1$queryWANN(w2$.CppObject, 1, 0), nn2(p1, p2, k=1))
})

test_that("Build and delete trees explicitly", {
  p1=kcpoints[[1]]
  w1=WANN(p1)
  p2=kcpoints[[2]]
  w2=WANN(p2)
  expect_is(w1$querySelf(1, 0), 'list')
  expect_is(w1$query(p2, 1, 0), 'list')
  expect_is(w1$queryWANN(w2$.CppObject, 1, 0), 'list')
})
