testthat::test_that("conversion works", {
  testthat::expect_equal(MiMeMo.tools::micro_to_full(1e6), 1)
})
