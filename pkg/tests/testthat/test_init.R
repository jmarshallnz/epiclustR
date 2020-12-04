context("Initialisation")
library(epiclustR)

test_that("init_priors returns a named list", {
  def_priors <- init_priors()
  expect_type(def_priors, "list")
  expect_named(def_priors)
})
