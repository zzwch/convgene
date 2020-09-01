library(testthat)

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("ConvertHomoGene", {
  expect_error(ConvertHomoGene(c("PECAM1", "GYPA", "GZMB", "GZMH", "GZMK"),
                               orgA = "human",
                               orgB = "mouse"), NA)
})

test_that("str_to_gene", {
  expect_identical(str_to_gene("Etv2, Pecam1", tolower = T, capitalize = T), c("Etv2", "Pecam1"))
})
