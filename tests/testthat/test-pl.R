context("test error")

#print(load("tests/testthat/testdata.Rdata"))
print(load("testdata.Rdata"))
test_that("pl_dotplot", {
  expect_error(pl_dotplot(data, features, group_by, groups = c("brEC","unEC",'liEC')), NA)
})

test_that("pl_vioboxplot", {
  expect_error(pl_vioboxplot(cbind(meta_data, t(as.matrix(data[features,]))), "Cluster", "PECAM1"), NA)
})

test_that("pl_averageHeatmap", {
  expect_error(pl_averageHeatmap(data, features, group_by, groups = c("brEC","unEC",'liEC'), show_rownames = T), NA)
})

test_that("pl_stackedViolinPlot", {
  expect_error(pl_stackedViolinPlot(data, features, group_by, groups = c("brEC","unEC",'liEC')), NA)
})

test_that("pl_tableHeatmap", {
  expect_error(pl_tableHeatmap(table(meta_data$Cluster, meta_data$Batch)), NA)
})


test_that("tl_crossTableEnrichment", {
  expect_error(tl_crossTableEnrichment(table(meta_data$Cluster, meta_data$Batch)), NA)
})
