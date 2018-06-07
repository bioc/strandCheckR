
test_that("GetWinFromBamFile Can Be Run", {
  file <- system.file("extdata","s2.sorted.bam",package = "strandCheckR")
  win <- getWinFromBamFile(file)
  require_cols <- c("Seq","Start","End","NbPositive","NbNegative",
                    "CovPositive","CovNegative","MaxCoverage")
  expect_true(all(require_cols %in% colnames(win)))
})
