
test_that("GetWinFromBamFile Can Be Run", {
  file <- system.file("extdata","paired.bam",package = "strandCheckR")
  win <- getWinFromBamFile(file,sequences="10")
  require_cols <- c("Seq","Start","End","NbPositive","NbNegative",
                    "CovPositive","CovNegative","MaxCoverage")
  expect_true(all(require_cols %in% colnames(win)))
})
