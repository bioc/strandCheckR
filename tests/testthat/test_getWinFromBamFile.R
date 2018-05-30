
test_that("GetWinFromBamFile Can Be Run", {
  file <- system.file("extdata","s2.sorted.bam",package = "strandCheckR")
  win <- getWinFromBamFile(file)
  expect_equivalent(win$End - win$Start + 1,rep(1000,nrow(win)))
})
