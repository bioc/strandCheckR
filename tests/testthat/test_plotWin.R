test_that("plotWin gives a ggplot object", {
    file <- system.file("extdata", "ex1.bam", package="Rsamtools")
    win <- getWinFromBamFile(file)
    g <- plotWin(win, group_by = "Type")
    expect_is(g,"ggplot")
})

