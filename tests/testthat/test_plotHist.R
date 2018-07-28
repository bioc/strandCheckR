test_that("plotHist gives a ggplot object with correct data to plot", {
    file <- system.file("extdata", "ex1.bam", package="Rsamtools")
    win <- getWinFromBamFile(file)
    g <- plotHist(win, group_by = "Type")
    expect_is(g,"ggplot")
    expect_equal(names(g$data),c("Type","PosProp","Coverage","ReadCountProp"))
})
