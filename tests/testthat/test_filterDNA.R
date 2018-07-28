test_that("filterDNA returns proper DataFrame", {
    file <- system.file("extdata", "ex1.bam", package="Rsamtools")
    destination <- tempfile()
    win <- filterDNA(file,destination = destination, getWin = TRUE)
    require_cols <- c("Seq","Start","End","NbPos","NbNeg",
                      "CovPos","CovNeg","MaxCoverage","File")
    expect_true(all(require_cols %in% colnames(win)))
})

test_that("filterDNA creates destination file", {
    file <- system.file("extdata", "ex1.bam", package="Rsamtools")
    destination <- tempfile()
    filterDNA(file,destination = destination)
    expect_true(file.exists(paste0(destination,".bam")))
})

