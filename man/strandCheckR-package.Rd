\name{strandCheckR-package}
\alias{strandCheckR-package}
\alias{strandCheckR}
\docType{package}
\title{
  Quantify and Filter putative double strand DNA from strand-specific RNA bam file
}
\description{
  This package aims to quantify and remove putative double strand DNA from a strand-specific RNA sample. There are also options and methods to plot the positive/negative proportions of all sliding windows, which allow users to have an idea of how much the sample was contaminated and the appropriate threshold to be used for filtering.
}

\details{
 The package has some following main functions:
 
 - getWinFromBamFile: calculate positive/negative proprortion and sum of reads over all sliding windows from a single end bam file
 
 - getWinFromPairedBamFile: calculate positive/negative proprortion and sum of reads over all sliding windows from a paried end bam file
 
 - filterDNA: filter a bam file
}
\author{
Thu-Hien To & Steve Pederson

Maintainer: Hien To <hien.to@adelaide.edu.au>
}
%\references{
%  This optional section can contain literature or other references for
%  background information.
%}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
%\keyword{ package }
%\seealso{
%  Optional links to other man pages
%}
\examples{
bamfilein <- system.file("extdata","s1.chr1.bam",package = "strandCheckR")
windows <- getWinFromBamFile(bamfilein)
plotWin(windows)
plotHist(windows)
filterDNA(file = bamfilein,fileout = "filter.bam")
}
