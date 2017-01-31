\name{rnaCleanR-package}
\alias{rnaCleanR-package}
\alias{rnaCleanR}
\docType{package}
\title{
  Quantify and Filter putative double strand DNA from strand-specific RNA bam file
}
\description{
  This package aims to quantify and remove putative double strand DNA from a strand-specific RNA sample. There are also options and methods to plot the positive/negative proportions of all sliding windows, which allow users to have an idea of how much the sample was contaminated and the appropriate threshold to be used for filtering.
}

\details{
 The package has three main functions:
 
 - getPlot: generate plots of positive/negative proprortion and sum of reads over all sliding windows.
 
 - filterOne: filter one bam file with a given threshold.
 
 - filterMulti: filter several bam files with a given threshold.
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
bamfilein <- system.file("data","s1.chr1.bam",package = "rnaCleanR")
getPlot(bamfilein,readLength=100)
filterOne(bamfilein,bamfileout="out.bam",statfile = "out.stat",histPlot = TRUE,histfile = "hist.pdf", readLength = 100, threshold = 0.7)
}
