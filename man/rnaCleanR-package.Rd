\name{rnaCleanR-package}
\alias{rnaCleanR-package}
\alias{rnaCleanR}
\docType{package}
\title{
  Filter contaminated double strand DNA from strand-specific RNA
}
\description{
  This package aims to quantify and remove contaminated double strand DNA from a strand-specific RNA sample. It uses a window sliding across each chromosome, and consider the proportion of positive/negative reads in each window comparing to a given threshold, to decide wheather a window is kept or removed. Then we assigne to each read fragment a probability of being kept based on the proportion of positive/negative of the kept windows that contain the fragment. There are also options and methods to plot the positive/negative proportions of all sliding windows, which allow users to have an idea of how much the sample was contaminated and the appropriate threshold to be used.
}

\details{
 The input bam file must be sorted before and an index file must be found at the same folder. The average read length is required to adjust the size of sliding window and its step. The functions filterOne removes potential contaminated DNA from one bam file, while filterMulti is used for several bam files. The later method uses the merged alignments from all input bam files to decide which windows are kept; then each individual bam file is filtered based on the sam set of kept windows.
 The function getPlot can be used to have an idea of how much the sample was contaminated and which threshold to be used. You can also have these plots while filtering the bam file using function filterOne with options histPlot=TRUE or winPlot=TRUE.
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
filterOne(bamfilein,bamfileout="out.bam",statfile = "out.stat",histPlot = TRUE,histfile = "hist.pdf", readLength = 100, threshold = 0.7)
}
