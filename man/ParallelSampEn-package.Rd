\name{ParallelSampEn-package}
\alias{ParallelSampEn-package}
\alias{ParallelSampEn}
\docType{package}
\title{
  \packageTitle{ParallelSampEn}
}
\description{
\packageDescription{ParallelSampEn}
}
\details{
sampen(data, m=2, r=0.15, s=1, opencl=TRUE)

Function to calculate sample entropy of a time series for a specific scale.

ARGUMENTS

data: A numeric vector representing a time series.

m: A integer, length of template vector, also known as embedding dimension.

r: coefficient used to calculate similarity thresholds by r*sd(data)

scale: A integer, scale factor for coarse graining of the time series.

opencl: A logical, TRUE for OpenCL, FALSE for basic.

VALUE

the sample entroy value of specified parameters: m, r, s.
}
\author{
\packageAuthor{ParallelSampEn}

Maintainer: \packageMaintainer{ParallelSampEn}
}
\references{
Costa M., Goldberger A.L., Peng C.-K. Multiscale entropy analysis of physiologic time series. Phys Rev Lett 2002; 89:062102.

Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG, Mietus JE, Moody GB, Peng C-K, Stanley HE. PhysioBank, PhysioToolkit, and PhysioNet: Components of a New Research Resource for Complex Physiologic Signals (2003). Circulation. 101(23):e215-e220.
}
\keyword{ package }

\examples{
\dontrun{
# install.packages("ParallelSampEn_1.0.tar.gz", repos = NULL)
library("ParallelSampEn")
file <- system.file("data", "rrtest.txt", package = "ParallelSampEn")
data <- read.table(file)
print(nrow(data))

sampen.diff.len <- function(isParallel){
  for(len in 1:10)
  {
    se <- sampen(data[1:(len*10000),1], m = 2,  r=0.15, s=1, opencl = isParallel)
    print(se)
  }
}
print("Parallel:")
print(system.time(sampen.diff.len (TRUE)))
print("Basic:")
print(system.time(sampen.diff.len (FALSE)))
}
}
