\name{sampen}
\alias{sampen}
\title{
 Function to calculate sample entropy of a time series for a specific scale.
}
\description{
 Function to calculate sample entropy of a time series for a specific scale.
}
\usage{
sampen(data, m=2, r=0.15, s=1, opencl=TRUE)	
}
\arguments{
  \item{data}{
A numeric vector representing a time series.
}
  \item{m}{
A integer, length of template vector, also known as embedding dimension.
}
  \item{r}{
coefficient used to calculate similarity thresholds by r*sd(data)
}
  \item{scale}{
A integer, scale factor for coarse graining of the time series.
}
  \item{opencl}{
A logical, TRUE for OpenCL, FALSE for basic.
}
}
\value{
The sample entroy value of specified parameters: m, r, s.
}
\examples{
\dontrun{
library("ParallelSampEn")
file <- system.file("data", "rrtest.txt", package = "ParallelSampEn")
data <- read.table(file)
print(nrow(data))
print(sampen(data[1:10000,1], m = 2,  r=0.15, s=1, opencl = TRUE))
}
}
