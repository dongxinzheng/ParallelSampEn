# ParallelSampEn
A parallel implementation of sample entropy based on OpenCL
# Install
```R
#install "remotes" if it was not installed before
#install.packages("remotes")
remotes::install_github("dongxinzheng/ParallelSampEn")
```
# Example
```R
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
```
