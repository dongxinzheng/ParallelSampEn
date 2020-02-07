# install.packages("ParallelSampEn_1.0.tar.gz", repos = NULL)
library("ParallelSampEn")
file <- system.file("data", "rrtest.txt", package = "ParallelSampEn")
data <- read.table(file)
print(nrow(data))
print(sampen(data[1:10000,1], m = 2,  r=0.15, s=1, opencl = TRUE))

sampen.diff.len <- function(isParallel){
  for(len in 1:2)
  {
    se <- sampen(data[1:(len*10000),1], m = 2,  r=0.15, s=1, opencl = isParallel)
    print(se)
  }
}
print("Parallel:")
print(system.time(sampen.diff.len (TRUE)))
print("Basic:")
print(system.time(sampen.diff.len (FALSE)))

