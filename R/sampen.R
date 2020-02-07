sampen <- function(data, m=2, r=0.15, s=1, opencl=TRUE)
{
  ############################################################################################
  # Function to calculate sample entropy of a time series for a specific scale.
  #
  # ARGUMENTS
  #  data: A numeric vector representing a time series.
  #  m: A integer, length of template vector, also known as embedding dimension.
  #  r: coefficient used to calculate similarity thresholds by r*sd(data)
  #  scale: A integer, scale factor for coarse graining of the time series.
  #  opencl: A logical, TRUE for OpenCL, FALSE for basic.
  # VALUE
  #  the sample entroy value of specified parameters: m, r, s.
  ############################################################################################
  print(paste0("len=", length(data), ", m=", m, ", r=", r, ", s=", s))
  stopifnot(length(data)>100, m>0, s>0, r>0)
  stopifnot(!anyNA(data))
  if(opencl)
    opencl_sample_entropy(data,m,r,s)
  else
    basic_sample_entropy(data,m,r,s)
}