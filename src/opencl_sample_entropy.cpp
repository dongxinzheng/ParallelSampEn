#include <Rcpp.h>
#include <vector>
#include <math.h>
#include <float.h>
#include <algorithm>
#include "CL/opencl.h"
#include "oclUtils.h"
#include "helper.h"
using std::vector;
using std::sort;
using namespace Rcpp;

// declarations
double* coarse_grain(double* data, int len, int s);
double do_opencl_sample_entropy(double* data, int len, int m, double r, double sd);
/*
############################################################################################
# Function to calculate sample entropy using OpenCL
#
# ARGUMENTS
#  data: A numeric vector representing a time series.
#  m: A integer value, length of template vector.
#  r: A double value, coefficient used to calculate similarity thresholds by r*sd(data)
#  s: A integer value, scale factor for coarse graining of the time series.
# VALUE
#  the sample entroy value of specified parameters: m, r, s.
############################################################################################
*/

// [[Rcpp::export]]
NumericVector opencl_sample_entropy(NumericVector data, IntegerVector m, NumericVector r, IntegerVector s) {
  double sd_value = sd(noNA(data));
  double* cg_data = coarse_grain(data.begin(),data.length(),s[0]);
  double se = do_opencl_sample_entropy(cg_data, data.length()/s[0],m[0], r[0], sd_value);
  delete[] cg_data; 
  return NumericVector::create(se);
}


// Name of the file with the source code for the computation kernel
const char* cSourceFile = "inst/entropy.c";
const char* kernalSource = "\
__kernel void SortedMatchCount(__global const float* in, __global int2* out, const float r_new)\
{\
  int iGID = get_global_id(0);\
  int nlin_j = LEN - M;\
  if (iGID >= nlin_j)\
  {\
    return;\
  }\
  int cont[M+2]={0};\
  for (int l = iGID+1; l < nlin_j; ++l)\
  {\
    int pos1 = iGID*(M+1);\
    int pos2 = l*(M+1);\
    if(in[pos2]-in[pos1] > r_new)\
      break;\
    int k = 0;\
    cont[++k]++;\
    while (k <= M && fabs(in[pos2+k] - in[pos1+k]) <= r_new)\
      cont[++k]++;\
  }\
  out[iGID].x = cont[M];\
  out[iGID].y = cont[M+1];\
}";

// global Vars
cl_context cxGPUContext=NULL;        // OpenCL context
cl_command_queue cqCommandQueue=NULL;// OpenCL command que
cl_program cpProgram=NULL;           // OpenCL program
cl_kernel ckKernel=NULL;             // OpenCL kernel
// char* cSourceCL = NULL;         // Buffer to hold source for compilation
bool hasContext = false;
cl_int ciErr1, ciErr2;			// Error code var
cl_device_id cdDevice;
size_t szKernelLength;

// [[Rcpp::export]]
bool createContext()
{
  // has created context before
  if(hasContext)
    return hasContext;
  
  // Get the device of specified type
  ciErr1 = oclGetDeviceID(&cdDevice, CL_DEVICE_TYPE_GPU);
  if(ciErr1 != CL_SUCCESS)
    return false;
  
  //Create the context
  cxGPUContext = clCreateContext(0, 1, &cdDevice, NULL, NULL, &ciErr1);
  if(ciErr1 != CL_SUCCESS)
    return false;
  
  //Create a command-queue
  cqCommandQueue = clCreateCommandQueue(cxGPUContext, cdDevice, 0, &ciErr1);
  if(ciErr1 != CL_SUCCESS)
    return false;
  
  // Read the OpenCL kernel in from source file
  // cSourceCL = oclLoadProgSource(cSourceFile, "", &szKernelLength);
  
  // create context successfully
  hasContext = true;
  return hasContext;
}

// [[Rcpp::export]]
void cleanup ()
{
  // cleanup allocated objects
  if(ckKernel)
  {
    clReleaseKernel(ckKernel);
    ckKernel=NULL;
  }
  if(cpProgram){
    clReleaseProgram(cpProgram);
    cpProgram=NULL;
  }
  if(cqCommandQueue){
    clReleaseCommandQueue(cqCommandQueue);
    cqCommandQueue=NULL;
  }
  if(cxGPUContext){
    clReleaseContext(cxGPUContext);
    cxGPUContext=NULL;
  }
  // if(cSourceCL)
  // {
  //   free(cSourceCL);
  //   cSourceCL=NULL;
  // }
  
  hasContext=false;
  
}
// Round Up Division function
size_t roundUp(int group_size, int global_size)
{
  int r = global_size % group_size;
  if(r == 0)
  {
    return global_size;
  } else
  {
    return global_size + group_size - r;
  }
}

double getEntropy(cl_int2* pfData, int iSize)
{
  int i;
  long long a=0,b=0;
  for (i = 0; i < iSize; ++i)
  {
    a = a + pfData[i].y;
    b = b + pfData[i].x;
  }
  
  double se;
  if (a == 0 || b == 0)
    se = -log((double)1/((iSize)*(iSize-1)));
  else
    se = -log((double)a/b);
  
  return se;
}

float* sortData(const double* data, int len, int m, double r_new)
{
  int i,j;
  int nlin_j = len - m;
  struct Item *sorted;
  if ((sorted = (Item*)calloc (nlin_j, sizeof (struct Item))) == NULL) {
    Rprintf("insufficient memory\n");
    exit(2);
  }
  
  float* sortedVector;
  if ((sortedVector = (float*)calloc (nlin_j*(m+1), sizeof (float))) == NULL) {
    Rprintf("insufficient memory\n");
    exit(2);
  }
  
  // get the total number of buckets(bucket width = r*sd)
  double minV = DBL_MAX;
  double maxV = -1.0;
  for(i=0;i<nlin_j;++i)
  {
    if(minV > data[i])
      minV = data[i];
    else if(maxV < data[i])
      maxV = data[i];
  }
  int bucketNum = ceil(2*(maxV-minV)/r_new);
  
  // create and fill buckets
  vector<Item> * bucketList = new vector<Item>[bucketNum]();
  for(int i=0;i<nlin_j;++i)
  {
    int index = (data[i]-minV)/(r_new/2);
    Item item = {data[i], i};
    bucketList[index].push_back(item);
  }
  
  // sort elements for each bucket, ascending order
  for(int i=0;i<bucketNum;++i)
  {
    sort(bucketList[i].begin(), bucketList[i].end());
  }
  
  int k=0;
  for(i=0;i<bucketNum;++i)
  {
    for(unsigned int j=0;j<bucketList[i].size();++j)
      sorted[k++] = bucketList[i][j];
  }
  
  for(i=0;i<nlin_j;i++)
  {
    for(j=0;j<=m;j++)
    {
      sortedVector[i*(m+1)+j] = data[sorted[i].pos+j];
    }
  }
  
  free(sorted);
  return sortedVector;
}
int getProgramBuildInfo(cl_program program,cl_device_id device)
{
  size_t log_size;
  char *program_log;
  /* Find size of log and print to std output */
  clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
                        0, NULL, &log_size);
  program_log = (char*) malloc(log_size+1);
  program_log[log_size] = '\0';
  clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
                        log_size+1, program_log, NULL);
  Rprintf("%s\n", program_log);
  free(program_log);
  return 0;
}


double do_opencl_sample_entropy(double* data, int len,  int m, double r, double sd)
{
  if(!createContext())
    return -1;
  
  int iNumElements = len;
  // set and log Global and Local work size dimensions
  size_t szLocalWorkSize = 128;
  size_t szGlobalWorkSize = roundUp((int)szLocalWorkSize, iNumElements);  // rounded up to the nearest multiple of the LocalWorkSize
  
  // Allocate the OpenCL buffer memory objects for source and result on the device GMEM
  cl_mem cmDevSrc = clCreateBuffer(cxGPUContext, CL_MEM_READ_ONLY, sizeof(cl_float) * szGlobalWorkSize*(m+1), NULL, &ciErr1);
  ciErr1 |= ciErr2;
  cl_mem cmDevDst = clCreateBuffer(cxGPUContext, CL_MEM_WRITE_ONLY, sizeof(cl_int2) * szGlobalWorkSize, NULL, &ciErr2);
  ciErr1 |= ciErr2;
  oclCheckError(ciErr1, CL_SUCCESS);
  
  // Create the program
  if(cpProgram){
    clReleaseProgram(cpProgram);
    cpProgram=NULL;
  }
  szKernelLength = strlen(kernalSource);
  cpProgram = clCreateProgramWithSource(cxGPUContext, 1, (const char **)&kernalSource, &szKernelLength, &ciErr1);
  
  // Build the program with 'mad' Optimization option
  //    char* flags = "-cl-fast-relaxed-math";
  char options[256];
  sprintf(options, "-D M=%d -D LEN=%d", m, iNumElements);
  ciErr1 = clBuildProgram(cpProgram, 0, NULL, options, NULL, NULL);
  if(ciErr1 != 0)
    getProgramBuildInfo(cpProgram ,cdDevice);
  oclCheckError(ciErr1, CL_SUCCESS);
  // Create the kernel
  if(ckKernel)
  {
    clReleaseKernel(ckKernel);
    ckKernel=NULL;
  }
  
  ckKernel = clCreateKernel(cpProgram, "SortedMatchCount", &ciErr1);
  oclCheckError(ciErr1, CL_SUCCESS);
  
  // Set the Argument values
  float r_new = r*sd;
  ciErr1 = clSetKernelArg(ckKernel, 0, sizeof(cl_mem), (void*)&cmDevSrc);
  ciErr1 |= clSetKernelArg(ckKernel, 1, sizeof(cl_mem), (void*)&cmDevDst);
  ciErr1 |= clSetKernelArg(ckKernel, 2, sizeof(cl_float), (void*)&r_new);
  oclCheckError(ciErr1, CL_SUCCESS);
  
  // Start Core sequence... copy input data to GPU, compute, copy results back
  cl_int2 *dst = (cl_int2 *)malloc(sizeof(cl_int2) * szGlobalWorkSize);
  float *sorted = sortData(data, len, m, r_new);
  // Asynchronous write of data to GPU device
  ciErr1 = clEnqueueWriteBuffer(cqCommandQueue, cmDevSrc, CL_FALSE, 0, sizeof(cl_float) * (len-m)*(m+1), sorted, 0, NULL, NULL);
  oclCheckError(ciErr1, CL_SUCCESS);
  // Launch kernel
  ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, ckKernel, 1, NULL, &szGlobalWorkSize, &szLocalWorkSize, 0, NULL, NULL);
  oclCheckError(ciErr1, CL_SUCCESS);
  
  // Synchronous/blocking read of results, and check accumulated errors
  ciErr1 = clEnqueueReadBuffer(cqCommandQueue, cmDevDst, CL_TRUE, 0, sizeof(cl_int2) * szGlobalWorkSize, dst, 0, NULL, NULL);
  
  double en = getEntropy((cl_int2*)dst, iNumElements-m);
  
  if(cmDevSrc)clReleaseMemObject(cmDevSrc);
  if(cmDevDst)clReleaseMemObject(cmDevDst);
  // Free host memory
  free(dst);
  free(sorted);
  
  return en;
}


