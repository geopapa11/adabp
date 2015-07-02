//
// Matlab MEX setdiff for sorted uint32 arrays
//

#include "mex.h"
#include <stdint.h>

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  if (nrhs != 2) mexErrMsgTxt("Takes 2 input arguments");
  int N1=mxGetN(prhs[0]); int N2=mxGetN(prhs[1]);
  if (!mxIsUint32(prhs[0]) && N1>0) mexErrMsgTxt("Arguments must be sorted uint32s");
  if (!mxIsUint32(prhs[1]) && N2>0) mexErrMsgTxt("Arguments must be sorted uint32s");

  uint32_t *src, *cmp, *dest;
  src = (uint32_t *) mxGetData(prhs[0]); 
  cmp = (uint32_t *) mxGetData(prhs[1]);
  
  plhs[0] = mxCreateNumericMatrix(1,N1,mxUINT32_CLASS,0); 
  dest=(uint32_t *) mxGetData(plhs[0]);

  unsigned int i,j,k;
  for (i=0,j=0,k=0; i<N1 && j<N2; ) {
    if (src[i] < cmp[j]) dest[k++]=src[i++]; 
    else if (src[i]==cmp[j]) { i++; }
    else if (src[i] > cmp[j] && j<N2) j++;
  }
  while (i<N1) dest[k++]=src[i++];
  mxSetN(plhs[0],k);
  
}
