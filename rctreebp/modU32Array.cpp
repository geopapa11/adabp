//////////////////////////////////////////////////////////////////////////
// modU32Array.cpp : Mex file for changing matlab uint32 matrix entries by reference
//
// modArray(C,i,A);  // C(i)=A, by reference
//
// WARNING! WARNING! WARNING!  This MEX file is not for the faint of heart.
//  It attempts to circumvent Matlab's memory-handling functions in order
//  to give more efficient code, and uses undocumented functions to do so.
//  It comes with no guarantees whatsoever.  It may not work on other
//  versions of Matlab, or other platforms, or even at all.  It may result
//  in memory leaks, segmentation faults, destroy data, set your computer on 
//  fire, or suck it into a black hole, for all I know.  If you do not accept 
//  these possible risks, do not use this MEX file.
//
//////////////////////////////////////////////////////////////////////////

#include "mex.h"
#include <stdint.h>

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  // verify arguments
  if (nrhs > 3 || nrhs < 3) mexErrMsgTxt("Takes 3 input arguments");
  if (nlhs != 0) mexErrMsgTxt("Outputs no results; modifies passed matrix by reference.");
 
  uint32_t* matrix = (uint32_t*) mxGetData(prhs[0]);
  uint32_t index = (uint32_t) mxGetScalar(prhs[1]);
  uint32_t val = (uint32_t) mxGetScalar(prhs[2]);
  matrix[index-1]=val;
}


