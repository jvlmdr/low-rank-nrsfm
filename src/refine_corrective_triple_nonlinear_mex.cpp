#include "mex.h"
#include "corrective-triple.hpp"

void checkDimensions(const mxArray* M, const mxArray* G) {
  bool valid = true;
  mwSize num_dims;
  const mwSize* dims;

  int F;
  int K;

  num_dims = mxGetNumberOfDimensions(M);
  if (num_dims != 3) {
    valid = false;
  } else {
    dims = mxGetDimensions(M);
    if (dims[0] != 2) {
      valid = false;
    } else {
      if (dims[1] % 3 != 0) {
        valid = false;
      }
      K = dims[1] / 3;
      F = dims[2];
    }
  }
  if (!valid) {
    mexErrMsgTxt("M should be 2 x 3K x F");
  }

  num_dims = mxGetNumberOfDimensions(G);
  if (num_dims != 2) {
    valid = false;
  } else {
    dims = mxGetDimensions(G);
    if (dims[0] != 3 * K) {
      valid = false;
    } else if (dims[1] != 3) {
      valid = false;
    }
  }
  if (!valid) {
    mexErrMsgTxt("G should be 3K x 3");
  }
}

void mexFunction(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
  const int NARGIN = 4;
  const int NARGOUT = 1;

  if (nrhs < NARGIN) {
    mexErrMsgTxt("Not enough input arguments");
  }
  if (nrhs > NARGIN) {
    mexErrMsgTxt("Too many input arguments");
  }
  if (nlhs < NARGOUT) {
    mexErrMsgTxt("Not enough output arguments");
  }
  if (nlhs > NARGOUT) {
    mexErrMsgTxt("Too many output arguments");
  }

  const mxArray* M = prhs[0];
  const mxArray* G_init = prhs[1];
  int max_iter = mxGetScalar(prhs[2]);
  double tol = mxGetScalar(prhs[3]);

  checkDimensions(M, G_init);
  int F = mxGetDimensions(M)[2];
  int K = mxGetDimensions(M)[1] / 3;

  mxArray* G = plhs[0] = mxDuplicateArray(G_init);

  refineCorrectiveTriple(mxGetPr(M), mxGetPr(G), F, K, max_iter, tol, true,
      false);
}
