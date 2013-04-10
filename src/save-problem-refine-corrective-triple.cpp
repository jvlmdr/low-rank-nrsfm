#include "mex.h"
#include <iostream>
#include <string>
#include <fstream>
#include <google/protobuf/repeated_field.h>
#include "nrsfm.pb.h"

using std::string;
using google::protobuf::RepeatedFieldBackInserter;

// M -- 2 x 3K x F column major
// G -- 3K x 3 column major
bool save(const double* M,
          const double* G,
          int F,
          int K,
          const string& filename) {
  nrsfm::CorrectiveTripleProblem problem;

  problem.set_num_frames(F);
  problem.set_num_bases(K);

  // Copy data.
  std::copy(M, M + 6 * K * F,
      RepeatedFieldBackInserter(problem.mutable_uncorrected_motion()));
  std::copy(G, G + 9 * K,
      RepeatedFieldBackInserter(problem.mutable_corrective_triple()));

  // Save.
  std::ofstream output(filename.c_str(), std::ios::trunc | std::ios::binary);
  bool ok = problem.SerializeToOstream(&output);

  return ok;
}

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
  const int NARGIN = 3;

  if (nrhs < NARGIN) {
    mexErrMsgTxt("Not enough input arguments");
  }
  if (nrhs > NARGIN) {
    mexErrMsgTxt("Too many input arguments");
  }

  const mxArray* M = prhs[0];
  const mxArray* G = prhs[1];
  const mxArray* filename_arg = prhs[2];

  checkDimensions(M, G);
  int F = mxGetDimensions(M)[2];
  int K = mxGetDimensions(M)[1] / 3;

  if (mxIsChar(filename_arg) != 1) {
    mexErrMsgTxt("Filename must be a string");
  }
  if (mxGetM(filename_arg) != 1) {
    mexErrMsgTxt("Filename must be a row vector");
  }
  string filename(mxArrayToString(filename_arg));

  bool ok = save(mxGetPr(M), mxGetPr(G), F, K, filename);
  if (!ok) {
    mexErrMsgTxt("Failed to write to file");
  }
}
