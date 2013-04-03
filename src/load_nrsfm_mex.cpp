#include "mex.h"
#include <iostream>
#include <string>
#include <fstream>
#include <google/protobuf/repeated_field.h>
#include "nrsfm.pb.h"

using std::string;
using std::fstream;
using google::protobuf::RepeatedFieldBackInserter;

// Q -- 4 x F column major
// B -- 3 x K x P column major
// C -- K x F column major
void loadNrsfm(const string& filename,
               mxArray*& quaternions,
               mxArray*& basis,
               mxArray*& coeff) {
  fstream file(filename.c_str(), std::ios::in | std::ios::binary);
  if (!file) {
    mexErrMsgTxt("File not found");
  }

  nrsfm::Problem problem;
  bool ok = problem.ParseFromIstream(&file);
  if (!ok) {
    mexErrMsgTxt("Could not parse problem from file");
  }

  int F = problem.num_frames();
  int P = problem.num_points();
  int K = problem.num_bases();

  // Check dimensions.
  if (problem.quaternions_size() != 4 * F) {
    mexErrMsgTxt("Incorrect number of elements in quaternions");
  }
  if (problem.basis_size() != 3 * K * P) {
    mexErrMsgTxt("Incorrect number of elements in basis");
  }
  if (problem.coeff_size() != F * K) {
    mexErrMsgTxt("Incorrect number of elements in coeff");
  }

  quaternions = mxCreateDoubleMatrix(4, F, mxREAL);
  mwSize dims[3] = { 3, K, P };
  basis = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  coeff = mxCreateDoubleMatrix(K, F, mxREAL);

  // Copy data.
  std::copy(problem.quaternions().begin(), problem.quaternions().end(),
      mxGetPr(quaternions));
  std::copy(problem.basis().begin(), problem.basis().end(), mxGetPr(basis));
  std::copy(problem.coeff().begin(), problem.coeff().end(), mxGetPr(coeff));
}

void mexFunction(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
  const int NARGIN = 1;
  const int NARGOUT = 3;

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

  const mxArray* filename_arg = prhs[0];
  if (mxIsChar(filename_arg) != 1) {
    mexErrMsgTxt("Filename must be a string");
  }
  if (mxGetM(filename_arg) != 1) {
    mexErrMsgTxt("Filename must be a row vector");
  }

  string filename(mxArrayToString(filename_arg));

  loadNrsfm(filename, plhs[0], plhs[1], plhs[2]);
}
