#include "mex.h"
#include <iostream>
#include <string>
#include <fstream>
#include <google/protobuf/repeated_field.h>
#include "nrsfm.pb.h"

using std::string;
using std::fstream;
using google::protobuf::RepeatedFieldBackInserter;

// W -- 2 x P x F column major
// Q -- 4 x F column major
// B -- 3 x K x P column major
// C -- K x F column major
bool saveNrsfm(const double* W,
               const double* Q,
               const double* B,
               const double* C,
               int F,
               int P,
               int K,
               const string& filename) {
  nrsfm::Problem problem;

  problem.set_num_frames(F);
  problem.set_num_points(P);
  problem.set_num_bases(K);

  // Copy data.
  std::copy(W, W + 2 * F * P,
      RepeatedFieldBackInserter(problem.mutable_projections()));
  std::copy(Q, Q + 4 * F,
      RepeatedFieldBackInserter(problem.mutable_quaternions()));
  std::copy(B, B + 3 * K * P,
      RepeatedFieldBackInserter(problem.mutable_basis()));
  std::copy(C, C + K * F, RepeatedFieldBackInserter(problem.mutable_coeff()));

  // Save.
  fstream output(filename.c_str(),
      std::ios::out | std::ios::trunc | std::ios::binary);
  bool ok = problem.SerializeToOstream(&output);

  return ok;
}

void checkDimensions(const mxArray* projections,
                     const mxArray* quaternions,
                     const mxArray* bases,
                     const mxArray* coeffs) {
  bool valid = true;
  mwSize num_dims;
  const mwSize* dims;

  int F;
  int P;
  int K;

  num_dims = mxGetNumberOfDimensions(projections);
  if (num_dims != 3) {
    valid = false;
  } else {
    dims = mxGetDimensions(projections);
    if (dims[0] != 2) {
      valid = false;
    } else {
      P = dims[1];
      F = dims[2];
    }
  }
  if (!valid) {
    mexErrMsgTxt("projections should be 2 x P x F");
  }

  num_dims = mxGetNumberOfDimensions(quaternions);
  if (num_dims != 2) {
    valid = false;
  } else {
    dims = mxGetDimensions(quaternions);
    if (dims[0] != 4) {
      valid = false;
    } else if (dims[1] != F) {
      valid = false;
    }
  }
  if (!valid) {
    mexErrMsgTxt("quaternions should be 4 x F");
  }

  num_dims = mxGetNumberOfDimensions(bases);
  if (num_dims != 3) {
    valid = false;
  } else {
    dims = mxGetDimensions(bases);
    if (dims[0] != 3) {
      valid = false;
    } else if (dims[2] != P) {
      valid = false;
    } else {
      K = dims[1];
    }
  }
  if (!valid) {
    mexErrMsgTxt("bases should be 3 x K x P");
  }

  num_dims = mxGetNumberOfDimensions(coeffs);
  if (num_dims != 2) {
    valid = false;
  } else {
    dims = mxGetDimensions(coeffs);
    if (dims[0] != K) {
      valid = false;
    } else if (dims[1] != F) {
      valid = false;
    }
  }
  if (!valid) {
    mexErrMsgTxt("coeffs should be K x F");
  }
}

void mexFunction(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
  const int NARGIN = 5;

  if (nrhs < NARGIN) {
    mexErrMsgTxt("Not enough input arguments");
  }
  if (nrhs > NARGIN) {
    mexErrMsgTxt("Too many input arguments");
  }

  const mxArray* projections = prhs[0];
  const mxArray* quaternions = prhs[1];
  const mxArray* bases = prhs[2];
  const mxArray* coeffs = prhs[3];
  checkDimensions(projections, quaternions, bases, coeffs);

  const mxArray* filename_arg = prhs[4];

  if (mxIsChar(filename_arg) != 1) {
    mexErrMsgTxt("Filename must be a string");
  }
  if (mxGetM(filename_arg) != 1) {
    mexErrMsgTxt("Filename must be a row vector");
  }

  string filename(mxArrayToString(filename_arg));

  int num_frames = mxGetDimensions(projections)[2];
  int num_points = mxGetDimensions(projections)[1];
  int num_bases = mxGetDimensions(bases)[1];

  bool ok = saveNrsfm(mxGetPr(projections), mxGetPr(quaternions),
      mxGetPr(bases), mxGetPr(coeffs), num_frames, num_points, num_bases,
      filename);

  if (!ok) {
    mexErrMsgTxt("Failed to write to file");
  }
}
