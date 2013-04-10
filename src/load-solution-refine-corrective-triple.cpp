#include "mex.h"
#include <string>
#include <fstream>
#include <google/protobuf/repeated_field.h>
#include "nrsfm.pb.h"

using std::string;
using google::protobuf::RepeatedFieldBackInserter;

// G -- 3K x 3 column major
void load(const string& filename, mxArray*& G) {
  std::ifstream file(filename.c_str(), std::ios::binary);
  if (!file) {
    mexErrMsgTxt("File not found");
  }

  nrsfm::CorrectiveTripleSolution solution;
  bool ok = solution.ParseFromIstream(&file);
  if (!ok) {
    mexErrMsgTxt("Could not parse solution from file");
  }

  int K = solution.num_bases();

  // Check dimensions.
  if (solution.corrective_triple_size() != 9 * K) {
    mexErrMsgTxt("Incorrect number of elements in corrective triple");
  }

  G = mxCreateDoubleMatrix(3 * K, 3, mxREAL);

  // Copy data.
  std::copy(solution.corrective_triple().begin(),
      solution.corrective_triple().end(), mxGetPr(G));
}

void mexFunction(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
  const int NARGIN = 1;
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

  const mxArray* filename_arg = prhs[0];
  if (mxIsChar(filename_arg) != 1) {
    mexErrMsgTxt("Filename must be a string");
  }
  if (mxGetM(filename_arg) != 1) {
    mexErrMsgTxt("Filename must be a row vector");
  }
  string filename(mxArrayToString(filename_arg));

  load(filename, plhs[0]);
}
