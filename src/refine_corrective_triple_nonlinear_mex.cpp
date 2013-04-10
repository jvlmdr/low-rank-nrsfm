#include "mex.h"
#include <iostream>
#include <algorithm>
#include <ceres/ceres.h>
#include "chain.hpp"
#include "unit-vector-parameterization.hpp"

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::RowMajor;
using Eigen::ColMajor;

namespace {

class GramOfCorrectedRowPair : public ceres::CostFunction {
  public:
    GramOfCorrectedRowPair(int K, const double* row_pair);

    bool Evaluate(const double* const* parameters,
                  double* residuals,
                  double** jacobians) const;

  private:
    int K_;
    const double* row_pair_;
};

GramOfCorrectedRowPair::GramOfCorrectedRowPair(int K, const double* row_pair)
    : K_(K), row_pair_(row_pair) {
  // Corrective triple (3K x 3)
  mutable_parameter_block_sizes()->push_back(9 * K_);
  // Outputs 2x2 Gram matrix
  set_num_residuals(4);
}

bool GramOfCorrectedRowPair::Evaluate(const double* const* parameters,
                                      double* residuals,
                                      double** jacobians) const {
  typedef Matrix<double, Dynamic, Dynamic, ColMajor> MatlabMatrix;

  Eigen::Map<const MatlabMatrix> M(row_pair_, 2, 3 * K_);
  Eigen::Map<const MatlabMatrix> G(parameters[0], 3 * K_, 3);

  // Choice of X and Y being col-major or row-major is arbitrary.
  // Ceres uses row major (for Jacobians), Matlab column major.
  Matrix<double, 2, 3, ColMajor> X = M * G;

  Eigen::Map<Matrix<double, 2, 2, ColMajor> > Y(residuals);
  Y = X * X.transpose();

  if (jacobians != NULL) {
    typedef Matrix<double, 4, Dynamic, RowMajor> JacobianMatrix;

    if (jacobians[0] != NULL) {
      // (3K x 3) -> (2 x 2)
      Eigen::Map<JacobianMatrix> dYdG(jacobians[0], 4, 9 * K_);
      dYdG.setZero();

      for (int i = 0; i < 2; i += 1) {
        for (int j = 0; j < 2; j += 1) {
          // Y is column major 2x2
          int ij = i + j * 2;
          for (int u = 0; u < 3 * K_; u += 1) {
            for (int v = 0; v < 3; v += 1) {
              // G is column major 3K x 3
              int uv = u + v * (3 * K_);
              dYdG(ij, uv) = M(i, u) * X(j, v) + M(j, u) * X(i, v);
            }
          }
        }
      }
    }
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////

class ScaledOrthographicErrorFunction {
  public:
    template<class T>
    bool operator()(const T* const M, T* residual) const {
      // Difference between diagonals
      residual[0] = M[0] - M[3];
      // Off-diagonal
      residual[1] = M[1];

      return true;
    }
};

////////////////////////////////////////////////////////////////////////////////

// M -- 2 x 3K x F column major
// G -- 3K x 3 column major
void refineCorrectiveTriple(const double* M,
                            double* G,
                            int F,
                            int K,
                            int max_iter,
                            double tol,
                            bool verbose,
                            bool check_gradients) {
  ceres::Problem problem;

  for (int t = 0; t < F; t += 1) {
    const double* M_t = &M[t * (6 * K)];

    chain::ComposedCostFunction* function =
        new chain::ComposedCostFunction(
          new ceres::AutoDiffCostFunction<ScaledOrthographicErrorFunction, 2, 4>(
            new ScaledOrthographicErrorFunction()));
    function->SetInput(0, new GramOfCorrectedRowPair(K, M_t));

    problem.AddResidualBlock(function, NULL, G);
  }

  problem.SetParameterization(G, new UnitVectorParameterization(9 * K));

  if (verbose) {
    std::cout << "Problem has " << problem.NumParameters() <<
        " parameters in " << problem.NumParameterBlocks() << " blocks" <<
        std::endl;
    std::cout << "Problem has " << problem.NumResiduals() << " residuals in " <<
        problem.NumResidualBlocks() << " blocks" << std::endl;
  }

  ceres::Solver::Options options;
  options.max_num_iterations = max_iter;
  options.function_tolerance = tol;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = verbose;
  options.check_gradients = check_gradients;
  options.gradient_check_relative_precision = 1e-3;

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  std::cout << summary.BriefReport() << std::endl;
}

}

////////////////////////////////////////////////////////////////////////////////

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
