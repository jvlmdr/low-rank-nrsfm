#include "mex.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <boost/scoped_array.hpp>
#include "chain.hpp"

using std::vector;
using boost::scoped_array;

class BasisFunction : public ceres::CostFunction {
  public:
    BasisFunction(int K) : K_(K) {
      // Basis
      mutable_parameter_block_sizes()->push_back(3 * K_);
      // Coefficients
      mutable_parameter_block_sizes()->push_back(K_);

      // Residual is 2D point error.
      set_num_residuals(3);
    }

    bool Evaluate(const double* const* parameters,
                  double* residuals,
                  double** jacobians) const {
      const double* basis = parameters[0];
      const double* coeff = parameters[1];

      // Matrix times vector.
      for (int k = 0; k < K_; k += 1) {
        for (int d = 0; d < 3; d += 1) {
          // 3 x K column major, K x 3 row major
          int dk = k * 3 + d;
          residuals[d] += basis[dk] * coeff[k];
        }
      }

      if (jacobians != NULL) {
        double* jac_basis = jacobians[0];
        double* jac_coeff = jacobians[1];

        if (jac_basis != NULL) {
          int n = 3 * 3 * K_;
          std::fill(jac_basis, jac_basis + n, 0.);

          for (int k = 0; k < K_; k += 1) {
            for (int d = 0; d < 3; d += 1) {
              // 3 x (K x 3) row major
              int dkd = (d * K_ + k) * 3 + d;
              jac_basis[dkd] += coeff[k];
            }
          }
        }

        if (jac_coeff != NULL) {
          int n = 3 * K_;
          std::fill(jac_coeff, jac_coeff + n, 0.);

          for (int k = 0; k < K_; k += 1) {
            for (int d = 0; d < 3; d += 1) {
              // Jacobian is 3 x K row major
              int dk = d * K_ + k;
              // Basis is 3 x K column major
              int kd = k * 3 + d;
              jac_coeff[dk] += basis[kd];
            }
          }
        }
      }

      return true;
    }

  private:
    int K_;
};

class ProjectionErrorFunction {
  public:
    template<class T>
    bool operator()(const T* const camera,
                    const T* const object_point,
                    T* error) const {
      // Rotate into camera co-ordinates.
      T camera_point[3];
      ceres::UnitQuaternionRotatePoint(camera, object_point, camera_point);

      // Subtract from point.
      for (int d = 0; d < 2; d += 1) {
        error[d] = w_[d] - camera_point[d];
      }

      return true;
    }

    ProjectionErrorFunction(double x, double y) : w_() {
      w_[0] = x;
      w_[1] = y;
    }

  private:
    double w_[2];
};

// W -- 2 x P x F column major
// Q -- 4 x F column major
// B -- 3 x K x P column major
// C -- K x F column major
void nrsfmNonlinear(const double* W,
                    double* Q,
                    double* B,
                    double* C,
                    int F,
                    int P,
                    int K,
                    int max_iter,
                    double tol,
                    bool verbose,
                    bool check_gradients) {
  ceres::Problem problem;

  for (int t = 0; t < F; t += 1) {
    double* q = &Q[t * 4];
    double* c = &C[t * K];

    for (int i = 0; i < P; i += 1) {
      int it = t * P + i;
      const double* w = &W[it * 2];
      double* b = &B[i * (3 * K)];

      chain::ComposedCostFunction* function = new chain::ComposedCostFunction(
          new ceres::AutoDiffCostFunction<ProjectionErrorFunction, 2, 4, 3>(
            new ProjectionErrorFunction(w[0], w[1])));
      function->SetInput(1, new BasisFunction(K));

      problem.AddResidualBlock(function, NULL, q, b, c);
    }

    problem.SetParameterization(q, new ceres::QuaternionParameterization());
  }

  std::cout << "Problem has " << problem.NumParameters() << " parameters in " <<
      problem.NumParameterBlocks() << " blocks" << std::endl;
  std::cout << "Problem has " << problem.NumResiduals() << " residuals in " <<
      problem.NumResidualBlocks() << " blocks" << std::endl;

  ceres::Solver::Options options;
  options.max_num_iterations = max_iter;
  options.function_tolerance = tol;
  options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  options.minimizer_progress_to_stdout = verbose;
  options.check_gradients = check_gradients;
  options.gradient_check_relative_precision = 1e-3;

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  std::cout << summary.BriefReport() << std::endl;
}

void checkDimensions(const mxArray* projections,
                     const mxArray* cameras,
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

  num_dims = mxGetNumberOfDimensions(cameras);
  if (num_dims != 2) {
    valid = false;
  } else {
    dims = mxGetDimensions(cameras);
    if (dims[0] != 4) {
      valid = false;
    } else if (dims[1] != F) {
      valid = false;
    }
  }
  if (!valid) {
    mexErrMsgTxt("cameras should be 4 x F");
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
  const int NARGIN = 6;
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

  const mxArray* projections = prhs[0];
  const mxArray* cameras_init = prhs[1];
  const mxArray* bases_init = prhs[2];
  const mxArray* coeffs_init = prhs[3];
  int max_iter = mxGetScalar(prhs[4]);
  double tol = mxGetScalar(prhs[5]);

  checkDimensions(projections, cameras_init, bases_init, coeffs_init);

  int F = mxGetDimensions(projections)[2];
  int P = mxGetDimensions(projections)[1];
  int K = mxGetDimensions(bases_init)[1];

  mxArray* cameras = plhs[0] = mxDuplicateArray(cameras_init);
  mxArray* bases = plhs[1] = mxDuplicateArray(bases_init);
  mxArray* coeffs = plhs[2] = mxDuplicateArray(coeffs_init);

  const double* W = mxGetPr(projections);
  double* Q = mxGetPr(cameras);
  double* B = mxGetPr(bases);
  double* C = mxGetPr(coeffs);

  nrsfmNonlinear(W, Q, B, C, F, P, K, max_iter, tol, true, false);
}
