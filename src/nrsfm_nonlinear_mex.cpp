#include "mex.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <boost/scoped_array.hpp>

using std::vector;
using boost::scoped_array;

//class UnitVectorParameterization : public ceres::LocalParameterization {
//  public:
//    UnitVectorParameterization(int n) : n_(n), x_(), A_() {}
//    ~UnitVectorParameterization() {}
//
//    bool Plus(const double* x, const double* delta, double* x_prime) const {
//      update(x);
//
//      // Matrix multiplication.
//      for (int i = 0; i < n_; i += 1) {
//        x_prime[i] = x[i];
//        for (int j = 0; j < n_ - 1; j += 1) {
//          x_prime[i] += A_[(n_ - 1) * i + j] * delta[j];
//        }
//      }
//
//      // Compute norm.
//      double norm = 0;
//      for (int i = 0; i < n_; i += 1) {
//        norm += x_prime[i] * x_prime[i];
//      }
//      norm = std::sqrt(norm);
//
//      // Normalize.
//      for (int i = 0; i < n_; i += 1) {
//        x_prime[i] = x_prime[i] / norm;
//      }
//    }
//
//    bool ComputeJacobian(const double* x, double* jacobian) const {
//      update(x);
//      std::copy(A_.begin(), A_.end(), jacobian);
//    }
//
//    int GlobalSize() const {
//      return n_;
//    }
//
//    int LocalSize() const {
//      return n_ - 1;
//    }
//
//  private:
//    void recompute() const {
//      const int nlhs = 3;
//      mxArray* U = mxCreateDoubleMatrix(1, 1, mxREAL);
//      mxArray* S = mxCreateDoubleMatrix(1, n_, mxREAL);
//      mxArray* V = mxCreateDoubleMatrix(n_, n_, mxREAL);
//      mxArray* plhs[nlhs] = { U, S, V };
//
//      const int nrhs = 1;
//      mxArray* X = mxCreateDoubleMatrix(1, n_, mxREAL);
//      std::copy(x_.begin(), x_.end(), mxGetPr(X));
//      mxArray* prhs[nrhs] = { X };
//
//      // Compute SVD.
//      mexCallMATLAB(nlhs, plhs, nrhs, prhs, "svd");
//
//      // Copy out n_ - 1 smallest right singular vectors.
//      A_.assign(mxGetPr(V) + n_, mxGetPr(V) + n_ * n_);
//
//      mxDestroyArray(U);
//      mxDestroyArray(S);
//      mxDestroyArray(V);
//    }
//
//    void update(const double* x) const {
//      bool need_recompute;
//
//      if (x_.empty()) {
//        x_.assign(n_, 0);
//        need_recompute = true;
//      } else {
//        need_recompute = false;
//
//        for (int i = 0; i < n_ && !need_recompute; i += 1) {
//          if (x_[i] != x[i]) {
//            need_recompute = true;
//          }
//        }
//      }
//
//      if (need_recompute) {
//        std::cerr << "Recomputing tangent plane" << std::endl;
//        recompute();
//      } else {
//        std::cerr << "Using cached tangent plane" << std::endl;
//      }
//    }
//
//    int n_;
//
//    // Cache previous tangent plane.
//    mutable vector<double> x_;
//    mutable vector<double> A_;
//};

class ProjectionResidual {
  public:
    // rotation -- Camera rotation as normalized quaternion, 4 vector
    // point -- World point, 3 vector
    template<class T>
    bool operator()(const T* const q,
                    const T* const x,
                    T* residual) {
      // Rotate point.
      T p[3];
      QuaternionRotatePoint(q, x, p);

      // Compute residual.
      for (int d = 0; d < 2; d += 1) {
        residual[d] = w_[d] - x[d];
      }

      return true;
    }

    ProjectionResidual(double x, double y) : w_() {
      w_[0] = x;
      w_[1] = y;
    }

  private:
    double w_[2];
};

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

      // Matrix multiplication.
      for (int k = 0; k < K_; k += 1) {
        for (int d = 0; d < 3; d += 1) {
          // Row major
          int dk = d * K_ + k;
          residuals[d] += basis[dk] * coeff[k];
        }
      }

      if (jacobians != NULL) {
        double* jac_basis = jacobians[0];
        double* jac_coeff = jacobians[1];

        if (jac_basis != NULL) {
          int n = 3 * 3 * K_;
          std::fill(jac_basis, jac_basis + n, 0.);
        }

        if (jac_coeff != NULL) {
          int n = 3 * K_;
          std::fill(jac_coeff, jac_coeff + n, 0.);
        }

        for (int k = 0; k < K_; k += 1) {
          for (int d = 0; d < 3; d += 1) {
            if (jac_basis != NULL) {
              // 3 x (3 x K), row major
              int ddk = (d * 3 + d) * K_ + k;
              jac_basis[ddk] = coeff[k];
            }

            if (jac_coeff != NULL) {
              // 3 x K, row major
              int dk = d * K_ + k;
              jac_coeff[dk] = basis[dk];
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

class ProjectionCostFunction : public ceres::CostFunction {
  public:
    ProjectionCostFunction(double x, double y, int K)
        : K_(K),
          point_func_(K),
          error_func_(new ProjectionErrorFunction(x, y)) {
      // Cameras
      mutable_parameter_block_sizes()->push_back(4);
      // Basis
      mutable_parameter_block_sizes()->push_back(3 * K_);
      // Coefficients
      mutable_parameter_block_sizes()->push_back(K_);

      // Residual is 2D point error.
      set_num_residuals(2);
    }

    bool Evaluate(const double* const* params,
                  double* error,
                  double** jacobians) const {
      const double* camera = params[0];
      const double* basis = params[1];
      const double* coeff = params[2];

      double* jac_error_camera = NULL;
      double* jac_error_basis = NULL;
      double* jac_error_coeff = NULL;

      if (jacobians != NULL) {
        jac_error_camera = jacobians[0];
        jac_error_basis = jacobians[1];
        jac_error_coeff = jacobians[2];
      }

      // Compute point from basis and coefficients.
      vector<double> point(3);
      scoped_array<double> jac_point_basis;
      if (jac_error_basis != NULL) {
        jac_point_basis.reset(new double[3 * 3 * K_]);
      }
      scoped_array<double> jac_point_coeff;
      if (jac_error_coeff != NULL) {
        jac_point_coeff.reset(new double[3 * K_]);
      }
      {
        const double* func_params[] = { basis, coeff };
        double* func_jacobians[2] = { jac_point_basis.get(),
            jac_point_coeff.get() };
        point_func_.Evaluate(func_params, point.data(), func_jacobians);
      }

      // Evaluate projection error of this point.
      scoped_array<double> jac_error_point;
      if (jac_error_basis != NULL || jac_error_coeff != NULL) {
        jac_error_point.reset(new double[2 * 3]);
      }
      {
        const double* func_params[] = { camera, point.data() };
        double* func_jacobians[] = { jac_error_camera, jac_error_point.get() };
        error_func_.Evaluate(func_params, error, func_jacobians);
      }

      // Compose Jacobians using the chain rule.
      if (jac_error_basis != NULL) {
        for (int i = 0; i < 2; i += 1) {
          for (int j = 0; j < 3; j += 1) {
            for (int k = 0; k < 3 * K_; k += 1) {
              // Row major
              int ik = (i * 3 * K_) + k;
              int ij = (i * 3) + j;
              int jk = (j * 3 * K_) + k;
              jac_error_basis[ik] = jac_error_point[ij] * jac_point_basis[jk];
            }
          }
        }
      }
      if (jac_error_coeff != NULL) {
        for (int i = 0; i < 2; i += 1) {
          for (int j = 0; j < 3; j += 1) {
            for (int k = 0; k < K_; k += 1) {
              // Row major
              int ik = (i * K_) + k;
              int ij = (i * 3) + j;
              int jk = (j * K_) + k;
              jac_error_coeff[ik] = jac_error_point[ij] * jac_point_coeff[jk];
            }
          }
        }
      }

      return true;
    }

  private:
    typedef ceres::AutoDiffCostFunction<ProjectionErrorFunction, 2, 4, 3>
        AutoDiffProjectionError;

    int K_;
    BasisFunction point_func_;
    AutoDiffProjectionError error_func_;
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
                    int K) {
  ceres::Problem problem;

  for (int t = 0; t < F; t += 1) {
    for (int i = 0; i < P; i += 1) {
      int it = t * P + i;
      const double* w = &W[it * 2];
      double* q = &Q[t * 4];
      double* b = &B[i * (3 * K)];
      double* c = &C[t * K];

      ceres::CostFunction* function = new ProjectionCostFunction(w[0], w[1], K);
      problem.AddResidualBlock(function, NULL, q, b, c);
    }
  }

  std::cout << "Problem has " << problem.NumParameters() << " parameters in " <<
      problem.NumParameterBlocks() << " blocks" << std::endl;
  std::cout << "Problem has " << problem.NumResiduals() << " residuals in " <<
      problem.NumResidualBlocks() << " blocks" << std::endl;

  ceres::Solver::Options options;
  options.max_num_iterations = 1000;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = true;

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
  if (nrhs < 4) {
    mexErrMsgTxt("Not enough input arguments");
  }
  if (nrhs > 4) {
    mexErrMsgTxt("Too many input arguments");
  }
  if (nlhs < 3) {
    mexErrMsgTxt("Not enough output arguments");
  }
  if (nlhs > 3) {
    mexErrMsgTxt("Too many output arguments");
  }

  const mxArray* projections = prhs[0];
  const mxArray* cameras_init = prhs[1];
  const mxArray* bases_init = prhs[2];
  const mxArray* coeffs_init = prhs[3];

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

  nrsfmNonlinear(W, Q, B, C, F, P, K);
}
