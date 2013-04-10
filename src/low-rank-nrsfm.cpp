#include "low-rank-nrsfm.hpp"
#include <iostream>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include "chain.hpp"

namespace {

class ComposeLowRankStructure : public ceres::CostFunction {
  public:
    ComposeLowRankStructure(int K);

    bool Evaluate(const double* const* parameters,
                  double* residuals,
                  double** jacobians) const;

  private:
    int K_;
};

ComposeLowRankStructure::ComposeLowRankStructure(int K) : K_(K) {
  // Basis
  mutable_parameter_block_sizes()->push_back(3 * K_);
  // Coefficients
  mutable_parameter_block_sizes()->push_back(K_);

  // Residual is 2D point error.
  set_num_residuals(3);
}

bool ComposeLowRankStructure::Evaluate(const double* const* parameters,
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

////////////////////////////////////////////////////////////////////////////////

class ProjectionErrorFunction {
  public:
    ProjectionErrorFunction(double x, double y) : w_() {
      w_[0] = x;
      w_[1] = y;
    }

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

  private:
    double w_[2];
};

////////////////////////////////////////////////////////////////////////////////

class DifferenceFunction {
  public:
    DifferenceFunction(int D) : D_(D) {}

    template<class T>
    bool operator()(const T* const x, const T* const y, T* r) const {
      for (int i = 0; i < 3; i += 1) {
        r[i] = x[i] - y[i];
      }

      return true;
    }

  private:
    int D_;
};

} // namespace

////////////////////////////////////////////////////////////////////////////////

// W -- 4 x P x F column major
// Q -- 4 x F column major
// B -- 3 x K x P column major
// C -- K x F column major
void nrsfm(const double* W,
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
      function->SetInput(1, new ComposeLowRankStructure(K));

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

// W -- 4 x P x F column major
// Q -- 4 x F column major
// B -- 3 x K x P column major
// C -- K x F column major
void findStructure(const double* W,
                   const double* Q,
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
    // :/
    double* q = const_cast<double*>(&Q[t * 4]);
    double* c = &C[t * K];

    for (int i = 0; i < P; i += 1) {
      int it = t * P + i;
      const double* w = &W[it * 2];
      double* b = &B[i * (3 * K)];

      chain::ComposedCostFunction* function = new chain::ComposedCostFunction(
          new ceres::AutoDiffCostFunction<ProjectionErrorFunction, 2, 4, 3>(
            new ProjectionErrorFunction(w[0], w[1])));
      function->SetInput(1, new ComposeLowRankStructure(K));

      problem.AddResidualBlock(function, NULL, q, b, c);
    }

    problem.SetParameterBlockConstant(q);
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

// W -- 4 x P x F column major
// Q -- 4 x F column major
// S -- 3 x P x F column major
// B -- 3 x K x P column major
// C -- K x F column major
void findApproxLowRankStructure(const double* W,
                                const double* Q,
                                double* S,
                                double* B,
                                double* C,
                                double lambda,
                                int F,
                                int P,
                                int K,
                                int max_iter,
                                double tol,
                                bool verbose,
                                bool check_gradients) {
  ceres::Problem problem;

  for (int t = 0; t < F; t += 1) {
    // :/
    double* q = const_cast<double*>(&Q[t * 4]);
    double* c = &C[t * K];

    for (int i = 0; i < P; i += 1) {
      int it = t * P + i;
      const double* w = &W[it * 2];
      double* s = &S[it * 3];
      double* b = &B[i * (3 * K)];

      ceres::CostFunction* projection_residual =
          new ceres::AutoDiffCostFunction<ProjectionErrorFunction, 2, 4, 3>(
            new ProjectionErrorFunction(w[0], w[1]));
      problem.AddResidualBlock(projection_residual, NULL, q, s);

      // Difference between projected structure and low-rank structure.
      chain::ComposedCostFunction* rank_residual =
          new chain::ComposedCostFunction(
            new ceres::AutoDiffCostFunction<DifferenceFunction, 3, 3, 3>(
              new DifferenceFunction(3)));
      rank_residual->SetInput(1, new ComposeLowRankStructure(K));
      // Scale by lambda.
      ceres::LossFunction* loss = new ceres::ScaledLoss(NULL, lambda,
          ceres::TAKE_OWNERSHIP);
      problem.AddResidualBlock(rank_residual, loss, s, b, c);
    }

    problem.SetParameterBlockConstant(q);
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
