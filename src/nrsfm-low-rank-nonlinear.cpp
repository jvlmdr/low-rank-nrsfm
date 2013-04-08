#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <google/protobuf/repeated_field.h>
#include "nrsfm.pb.h"
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include "chain.hpp"

using std::vector;
using std::string;
using std::fstream;

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

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Not enough input arguments" << std::endl;
    return 1;
  }
  if (argc > 3) {
    std::cerr << "Too many input arguments" << std::endl;
    return 1;
  }

  string infile(argv[1]);
  string outfile(argv[2]);

  fstream file(infile.c_str(), std::ios::in | std::ios::binary);
  if (!file) {
    std::cerr << "File not found" << std::endl;
    return 1;
  }

  nrsfm::Problem problem;
  bool ok = problem.ParseFromIstream(&file);
  if (!ok) {
    std::cerr << "Could not parse problem from file" << std::endl;
    return 1;
  }

  int F = problem.num_frames();
  int P = problem.num_points();
  int K = problem.num_bases();

  // Check dimensions.
  if (problem.quaternions_size() != 4 * F) {
    std::cerr << "Incorrect number of elements in quaternions" << std::endl;
    return 1;
  }
  if (problem.basis_size() != 3 * K * P) {
    std::cerr << "Incorrect number of elements in basis" << std::endl;
    return 1;
  }
  if (problem.coeff_size() != F * K) {
    std::cerr << "Incorrect number of elements in coeff" << std::endl;
    return 1;
  }

  vector<double> projections;
  vector<double> quaternions;
  vector<double> basis;
  vector<double> coeff;

  // Copy data.
  std::copy(problem.projections().begin(), problem.projections().end(),
      std::back_inserter(projections));
  std::copy(problem.quaternions().begin(), problem.quaternions().end(),
      std::back_inserter(quaternions));
  std::copy(problem.basis().begin(), problem.basis().end(),
      std::back_inserter(basis));
  std::copy(problem.coeff().begin(), problem.coeff().end(),
      std::back_inserter(coeff));

  nrsfmNonlinear(&projections.front(), &quaternions.front(), &basis.front(),
      &coeff.front(), F, P, K, 200, 1e-4,  true, false);

  nrsfm::Problem solution;
  solution.set_num_frames(F);
  solution.set_num_points(P);
  solution.set_num_bases(K);

  // Copy data.
  std::copy(quaternions.begin(), quaternions.end(),
      RepeatedFieldBackInserter(solution.mutable_quaternions()));
  std::copy(basis.begin(), basis.end(),
      RepeatedFieldBackInserter(solution.mutable_basis()));
  std::copy(coeff.begin(), coeff.end(),
      RepeatedFieldBackInserter(solution.mutable_coeff()));

  // Save.
  fstream output(outfile.c_str(),
      std::ios::out | std::ios::trunc | std::ios::binary);
  ok = solution.SerializeToOstream(&output);
  if (!ok) {
    std::cerr << "Could not save solution to file" << std::endl;
    return 1;
  }

  return 0;
}
