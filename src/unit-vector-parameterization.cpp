#include "unit-vector-parameterization.hpp"

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::RowMajor;
using Eigen::VectorXd;

UnitVectorParameterization::UnitVectorParameterization(int n) : n_(n) {}

bool UnitVectorParameterization::Plus(const double* x_data,
                                      const double* p_data,
                                      double* y_data) const {
  Eigen::Map<const VectorXd> x(x_data, n_);
  Eigen::Map<const VectorXd> p(p_data, n_);
  Eigen::Map<VectorXd> y(y_data, n_);

  VectorXd z = x + p;
  y = 1. / z.norm() * z;

  return true;
}

bool UnitVectorParameterization::ComputeJacobian(const double* x_data,
                                                 double* dy_dx_data) const {
  // Ceres uses row-major for Jacobians.
  typedef Matrix<double, Dynamic, Dynamic, RowMajor> Jacobian;
  Eigen::Map<const VectorXd> x(x_data, n_);
  // Jacobian is n x n.
  Eigen::Map<Jacobian> dy_dx(dy_dx_data, n_, n_);
  dy_dx.setZero();

  double norm_x = x.norm();
  double norm_x3 = norm_x * norm_x * norm_x;

  for (int i = 0; i < n_; i += 1) {
    for (int j = 0; j < n_; j += 1) {
      if (i == j) {
        dy_dx(i, j) = 1. / norm_x - x(i) * x(i) / norm_x3;
      } else {
        dy_dx(i, j) = -x(i) * x(j) / norm_x3;
      }
    }
  }

  return true;
}

int UnitVectorParameterization::GlobalSize() const {
  return n_;
}

int UnitVectorParameterization::LocalSize() const {
  return n_;
}
