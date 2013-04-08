#ifndef LOW_RANK_NRSFM_HPP_
#define LOW_RANK_NRSFM_HPP_

#include <ceres/ceres.h>

class ComposeLowRankStructure : public ceres::CostFunction {
  public:
    ComposeLowRankStructure(int K);

    bool Evaluate(const double* const* parameters,
                  double* residuals,
                  double** jacobians) const;

  private:
    int K_;
};

class ProjectionErrorFunction {
  public:
    ProjectionErrorFunction(double x, double y);

    template<class T>
    bool operator()(const T* const camera,
                    const T* const object_point,
                    T* error) const;

  private:
    double w_[2];
};

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
           bool check_gradients);

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
                   bool check_gradients);

#include "low-rank-nrsfm.inl"

#endif
