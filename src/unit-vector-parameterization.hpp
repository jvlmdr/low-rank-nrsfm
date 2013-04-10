#include <ceres/ceres.h>

class UnitVectorParameterization : public ceres::LocalParameterization {
  public:
    UnitVectorParameterization(int n);
    bool Plus(const double* x, const double* p, double* y) const;
    bool ComputeJacobian(const double* x, double* dy_dx) const;
    int GlobalSize() const;
    int LocalSize() const;

  private:
    int n_;
};
