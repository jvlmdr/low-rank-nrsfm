#ifndef CHAIN_HPP_
#define CHAIN_HPP_

#include <vector>
#include <deque>
#include <map>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <ceres/ceres.h>

namespace chain {

using std::vector;
using std::deque;
using std::map;
using boost::scoped_ptr;
using boost::shared_ptr;

// h(x) = f(z), z = g(x)
class ComposedCostFunction : public ceres::CostFunction {
  public:
    // Takes ownership of f.
    ComposedCostFunction(ceres::CostFunction* f);
    ~ComposedCostFunction();

    bool Evaluate(const double* const* x, double* y, double** dy_dx) const;

    // Compute z(i) using gi. Takes ownership of gi.
    // Output dimensionality of gi must match that of block i of f.
    void SetInput(int i, ceres::CostFunction* gi);

  private:
    // Parent function.
    scoped_ptr<ceres::CostFunction> f_;
    // Child functions (NULL to pass through).
    vector<shared_ptr<ceres::CostFunction> > g_;

    // Re-populates list of parameter block sizes.
    void updateParameterBlocks();
};

}

#endif
