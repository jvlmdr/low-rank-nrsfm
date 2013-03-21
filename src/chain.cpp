#include "chain.hpp"
#include <glog/logging.h>

namespace chain {

namespace {

typedef vector<double> Vector;

// Multiplies A (m x p) by B (p x n) and stores the result in C (m x n).
// Assumes row major (Jacobian format of Ceres).
void multiply(const double* A,
              const double* B,
              double* C,
              int m,
              int p,
              int n) {
  std::fill(C, C + m * n, 0.);

  int ik = 0;
  for (int i = 0; i < m; i += 1) {
    int kj = 0;
    for (int k = 0; k < p; k += 1, ik += 1) {
      int ij = i * n;
      for (int j = 0; j < n; j += 1, ij += 1, kj += 1) {
        C[ij] += A[ik] * B[kj];
      }
    }
  }
}

struct SetAddressInTable {
  vector<double*>* table;

  void operator()(map<int, Vector>::reference e) {
    table->at(e.first) = &e.second.front();
  }

  SetAddressInTable(vector<double*>& table) : table(&table) {}
};

void mapToTable(map<int, Vector>& subset, vector<double*>& table) {
  std::for_each(subset.begin(), subset.end(), SetAddressInTable(table));
}

// Intermediate result of evaluating g(x).
//
// Contains the value z = g(x) and optionally the Jacobians dz_i/dx_ij with
// respect to some input values and the Jacobian dy/dz_i for the output of the
// parent function.
struct Intermediate {
  Vector z;
  Vector dy_dz;
  map<int, Vector> dz_dx;

  // g(x) maps to R^p. Note that Jacobians are empty.
  Intermediate(int p) : z(p, 0), dy_dz(), dz_dx() {}

  Intermediate() : z(), dy_dz(), dz_dx() {}

  void swap(Intermediate& other) {
    z.swap(other.z);
    dy_dz.swap(other.dy_dz);
    dz_dx.swap(other.dz_dx);
  }
};

// Functor to compose dy/dx = dy/dz dz/dx.
struct ComposeJacobian {
  double** dy_dx;
  const double* dy_dz;
  int m;
  int p;

  void operator()(map<int, Vector>::const_reference e) {
    int j = e.first;
    const Vector& dz_dx = e.second;

    CHECK_EQ(dz_dx.size() % p, 0);
    int n = dz_dx.size() / p;
    multiply(dy_dz, &dz_dx.front(), dy_dx[j], m, p, n);
  }

  ComposeJacobian(double** dy_dx, const double* dy_dz, int m, int p)
      : dy_dx(dy_dx), dy_dz(dy_dz), m(m), p(p) {}
};

// Functor to compute dy/dx_ij = dy/dz_i dz_i/dx_ij for all i, j requested.
struct ComposeJacobians {
  double** dy_dx;
  const vector<int>* offsets;

  void operator()(map<int, Intermediate>::const_reference e) {
    int i = e.first;
    const Intermediate& intermediate = e.second;
    const Vector& z = intermediate.z;
    const Vector& dy_dz = intermediate.dy_dz;
    const map<int, Vector>& dz_dx = intermediate.dz_dx;

    if (!dy_dz.empty()) {
      int p = z.size();
      CHECK_EQ(dy_dz.size() % p, 0);
      int m = dy_dz.size() / p;

      ComposeJacobian compose(dy_dx + (*offsets)[i], &dy_dz.front(), m, p);
      std::for_each(dz_dx.begin(), dz_dx.end(), compose);
    }
  }

  ComposeJacobians(double** dy_dx, const vector<int>& offsets)
      : dy_dx(dy_dx), offsets(&offsets) {}
};

}

////////////////////////////////////////////////////////////////////////////////

ComposedCostFunction::ComposedCostFunction(ceres::CostFunction* function)
    : f_(function), g_() {
  set_num_residuals(f_->num_residuals());
  // Default to no input functions (all NULL).
  int num_blocks = f_->parameter_block_sizes().size();
  g_.assign(num_blocks, shared_ptr<ceres::CostFunction>());
  updateParameterBlocks();
}

ComposedCostFunction::~ComposedCostFunction() {}

bool ComposedCostFunction::Evaluate(const double* const* x,
                                    double* y,
                                    double** dy_dx) const {
  typedef vector<ceres::int16> SizeList;

  // Some z are passed through and some are computed.
  vector<const double*> z;
  map<int, Intermediate> intermediates;
  // x[offsets[i] + j] = x[ij]
  vector<int> offsets;

  // Pointers to Jacobian outputs.
  vector<double*> dy_dz;
  if (dy_dx != NULL) {
    int num_blocks = f_->parameter_block_sizes().size();
    dy_dz.assign(num_blocks, NULL);
  }

  // Global parameter index.
  int ij = 0;

  for (int i = 0; i < int(g_.size()); i += 1) {
    CHECK_LT(ij, parameter_block_sizes().size());
    offsets.push_back(ij);

    if (!g_[i]) {
      // This parameter is not computed using a function. Pass through.
      z.push_back(x[ij]);
      // Likewise for Jacobian.
      if (dy_dx != NULL) {
        dy_dz.at(i) = dy_dx[ij];
      }
      // Move to next global parameter.
      ij += 1;
    } else {
      ceres::CostFunction& gi = *g_[i];
      const double* const* xi = &x[ij];

      // Number of outputs.
      int p = gi.num_residuals();
      // List of parameter block sizes for this function.
      const SizeList& blocks = gi.parameter_block_sizes();

      Intermediate intermediate(p);
      // For convenience.
      Vector& zi = intermediate.z;
      Vector& dy_dzi = intermediate.dy_dz;
      map<int, Vector>& dzi_dx = intermediate.dz_dx;

      // Add output of this function as input to parent function.
      z.push_back(&zi.front());

      // Allocate memory for Jacobians where required.
      for (int j = 0; j < int(blocks.size()); j += 1) {
        // Number of inputs.
        int n = blocks.at(j);

        if (dy_dx != NULL) {
          if (dy_dx[ij] != NULL) {
            // Computing this Jacobian. Allocate memory for p x n matrix.
            Vector dzi_dxj(p * n, 0);
            dzi_dx[j].swap(dzi_dxj);
          }
        }

        // Move to next global parameter.
        ij += 1;
      }

      // Obtain list of Jacobian addresses.
      vector<double*> dzi_dx_ptrs;
      if (dy_dx != NULL) {
        if (!dzi_dx.empty()) {
          int num_blocks = gi.parameter_block_sizes().size();
          dzi_dx_ptrs.assign(num_blocks, NULL);
          mapToTable(dzi_dx, dzi_dx_ptrs);
        }
      }

      // Evaluate input function.
      double** dzi_dx_ptrs_ptr = NULL;
      if (!dzi_dx_ptrs.empty()) {
        dzi_dx_ptrs_ptr = &dzi_dx_ptrs.front();
      }
      gi.Evaluate(xi, &zi.front(), dzi_dx_ptrs_ptr);

      // If we need to compute any Jacobians for this input function, then
      // also find the Jacobian of the parent function with respect to it.
      if (dy_dx != NULL) {
        if (!dzi_dx.empty()) {
          dy_dzi.assign(num_residuals() * p, 0);
          dy_dz.at(i) = &dy_dzi.front();
        }
      }

      // Add to list of intermediate results.
      intermediates[i].swap(intermediate);
    }
  }

  // Evaluate parent function.
  double** dy_dz_ptr = NULL;
  if (!dy_dz.empty()) {
    dy_dz_ptr = &dy_dz.front();
  }
  f_->Evaluate(&z.front(), y, dy_dz_ptr);

  // Combine Jacobians using chain rule.
  if (dy_dx != NULL) {
    std::for_each(intermediates.begin(), intermediates.end(),
        ComposeJacobians(dy_dx, offsets));
  }

  return true;
}

void ComposedCostFunction::SetInput(int i, ceres::CostFunction* gi) {
  if (gi != NULL) {
    // Ensure number of inputs equals number of outputs.
    CHECK_EQ(f_->parameter_block_sizes()[i], gi->num_residuals());
  }

  // Replace old function.
  shared_ptr<ceres::CostFunction> ptr(gi);
  g_.at(i).swap(ptr);

  updateParameterBlocks();
}

// Re-populates list of parameter block sizes.
void ComposedCostFunction::updateParameterBlocks() {
  typedef vector<ceres::int16> SizeList;

  const SizeList& parent_blocks = f_->parameter_block_sizes();
  SizeList& blocks = *mutable_parameter_block_sizes();
  blocks.clear();

  for (int i = 0; i < int(g_.size()); i += 1) {
    if (g_[i]) {
      // This parameter is computed using a function. Use child arguments.
      const SizeList& child_blocks = g_[i]->parameter_block_sizes();
      std::copy(child_blocks.begin(), child_blocks.end(),
          std::back_inserter(blocks));
    } else {
      // This parameter is fed through.
      blocks.push_back(parent_blocks.at(i));
    }
  }
}

} // namespace
