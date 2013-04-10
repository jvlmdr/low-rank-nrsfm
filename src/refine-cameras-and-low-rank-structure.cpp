#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <google/protobuf/repeated_field.h>
#include "nrsfm.pb.h"
#include "low-rank-nrsfm.hpp"

using std::vector;
using std::string;

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

  std::ifstream file(infile.c_str(), std::ios::binary);
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
  if (problem.cameras_size() != 4 * F) {
    std::cerr << "Incorrect number of elements in cameras" << std::endl;
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
  vector<double> cameras;
  vector<double> basis;
  vector<double> coeff;

  // Copy data.
  std::copy(problem.projections().begin(), problem.projections().end(),
      std::back_inserter(projections));
  std::copy(problem.cameras().begin(), problem.cameras().end(),
      std::back_inserter(cameras));
  std::copy(problem.basis().begin(), problem.basis().end(),
      std::back_inserter(basis));
  std::copy(problem.coeff().begin(), problem.coeff().end(),
      std::back_inserter(coeff));

  refineCamerasAndLowRankStructure(&projections.front(), &cameras.front(),
      &basis.front(), &coeff.front(), F, P, K, 200, 1e-4, true, false);

  nrsfm::Problem solution;
  solution.set_num_frames(F);
  solution.set_num_points(P);
  solution.set_num_bases(K);

  // Copy data.
  std::copy(cameras.begin(), cameras.end(),
      RepeatedFieldBackInserter(solution.mutable_cameras()));
  std::copy(basis.begin(), basis.end(),
      RepeatedFieldBackInserter(solution.mutable_basis()));
  std::copy(coeff.begin(), coeff.end(),
      RepeatedFieldBackInserter(solution.mutable_coeff()));

  // Save.
  std::ofstream output(outfile.c_str(), std::ios::trunc | std::ios::binary);
  ok = solution.SerializeToOstream(&output);
  if (!ok) {
    std::cerr << "Could not save solution to file" << std::endl;
    return 1;
  }

  return 0;
}
