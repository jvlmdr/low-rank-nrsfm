#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <google/protobuf/repeated_field.h>
#include "nrsfm.pb.h"
#include "corrective-triple.hpp"

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

  nrsfm::CorrectiveTripleProblem problem;
  bool ok = problem.ParseFromIstream(&file);
  if (!ok) {
    std::cerr << "Could not parse problem from file" << std::endl;
    return 1;
  }

  int F = problem.num_frames();
  int K = problem.num_bases();

  // Check dimensions.
  if (problem.uncorrected_motion_size() != 6 * K * F) {
    std::cerr << "Incorrect number of elements in uncorrected motion" <<
        std::endl;
    return 1;
  }
  if (problem.corrective_triple_size() != 9 * K) {
    std::cerr << "Incorrect number of elements in corrective triple" <<
        std::endl;
    return 1;
  }

  // Copy data.
  vector<double> motion;
  vector<double> triple;
  std::copy(problem.uncorrected_motion().begin(),
      problem.uncorrected_motion().end(), std::back_inserter(motion));
  std::copy(problem.corrective_triple().begin(),
      problem.corrective_triple().end(), std::back_inserter(triple));

  // Solve!
  refineCorrectiveTriple(&motion.front(), &triple.front(), F, K, 1000, 1e-4,
      true, false);

  nrsfm::CorrectiveTripleSolution solution;
  solution.set_num_bases(K);
  // Copy data.
  std::copy(triple.begin(), triple.end(),
      RepeatedFieldBackInserter(solution.mutable_corrective_triple()));

  // Save.
  std::ofstream output(outfile.c_str(), std::ios::trunc | std::ios::binary);
  ok = solution.SerializeToOstream(&output);
  if (!ok) {
    std::cerr << "Could not save solution to file" << std::endl;
    return 1;
  }

  return 0;
}
