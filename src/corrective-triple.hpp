#ifndef CORRECTIVE_TRIPLE_HPP_
#define CORRECTIVE_TRIPLE_HPP_

// M -- 2 x 3K x F column major
// G -- 3K x 3 column major
void refineCorrectiveTriple(const double* M,
                            double* G,
                            int F,
                            int K,
                            int max_iter,
                            double tol,
                            bool verbose,
                            bool check_gradients);

#endif
