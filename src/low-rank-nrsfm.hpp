#ifndef LOW_RANK_NRSFM_HPP_
#define LOW_RANK_NRSFM_HPP_

void refineCamerasAndLowRankStructure(const double* W,
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

void refineLowRankStructure(const double* W,
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

void refineApproxLowRankStructure(const double* W,
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
                                  bool check_gradients);

#endif
