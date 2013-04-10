#!/bin/bash

protoc --cpp_out=. nrsfm.proto
g++ -I/nwdata/val064/local/include \
  -fPIC -shared -o nrsfm.pb.so nrsfm.pb.cc

mex -I/nwdata/val064/local/include -L/nwdata/val064/local/lib \
  -o save_problem_refine_cameras_and_low_rank_structure \
  save-problem-refine-cameras-and-low-rank-structure.cpp nrsfm.pb.so \
  -lprotobuf

mex -I/nwdata/val064/local/include -L/nwdata/val064/local/lib \
  -o load_solution_refine_cameras_and_low_rank_structure \
  load-solution-refine-cameras-and-low-rank-structure.cpp nrsfm.pb.so \
  -lprotobuf

g++ -o refine-cameras-and-low-rank-structure \
  -I/nwdata/val064/local/include -I/nwdata/val064/local/include/eigen3 \
  -L/nwdata/val064/local/lib -L/nwdata/val064/src/ceres-solver/release/lib \
  refine-cameras-and-low-rank-structure.cpp nrsfm.pb.cc \
  low-rank-nrsfm.cpp chain.cpp \
  -lprotobuf -lceres -pthread -lglog -lcholmod -lamd -lcamd -lccolamd -lcolamd \
  -lmetis -lcxsparse -lsuitesparseconfig -lblas -llapack -lrt

mex -I/nwdata/val064/local/include -L/nwdata/val064/local/lib \
  -o load_solution_refine_low_rank_structure \
  load-solution-refine-low-rank-structure.cpp nrsfm.pb.so \
  -lprotobuf

g++ -o refine-low-rank-structure \
  -I/nwdata/val064/local/include -I/nwdata/val064/local/include/eigen3 \
  -L/nwdata/val064/local/lib -L/nwdata/val064/src/ceres-solver/release/lib \
  refine-low-rank-structure.cpp nrsfm.pb.cc low-rank-nrsfm.cpp chain.cpp \
  -lprotobuf -lceres -pthread -lglog -lcholmod -lamd -lcamd -lccolamd -lcolamd \
  -lmetis -lcxsparse -lsuitesparseconfig -lblas -llapack -lrt

mex -I/nwdata/val064/local/include -L/nwdata/val064/local/lib \
  -o save_problem_refine_corrective_triple \
  save-problem-refine-corrective-triple.cpp nrsfm.pb.so \
  -lprotobuf

mex -I/nwdata/val064/local/include -L/nwdata/val064/local/lib \
  -o load_solution_refine_corrective_triple \
  load-solution-refine-corrective-triple.cpp nrsfm.pb.so \
  -lprotobuf

g++ -o refine-corrective-triple \
  -I/nwdata/val064/local/include -I/nwdata/val064/local/include/eigen3 \
  -L/nwdata/val064/local/lib -L/nwdata/val064/src/ceres-solver/release/lib \
  refine-corrective-triple.cpp nrsfm.pb.cc \
  corrective-triple.cpp chain.cpp unit-vector-parameterization.cpp \
  -lprotobuf -lceres -pthread -lglog -lcholmod -lamd -lcamd -lccolamd -lcolamd \
  -lmetis -lcxsparse -lsuitesparseconfig -lblas -llapack -lrt
