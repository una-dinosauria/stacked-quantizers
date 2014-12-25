
% Compile utilities of OPQ by @norouzi.
cd OPQ/search/
mex linscan_aqd_knn_mex.cc linscan_aqd.cc -Iinclude CXXFLAGS="\$CXXFLAGS ...
    -fopenmp -Wall" LDFLAGS="\$LDFLAGS -fopenmp";

cd ../utils/
mex euc_nn_mex.cc CXXFLAGS="\$CXXFLAGS -fopenmp -Wall" LDFLAGS="\$LDFLAGS -fopenmp";
mex kmeans_iter_mex.cc CXXFLAGS="\$CXXFLAGS -Wall"

% Compile search for SQ/AQ.
cd ../../SQ/util/
mex SQ_search.cpp CXXFLAGS="\$CXXFLAGS -fopenmp -Wall" LDFLAGS="\$LDFLAGS -fopenmp"

% Get back to the top directory.
cd ../../