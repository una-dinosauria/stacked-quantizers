function D = sqdist(A, B)
% SQDIST Computes squared Euclidean distance matrix.
% D = SQDIST(A, B) Computes a rectangular matrix of pairwise distances
% between vectors in A and B, both given in columns.

% Very fast implementation taken from Roland Bunschoten.

D = bsxfun(@plus, -2*A'*B, sum(B.*B,1));
D = abs(bsxfun(@plus, D', sum(A.*A,1)))';
