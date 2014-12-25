function [ CODES, CODEBOOKS ] = SQ_initialize_kmeans( X, NWORDS, NLEVELS, NITS, V )
%SQ_INITIALIZE_KMEANS Initialize a series of stacked quantizers.
%
% [ CODES, CODEBOOKS ] = SQ_INITIALIZE_KMEANS( X, NWORDS, NLEVELS, V )
%
% Input
%   X       : d-by-n matrix. Each row is a datapoint.
%   NWORDS  : Integer. Number of words per codebook.
%   NLEVELS : Integer. Number of codebooks.
%   V       : Boolean. Whether to print progress. Optional, defaults to
%               false.
%
% Output
%   CODES     : nlevels-by-n matrix. Each row is the encoding at the jth
%                 level.
%   CODEBOOKS : nlevels-long cell array. Each entry contains an 
%                 d-by-nwords matrix of codewords used to encode the 
%                 database.

% --
% Julieta

[~, n] = size( X );

if nargin < 5,
    V = false;
end

% Save the codewords here.
CODEBOOKS = cell( NLEVELS, 1);
CODES     = zeros( NLEVELS, n);

% Do k-means, substract, and do k-means again.
for j = 1:NLEVELS,

    % Do k-means.
    if V, fprintf('Running kmeans on level %d... ', j); tic; end

    % Using quick-and-dirty k-means implementation here.
    [centers, idx] = jlt_kmeans( single(X), NWORDS, NITS );
    
    CODES(j, :) = idx;

    if V, fprintf('done in %.2f seconds. ', toc); end

    % Save the ith codebook.
    CODEBOOKS{j} = centers;

    % Substract the encoding to update the residual.
    X = X - centers(:,idx);

    if V, fprintf('Qerror is %e\n', mean(sum(X.^2))); end

end

end

function [centers, idx] = jlt_kmeans( X, k, niters )
% Quick and dirty k-means with random initialization by Julieta.
[~, n] = size( X );

% Initialize codes at random.
idx = int32( ceil( rand( n, 1 ) * k ));

% Update centers and assignments.
for i = 1:niters,
    centers = kmeans_iter_mex( X, idx, k );
    idx     = euc_nn_mex(centers, X);
end

end

