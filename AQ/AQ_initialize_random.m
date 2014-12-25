function [CODES, CODEBOOKS] = AQ_initialize_random( X, NWORDS )
%AQ_INITIALIZE_RANDOM Random initialization for Additive Quantization 
% proposed by Babenko and Lempitsky.
%
% CODES = AQ_INITIALIZE_RANDOM(X, NWORDS) creates random 
% codes for the passed data points in X.
%
% [CODES, CODEBOOK] = AQ_INITIALIZE_RANDOM(X, NWORDS) creates
% random codes and a corresponding codebook for the points in X.
%
% Input
%   X      : d-by-n matrix. Each column is a datapoint to encode.
%   NWORDS : M-long vector. It has as many entries as codebooks are
%            desired. The entries specify the number of codewords to use 
%            for each codebook.
%
% Output
%   CODES     : M-by-n matrix. X encoded.
%   CODEBOOKS : M-long cell array. The ith entry is a d-by-NWORDS(i)
%               codebook.
% 
%   E.g. CODES = AQ_INITIALIZE_RANDOM( X, [256, 256, 256, 256] )
%   Creates codes for 4 codebooks, each with 256 codewords.

% --
% Julieta

if isrow( NWORDS ),
    NWORDS = NWORDS';
end

m = numel( NWORDS );

[d, n] = size( X );

% Create random codes and scale by the number of words.
CODES = bsxfun(@times, rand( m, n ), NWORDS);
CODES = ceil( CODES );

if nargout > 1
    
    % Create dummy codebooks to call the standard codebook update function.
    dummy_codebooks = cell(m, 1);
    for i = 1:m,
        dummy_codebooks{i} = zeros( d, NWORDS(i) );
    end
    
    CODEBOOKS = AQ_update_codebooks(X, CODES, dummy_codebooks);
end

end