function [ CODES, DISTORTIONS ] = AQPQ_encoding( X, CODEBOOKS, N, V )
%AQPQ_ENCODING Encoding for Additive Product Quantization proposed by 
% Babenko and Lempitsky. This is very similar to AQ_ENCODING, but is
% designed to handle a large (>4) number of codebooks. In this case,
% independent subspaces are found first, and then each encoded separately
% by AQ_encoding.
%
% [ CODES, DISTORTIONS ] = AQPQ_ENCODING( X, CODEBOOKS, N )
%
% Input
%   X         : d-by-n matrix. Each column is a datapoint to encode.
%   CODEBOOKS : n-long cell array. Each entry is a d-by-k codebook.
%   N         : Integer. Beam search depth.
%   V         : Boolean. Whether to print progrees. Defaults to false.
%
% Output
%   CODES       : numel(CODEBOOKS)-by-n matrix. X encoded.
%   DISTORTIONS : n-long vector. The squared distortion of each point.
%
% This code can handle d-dimensional codeboks, and is aware that only a
% portion of them should be used (as in PQ).

% --
% Julieta

if nargin < 4,
    V = 1;
end

M = numel( CODEBOOKS );
assert( M > 4, 'AQ/PQ encoding called with less than 4 codebooks');

[~, n] = size( X );

SUBAQ = 4;         % Every 4 subgroups are processed with AQ.
nsub  = M / SUBAQ; % Number of subspaces of the codebooks.

% Convert the passed (fully-dimensional) codebooks to codebook to subgroups
% of 4 codebooks that can be encoded independently.
[ APQ_CODEBOOKS, KEEP ] = AQ2APQcodebooks( CODEBOOKS, SUBAQ );

assert( numel(KEEP) == nsub, ...
    'Mismatch between the number of subcodebooks and keep indices' );

% Create space for the outputs.
CODES = cell(nsub, 1);

if nargout > 2,
    DISTORTIONS = zeros( nsub, n );
end

% Now we quantize each subspace with AQ_ENCODING.

% Loop through the subspaces.
for i = 1:nsub,
    
    if V, fprintf('  Working on subspace %d / %d... ', i, nsub); tic; end
    
    % Single out the codebooks to use.
    subaq_idx       = (i-1)*SUBAQ+1 : i*SUBAQ;
    subaq_codebooks = APQ_CODEBOOKS( subaq_idx );
    
    % Quantize with AQ.
    if nargout > 1
        [CODES{i}, DISTORTIONS(i, :)] = AQ_encoding( X(KEEP{i}, :), subaq_codebooks, N, V );
    else
        CODES{i} = AQ_encoding( X(KEEP{i}, :), subaq_codebooks, N, V );
    end
    
    if V, fprintf('done in %.2f seconds.\n', toc); end
    
end

% Join the codes!
CODES = vertcat( CODES{:} );


end

