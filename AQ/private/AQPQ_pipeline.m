function [CODES, CODEBOOKS, DISTORTIONS, R] = AQPQ_pipeline( ...
    X, NWORDS, NITQ_OPQ, NITS_AQ, N, V )
%AQPQ_PIPELINE Quantizes according to the hybrid APQ algorithm proposed by 
% Babenko and Lempitsky at CVPR 14.
%
% [CODES, CODEBOOKS, DISTORTIONS, R] = AQPQ_pipeline( ...
%   X, NWORDS, NITQ_OPQ, NITS_AQ, N, V )
%   
% Input
%   X      : d-by-n matrix. Each column is a datapoint to encode.
%   NWORDS : M-long vector. It has as many entries as codebooks are
%            desired. The entries specify the number of codewords to use 
%            for each codebook.
%   NITS_OPQ : Integer. Number of iterations to run OPQ for.
%   NITS_AQ  : Integer. Number of iterations to run AQ for.
%   N        : Integer. Beam search depth.
%   V        : Booleam. Whether to print progress. Defaults to false.
%
% Output
%   CODES       : NITS+1-long cell array. The first entry contains the
%                 initialization (random) codes. Entry i+1 contains the
%                 codes after the ith iteration.
%   CODEBOOKS   : NITS+1-long cell array. the first entry contains the
%                 initialization codebooks. Entry i+1 contains the
%                 codebooks after the ith iteration.
%   DISTORTIONS : NITS+1-long vector. The first entry contains the
%                 distorition after initialization. Entry i+1 contains the
%                 distortion after iteration i;
%   R           : d-by-d orthogonal matrix. Learned Rotation from OPQ.

% --
% Julieta

SUBAQ = 4;

M     = numel( NWORDS );
nsub  = M / SUBAQ;

% Create the output.
CODES       = cell(  NITS_AQ+1, 1 );
CODEBOOKS   = cell(  NITS_AQ+1, 1 );
DISTORTIONS = zeros( NITS_AQ+1, 1 );

%% === Initialize with OPQ ===
[model, B, ~] = ckmeans(X, M, NWORDS, NITQ_OPQ, 'natural', false);

R               = model.R;
codesopq        = B;
codebooks_opqnp = model.centers';

% Convert the OPQ codebooks to groups processable in parallel by AQ.
AQ_CODEBOOKS  = OPQ2AQcodebooks( codebooks_opqnp );

% Apply the rotation to the training data.
X   = R' * X;

%% === Update codebooks ===
CB = AQPQ_update_codebooks( X, codesopq, AQ_CODEBOOKS, SUBAQ );

%% === Save initial values ===
CODES{1}       = codesopq;
CODEBOOKS{1}   = CB;
DISTORTIONS(1) = get_qerror( X, codesopq, CB );
if V, fprintf('=== Initialization. Qerror is %e. ===\n', DISTORTIONS(1)); end

%% === Iterate codebook update and codes update ===
for i = 1:NITS_AQ,
    
    if V, fprintf('=== Iteration %d / %d. ===\n', i, NITS_AQ); end
    
    if V, fprintf('Encoding... '); tic; end
    C      = AQPQ_encoding( X, CB, N );
    qerror = get_qerror( X, C, CB );
    if V, fprintf('done in %.2f seconds. Qerror is %e.\n', toc, qerror); end
    
    if V, fprintf('Updating codebooks... '); tic; end
    CB     = AQPQ_update_codebooks( X, C, CB, SUBAQ );
    qerror = get_qerror( X, C, CB );
    if V, fprintf('done in %.2f seconds. Qerror is %e.\n', toc, qerror); end
    
    CODES{ i+1 }       = C;
    CODEBOOKS{ i+1 }   = CB;
    
    DISTORTIONS( i+1 ) = qerror;
    
end

end