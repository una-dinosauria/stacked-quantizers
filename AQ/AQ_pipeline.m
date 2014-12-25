function [CODES, CODEBOOKS, DISTORTIONS, R] = AQ_pipeline( ...
    X, NWORDS, NITS, N, V )
%AQ_PIPELINE Quantizes X according to the algorithm proposed by Babenko and
%Lempitsky at CVPR 14.
%
% [CODES, CODEBOOKS, DISTORTIONS] = AQ_PIPELINE( X, NWORDS, NITS, N, V )
%
% Input
%   X      : d-by-n matrix. Each column is a datapoint to encode.
%   NWORDS : M-long vector. It has as many entries as codebooks are
%            desired. The entries specify the number of codewords to use 
%            for each codebook.
%   NITS   : Integer. Number of iterations to run the optimization for.
%   N      : Integer. Beam search depth.
%   V      : Booleam. Whether to print progress. Defaults to false.
%
% Output
%   CODES       : numel(NWORDS)-by-n matrix. X encoded.
%   CODEBOOKS   : numel(NWORDS)-long cell array. Each entry is a learned
%                 subcodebook.
%   DISTORTIONS : NITS+1-long vector. The first entry contains the
%                 distorition after initialization. Entry i+1 contains the
%                 distortion after iteration i;
%   R           : d-by-d orthogonal matrix. The rotation of the data. 
%                 If the number of codebooks, M <= 4, then this will be the
%                 identity. If M > 4 then this will correspond to the
%                 preprocessing of OPQ.

% --
% Julieta.

if nargin < 5,
    V = false;
end

NITQ_OPQ = 100;

% If we are asked for more than 4 codebooks then do APQ.
if numel( NWORDS ) > 4,
    [CODES, CODEBOOKS, DISTORTIONS, R] = AQPQ_pipeline( ...
        X, NWORDS, NITQ_OPQ, NITS, N, V );
    return
else
    R = eye( size(X,1) );
end

%% Create space for the output.
DISTORTIONS = zeros( NITS+1, 1 );

%% Initialize at random.
[CODES, CODEBOOKS] = AQ_initialize_random( X, NWORDS );

qerror         = get_qerror( X, CODES, CODEBOOKS );
DISTORTIONS(1) = qerror;

if V, fprintf('Error after initialization  is %e.\n', qerror); end;

%% Iterate between encoding and updating codebooks.
for i = 1:NITS,
    
    if V, fprintf('=== Iteration %d / %d. ===\n', i, NITS); end
    
    if V, fprintf('Encoding... '); tic; end
    CODES  = AQ_encoding(X, CODEBOOKS, N, V > 1);
    qerror = get_qerror( X, CODES, CODEBOOKS );
    if V, fprintf('done in %.2f seconds. Qerror is %e\n', toc, qerror); end
    
    if V, fprintf('Updating codebooks... '); tic; end
    CODEBOOKS = AQ_update_codebooks( X, CODES, CODEBOOKS );
    qerror    = get_qerror( X, CODES, CODEBOOKS );
    if V, fprintf('done in %.2f seconds. Qerror is %e\n', toc, qerror); end
    
    DISTORTIONS( i+1 ) = qerror;
    
end


end

