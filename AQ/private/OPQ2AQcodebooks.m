function [ AQ_CODEBOOKS ] = OPQ2AQcodebooks( OPQ_CODEBOOKS )
%OPQ2AQCODEBOOKS Convert codebooks from OPQ to codebooks for AQ.
%
% AQ_CODEBOOKS = AQ2APQCODEBOOKS( X, OPQ_CODEBOOKS )
%
% Input
%   OPQ_CODEBOOKS : M-long cell array. Each entry is a k-by-(d/M) matrix
%                   that represents a codebooks from OPQ.
%
% Output
%   AQ_CODEBOOKS : M-long cell array. Each entry is a d-by-k matrix that
%                  represents a d-dimensional codebook for AQ.

% --
% Julieta

% Get the total number of codebooks.
M     = numel( OPQ_CODEBOOKS );
% And the total dimensionality.
d     = sum(cellfun(@(X) size(X,1), OPQ_CODEBOOKS ));

% Create space for the output.
AQ_CODEBOOKS = cell( size(OPQ_CODEBOOKS) );

for i = 1:numel( OPQ_CODEBOOKS ),
    
    % Single out the ith codebook.
    cbi  = OPQ_CODEBOOKS{i};
    
    % Get the number of elements in it.
    [~, ncbi] = size( cbi );
    
    % Create a full-dimensional matrix to replace this codebook.
    full_codebook = zeros(d, ncbi);
    
    % Single out the dimensions where this codebook should be embedded.
    dims_idx = (i-1)*(d/M)+1:i*(d/M);
    
    % Embed the codebook.
    full_codebook( dims_idx , : ) = cbi;
    
    % And save it for later.
    AQ_CODEBOOKS{i} = full_codebook;
    
end

end

