function [ AQ_CODEBOOKS ] = APQ2AQcodebooks( APQ_CODEBOOKS, KEEP, SUBAQ )
%APQ2AQCODEBOOKS Additive Quantization to Additive Product Quantization
% codebook converter. This function gives you full-dimensional vectors.
%
% [ AQ_CODEBOOKS ] = APQ2AQCODEBOOKS( APQ_CODEBOOKS, KEEP, SUBAQ )
%
% Input
%   APQ_CODEBOOKS : M-long cell array. Each entry is a (d/M)*SUBAQ-by-k 
%                   matrix that represents a codebook for APQ.
%   KEEP          : M/SUBAQ-long cell array. Each entry is a boolean array
%                   that contains the dimensions that are non-zero for the
%                   ith subgroup of codebooks.
%   SUBAQ         : Integer. APQ was run on groups of this many codebooks.
%                   This is usually 4.
%
% Output
%   AQ_CODEBOOKS : M-long cell array. Each entry is a d-by-k matrix that 
%                   represents a codebook for APQ.

% --
% Julieta

M    = numel( APQ_CODEBOOKS );
nsub = M / SUBAQ;

% The full dimensionality is the size of the first one, times the number of
% subspaces.
d = size( APQ_CODEBOOKS{1}, 1) * nsub;

AQ_CODEBOOKS = cell( size(APQ_CODEBOOKS) );

% Convert the codebooks back to full-dimensional vectors.
for i = 1:nsub,
    
    % The indices of this subgroup.
    dimidx = (i-1)*SUBAQ+1 : i*SUBAQ;
    
    % Single out the codebooks of this subgroup.
    subpaq_codebooks = APQ_CODEBOOKS( dimidx );
    
    % Embed each of them in a full-dimensional codebook.
    for j = 1:numel( subpaq_codebooks )
        % Get the number of vectors.
        ncbi = size( subpaq_codebooks{j}, 2 );
        
        % Embed in the full-dimensional codebook.
        full_codebook = zeros(d, ncbi);
        full_codebook( KEEP{i}, :) = subpaq_codebooks{j};
        subpaq_codebooks{j} = full_codebook;
    end
    
    % Save to the output.
    AQ_CODEBOOKS( dimidx ) = subpaq_codebooks;
end

end

