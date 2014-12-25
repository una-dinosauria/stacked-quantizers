function [ NEW_CODEBOOKS ] = AQPQ_update_codebooks( X, CODES, AQ_CODEBOOKS, SUBAQ )
%AQPQ_UPDATE_CODEBOOKS Update codebooks via least-squares in each PQ
% subspace. Return the updated codebooks.
%
% NEW_CODEBOOKS = AQPQ_update_codebooks( X, CODES, APQ_CODEBOOKS, SUBAQ, KEEP )
%
% Input
%   X              : d-by-n matrix. Each column is a datapoint.
%   CODES          : numel( AQ_CODEBOOOKS )-by-n matrix. X encoded.
%   AQ_CODEBOOOKS  : M-long cell array. Each entry is a d-by-k matrix that
%                    represents a d-dimensional codebook for AQ.
%   SUBAQ          : Integer. The number of subcodebooks to be
%                    independently quantized by AQ when using Additive 
%                    Product Quantization.
%
% Output
%   NEW_CODEBOOKS  : M-long cell arry. Follows the format of AQ_CODEBOOKS,
%                    but the codebooks have been updated via least-squares.

% --
% Julieta

% Convert to a suitable dimensionality for AQ processing.
[APQ_CODEBOOKS, KEEP] = AQ2APQcodebooks( AQ_CODEBOOKS );

% Number of subgroups to work on.
nsub = numel( KEEP );

assert( nsub == numel(APQ_CODEBOOKS)/SUBAQ, ...
    'The number of passed keep indices does not match the number of codebooks.')

NEW_CODEBOOKS = cell(size( APQ_CODEBOOKS ));

for i = 1:nsub,
    
    cbidx = (i-1)*SUBAQ+1 : i*SUBAQ;
    
    % Single out the codebooks of this subgroup.
    subcodebooks = APQ_CODEBOOKS( cbidx );
    
    % Single out the codes.
    subcodes = CODES( cbidx, : );
    
    % Single out the data to encode.
    subX     = X( KEEP{i}, : );
    
    % Solve this with regular AQ codebook update.
    NEW_CODEBOOKS( cbidx ) = AQ_update_codebooks( subX, subcodes, subcodebooks );
    
end

% Convert to full-dimensional codebooks again.
NEW_CODEBOOKS = APQ2AQcodebooks( NEW_CODEBOOKS, KEEP, SUBAQ );

end

