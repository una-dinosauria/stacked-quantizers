function [ APQ_CODEBOOKS, KEEP ] = AQ2APQcodebooks( AQ_CODEBOOKS, SUBAQ )
%AQ2APQCODEBOOKS Convertes AQ codebooks to APQ codebooks.
% 
% [ APQ_CODEBOOKS, KEEP ] = AQ2APQCODEBOOKS( AQ_CODEBOOKS, SUBAQ )
% 
% AQ is Additive Quantization, and APQ is Additive Product Quantization.
% The idea is that after calling this function, each subgroup of SUBAQ
% codebooks can be optimized independently by AQ.
%
% Input
%   AQ_CODEBOOKS : M-long cell array. Each entry is a d-by-k matrix that
%                  represents a d-dimensional codebook for AQ.
%   SUBAQ        : Integer. AQ will receive this many codebooks. The
%                  default is 4.
%
% Output
%   APQ_CODEBOOKS : M-long cell array. Each entry is a (d/M)*SUBAQ-by-k
%                   matrix that represents a codebook for APQ.
%   KEEP          : M/SUBAQ-long cell array. Each entry is a boolean array
%                   that contains the dimensions that are non-zero for the
%                   ith subgroup of codebooks.

% --
% Julieta

if nargin < 2,
    SUBAQ = 4;
end

% Total number of codebooks.
M = numel( AQ_CODEBOOKS );

% Number of groups to create.SUBAQ
nsub = M / SUBAQ;

APQ_CODEBOOKS = cell( size(AQ_CODEBOOKS) );
KEEP          = cell( nsub, 1 );

for i = 1:nsub,
    % Single out this group of codebooks.
    subcodebooks = AQ_CODEBOOKS( (i-1)*SUBAQ+1 : i*SUBAQ );
    
    % Find the dimensions where these codebooks are non-zero.
    KEEP{i} = sum( cell2mat( subcodebooks' ), 2 ) ~= 0;
    
    for j = 1:numel( subcodebooks ),
        % Filter out the undesired dimensions.
        subcodebooks{j} = subcodebooks{j}( KEEP{i}, : );
    end
    
    % Save the codebooks to the output.
    APQ_CODEBOOKS( (i-1)*SUBAQ+1 : i*SUBAQ ) = subcodebooks;
end

end

