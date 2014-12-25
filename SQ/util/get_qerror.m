function qerror = get_qerror( X, CODES, CODEBOOKS, R)
%GET_QERROR Compute quantization error.
% 
% QERROR = GET_QERROR( X, CODES, CODEBOOKS ) Computes the quantization
% error on X, as encoded by CODES given the passed CODEBOOKS.
%
% QERROR = GET_QERROR( X, CODES, CODEBOOKS, R ) Rotates the data by R
% previous to computing the quantization error.

% --
% Julieta.

if nargin < 4,
    % No rotation :(
    qerror = mean( sum( ( X - SQ_decode(CODES, CODEBOOKS) ).^2 ) );
else
    qerror = mean( sum( ( (X'*R)' - SQ_decode(CODES, CODEBOOKS) ).^2 ) );
end
    
end