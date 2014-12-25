function [ X ] = SQ_decode( CODES, CODEBOOKS )
%SQ_DECODE Reconstruct vectors encoded by stacked quantizers.
%
% X = SQ_decode( CODES, CODEBOOKS ) returns the reconstructed database X
% from the passed CODES and CODEBOOKS.
% 
% Input
%   CODES     : nlevels-by-n matrix. The encoded database.
%   CODEBOOKS : nlevels-long cell array. Each entry contains a 
%               d-by-nwords matrix of codewords used to encode the 
%               database.
% Output
%   X         : d-by-n matrix. The reconstructed vectors.

% --
% Julieta

% Reconstruct with the first codebook.
X =  CODEBOOKS{1}(:, CODES(1,:));

% Add the rest of the codebooks.
for i = 2:size(CODES,1),
     X = X + CODEBOOKS{i}(:, CODES(i,:));
end

end
