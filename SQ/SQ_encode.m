function [CODES, SINGLETONS] = SQ_encode( X, CODEBOOKS, V )
%SQ_ENCODE Encode a database using stacked quantizers.
%
% [CODES, SINGLETONS] = SQ_encode( X, CODEBOOKS, V )
%
% Input
%   X         : d-by-n matrix. Each column is a datapoint to encode.
%   CODEBOOKS : nlevels-long cell array. Each entry contains an
%                d-by-nwords matrix of codewords used to encode the 
%                database.
%   V         : Boolean. Whether to print progress.
%
% Output
%   CODES: nlevels-by-n matrix. X encoded.

% --
% Julieta

% Get the number of datapoints.
[~, n] = size( X );

% Get the number of codebooks.
nlevels = numel( CODEBOOKS );

if nargin < 4,
    V = false;
end

% Make space for the output.
CODES = zeros( nlevels, n );

% Return possible singletons.
SINGLETONS = cell( nlevels, 1 );

xx = X;

% Compute the codes for each level.
for i = 1:numel( CODEBOOKS ),
    
    yy = CODEBOOKS{i};
    
    % Encode at this level.
    if V, fprintf('Encoding level %d... ', i); tic; end
    
    % Memory-efficient distance computation by @norouzi.
    yy = single(yy);
    xx = single(xx);
    [idx, mindists] = euc_nn_mex( yy, xx );
    
    if V, fprintf('Done in %.2f seconds.\n', toc); end
    
    % Are there empty clusters? If so make singletons out of the points
    % farthest from the cluster centres.
    if nargout > 1,
    
        nwords_expected = size( CODEBOOKS{i}, 2 );
        nwords_assigned = numel( unique( idx ));
        
        if nwords_assigned < nwords_expected,
            [~, sdidx] = sort( mindists );
            sdidx = gather( sdidx);
            SINGLETONS{i} = X(:, sdidx( end-(nwords_expected - nwords_assigned)+1:end ) );
        end
    
    end
    
    % Save the produced codes.
    CODES( i, : ) = idx;
    
    % Subtract the encoded part.
    xx = xx - yy(:, idx);
    
end

end

