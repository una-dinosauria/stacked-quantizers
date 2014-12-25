function [ CODES, CODEBOOOKS ] = SQ_pipeline( X, NWORDS, NLEVELS, NITS, V )
%SQ_PIPELINE Training pipeline for Stacked Quantizers.
% 
% [ CODES, CODEBOOOKS ] = SQ_PIPELINE( X, NWORDS, NLEVELS, NITS, V )
%
% Input
%   X       : d-by-n matrix. Contains n d-dimensional points to use for
%               training.
%   NWORDS  : Integer. Number of words to use in each subcodebook.
%   NLEVELS : Integer. Number of subcodebooks to use.
%   NITS    : Integer. Number of iterations for optimization after
%               initialization.
%   V       : Boolean. Whether to print progress.
%
% Output
%   CODES     : NLEVELS-by-n matrix. The encoding of X.
%   CODEBOOKS : NLEVELS-long cell array. Each entry has a d-by-NWORDS
%                 matrix that represents a codebook.

% --
% Julieta

if nargin < 4,
    V = false;
end

tic;
%% Initialize with k-means.
[ CODES, CODEBOOOKS ] = SQ_initialize_kmeans( X, NWORDS, NLEVELS, 100, V );

qerror = get_qerror( X, CODES, CODEBOOOKS );
if V, fprintf('Error after initialization  is %e. Done in %.2f seconds.\n', qerror, toc); end

%% Iterate between encoding and updating codebooks.
for i = 1:NITS,
    
    if V, fprintf('=== Iteration %d / %d. ===\n', i, NITS); tic; end
    
    % Dummy singletons.
    singletons = cell(2,1);
    
    % Xr keeps what we have to encode.
    Xr = X;
    Xd = X - SQ_decode( CODES(2:end,:), CODEBOOOKS(2:end) );
    
    for j = 1:NLEVELS,
    
        if V, fprintf('Updating codebook %02d... ', j); tic; end
        
        % Compute Xd as the reconstruction with all the codebooks *except*
        % the jth one, which we are going to update in this iteration.
        if j == NLEVELS,
            Xd = Xr - SQ_decode( CODES(j-1,:), CODEBOOOKS(j-1) );
        elseif j > 1,
            Xd = Xr - SQ_decode( CODES([j-1, j+1:end],:), CODEBOOOKS([j-1, j+1:end]) );
        end

        codes2update = CODES( j, :);
        
        % Update codebook using the implementation of @norouzi.
        k = size( CODEBOOOKS{j}, 2 );
        new_codebook = kmeans_iter_mex( Xd, int32( codes2update' ), k );

        % Add singletons to the codebook.
        aa = true(k, 1);
        aa( codes2update ) = false;
        new_codebook( :, aa ) = singletons{2};
        
        CODEBOOOKS{j} = new_codebook;
        
        if V, fprintf('done in %.2f seconds.\n', toc); end
        % Update the residual.
        if j > 1,
            Xr = Xr - CODEBOOOKS{j-1}(:, CODES(j-1,:));
        end
        
        if V, fprintf('Updating codes... '); tic; end
        [CODES(j:end,:), singletons] = SQ_encode( Xr, CODEBOOOKS(j:end), V > 1 );
        if V, fprintf('done in %.2f seconds.\n', toc); end
    
    end
    
    if V,
        qerror = get_qerror( X, CODES, CODEBOOOKS );
        fprintf('Iteration %d / %d done in %.2f seconds. Qerror is %e.\n', i, NITS, toc, qerror);
    end
    
end

end