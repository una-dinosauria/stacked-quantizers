function [ CODES, DISTORTIONS ] = AQ_encoding( X, CODEBOOKS, N, V )
%AQ_ENCODING Encoding for additive quantization proposed by Babenko and 
% Lempitsky based on beam search. This is a *highly* optimized
% implementation, with lots of precomputation and vectorization in the
% loop.
%
% [CODES, DISTORTIONS] = AQ_ENCODING( X, CODEBOOKS, N, V )
%
% Input
%   X         : d-by-n matrix. Each column is a datapoint to encode.
%   CODEBOOKS : n-long cell array. Each entry is a d-by-k codebook.
%   N         : Integer. Beam search depth.
%   V         : Boolean. Whether to print progrees. Defaults to false.
%
% Output
%   CODES       : numel(codebooks)-by-n matrix. X encoded.
%   DISTORTIONS : n-long vector. The squared distortion of each point.

% --
% Julieta

if nargin < 4,
    V = false;
end

[~, n] = size( X );

M   = numel( CODEBOOKS );
assert( mod(M,2) == 0, 'The number of codebooks must be a multiple of 4' );

% If we have more than 4 codebooks, run AQ/PQ encoding.
if M > 4,
    [ CODES, DISTORTIONS ] = AQPQ_encoding( X, CODEBOOKS, N, V );
    return;
end


% Assuming that all the codebooks have the same size.
cbsz = size( CODEBOOKS{1}, 2);

% Create output variables.
CODES       = zeros(M, n);
if nargout > 1,
    DISTORTIONS = zeros(n,1);
end

% Convert the codebooks to a single matrix.
cbmat = cell2mat( CODEBOOKS' );

% Create an index from the matrix to the codebook number.
cbmatidx = [];
full2ordered = [];

for i = 1:M,
    cb_sz = size( CODEBOOKS{i},2 );
    cbmatidx = [cbmatidx, i(ones(1, cb_sz))];
    
    full2ordered = [full2ordered, 1:cb_sz];
end

% 1:total_number_of_codewords.
full_codes = 1:numel( cbmatidx );

% Precompute lookup table.
booltable    = true( M, numel( cbmatidx ) );
btidx = 1;
for i = 1:M,
    cb_sz = size( CODEBOOKS{i},2 );
    booltable(i, btidx:btidx+cb_sz-1) = false;
    btidx = btidx + cb_sz;
end

% Precompute dot products of all the codebooks.
prec_cbmat = sum(cbmat.*cbmat,1);

% Seeds expansion index: [1, 1, 1, 1, ....., N, N, N, N, N]
Sexpidx = [];
for i = 1:N,
    Sexpidx = [Sexpidx, i(ones(1,N))];
end


% Loop through data to encode...
% for i = 1:n,
parfor i = 1:n,
    
    % The vector to encode.
    q = X(:, i);
    
    % Codebooks used at iteration M.
    taboo_list = zeros(M, N*N);
    % Encodings at iteration M (in the original matrix).
    best_encoding = zeros(M, N*N);
    
    % Find N seeds.
    d = bsxfun(@plus, -2*q'*cbmat, prec_cbmat);
    [~, seeds_idx] = sort( d );
    
    seeds    = cbmat(:,  seeds_idx(1:N));
    
    % Update taboo and encoding lists.
    taboo_list( 1, : )    = cbmatidx(  seeds_idx( Sexpidx ));
    best_encoding( 1, : ) = full_codes( seeds_idx( Sexpidx ));
    
    % Start iterations.
    for j = 2:M,
        
        % Update the query.
        res = bsxfun(@minus, q, seeds );
        
        % Find N best continuations for this seed.
        d = bsxfun(@plus, -2*(res)'*cbmat, prec_cbmat);
        
        seeds_idx      = zeros( N, (M-j+1)*cbsz, 'single' );
        allscbmatidx   = zeros( N, (M-j+1)*cbsz, 'single' );
        allsfull_codes = zeros( N, (M-j+1)*cbsz, 'single' );
        
        cb_idxes = taboo_list(1:j-1, 1:N:end );
        
        % Loop through the seeds.
        for k = 1:N,
            
            % Codebooks already used by this seed.
            cb_idx = cb_idxes(:, k );
            
            % Remove the codebooks that the seed is using. This is actually
            % faster than a containers.Map.
            ismidx = all(booltable( cb_idx , : ), 1);
            
            allscbmatidx(k,:)   = cbmatidx(   ismidx );
            allsfull_codes(k,:) = full_codes( ismidx );
            
            seeds_idx(k,:) = d(k, ismidx);
            
        end
        
        [~, seeds_idx] = sort( seeds_idx, 2 );
        
        seeds_idx = seeds_idx';
        seeds_idx = bsxfun(@plus, seeds_idx, size(allscbmatidx,2)*(0:size(allscbmatidx,1)-1) );
        
        allscbmatidx   = allscbmatidx';
        allscbmatidx   = allscbmatidx( seeds_idx );
        
        allsfull_codes = allsfull_codes';
        allsfull_codes = allsfull_codes( seeds_idx );
        
        allscbmatidx   = allscbmatidx(1:N, :);
        allsfull_codes = allsfull_codes(1:N, :);
        
        nbest = cbmat( :, allsfull_codes(:) );

        % Update the candidate list.
        candidates = nbest + seeds(:, Sexpidx);

        % And remember the codebooks where the candidates came from, as 
        % well as their encodings.
        taboo_list(    j, : ) = allscbmatidx( : );
        best_encoding( j, : ) = allsfull_codes( : );
        
        % Find the closest elements from all the candidates.
        d = bsxfun(@plus, -2*q'*candidates, sum(candidates.*candidates,1));
        [didx, seeds_idx] = sort( d );
        seeds = candidates( :, seeds_idx(1:N) );
        
        seeds_codebooks = taboo_list(1:j, seeds_idx(1:N) );
        seeds_encoding  = best_encoding(1:j, seeds_idx(1:N) );
        
        % Update the taboo and encoding lists.
        taboo_list   ( 1:j, :) = seeds_codebooks( :, Sexpidx );
        best_encoding( 1:j, :) =  seeds_encoding( :, Sexpidx );
        
    end
    
    % The best encoding found during the last iteration.
    CODES(:,i) = full2ordered( sort(best_encoding(:,1)) );
    % The distortion of such encoding.
    if nargout > 1,
        DISTORTIONS(i) = didx(1) + q'*q;
    end
    
    if V && ~mod(i,10000),
        fprintf('Done with %d / %d; %.2f%%.\n', i, n, 100*i/n);
        %toc;
        %tic;
    end
end


end