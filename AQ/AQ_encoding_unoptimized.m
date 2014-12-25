function [ CODES, DISTORTIONS ] = AQ_encoding_unoptimized( X, CODEBOOKS, N )
%AQ_ENCODING Encoding for additive quantization proposed by Babenko and 
% Lempitsky based on beam search.
%
% This is a version written with clarity in mind, just to make sure that
% the optimized code produces the same output. Do no use this in
% experiments.
 
% --
% Julieta

[dims, n] = size( X );

M = numel( CODEBOOKS );

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
for i = 1:numel( CODEBOOKS ),
    cb_sz = size( CODEBOOKS{i},2 );
    cbmatidx = [cbmatidx, i(ones(1, cb_sz))];
    
    full2ordered = [full2ordered, 1:cb_sz];
end

% 1:total_number_of_codewords.
full_codes = 1:numel( cbmatidx );

% Precompute dot products of all the codebooks.
prec_cbmat = sum(cbmat.*cbmat,1);

% Loop through data to encode...
for i = 1:n,
    
    % The vector to encode.
    q = X(:, i);
    
    % Approximations at iteration m.
    candidates = zeros(dims, N*N);
    
    % Codebooks used at iteration M.
    taboo_list = zeros(M, N*N);
    % Encodings at iteration M (in the original matrix).
    best_encoding = zeros(M, N*N);
    
    % Find N seeds.
    d = bsxfun(@plus, -2*q'*cbmat, prec_cbmat);
    [~, seeds_idx] = sort( d );
    
    seeds    = cbmat(:,  seeds_idx(1:N));
    
    % Update taboo list.
    for j = 1:N,
        taboo_list   ( 1, (j-1)*N+1:j*N ) = cbmatidx( seeds_idx(j) );
        best_encoding( 1, (j-1)*N+1:j*N ) = full_codes( seeds_idx(j) );
    end
    
    % Start iterations.
    for j = 2:M,
        
        % Loop through the seeds.
        for k=1:N,
            
            % The seed.
            s = seeds(:, k);
            % Codebooks already used by this seed.
            cb_idx = taboo_list( 1:j-1, (k-1)*N+1 );
            
            % Remove the codebooks that the seed is using.
            
            %ismidx = ~ismember(cbmatidx, cb_idx);
            sorted_cb_idx =  sort(cb_idx);
            ismidx = ~ismembc(cbmatidx, sorted_cb_idx);
            
            scbmat      = cbmat( :,   ismidx );
            scbmatidx   = cbmatidx(   ismidx );
            sfull_codes = full_codes( ismidx );
            
            % Find N best continuations for this seed.
            % d = bsxfun(@plus, -2*(q-s)'*scbmat, sum(scbmat.*scbmat,1));
            d = bsxfun(@plus, -2*(q-s)'*scbmat, prec_cbmat( :, ismidx ));
            [~, seeds_idx] = sort( d );
            nbest = scbmat(:,  seeds_idx(1:N));
            
            % Update the candidate list.
            candidates( :, (k-1)*N+1 : k*N ) = bsxfun(@plus, nbest, s );
            
            % And remember the codebooks where they came from, as well as
            % their encodings.
            taboo_list(    j, (k-1)*N+1 : k*N ) = scbmatidx( seeds_idx(1:N) );
            best_encoding( j, (k-1)*N+1 : k*N ) = sfull_codes( seeds_idx(1:N) );
            
            
        end
        
        % Find the closest elements from all the candidates.
        d = bsxfun(@plus, -2*q'*candidates, sum(candidates.*candidates,1));
        [didx, seeds_idx] = sort( d );
        seeds = candidates( :, seeds_idx(1:N) );
        
        seeds_codebooks = taboo_list(1:j, seeds_idx(1:N) );
        seeds_encoding  = best_encoding(1:j, seeds_idx(1:N) );
        
        % Update the taboo list.
        for k = 1:N,
            will_settle = seeds_codebooks( :, k );
            taboo_list( 1:j, (k-1)*N+1:k*N ) = will_settle(:, ones(1,N));
            
            will_settle = seeds_encoding( :, k );
            best_encoding( 1:j, (k-1)*N+1:k*N ) = will_settle(:, ones(1,N));
        end
        
    end
    
    % The best encoding found during the last iteration.
    CODES(:,i) = full2ordered( sort(best_encoding(:,1)) );
    % The distortion of such encoding.
    if nargout > 1,
        DISTORTIONS(i) = didx(1) + q'*q;
    end
    
    if ~mod(i,1000),
        fprintf('Done with %d / %d; %.2f%%.\n', i, n, 100*i/n);
    end
end