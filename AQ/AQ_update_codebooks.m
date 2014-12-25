function [ NEW_CODEBOOKS ] = AQ_update_codebooks( X, CODES, CODEBOOKS )
%AQ_UPDATE_CODEBOOKS Update codebooks by solving an over-constrained
%least squares problem, as proposed by Babenko and Lempitsky.
%
%   NEW_CODEBOOKS = AQ_UPDATE_CODEBOOKS(X, CODES, CODEBOOKS);
%
%   Input
%       X         : d-by-n matrix. The original vectors.
%       CODES     : nlevels-by-n matrix. X encoded.
%       CODEBOOKS : nlevels-long cell array. Entry i contains a d-by-k
%                   matrix, where each column is a codeword.
%   Output
%       NEW_CODEBOOKS : nlevels-long cell array. The updated codewords.

% --
% Julieta

d = size( X, 1 );
[ncodebooks, n] = size( CODES );

% First, we need to convert the codes to a sparse matrix. The following
% steps can surely be optimized to avoid the creation of the full matrix in
% the first place. However, the code as is is very readable and should not
% be a huge problem if you have a fair amount of memory.

C = cell( size(CODEBOOKS) );
for i = 1:ncodebooks,

    % Get the number of words in this codebook.
    [~, nwords] = size( CODEBOOKS{i} );
    
    % Create the sparse matrix per codebook (though it is defined as full now).
    spmatrix    = zeros( nwords, n );

    % Convert the codes to linear indices and set to code values to 1s in
    % the sparse matrix.
    linidx             = sub2ind( [nwords, n], int32(CODES(i,:)), int32(1:n) );
    spmatrix( linidx ) = ones(1,n);

    % Add to the matrix of all codebooks.
    C{i} = spmatrix;
end

% Concatenate all the code matrices and convert to (actually) sparse.
C  = cell2mat( C );
CC = sparse(C');

% Solve d independent problems separately.
K = zeros(d, size(C,1));

for i = 1:size(X, 1 ),
    [K(i,:), ~] = lsqr( CC, double(X(i,:)') );
end

% Transform back to cell format.
NEW_CODEBOOKS = cell( size(CODEBOOKS ));

new_codebooks_idx = 1;
for i = 1:ncodebooks
    [~, nwords] = size( CODEBOOKS{i} );

    last_idx = new_codebooks_idx + nwords - 1;

    NEW_CODEBOOKS{ i } = K( :, new_codebooks_idx:last_idx );

    new_codebooks_idx = last_idx + 1;
end

end
