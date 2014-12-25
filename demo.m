% Approximate nearest neighbour search demo.

% --
% Julieta, 2014

clc
clear

% Load training, query and base datasets, and the query ground truth.
X = load( 'data/features_m_128.mat', 'feats_m_128_train', 'feats_m_128_test', 'feats_m_128_base');
X_train = X.feats_m_128_train;
X_test  = X.feats_m_128_test;
X_base  = X.feats_m_128_base;
nquery  = size( X_test, 2 ); % Number of query vectors.

gt      = load( 'data/features_m_128_gt.mat', 'gt' ); gt = gt.gt;
K       = 1; % Number of nearest neighbours to search for.
gt      = gt(1:nquery, 1:K)';

%% Set train and search parameters
m       = 4;   % Number of subcodebooks.
h       = 256; % Number of cluster centres per subcodebook.
nbits   = log2(h) * m;  % Number of bits in the final code.

% === The following values are here for a quick demo. 
nitsOPQ = 10; % Iterations in PQ / OPQ. (100 in paper)
nitsSQ  = 10; % Iterations in SQ.       (100 in paper)
nitsAQ  = 10; % Iterations in AQ.       ( 10 in paper)

verbose = 1; % Print progress.

N_train = 16; % Pool size for beam search in AQ during training.
N_base  = 64; % Pool size for beam search in AQ during database encoding.

selectivity = 10000; % Number of nearest-neighbours to retrieve.

%% === PQ (no preprocessing) ===
fprintf('=== PQ: %d codebooks. ===\n', m);

% Train
[model, ~] = product_quantization( X_train, m, h, nitsOPQ );

% Quantize the database
cbase = uint8( quantize_by_ckmeans(X_base, model) -1 );

% Search
centers = double(cat(3, model.centers{:}));
npoints = size(cbase, 2);

fprintf('Searching... '); tic;
queryR       = double( model.R' * X_test );
[ids_aqd, ~] = linscan_aqd_knn_mex( cbase, queryR, npoints, nbits, selectivity, centers);
fprintf('done in %.2f seconds\n', toc);

% Plot recall@N curve
recall_at_k_aqd_pq = eval_recall_vs_sel( double(ids_aqd'), nquery, double(gt'), K, selectivity );
semilogx( recall_at_k_aqd_pq, 'b-', 'linewidth', 2 ); 
grid on; hold on; xlabel('N'); ylabel('Recall@N');
legend('PQ', 'location', 'northwest');
pause(0.5);

%% === OPQ ===
fprintf('=== OPQ: %d codebooks. ===\n', m);

% Train
[model, ~] = ckmeans(X_train, m, h, nitsOPQ, 'natural');

% Quantize the database
cbase = uint8( quantize_by_ckmeans(X_base, model) -1 );

% Search
centers = double(cat(3, model.centers{:}));
npoints = size(cbase, 2);

fprintf('Searching... '); tic;
queryR       = double( model.R' * X_test );
[ids_aqd, ~] = linscan_aqd_knn_mex( cbase, queryR, npoints, nbits, selectivity, centers);
fprintf('done in %.2f seconds\n', toc);

% Plot recall@N curve
recall_at_k_aqd_opq = eval_recall_vs_sel( double(ids_aqd'), nquery, double(gt'), K, selectivity );
semilogx( recall_at_k_aqd_opq, 'r-', 'linewidth', 2 );
legend('PQ', 'OPQ', 'location', 'northwest');
pause(0.5);

%% === SQ ===
fprintf('=== SQ ncodebooks %d ===.\n', m);

% Train
[~, codebooks] = SQ_pipeline( X_train, h, m, nitsSQ, verbose );

% Quantize the database
cbase = SQ_encode( X_base, codebooks, verbose );

% Compute database l2 norms.
dbnorms = single( sum( SQ_decode( cbase,  codebooks ).^2, 1 ) );

% Convert cbase to uint8
cbase = uint8( cbase -1 );

fprintf('Searching... '); tic;
[~, idx] = SQ_search( cbase, codebooks, X_test, dbnorms, selectivity);
fprintf('done in %.2f seconds\n', toc);

% Plot recall@N curve
recall_at_k_aqd_sq = eval_recall_vs_sel( double(idx'), nquery, double(gt'), K, 10000 );
semilogx( recall_at_k_aqd_sq, 'm-', 'linewidth', 2 );
legend('PQ', 'OPQ', 'SQ', 'location', 'northwest');
pause(0.5);

%% === AQ === 
fprintf('=== AQ ncodebooks %d ===.\n', m);

% Train
[codes, codebooks, distortions, R] = AQ_pipeline( X_train, h(ones(1,m)), nitsAQ, N_train, verbose );

% Quantize the database
fprintf('Encoding the database; this is gonna take a while... '); tic;
cbase = AQ_encoding( R'*X_base, codebooks, N_base, verbose );
fprintf('done in %.2f seconds\n', toc);

% Compute database l2 norms.
dbnorms = single( sum( SQ_decode( cbase,  codebooks ).^2, 1 ) );

% Convert cbase to uint8
cbase = uint8( cbase -1 );

fprintf('Searching... '); tic;
for i = 1:numel(codebooks), codebooks{i} = single( codebooks{i} ); end
[~, idx] = SQ_search( cbase, codebooks, X_test, dbnorms, selectivity);
fprintf('done in %.2f seconds\n', toc);

% Plot recall@N curve
recall_at_k_aqd_aq = eval_recall_vs_sel( double(idx'), nquery, double(gt'), K, 10000 );
semilogx( recall_at_k_aqd_aq, 'k-', 'linewidth', 2 );
legend('PQ', 'OPQ', 'SQ', 'AQ', 'location', 'northwest');

