% Runs Cartesian k-means (ck-means) on the dataset X (n p-dimensional
% data points) and generates m subspaces with h centers in each. As a
% result a total of h^m centers can be represented.

% =====
% 
% Slight modification of ckmeans.m (@norouzi) by Julieta (@jltmtz).
%
% =====

function [model, B] = product_quantization( X, m, h, niter )

% Inputs:
%       X: p*n -- n p-dimensional input points used for training the quantizer.
%       m: number of subspaces in the encoding.
%       h: number of centers per subspace.
%          if h is an m-dim array, h(i) is the number of centers in
%          the i-th subspace
%       niter: number of iterations for the k-means optimization.

% Outputs:
%       model: the learned quantization model and its parameters.
%       B: the leaned codes.

% obj: the quantization error objective.
obj = Inf;

model.type = 'pq';
n = size(X, 2);
p = size(X, 1);
model.m = m;
model.p = p;

if (length(h) == 1)
  h = ones(m, 1) * h;
end
nbits = sum(log2(h));
model.h = h;
model.nbits = nbits;

% len: an m-dim array the i-th element of which, i.e. len(i), stores
%      the dimensionality of the i-th subspace.
len = ones(m, 1) * floor(p / m);
len(1:mod(p, m)) = len(1:mod(p, m)) + 1;  % p = m * floor(p / m) + mod(p, m)
model.len = len;

len0 = 1 + cumsum([0; len(1:end-1)]);
len1 = cumsum(len);

DB = zeros(size(X), 'single');  % DB stores D*B
 
% Rotation is the identity (no rotation)
R       = eye(p, p, 'single');
model.R = R;
RX      = X;

% initialize D
D = cell(m, 1);

% inializing D by random selection of subspace centers (after rotation).
for i=1:m,
  perm = randperm(n, h(i));
  D{i} = RX(len0(i):len1(i), perm);
end

% initialize B
B = zeros(n, m, 'int32');
for i=1:m,
  B(:, i) = euc_nn_mex(D{i}, RX(len0(i):len1(i), :));
  DB(len0(i):len1(i), :) = D{i}(:, B(:, i));
end

for iter=0:niter,
  if mod(iter, 10) == 0,
    objlast = obj;
    tmp = R * DB;
    tmp = tmp - X;
    tmp = tmp.^2;
    obj = mean(sum(tmp, 'double'));
    clear tmp;
    fprintf('%3d %e   \n', iter, obj);  
    model.obj(iter+1) = obj;
  else
    fprintf('%3d\r', iter);
  end
  
  if (objlast - obj < model.obj(1) * 1e-5)
    fprintf('not enough improvement in the objective... breaking.\n')
    break;
  end
  
  for i=1:m,
    % update D
    D{i} = kmeans_iter_mex(RX(len0(i):len1(i), :), B(:, i), h(i));
    
    % update B
    B(:, i) = euc_nn_mex(D{i}, RX(len0(i):len1(i), :));
    
    % update D*B
    DB(len0(i):len1(i), :) = D{i}(:, B(:, i));
  end
end

for i=1:m,
  model.centers{i} = D{i};
end

B = B';