function recall_vs_sel = eval_recall_vs_sel(index, num_test, gt, K, selectivity)
%RECALL_VS_SEL Make recal vs. selectivity curves.

% Modified from Kaiming He's K-means hashing code.

index = index(:,1:selectivity);

% find the sucess items
recall_vs_sel = zeros(num_test, selectivity);

for i=1:num_test
    index_one_query = index(i, :);
    gt_one_query = gt(i, 1:K);
    
    % find correct items
    hit = (sqdist(gt_one_query, index_one_query) == 0);
    recall_vs_sel(i, :) = cumsum(sum(hit, 1)) / K;
end
recall_vs_sel = mean(recall_vs_sel, 1);

end
