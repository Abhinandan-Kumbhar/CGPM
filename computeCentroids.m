function centroids = computeCentroids(X, idx, K)

% Useful variables
[m n] = size(X);

centroids = zeros(K, n);
for i = 1:K
  g = find(idx==i);
  centroids(i,:) = 1/size(g,1).*sum(X(g,:));
endfor

% =============================================================
end