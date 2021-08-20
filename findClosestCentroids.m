function idx = findClosestCentroids(X, centroids)

% Set K
K = size(centroids, 1);
m = size(X, 1);
% You need to return the following variables correctly.
idx = zeros(size(X,1), 1);

for i = 1: m
  for j = 1:K
    dist = (X(i,:)-centroids(j,:));
    distance(j) = (dist*dist');
  endfor  
  idx(i) = min(find(distance == min(distance)));
endfor
end