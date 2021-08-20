function centroids = kMeansInitCentroids(X, K)

% You should return this values correctly
centroids = zeros(K, size(X, 2));

% ====================== YOUR CODE HERE ======================
% Instructions: You should set centroids to randomly chosen examples from
%               the dataset X
indices = randperm(size(X,1));
Kindices = indices(1:K);
centroids = X(Kindices,:);







% =============================================================

end

