function [centroids, idx] = runkMeans(X, initial_centroids, ...
                                      max_iters)

% Initialize values
[m n] = size(X);
K = size(initial_centroids, 1);
centroids = initial_centroids;
previous_centroids = centroids;
idx = zeros(m, 1);

% Run K-Means
for i=1:max_iters
    idx = findClosestCentroids(X, centroids);
    % Given the memberships, compute new centroids
    centroids = computeCentroids(X, idx, K);    
end
printf('\nRunning K means iteration\n')
end