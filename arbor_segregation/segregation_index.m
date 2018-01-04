function seg_idx = segregation_index(stat,n_cluster)

if ~exist('n_cluster','var') || isempty(n_cluster)
    n_cluster = 2;
end
[cluster, centroid] = kmeans_iter(stat,n_cluster,1000);

% Std
den = 0;
for i = 1:n_cluster
    n_elem = sum(cluster==i);
    den = den + sum((stat(cluster==i) - centroid(i)).^2)/n_elem;
end

den = (den/n_cluster)^0.5;
num = centroid(2)-centroid(1);

seg_idx = num/den; 

end