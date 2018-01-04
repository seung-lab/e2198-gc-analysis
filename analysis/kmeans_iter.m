function [cluster, centroid] = kmeans_iter(x,k,n_iter)

cluster = zeros(size(x,1),1);
centroid = zeros(size(x,1),k);
for i = 1:n_iter
    [cluster_i,centroid_i] = kmeans(x,k);
    
    [cluster(:,i),centroid(i,:)] = sort_kmeans(cluster_i,centroid_i);
end

cluster = round(mean(cluster,2));
centroid = mean(centroid,1);

end