%% path_len 
% Compute path length of the skeletons
% [mm]

function l_list = path_len(cell_list)
n_cell = length(cell_list);

for i = 1:n_cell
    cell = cell_list(i);
    
    load(strcat('../skeletons/skel_',num2str(cell)));
    n(:,1) = n(:,1)*23/16.5;
    
    l = 0;
    for j = 1:size(e,1)
        l = l + sum((n(e(j,2),:) - n(e(j,1),:)).^2)^0.5 * 66*10^-6;
    end
    
    l_list(i) = l;
end

end
