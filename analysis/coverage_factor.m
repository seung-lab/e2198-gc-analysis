function cf = coverage_factor(type)
    load('./cell_hull');
    
    size_x = 5376;
    size_y = round(3456*23/16.5);

    cut_x = [1000,1000];
    cut_y = [1000,1000];
%     
%     cut_x = [0,0];
%     cut_y = [0,0];
    
    size_x_crop = size_x - sum(cut_x);
    size_y_crop = size_y - sum(cut_y);
    
    gc = cell_info_typedef_gc();
    cell_list = gc(strcmp({gc.name},type)).cells;
    n_cell = length(cell_list);
    
    area_list = zeros(n_cell,1);
    
    patch_all = zeros(size_x_crop,size_y_crop);
    for i = 1:n_cell
        cell = cell_list(i);
        hull_cell = cell_hull{cell};
        
        in_hull = poly2mask(hull_cell(:,1),hull_cell(:,2),size_y,size_x);
        patch_hull = transpose(in_hull);
        
        patch_hull = cutedge(patch_hull,cut_x,cut_y);
        
        patch_all = patch_all + patch_hull;
        area_list(i) = sum(sum(patch_hull));
    end
    
    cf = sum(area_list)/sum(patch_all(:)>=1);
    
end