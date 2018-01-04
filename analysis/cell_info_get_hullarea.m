%% cell_info_get_hullarea
% Update convex hull area (um^2) on cell_info 

function cell_info = cell_info_get_hullarea(cell_info)
    load('cell_hull.mat'); 
    load('gc_list.mat');
    
    for cell = gc_list
        hull_cell = cell_hull{cell};
        cell_info([cell_info.cell_id]==cell).hull_area = polyarea(hull_cell(:,1),hull_cell(:,2))*0.066*0.066;
    end
 
end