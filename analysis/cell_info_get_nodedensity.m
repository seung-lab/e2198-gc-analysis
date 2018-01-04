%% cell_info_get_nodedensity
% Update arbor density (path length per hull area (um^-1)) on cell_info

function cell_info = cell_info_get_nodedensity(cell_info)
    load('cell_hull.mat');
    load('numnodes.mat');
    load('gc_list.mat')
    
    for cell = gc_list
        if isempty(numnodes(cell))
            cell_info([cell_info.cell_id]==cell).node_density = 0;
        else
            hull_cell = cell_hull{cell};
            area = polyarea(hull_cell(:,1),hull_cell(:,2))*0.066*0.066;
            cell_info([cell_info.cell_id]==cell).node_density = numnodes(cell)*0.0957/area;
        end
    end
 
end