%% cell_info_get_nodedensity
% Update arbor density (number of nodes per hull area) on cell_info

function cell_info = cell_info_get_nodedensity(cell_info)
    load('area_hull.mat');
    load('numnodes.mat');
    load('gc_list.mat');
    
    for cell = gc_list
        if isempty(numnodes(cell))
            cell_info([cell_info.cell_id]==cell).node_density = 0;
        else
            cell_info([cell_info.cell_id]==cell).node_density = numnodes(cell)/(area_hull{cell}*0.066*0.066);
        end
    end
 
end