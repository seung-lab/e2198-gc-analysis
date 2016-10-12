%% cell_info_get_branch
% Update arbor complexity (branch point density) on cell_info 

function cell_info = cell_info_get_branch(cell_info)
    load('numbranch.mat');
%     load('numnodes.mat');
    load('gc_list.mat');
    load('path_length'); 
    
    for cell = gc_list
        if isempty(numbranch(cell))
            cell_info([cell_info.cell_id]==cell).branch = 0;
        else
            cell_info([cell_info.cell_id]==cell).branch = numbranch(cell)/path_length(cell);
        end
    end
 
end