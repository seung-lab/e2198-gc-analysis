%% cell_info_update_types

function cell_info = cell_info_update_types(cell_info)
    load('cell_types');
    load('clusters');
    
    for i = 1:length(cell_types)
        type = cell_types{i};
        cell_list = clusters(type);
        for cell = cell_list'
            cell_info([cell_info.cell_id]==cell).type = type;
        end
    end
end