%% cell_types_order

function types_ordered = cell_types_order(types)
    load cell_types;
    
    idx_list = [];
    for i = 1:length(types)
        idx_list(i) = find(strcmp(cell_types,types{i}));
    end
    
    idx_list_sort = sort(idx_list);
    
    for i = 1:length(idx_list_sort)
        types_ordered{i} = cell_types{idx_list_sort(i)};
    end
    
end