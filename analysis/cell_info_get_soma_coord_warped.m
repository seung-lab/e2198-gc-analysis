%% cell_info_get_soma_coord_warped

function cell_info = cell_info_get_soma_coord_warped(cell_info)
    gc = cell_info_typedef_gc();

    for i = 1:length(gc)
        for c = gc(i).cells
            cell_info([cell_info.cell_id]==c).soma_coord_warped = get_m2_warped_omni_coords(cell_info([cell_info.cell_id]==c).soma_coord);
        end
    end
end