%% cell_info_get_somasize
% Update soma size data on cell_info 

function cell_info = cell_info_get_somasize(cell_info)
    load('somasize.mat');
    load('gc_list.mat');
    
    for cell = gc_list
        if isempty(somasize(cell))
            cell_info([cell_info.cell_id]==cell).soma_size = 0;
        else
            cell_info([cell_info.cell_id]==cell).soma_size = somasize(cell)*165*165*230/10^9;
        end
    end
 
end