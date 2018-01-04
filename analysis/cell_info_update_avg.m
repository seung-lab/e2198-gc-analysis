%% cell_info_update_avg()

function cell_info = cell_info_update_avg(cell_info)

load cell_types;

% 99901 ~ 99947
for t = 1:length(cell_types)
    cell_idx = find(strcmp({cell_info.type},cell_types{t}));
    
    strat = zeros(723,1);
    branch = 0;
    nodearea = 0;
    soma = 0;
    diameter = 0;
    
    for i = 1:length(cell_idx)
        strat = strat + cell_info(cell_idx(i)).strat_nrml(:,2);
        branch = branch + cell_info(cell_idx(i)).branch;
        nodearea = nodearea + cell_info(cell_idx(i)).node_density;
        soma = soma + cell_info(cell_idx(i)).soma_size;
        diameter = diameter + cell_info(cell_idx(i)).max_diameter;
    end
    depth = cell_info(cell_idx(i)).strat_nrml(:,1);
    strat = strat/length(cell_idx);
    branch = branch/length(cell_idx);
    nodearea = nodearea/length(cell_idx);
    soma = soma/length(cell_idx);
    diameter = diameter/length(cell_idx);
    
    cell_info(1742+t).cell_id = 99900 + t;
    cell_info(1742+t).class = 99;
    cell_info(1742+t).soma = 99;
    cell_info(1742+t).zone = 99;
    cell_info(1742+t).ctn = 99;
    cell_info(1742+t).complete = 99;
    cell_info(1742+t).type = strcat(cell_types{t},'avg');
    cell_info(1742+t).annotation = 'average';
    cell_info(1742+t).strat_nrml(:,1) = depth;
    cell_info(1742+t).strat_nrml(:,2) = strat;
    cell_info(1742+t).soma_size = soma;
    cell_info(1742+t).branch = branch;
    cell_info(1742+t).node_density = nodearea;
    cell_info(1742+t).max_diameter = diameter;
    
end