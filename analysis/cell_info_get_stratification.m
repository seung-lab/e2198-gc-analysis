%% cell_info_get_stratification
% Update strat_nrml and strat_unrml fields of cell_info

function [cell_info, cell_info_off, cell_info_on, cell_info_sus, cell_info_susoff,...
    cell_info_suson, cell_info_tran, cell_info_tranoff, cell_info_tranon] = cell_info_get_stratification(cell_info)

load('cell_list.mat');

load('strat_all.mat');
skel_strat = strat;

center = 47;
onsac  = 62;
offsac = 28;

cell_info_off = cell_info;
cell_info_on = cell_info;
cell_info_sus = cell_info;
cell_info_susoff = cell_info;
cell_info_suson = cell_info;
cell_info_tran = cell_info;
cell_info_tranoff = cell_info;
cell_info_tranon = cell_info;

for cell = cell_list
    if isempty(skel_strat{cell})
        cell_info([cell_info.cell_id]==cell).strat_nrml = zeros(723,2);
        cell_info([cell_info.cell_id]==cell).strat_unrml = zeros(723,2);
    else
        depth = skel_strat{cell}(:,1);
        raw = skel_strat{cell}(:,2);
        normalized = skel_strat{cell}(:,3);
        
        cell_info([cell_info.cell_id]==cell).strat_unrml(:,1) = depth;
        cell_info([cell_info.cell_id]==cell).strat_unrml(:,2) = raw;
        cell_info([cell_info.cell_id]==cell).strat_nrml(:,1) = depth;
        cell_info([cell_info.cell_id]==cell).strat_nrml(:,2) = normalized;
        
        raw_region = zeros(723,1);
        raw_region(depth<center) = raw(depth<center);
        cell_info_off([cell_info.cell_id]==cell).strat_unrml(:,1) = depth;
        cell_info_off([cell_info.cell_id]==cell).strat_unrml(:,2) = raw_region;
        cell_info_off([cell_info.cell_id]==cell).strat_nrml(:,1) = depth;
        cell_info_off([cell_info.cell_id]==cell).strat_nrml(:,2) = raw_region/sum(raw_region);
        
        raw_region = zeros(723,1);
        raw_region(depth>center) = raw(depth>center);
        cell_info_on([cell_info.cell_id]==cell).strat_unrml(:,1) = depth;
        cell_info_on([cell_info.cell_id]==cell).strat_unrml(:,2) = raw_region;
        cell_info_on([cell_info.cell_id]==cell).strat_nrml(:,1) = depth;
        cell_info_on([cell_info.cell_id]==cell).strat_nrml(:,2) = raw_region/sum(raw_region);
        
        raw_region = zeros(723,1);
        raw_region(depth<offsac | depth>onsac) = raw(depth<offsac | depth>onsac);
        cell_info_sus([cell_info.cell_id]==cell).strat_unrml(:,1) = depth;
        cell_info_sus([cell_info.cell_id]==cell).strat_unrml(:,2) = raw_region;
        cell_info_sus([cell_info.cell_id]==cell).strat_nrml(:,1) = depth;
        cell_info_sus([cell_info.cell_id]==cell).strat_nrml(:,2) = raw_region/sum(raw_region);
        
        raw_region = zeros(723,1);
        raw_region(depth<offsac) = raw(depth<offsac);
        cell_info_susoff([cell_info.cell_id]==cell).strat_unrml(:,1) = depth;
        cell_info_susoff([cell_info.cell_id]==cell).strat_unrml(:,2) = raw_region;
        cell_info_susoff([cell_info.cell_id]==cell).strat_nrml(:,1) = depth;
        cell_info_susoff([cell_info.cell_id]==cell).strat_nrml(:,2) = raw_region/sum(raw_region);
        
        raw_region = zeros(723,1);
        raw_region(depth>onsac) = raw(depth>onsac);
        cell_info_suson([cell_info.cell_id]==cell).strat_unrml(:,1) = depth;
        cell_info_suson([cell_info.cell_id]==cell).strat_unrml(:,2) = raw_region;
        cell_info_suson([cell_info.cell_id]==cell).strat_nrml(:,1) = depth;
        cell_info_suson([cell_info.cell_id]==cell).strat_nrml(:,2) = raw_region/sum(raw_region);
        
        raw_region = zeros(723,1);
        raw_region(depth>offsac & depth<onsac) = raw(depth>offsac & depth<onsac);
        cell_info_tran([cell_info.cell_id]==cell).strat_unrml(:,1) = depth;
        cell_info_tran([cell_info.cell_id]==cell).strat_unrml(:,2) = raw_region;
        cell_info_tran([cell_info.cell_id]==cell).strat_nrml(:,1) = depth;
        cell_info_tran([cell_info.cell_id]==cell).strat_nrml(:,2) = raw_region/sum(raw_region);
        
        raw_region = zeros(723,1);
        raw_region(depth>offsac & depth<center) = raw(depth>offsac & depth<center);
        cell_info_tranoff([cell_info.cell_id]==cell).strat_unrml(:,1) = depth;
        cell_info_tranoff([cell_info.cell_id]==cell).strat_unrml(:,2) = raw_region;
        cell_info_tranoff([cell_info.cell_id]==cell).strat_nrml(:,1) = depth;
        cell_info_tranoff([cell_info.cell_id]==cell).strat_nrml(:,2) = raw_region/sum(raw_region);
        
        raw_region = zeros(723,1);
        raw_region(depth>center & depth<onsac) = raw(depth>center & depth<onsac);
        cell_info_tranon([cell_info.cell_id]==cell).strat_unrml(:,1) = depth;
        cell_info_tranon([cell_info.cell_id]==cell).strat_unrml(:,2) = raw_region;
        cell_info_tranon([cell_info.cell_id]==cell).strat_nrml(:,1) = depth;
        cell_info_tranon([cell_info.cell_id]==cell).strat_nrml(:,2) = raw_region/sum(raw_region);
        
    end
    
end

end