%% cell_info_update

load cell_info;

cell_info = cell_info_update_types(cell_info);
cell_info_off = cell_info_update_types(cell_info_off);
cell_info_on = cell_info_update_types(cell_info_on);
cell_info_sus = cell_info_update_types(cell_info_sus);
cell_info_susoff = cell_info_update_types(cell_info_susoff);
cell_info_suson = cell_info_update_types(cell_info_suson);
cell_info_tran = cell_info_update_types(cell_info_tran);
cell_info_tranoff = cell_info_update_types(cell_info_tranoff);
cell_info_tranon = cell_info_update_types(cell_info_tranon);

cell_info = cell_info_update_avg(cell_info);
save('cell_info','cell_info','cell_info_off','cell_info_on','cell_info_sus','cell_info_susoff','cell_info_suson','cell_info_tran','cell_info_tranoff','cell_info_tranon');

get_ca_avg;
get_tempresp_avg;

