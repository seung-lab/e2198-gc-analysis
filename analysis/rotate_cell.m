function nodes_rotate = rotate_cell(nodes_cell,deg,soma_loc)
    rot_mat = [cos(deg),-sin(deg);sin(deg),cos(deg)];
    
    nodes_rotate = nodes_cell*rot_mat + repmat(soma_loc - soma_loc*rot_mat,[size(nodes_cell,1),1]);
    nodes_rotate = round(nodes_rotate);
end