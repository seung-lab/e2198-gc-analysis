function nodes_reflect = reflect_cell(nodes_cell,direction,soma_loc)
    % direction - 1: |
    %             2: ㅡ
    %             3: \
    %             4: /
    
    if direction == 1
        ref_mat = [-1,0;0,1];
        ref_bias = soma_loc*[1,0;0,0];
    
        nodes_reflect = (nodes_cell-repmat(ref_bias,[size(nodes_cell,1),1]))*ref_mat + repmat(ref_bias,[size(nodes_cell,1),1]);
        nodes_reflect = round(nodes_reflect);
        
    elseif direction == 2
        ref_mat = [1,0;0,-1];
        ref_bias = soma_loc*[0,0;0,1];
        
        nodes_reflect = (nodes_cell-repmat(ref_bias,[size(nodes_cell,1),1]))*ref_mat + repmat(ref_bias,[size(nodes_cell,1),1]);
        nodes_reflect = round(nodes_reflect);
        
    elseif direction == 3
        nodes_cell = rotate_cell(nodes_cell,pi/4,soma_loc);
        ref_mat = [-1,0;0,1];
        ref_bias = soma_loc*[1,0;0,0];
    
        nodes_reflect = (nodes_cell-repmat(ref_bias,[size(nodes_cell,1),1]))*ref_mat + repmat(ref_bias,[size(nodes_cell,1),1]);
        nodes_reflect = rotate_cell(nodes_reflect,-pi/4,soma_loc);
      
        nodes_reflect = round(nodes_reflect);
        
    elseif direction == 4
        nodes_cell = rotate_cell(nodes_cell,-pi/4,soma_loc);
        ref_mat = [-1,0;0,1];
        ref_bias = soma_loc*[1,0;0,0];
    
        nodes_reflect = (nodes_cell-repmat(ref_bias,[size(nodes_cell,1),1]))*ref_mat + repmat(ref_bias,[size(nodes_cell,1),1]);
        nodes_reflect = rotate_cell(nodes_reflect,pi/4,soma_loc);
      
        nodes_reflect = round(nodes_reflect);
        
    else
        ref_mat = [1,0;0,1];
        ref_bias = [0,0];
        
        nodes_reflect = (nodes_cell-repmat(ref_bias,[size(nodes_cell,1),1]))*ref_mat + repmat(ref_bias,[size(nodes_cell,1),1]);
        nodes_reflect = round(nodes_reflect);
    end
    
    
end