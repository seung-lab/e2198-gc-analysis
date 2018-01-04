function [density_patch, patch_node, nodes, soma_coord, winner_patch, winner_id] = density_grid(cell_info,cell_list,in_thr,out_thr,grid_size,permute,crop)
load('./cell_hull');
load('./gc_list');
load('./clusters');

visual = 0;

tic;
if isstr(cell_list)
    type = cell_list;
    gc = cell_info_typedef_gc();
    cell_list = gc(strcmp({gc.name},type)).cells;
end

if ~exist('crop','var') || isempty(crop)
    crop = 1;
end

n_cell = length(cell_list);

if length(in_thr) == 1
    in_thr = repmat(in_thr,[1,n_cell]);
end
if length(out_thr) == 1
    out_thr = repmat(out_thr,[1,n_cell]);
end


% Size of the patch
size_x = 5376;
size_y = round(3456*23/16.5);
% size_x = 5280;
% size_y = 4570;


if crop == 0
    cut_x = [0,0]; cut_y = [0,0];
else
    cut_x = [1000,1000]; cut_y = [1000,1000];
end
size_x_cut = size_x-sum(cut_x);
size_y_cut = size_y-sum(cut_y);

% Patch of skeleton nodes
patch_node = sparse(size_x,size_y);
patch_cells = sparse(size_x_cut*size_y_cut,n_cell);

if strcmp(permute,'random_rot')  
    % Random soma, rotation/reflection
    rotate_angle = permute_rotation(type,0);
    reflect_direction = permute_reflection(type,0);
    
    soma_whole = zeros(n_whole,2);
    for i = 1:n_whole
        cell = whole_list(i);
        nodes_cell = trunkate_cell(cell,in_thr(i),out_thr(i));
        nodes_whole{i} = nodes_cell;
        soma = cell_info([cell_info.cell_id]==cell).soma_coord_warped;
        soma_plane = soma(2:3);
        soma_plane(2) = soma_plane(2)*(23/16.5);
        soma_whole(i,:) = soma_plane;
    end
    
    soma_random(:,1) = randi(size_x,n_cell,1);
    soma_random(:,2) = randi(size_y,n_cell,1);
    soma_coord = soma_random;
    for i = 1:n_cell
        sample = randi(n_whole);
        
        nodes_cell = nodes_whole{sample};
        soma_cell = soma_whole(sample,:);
        nodes_cell = rotate_cell(nodes_cell,rotate_angle(i),soma_cell);
        nodes_cell = reflect_cell(nodes_cell,reflect_direction(i),soma_cell);
        
        shift = soma_random(i,:) - soma_whole(sample,:);
        nodes_cell = nodes_cell + repmat(shift,[size(nodes_cell,1),1]);
        
        nodes_cell(nodes_cell(:,1)<1 | nodes_cell(:,2)<1 | ...
            nodes_cell(:,1)>size_x | nodes_cell(:,2)>size_y,:) = [];
        
        nodes{i} = nodes_cell;
        
        % Patch of skeleton nodes
        patch_node = patch_node + sparse(nodes_cell(:,1),round(nodes_cell(:,2)),1,size_x,size_y);
    end
    
elseif strcmp(permute,'soma_rot')
     % Permute soma, rotation/reflection
    rotate_angle = permute_rotation(type,0);
    reflect_direction = permute_reflection(type,0);
    
    soma_whole = zeros(n_whole,2);
    for i = 1:n_whole
        cell = whole_list(i);
        nodes_cell = trunkate_cell(cell,in_thr(i),out_thr(i));
        nodes_whole{i} = nodes_cell;
        soma = cell_info([cell_info.cell_id]==cell).soma_coord_warped;
        soma_plane = soma(2:3);
        soma_plane(2) = soma_plane(2)*(23/16.5);
        soma_whole(i,:) = soma_plane;
    end
    
    for i = 1:n_cell
        cell = cell_list(i);
        soma = cell_info([cell_info.cell_id]==cell).soma_coord_warped;
        soma_plane = soma(2:3);
        soma_plane(2) = soma_plane(2)*(23/16.5); 
        soma_coord(i,:) = soma_plane;
        
        sample = randi(n_whole);
        
        nodes_cell = nodes_whole{sample};
        soma_cell = soma_whole(sample,:);
        nodes_cell = rotate_cell(nodes_cell,rotate_angle(i),soma_cell);
        nodes_cell = reflect_cell(nodes_cell,reflect_direction(i),soma_cell);
        
        shift = soma_coord(i,:) - soma_whole(sample,:);
        nodes_cell = nodes_cell + repmat(shift,[size(nodes_cell,1),1]);
        
        nodes_cell(nodes_cell(:,1)<1 | nodes_cell(:,2)<1 | ...
            nodes_cell(:,1)>size_x | nodes_cell(:,2)>size_y,:) = [];
        
        nodes{i} = nodes_cell;
        
        % Patch of skeleton nodes
        patch_node = patch_node + sparse(nodes_cell(:,1),round(nodes_cell(:,2)),1,size_x,size_y);
    end
    
    
elseif strcmp(permute,'orbit_rot_all')
    % orbit, rotation according to location, all cells
%     rotate_angle = permute_rotation(type,1);
    
    for i = 1:n_cell
        cell = cell_list(i);
        nodes_cell = trunkate_cell(cell,in_thr(i),out_thr(i));
        
        soma = cell_info([cell_info.cell_id]==cell).soma_coord_warped;
        soma_plane = soma(2:3);
        soma_plane(2) = soma_plane(2)*(23/16.5);
        soma_cell = soma_plane;
        
        edge_order = [1,3,4,2];
        side_xy = zeros(1,3);
        [side_x_min,side_xy(1)] = min([soma_plane(1),size_x-soma_plane(1)+1]);
        [side_y_min,side_xy(2)] = min([soma_plane(2),size_y-soma_plane(2)+1]);
        [dist_edge,side_xy(3)] = min([side_x_min,side_y_min]);
        side = edge_order((side_xy(3)-1)*2+side_xy(side_xy(3)));
        dist_edge = round(dist_edge);
        
        n_range_edge = [size_x-2*dist_edge+1,size_y-2*dist_edge+1];
        n_range = 2*sum(n_range_edge);
        soma_range = [];
        soma_range = [soma_range;ones(n_range_edge(2),1)*dist_edge,(dist_edge:size_y-dist_edge)'];
        soma_range = [soma_range;(dist_edge+1:size_x-dist_edge+1)',ones(n_range_edge(1),1)*(size_y-dist_edge+1)];
        soma_range = [soma_range;ones(n_range_edge(2),1)*(size_x-dist_edge+1),(dist_edge+1:size_y-dist_edge+1)'];
        soma_range = [soma_range;(dist_edge+1:size_x-dist_edge+1)',ones(n_range_edge(1),1)*dist_edge];
        
        orbit = [dist_edge,dist_edge;dist_edge,size_y-dist_edge+1;...
            size_x-dist_edge+1,size_y-dist_edge+1;size_x-dist_edge+1,dist_edge;...
            dist_edge,dist_edge];
        
        rotate_angle_range = [];
        rotate_angle_range = [rotate_angle_range;ones(n_range_edge(2),1)];
        rotate_angle_range = [rotate_angle_range;ones(n_range_edge(1),1)*2];
        rotate_angle_range = [rotate_angle_range;ones(n_range_edge(2),1)*3];
        rotate_angle_range = [rotate_angle_range;ones(n_range_edge(1),1)*4];
        rotate_angle_range = rem(rotate_angle_range + (4-side),4)*(pi/2);
        
      
        sample = randi(n_range);
        soma_random(i,:) = soma_range(sample,:);
        shift = soma_random(i,:) - soma_cell;
        
        rotate_angle(i) = rotate_angle_range(sample);
        
        nodes_cell = rotate_cell(nodes_cell,rotate_angle(i),soma_cell);
        nodes_cell = nodes_cell + repmat(shift,[size(nodes_cell,1),1]);
        
        nodes_cell(nodes_cell(:,1)<1 | nodes_cell(:,2)<1 | ...
            nodes_cell(:,1)>size_x | nodes_cell(:,2)>size_y,:) = [];
        
        orbits{i} = orbit;
        nodes{i} = nodes_cell;
        
        % Patch of skeleton nodes
        patch_cell = sparse(nodes_cell(:,1),round(nodes_cell(:,2)),1,size_x,size_y);
        patch_node = patch_node + patch_cell;
        patch_cell = cutedge(patch_cell,cut_x,cut_y);
        patch_cells(:,i) = reshape(patch_cell,size_x_cut*size_y_cut,1);

    end
    soma_coord = soma_random;
    
    
elseif strcmp(permute,'ref_all')
    % reflection according to location, all cells
    soma_coord = zeros(n_cell,2);
    reflect_direction = zeros(n_cell,1);
    for i = 1:n_cell
        cell = cell_list(i);
        nodes_cell = trunkate_cell(cell,in_thr(i),out_thr(i));
        
        soma = cell_info([cell_info.cell_id]==cell).soma_coord_warped;
        soma_plane = soma(2:3);
        soma_plane(2) = soma_plane(2)*(23/16.5);
        soma_cell = round(soma_plane);
        soma_coord(i,:) = soma_cell; 
        
        edge_order = [1,3,4,2];
        side_xy = zeros(1,3);
        [side_x_min,side_xy(1)] = min([soma_plane(1),size_x-soma_plane(1)]);
        [side_y_min,side_xy(2)] = min([soma_plane(2),size_y-soma_plane(2)]);
        [dist_edge,side_xy(3)] = min([side_x_min,side_y_min]);
        side = edge_order((side_xy(3)-1)*2+side_xy(side_xy(3)));
        dist_edge = round(dist_edge);    
        
        n_range_edge = [size_x-2*dist_edge-1,size_y-2*dist_edge-1];
        soma_orbit = soma_cell - [dist_edge-1,dist_edge-1];
        if soma_orbit(:,1)<0.25*n_range_edge(1) && soma_orbit(:,2)<0.25*n_range_edge(2)
            sample = randi([0,1]);
            reflect_direction(i) = sample*4;
   
        elseif soma_orbit(:,1)<0.25*n_range_edge(1) && soma_orbit(:,2)>=0.25*n_range_edge(2) ...
            && soma_orbit(:,2)<=0.75*n_range_edge(2)
            sample = randi([0,1]);
            reflect_direction(i) = sample*2;
   
        elseif soma_orbit(:,1)<0.25*n_range_edge(1) && soma_orbit(:,2)>0.75*n_range_edge(2)
            sample = randi([0,1]);
            reflect_direction(i) = sample*3;
  
        elseif soma_orbit(:,1)>=0.25*n_range_edge(1) && soma_orbit(:,1)<=0.75*n_range_edge(1) ...
            && soma_orbit(:,2)>0.75*n_range_edge(2)
            sample = randi([0,1]);
            reflect_direction(i) = sample*1;
     
        elseif soma_orbit(:,1)>0.75*n_range_edge(1) && soma_orbit(:,2)>0.75*n_range_edge(2)
            sample = randi([0,1]);
            reflect_direction(i) = sample*4;

        elseif soma_orbit(:,1)>0.75*n_range_edge(1) && soma_orbit(:,2)>=0.25*n_range_edge(2) ...
            && soma_orbit(:,2)<=0.75*n_range_edge(2)
            sample = randi([0,1]);
            reflect_direction(i) = sample*2;

        elseif soma_orbit(:,1)>0.75*n_range_edge(1) && soma_orbit(:,2)<0.25*n_range_edge(2)
            sample = randi([0,1]);
            reflect_direction(i) = sample*3;

        else
            sample = randi([0,1]);
            reflect_direction(i) = sample*1;
  
        end
                
        nodes_cell = reflect_cell(nodes_cell,reflect_direction(i),soma_cell);
        
        nodes_cell(nodes_cell(:,1)<1 | nodes_cell(:,2)<1 | ...
            nodes_cell(:,1)>size_x | nodes_cell(:,2)>size_y,:) = [];
        
        nodes{i} = nodes_cell;
        
        % Patch of skeleton nodes
        patch_node = patch_node + sparse(nodes_cell(:,1),round(nodes_cell(:,2)),1,size_x,size_y);
        patch_cell = sparse(nodes_cell(:,1),round(nodes_cell(:,2)),1,size_x,size_y);
        patch_cell = cutedge(patch_cell,cut_x,cut_y);
        patch_cells(:,i) = reshape(patch_cell,size_x_cut*size_y_cut,1);
        
    end
    
elseif permute == 0 
    soma_coord = zeros(n_cell,2);
    for i = 1:n_cell
        cell = cell_list(i);
        nodes_cell = trunkate_cell(cell,in_thr(i),out_thr(i));
        nodes{i} = nodes_cell;

        soma = cell_info([cell_info.cell_id]==cell).soma_coord_warped;
        soma_plane = soma(2:3);
        soma_plane(2) = soma_plane(2)*(23/16.5);   
        soma_coord(i,:) = soma_plane;
        
        % Patch of skeleton nodes
        patch_cell = sparse(nodes_cell(:,1),round(nodes_cell(:,2)),1,size_x,size_y);
        patch_node = patch_node + patch_cell;
        patch_cell = cutedge(patch_cell,cut_x,cut_y);
        patch_cells(:,i) = reshape(patch_cell,size_x_cut*size_y_cut,1);

    end
end

%% Plot
if visual == 1
color_list = distinguishable_colors(n_cell);

figure();
for i = 1:length(cell_list)
    nodes_cell = nodes{i};
    soma_cell = soma_coord(i,:);
    
    hold on;
    plot(nodes_cell(:,2),nodes_cell(:,1),'.','Color',color_list(i,:),'MarkerSize',3);

end

for i = 1:n_cell
    hold on;
    plot(soma_coord(i,2),soma_coord(i,1),'ko','MarkerFaceColor',color_list(i,:),'MarkerEdgeColor',color_list(i,:),'MarkerSize',12);
end

xlim([0,size_y]);
ylim([0,size_x]);

fill([0,0,size_y,size_y],[0,cut_x(1),cut_x(1),0],[192,192,192]/255,'FaceAlpha',0.7,'EdgeColor','none')
fill([0,0,size_y,size_y],[size_x-cut_x(2),size_x,size_x,size_x-cut_x(2)],[192,192,192]/255,'FaceAlpha',0.7,'EdgeColor','none')
fill([0,0,cut_y(1),cut_y(1)],[cut_x(1),size_x-cut_x(2),size_x-cut_x(2),cut_x(1)],[192,192,192]/255,'FaceAlpha',0.7,'EdgeColor','none')
fill([size_y-cut_y(2),size_y-cut_y(2),size_y,size_y],[cut_x(1),size_x-cut_x(2),size_x-cut_x(2),cut_x(1)],[192,192,192]/255,'FaceAlpha',0.7,'EdgeColor','none')

hold on;
plot([cut_y(1),size_y-cut_y(2),size_y-cut_y(2),cut_y(1),cut_y(1)],[cut_x(1),cut_x(1)...
    ,size_x-cut_x(2),size_x-cut_x(2),cut_x(1)],'LineWidth',2,'LineStyle','-','Color',[192,192,192]/255);
hold on;
plot([0,0,size_y,size_y,0],[0,size_x...
    ,size_x,0,0],'LineWidth',2,'LineStyle','-','Color',[192,192,192]/255);

axis equal;
axis off;
end

%%

% Crop patch
patch_node = cutedge(patch_node,cut_x,cut_y);

% Divide into grid squares
grid_x = round(size_x_cut/grid_size);
grid_y = round(size_y_cut/grid_size);

% Calculate density
density_patch = zeros(grid_x,grid_y);
winner_patch = zeros(grid_x,grid_y);
winner_id = zeros(grid_x,grid_y);
for m = 1:grid_x
    for n = 1:grid_y
        if m == grid_x && n ~= grid_y
            grid_square = patch_node((m-1)*grid_size+1:end,(n-1)*grid_size+1:n*grid_size);
            density_square = sum(sum(grid_square))/numel(grid_square);
            density_patch(m,n) = density_square;
            
            density_cell = zeros(1,n_cell);
            for i = 1:n_cell
                patch_cell = reshape(patch_cells(:,i),size_x_cut,size_y_cut);
                grid_square = patch_cell((m-1)*grid_size+1:end,(n-1)*grid_size+1:n*grid_size);
                density_square = sum(sum(grid_square))/numel(grid_square);
                density_cell(i) = density_square;
            end
            [winner_density,winner_idx] = max(density_cell);
            winner_patch(m,n) = winner_density/sum(density_cell);
            winner_id(m,n) = (winner_density>=0.001)*(winner_patch(m,n)>0.5)*winner_idx;
%             winner_id(m,n) = (winner_density>=0.001)*winner_idx;
           
        elseif m ~= grid_x && n == grid_y
            grid_square = patch_node((m-1)*grid_size+1:m*grid_size,(n-1)*grid_size+1:end);
            density_square = sum(sum(grid_square))/numel(grid_square);
            density_patch(m,n) = density_square;
            
            density_cell = zeros(1,n_cell);
            for i = 1:n_cell
                patch_cell = reshape(patch_cells(:,i),size_x_cut,size_y_cut);
                grid_square = patch_cell((m-1)*grid_size+1:m*grid_size,(n-1)*grid_size+1:end);
                density_square = sum(sum(grid_square))/numel(grid_square);
                density_cell(i) = density_square;
            end
            [winner_density,winner_idx] = max(density_cell);
            winner_patch(m,n) = winner_density/sum(density_cell);
            winner_id(m,n) = (winner_density>=0.001)*(winner_patch(m,n)>0.5)*winner_idx;
%             winner_id(m,n) = (winner_density>=0.001)*winner_idx;
            
        elseif m == grid_x && n == grid_y
            grid_square = patch_node((m-1)*grid_size+1:end,(n-1)*grid_size+1:end);
            density_square = sum(sum(grid_square))/numel(grid_square);
            density_patch(m,n) = density_square;
            
            density_cell = zeros(1,n_cell);
            for i = 1:n_cell
                patch_cell = reshape(patch_cells(:,i),size_x_cut,size_y_cut);
                grid_square = patch_cell((m-1)*grid_size+1:end,(n-1)*grid_size+1:end);
                density_square = sum(sum(grid_square))/numel(grid_square);
                density_cell(i) = density_square;
            end
            [winner_density,winner_idx] = max(density_cell);
            winner_patch(m,n) = winner_density/sum(density_cell);
            winner_id(m,n) = (winner_density>=0.001)*(winner_patch(m,n)>0.5)*winner_idx;
%             winner_id(m,n) = (winner_density>=0.001)*winner_idx;
            
        else
            grid_square = patch_node((m-1)*grid_size+1:m*grid_size,(n-1)*grid_size+1:n*grid_size);
            density_square = sum(sum(grid_square))/numel(grid_square);
            density_patch(m,n) = density_square;
            
            density_cell = zeros(1,n_cell);
            for i = 1:n_cell
                patch_cell = reshape(patch_cells(:,i),size_x_cut,size_y_cut);
                grid_square = patch_cell((m-1)*grid_size+1:m*grid_size,(n-1)*grid_size+1:n*grid_size);
                density_square = sum(sum(grid_square))/numel(grid_square);
                density_cell(i) = density_square;
            end
            [winner_density,winner_idx] = max(density_cell);
            winner_patch(m,n) = winner_density/sum(density_cell);
            winner_id(m,n) = (winner_density>=0.001)*(winner_patch(m,n)>0.5)*winner_idx;
%             winner_id(m,n) = (winner_density>=0.001)*winner_idx;
            
        end
    end
end

toc;

end


