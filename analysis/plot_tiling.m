function [p lims] = plot_tiling( cell_ids, downsampling, yz_not_yx, color_arr)

if ~exist('color_arr','var') || isempty(color_arr)
    color_arr = default_color_arr();
end

skel_point_size = 1;


if yz_not_yx
    soma_point_size = 75;
else
    soma_point_size = 75;
end

cell_info_filename = 'cell_info_clustering.20160623warpedSomaCorrection.mat';

ci = load(cell_info_filename, 'cell_info');
ci = ci.cell_info;


if yz_not_yx
    lims = [0 5376 ... %y limits
        0 (3456*23/16.5)]; %z limits
else
    lims = [0 5376 ... %y limits
        445 1166]; %x limits (
end

for i=1:length(cell_ids)
    cell_id = cell_ids(i);
    
    %Plotting skeleton
    skel_nodes = load_skeleton_points( cell_id );
    
    skel_nodes = skel_nodes(1:downsampling:end, :);
    size_column = repmat(skel_point_size, [size(skel_nodes,1) 1]);
    
   
    if yz_not_yx
        hold on;
        plot(skel_nodes(:,2), skel_nodes(:,3), '.', ...
            'MarkerSize',skel_point_size,'Color',color_arr(i,:));
        xlim([0 5376]);
        ylim([0 (3456*23/16.5)]);
        
    else
        hold on;
        plot(skel_nodes(:,2), skel_nodes(:,1), '.', ...
            'MarkerSize',skel_point_size,'Color',color_arr(i,:));
        xlim([0 5376]);
        ylim([200 1166]);
        
    end
    
end

for i=1:length(cell_ids)
    cell_id = cell_ids(i);
    
    cell_info_filename = 'cell_info_clustering.20160623warpedSomaCorrection.mat';
    
    ci = load(cell_info_filename, 'cell_info');
    ci = ci.cell_info;
    
    %Plotting soma
    ci_id = [ci.cell_id] == cell_id;
    cell_soma_coords = ci(ci_id).soma_coords_warped_mip2_zscaled;
    cell_soma_coords = cell_soma_coords([3 1 2]);
    
    if yz_not_yx
        hold on;
        plot(cell_soma_coords(:,2), cell_soma_coords(:,3), '.', ...
            'MarkerSize', soma_point_size, 'Color',color_arr(i,:))
        xlim([0 5376]);
        ylim([0 (3456*23/16.5)]);
        
    else
        hold on;
        plot(cell_soma_coords(:,2), cell_soma_coords(:,1), '.', ...
            'MarkerSize', soma_point_size, 'Color',color_arr(i,:))
        xlim([0 5376]);
        ylim([200 1166]);
        
    end
end

set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% set(gca,'visible','off');
%   axis(lims)
axis off;

end


function points = load_skeleton_points( cell_id )

skel_fmt = '../skeletons/skel_%d.mat';

load_obj = load(sprintf(skel_fmt, cell_id), 'n');

points = load_obj.n(:,[3,2,1]);
points(:,3) = 1 * points(:,3) * 23 / 16.5;
end

function color_arr = default_color_arr()

color_arr =  [ 0.0         0.0       0.0 ; ...
    1.0         1.0       0.384314; ...
    1.0         0.623529  1.0     ; ...
    0.0         0.839216  1.0     ; ...
    0.843137    0.266667  0.0     ; ...
    0.0         0.501961  0.160784; ...
    0.0         0.372549  0.835294; ...
    0.576471    0.0       0.407843; ...
    1.0         0.796078  0.709804; ...
    0.654902    0.52549   0.0     ; ...
    0.0         1.0       0.772549; ...
    0.0         0.494118  0.52549 ; ...
    0.32549     0.258824  0.0     ; ...
    0.54902     0.529412  0.576471; ...
    0.407843    0.0       0.0     ; ...
    0.0         0.0       0.352941; ...
    1.0         0.172549  0.435294; ...
    0.968627    1.0       0.964706; ...
    0.0         0.239216  0.145098; ...
    0.568627    0.670588  1.0     ; ...
    0.631373    0.713725  0.517647; ...
    0.376471    0.313725  0.317647; ...
    0.756863    0.160784  0.886275; ...
    1.0         0.556863  0.117647; ...
    0.0         0.247059  0.317647; ...
    0.47451     0.470588  0.376471; ...
    0.180392    0.0       0.109804; ...
    0.00392157  0.686275  0.6     ; ...
    1.0         0.901961  1.0     ; ...
    0.466667    0.854902  0.0];
end
