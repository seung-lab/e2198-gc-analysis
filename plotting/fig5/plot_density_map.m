function plot_density_map(density_patch,nodes,soma_coord)
    %% Plot 3D density map
    
    figure();
    surf(density_patch,'EdgeColor','interp','FaceColor','interp','FaceAlpha',0.5);
    colormap winter;
    colorbar();
    
    % Axis control
    size_x = 5376;
    size_y = round(3456*23/16.5);
    cut_x = [1000,1000]; cut_y = [1000,1000];
    size_x_cut = size_x-sum(cut_x);
    size_y_cut = size_y-sum(cut_y);
    grid_size = 600;
    
    xlim([0.5,size(density_patch,2)+0.5]);
    xlim([0.5,size_y_cut/grid_size + 0.5]);
    ylim([0.5,size(density_patch,1)+0.5]);
    ylim([0.5,size_x_cut/grid_size + 0.5]);
    zaxis = get(gca,'ZLim');
    zlim([0,0.9]);
    caxis([0,0.9]);
    
    % Axis equal
    dar = get(gca,'DataAspectRatio');
    set(gca,'Color','none','DataAspectRatio',[1,1,1/max(dar(1:2))],...
        'XTick',1:size(density_patch,2),'XTickLabel',{},...
        'YTick',1:size(density_patch,1),'YTickLabel',{},...
        'FontName','Arial');
%     set(gcf,'Color','none','inverthardcopy','off');
%     
    
    % Plot skeletons    
%     n_cell = size(soma_coord,1);
%     
%     color_list = [176,178,218;249,185,187;179,135,188;68,182,73;197,23,137]/255;
%     soma_coord(:,1) = soma_coord(:,1) - cut_x(1);
%     soma_coord(:,2) = soma_coord(:,2) - cut_y(1);
%     for i = 1:n_cell
%         nodes_cell = nodes{i};
%         soma_cell = soma_coord(i,:);
%         
%         nodes_cell(nodes_cell(:,1)<=cut_x(1) | nodes_cell(:,2)<=cut_y(1) | ...
%             nodes_cell(:,1)>size_x-cut_x(2) | nodes_cell(:,2)>size_y-cut_y(2),:) = [];
%         
%         nodes_cell(:,1) = nodes_cell(:,1) - cut_x(1);
%         nodes_cell(:,2) = nodes_cell(:,2) - cut_y(1);
%         hold on;
%         plot(nodes_cell(:,2)/600 + 0.5,nodes_cell(:,1)/600 + 0.5,'.','Color',color_list(i,:),'MarkerSize',3);
%         hold on;
%         plot(soma_coord(i,2)/600 + 0.5,soma_coord(i,1)/600 + 0.5,'ko','MarkerFaceColor',color_list(i,:),'MarkerEdgeColor',color_list(i,:),'MarkerSize',12);
%     end
        
end