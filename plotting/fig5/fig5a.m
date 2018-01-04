load('./cell_hull');

type = '25';

gc = cell_info_typedef_gc();
cell_list = gc(strcmp({gc.name},type)).cells;
n_cell = length(cell_list);

size_x = 5376;
size_y = round(3456*23/16.5);
cut_x = [1000,1000]; cut_y = [1000,1000];

patch_overlap = zeros(size_x,size_y);
for i = 1:n_cell
    cell = cell_list(i);
    hull_cell = cell_hull{cell};
    
    in_hull = poly2mask(hull_cell(:,1),hull_cell(:,2),size_y,size_x);
    patch_overlap = patch_overlap + transpose(in_hull);
end

figure();
i = imagesc(patch_overlap);
% hold on;
% plot([cut_y(1),size_y-cut_y(2),size_y-cut_y(2),cu6t_y(1),cut_y(1)],[cut_x(1),cut_x(1)...
%         ,size_x-cut_x(2),size_x-cut_x(2),cut_x(1)],'LineWidth',2,'LineStyle','-','Color',[192,192,192]/255);
% hold on;
% fill([0,0,size_y,size_y],[0,cut_x(1),cut_x(1),0],[192,192,192]/255,'FaceAlpha',0.7,'EdgeColor','none')
% fill([0,0,size_y,size_y],[size_x-cut_x(2),size_x,size_x,size_x-cut_x(2)],[192,192,192]/255,'FaceAlpha',0.7,'EdgeColor','none')
% fill([0,0,cut_y(1),cut_y(1)],[cut_x(1),size_x-cut_x(2),size_x-cut_x(2),cut_x(1)],[192,192,192]/255,'FaceAlpha',0.7,'EdgeColor','none')
% fill([size_y-cut_y(2),size_y-cut_y(2),size_y,size_y],[cut_x(1),size_x-cut_x(2),size_x-cut_x(2),cut_x(1)],[192,192,192]/255,'FaceAlpha',0.7,'EdgeColor','none')

axis equal;
set(gca,'XLim',[0,size_y],'YLim',[0,size_x]);
c = colorbar();
c.Ticks = 0:5;
c.FontName = 'Arial';
c.FontSize = 25;
axis off;
