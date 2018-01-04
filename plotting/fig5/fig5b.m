type = '25';
load('./cell_info');

colormap('default');
c = colormap;

trunk_thr = get_trunkthr(type);
[density_overlap, area_overlap] = density_hull_overlap(cell_info,type,trunk_thr,-0.2,0,0);

n_overlap = length(density_overlap)-1;
area_coverage = zeros(n_overlap,1);

figure();
for i = 1:n_overlap
    area_list = area_overlap{i+1};
    area_coverage(i) = sum(area_list)*66*66/10^6;
    
    hold on;
    bar(i,area_coverage(i),'FaceColor',c(round(1+(63/5)*i),:));
end


xlim([0.5,n_overlap+0.5]);
set(gca,'XTick',1:n_overlap,'FontName','Arial','FontSize',15);
xlabel('Coverage','FontName','Arial','FontSize',20);
ylabel('Area (um^2)','FontName','Arial','FontSize',20);
axis square;