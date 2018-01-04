type = '25';
load('./cell_info');

colormap('default');
c = colormap;

trunk_thr = get_trunkthr(type);
[density_overlap, area_overlap] = density_hull_overlap(cell_info,type,trunk_thr,-0.2,0,0);

n_overlap = length(density_overlap)-1;
density_coverage = zeros(n_overlap,1);

figure();
for i = 1:n_overlap
    density_list = density_overlap{i+1};
    density_list = density_list*0.0957*10^6/(66*66);
    
    density_coverage(i) = mean(density_list);
    
%     ptile = prctile(density_list,[5,95]);
%     neg = density_coverage(i) - ptile(1);
%     pos = ptile(2) - density_coverage(i);
    neg = std(density_list)/length(density_list)^0.5;
    pos = std(density_list)/length(density_list)^0.5; 
    length(density_list)
    
    hold on;
    bar(i,density_coverage(i),'FaceColor',c(round(1+(63/5)*i),:));
    hold on;
    errorbar(i,density_coverage(i),neg,pos,'Color','k','LineStyle','none','LineWidth',1.5,'CapSize',8);
end

xlim([0.5,n_overlap+0.5]);
% ylim([0,0.08]);
set(gca,'XTick',1:n_overlap,'FontName','Arial','FontSize',20);
xlabel('Coverage','FontName','Arial','FontSize',25);
ylabel('Density','FontName','Arial','FontSize',25);
% legend({'Mean Density','\pm S.E.M.'},'FontName','Arial','FontSize',15)
axis square;