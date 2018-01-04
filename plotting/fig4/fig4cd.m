%% Figure 4 histograms
load('./cell_types');
load('./group_idx');

figure();
% Number of cells (bar)
load('./numcell');

numcell_bar_axes = axes('Position',[0.3 0.74 0.31 0.25]);

cell_types_sort = cell_types(sort_idx);

bar(find(masland_idx(sort_idx)==1),numcell_sort(masland_idx(sort_idx)),1,'FaceColor',[183 36 103]/255);
hold on;
bar(find(seminovel_idx(sort_idx)==1),numcell_sort(seminovel_idx(sort_idx)),1,'FaceColor',[82 79 161]/255);
hold on;
bar(find(novel_idx(sort_idx)==1),numcell_sort(novel_idx(sort_idx)),1,'FaceColor',[109 200 191]/255);
hold on;
bar(find(invalid_idx(sort_idx)==1),numcell_sort(invalid_idx(sort_idx)),1,'FaceColor',[131 108 121]/255);

xlim([0.5,47.5]);
ylim([0,30.5]);
set(gca,'XTick',1:47,'XTickLabel',cell_types_sort,'TickLength',[0.005 0],'FontSize',8,'FontName','Arial');
ylabel('Number of Cells','FontSize',15,'FontName','Arial');
box off;
ax = gca;
ax.XTickLabelRotation=90;
legend({'Sanes & Masland, 2015','Semi-novel types','Novel types','Others'},'FontName','Arial','FontSize',9,'EdgeColor','none');

% Number of cells (histogram)
numcell_hist_axes = axes('Position',[0.63 0.74 0.07 0.25]);

[counts,bins] = hist(numcell_sort,0:1:30);
h = barh(bins,counts,1,'FaceColor',[119 136 153]/255);
ylim([0,30.5]);
xlim([0,7]);
set(gca,'FontName','Arial','FontSize',8);
xlabel('Number of Clusters','FontSize',10,'FontName','Arial');
box off; 


% Coverage factor (bar)
load('./coverage_factor.mat');

cf_bar_axes = axes('Position',[0.3 0.4 0.31 0.25]);
    
hull_cf_sort = cf(sort_idx);

bar(find(masland_idx(sort_idx)==1),hull_cf_sort(masland_idx(sort_idx)),1,'FaceColor',[183 36 103]/255);
hold on;
bar(find(seminovel_idx(sort_idx)==1),hull_cf_sort(seminovel_idx(sort_idx)),1,'FaceColor',[82 79 161]/255);
hold on;
bar(find(novel_idx(sort_idx)==1),hull_cf_sort(novel_idx(sort_idx)),1,'FaceColor',[109 200 191]/255);
hold on;
bar(find(invalid_idx(sort_idx)==1),hull_cf_sort(invalid_idx(sort_idx)),1,'FaceColor',[131 108 121]/255);
xlim([0.5,47.5]);
ylim([0,10]);
set(gca,'XTick',1:47,'XTickLabel',cell_types_sort,'YTick',0:2:11,'TickLength',[0.005 0],'FontSize',8,'FontName','Arial');
ylabel('Coverage Factor','FontSize',15,'FontName','Arial');
box off;
ax = gca;
ax.XTickLabelRotation=90;

% Coverage factor (histogram)
cf_hist_axes = axes('Position',[0.63 0.4 0.07 0.25]);

[counts,bins] = hist(hull_cf_sort,0:0.4:11);
barh(bins,counts,1,'FaceColor',[119 136 153]/255)
xlim([0,12]);
ylim([0,10]);
set(gca,'FontName','Arial','FontSize',8);
xlabel('Number of Clusters','FontSize',10,'FontName','Arial');
box off; 


% % Winner fraction (bar)
% load('./data/winner_fraction');
% 
% wf_bar_axes = axes('Position',[0.3 0.06 0.31 0.25]);
% 
% winner_fraction = winner_fraction(sort_idx,:);
% winner_cf_sort = mean(winner_fraction,2,'omitnan');
% 
% bar(find(masland_idx(sort_idx)==1),winner_cf_sort(masland_idx(sort_idx)),1,'FaceColor',[183 36 103]/255);
% hold on;
% bar(find(seminovel_idx(sort_idx)==1),winner_cf_sort(seminovel_idx(sort_idx)),1,'FaceColor',[82 79 161]/255);
% hold on;
% bar(find(novel_idx(sort_idx)==1),winner_cf_sort(novel_idx(sort_idx)),1,'FaceColor',[109 200 191]/255);
% hold on;
% bar(find(invalid_idx(sort_idx)==1),winner_cf_sort(invalid_idx(sort_idx)),1,'FaceColor',[131 108 121]/255);
% xlim([0.5,47.5]);
% hold on;
% plot([0.5,47.5],[0.5,0.5],'--k');
% 
% mean_wf = mean(winner_fraction,2,'omitnan');
% stdev = std(winner_fraction');
% sem = stdev/size(winner_fraction,2)^0.5;
% 
% hold on;
% errorbar(1:47,mean_wf,sem,sem,'Color','k','LineStyle','none','LineWidth',1,'CapSize',5);
% set(gca,'XTick',1:47,'XTickLabel',cell_types_sort,'YTick',0:0.2:1,'TickLength',[0.005 0],'FontSize',8,'FontName','Arial');
% ylabel('Winner Fraction','FontName','Arial','FontSize',15);
% ylim([0,1.025]);
% box off;
% ax = gca;
% ax.XTickLabelRotation=90;
% 
% % Winner fraction (histogram)
% wf_hist_axes = axes('Position',[0.63 0.06 0.07 0.25]);
% 
% [counts,bins] = hist(winner_cf_sort,0:0.05:1);
% barh(bins,counts,1,'FaceColor',[119 136 153]/255)
% 
% hold on;
% plot([0,12],[0.5,0.5],'--k');
% set(gca,'FontName','Arial','FontSize',8,'YTick',0:0.2:1,'XTick',0:2:12);
% xlim([0,12]);
% ylim([0,1.025]);
% xlabel('Number of Clusters','FontSize',10,'FontName','Arial');
% box off; 
