%% Figure 2h

load('seg_idx_marginal-central_gc_wods');
% load('./data/sep_idx/sep_idx_marginal-central_bc.mat');

figure(); 
imagesc(seg_idx(:,1:36));
c = colorbar();
c.Limits(1) = 0;

set(gca,'XTick',1:5:36,'XTickLabels',0.5:0.05:0.9,...
    'YTick',1:5:51,'YTickLabels',0:0.05:0.5,...
    'TickLength',[0.015,0.015],'FontName','Arial','FontSize',15);
xlabel('Inner Boundary (IPL Depth)','FontName','Arial','FontSize',20);
ylabel('Outer Boundary (IPL Depth)','FontName','Arial','FontSize',20);
xlim([1,36]);
ylim([1,51]);

hold on;
plot([13,13],[1,51],'--k');
hold on;
plot([1,36],[29,29],'--k');
axis square;



