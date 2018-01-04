load('seg_idx_inner-outer_bc_0_100');

plot(0.05:0.01:0.95,seg_idx(6:96),'LineWidth',1.5,'Color','k');
y_lim = get(gca,'YLim');
y_lim(1) = 0;

hold on;
plot([0.47,0.47],y_lim,':k');
hold on;
plot([0.28,0.28],y_lim,'--k');
hold on;
plot([0.62,0.62],y_lim,'--k');


set(gca,'XTick',0.1:0.1:0.9,'YTick',0:2:y_lim(2),'FontName','Arial','FontSize',15);
xlim([0.05,0.95]);
xlabel({'Inner/Outer Boundary','(IPL Depth)'},'FontName','Arial','FontSize',25);
ylabel('Segregation Index','FontName','Arial','FontSize',25);
axis square;
