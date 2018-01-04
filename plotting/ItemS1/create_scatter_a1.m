%% create 0_a scatter plot


[onsac, offsac] = get_sac_strat(cell_info);
off_sac = 28;
on_sac = 62;
center = 47;

cell_types = {'7o','7id','7ir','7iv','37r','37v','37d','37c','63','1ws','1wt','1no','1ni','2an','2aw','2o','2i','3o','3i','4on','4ow','4i','5to','5ti','5so','5si','6sn','6sw','6t','8w','8n','9n','9w','51','25','85','73','72','27','81i','81o','82wo','82wi','82n','915','28','91'};

figure();
for t = 1:length(cell_types)
    idx = find(strcmp({cell_info.type},cell_types{t}));

    cell_list = zeros(1,length(idx));
    cell_stat_on = zeros(1,length(idx));
    cell_stat_off = zeros(1,length(idx));
    
    for i = 1:length(idx)
        cell_list(i) = cell_info(idx(i)).cell_id;
        
        cell_info_elem = get_cell_info(cell_info,cell_list(i));
        cell_stat_on(i) = cell_info_get_strat_property(cell_info_elem, 'corr', off_sac, on_sac, center, true, onsac);
        cell_stat_off(i) = cell_info_get_strat_property(cell_info_elem, 'corr', off_sac, on_sac, center, true, offsac);
    end
    
    hold on;
    if t <= 9
        scatter(cell_stat_on, cell_stat_off,80,'filled');
        l{t} = cell_types{t};
    else
        scatter(cell_stat_on, cell_stat_off,80,[119 136 153]/255,'filled');
    end
    
end

rad = 0:0.01:pi/2;
x = 0.835*cos(rad);
y = 0.835*sin(rad);
hold on;
plot(x,y,'linewidth',1,'Color',[0 0 0]);

legend(l,'FontName','Arial','FontSize',20);
set(gca,'XTick',0:0.2:1,'YTick',0:0.2:1,'FontSize',20,'FontName','Arial');
title('Split a-1','FontSize',25,'FontName','Arial');
xlabel('On SAC similarity','FontSize',25,'FontName','Arial');
ylabel('Off SAC similarity','FontSize',25,'FontName','Arial');
axis square;

% print('-dpng','../Classification/histograms/a-1.png','-r300');
% print('-depsc','../Classification/histograms/a-1.eps','-r500');
% print('-dsvg','../Classification/histograms/a-1.svg','-r300');

   