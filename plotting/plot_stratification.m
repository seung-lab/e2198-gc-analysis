function plot_stratification(cell_list,average,bin_resolution,var,vol)
all_types = cell_info_typedef();

if ~exist('bin_resolution','var') || isempty(bin_resolution)
    bin_resolution = 1;
end

if ~exist('vol','var') || isempty(vol)
    vol = 0;
end
legend_on = 0; 

figure();
if iscell(cell_list)
    group_list = cell_list;
    n_group = length(group_list);
%     color_list = distinguishable_colors(n_group);
    color_list = get(gca,'colororder');
    
    for i = 1:n_group
        group = group_list{i};
        
        if isstr(group)
            legend_on = 1;
            cell_list = all_types(strcmp({all_types.name},group)).cells;
            
        else
            cell_list = group;
        end
        
        n_cell = length(cell_list);
        x_idx = 1:bin_resolution:723;
        
        strat_avg = zeros(1,length(x_idx));
        for j = 1:n_cell
            cell = cell_list(j);
            load(strcat('../skeletons/skel_',num2str(cell)));
            
            if vol == 1
                n = p;
            end
            
            [x strat_y] = strat_skel(n,1);
            strat_y = strat_y/sum(strat_y)*514.71;
            strat = strat_y;
            
            x = x(x_idx);
            strat_bin = strat(x_idx);
            for b = 1:bin_resolution-1
                idx_next = x_idx+b;
                if idx_next(end) > 723
                    idx_next(end) = [];
                    temp = strat(idx_next);
                    temp(end+1) = 0;
                else
                    temp = strat(idx_next);
                end
                
                strat_bin = strat_bin + temp;
            end
            
            strat_bin = strat_bin/bin_resolution;
            
            strat_type(j,:) = strat_bin;
         
            if average
                strat_avg = strat_avg + strat_bin;
            else
                hold on;
                plot(x/100,strat_bin,'LineWidth',1.5);
                
            end
        end
        
        if average
            strat_avg = strat_avg/n_cell;
            
            sem = std(strat_type)/n_cell^0.5;
            sem = std(strat_type);
            
            if var == 1
                shade_bound = [(strat_avg - sem)',sem',sem'];
                shade = area(x/100,shade_bound,'LineStyle','none');
                shade(1).FaceAlpha = 0;
                shade(2).FaceColor = color_list(i,:);
                shade(2).FaceAlpha = 0.2;
                shade(3).FaceColor = color_list(i,:);
                shade(3).FaceAlpha = 0.2;
            end
            
            hold on;
            pl(i) = plot(x/100,strat_avg,'LineWidth',1.5);
           
%             pl(i) = plot(x/100,strat_avg,'LineWidth',1.5,'Color',color_list(i,:));
%             save('strat_avg','x','strat_avg');
            
           
        end
        
    end
    
% cell_list not in cell format
else
    if isstr(cell_list)
        type = cell_list;
        cell_list = all_types(strcmp({all_types.name},type)).cells;
    end
    
    n_cell = length(cell_list);
    
    strat_avg = zeros(1,length(1:bin_resolution:723));
    for i = 1:length(cell_list)
        cell = cell_list(i);
        load(strcat('../skeletons/skel_',num2str(cell)));
        if vol == 1
            n = p;
        end
        
        [x strat_y] = strat_skel(n,1);
        strat_y = strat_y/sum(strat_y)*514.71;
        
        strat = strat_y;
            
        x_idx = 1:bin_resolution:723;
        x = x(x_idx);
        strat_bin = strat(x_idx);
        for b = 1:bin_resolution-1
            idx_next = x_idx+b;
            if idx_next(end) > 723
                idx_next(end) = [];
                temp = strat(idx_next);
                temp(end+1) = 0;
            else
                temp = strat(idx_next);
            end
            
            strat_bin = strat_bin + temp;
        end
        
        strat_bin = strat_bin/bin_resolution;
        
        if average
            hold on;
            strat_avg = strat_avg + strat_bin;
        
        else
            hold on;
            plot(x/100,strat_bin,'LineWidth',1.5);
            
        end
    end
    
    if average
          strat_avg = strat_avg/n_cell;
            
            sem = std(strat_type)/n_cell^0.5;
            
            if var == 1
                shade_bound = [(strat_avg - sem)',sem',sem'];
                shade = area(x/100,shade_bound,'LineStyle','none');
                shade(1).FaceAlpha = 0;
                shade(2).FaceColor = [192,192,192]/255;
                shade(2).FaceAlpha = 0.3;
                shade(3).FaceColor = [192,192,192]/255;
                shade(3).FaceAlpha = 0.3;
            end
            
            hold on;
            plot(x/100,strat_avg,'LineWidth',1.5);
    end
    
end
strat_ylim = get(gca,'YLim');
line([0.28 0.28],strat_ylim,'LineStyle','--','Color','k');
line([0.62 0.62],strat_ylim,'LineStyle','--','Color','k');
line([0.47 0.47],strat_ylim,'LineStyle',':','Color','k');
% xlim([-0.1, 1.1]);
xlim([0,1]);
ylim([0,strat_ylim(2)]);
set(gca,'FontName','Arial','FontSize',25,'XTick',[0,0.28,0.47,0.62,1],'Box','on')
xlabel('IPL Depth','FontName','Arial','FontSize',30);
ylabel('Skeleton Density','FontName','Arial','FontSize',30);

% legend(pl,{'outer central','inner-outer central','inner central','outer marginal','inner marginal'},'FontName','Arial','FontSize',15,'Box','off')

if legend_on == 1
    legend(pl,group_list,'FontName','Arial','FontSize',15,'Box','off');
end

end