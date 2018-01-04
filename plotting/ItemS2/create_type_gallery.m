%% create_type_gallery

load cell_types;
load cell_info;
load temporal_response;
load tempresp_avg;
load ca_export_museum;
load ca_cells;

gc = cell_info_typedef_gc();

colors = distinguishable_colors(30);
types_to_plot = cell_types;
% types_to_plot = {'1ni'};

for t = 1:length(types_to_plot)
    type_idx = find(strcmp(cell_types,types_to_plot{t}));
    cells = gc(type_idx).cells;
    
    figure();
    type_name_axes = axes('Position',[0.05 0.9 0.4 0.1]);
    text(0.02,0.5,types_to_plot{t},'FontName','Arial','FontSize',45);
    axis off;
    
    colors = distinguishable_colors(30);
    
    % Tiling planar view panel
    tile_plane_axes = axes('Position',[0.05 0.33 0.4 0.6]);
    
    plot_tiling(cells,1,1,colors);
    
    labela_axes = axes('Position',[0.03 0.89 0.02 0.02]);
    text(0.5,0.5,'a','Units','normalized','FontName','Arial','FontSize',25,'FontWeight','bold');
    axis off;
    
    
    % Tiling side view panel
    tile_side_axes = axes('Position',[0.05 0.1 0.4 0.2]);
    
    plot_tiling(cells,1,0,colors);
    
    line([0 5376],[744 744],'LineStyle','--','LineWidth',1.5,'Color',[192 192 192]/255);
    line([0 5376],[919 919],'LineStyle','--','LineWidth',1.5,'Color',[192 192 192]/255);
    
    labelb_axes = axes('Position',[0.03 0.28 0.02 0.02]);
    text(0.5,0.5,'b','Units','normalized','FontName','Arial','FontSize',25,'FontWeight','bold');
    axis off;


    % Stratification panel
    strat_axes = axes('Position',[0.5 0.55 0.45 0.4]);
        
    for i = 1:length(cells)
        strat_id = cells(i);
        
        bin_res = 4;
        x = cell_info([cell_info.cell_id]==strat_id).strat_nrml(:,1)/100;
        strat = cell_info([cell_info.cell_id]==strat_id).strat_nrml(:,2)*514.71;
        
        x_idx = 1:bin_res:723;
        x = x(x_idx);
        strat_bin = strat(x_idx);
        for b = 1:bin_res-1
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
        
        strat_bin = strat_bin/bin_res;
        
        hold on;
        plot(x,strat_bin,'Color',colors(i,:));
    end
    
    % Plot average
    strat_id = type_idx + 99900;
    
    bin_res = 4;
    x = cell_info([cell_info.cell_id]==strat_id).strat_nrml(:,1)/100;
    strat = cell_info([cell_info.cell_id]==strat_id).strat_nrml(:,2)*514.71;
    
    x_idx = 1:bin_res:723;
    x = x(x_idx);
    strat_bin = strat(x_idx);
    for b = 1:bin_res-1
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
    
    strat_bin = strat_bin/bin_res;
    hold on;
    plot(x,strat_bin,'Color',[119 136 153]/255,'LineWidth',4);
    
    
    strat_ylim = get(gca,'YLim');
    line([0.28 0.28],strat_ylim,'LineStyle','--','LineWidth',1.5,'Color',[192 192 192]/255);
    line([0.62 0.62],strat_ylim,'LineStyle','--','LineWidth',1.5,'Color',[192 192 192]/255);
    
    xlim([-0.1 1.1]);
    set(gca,'XTick',[0 0.28 0.47 0.62 1],'XTickLabel',{'0','Off SAC','0.47','On SAC','1'},'FontName','Arial','FontSize',15);
    xlabel('IPL depth','FontName','Arial','FontSize',20);
    ylabel('Skeleton density','FontName','Arial','FontSize',20);
    
    labelc_axes = axes('Position',[0.45 0.95 0.02 0.02]);
    text(0.5,0.5,'c','Units','normalized','FontName','Arial','FontSize',25,'FontWeight','bold');
    axis off;
    
    % Tuning curve panel
    tuning_axes = polaraxes('Position',[0.48 0.1 0.2 0.3]);
 
    theta = pi/4 * [0:7].';
    for i = 1:length(cells)
        cell = cells(i);
        ca_id = ca_cells(find(ca_cells(:,2)==cell),1);

        if isempty(ca_id)
            continue;
        end
        
        hold on;
        rho = tuning_ordered_unified_coord_base0(3,:,ca_id);
        polarplot([theta(:);theta(1)],[rho(:);rho(1)],'LineWidth',1,'Color',colors(i,:));
    end
    set(gca,'FontName','Arial','FontSize',15);
    
    labeld_axes = axes('Position',[0.45 0.45 0.02 0.02]);
    text(0.5,0.5,'d','Units','normalized','FontName','Arial','FontSize',25,'FontWeight','bold');
    axis off;
    
    
    % Temporal response panel
    temp_axes = axes('Position',[0.75 0.07 0.2 0.38]);
    
    x = 0:4/30:4;
    
    for i = 1:length(cells)
        temp_resp = temporal_response{cells(i)};
        if ~isempty(temp_resp)
            hold on;
            plot(x,temp_resp,'Color',colors(i,:));
        end
    end
    
    temp_resp = tempresp_avg(:,type_idx);
    
    hold on;
    plot(x,temp_resp,'Color',[119 136 153]/255,'LineWidth',4);
    temp_ylim = get(gca,'YLim');
    line([x(8) x(8)],temp_ylim,'LineStyle',':','LineWidth',1.5,'Color',[192 192 192]/255);
    line([x(16) x(16)],temp_ylim,'LineStyle',':','LineWidth',1.5,'Color',[192 192 192]/255);
    
    xlim([0 4]);
    set(gca,'FontName','Arial','FontSize',15);
    xlabel('Time (s)','FontName','Arial','FontSize',20);
    ylabel('\DeltaF/F (%)','Interpreter','tex','FontName','Arial','FontSize',20);
    
    labele_axes = axes('Position',[0.70 0.45 0.02 0.02]);
    text(0.5,0.5,'e','Units','normalized','FontName','Arial','FontSize',25,'FontWeight','bold');
    axis off;
    
end

