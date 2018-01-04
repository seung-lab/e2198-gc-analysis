function plot_10barstrat(types)

load cell_info;
load cell_types;

if isempty(types)
    figure();
    
    for t = 1:length(cell_types)
        type_idx = find(strcmp(cell_types,cell_types{t}));
        x = cell_info([cell_info.cell_id]==99900+type_idx).strat_nrml(:,1);
        strat = cell_info([cell_info.cell_id]==99900+type_idx).strat_nrml(:,2);
        
        for i = 1:10
            bar_strat(i) = sum(strat(x>=(i-1)*10 & x<i*10))*5.1471;
        end
        
        subplot(7,7,t);
        bar(bar_strat);
        ylabel(cell_types{t});
        xlim([0.5,10.5]);
        
    end
    
else
    
    for t = 1:length(types)
        if isstr(types{t})
            type_idx = find(strcmp(cell_types,types{t}));
            x = cell_info([cell_info.cell_id]==99900+type_idx).strat_nrml(:,1);
            strat = cell_info([cell_info.cell_id]==99900+type_idx).strat_nrml(:,2);
        else
            x = cell_info([cell_info.cell_id]==types{t}).strat_nrml(:,1);
            strat = cell_info([cell_info.cell_id]==types{t}).strat_nrml(:,2);
        end
        
        for i = 1:10
            bar_strat(i) = sum(strat(x>=(i-1)*10 & x<i*10))*5.1471;
        end
        
        figure();
        bar(bar_strat);
        xlabel('Layer');
        ylabel('Skeleton Density');
        xlim([0.5,10.5]);
        l = types{t};
        legend(l);
    end
end

end

