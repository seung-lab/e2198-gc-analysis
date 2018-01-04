function cell_info_scatter_plot(cell_info,type_names, stat_types)

washold = ishold();
if ~washold
    cla;
end
hold on;

N = numel(type_names);
D = size(stat_types, 1); % D = numel(stat_types);
if D<2
    D = numel(stat_types);
    stat_types = stat_types.';
    if D>3 || D<2
        error('unsupported number of stats')
    end
end

colors=distinguishable_colors(N);

if strncmp(stat_types{1}, 'sac_corr', 8)
    [onsac, offsac] = get_sac_strat(cell_info);
end
if strcmp(stat_types{1}, 'corr') || strcmp(stat_types{1}, 'corr_unrml') || strcmp(stat_types{1}, 'corr_u') || strcmp(stat_types{1}, 'corr_uu')
    %{
    bcs = get_cell_info(cell_info, 3);  % class 3 = bc
    bctypes = cellstr(unique(char(bcs.type), 'rows'));
    bctypes = bctypes(2:end);  % bctypes{1} == ''
    %}
    bctypes = stat_types(:, 2);
    for d = 1:D
        strat_bc{d} = get_avg_strat(cell_info, bctypes{d});
    end
    %{
    [bc7] = get_avg_strat(cell_info, 'bc7');
    [bc89] = get_avg_strat(cell_info, 'bc8/9');
    %}
end

for k=1:N

    cells = get_cell_info(cell_info, type_names{k});

    %x = (cell_info_get_strat_property(cells, 'cum_sus'));
    %x(isinf(x)) = 0;
    %y = (cell_info_get_strat_property(cells, 'cum_trans'));
    %y(isinf(y)) = 0;
    x = cell_info_get_strat_property(cells, 'cum_sus') - cell_info_get_strat_property(cells, 'cum_trans');
    y = cell_info_get_strat_property(cells, 'cum_on') - cell_info_get_strat_property(cells, 'cum_off');
    y = cell_info_get_strat_property(cells, 'sus_on-trans_on');

    xyz = [];
    if strncmp(stat_types{1}, 'sac_corr', 8)
        x = cell_info_get_strat_property(cells, 'sac_corr', true, onsac);
        y = cell_info_get_strat_property(cells, 'sac_corr', true, offsac);
        xyz = [x y];
    elseif strcmp(stat_types{1}, 'corr') || strcmp(stat_types{1}, 'corr_unrml') || strcmp(stat_types{1}, 'corr_u')
        for d = 1:D
            xyz(:,d) = cell_info_get_strat_property(cells, stat_types{d}, true, strat_bc{d});
        end
        %{
        x = cell_info_get_strat_property(cells, 'sac_corr', true, bc7);
        y = cell_info_get_strat_property(cells, 'sac_corr', true, bc89);
        %}
    elseif strcmp(stat_types{1}, 'corr_uu')
        for d = 1:D
            xyz(:,d) = cell_info_get_strat_property(cells, stat_types{d}, false, strat_bc{d});
        end
    else
        for d = 1:D
            xyz(:,d) = cell_info_get_strat_property(cells, stat_types(d,:));
        end
        %{
        x = cell_info_get_strat_property(cells, stat_types{1});
        y = cell_info_get_strat_property(cells, stat_types{2});
        if D>2
            z = cell_info_get_strat_property(cells, stat_types{3});
        end
        %}
    end

    if D<3
        scatter(xyz(:,1),xyz(:,2), 'filled','MarkerFaceColor',colors(k,:),'MarkerEdgeColor',colors(k,:));
        for j = 1:size(xyz, 1)
            text(xyz(:,1),xyz(:,2), type_names{k});
        end
    else
        scatter3(xyz(:,1),xyz(:,2),xyz(:,3),'filled','MarkerFaceColor',colors(k,:),'MarkerEdgeColor',colors(k,:));
        for j = 1:size(xyz, 1)
            text(xyz(:,1),xyz(:,2),xyz(:,3), type_names{k});
        end
        grid on
    end

end
%axis equal;
legend(type_names);
xlabel(stat_types(1,:))
ylabel(stat_types(2,:))
if D>2
    zlabel(stat_types(3,:))
end
title(strjoin(stat_types, '   '), 'Interpreter', 'none')

if washold ~= ishold()
    hold();
end
