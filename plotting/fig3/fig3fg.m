%% Figure 2 panel f,g
cell_list = [17160,20181,17077,20166,20132];

figure();
% Plot GC inner-outer histogram
histogram_axes = axes('Position',[0.3,0.25,0.4,0.63]);

load cell_info;
load clusters;

offsac = 28; 
onsac = 62;
center = 47;

stat_type_arg = {'sus-trans'};
stat_type = stat_type_arg{1};
stat_arg = stat_type_arg(2:end);

type_names = clusters('a2');
binstep=0;
x_lim = [-1,1];
binsize = [];
cutoff = Inf;

cells=[];
if isnumeric(type_names)
    cells = type_names;
else
    for j=1:numel(type_names)
        %     idx=strncmp({cell_info.type},type_names{j}, length(type_names{j}));
        idx=strcmp({cell_info.type},type_names{j});
        if isempty(find(idx))
            error(sprintf('Unrecognized type "%s"', type_names{j}));
        end
        cells=[cells; [cell_info(idx).cell_id]'];
    end
end
N=numel(cells);

strat=cell_info_bin_strat(cell_info,binstep);
cell_stat=zeros(size(cells));

jj=1;
type_names = {};
for j=1:N
    cell_info_elem = get_cell_info(cell_info, cells(j));
        

    ctype{j}=cell_info([cell_info.cell_id]==cells(j)).type;
    if sum(strcmp(type_names,ctype{j}))==0
        type_names{jj} = ctype{j};
        jj = jj + 1;
    end
    s=strat{cells(j)}(:,2);
    x=strat{cells(j)}(:,1)/100;
    s=s(x<cutoff);
    x=x(x<cutoff);
    
    cell_stat(j) = cell_info_get_strat_property(cell_info_elem, stat_type, offsac, onsac, center);

end


range = x_lim(2)-x_lim(1);
if ~isempty(binsize) && ( binsize >= range || binsize < range/1000 )
    binsize = [];
    warning('Auto resetting improper binsize');
end
if isempty(binsize)
    if range < 100 && range > 5
        binsize = 1;
    else
        binsize = range / 30;
    end
end


lb = x_lim(1);
ub = x_lim(2);
if mod(range, binsize) ~= 0
    % fix missing last bin % TODO: maybe I should use the native support for automatic binning in the histcounts().
    ub = ub + binsize;
    lb = lb - binsize;
end
binranges=x_lim(1):binsize:ub;
cnts = histcounts(cell_stat,binranges);
bin = discretize(cell_stat,binranges);
bargraph = bar(binranges(1:end-1),cnts,'histc');
set(bargraph,'FaceColor',[119 136 153]/255);
set(gca,'FontSize',15,'FontName','Arial','XAxisLocation','top');
xlabel('marginal - central arbor','FontSize',25,'FontName','Arial');
ylabel('Number of GC','FontSize',25,'FontName','Arial');
xlim([lb,ub]);

% Plot cell stats
stat = zeros(length(cell_list),1);
for i = 1:length(cell_list)
    cell = cell_list(i);
    stat(i) = cell_stat(cells == cell);
%     hold on;
%     plot(stat(i),7,'.k','MarkerSize',20);
end

% Separation index illustration
[cluster,centroid] = kmeans_iter(cell_stat,2,1000);

hold on;
plot(centroid,[7,7],'xk','MarkerSize',15,'LineWidth',2);
% hold on;
% plot([centroid(1),centroid(1)],[0,35],'--k');
% hold on;
% plot([centroid(2),centroid(2)],[0,35],'--k');


% Plot GC inner/outer diagram
diagram_axes = axes('Position',[0.3,0.01,0.4,0.25]);

cell_info_filename = 'cell_info_clustering.20160623warpedSomaCorrection.mat';

position_list = (stat + 1).*1500;
position_list(4) = position_list(4)-100;
compress_ratio = [4,5,5,9,7];
% compress_ratio = [7,7,7,7,7];

boundary_list = [0.62,0.28];
boundary_list = -(boundary_list*100*175/34 - 1063.12);

if length(boundary_list) == 1
    color_arr = [106,189,69;...
        33,78,37];
else
    color_arr = [33,78,37;...
        170,18,20;...
        33,78,37];
end
color_arr = color_arr/255;


for i = 1:length(cell_list)
    cell = cell_list(i);
    [skel_nodes, cell_soma_coords] = compress_cell(cell,compress_ratio(i),70);
    
    skel_nodes(:,2) = skel_nodes(:,2) + position_list(i) - cell_soma_coords(2);
    cell_soma_coords(2) = position_list(i);
    
    for j = 1:length(boundary_list)+1
        if j == 1
            valid = skel_nodes(:,1) < boundary_list(1);
        elseif j == length(boundary_list)+1
            valid = skel_nodes(:,1) > boundary_list(end);
        else
            valid = skel_nodes(:,1) < boundary_list(j) & skel_nodes(:,1) > boundary_list(j-1);
        end
        
        skel_nodes_plot = skel_nodes(valid,:);

        hold on;
        plot(skel_nodes_plot(:,2), skel_nodes_plot(:,1), '.', ...
            'MarkerSize',1,'Color',color_arr(j,:));
    end
    
    hold on;
    plot(cell_soma_coords(:,2), 570, '.', ...
        'MarkerSize', 50, 'Color','k')
end

axis equal;
ylim([510,1146]); 
xlim([0,3000]);

x_lim = get(gca,'XLim');
for i = 1:length(boundary_list)
    boundary = boundary_list(i);
    hold on;
    plot(x_lim,[boundary,boundary],'--k','LineWidth',1.5);
end

set(gca,'XTick',[],'XColor',[1,1,1],'YTick',[548.41,1063.12],...
    'YTickLabel',[1,0],'FontName','Arial','FontSize',15);
xlabel('','FontName','Arial','FontSize',25);
ylabel('','FontName','Arial','FontSize',25);
