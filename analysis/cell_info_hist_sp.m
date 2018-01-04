function [cells,cell_stat,ctype,bin]=cell_info_hist_sp(cell_info,type_names, stat_type_arg, x_lim, e1, e2, graph_title, graph_xlb, varargin) %,  p, pminus, printfigure

% stat_type: prcntile, prcntileDiff, peakWidth

nvarargin = length(varargin);
optargs = {0.2, [], [], '', [], Inf, false, false, [], [], []};
optargs(1:nvarargin) = varargin;
[p, pminus, divisions, printfigure, binsize, cutoff, printcells, showstrat, offsac, onsac, center] = optargs{:};
if isempty(cutoff)
    cutoff = Inf;
end

if ~exist('offsac','var') || isempty(offsac)
    offsac = 28;
end
if ~exist('onsac','var') || isempty(onsac)
    onsac = 62;
end
if ~exist('center','var') || isempty(center)
    center = 47;
end

if ~iscell(stat_type_arg)
    stat_type_arg = {stat_type_arg};
end
stat_type = stat_type_arg{1};
stat_arg = stat_type_arg(2:end);

titletext2 = '';
switch stat_type
case {'corr', 'corr_unrml', 'corr_u', 'corr_u-corr_u', 'corr_uu'}
    bctype = stat_type_arg{2};
    %{
    switch bctype  % cut off cell body for SACs
        case 'ON SAC'
            [onsac, offsac] = get_sac_strat(cell_info);
            corr_against = onsac;
        case 'OFF SAC'
            [onsac, offsac] = get_sac_strat(cell_info);
            corr_against = offsac;
        otherwise
            corr_against = get_avg_strat(cell_info, bctype);
    end
    %}
    corr_against = get_avg_strat(cell_info, bctype);
    %titletext2 = strjoin(stat_arg);
    titletext2 = strsplit(evalc('disp(cell2table(stat_arg))'), '\n');
    titletext2 = titletext2{end-1};
case {'sac_corr'}  % cut off cell body for SACs
    [onsac, offsac] = get_sac_strat(cell_info);
    switch stat_type_arg{2}
        case 'ON SAC'
            corr_against = onsac;
        case 'OFF SAC'
            corr_against = offsac;
        otherwise
            error('not recognized SAC type')
    end
    titletext2 = strjoin(stat_arg);
case {'sac_sac_2'}
    [onsac, offsac] = get_sac_strat(cell_info);
    corr_against = {onsac, offsac};
    titletext2 = ': sac^2 + sac^2';
otherwise
%TODO: REFACTOR: move stat specific arguments (p, pminus) into this argument
    titletext2 = sprintf('%g %g', p, pminus);
    % do nothing
end

binstep=0;

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

if N==0
    error('no cells found')
end

strat=cell_info_bin_strat(cell_info,binstep);
cell_stat=zeros(size(cells));
ctype=cell(size(cells));

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
    %%{
    s=s(x<cutoff);
    x=x(x<cutoff);
    %}

    switch stat_type
    case {'prcntile', 'ptile', 'ptileDiff', 'percntileDiff', 'ptile-'}
        %[locm,locl,locr]=cell_info_get_cell_height_from_prcntile([x s],p);
        %cell_stat(j)=locr;  %20
        %cell_stat(j)=locl; %80
        %cell_stat(j)=locl-locr;
        %tmp = sortrows([x s], 1);   % put into ascending order
        if ~isnan(s(1))
            cell_stat(j) = get_percentile([x s],p);
            if ~isempty(pminus) || strcmp(stat_type, 'percntileDiff') || strcmp(stat_type, 'ptileDiff')	...
                        || strcmp(stat_type, 'ptile-')  % compute height between p and pminus
                %[~,~,minus]=cell_info_get_cell_height_from_prcntile([x s],pminus);
                minus = get_percentile([x s],pminus);
                cell_stat(j)=cell_stat(j) - minus;
                if cell_stat(j) < 0
                    cell_stat(j) = -cell_stat(j);
                end
            end
        
        else
            cell_stat(j) = 0;
        end
    case {'ptile+'}
        cell_stat(j) = get_percentile([x s],p);
        if ~isempty(pminus) || strcmp(stat_type, 'percntileDiff') || strcmp(stat_type, 'ptileDiff') % compute height between p and pminus
            %[~,~,minus]=cell_info_get_cell_height_from_prcntile([x s],pminus);
            minus = get_percentile([x s],pminus);
            cell_stat(j)=cell_stat(j) + minus;
        end
    case {'part-'}        
        s_minus = s;
        s_minus(1:475) = 0;     
        if sum(s_minus) == 0
            minus = 0;
        else
            s_minus = s_minus/sum(s_minus);
            minus = get_percentile([x s_minus],pminus);
        end
        
        s(300:end) = 0;
        if sum(s) == 0
            perc = 0;
        else
            s = s/sum(s);
            perc = get_percentile([x s],p);
        end
        
        cell_stat(j) = perc - minus;
        
   case {'part+'}        
        s_minus = s;
        s_minus(1:475) = 0;     
        if sum(s_minus) == 0
            minus = 0;
        else
            s_minus = s_minus/sum(s_minus);
            minus = get_percentile([x s_minus],pmineus);
        end
        
        s(300:end) = 0;
        if sum(s) == 0
            perc = 0;
        else
            s = s/sum(s);
            perc = get_percentile([x s],p);
        end
        
        cell_stat(j) = perc + minus;

    case {'peak_no_split', 'PNS'}
        [x1, x2]=cell_info_get_cell_height_from_peak_not_allow_split([x s],p);
        cell_stat(j) = x2 - x1;

    case {'peaks', 'PS'}
        [x1, x2]=cell_info_get_cell_height_from_peak_allow_split([x s],p);
        %cell_stat(j) = x2 - x1;
        xx1=min([x1; x2]);
        xx2=max([x1; x2]);
        cell_stat(j) = xx2 - xx1;

    case {'asymmetry', 'asym'}
        cell_stat(j) = norm(cell_info_elem.asymm_index);

    case {'asym_2an', 'asym2'}
        cell_stat(j) = log(norm(cell_info_elem.asymm_2an));
        if (cell_info_elem.is_cutoff)
            cell_stat(j) = NaN;
        end
    case {'asym_2an_p', 'asym2_p'}
        cell_stat(j) = log(norm(cell_info_elem.asymm_2an_prj));
        if (cell_info_elem.is_cutoff)
            cell_stat(j) = NaN;
        end

    case {'density'}
        cell_stat(j) = cell_info_elem.area_projection / cell_info_elem.area_hull;

    case {'size'}
        cell_stat(j) = cell_info_elem.area_hull;

    case {'maxdiameter'}
        cell_stat(j) = cell_info_elem.max_diameter;

    case {'radius'}
%         cell_stat(j) = cell_info_elem.radii_hull;
        cell_stat(j) = cell_info_elem.radius;
    case {'area'}
%         cell_stat(j) = cell_info_elem.area_projection;
        cell_stat(j) = cell_info_elem.hull_area;
        
    case {'somasize'}
        cell_stat(j) = cell_info_elem.soma_size;

    case {'branch'}
        cell_stat(j) = cell_info_elem.branch;     
        
    case {'node_density'}
        cell_stat(j) = cell_info_elem.node_density;
        
    case {'value'}
        cell_stat(j) = mean( s(find(x>p & x<p+1)) );


    case {'s:t.stratratio'}
        as = sum(s(x<28)) + sum(s(x>62));
        at = sum(s(x>28 & x<62));
        if at
            cell_stat(j) = log(as/at);
        else
            cell_stat(j) = 20;
        end

    case {'corr', 'corr_unrml', 'corr_u', 'sac_sac_2'}
        cell_stat(j) = cell_info_get_strat_property(cell_info_elem, stat_type, offsac, onsac, center, true, corr_against);

    case {'corr_uu'}
        cell_stat(j) = cell_info_get_strat_property(cell_info_elem, stat_type, offsac, onsac, center, false, corr_against);

    case {'corr-corr'}
        corr_against1 = get_avg_strat(cell_info, 'bc2');
        corr_against2 = get_avg_strat(cell_info, 'bc4');
        cell_stat(j) = cell_info_get_strat_property(cell_info_elem, 'corr', offsac, onsac, center, true, corr_against1) ...
            - cell_info_get_strat_property(cell_info_elem, 'corr', offsac, onsac, center, true, corr_against2);
    
    otherwise
        %try
            cell_stat(j) = cell_info_get_strat_property(cell_info_elem, stat_type, offsac, onsac, center);
        %catch
        %error('not recognized stat name') 
        %end

    end
end
%cell_stat


autoxlim = isempty(x_lim);
if ~autoxlim && ( x_lim(2) > 20 * max(cell_stat) || x_lim(2) < min(cell_stat) )
    autoxlim = true;
    warning('Auto resetting improper limits');
end
if autoxlim
    x_lim(1) = min(cell_stat);
    x_lim(2) = max(cell_stat);
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

ub = x_lim(2);
if mod(range, binsize) ~= 0
    % fix missing last bin % TODO: maybe I should use the native support for automatic binning in the histcounts().
    ub = ub + binsize;
end
binranges=x_lim(1):binsize:ub;
%[cnts,bin]=histc(cell_stat,binranges);
%bar(binranges,cnts,'histc');
cnts = histcounts(cell_stat,binranges);
bin = discretize(cell_stat,binranges);
%binranges
%cnts = histcounts(cell_stat(~isnan(cell_stat)), binranges);
%bin = discretize(cell_stat(~isnan(cell_stat)), binranges);
figure;
edge1 = e1;
edge2 = e2;
if isempty(e2)
    bargraph1 = bar(binranges(1:edge1),cnts(1:edge1),'histc');
    hold on;
    bargraph2 = bar(binranges(edge1:end-1),cnts(edge1:end),'histc');

%     set(bargraph1,'FaceColor',[255 0 0]/255);
%     set(bargraph2,'FaceColor',[0 0 255]/255);
    set(bargraph1,'FaceColor',[213 82 115]/255);
    set(bargraph2,'FaceColor',[0 97 187]/255);
else
    bargraph1 = bar(binranges(1:edge1),cnts(1:edge1),'histc');
    hold on;
    bargraph2 = bar(binranges(edge1:end-1),cnts(edge1:end),'histc');
    hold on;
    bargraph3 = bar(binranges(edge2:end-1),cnts(edge2:end),'histc');
    
    set(bargraph1,'FaceColor',[213 82 115]/255);
    set(bargraph2,'FaceColor',[115 115 46]/255);
    set(bargraph3,'FaceColor',[0 97 187]/255);
%     set(bargraph1,'FaceColor',[255 0 0]/255);
%     set(bargraph2,'FaceColor',[0 255 0]/255);
%     set(bargraph3,'FaceColor',[0 0 255]/255);
end
ylim([0 max(cnts)+1])

strat_mean = [];
%%{
% order by stat within each type
for type = type_names(:).'
    type = type{1};
%     idx = find(strncmp(ctype, type, length(type)));
    idx = find(strcmp(ctype, type));
    [B,I] = sort(cell_stat(idx));
    cells(idx)     = cells(idx(I));
    cell_stat(idx) = cell_stat(idx(I));
    bin(idx)       = bin(idx(I));
    if showstrat
        tmp = cat(3, strat{cells(idx)});
        strat_mean(:, end+1) = squeeze(mean(tmp(:,2,:), 3));
    end
end
%}



ax=gca();
if ~autoxlim
    ax.XLim=x_lim;
end


if strcmp(stat_type,'ptile') || strcmp(stat_type,'ptile-') || strcmp(stat_type,'ptile+')
    stat_type = 'IPL depth';
end

xlabel(graph_xlb,'FontName','Arial','FontSize',45);
set(gca,'XTick',[],'YTick',[],'FontSize',25,'FontName','Arial');


if ~isempty(divisions) && max(divisions) < x_lim(2) && min(divisions) > x_lim(1)
    n = length(divisions);
    divisions = divisions(:).'; %make row vec
    %line([divisions; divisions], repmat([0; 10], 1, n), 'Color', [0 .6 .7])
    line([divisions; divisions], repmat([0; ax.YLim(2)], 1, n), 'Color', [0 0 0],'LineWidth',1)
%     xtick = ax.XTick;
%     ax.XTick = unique([xtick divisions]);
end

if showstrat
    new_axis = axes('position',[0.6 0.7 0.4 0.3]);
    x=strat{cells(1)}(:,1);
    plot(x, strat_mean);
    legends = type_names;
    switch stat_type
    case {'corr', 'corr_unrml', 'corr_u', 'corr_u-corr_u'}
        %strat_mean(:, end+1) = corr_against;
        hold on
        plot(x, corr_against, ':', 'LineWidth', 2)
        legends(end+1) = {bctype};
    otherwise
        % do nothing
    end
    set(new_axis,'color','none');
    set(new_axis,'visible','off');
    legend(legends);
end



if ~isempty(printfigure) && ischar(printfigure)
    v = strsplit(version());
    v = v{1};
    idx = find(v=='.');
    v(idx(1:2)) = [];
    if str2num(v) >= 900.3413 % >= 9.0.0.341360 (R2016a)
        h = gcf();
        pos = h.Position;
        pos(3:4) = pos(3:4) * 2;
        h.Position = pos;
    end

    print(gcf, '-r300', printfigure, '-dpng');
    %print(gcf, '-r300', printfigure, '-depsc');
    print(gcf, '-r300', printfigure, '-dsvg');
    fprintf('saved file %s \n', printfigure)
    %close(summary_fig_h);
end

if printcells
    tab = repmat(sprintf('\t'), N, 1);
    [num2str(cells) tab num2str(cell_stat) tab char(ctype)]
    %{
    for j=1:N
        fprintf(stat)
    end
    %}
end

end