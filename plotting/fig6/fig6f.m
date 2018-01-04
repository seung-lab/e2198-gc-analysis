%% Soma size
load cell_info;
load clusters;

offsac = 28; 
onsac = 62;
center = 47;

stat_type_arg = {'somasize'};
stat_type = stat_type_arg{1};
stat_arg = stat_type_arg(2:end);

type_names = clusters('e1');

binstep=0;
x_lim = [];
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
    
    cell_stat(j) = cell_info_elem.soma_size;
end

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
figure;
bargraph = bar(binranges(1:end-1),cnts,'histc');
set(bargraph,'FaceColor',[119 136 153]/255);
set(gca,'FontSize',25,'FontName','Arial');
xlabel('inner - outer','FontSize',25,'FontName','Arial');
ylabel('Number of BC','FontSize',25,'FontName','Arial');
xlim([lb,ub]);

    
