load('./cell_types.mat');
load('./numcell.mat');
load('./cell_info');
n_types = length(cell_types);

dir = './perm_orbit_all_original_600_10000/';
cell_types = cell_types(sort_idx);
iter = 10000;

disp_perm_mat = zeros(iter,n_types);
disp_type = zeros(2,n_types);
med_perm = zeros(1,n_types);
neg = zeros(1,n_types);
pos = zeros(1,n_types);
mean_shift = zeros(1,n_types);
for t = 1:n_types
    type = cell_types{t}

    load(strcat(dir,type,'_perm'));
    disp_perm_mat(:,t) = disp_perm;
%     load(strcat('./figures/grid_density/',type));
    trunk_thr = get_trunkthr(type);
    density_type = density_grid(cell_info,type,trunk_thr,-0.2,600,0);
    [disp_patch,mean_patch,std_patch] = coeffvar(density_type);
    
    disp_type(1,t) = disp_patch;
    disp_type(2,t) = pvalue(disp_perm_mat(:,t),disp_patch);
    
    med_perm(t) = median(disp_perm,1);
    prc = prctile(disp_perm,[1,99]);
    neg(t) = med_perm(t) - prc(1);
    pos(t) = prc(2) - med_perm(t);
    
    mean_shift(t) = (mean(mean_perm) - mean_patch)*100/mean_patch;
end

% Plot box plot
figure();
errorbar(1:n_types,med_perm,neg,pos,'Color','k','LineStyle','none','LineWidth',1,'CapSize',5);
hold on;
boxplot(disp_perm_mat,'Labels',cell_types,'Color','k','Symbol','','PlotStyle','traditional','Whisker',0);
hold on;
plot(1:n_types,disp_type(1,:),'.r','MarkerSize',30,'Linewidth',2)

% save(strcat(dir,'dispersion'),'disp_type');


fail_idx = find(disp_type(2,:)>0.01);
for i = 1:length(fail_idx)
    idx = fail_idx(i);
    fill([0.5+(idx-1),0.5+idx,0.5+idx,0.5+(idx-1)],[0,0,2.3,2.3],[192,192,192]/255,...
        'FaceAlpha',0.3,'EdgeColor','none');
end

ylim([0,2]);
set(gca,'XTickLabelRotation',90,'FontName','Arial','FontSize',20,'TickLength',[0.005,0.005]);
ylabel('Coefficient of variation','FontName','Arial','FontSize',20);
% legend({'Randomized 99/1 percentile','Type'},'FontName','Arial','FontSize',15)
% DataAspectRatio - [1,0.15,1]







