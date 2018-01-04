% Permutation test for grid density analysis
load('cell_types');
load('cell_info');
n_types = length(cell_types);

% Parameters
iter = 10000;
grid_size = 600;
bin_size = 600/grid_size^2;
bin_end = 1/bin_size;

%% ACTION NEEDED
% Mode
mode = 'orbit_rot_all'; 

% Output directory
out_dir = 'output_directory/';
%%

for t = 1:47
    type = cell_types{t}
    trunk_thr = get_trunkthr(type);
    
    density_type = density_grid(cell_info,type,trunk_thr,-0.2,grid_size,0);
    hist_type = histogram_grid_density(density_type,bin_size);
    hist_100 = hist_type.Values;
    close;
    
    density_perm = zeros(iter,numel(density_type));
    winner_perm = zeros(iter,numel(density_type));
    disp_perm = zeros(iter,1);
    mean_perm = zeros(iter,1);
    std_perm = zeros(iter,1);
    hist_perm = zeros(iter,bin_end);
    
    for i = 1:iter
        if rem(i,100) == 0 
            disp([type,'   ',num2str(i)]);
        end
        [density_patch,~,nodes,soma,winner_patch,winner_id] = density_grid(cell_info,type,trunk_thr,-0.2,grid_size,mode);
       
        if max(density_patch(:)) < bin_size
            continue;
        end
        h = histogram(density_patch,[0:bin_size:ceil(max(density_patch(:))/bin_size)*bin_size]);
        hist = h.Values;
        hist(bin_end) = 0; 
        hist_perm(i,:) = hist;
        density_perm(i,:) = reshape(density_patch,[numel(density_patch),1]);
        winner_perm(i,:) = reshape(winner_patch,[numel(winner_patch),1]);
        
        close;
        
        [disp_patch,mean_patch,std_patch] = dispersion(density_patch);
        disp_perm(i) = disp_patch;
        mean_perm(i) = mean_patch;
        std_perm(i) = std_patch;
    end
    
    
    x_range = find(sum(hist_perm,1)>0);
    x_end = x_range(end);
    if length(hist_100) > x_end
        x_end = length(hist_100);
    else
        hist_100(x_end) = 0;
    end
    hist_perm(:,x_end+1:end) = [];
    
    med_hist = median(hist_perm,1);
    ptile = prctile(hist_perm,[5,95]);
    neg = med_hist - ptile(1,:);
    pos = ptile(2,:) - med_hist;
    
    figure();
    bar(bin_size/2:bin_size:x_end*bin_size-bin_size/2,med_hist,1);
    hold on;
    errorbar(bin_size/2:bin_size:x_end*bin_size-bin_size/2,med_hist,neg,pos,...
        'k.','LineWidth',1.5);
    hold on;
    plot(bin_size/2:bin_size:x_end*bin_size-bin_size/2,hist_100,'Color','r','LineWidth',2);
    
    set(gca,'XTick',0:bin_size*length(hist_100)/4:x_end*bin_size);
    xlabel('Density','FontSize',15);
    ylabel('Count','FontSize',15);

    save(strcat(out_dir,type,'_perm'),'density_perm','disp_perm','mean_perm','std_perm');
    print(strcat(out_dir,type,'_perm_hist'),'-dpng','-r0');
    close;
    
end

