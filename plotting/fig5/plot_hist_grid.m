function plot_hist_grid(density_patch)
%% Plot histogram for 3D density map

h = histogram_grid_density(density_patch,0.06);
hist = h.Values;
box off;
% set(gcf,'Color','none','InvertHardcopy','off');
% set(gca,'Color','none')

xlim([0,0.9]);
yl = [0,13];
ylim(yl);

mu = mean(density_patch(:));
sd = std(density_patch(:));
hold on;
plot([mu,mu],[0,yl(2)],'--k');
hold on;
plot([mu-sd,mu-sd],[0,yl(2)],'-.k');
hold on;
plot([mu+sd,mu+sd],[0,yl(2)],'-.k');

end