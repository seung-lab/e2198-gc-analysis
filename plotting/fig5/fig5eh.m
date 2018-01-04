load('./fig5_data');

density_type = density_type*0.0957*10^6/(66*66);
density_rand = density_rand*0.0957*10^6/(66*66);

%% 6sw density map
plot_density_map(density_type,nodes,soma);

plot_hist_grid(density_type);
set(gca,'XTickLabel',[],'YTickLabel',[],'YColor','none');

%% Randomized density map
plot_density_map(density_rand,nodes_rand,soma_rand);

plot_hist_grid(density_rand);
set(gca,'XTickLabel',[],'YTickLabel',[],'YColor','none');
