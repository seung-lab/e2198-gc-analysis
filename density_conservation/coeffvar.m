function [disp_patch,mean_patch,std_patch] = coeffvar(density_patch)
%% Measure the dispersion value of density patch
% disp = std/mean

density = reshape(density_patch,[1,numel(density_patch)]);

mean_patch = mean(density);
std_patch = std(density);
disp_patch = std_patch/mean_patch;
end