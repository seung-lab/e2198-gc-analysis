function [s] = get_avg_strat(cell_info, typename, normalization_method)

stratname = 'strat_nrml';

sac = get_cell_info(cell_info, typename);
strat = cat(3, sac.(stratname));

if exist('normalization_method', 'var') && strcmp(normalization_method, 'vectorlength')
	norms = sqrt(sum(strat(:,2,:).^2,1));
	strat(:,2,:) = strat(:,2,:) ./ repmat(norms,size(strat,1),1,1);
end

strat = mean(strat, 3);
s=strat(:,2);
x=strat(:,1);
onsac = s;
