function cell_stat = cell_info_get_strat_property(cell_info, property_type_arg, offsac, onsac, center, use_normalized_strat, varargin)
% Note the cell_info input here is expected to be the already filtered subset of interest, rather than the full array of all cells.
%   	# TODO: perhaps I should just add the celltype argument and do filtering here.
% varargin: used for certain properties, such as 'corr' against a certain type of known strat profile.

self = @cell_info_get_strat_property;

N = length(cell_info);
cell_stat=zeros(N, 1);

if ~exist('use_normalized_strat', 'var')
	use_normalized_strat = true;
end

%strat=cell_info_bin_strat(cell_info,binstep);

if use_normalized_strat
	stratname = 'strat_nrml';
else
	stratname = 'strat_unrml';
end

if ~iscell(property_type_arg)
    property_type_arg = {property_type_arg};
end
property_type = property_type_arg{1};
property_arg = property_type_arg(2:end);

%TODO: this for loop can be made into an array operation instead

for k = 1:N
	if use_normalized_strat
		strat = cell_info(k).strat_nrml;
    	binwidth = abs(strat(2,1) - strat(1,1));
%     	strat(:,2) = strat(:,2) * binwidth;
	else
		strat = cell_info(k).strat_unrml;
	end
	s=strat(:,2);
    x=strat(:,1);
    
    cell_info_elem = cell_info(k);

	switch property_type
	case {'sus_strat_vol', 'cum_sus', 'sus'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'sus_on', offsac, onsac, center) + cell_info_get_strat_property(cell_info(k), 'sus_off', offsac, onsac, center);

	case {'trans_strat_vol', 'cum_trans', 'trans'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'trans_on', offsac, onsac, center) + cell_info_get_strat_property(cell_info(k), 'trans_off', offsac, onsac, center);
        % cell_stat(k) = cell_stat(k) + sum(s(x==28 | x==62));

	case {'cum_on', 'on'}
        cell_stat(k) = sum(s(x>center));

	case {'cum_off', 'off'}
        cell_stat(k) = sum(s(x<center));

	case {'cum_trans_on', 'trans_on'}
        cell_stat(k) = sum(s(x>center & x<onsac));

	case {'cum_sus_on', 'sus_on'}
        cell_stat(k) = sum(s(x>onsac));

	case {'cum_trans_off', 'trans_off'}
        cell_stat(k) = sum(s(x<center & x>offsac));

	case {'cum_sus_off', 'sus_off'}
        cell_stat(k) = sum(s(x<offsac));

   	case {'on-off'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'on', offsac, onsac, center) - cell_info_get_strat_property(cell_info(k), 'off', offsac, onsac, center);

   	case {'sus-trans'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'sus', offsac, onsac, center) - cell_info_get_strat_property(cell_info(k), 'trans', offsac, onsac, center);

   	case {'sus_on-trans_on'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'sus_on', offsac, onsac, center) - cell_info_get_strat_property(cell_info(k), 'trans_on', offsac, onsac, center);

   	case {'sus_off-trans_off'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'sus_off', offsac, onsac, center) - cell_info_get_strat_property(cell_info(k), 'trans_off', offsac, onsac, center);

   	case {'sus_on-sus_off'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'sus_on', offsac, onsac, center) - cell_info_get_strat_property(cell_info(k), 'sus_off', offsac, onsac, center);

    case {'trans_on-trans_off'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'trans_on', offsac, onsac, center) - cell_info_get_strat_property(cell_info(k), 'trans_off', offsac, onsac, center);

    case {'sus_off-trans_on'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'sus_off', offsac, onsac, center) - cell_info_get_strat_property(cell_info(k), 'trans_on', offsac, onsac, center);

	case {'corr', 'sac_corr'}
		sac = varargin{1};
		cell_stat(k) = s.'*sac / sqrt(sum(s.^2) * sum(sac.^2));

	case {'corr_unrml', 'corr_u', 'corr_uu'}
		sac = varargin{1};
		cell_stat(k) = s.'*sac;% / sqrt(sum(s.^2) * sum(sac.^2));

    case {'corr_sac_sac', 'sac_sac_2'}
        sac = varargin{1};
        cell_stat(k) = self(cell_info(k), 'sac_corr', use_normalized_strat, sac{1}).^2 ...
                     + self(cell_info(k), 'sac_corr', use_normalized_strat, sac{2}).^2;

    case {'prcntile', 'ptile'}
        cell_stat(k) = get_percentile([x s], property_arg{1});

    case {'asym_2an_p', 'asym2_p'}
        cell_stat(k) = log(norm(cell_info_elem.asymm_2an_prj));
        if (cell_info_elem.is_cutoff)
            cell_stat(k) = NaN;
        end

    case {'density'}
        cell_stat(k) = cell_info_elem.area_projection / cell_info_elem.area_hull;

    case {'size'}
        cell_stat(k) = cell_info_elem.area_hull;

    case {'maxdiameter'}
        cell_stat(k) = cell_info_elem.max_diameter;


    % case {'sus/tran'}
    % 	cell_stat(k) = cell_info_get_strat_property(cell_info())
    %     if at
    %         cell_stat(k) = log(as/at);
    %     else
    %         cell_stat(k) = 20;
    %     end

    otherwise
        error('not recognized stat name') 

    end
end
