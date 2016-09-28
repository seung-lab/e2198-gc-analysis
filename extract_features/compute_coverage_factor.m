%% compute_coverage_factor
% Code to compute coverage factor

load cell_types;
load cell_hull;

% figure();
hull_all = [];
hull_area_sum = zeros(1,length(cell_types));
for t = 1:length(cell_types)
    cell_idx = find(strcmp({cell_info.type},cell_types{t}));
    
    hull_type = [];
    for i = 1:length(cell_idx)
        cell = cell_info(cell_idx(i)).cell_id;
        
        hull = cell_hull{cell};
        temp = hull_type;
        hull_type = [];
        if isempty(temp)
            hull_type = hull;
        else
            [hull_type(:,1) hull_type(:,2)] = polybool('union',hull(:,1),hull(:,2),temp(:,1),temp(:,2));
        end
        hull_area_sum(t) = hull_area_sum(t) + polyarea(hull(:,1),hull(:,2));
    end
   
    hold on;
%     plot(hull_type(:,1),hull_type(:,2));
    
    temp_all = hull_all;
    hull_all = [];
    if isempty(temp_all)
        hull_all = hull_type;
    else
        [hull_all(:,1) hull_all(:,2)] = polybool('union',hull_type(:,1),hull_type(:,2),temp_all(:,1),temp_all(:,2));
    end
    hull_list{t} = hull_type;
    hull_area(t) = polyarea(hull_type(:,1),hull_type(:,2));
end

hold on;
% plot(hull_all(:,1),hull_all(:,2),'LineWidth',3);
hull_total_area = polyarea(hull_all(:,1),hull_all(:,2));
hull_area(21) = 1.9988*10^7; % This one has a hole in the middle so manually modified the area

coverfactor = hull_area_sum./hull_area;

% save('coverfactor','coverfactor','hull_area_sum','hull_area','hull_total_area');


