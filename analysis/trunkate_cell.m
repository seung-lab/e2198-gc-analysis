function n_trunkated = trunkate_cell(cell_id,in_thr,out_thr)
%% Trunkate cell

skel_dir = '../skeletons/';

load(strcat(skel_dir,'skel_',num2str(cell_id),'.mat'));

if ~exist('in_thr','var') || isempty(in_thr)   
    in_thr = 1.2;
end
if ~exist('out_thr','var') || isempty(out_thr)
    out_thr = -0.2;
end

% Convert IPL depth into indices
in_thr = 1063-100*in_thr*(175/34);
out_thr = 1063-100*out_thr*(175/34);

n = n(:,[3 2 1]);
n = n((n(:,1)>in_thr)&(n(:,1)<out_thr),:);
%     n = unique(n,'rows');
n_trunkated = n(:,2:3);
n_trunkated(:,2) = n_trunkated(:,2) * (23/16.5);

end