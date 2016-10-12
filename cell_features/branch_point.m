%% branch_point
% Function to acquire branch nodes and compute the number of branch points

function [branch_node, num_branch] = branch_point(n,e)
    num_nodes = size(n,1);
    num_branch = 0;
    branch_node = [0, 0, 0];
    
    for i = 1:num_nodes
        if length(find(e==i)) >= 3
            num_branch = num_branch + 1;
            branch_node(num_branch,:) = n(i,:);
        end
    end
    
end