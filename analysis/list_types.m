function [gc, ac, bc, alltypes] = list_types(cell_info)

% cantknow=0, GC=1, AC=2, BC=3, uncertain=4, etc=5

for classId = 1:3
	gc = get_cell_info(cell_info, classId);
	types{classId} = cellstr(unique(char(gc.type), 'rows'));
end

[gc, ac, bc] = types{:};

alltypes = {};
for subgroup = types(:).'
	subgroup = subgroup{1};
	if isempty(subgroup{1})     % remove the empty type
    	subgroup = subgroup(2:end);
    end
    alltypes = [alltypes; subgroup];
end

assert( length(unique(alltypes)) == length(alltypes) );

end
