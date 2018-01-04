function tmp_print_division(ctype, clocation, dip1, dip2)

idx=clocation;
if ~exist('dip2', 'var')  %print saddle
	div = dip1;
	C = intersect(ctype(idx > div), ctype(idx <=div))
	for c = C.'
		below = idx <= div & strcmp(ctype, c);
		above = idx > div & strcmp(ctype, c);
		fprintf('%s \t %d \t %d \n', c{1}, sum(below), sum(above))
	end
else
	if dip1 > dip2
		[dip1, dip2] = deal(dip2, dip1);
	end
	indices = idx > dip1 & idx <= dip2;
	C = unique(ctype(indices))
	for c = C.'
		inside = indices & strcmp(ctype, c);
		below = idx <= dip1 & strcmp(ctype, c);
		above = idx > dip2 & strcmp(ctype, c);
		%outside = ~indices & strcmp(ctype, c);
		if sum(below) || sum(above)
			%fprintf('%s \t %d \t\t %d \t %d \n', c{1}, sum(inside), sum(below), sum(above))
			fprintf('%s \t %d \t\t %d \t %d \t %d \n', c{1}, sum(inside), sum(below), sum(inside), sum(above))
		else
			fprintf('%s \t %d \t\t %s \t %d \t %s \n', c{1}, sum(inside), ' ', sum(inside), ' ')
			%fprintf('%s \t %d \n', c{1}, sum(inside))
		end
	end
end