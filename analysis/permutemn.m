function perm_list = permutemn(m,n,unique)
i = 2;
perm_list = 1:m;

if unique
    while i~=n+2
        perm = randperm(m);
        if i == 1
            perm_list(i,:) = perm;
            i = i + 1;
        else
            overlap_row = sum(perm_list == repmat(perm,[size(perm_list,1),1]),2);
            if isempty(find(overlap_row==m))
                perm_list(i,:) = perm;
                i = i + 1;
            end
        end
    end
    
else
    while i~=n+2
        perm = randperm(m);
        perm_list(i,:) = perm;
        i = i + 1;
    end
end

perm_list = perm_list(2:end,:);

end