function p = pvalue(list,value)
    res = 1000;
    
    percent_list = 0:0.001:100;
    prc = prctile(list,percent_list);
    [~, idx] = min(abs(prc-value));
    
    p = percent_list(idx)/100;
end