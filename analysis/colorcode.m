colors = [255 255 255; 0 0 143;0 0 255;0 111 255;0 239 255;64 255 191;143 255 111;191 255 64;255 255 0; 255 207 0;255 159 0; 255 80 0;255 16 0; 143 0 0];  
values = unique(dbf_cyl);

figure();
for i = 1:11
    for j = 1:11
        
         idx = find(values==dbf_cyl(i,j));
            
         hold on;
         scatter3(i*ones(1,100),j*ones(1,100),1:100, '.', 'MarkerEdgeColor',colors(idx,:)/255,'MarkerFaceColor',colors(idx,:)/255);
           
       
    end
end