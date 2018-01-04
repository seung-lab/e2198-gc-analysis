load('cell_info');

clusters = containers.Map;

ds_types = {'37c','37d','37r','37v','7id','7ir','7iv','7o'};

for i = 1:length(ds_types)
    type = ds_types{i};
    
    idx = find(strcmp({cell_info.type},type));
    
    cell_list = zeros(length(idx),1);
    for j = 1:length(idx)
        cell_list(j) = cell_info(idx(j)).cell_id;
    end
    
    clusters(type) = cell_list;
end

cell_types_wods = {'1ws','1wt','1no','1ni','2an','2aw','2o','2i','3o','3i','4on','4ow','4i','5to','5ti','5so','5si','6sn','6sw','6t','8w','8n','9n','9w','51','25','85','63','73','72','27','81i','81o','82wo','82wi','82n','915','28','91'};
cell_types_wods = cell_types_order(cell_types_wods);


% a-1
create_scatter_a1;
close;

% a-2
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,cell_types_wods,{'sus-trans'},[-1,1],[3],[5],'Split a-2','Marginal - central',[],[],[0.05],'',[],[],[],[],[]);
ylim([0,160]);
clusters('a2') = cell_ids;
close;

split = kmeans_iter(stat,2,1000);
clusters('a3') = cell_ids(split==1);
clusters('a4') = cell_ids(split==2);

% a-3
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info_tran,clusters('a3'),{'on-off'},[-1,1],[1],[5],'Split a-3','Central inner - outer',[],[],[-0.5,0.5],'',[],[],[],[]);
% [cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('a3'),{'trans_off'},[],[1],[5],'Split a-3','outer intra-SAC',[],[],[0.098],'',[],[],[],[]);
close;

split = kmeans_iter(stat,3,1000);
clusters('b1') = cell_ids(split==1);
clusters('c1') = cell_ids(split==2);
clusters('d1') = cell_ids(split==3);

% a-4
% [cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('a4'),{'sus_on-sus_off'},[],[2],[5],'Split a-4','marginal inner - marginal outer',[],[],[-0.04],'',[0.05],[],[],[]);
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info_sus,clusters('a4'),{'on-off'},[-1,1],[2],[5],'Split a-4','Marginal inner - outer',[],[],[-0.15],'',[],[],[],[]);
ylim([0,110]);
close;

split = kmeans_iter(stat,2,1000);
clusters('e1') = cell_ids(split==1);
clusters('f1') = cell_ids(split==2);

% b-1
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('b1'),{'ptile'},[0.21 0.36],[0.5],[1],'Split b-1','10th percentile',[0.1],[],[0.28],'',[0.01],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('5to') = cell_ids(split==1);
clusters('b2') = cell_ids(split==2);

% b-2
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('b2'),{'somasize'},[700 3100],[0.5],[1],'Split b-2','Soma size (\mum^3)',[],[],[300000*0.165*0.165*0.23],'',[100],[],[]);
ylim([0,9]);
close;

split = kmeans_iter(stat,2,1000);
clusters('b3') = cell_ids(split==1);
clusters('4ow') = cell_ids(split==2);

% b-3
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('b3'),{'ptile-'},[0.05 0.11],[0.5],[1],'Split b-3','70th percentile - 15th percentile',[0.7],[0.15],[0.075],'',[0.004],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('4on') = cell_ids(split==1);
clusters('4i') = cell_ids(split==2);

% c-1
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('c1'),{'corr','OFF SAC'},[0 0.8],[1],[5],'Split c-1','Off SAC similarity',[],[],[0.42],'',[0.05],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('c2') = cell_ids(split==1);
clusters('63') = cell_ids(split==2);

% c-2
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('c2'),{'ptile'},[-0.05 0.45],[1],[5],'Split c-2','5th percentile',[0.05],[],[0.285],'',[0.02],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('c3') = cell_ids(split==1);
clusters('c4') = cell_ids(split==2);

% c-3
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('c3'),{'branch'},[0.055 0.093],[1],[5],'Split c-3','Arbor complexity (\mum^{-1})',[0.15],[],[0.0745],'',[0.002],[],[]);
close;

split = kmeans_iter(stat,2,1000);
clusters('51') = cell_ids(split==1);
clusters('5ti') = cell_ids(split==2);

% c-4
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('c4'),{'ptile'},[0.5 0.6],[0.7],[1],'Split c-4','80th percentile',[0.8],[],[0.542],'',[0.005],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('5so') = cell_ids(split==1);
clusters('5si') = cell_ids(split==2);

% d-1
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('d1'),{'ptile'},[0.642 0.85],[0.5],[1],'Split d-1','95th percentile',[0.95],[],[0.75],'',[0.02],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('d2') = cell_ids(split==1);
clusters('6t') = cell_ids(split==2);

% d-2
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('d2'),{'somasize'},[700 1700],[0.5],[1],'Split d-2','Soma size (\mum^3)',[],[],[180000*0.165*0.165*0.23],'',[100],[],[]);
close;

split = kmeans_iter(stat,2,1000);
clusters('6sn') = cell_ids(split==1);
clusters('6sw') = cell_ids(split==2);


%% outer marginal
% e-1
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info_sus,clusters('e1'),{'ptile'},[0 0.9],[3],[5],'Split e-1','Marginal 85th percentile',[0.85],[],[0.45],'',[],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('e2') = cell_ids(split==1);
clusters('e12') = cell_ids(split==2);

% e-2
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('e2'),{'somasize'},[450 2550],[2],[5],'Split e-2','Soma size (\mum^3)',[],[],[195000*0.165*0.165*0.23],'',[],[],[]);
close;

split = kmeans_iter(stat,2,1000);
clusters('e3') = cell_ids(split==1);
clusters('e11') = cell_ids(split==2);

% e-3
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('e3'),{'ptile'},[0.03 0.275],[1],[5],'Split e-3','50th percentile',[0.5],[],[0.11],'',[],[],[],[2]);
close;

split = kmeans_iter(stat,2,1000);
clusters('e5') = cell_ids(split==1);
clusters('e4') = cell_ids(split==2);

% e-4
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('e4'),{'branch'},[0.01 0.09],[1],[5],'Split e-4','Arbor complexity (\mum^{-1})',[],[],[0.05],'',[],[],[]);
close;

split = kmeans_iter(stat,2,1000);
clusters('e6') = cell_ids(split==1);
clusters('e10') = cell_ids(split==2);

% e-5
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('e5'),{'node_density'},[0.005 0.19],[0.5],[1],'Split e-5','Arbor density (\mum^{-1})',[],[],[0.07],'',[0.01],[],[],[]);
close;

split = kmeans_iter(stat,2,1000);
clusters('1ws') = cell_ids(split==1);
clusters('e7') = cell_ids(split==2);

% e-6
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info_off,clusters('e6'),{'ptile+'},[0.25 0.46],[1.2],[5],'Split e-6','Outer 80th percentile + outer 10th percentile',[0.8],[0.1],[0.35],'',[],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('e8') = cell_ids(split==1);
clusters('e9') = cell_ids(split==2);

% e-7
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('e7'),{'ptile'},[0.06 0.3],[0.5],[1],'Split e-7','85th percentile',[0.85],[],[0.15],'',[0.02],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('1no') = cell_ids(split==1);
clusters('1ni') = cell_ids(split==2);

% e-8
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('e8'),{'node_density'},[0.1 0.255],[0.5],[1],'Split e-8','Arbor density (\mum^{-1})',[],[],[0.152],'',[0.007],[],[]);
close;

split = kmeans_iter(stat,2,1000);
clusters('2aw') = cell_ids(split==1);
clusters('2i') = cell_ids(split==2);

% e-9
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info_off,clusters('e9'),{'ptile-'},[0.205 0.282],[0.5],[1],'Split e-9','Outer 90th percentile - outer 10th percentile',[0.9],[0.1],[0.238],'',[0.004],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('3o') = cell_ids(split==1);
clusters('3i') = cell_ids(split==2);

% e-10
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('e10'),{'ptile-'},[0.135 0.38],[0.5],[1],'Split e-10','90th percentile - 45th percentile',[0.9],[0.45],[0.3], '', [0.01],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('2an') = cell_ids(split==1);
clusters('25') = cell_ids(split==2);

% e-11
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('e11'),{'branch'},[0.008 0.06],[0.5],[1],'Split e-11','Arbor complexity (\mum^{-1})',[],[],[0.03],'',[],[],[]);
close;

split = kmeans_iter(stat,2,1000);
clusters('1wt') = cell_ids(split==1);
clusters('2o') = cell_ids(split==2);

% e-12
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('e12'),{'ptile'},[0.67 0.89],[0.7],[1],'Split e-12','95th percentile',[0.95],[],[0.78],'',[0.02],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('27') = cell_ids(split==1);
clusters('28') = cell_ids(split==2);


%% inner-marginal
% f-1
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('f1'),{'somasize'},[400 2900],[1],[5],'Split f-1','Soma size (\mum^3)',[],[],[290000*0.165*0.165*0.23],'',[],[],[]);
close;

split = kmeans_iter(stat,2,1000);
clusters('f2') = cell_ids(split==1);
clusters('8w') = cell_ids(split==2);

% f-2
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('f2'),{'ptile'},[0.71 0.97],[1],[5],'Split f-2','95th percentile',[0.95],[],[0.81],'',[0.015],[],[],[2]);
close;

split = kmeans_iter(stat,2,1000);
clusters('f6') = cell_ids(split==1);
clusters('f3') = cell_ids(split==2);

% f-3
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('f3'),{'trans_on'},[0 0.28],[1],[5],'Split f-3','Inner central',[],[],[0.15],'',[0.015],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('f4') = cell_ids(split==1);
clusters('85') = cell_ids(split==2);

% f-4
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('f4'),{'ptile'},[0 0.8],[1],[5],'Split f-4','5th percentile',[0.05],[],[0.45],'',[0.05],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('f5') = cell_ids(split==1);
clusters('f8') = cell_ids(split==2);

% f-5
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('f5'),{'trans'},[0 0.22],[1],[5],'Split f-5','Central',[],[],[0.09],'',[0.012],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('f7') = cell_ids(split==1);
clusters('f12') = cell_ids(split==2);

% f-6
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('f6'),{'ptile'},[0.67 0.735],[1],[5],'Split f-6','80th percentile',[0.8],[],[0.71],'',[0.005],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('f9') = cell_ids(split==1);
clusters('f10') = cell_ids(split==2);

% f-7
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info_on,clusters('f7'),{'ptile'},[0.7 0.85],[1],[5],'Split f-7','Inner 50th percentile',[0.5],[],[0.775],'',[0.01],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('f11') = cell_ids(split==1);
clusters('91') = cell_ids(split==2);

% f-8
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('f8'),{'ptile'},[0.77 0.85],[1],[5],'Split f-8','50th percentile',[0.5],[],[0.8],'',[0.005],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('8n') = cell_ids(split==1);
clusters('f13') = cell_ids(split==2);

% f-9
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('f9'),{'trans'},[0.11 0.41],[1],[5],'Split f-9','Central',[],[],[0.23],'',[0.03],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('72') = cell_ids(split==1);
clusters('73') = cell_ids(split==2);

% f-10
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('f10'),{'ptile'},[0.1 0.5],[1],[5],'Split f-10','25th percentile',[0.25],[],[0.35],'',[0.05],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('82wo') = cell_ids(split==1);
clusters('81o') = cell_ids(split==2);

% f-11
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('f11'),{'ptile'},[0.05 0.72],[1],[5],'Split f-11','25th percentile',[0.25],[],[0.3],'',[0.04],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('81i') = cell_ids(split==1);
clusters('82wi') = cell_ids(split==2);

% f-12 
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('f12'),{'sus-trans'},[0.58 0.81],[1],[5],'Split f-12','Marginal - central',[],[],[0.7],'',[0.01],[],[],[1]);
close;

split = kmeans_iter(stat,2,1000);
clusters('915') = cell_ids(split==1);
clusters('82n') = cell_ids(split==2);

% f-13
[cell_ids,stat,ctype,bin]=cell_info_hist(cell_info,clusters('f13'),{'area'},[8000 88000],[1],[5],'Split f-13','Dendritic field area (\mum^2)',[],[],[60000],'',[8000],[],[],[]);
close;

split = kmeans_iter(stat,2,1000);
clusters('9n') = cell_ids(split==1);
clusters('9w') = cell_ids(split==2);


%% BC
bc = cell_info_typedef_bc();

for i = 1:length(bc)
    clusters(bc(i).name) = bc(i).cells'; 
end

% Save file
save('clusters','clusters');