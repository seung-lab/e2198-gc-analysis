%% create_table

load cell_info;
load tempresp_avg;
load cell_types;
gc = cell_info_typedef_gc();

figure();
h = 2;
w = 1.2;
xmax = 5.5;
ymax = 98;
xlim([0 16]);
ylim([0 ymax]);

% rectangle('Position',[1 ymax-29 4 26],'FaceColor',[0 0.4 0],'LineStyle','none');
% rectangle('Position',[1 ymax-37 4 8],'FaceColor',[1,1,1],'LineStyle','none');
% rectangle('Position',[1 ymax-43 4 6],'FaceColor',[0.6 0 0],'LineStyle','none');
% rectangle('Position',[1 ymax-55 4 12],'FaceColor',[0.8 0 0],'LineStyle','none');
% rectangle('Position',[1 ymax-61 4 6],'FaceColor',[1 0 0],'LineStyle','none');
% rectangle('Position',[1 ymax-69 4 8],'FaceColor',[1 1 1],'LineStyle','none');
% rectangle('Position',[1 ymax-95 4 26],'FaceColor',[0 1 0],'LineStyle','none');
rectangle('Position',[1 0 1+2*w ymax],'FaceColor',[1 1 1],'LineStyle','none');

line([1.3 xmax],[ymax-3 ymax-3],'LineWidth',1,'Color',[0 0 0]);

font_size = 8;
text(1.9,ymax-1.5,'Cluster','FontSize',font_size,'HorizontalAlignment','right','FontName','Arial');
% text(2+w/2,ymax+1,'Strat','FontSize',11,'HorizontalAlignment','center');
text(2+w/2,ymax-1.5,'Stratification','FontSize',font_size,'HorizontalAlignment','center','FontName','Arial');
text(2+3*w/2,ymax-1.5,'Visual response','FontSize',font_size,'HorizontalAlignment','center','FontName','Arial');
text(4.5,ymax-1.5,'Other name','FontSize',font_size,'HorizontalAlignment','left','FontName','Arial');

for i = 1:length(cell_types)
    type = cell_types{i};
    other_name = gc(strcmp({gc.name},type)).annotation;
    text(1.9,ymax-(2+2*i),type,'FontSize',font_size,'HorizontalAlignment','right','FontName','Arial');
    text(4.5,ymax-(2+2*i),'name','FontSize',font_size,'HorizontalAlignment','left','FontName','Arial');
    
    type_id = 99900+i;

    % plot strat
    bin_res = 4;
    x = cell_info([cell_info.cell_id]==type_id).strat_nrml(:,1)/100;
    strat = cell_info([cell_info.cell_id]==type_id).strat_nrml(:,2);
    
    x_idx = 1:bin_res:723;
    x = x(x_idx);
    strat_bin = strat(x_idx);
    for b = 1:bin_res-1
        idx_next = x_idx+b;
        if idx_next(end) > 723
            idx_next(end) = [];
            temp = strat(idx_next);
            temp(end+1) = 0;
        else
            temp = strat(idx_next);
        end
        
        strat_bin = strat_bin + temp;
    end
        
    strat_bin = strat_bin(x>-0.1 & x<1.1);
    x = x(x>-0.1 & x<1.1);
    x1_idx = find(x<=0.28);
    x2_idx = find(x>=0.28 & x<=0.475);
    x3_idx = find(x>=0.28 & x<=0.62);
    x4_idx = find(x>=0.47 & x<=0.62);
    x5_idx = find(x>=0.62);
    x = w*(x+0.1)/1.2;
    strat_bin = h*strat_bin/(1.3*max(strat_bin));

    hold on;
    plot(x+2,strat_bin+ymax-(2.98+2*i),'LineWidth',0.7,'Color',[0 0 0]);
    if i<=13
        area(x(x1_idx)+2,strat_bin(x1_idx)+ymax-(2.98+2*i),'FaceColor',[0,0.3,0],'EdgeAlpha',0);
%         area(x([x5_idx;x4_idx;x2_idx])+2,strat_bin([x5_idx;x4_idx;x2_idx])+ymax-(2.98+2*i),'FaceColor',[1,1,1],'EdgeAlpha',0);
%         area(x+2,ymax-(3+2*i)*ones(length(x),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
        
        % plot ca
        x_ca = w/31:w/31:w;
        ca = tempresp_avg(:,i);
        ca = (ca-min(ca))/(max(ca)-min(ca));
        ca = h*ca/(1.3*max(ca));
        
        hold on;
        plot(x_ca+2+w,ca+ymax-(2.98+2*i),'LineWidth',1,'Color',[0 0 0]);
        
        area(x_ca+2+w,ca+ymax-(2.98+2*i),'FaceColor',[0.75,0.75,0.75],'EdgeAlpha',0);
        area(x_ca+2+w,ymax-(3+2*i)*ones(length(x_ca),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
        
        area([1 2+2*w],ymax-(3+2*i)*ones(2,1),'FaceColor',[1 1 1],'EdgeAlpha',0);
    elseif i>=14 & i<=17
        area(x+2,strat_bin+ymax-(2.98+2*i),'FaceColor',[1,1,0],'EdgeAlpha',0);
%         area(x+2,ymax-(3+2*i)*ones(length(x),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
        
        % plot ca
        x_ca = w/31:w/31:w;
        ca = tempresp_avg(:,i);
        ca = (ca-min(ca))/(max(ca)-min(ca));
        ca = h*ca/(1.3*max(ca));
        
        hold on;
        plot(x_ca+2+w,ca+ymax-(2.98+2*i),'LineWidth',1,'Color',[0 0 0]);
        
        area(x_ca+2+w,ca+ymax-(2.98+2*i),'FaceColor',[0.75,0.75,0.75],'EdgeAlpha',0);
        area(x_ca+2+w,ymax-(3+2*i)*ones(length(x_ca),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
        
        area([1 2+2*w],ymax-(3+2*i)*ones(2,1),'FaceColor',[1 1 1],'EdgeAlpha',0);
    elseif i>=18 & i<=21
        area(x(x2_idx)+2,strat_bin(x2_idx)+ymax-(2.98+2*i),'FaceColor',[0.5,0,0],'EdgeAlpha',0);
%         area(x([x5_idx;x4_idx;x1_idx])+2,strat_bin([x5_idx;x4_idx;x1_idx])+ymax-(2.98+2*i),'FaceColor',[1,1,1],'EdgeAlpha',0);
%         area(x+2,ymax-(3+2*i)*ones(length(x),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
        
        % plot ca
        x_ca = w/31:w/31:w;
        ca = tempresp_avg(:,i);
        ca = (ca-min(ca))/(max(ca)-min(ca));
        ca = h*ca/(1.3*max(ca));
        
        hold on;
        plot(x_ca+2+w,ca+ymax-(2.98+2*i),'LineWidth',1,'Color',[0 0 0]);
        
        area(x_ca+2+w,ca+ymax-(2.98+2*i),'FaceColor',[0.75,0.75,0.75],'EdgeAlpha',0);
        area(x_ca+2+w,ymax-(3+2*i)*ones(length(x_ca),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
        
        area([1 2+2*w],ymax-(3+2*i)*ones(2,1),'FaceColor',[1 1 1],'EdgeAlpha',0);
    elseif i>=22 & i<=26
        area(x(x2_idx)+2,strat_bin(x2_idx)+ymax-(2.98+2*i),'FaceColor',[0.5,0,0],'FaceAlpha',1,'EdgeAlpha',0);
        area(x(x4_idx)+2,strat_bin(x4_idx)+ymax-(2.98+2*i),'FaceColor',[1,0,0],'FaceAlpha',0.5,'EdgeAlpha',0);
%         area(x([x5_idx;x1_idx])+2,strat_bin([x5_idx;x1_idx])+ymax-(2.98+2*i),'FaceColor',[1,1,1],'EdgeAlpha',0);
%         area(x+2,ymax-(3+2*i)*ones(length(x),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
        
        % plot ca
        x_ca = w/31:w/31:w;
        ca = tempresp_avg(:,i);
        ca = (ca-min(ca))/(max(ca)-min(ca));
        ca = h*ca/(1.3*max(ca));
        
        hold on;
        plot(x_ca+2+w,ca+ymax-(2.98+2*i),'LineWidth',1,'Color',[0 0 0]);
        
        area(x_ca+2+w,ca+ymax-(2.98+2*i),'FaceColor',[0.75,0.75,0.75],'EdgeAlpha',0);
        area(x_ca+2+w,ymax-(3+2*i)*ones(length(x_ca),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
        
        area([1 2+2*w],ymax-(3+2*i)*ones(2,1),'FaceColor',[1 1 1],'EdgeAlpha',0);
    elseif i>=27 & i<=29
        area(x(x4_idx)+2,strat_bin(x4_idx)+ymax-(2.98+2*i),'FaceColor',[1,0,0],'FaceAlpha',0.5,'EdgeAlpha',0);
%         area(x([x5_idx;x2_idx;x1_idx])+2,strat_bin([x5_idx;x2_idx;x1_idx])+ymax-(2.98+2*i),'FaceColor',[1,1,1],'EdgeAlpha',0);
%         area(x+2,ymax-(3+2*i)*ones(length(x),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
        
        % plot ca
        x_ca = w/31:w/31:w;
        ca = tempresp_avg(:,i);
        ca = (ca-min(ca))/(max(ca)-min(ca));
        ca = h*ca/(1.3*max(ca));
        
        hold on;
        plot(x_ca+2+w,ca+ymax-(2.98+2*i),'LineWidth',1,'Color',[0 0 0]);
        
        area(x_ca+2+w,ca+ymax-(2.98+2*i),'FaceColor',[0.75,0.75,0.75],'EdgeAlpha',0);
        area(x_ca+2+w,ymax-(3+2*i)*ones(length(x_ca),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
        
        area([1 2+2*w],ymax-(3+2*i)*ones(2,1),'FaceColor',[1 1 1],'EdgeAlpha',0);
    elseif i>=30 & i<=33
        area(x+2,strat_bin+ymax-(2.98+2*i),'FaceColor',[1,1,0],'EdgeAlpha',0);
%         area(x+2,ymax-(3+2*i)*ones(length(x),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
        
         % plot ca
         x_ca = w/31:w/31:w;
         ca = tempresp_avg(:,i);
         ca = (ca-min(ca))/(max(ca)-min(ca));
         ca = h*ca/(1.3*max(ca));
         
         hold on;
         plot(x_ca+2+w,ca+ymax-(2.98+2*i),'LineWidth',1,'Color',[0 0 0]);
         
         area(x_ca+2+w,ca+ymax-(2.98+2*i),'FaceColor',[0.75,0.75,0.75],'EdgeAlpha',0);
         area(x_ca+2+w,ymax-(3+2*i)*ones(length(x_ca),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
         
         area([1 2+2*w],ymax-(3+2*i)*ones(2,1),'FaceColor',[1 1 1],'EdgeAlpha',0);
    elseif i>=34 & i<=47
        area(x(x5_idx)+2,strat_bin(x5_idx)+ymax-(2.98+2*i),'FaceColor',[0,1,0],'FaceAlpha',0.5,'EdgeAlpha',0);
%         area(x([x4_idx;x2_idx;x1_idx])+2,strat_bin([x4_idx;x2_idx;x1_idx])+ymax-(2.98+2*i),'FaceColor',[1,1,1],'EdgeAlpha',0);
        
        % plot ca
        x_ca = w/31:w/31:w;
        ca = tempresp_avg(:,i);
        ca = (ca-min(ca))/(max(ca)-min(ca));
        ca = h*ca/(1.3*max(ca));
        
        hold on;
        plot(x_ca+2+w,ca+ymax-(2.98+2*i),'LineWidth',1,'Color',[0 0 0]);
        
        area(x_ca+2+w,ca+ymax-(2.98+2*i),'FaceColor',[0.75,0.75,0.75],'EdgeAlpha',0);
        area(x_ca+2+w,ymax-(3+2*i)*ones(length(x_ca),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
        
        a = area([1 2+2*w],ymax-(3+2*i)*ones(2,1),'FaceColor',[1 1 1],'EdgeAlpha',0);
    end
    
   
%     if i<=13
%         area(x_ca(2:end-1)+3.5,ca(2:end-1)+ymax-(2.98+2*i),'FaceColor',[0,0.4,0],'EdgeAlpha',0);
%         area(x_ca+3.5,ymax-(3+2*i)*ones(length(x_ca),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
%     elseif i>=14 & i<=17
%         area(x_ca(2:end-1)+3.5,ca(2:end-1)+ymax-(2.98+2*i),'FaceColor',[1,1,1],'EdgeAlpha',0);
%         area(x_ca+3.5,ymax-(3+2*i)*ones(length(x_ca),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
%     elseif i>=18 & i<=20
%         area(x_ca(2:end-1)+3.5,ca(2:end-1)+ymax-(2.98+2*i),'FaceColor',[0.6,0,0],'EdgeAlpha',0);
%         area(x_ca+3.5,ymax-(3+2*i)*ones(length(x_ca),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
%     elseif i>=21 & i<=26
%         area(x_ca(2:end-1)+3.5,ca(2:end-1)+ymax-(2.98+2*i),'FaceColor',[0.8,0,0],'EdgeAlpha',0);
%         area(x_ca+3.5,ymax-(3+2*i)*ones(length(x_ca),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
%     elseif i>=27 & i<=29
%         area(x_ca(2:end-1)+3.5,ca(2:end-1)+ymax-(2.98+2*i),'FaceColor',[1,0,0],'EdgeAlpha',0);
%         area(x_ca+3.5,ymax-(3+2*i)*ones(length(x_ca),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
%     elseif i>=30 & i<=33
%         area(x_ca(2:end-1)+3.5,ca(2:end-1)+ymax-(2.98+2*i),'FaceColor',[1,1,1],'EdgeAlpha',0);
%         area(x_ca+3.5,ymax-(3+2*i)*ones(length(x_ca),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
%     elseif i>=34 & i<=46
%         area(x_ca(2:end-1)+3.5,ca(2:end-1)+ymax-(2.98+2*i),'FaceColor',[0,1,0],'EdgeAlpha',0);
%         area(x_ca+3.5,ymax-(3+2*i)*ones(length(x_ca),1),'FaceColor',[1 1 1],'EdgeAlpha',0);
%     end
   
end
base = a.BaseLine;
base.Visible = 'off';
axis off;

line([2.38 2.38],[1 ymax-3],'LineStyle',':','LineWidth',0.5,'Color',[0 0 0]);
% line([2.55 2.55],[1 ymax-3],'LineStyle',':','LineWidth',1,'Color',[0 0 0]);
line([2.72 2.72],[1 ymax-3],'LineStyle',':','LineWidth',0.5,'Color',[0 0 0]);
% text(2.38,-2,'Off SAC','HorizontalAlignment','center','Rotation',-90);
% text(2.72,-2,'On SAC','HorizontalAlignment','center','Rotation',-90);

line([2+w+w*8/31 2+w+w*8/31],[1 ymax-3],'LineStyle',':','LineWidth',0.5,'Color',[0 0 0]);
line([2+w+w*16/31 2+w+w*16/31],[1 ymax-3],'LineStyle',':','LineWidth',0.5,'Color',[0 0 0]);

line([1.7485+w 2+w],[0.2 0.2],'LineWidth',2,'Color',[0 0 0]);
line([2+w+w*23.1875/31 2+2*w],[0.2 0.2],'LineWidth',2,'Color',[0 0 0]);
text(1.8742+w,-1,'10um','HorizontalAlignment','center','FontName','Arial','FontSize',font_size);
text((4+3*w+w*23.1875/31)/2,-1,'1s','interpreter','tex','HorizontalAlignment','center','FontName','Arial','FontSize',font_size);





