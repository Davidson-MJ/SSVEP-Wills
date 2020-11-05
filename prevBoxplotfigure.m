figure(6); clf % arrange columns as f1, IM , f2 (1,3,2)
set(gcf, 'units', 'normalized', 'position', [0 .55 .4 .4]);
tmp = LOO_firstSIG;
tmp(:,2) = LOO_firstSIG(:,3);
tmp(:,3) = LOO_firstSIG(:,2);


%add boxplot and reorient.
 boxplot(tmp, 'Colors', 'k', 'Symbol', ['']);
 boxcol = {'k', 'm', 'b'};   
 SCcol = {'b', 'm', 'k'};
for ibox = 1:size(tmp,2)
 h=findobj(gca, 'tag','Box');
 h(ibox).LineWidth = 1;
 h(ibox).Color = boxcol{ibox};
 %shift box.
%  h(ibox).Xdata = h(ibox).Xdata +.1;
 h=findobj(gca, 'tag','Median');
 h(ibox).LineWidth = 3;
 h(ibox).Color = boxcol{ibox};
 h=findobj(gca, 'tag','Upper Whisker');
 h(ibox).LineWidth = 1;
 h(ibox).Color = boxcol{ibox};
 h=findobj(gca, 'tag','Lower Whisker');
 h(ibox).LineWidth = 1; 
 h(ibox).Color = boxcol{ibox};
 h=findobj(gca, 'tag','Lower Adjacent Value');
 h(ibox).LineWidth = 1;
 h(ibox).Color = boxcol{ibox};
 h=findobj(gca, 'tag','Upper Adjacent Value');
 h(ibox).LineWidth = 1;
 h(ibox).Color = boxcol{ibox};
 
 hold on
 % overlay individual data points.
 
     useD =squeeze(tmp(:,ibox));
     
     ph = repmat(ibox, [length(useD),1]);
     if ibox==1
         ph = ph+.35; % offset to centre (shifting right).
     elseif ibox==3         
         %offset to the middle (shift left)
         ph = ph-.35;
     end
 sc=scatter(ph, useD);
 sc.LineWidth = 2;
 sc.SizeData = 200;
 sc.MarkerFaceColor = [.9 .9 .9];
 sc.MarkerEdgeColor = SCcol{ibox};
 sc.MarkerFaceAlpha = 0;
 %
 
 end
 set(gca, 'Xticklabels', {'f1' 'IM', 'f2'})
 hold on; plot(xlim, [0, 0], ['k',':'], 'linew', 2);
 %%
 % add connecting lines between f2 and f1.
for isubs = 1:size(tmp,1)    
    hold on
    if ~isnan(tmp(isubs,:));
    xp =[1+.35, 2, 3-.35];
    yp =[tmp(isubs,1), tmp(isubs,2),  tmp(isubs,3)];
        
        %xpoints
    else % skip the IM.
    xp =[1+.35,3-.35];
    yp =[tmp(isubs,1),  tmp(isubs,3)];
    end
     % sort by temporal order:
    [ypn, id]=sort(yp, 'ascend');
    line(xp(id),ypn, 'Color', [.8 .8 .8],'linestyle', '-');
%     line(xp,yp, 'Color', [.8 .8 .8],'linestyle', '-');
    
end
%
view([90 -90]);
    axis ij    
    set(gca, 'ydir', 'normal', 'fontsize', 25);
 shg
 ylabel('time from button press (s)')
 ylim([-1.5 0.5]);
 set(gcf, 'color', 'w')
 