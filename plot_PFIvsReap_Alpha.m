 
     normON=1;
     ntype=1; % 1 for division, 2 for baselinesub.
     
     legH=[];
     figure(1); clf
     for iPFI = 1%:2
         
         switch iPFI
             case 1
                 dtype = allPFI;
%                  dtype= GFX_rightPFI;
                 colis= 'b';
%                  Subchans = Lchans;
             case 2
                 dtype = all_leftPFI;
                 colis='k';
%                  Subchans = Rchans;
         end
         
         %first plot the alpha power, baseline window (across channels),
         %before normalization.
         
         %show topoplot of certain periods, and marked electrodes.
         if iPFI==1
             clf
         subplot(3,3,1);
         %which time windows for the topoplots?
         t1 = dsearchn([timeax'],[-2.5 -1.5]');

         tp = squeeze(nanmean(nanmean(dtype(:,:,t1(1):t1(2)),3),1));
         topoplot(tp, elocs(1:64), 'conv', 'on','emarker2', {[Subchans], 'o' 'w', 2}); shg;
         c=colorbar; 
         caxis([5 35 ])
         ylabel(c, 'Alpha amp. (\rm\muV)', 'fontsize', 15)
         set(gca, 'fontsize', 15)
         end
         %
         
         if normON==1 % normalize within ppant.
         tmp = zeros(size(dtype));
         %divide by baseline per subject.
         baselineis = dsearchn(timeax', [-2.5 -2]');
         
         for isub = 1:size(dtype,1)
                for ichan = 1:64
                 tmpchan = squeeze(dtype(isub, ichan, :));
                 tmpbase = nanmean(tmpchan(baselineis(1):baselineis(2)));
                 tmp(isub,ichan,:)= (tmpchan./ tmpbase) -1;
                end
             end
             dtype=tmp;
         end
         %
%              if iPFI==1
         subplot(3,3,iPFI+1);
         %which time windows for the topoplots?
         t2 = dsearchn([timeax'],[-1 2]');
         
         tp = squeeze(nanmean(nanmean(dtype(:,:,t2(1):t2(2)),3),1));
         topoplot(tp, elocs(1:64), 'conv', 'on','emarker2', {[Subchans], 'o' 'w', 2});
         c=colorbar;
         ylabel(c, '% change', 'fontsize', 15)
         caxis([-.1 .1])
         set(gca, 'fontsize', 15)
%              end
         %now hilbs.   
          
         %chansubset
         plotd = squeeze(nanmean(dtype(:,Subchans,:),2));       
         plotMean = squeeze(nanmean(plotd,1));
         stE = CousineauSEM(plotd);
         
        stE = std(plotd)./ sqrt(size(plotd,1));
         subplot(3,3,4:9); hold on
         sh=shadedErrorBar(timeax, plotMean, stE, colis,1);
         sh.mainLine.LineWidth=3;
         legH(iPFI) = sh.mainLine;
         set(gca, 'fontsize', 25)         
        hold on 
        
        %also add patches
        pch = patch([-2.5 -2.5 -1.5 -1.5], [-.1 .1 .1 -.1], [.9 .9 .9]);
        pch.FaceAlpha = .3;
        pch.LineStyle = 'none';
        
        pch = patch([-1 -1 2 2], [-.1 .1 .1 -.1], [1 .9 .9]);
        pch.FaceAlpha = .3;
        pch.LineStyle = 'none';
     end
     %
     
% create cold map 

     colormap('inferno');
     hold on;     
     plot([0 0], ylim, ['k-']); axis tight
     plot(xlim, [0 0], ['k-']); axis tight
     % test for sig:
     ylim([-.11 .05])
     axis  tight
     
     pvals= [];
     for itime = 1:size(plotd,2)
         
        [~, pvals(itime)] = ttest(plotd(:,itime),0);%ttestdata(2,:,itime));
        
        if pvals(itime)<.05
            text(timeax(itime), -.1, '*','fontsize', 55)
        end
     end
         ylabel({['Alpha amplitude'];['(% change)']})
     set(gcf, 'color', 'w')
%      set(gca, 'fontsize', 15)
     
%      legend(legH(1), legH(2), 'Left PFI', 'right PFI') 
xlim([-2.5 2.5])
     xlabel('Time from button press')
     %% when all done, change topoplot crange
     
     tmp= cbrewer('div', 'RdBu', 100);
     colormap(flipud(tmp));
     %%
