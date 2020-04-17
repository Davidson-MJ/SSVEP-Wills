% plot the changes in normalized alpha power over time, for a subset of
% elecctrodes.


%if alpha power decrease contralateral to attended target, then we should
%compare PFI on the LEFT, at RIGHT electrodes.

% as per Handal & Jansen, (2014). Normalise power diff:
%[L-R/L+R (left sensors)] - [ L-/L+R (right sensors)]. 

%%
% Contra - Ipsi
% we want to compare the contralateral response, for Left and Right target
% disappearances.
% for Right electrodes, focus on Left Disap.
PFI_contraL1= squeeze(mean(GFX_leftPFI(:, Rchans,:),2)) ;
PFI_contraL2= squeeze(mean(GFX_rightPFI(:, Rchans,:),2)) ;

% for the topoplots:

PFI_contraR1 = squeeze(mean(GFX_rightPFI(:, Lchans,:),2)); 
PFI_contraR2 = squeeze(mean(GFX_leftPFI(:, Lchans,:),2));% -squeeze(mean(GFX_rightReap(:, Rchans,:),2));
%% take the mean
% tmp=cat(3,PFI_contraL1, PFI_contraL2);
% CONTRAL = mean(tmp,3);
% tmp=cat(3,PFI_contraR1, PFI_contraR2);
% CONTRAR= mean(tmp,3);
%% take the difference
CONTRAL = PFI_contraL1- PFI_contraL2;
CONTRAR = PFI_contraR1- PFI_contraR2;

%% 

 
 
figure(46); clf
     normON=1;
     legH=[];
     for iPFI = 1:2
         switch iPFI
             case 1
                 plotd = CONTRAL;
                 topod = GFX_leftPFI
                 colis= 'k';
                 Subchans = Rchans;
             case 2
                plotd= CONTRAR;
                 colis= 'm';
                 Subchans = Lchans;
         end
         
         %chansubset
%          plotd = squeeze(mean(dtype(:,Subchans,:),2));
         
         
         if normON==1 % normalize within ppant.
         tmp = zeros(size(plotd));
             for isub = 1:size(plotd,1)
                 %divide by baseline per subject.
                 baselineis = dsearchn(timeax', [-2.5 -2]');
                 tmpbase = mean(plotd(isub, baselineis(1):baselineis(2)));
                 tmp(isub,:)= (plotd(isub,:)./ tmpbase) -1;
             end
             plotd=tmp;
         end
             
         plotMean = squeeze(mean(plotd,1));
         stE = CousineauSEM(plotd);
         
         
         %show topoplot of entire period, and marked electrodes.
%          subplot(3,3,iPFI);
%          tp = squeeze(mean(mean(dtype,3),1));
%          topoplot(tp, elocs(1:64), 'conv', 'on','emarker2', {[Subchans], '*' 'w'});
     
         %now hilbs.     
         subplot(3,3,4:9);
         sh=shadedErrorBar(timeax, plotMean, stE, colis,1);
         sh.mainLine.LineWidth=3;
         legH(iPFI) = sh.mainLine;
         
         
        hold on 
%      ttestdata(iPFI,:,:) = plotd;
     
%      legend([legH(1), legH(2)], 'Contra PFI left', 'Contra PFI right', 'autoupdate', 'off')
     hold on;
     %
     plot([0 0], ylim, ['k-']); axis tight
     % test for sig:
     pvals= [];
     for itime = 1:size(plotd,2)
         
        [~, pvals(itime)] = ttest(plotd(:,itime),0); %ttestdata(2,:,itime));
        
        if pvals(itime)<.05
            text(timeax(itime), -.09, '*','fontsize', 15)
        end
        
     end
     hold on
     plot(xlim, [0 0], ['k-'])
     end
         ylabel({['Alpha amplitude'];['(% change)']});
     set(gcf, 'color', 'w')
     set(gca, 'fontsize', 15)
     
%      legend('PFI', 'reap') 
xlim([-2 2])
     xlabel('Time from button press')
     %%