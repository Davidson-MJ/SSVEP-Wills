%
addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')


basefol=pwd;
clearvars -except basefol allppants
dbstop if error

cd('EEG');
pdirs = dir([pwd filesep '*EEG']);


job.concatTOPOacrossppanst=1;

job.plotTOPOsacrossppants=1;
normON=0; %normalize across ppants?

job.plotTOPOsacrossppants_crunched=0;




getelocs;
allppants=[1,2,4,6,9:16,18]; %
% window=[-2 4];
window=[-3 3];

srate=250;

epochdur = sum(abs(window))*srate;

timeid = [0:1/srate:epochdur];
timeid= timeid-3;

onsetc = ceil(epochdur)/2;
% peakfreqsare=[20,40]; %hz
%timing
tt = 0:1/srate:60;



% which frequencies to analyze?/apply?
peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; % don't change!      so that correct RESS filter gets applied.

if job.concatTOPOacrossppanst==1
    
    %%
    
    acrossRESSSNRTOPO=zeros(length(allppants), length(peakfreqsare),64); %now 10 freqs f1 TGs x4, BG, 2BG, 2fTGsx4


    icounter=1;
    for ippant=allppants
        cd(basefol)
        cd('EEG')
        cd(pdirs(ippant).name)
        
        load('RESSfilterdata')
        acrossRESSSNRTOPO(icounter,:,:)= ressEVEC_byHz_MAPS;
        
        
        icounter=icounter+1;
    end
    %%
    cd(basefol)
    cd('EEG')
    cd('GFX_EEG-RESS')
    %%
    save('GFX_RESSmaps', 'acrossRESSSNRTOPO','peakfreqsare')
    
    
    
end
%%

if job.plotTOPOsacrossppants==1
    %%
    cd(basefol)
    cd('EEG')
    cd('GFX_EEG-RESS')
    load('GFX_RESSmaps')
    %%
    normON=0;
    %norm within ppants?
    if normON==1;
    %%
    
    % normalize each topo for plotting?
    for ippant = 1:size(acrossRESSSNRTOPO,1)
        for iHz = 1:size(acrossRESSSNRTOPO,2)
            
                acrossRESSSNRTOPO(ippant, iHz,:) = acrossRESSSNRTOPO(ippant, iHz,:)./ max(acrossRESSSNRTOPO(ippant, iHz, :));
            
            
        end
        
    end
    end
  
    
    meanTOPOs=squeeze(mean(acrossRESSSNRTOPO,1));
    %%
%     peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; 
    plotORDER = [1,4,2,5,3,6,7,8,9];
    titles={'(TG) f1', '(TG) 2f1', '(TG) 3f1',...
        '(BG) f2', '(BG) 2f2', '(BG) 3f2',...
        '(IM) f2-f1', '(IM) 2f2-f1', '(IM) f1 +f2'};
    %%
        
        for ipl=1:9
            subplot(3,3,ipl)
            % find the correct index.
            findx = find(plotORDER==ipl);
            plotme= squeeze(meanTOPOs(findx,:));
            
            % plot
            topoplot(plotme, elocs(1:64), 'conv', 'on');            
            c=colorbar;
            
%             title([num2str(peakfreqsare(findx)) ' Hz'])
                title(titles{ipl})
            ylabel(c, 'RESS filter weights')
            caxis([0 .5])
            set(gca, 'fontsize', 25)
        end
          set(gcf, 'color', 'w')
        shg
    %%
      
        cd(basefol)        
        cd('Figures')
    %%
    print('-dpng', ['GFX_RESS TOPOs by Background frequencies.png'])
    %
end

if job.plotTOPOsacrossppants_crunched==1
    %%
    cd(basefol)
    cd('newplots-MD')
    load('GFX_RESSmaps.mat');
    %%
    figure(1)
    clf
    % normalize each topo for plotting.
    for ippant = 1:size(acrossRESSSNRTOPO,1)
        for iHz = 1:6%size(acrossRESSSNRTOPO,2)
            for iloc=1:size(acrossRESSSNRTOPO,3)
                acrossRESSSNRTOPO(ippant, iHz, iloc,:) = acrossRESSSNRTOPO(ippant, iHz, iloc,:)./ max(acrossRESSSNRTOPO(ippant, iHz, iloc,:));
            end
            
        end
%         subplot(7,3,ippant)
%         topoplot(squeeze(mean(acrossRESSSNRTOPO(ippant,4,:,:),3)), elocs(1:64));
%         title(num2str(ippant))
    end
    %%
%     figure(2);
%     topoplot(squeeze(mean(acrossRESSSNRTOPO(7,:,iloc,:),2)), elocs(1:64));
%     c=colorbar;
%     ylabel(c, 'RESS filter weights')
%     caxis([-.4 .4])
%             set(gca, 'fontsize', 25)
%     set(gcf, 'color', 'w'); 
    %%
    meanTOPOs=squeeze(mean(acrossRESSSNRTOPO,1));
    %%
    %mean over LOC
    figure(2);
    clf
    for it=1:2
        
        switch it
            case 1
                % Tg freqs, all locs,
                dis = squeeze(mean(squeeze(meanTOPOs(1:4,:,:)),1));
                titleis={['Target flicker'];['1st harmonic']};
            case 2
                dis = squeeze(mean(squeeze(meanTOPOs(7:10,:,:)),1));
                titleis={['Target flicker'];['2nd harmonic']};
            case 3
                dis = squeeze(mean(squeeze(meanTOPOs(5,:,:)),1));                
                titleis={['Background flicker'];['1st harmonic']};
             case 4
                dis = squeeze(mean(squeeze(meanTOPOs(6,:,:)),1));                
                titleis={['Background flicker'];['2nd harmonic']};
            
        end
        %then mean over LOCS.
        dplot = squeeze(mean(dis,1));
        %plot the mean
        subplot(2,2,it);
        topoplot(dplot, elocs(1:64), 'conv', 'on')
        c=colorbar;
            title(titleis)
            ylabel(c, 'RESS filter weights')
            if it<3
                caxis([-.1 .1])
            else
            caxis([-.2 .2])
            end
            set(gca, 'fontsize', 20)
        end
        %%
        set(gcf, 'color', 'w')
        cd(basefol)
        cd('newplots-MD')
        cd('Figures')
%         %%
%         if it==1
%             print('-dpng', ['RESS TOPOs by target Locs.png'])
%         else
            print('-dpng', ['RESS TOPOs by Harmonics.png'])
%         end
    
    %%
    figure(2)
    clf
    titlesare={'20 Hz' , '40 Hz'};
    for ipl=1:2
        subplot(2,1,ipl)
        
        topoplot(squeeze(mean(meanTOPOs(ipl+4,:,:),2)), elocs(1:64), 'emarker2', {30, '*', 'w', 20,5});
        hold on
        topoplot(squeeze(mean(meanTOPOs(ipl+4,:,:),2)), elocs(1:64), 'emarker2', {30, 'o', 'k', 25, 5});
         c=colorbar;
        title(titlesare{ipl})
        ylabel(c, 'RESS filter weights')
        caxis([-.2 .2])
        set(gca, 'fontsize', 25)
    end
    set(gcf, 'color', 'w')
    %
    print('-dpng', ['RESS TOPOs by Background frequencies.png'])
    %
end
