
clear all
addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/EEG')
basefol=pwd;
%%
pdirs = dir([pwd filesep '*_*' 'EEG']);

% allppants=[1,2,4,6,9:16,18]; %
 allppants=[1,2,4,6,7,9:19]; %

    %% %% %
    job2.plotSpacedTimetopo=0; %1 x 4 topos
    job2.plotComparativeTimetopo = 0; % new version, settled on - 0.5s, compares two topos at Disap and Reap in 2x2
    job2.plotMeanTIMEtopo_andtvals=0;
    job2.plotSpatialCorrelation_overtime=1; % this is for 
    job2.plotgroupedSNRovertime=0;
    
    
%
%
%
%

%subtract mean per channel first?
normEEGcorr=0;

%
%
%


 %%%%% which data type to plot?
%     usePFIorCatch=2; % 1 for PFI, 2 for Catch.
%
getelocs;

if job2.plotSpacedTimetopo==1
    cd(basefol)
    cd('GFX_Pre-RESS')
    %%
    icount=1;
    clf
    
    
    pochans= [17:32,51:64];% chans to plot;
%     pochans=1:64;
    
    %     colormap('parula')
    peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; % don't change!
    for usePFIorCatch=2%1:2
        
        if usePFIorCatch==1
            useD='PFI';
        else
            useD='Catch';
        end
        
        
%         cmaxbyHz = [.5, .5, .1,.1,.1,.1, .025]; %PFI
        cmaxbyHz = [.5, .75, .1,.1,.1,.1, .025]; %PFI
        
        for ihz=[1,2,3,4,7]
            usehz = peakfreqsare(ihz);
            cd(basefol)
            cd('GFX_Pre-RESS')
            clf
            
            load(['GFX_' useD 'performance_withSNR_' num2str(usehz) '_allchan'])
            icount=1;
            for itimezero=1:2%onset and offset
                
                
                
                switch itimezero
                    case 1
                        dPLOT = storeacrossPpant_onsetSNR_chans;
                        dirtype = 'disappearance';
                        %                     titlep=
                    case 2
                        dPLOT = storeacrossPpant_offsetSNR_chans;
                        dirtype = 'reappearancee';
                end
                TIMING = [ -.5 0.5 ];
                
                
                
                %adjust topo values by mean subtraction per eeg channel?
                if normEEGcorr==1;
                    dPLOTtmp=zeros(size(dPLOT));
                    
                    for ippant=1:size(dPLOT,1)
                        for ichan = 1:size(dPLOT,2)
                            prevsnr=squeeze(dPLOT(ippant,ichan,:));
                            thism=squeeze(mean(dPLOT(ippant,ichan,:)));
                            tmpsub= repmat(thism, [1, size(prevsnr,2)]);
                            %subtract channel mean from snr data.
                            dPLOTtmp(ippant,ichan,:) = prevsnr-tmpsub;
                        end
                        
                    end
                    dPLOT=dPLOTtmp;
                    
                    
                    
                end
                
                
                
                
                
                figure(1)
                
                timeIND = dsearchn([tgrm-3]', TIMING');
                %%
                %     icount=1
                
                hold on
                pvals=nan(64,length(timeIND));
                tvals=nan(64,length(timeIND));
                ipl=1;
                for itime=1:length(timeIND)
                    
                    subplot(1,4,icount)
                    
                    tid= timeIND(itime);
                    ipl= ipl + (length(TIMING)*(ihz-1));
                    
                    
                    
                    %remove mean per channel?
                    
                    
                    
                    
                    % plot sig.
                    timeTOPO=squeeze(dPLOT(:,:,tid));
                    
                    for ichan=1:64
                        [~,pvals(ichan,itime),~, stat]= ttest(timeTOPO(:,ichan), 0, 'tail', 'both');
                        tvals(ichan,itime)=stat.tstat;
                        
                    end
                    q=fdr(pvals(:),.05);
                    pmask = (pvals(:,itime)<q);
                    
                    
                    %         subplot(2,length(timeIND),ipl)
                    %         subplot(1,2,ipl)
                    %        tp= topoplot(squeeze(mean(dPLOT(:,:,tid),1)), elocs(1:64), 'pmask', pmask, 'conv', 'on');
                    tp= topoplot(squeeze(mean(dPLOT(:,:,tid),1)), elocs([1:64]), 'conv', 'on');
%                     tp= topoplot(squeeze(mean(dPLOT(:,:,tid),1)), elocs([1:64]), 'plotchans', pochans);
                    
                    %else plot tscores
                    %        tp= topoplot(tvals, elocs(1:64), 'pmask', pmask, 'conv', 'on');
                    
                    set(findobj(gca,'type','patch'),'facecolor',get(gcf,'color'))
                    if normEEGcorr==1
                        caxis([-.15 .15]);
%                         caxis([-.05 .05]);
                    else
                        %                         caxis([0 max(squeeze(mean(dPLOT(:,:,tid),1)))])
                        caxis([0 cmaxbyHz(ihz)])
                    end
                    %plot colorbar?
                    if icount==5
%%
                        c= colorbar;
                        ylabel(c, [ num2str(usehz) ' Hz, log(SNR)']);
%                         ylabel(c, ['\Delta log(SNR)']);
%                         ylabel(c, ['\Delta log(SNR)']);
                        
                    end
                    
                    
%                     title({[ num2str(TIMING(itime)) 's']})
                    set(gca, 'fontsize', 25)
                    icount=icount+1;
                end
                
                %%
                
                c=colormap('viridis');        %
                %%
%                 c(1,:)=[ 0 0 0];
                colormap(c)
                %         subplot(3,1,3)
                %         colorbar; caxis([-1 1])
                set(gcf, 'color', 'w')
                
                %%
            end
                %%
                cd(basefol)
                cd ../
                cd(['Figures' filesep 'GFX spatial correlations'])
                %%
%                 suptitle(['Top row ' useD ' Disappearance, bottom row  ' useD ' reappearance'])
                shg
                print('-dpng', ['SNR Topo time for ' num2str(usehz) ' during ' useD ])
            
        end
    end
end


if job2.plotComparativeTimetopo==1
    cd(basefol)
    cd('GFX_Pre-RESS')
    %%
    icount=1;
    
    %     colormap('parula')
    peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; % don't change!
    ymaxperfreq = [1, 1, .5, 1, .5 , .5, .1, .5, .5];
    for usePFIorCatch=1:2
        
        if usePFIorCatch==1
            useD='PFI';
        else
            useD='Catch';
        end
        
        
    onsetChans_both= zeros(2, 16, 64, 14);
    offsetChans_both= zeros(2, 16, 64, 14);
    
    
    %cycle through them all.
    hzcompareAll = [1, 2; 1,3;1,4;1,7;...
                 2, 3; 2,4;2,7;...
                 3,4;3,7;4,7];
    
             
    for allhzcombos = 1:size(hzcompareAll,1)
    hzcounter=1;
    hzcompare= hzcompareAll(allhzcombos,:);
    
    
    cd(basefol)
    cd('GFX_Pre-RESS')
    figure(1)
    clf
    
    
    for ihz = hzcompare
        usehz= peakfreqsare(ihz);
        load(['GFX_' useD 'performance_withSNR_' num2str(usehz) '_allchan'])
        onsetChans_both(hzcounter,:,:,:) = storeacrossPpant_onsetSNR_chans;
        offsetChans_both(hzcounter,:,:,:) = storeacrossPpant_offsetSNR_chans;        
        hzcounter=hzcounter+1;
    end
    
    %
    onsetChans_20=squeeze(onsetChans_both(1,:,:,:));
    offsetChans_20=squeeze(offsetChans_both(1,:,:,:));
    onsetChans_40=squeeze(onsetChans_both(2,:,:,:));
    offsetChans_40=squeeze(offsetChans_both(2,:,:,:));
    
    figure(1); 
    leg=[];
%     ttestdata= zeros(2, size(onsetChans_20,1),size(onsetChans_20,3));
    
for isubplot = 1:4
    switch isubplot
        case 1
            dPLOT=onsetChans_20; 
            hzused = peakfreqsare(hzcompare(1));
        case 2
            dPLOT=onsetChans_40;
            hzused = peakfreqsare(hzcompare(2));
        case 3
            dPLOT=offsetChans_20; 
            hzused = peakfreqsare(hzcompare(1));
        case 4
            dPLOT=offsetChans_40;
            hzused = peakfreqsare(hzcompare(2));
    end
    % this orders the topos side by side ina  2 x2 .
            
    
    %set ymx
    hid = find(peakfreqsare==hzused);
    ymax = ymaxperfreq(hid);
    
        
                % single TIME.
                TIMING = [-.5 ];
                
                figure(1)
                
                timeIND = dsearchn([tgrm-3]', TIMING');   
                
                %adjust topo values by mean subtraction per eeg channel?
                if normEEGcorr==1;
                    dPLOTtmp=zeros(size(dPLOT));
                    
                    for ippant=1:size(dPLOT,1)
                        for ichan = 1:size(dPLOT,2)
                            prevsnr=squeeze(dPLOT(ippant,ichan,:));
                            thism=squeeze(mean(dPLOT(ippant,ichan,:)));
                            tmpsub= repmat(thism, [1, size(prevsnr,2)]);
                            dPLOTtmp(ippant,ichan,:) = prevsnr-tmpsub;
                        end
                        
                    end
                    dPLOT=dPLOTtmp;
                    
                    
                    
                end
                
              % plot now
              subplot(2,2,isubplot)
              tp= topoplot(squeeze(mean(dPLOT(:,:,timeIND),1)), elocs(1:64), 'conv', 'on');      
              %else plot tscores
              %        tp= topoplot(tvals, elocs(1:64), 'pmask', pmask, 'conv', 'on');
              
              set(findobj(gca,'type','patch'),'facecolor',get(gcf,'color'))
              
              %         topoplot(squeeze(mean(dPLOT(:,:,tid),1)), elocs(1:64), 'emarker2', {find(pmask), 'o', 'w', 2}); caxis([-1 1])
              %         if itime==length(timeIND)
              c= colorbar;
              ylabel(c, [ num2str(hzused) ' Hz, log(SNR)']);
              
              if normEEGcorr==1
                  caxis([-.15 .15]);
              else
                  
                 
                  caxis([0 ymax])

              end
              
              
              title({[ num2str(TIMING)   's']})
              set(gca, 'fontsize', 25)
              icount=icount+1;
     end
                
                %%
                c=colormap('viridis');        %
                %%
%                 c(1,:)=[ 0 0 0];
                colormap(c)
                %         subplot(3,1,3)
                %         colorbar; caxis([-1 1])
                set(gcf, 'color', 'w')
                
                %%
            
%            suptitle(['Top row ' useD ' Disappearance, bottom row  ' useD ' reappearance'])
%          
         cd(basefol)
         cd ../
         cd(['Figures' filesep 'GFX spatial correlations'])
    
                
                print('-dpng', ['SNR grouped Topo at -0.5s,' num2str(peakfreqsare(hzcompare(1))) 'and ' num2str(peakfreqsare(hzcompare(2))) ' during ' useD  '.png'])
    
    
    end
                %%
               
    end
    
                %%
            
end





if job2.plotMeanTIMEtopo_andtvals==1
    %%
    cd(basefol)
    cd('newplots-MD')
    %looking at PFIINC>PFIdec
    %         timewinds = [-3 -1.1]; %'20Hz early'
    %         timewinds = [-.75 .25];%'20Hz late'
    
    %         timewinds = [-3 -1.5]; %'40Hz early'
    %         timewinds = [-.5 0];%'40Hz late'
    
    %looking at baseline:
    % timewinds=[-1 0.24]; %PFI increase 20Hz from baseline
    % timewinds=[-1.28 0.1]; %PFI increase 40Hz from baseline
    %
    % timewinds=[-.52 1.8]; %PFI decrease 20Hz from baseline
    % timewinds=[-.5 .5]; %PFI decrease 40Hz from baseline
    ipl=1;
    clf
    
    %only use for offsets.
    rmvbase=0;
    
    
    
    for ihz=1%1:2%:2;
        cd(basefol)
        cd('newplots-MD')
        switch ihz
            case 1
                load('GFX_PFIperformance_withSNR_20_allchan')
            case 2
                load('GFX_PFIperformance_withSNR_40_allchan')
        end
        dPLOT=[];
        dPLOT(1,:,:,:) = storeacrossPpant_onsetSNR_chans;
        
        
        if rmvbase==1
            bsrem=[-3 -1];
            tbase=tgrm-3;
            tmp = zeros(size(storeacrossPpant_offsetSNR_chans));
            tIND = dsearchn(tbase', [bsrem]');
            for ippant=1:size(storeacrossPpant_offsetSNR_chans,1)
                for itrial = 1:size(storeacrossPpant_offsetSNR_chans,2)
                    tbase = squeeze(nanmean(storeacrossPpant_offsetSNR_chans(ippant,itrial,tIND(1):tIND(2)),3));
                    rbase = repmat(tbase, [1  size(storeacrossPpant_offsetSNR_chans,3)]);
                    
                    tmp(ippant,itrial,:) =squeeze(storeacrossPpant_offsetSNR_chans(ippant,itrial,:))' - rbase;
                end
            end
            storeacrossPpant_offsetSNR_chans=tmp;
        end
        
        dPLOT(2,:,:,:) = storeacrossPpant_offsetSNR_chans;
        
        
        
        %
        %     TIMING = [  -.5 -.4 -.2 0 .2 .4 .5];
        
        timeIND = dsearchn([tgrm-3]', timewinds');
        %
        hold on
        %take mean over timewindow
        dTest= mean(dPLOT(:,:,:,timeIND(1):timeIND(2)),4);
        
        pvals=[];
        tvals=[];
        
        
        for ichan=1:64
            [~,pvals(ichan),~, stat]= ttest(squeeze(dTest(2,:,ichan)));
            %             [~,pvals(ichan),~, stat]= ttest(squeeze(dTest(1,:,ichan)), squeeze(dTest(2,:,ichan)));
            tvals(ichan)=stat.tstat;
        end
        pmask = (pvals<.05);
        %
        topoplot(tvals, elocs(1:64), 'pmask', pmask)
        shg
        checkchans=find(pmask);
        % perform spatial cluster FDR:
        %% compare the increase in ASH to zero (since a diff already)
        % extract tvalues at spatially coincident p<05 electrodes,
        %(uncorrected)
        
        %         checkchans= [25,62,30,61,29,63,31]; % PFI 20Hz, early.
        %         checkchans= [29,30,61:63]; % PFI 20Hz, late.
        % checkchans= [4 9 18 23 24 25 26 29 30 47 54 56 57 58 60 62 63];
        % checkchans=[10,44,53,58,59,26,60:63,29:31];
        tvalsperchan = tvals(checkchans);
        
        % this the accrued test statistic (total), we have to
        % check against, after shuffling the data (Maris&Oostenveld)
        observedCV = sum(abs(tvalsperchan));
        %
        % now repeat the above process, but create random samples:
        sumTestStatsShuff = zeros(1,2000);
        for irand = 1:2000
            %testing the null that it isn't mismatched - matched at time 2
            % which creates a diff. so select from either!
            shD= zeros(2,length(tvalsperchan),21);
            for ipartition = 1:2
                for ippant = 1:21
                    for chan=1:length(checkchans)
                        if mod(randi(100),2)==0 %if random even number
                            pdata = dTest(1,randi(21), checkchans(chan)); %select both chans
                            
                        else %
                            pdata = dTest(2,randi(21), checkchans(chan)); %select both chans
                            
                        end
                        
                        shD(ipartition,chan,ippant) = pdata;
                    end
                end
            end
            
            %now compute difference between out hypothetical topoplots,
            % and test for sig, checking the accumulated test statistic at our
            % chans of interest
            tvalsperchan = zeros(1,length(checkchans));
            if length(checkchans>1)
                testdata = squeeze(shD(1,:,:)) - squeeze(shD(2,:,:));
            else
                testdata = (shD(1,:,:)) - (shD(2,:,:));
                testdata=testdata';
            end
            
            for itest = [1:length(checkchans)] %test each channel
                
                if length(checkchans)>1
                    [~, p, ~,stat]= ttest(testdata(itest,:));
                else
                    [~, p, ~,stat]= ttest(testdata);
                end
                
                tvalsperchan(1,itest) = stat.tstat;
            end
            
            sumTestStatsShuff(1,irand) = sum(abs(tvalsperchan));
        end %repeat nshuff times
        %
        subplot(211)
        %plot histogram:
        H=histogram(abs(sort(sumTestStatsShuff)));
        % fit CDF
        cdf= cumsum(H.Data)/ sum(H.Data);
        %the X values (actual CV) corresponding to .01
        [~,cv05uncorr] = (min(abs(cdf-.95)));
        [~,cv01uncorr] = (min(abs(cdf-.99)));
        [~,cv001uncorr] = (min(abs(cdf-.999)));
        %THE Q VALUE FOR OUR OBSERVED DATA:
        
        
        hold on
        pCV=plot([observedCV observedCV], ylim, ['r-']);
        p05=plot([H.Data(cv05uncorr) H.Data(cv05uncorr)], ylim, ['k:']);
        plot([H.Data(cv01uncorr) H.Data(cv01uncorr)], ylim, ['k:']);
        plot([H.Data(cv001uncorr) H.Data(cv001uncorr)], ylim, ['k:']);
        legend([pCV p05], {['observed'] ['p005'] })
        % what is the Pvalue?
        
        %
        if observedCV>H.Data(cv05uncorr)
            %observed pvalue in distribution=
            [~, c2] = min(abs(H.Data-observedCV)); %closest value.
            pvalis= 1-cdf(c2);
            title({['Spatial Cluster  Significant!'];['Tvals = ' num2str(observedCV) ', p =' num2str(pvalis)]})
            %display Q
            %              checkchans= [25,62,30,61,29,63,31];
            pmaskchecked = zeros(1,64);
            pmaskchecked(checkchans)=1;
        else
            title('Spatial Cluster  ns!')
        end
        
        
        
        
        
        
        
        
        %%
        subplot(2, 1, 2) ;
        % pmaskchecked(9)=[0]; PFI in 20 vs 0
        % pmaskchecked(39)=[0]; %PIF in 40 vs 0
        
        % pmaskchecked(45)=[0];pmaskchecked(19)=[0]; %PIF in 40 vs 0
        
        tp=topoplot(tvals, elocs(1:64), 'pmask', pmaskchecked, 'conv', 'on');
        %         topoplot(tvals, elocs(1:64), 'emarker2', {find(pmask), 'o', 'w', 20});
        
        c= colorbar;
        ylabel(c, 't-value')
        caxis([-4 4])
        set(c, 'location', 'Eastoutside')
        set(gca, 'fontsize', 25)
        %%
        title({['\it p \rm\bf < .001 for'];[sprintf('%.2f',tgrm(timeIND(1))-3) ' to ' sprintf('%.2f',(tgrm(timeIND(2))-3)) 's']})
        %%
        set(gcf, 'color', 'w')
        cd(basefol)
        cd('newplots-MD')
        cd('figures')
        
        shg
    end
    %%
    shg
    set(gcf, 'color', 'w')
    print('-dpng', 'Topo PreBP 20Hz late PFIin>PFIde.png')
end



    


if  job2.plotSpatialCorrelation_overtime==1
    cd(basefol)
    cd('GFX_Pre-RESS')
    figure(1)
    clf
    
                usenormalPFIordownsampled=2; % 1 for normal, 2 for downsample

    % can plot the mean correlation across subjs (corr1or2=1),
    %or take the across subj correlation(corr1or2=2).
    corr1or2=1;
    
    selectchans=[1:64]; %use whole head
    %             selectchans = [23:32,56:64]; % parieto-occipital only
    
     peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; % don't change!        

   
for usePFIorCatch=2%1:2
    
    if usePFIorCatch==1
            useD='PFI';
            falpha = .5;
    else
            useD='Catch';
            falpha = .15;
            usenormalPFIordownsampled=1;
    end
        
   
    
    onsetChans_both=[]; %zeros(2, 16, 64, 14);
    offsetChans_both=[];% zeros(2, 16, 64, 14);
    
    
    %cycle through them all.
    hzcompareAll = [1, 2; 1,3;1,4;1,7;...
                 2, 3; 2,4;2,7;...
                 3,4;3,7;4,7];
    
             % completed for hz 1,2,7. 

             
    for allhzcombos = [1,4,7] % these are the comparisons for tg,bg, and im.
    hzcounter=1;
    hzcompare= hzcompareAll(allhzcombos,:);
    
    
    cd(basefol)
    cd('GFX_Pre-RESS')
    figure(1)
    clf
    
    
    for ihz = hzcompare
        usehz= peakfreqsare(ihz);
        load(['GFX_' useD 'performance_withSNR_' num2str(usehz) '_allchan'])
        
        if usenormalPFIordownsampled==1
        onsetChans_both(hzcounter,:,:,:) = storeacrossPpant_onsetSNR_chans;
        offsetChans_both(hzcounter,:,:,:) = storeacrossPpant_offsetSNR_chans;        
        else
          
            onsetChans_both(hzcounter,:,:,:,:) = storeacrossPpant_onsetSNR_chans_downsampled;
            offsetChans_both(hzcounter,:,:,:,:) = storeacrossPpant_offsetSNR_chans_downsampled;
            
        end
        hzcounter=hzcounter+1;
    end
    
    %
    if usenormalPFIordownsampled==1
    onsetChans_20=squeeze(onsetChans_both(1,:,:,:));
    offsetChans_20=squeeze(offsetChans_both(1,:,:,:));
    onsetChans_40=squeeze(onsetChans_both(2,:,:,:));
    offsetChans_40=squeeze(offsetChans_both(2,:,:,:));
    else % extra dimension.
        onsetChans_20=squeeze(onsetChans_both(1,:,:,:,:));
    offsetChans_20=squeeze(offsetChans_both(1,:,:,:,:));
    onsetChans_40=squeeze(onsetChans_both(2,:,:,:,:));
    offsetChans_40=squeeze(offsetChans_both(2,:,:,:,:));
    end
    
    figure(1); clf
    leg=[];
    if usenormalPFIordownsampled==1
        
        ttestdata= zeros(2, size(onsetChans_20,1),size(onsetChans_20,3));
    else
        ttestdata= zeros(2, size(onsetChans_20,1),size(onsetChans_20,4));
    end
    
    
    
    
    for iPFIdir=1:2%1:2 %onset and offset
        hold on
        switch iPFIdir
            case 1
                d1=onsetChans_20;
                d2=onsetChans_40;
                chis= 'target invisible';
                %calculate spatial correlation over time:
                colis='r';
                linestyle='--';
            case 2
                d1=offsetChans_20;
                d2=offsetChans_40;
                chis= 'target visible';
                %                            colis='k';
                linestyle='-';
        end
        
        
        
        
        if usenormalPFIordownsampled==1
            %if restricting the channels for comparison
            d1=d1(:,selectchans,:);
            d2=d2(:,selectchans,:);
        else
            d1=d1(:,selectchans,:,:);
            d2=d2(:,selectchans,:,:);
        end
        
        % normalize by within channel SNR if necessary.
        if normEEGcorr==1
            if usenormalPFIordownsampled==2
                error('code unfinished')
            end
            
            for id=1:2
                switch id
                    case 1
                        dPLOT=d1;
                    case 2
                        dPLOT=d2;
                end
                dPLOTtmp=zeros(size(dPLOT));
                
                for ippant=1:size(dPLOT,1)
                    for ichan = 1:size(dPLOT,2)
                        prevsnr=squeeze(dPLOT(ippant,ichan,:));
                        thism=squeeze(mean(dPLOT(ippant,ichan,:)));
                        tmpsub= repmat(thism, [1, size(prevsnr,2)]);
                        dPLOTtmp(ippant,ichan,:) = prevsnr-tmpsub;
                    end
                    
                end
                dPLOT=dPLOTtmp;
                
                switch id
                    case 1
                        d1=dPLOT;
                    case 2
                        d2=dPLOT;
                end
            end
            placesig=.33; % for placing results of ttests on figure
            
        else
            placesig=.66; % for placing results of ttests on figure
        end
        
        
        
        
        
        
     
                   if  corr1or2==1
                       % corr 1  = correlation across channels (within
                       % subject).
                       
                       if usenormalPFIordownsampled==1 % otherwise adjust for extra dimension, and plot across shuffles
                       corr_time = zeros(size(d1,1), size(d1,3));
                       p_time = zeros(size(d1,1), size(d1,3));
                       % per ppant, calculate corr over time
                       for ippant = 1:size(d1,1)
                           for itime= 1:size(d1,3)
                               
                               %                              scatter(squeeze(d1(ippant,:,itime))', squeeze(d2(ippant,:,itime))');
                               [r,p]= corr(squeeze(d1(ippant,:,itime))', squeeze(d2(ippant,:,itime))');
                               corr_time(ippant,itime)=r;
                               p_time(ippant,itime)=r;
                           end
                           
                       end
                       
                       
                       %within subj error bars
                       
                       
                       %adjust standard error as per COusineau(2005)
                       %confidence interval for within subj designs.
                       % y = x - mXsub + mXGroup,
                       x = corr_time;
                       
                       mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
                       mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
                       
                       %for each observation, subjtract the subj average, add
                       %the group average.
                       NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
                       
                       
                       corr_time=NEWdata;
                       
                       
                       %                    subplot(1,2,iPFIdir)
                       mP= squeeze(mean(corr_time,1));
                       
                       
                       stP = std(corr_time)/sqrt(size(corr_time,1));
                       
                       
                       
                       else % we are plotting median and shuffled CI.
                           
                           
                            corr_time = zeros(size(d1,1), size(d1,3), size(d1,4));
                            p_time = corr_time;
                       
                            % per ppant, calculate corr over time
                       for ippant = 1:size(d1,1)
                           for ishuff=1:size(d1,3)
                               for itime= 1:size(d1,4)
                                   
                                   [r,p]= corr(squeeze(d1(ippant,:,ishuff,itime))', squeeze(d2(ippant,:, ishuff,itime))');
                                   corr_time(ippant,ishuff,itime)=r;
%                                    p_time(ippant,itime)=p;
                               end
                           end
                           
                       end
                       
                       
                       %within subj error bars
                       %% sanity check of distribution at each time point:
%                        figure(2); clf
%                        for it=1:size(corr_time,3)
%                            
%                            subplot(6,6,it);
%                            histogram(squeeze(corr_time(3,:,it)));
%                            title([num2str(tgrm(it)-3) 's'])
%                            xlabel('Pearson r ')
%                        end
                       %%
                       corr_time=squeeze(mean(corr_time,2));
                       %adjust standard error as per COusineau(2005)
                       %confidence interval for within subj designs.
                       % y = x - mXsub + mXGroup,
                       x = corr_time;
                       
                       mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
                       mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
                       
                       %for each observation, subjtract the subj average, add
                       %the group average.
                       NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
                       
                       
                       corr_time=NEWdata;
                       
                       
                       %                    subplot(1,2,iPFIdir)
                       mP= squeeze(mean(corr_time,1));
                       
                       
                       stP = std(corr_time)/sqrt(size(corr_time,1));
                       
                       
                           
                           
                           
                           
                       end
                           
                       %%
                       figure(1);
                       st=shadedErrorBar(tgrm-3, mP, stP,[],1);
                       st.mainLine.Color= colis;
                       st.mainLine.LineStyle=linestyle;
                       st.mainLine.LineWidth=3;
                       st.patch.FaceColor=colis;
                       st.patch.FaceAlpha = falpha;
                       st.edge(1).Color=colis;
                       st.edge(2).Color=colis;
                       
                       leg(iPFIdir)=st.mainLine;
                           
                       

                       
                       %%
                       ttestdata(iPFIdir,:,:) =x;
                       
                       
                       
                       axis tight
                       figure(1); % adjust sig value placement?
                       
                       yl=get(gca, 'ylim');
                       
                       %extend by 10% (keeps sig points in relative space).
                       dyl=diff(yl)*.1;
                       ylim([yl(1)-dyl yl(2)+dyl])
                       
                       %% plot sig points if ver1 subj
                       if iPFIdir==2
                           %plot sig
                           pvals = zeros(1,size(ttestdata,3));
                           for itime=1:size(ttestdata,3)
                               
                               [h, pvals(itime) ]= ttest(ttestdata(1,:,itime), ttestdata(2,:,itime)) ;
                               
                           end
                           q=fdr(pvals, .05);
                           sigspots = find(pvals<=q);
                           
                           for itime = 1:length(sigspots)
                               
                               rtime = sigspots(itime);
                               hold on
                               
                               %                             plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', sh.mainLine.Color)
                               plot(tgrm(rtime)-3, yl(2)+dyl, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'k')
                           end
                           
                       end
                   else
                   %sanity check - compare to across subj correlation
                   da=squeeze(nanmean(d1,1));
                   db=squeeze(nanmean(d2,1));
                   corr_time_gfx=zeros(1,33);
                   p_time=corr_time_gfx;
                   
                   for itime= 1:size(d1,3)
                       
%                              scatter(squeeze(d1(ippant,:,itime))', squeeze(d2(ippant,:,itime))');
                          [r,p]= corr(squeeze(da(:,itime)), squeeze(db(:,itime)));
                   corr_time_gfx(itime)=r;
                       p_time(itime)=p;
                   
                   end
                   
                   
                   
                   plot(tgrm-3, corr_time_gfx, colis);
                   
                   sigspots= find(p_time<.05);
                   
                   axis tight
                   figure(1); % adjust sig value placement?
                   
                       yl=get(gca, 'ylim');
                   
                   %extend by 10% (keeps sig points in relative space).
                   dyl=diff(yl)*.1;
                   ylim([yl(1)-dyl yl(2)+dyl])

                   for itime = 1:length(sigspots)
                               
                               rtime = sigspots(itime);
                               hold on
                               
                               %                             plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', sh.mainLine.Color)
                               plot(tgrm(rtime)-3, yl(2), ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'k')
                   end
                   
                   
                   end
        
        
        %                    ylim([-.1 .3])
        axis tight
        xlim([-2.5 2.5])
        
        
        ylabel({[num2str(peakfreqsare(hzcompare(1))) ' vs ' num2str(peakfreqsare(hzcompare(2))) ' Hz '];['spatial correlation [\it r\rm ]']})
%         ylabel({['1f vs. 2f correlation [\itr\rm ]']})
%         xlabel(['Time from ' chis ])
        xlabel(['Time from ' useD ' report'])
        set(gca, 'fontsize', 25)
        set(gcf, 'color', 'w')
        
        
        
        %check distribution of r values at each time point? (for
        %Nao)
        %%
        %                    clf
        %                    for i=1:33;
        %                 %plot raw
        %                             subplot(7,5,i); hist(x(:,i),22); title([num2str(tgrm(i)-3) 's'])
        %
        %                        % or compare with logit transform
        %                     shd=x(:,i);
        %                     lp = log(shd)-log(1-shd);                    %
        %
        %               [lpn,id]=sort(lp);%convert to zscore
        %               %sort first.
        %             tz= zscore((lpn));
        %
        %                        subplot(7,5,i); hist(tz,22); title([num2str(tgrm(i)-3) 's'])
        %                    end
        %                    suptitle('During target disappearance')
        % %
        % %
    end
    
    
    %                legend([leg(1) leg(2)], {'target invisible' 'target visible'})
    %
 yl=get(gca, 'ylim');                   
                   %extend by 10% (keeps sig points in relative space).
                   dyl=diff(yl)*.1;
                   ylim([yl(1)-dyl yl(2)+dyl])

    
                   xlim([-1.75 1.75])
    %
    cd(basefol)
    cd ../
    cd('Figures')
    %
    cd('GFX spatial correlations')
    %
%     title(['During ' useD ])
    print('-dpng', ['SpatialCorrelation freqs ' num2str(peakfreqsare(hzcompare(1))) ' ' num2str(peakfreqsare(hzcompare(2))) ' ' useD])
    
    end % Hz combos.
end % PFI or Catch?
end
%%

         
       if  job2.plotgroupedSNRovertime==1;
                 cd(basefol)
        cd('newplots-MD')
             figure(2)
            clf
               
            
            % can plot the mean correlation across subjs (corr1or2=1), 
            %or take the across subj correlation(corr1or2=2).
            corr1or2=1; 
            
            selectchans=[1:64]; %use whole head
%             selectchans = [23:32,56:64]; % parieto-occipital only
            
                         load('GFX_PFIperformance_withSNR_20_allchan')
               onsetChans_20 = storeacrossPpant_onsetSNR_chans;
               offsetChans_20 = storeacrossPpant_offsetSNR_chans;
                         load('GFX_PFIperformance_withSNR_40_allchan')
                 onsetChans_40 = storeacrossPpant_onsetSNR_chans;
               offsetChans_40 = storeacrossPpant_offsetSNR_chans;
                 
               %%
               figure(1); clf
               leg=[];
               
               % also store for final plot.
               groupedharmonicCorrelations = zeros(2, length(allppants), 6, 33);
               for iPFIdir=1:2%1:2 %onset and offset
                  hold on
                   switch iPFIdir
                       case 1
                           d1all=onsetChans_20;
                           d2all=onsetChans_40;
                           chis= 'target invisible';
                   %calculate spatial correlation over time:
                   colis='r';
                   linestyle='-';
                       case 2
                           d1all=offsetChans_20;
                           d2all=offsetChans_40;
                           chis= 'target visible';
%                            colis='k';
linestyle=':';
                   end
                   
                   
                   legprint={};
                   for changroup=1:6
                       switch changroup
                           case 1
                               % Frontal 10 channels:
                               selectchans = [1:7, 33:40]; %n=15 chans
                               chang='Frontal';   
                               colis='r';
                           case 2
                               selectchans=[41:43, 47,48, 8,9,12,13];
                               chang='LeftTemp';
                               colis='k';
                           case 3
                               selectchans=[44:46, 49,50, 10,11,15,16];
                               chang='RightTemp';
                               colis= [.5 .5 .5];
                           case 4
                               selectchans=[17:19,23,24,25,51:53, 56,57];
                               chang='LeftParietal';
                               colis = 'b';
                           case 5
                               selectchans=[20:22,25:27,53:55,58,59];
                               chang='RightParietal';
                               colis = [0 0 .5];
                           case 6
                               selectchans=[60:64, 28:32];
                               chang='Occipital';
                               colis='m';
                               
                       end
                           
                   
                   %if restricting the channels for comparison
                   d1=d1all(:,selectchans,:);
                   d2=d2all(:,selectchans,:);
                   
                   
                   % normalize by within channel SNR if necessary.
                    if normEEGcorr==1;
                        for id=1:2
                            switch id
                                case 1
                                    dPLOT=d1;
                                case 2
                                    dPLOT=d2;
                            end
                                dPLOTtmp=zeros(size(dPLOT));
                
                                for ippant=1:size(dPLOT,1)
                                    for ichan = 1:size(dPLOT,2)
                                        prevsnr=squeeze(dPLOT(ippant,ichan,:));
                                        thism=squeeze(mean(dPLOT(ippant,ichan,:)));
                                        tmpsub= repmat(thism, [1, size(prevsnr,2)]);
                                        dPLOTtmp(ippant,ichan,:) = prevsnr-tmpsub;
                                    end
                                    
                                end
                                dPLOT=dPLOTtmp;
                                
                                switch id
                                    case 1
                                        d1=dPLOT;
                                    case 2
                                        d2=dPLOT;
                                end
                        end
                    placesig=.33; % for placing results of ttests on figure

                    else
                        placesig=.66; % for placing results of ttests on figure
                    end
            
                   
                   
                   
                   
                   
                   for ipl=1:3 % SNR then correlation
                       placeis = iPFIdir + 2*(ipl -1);
                       
                       subplot(3,2,placeis);
                       
                       hold on
                       switch ipl
                           case 1
                               datapl =squeeze(mean(d1,2)); %mean over chanss
                               tis = '1f raw SNR';
                           case 2
                               datapl =squeeze(mean(d2,2)); %mean over chans
                               tis = '2f raw SNR';
                           case 3
                               tis = '1f vs 2f [r]';
                       
                               corr_time = zeros(size(d1,1), size(d1,3));
                               p_time = zeros(size(d1,1), size(d1,3));
                               % per ppant, calculate corr over time
                               for ippant = 1:size(d1,1)
                                   for itime= 1:size(d1,3)
                                       
                                       %                              scatter(squeeze(d1(ippant,:,itime))', squeeze(d2(ippant,:,itime))');
                                       [r,p]= corr(squeeze(d1(ippant,:,itime))', squeeze(d2(ippant,:,itime))');
                                       corr_time(ippant,itime)=r;
                                       p_time(ippant,itime)=r;
                                   end
                                   
                               end
                               
                               datapl= corr_time;
                               
                       end
                       %within subj error bars
                       
                       
                       %adjust standard error as per COusineau(2005)
                       %confidence interval for within subj designs.
                       % y = x - mXsub + mXGroup,
                       x = datapl;
                       
                       mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
                       mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
                       
                       %for each observation, subjtract the subj average, add
                       %the group average.
                       NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
                       
                       
                       corr_time=NEWdata;
                       
                       
                       %                    subplot(1,2,iPFIdir)
                       mP= squeeze(mean(corr_time,1));
                       
                       
                       stP = std(corr_time)/sqrt(size(corr_time,1));
                       
                       
                       st=shadedErrorBar(tgrm-3, mP, stP,[],1);
                       st.mainLine.Color= colis;
                       st.mainLine.LineStyle=linestyle;
                       st.mainLine.LineWidth=3;
                       st.patch.FaceColor=colis;
                       st.edge(1).Color=colis;
                       st.edge(2).Color=colis;
                       
                       leg(changroup)=st.mainLine;
                      title(tis)     
                      set(gca, 'fontsize', 15)
                  axis tight
                  xlabel(['Time from reporting'])
                  ylabel(tis)
                   end
                   
%                    ylim([-.1 .3])
axis tight
                  
                   legprint=[legprint  {[chang ',n=' num2str(length(selectchans))]}];
                   
                   
                   ylabel(tis)
                   if changroup==6;
                   legend([leg], legprint);
                   end
                   
                   % store the correlation if that is what we are plotting,
                   if ipl==3
                   groupedharmonicCorrelations(iPFIdir, :, changroup,:)= datapl;
                   
                   end
                   end
                   
                   
                   
%                    xlabel(['Time from ' chis ])
               
%                    
               end
               %%
               figure(2)
               clf
               for idir=1:2
                   switch idir
                       case 1
                           linet=':'
                       case 2
                           linet='-'
                   end
                           hold on
              for igroup=6%:6
                  plot(squeeze(nanmean(groupedharmonicCorrelations(idir,:, igroup,:),2)), linet)
              end
               end
               shg
%                legend([leg(1) leg(2)], {'target invisible' 'target visible'})
               %%
               axis tight
%                ylim([.475 .675])
             %%
             cd('Figures')
             %%
             print('-dpng', 'grouped SNR BG freqs')
         end
    
         
    
     
% end