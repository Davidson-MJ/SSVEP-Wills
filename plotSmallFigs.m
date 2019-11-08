% plot new condensed figures.
%%%%% NEW FOR WILLS DATA!!!
clear all
clearvars -except basefol allppants

job.gradedchangesinERPimage_PFI=0;
job.BPandSSVEPtimecourseacrossppants_group_combinePFIandCATCH=0; %same format as previous script (epoch PFI etc.)

job.compareAMOUNT_PFIandPMD=0;
job.compareAMOUNT_PFIandPMD_raincloudver=1;


% Follows the format of s3_CE... only now applying to RESS timeseries

%having epoched tgs and catch, POST RESS. apply FFT and SNR to windows when
%all targets are present/ away from Catch stimulus onset.


% close all

addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/mDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')


basefol=pwd;

dbstop if error
%%



% also tr x tr of catch responses.
% for catch


%%%%% %%%%%%% %%%%%% %%%%%% %%%%%%
%%%%% %%%%%%% %%%%%% %%%%%% %%%%%%
% PARAMETERS FOR ALL SPECTROGRAM! KEEP CONSISTENT
%%%%% %%%%%%% %%%%%% %%%%%% %%%%%%
%%%%% %%%%%%% %%%%%% %%%%%% %%%%%%
param_spcgrm.tapers = [1 1];
param_spcgrm.Fs= [250];
param_spcgrm.Fpass= [0 50];
param_spcgrm.trialave=0;
param_spcgrm.pad = 2;
% movingwin=[2,.15]; % increase to 2.5?
movingwin=[2.5,.25]; % 

%%%%% %%%%%%% %%%%%% %%%%%% %%%%%%
%%%%% %%%%%%% %%%%%% %%%%%% %%%%%%
% what size will the output be? (taken from inside mtspecgramc)
        %frequency:
Nwin=round(250*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*250); % number of samples to step through
nfft=max(2^(nextpow2(Nwin)+param_spcgrm.pad),Nwin);
f=getfgrid(250,nfft,[param_spcgrm.Fpass]);
% use to automatically preallocate variable size.
Nf=length(f);
        %timesteps:
Nwin=round(250*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*250); % number of samples to step through
N=1501; %size EEG epochs
winstart=1:Nstep:N-Nwin+1;
winmid=winstart+round(Nwin/2);
% use to automatically preallocate variable size.
Nt=winmid/250;




% begin.
getelocs;

cd(basefol)
cd('EEG')
  
pdirs = dir([pwd filesep '*EEG']);

%% load data to determine physical catch timing.

%remaining participants after behavioral data analysis/exclusion
% allppants=[1,2,4,6,9:16,18]; %
allppants=[1,2,4,6,7,9:19]; % NEW ppants.
window=[-3 3];

srate=250;

epochdur = sum(abs(window));

timeidDYN = [0:1/srate:epochdur];
timeidDYN= timeidDYN-3;
%%
onsetc = ceil(epochdur)/2;
% peakfreqsare=[20,40]; %hz
%timing
tt = 0:1/srate:60;

% which frequencies to analyze?/apply?
peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; % don't change!        



if job.compareAMOUNT_PFIandPMD==1;
  
    %plot comparison betwen PFI and PMD on nTrags. 
    % as above, with some tweaks.
    
    %%
    plotstacked=1; % stack PFI and PMD ontop of each other
    clf
    %can plot BARS next to trial-by-trial image, or snr-timecourse.
for     plotBarsorTimeseries=1%

    
    getelocs
    rmvbase=0;
    
    counter2=1; % for plotting threeway PFI
    counter=1;
    hzcounter=1;
    
    %each Hz has a separate ylim range to optimize plots:
    ylimsareHz=[1.1,1.8;... %15 hz
             2.7,3.1;...%20 Hz
             .5, 1;... % 30 Hz
             2.8,3.3;... % 40 Hz
             0.4,.9;... % 45
             2.2,2.5;...
             .3, .6]; % 5 Hz
    
         hzlabel={'f1', 'f2', '2f1', '2f2', '3f1', '3f2*', 'f2-f1'};
         
         
    for hzis=2%:6%[1, 2, 7]
        cd(basefol)
        cd('EEG')
        cd('GFX_EEG-RESS')
        usehz= peakfreqsare(hzis);
        switch hzis
            case {1,3,5}
                                
                col='b';
                flickeris = 'Target';
            case {2,4}
%                 ylimsare=[2.5, 3.2];
                
                col=[.2 .2 .2];
                flickeris = 'Background';
            case 6
                col=[.2 1 .2];
                flickeris = 'Overlap';                
            case 7
%                 ylimsare=[2.5, 3.2];
                
                col='m';
                flickeris = 'Intermodulation';
        end
        
        ylimsare=ylimsareHz(hzis,:);
        
        
        for iPFInCatch=1:2
        switch iPFInCatch
            case 1
        
        load(['GFX_PFIperformance_withSNR_' num2str(usehz) '_min0_RESS'])
            case 2
                load(['GFX_Catch_performance_withSNR_' num2str(usehz) '_min0_RESS'])
                
%                 %use catch report onset
%           DIMSare{1} = 'catch onset';
%           DIMSare{2} = 'catch offset';
%           DIMSare{3} = 'reporting catch onset';
%           DIMSare{4} = 'reporting catch offset';
%           DIMSare{5} = 'invisible catch onset';
                
                storeacrossPpant_onsetSNR = squeeze(storeacrossPpant_catchEVENTS_SNR(:,3,:,:));
                storeacrossPpant_offsetSNR = squeeze(storeacrossPpant_catchEVENTS_SNR(:,4,:,:));
        end
%         tgrm = timeidDYN;
        %%
        for itimezero=1%:2
            
            
            switch itimezero
                case 1
                    useBP=storeacrossPpant_onsetBP;
                    useSNR=storeacrossPpant_onsetSNR;
%                     useRTs=storeacrossPpant_onsetRTs;
%                     useTOPO=storeacrossPpant_onsetTOPO;
                    
                    ctype= 'target invisible';
                    chtype='button press';
                    
                case 2
                    
                    useBP=storeacrossPpant_offsetBP;
                    useSNR=storeacrossPpant_offsetSNR;
%                     useRTs=storeacrossPpant_offsetRTs;
%                     useTOPO=storeacrossPpant_offsetTOPO;
                    
                    ctype= 'target visible';
                    chtype='button release';
                    
                    
            end
            
            %grouped nPFI
            
              %take mean across ppants for erp images..
            acrBP= squeeze(nanmean(useBP,1));
            
            acrSNR=squeeze(mean(useSNR,1));
            %%
            
            groupednPFISNR=nan(length(allppants), 3);  % for barchart
            groupednPFISNR_ts=nan(length(allppants), 3, size(acrSNR,2));  % for SNR time-seriestrace
            if itimezero==1 % interested in after BP                
            BPidx=181:361;
            SNRidx= 8:14;
            else
                BPidx=1:180;
                SNRidx= 1:7;
            end
                groupedBEHflpd = fliplr(squeeze(nanmean(acrBP(:,BPidx),2))'); %takes mean over all SNR timepoints.
            %%
            groupednPFISNR_thirds=nan(length(allppants),3);
            
            for ippant = 1:size(useBP,1); 
%                % per ppant, collect trial indices:
                ppant_BEH= squeeze(mean(useBP(ippant,:,BPidx),3)); %takes mean over all BP time points.
                ppant_SNR=squeeze(mean(useSNR(ippant,:,SNRidx),3)); %takes mean over only 3 second window.
%                 ppant_SNR=squeeze(mean(useSNR(ippant,:,:),3)); %takes mean over all SNR timepoints.
                ppant_SNR_ts=squeeze(useSNR(ippant,:,:)); %takes mean over all SNR timepoints.
               maxnPFI=ceil(max(ppant_BEH));
%                
               %flip BEH /SNR to correct orientation
               ppantBEHflpd=fliplr(ppant_BEH);
               ppantSNRflpd=fliplr(ppant_SNR);
               ppant_SNR_ts=flipud(ppant_SNR_ts);
               % only collect bin width.
%              
%                 nPFIidx = dsearchn(ppantBEHflpd', [0:3]');
               nPFIidx = dsearchn(groupedBEHflpd', [0:4]');
               
               nPFIidx=unique(nPFIidx);
               
               if nPFIidx(1) ~=1
                   nPFIidx(1)=1;
               end
               
               %take mean for these trials (actually rows in trialxtrial)
               for iPFI=1:length(nPFIidx)-1;
                   
                   try indices = nPFIidx(iPFI):nPFIidx(iPFI+1)-1;
                   catch %if max
                       indices = nPFIidx(iPFI):100;
                   end
                   
                   if iPFI==4;
                       shg
                   end
                   groupednPFISNR(ippant, iPFI)= squeeze(nanmean(ppantSNRflpd(indices)));
                   
                   
                   groupednPFISNR_ts(ippant, iPFI,:)= squeeze(nanmean(ppant_SNR_ts(indices,:)));
% %                    %
%                    figure(3)
%                    clf
%                    plot(ppantSNRflpd);
%                    hold on
%                    plot([indices(1) indices(1)], ylim, ['k:'])
%                    plot([indices(end) indices(end)], ylim, ['k:'])
%                    title([num2str(squeeze(nanmean(ppantSNRflpd(indices))))])
               end
               
% works well.  

groupednPFISNR_quarters(ippant,1)= squeeze(nanmean(ppant_SNR(1:24)));
groupednPFISNR_quarters(ippant,2)= squeeze(nanmean(ppant_SNR(25:50)));
groupednPFISNR_quarters(ippant,3)= squeeze(nanmean(ppant_SNR(51:75)));
groupednPFISNR_quarters(ippant,4)= squeeze(nanmean(ppant_SNR(76:100)));

%                
            end
            %% using raw thirds?
%             groupednPFISNR_thirds=fliplr(groupednPFISNR_thirds);
%                         groupednPFISNR=fliplr(groupednPFISNR_thirds);
%             groupednPFISNR
            
            %%
            figure(1); 
            if plotBarsorTimeseries==1
   %place Bar charts beneath SNR
            
            
         mbar=(nanmean(groupednPFISNR,1));
         
         
         %adjust errorbars:
          %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            x = groupednPFISNR;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant)
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
        stB=nanstd(NEWdata,1)./sqrt(length(allppants));
        
%             stB=nanstd(groupednPFISNR,1)./sqrt(length(allppants));
        
%plot single:
b1=bar(mbar);
b1.FaceColor=col;
hold on
errorbar(mbar, stB, 'color', 'k','LineStyle', 'none')

set(gca, 'fontsize', 15, 'xticklabels', {'0\leq1', '1\leq2', '2\leq3', '\geq3'})
xlabel('amount of PFI');
ylabel('RESS log(SNR)')

colormap('viridis')
ylim(ylimsare)
shg
%             title([num2str(usehz) ' Hz (' num2str(hzis) 'f)'])
% axis tight
xlim([.5 4.5])


            counter=counter+1;
%%%%%%% ANOVA and LME stats if need be.             
            nppants=length(allppants);
            nconds = size(groupednPFISNR,2);
%             %RMANOVA?
            dataan= reshape(groupednPFISNR, [nppants*nconds,1]);
            conds = [ones(nppants,1); ones(nppants,1)*2;ones(nppants,1)*3;ones(nppants,1)*4];
                subs = repmat([1:length(allppants)]', [nconds,1]);
%             
%             [anret]=rmanova(dataan, conds, subs);
% %%
%             %LME to be sure:            
          tblA=table(dataan, conds, subs);
%%             %%
             UNterm1='dataan ~  conds + (1|subs)';
%              UNterm2='dataan ~  conds';
             UNterm2='dataan ~  1+ (1|subs)';
%              %creating models:
             lmeUN =fitlme(tblA, UNterm1);
             lmeUN2 =fitlme(tblA, UNterm2);
%              
    LMEout=compare(lmeUN2, lmeUN)
            %%
% %%
            else % plot time series instead.
%place Bar charts beneath SNR
            
                placeis = hzcounter+ 9;
            
            subplot(4,3,placeis);
                tbase = timeidDYN;
    lgsav=[];
% can also plot time-series:
for iPFI=1:4
    x= squeeze(groupednPFISNR_ts(:,iPFI,:));
    %adjust errorbars.
    
       
            mXppant =squeeze( nanmean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = nanmean(nanmean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
            %             %compute new stErr %which version?
            stE = nanstd(NEWdata)/sqrt(size(x,1));
            %
            hold on
            sh=shadedErrorBar(tbase, nanmean(x,1),stE,['-'],[1]);
                
                sh.mainLine.Color = col;
                sh.mainLine.LineWidth = iPFI*1.5;
                sh.patch.FaceColor = col;
                sh.edge(1).Color = col;
                sh.edge(2).Color = col;
                
   lgsav(iPFI)= sh.mainLine; 
end
%%
% legend([lgsav], [{'nPFI 0 \leq 1'}, {'nPFI 1 \leq 2'}, {'nPFI 2 \leq 3 or 4'}])
%%
set(gca,'fontsize', 15); set(gcf, 'color', 'w')
xlabel(['time from ' chtype])
ylabel('RESS log(SNR)')
axis tight
ylim([ylimsare(1) ylimsare(2)])

            end
%%
counter2=counter2+1;
    
            %%
           
       
        end 
        
        %stash bars for stacked plots:
        stackedB(iPFInCatch,:) = mbar;
        stackedB_sd(iPFInCatch,:) = stB;
        
      
        
        end
        if plotstacked==1;
            
            clf
            
            b1=[];
            b1=bar(stackedB');
            b1(1).FaceColor=col;
            b1(2).FaceColor= 'r';
            b1(1).BarWidth = 2;
            hold on
            spots = 1:size(stackedB,2);
            places = [(spots)-.15; (spots)+.15];
            
            errorbar(places',stackedB', stackedB_sd', 'color', 'k','LineStyle', 'none')
            %%
%             set(gca, 'fontsize', 15, 'xticklabels', {'0\leq1', '1\leq2', '2\leq3', '\geq3'})            
            xlabel('amount of PFI or PMD');
            

            set(gca, 'fontsize', 15)            
            ylabel('RESS log(SNR)')
            
            colormap('viridis')
%             ylim(ylimsare)
            shg
            %             title([num2str(usehz) ' Hz (' num2str(hzis) 'f)'])
            % axis tight
            xlim([.5 4.5])
            if hzis==1
                ylim([0 1.95])
            elseif hzis==2
                ylim([1.5 3.3])
            end
                
            alpha(.6);
            lg=legend('PFI', 'PMD'); set(lg, 'location', 'NorthWest')
        end
            
        hzcounter=hzcounter+1;
    end
end
    %%
    shg
    %%
    
    colormap('viridis')
    cd(basefol)
    cd('Figures')
    cd('GFX PFI trial by trial')
    shg
     print('-dpng', ['PFI graded SNR summary graded, ' ctype '.png'])
    
    
end




if job.compareAMOUNT_PFIandPMD_raincloudver==1
  
    %plot comparison betwen PFI and PMD on nTrags. 
    % as above, with some tweaks.
    
    %%
    plotstacked=1; % stack PFI and PMD ontop of each other
    clf
    d=[];
    %can plot BARS next to trial-by-trial image, or snr-timecourse.
for     plotBarsorTimeseries=1%

    
    getelocs
    rmvbase=0;
    
    counter2=1; % for plotting threeway PFI
    counter=1;
    hzcounter=1;
    
    %each Hz has a separate ylim range to optimize plots:
    ylimsareHz=[1.1,1.8;... %15 hz
             2.7,3.1;...%20 Hz
             .5, 1;... % 30 Hz
             2.8,3.3;... % 40 Hz
             0.4,.9;... % 45
             2.2,2.5;...
             .3, .6]; % 5 Hz
    
         hzlabel={'f1', 'f2', '2f1', '2f2', '3f1', '3f2*', 'f2-f1'};
         
         
    for hzis=1%:6%[1, 2, 7]
        cd(basefol)
        cd('EEG')
        cd('GFX_EEG-RESS')
        usehz= peakfreqsare(hzis);
        switch hzis
            case {1,3,5}
                                
                col='b';
                flickeris = 'Target';
            case {2,4}
%                 ylimsare=[2.5, 3.2];
                
                col=[.2 .2 .2];
                flickeris = 'Background';
            case 6
                col=[.2 1 .2];
                flickeris = 'Overlap';                
            case 7
%                 ylimsare=[2.5, 3.2];
                
                col='m';
                flickeris = 'Intermodulation';
        end
        
        ylimsare=ylimsareHz(hzis,:);
        
        
        for iPFInCatch=1:2
        switch iPFInCatch
            case 1
        
        load(['GFX_PFIperformance_withSNR_' num2str(usehz) '_min0_RESS'])
            case 2
                load(['GFX_Catch_performance_withSNR_' num2str(usehz) '_min0_RESS'])
                
%                 %use catch report onset
%           DIMSare{1} = 'catch onset';
%           DIMSare{2} = 'catch offset';
%           DIMSare{3} = 'reporting catch onset';
%           DIMSare{4} = 'reporting catch offset';
%           DIMSare{5} = 'invisible catch onset';
                
                storeacrossPpant_onsetSNR = squeeze(storeacrossPpant_catchEVENTS_SNR(:,3,:,:));
                storeacrossPpant_offsetSNR = squeeze(storeacrossPpant_catchEVENTS_SNR(:,4,:,:));
        end
%         tgrm = timeidDYN;
        %%
        for itimezero=1%:2
            
            
            switch itimezero
                case 1
                    useBP=storeacrossPpant_onsetBP;
                    useSNR=storeacrossPpant_onsetSNR;
%                     useRTs=storeacrossPpant_onsetRTs;
%                     useTOPO=storeacrossPpant_onsetTOPO;
                    
                    ctype= 'target invisible';
                    chtype='button press';
                    
                case 2
                    
                    useBP=storeacrossPpant_offsetBP;
                    useSNR=storeacrossPpant_offsetSNR;
%                     useRTs=storeacrossPpant_offsetRTs;
%                     useTOPO=storeacrossPpant_offsetTOPO;
                    
                    ctype= 'target visible';
                    chtype='button release';
                    
                    
            end
            
            %grouped nPFI
            
              %take mean across ppants for erp images..
            acrBP= squeeze(nanmean(useBP,1));
            
            acrSNR=squeeze(mean(useSNR,1));
            %%
            
            groupednPFISNR=nan(length(allppants), 3);  % for barchart
            groupednPFISNR_ts=nan(length(allppants), 3, size(acrSNR,2));  % for SNR time-seriestrace
            if itimezero==1 % interested in after BP                
            BPidx=181:361;
            SNRidx= 8:14;
            else
                BPidx=1:180;
                SNRidx= 1:7;
            end
                groupedBEHflpd = fliplr(squeeze(nanmean(acrBP(:,BPidx),2))'); %takes mean over all SNR timepoints.
            %%
            groupednPFISNR_thirds=nan(length(allppants),3);
            
            for ippant = 1:size(useBP,1); 
%                % per ppant, collect trial indices:
                ppant_BEH= squeeze(mean(useBP(ippant,:,BPidx),3)); %takes mean over all BP time points.
                ppant_SNR=squeeze(mean(useSNR(ippant,:,SNRidx),3)); %takes mean over only 3 second window.
%                 ppant_SNR=squeeze(mean(useSNR(ippant,:,:),3)); %takes mean over all SNR timepoints.
                ppant_SNR_ts=squeeze(useSNR(ippant,:,:)); %takes mean over all SNR timepoints.
               maxnPFI=ceil(max(ppant_BEH));
%                
               %flip BEH /SNR to correct orientation
               ppantBEHflpd=fliplr(ppant_BEH);
               ppantSNRflpd=fliplr(ppant_SNR);
               ppant_SNR_ts=flipud(ppant_SNR_ts);
               % only collect bin width.
%              
%                 nPFIidx = dsearchn(ppantBEHflpd', [0:3]');
               nPFIidx = dsearchn(groupedBEHflpd', [0:4]');
               
               nPFIidx=unique(nPFIidx);
               
               if nPFIidx(1) ~=1
                   nPFIidx(1)=1;
               end
               
               %take mean for these trials (actually rows in trialxtrial)
               for iPFI=1:length(nPFIidx)-1;
                   
                   try indices = nPFIidx(iPFI):nPFIidx(iPFI+1)-1;
                   catch %if max
                       indices = nPFIidx(iPFI):100;
                   end
                   
                   if iPFI==4;
                       shg
                   end
                   groupednPFISNR(ippant, iPFI)= squeeze(nanmean(ppantSNRflpd(indices)));
                   
                   
                   groupednPFISNR_ts(ippant, iPFI,:)= squeeze(nanmean(ppant_SNR_ts(indices,:)));
% %                    %
%                    figure(3)
%                    clf
%                    plot(ppantSNRflpd);
%                    hold on
%                    plot([indices(1) indices(1)], ylim, ['k:'])
%                    plot([indices(end) indices(end)], ylim, ['k:'])
%                    title([num2str(squeeze(nanmean(ppantSNRflpd(indices))))])
               end
               
% works well.  

groupednPFISNR_quarters(ippant,1)= squeeze(nanmean(ppant_SNR(1:24)));
groupednPFISNR_quarters(ippant,2)= squeeze(nanmean(ppant_SNR(25:50)));
groupednPFISNR_quarters(ippant,3)= squeeze(nanmean(ppant_SNR(51:75)));
groupednPFISNR_quarters(ippant,4)= squeeze(nanmean(ppant_SNR(76:100)));

%                
            end
            %% using raw thirds?
%             groupednPFISNR_thirds=fliplr(groupednPFISNR_thirds);
%                         groupednPFISNR=fliplr(groupednPFISNR_thirds);
%             groupednPFISNR
            
            %%
            figure(1); 
            if plotBarsorTimeseries==1
   %place Bar charts beneath SNR
            
            
         mbar=(nanmean(groupednPFISNR,1));
         
         
         %adjust errorbars:
          %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            x = groupednPFISNR;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant)
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
        stB=nanstd(NEWdata,1)./sqrt(length(allppants));
        
%             stB=nanstd(groupednPFISNR,1)./sqrt(length(allppants));
        
%plot single:
b1=bar(mbar);
b1.FaceColor=col;
hold on
errorbar(mbar, stB, 'color', 'k','LineStyle', 'none')

set(gca, 'fontsize', 15, 'xticklabels', {'0\leq1', '1\leq2', '2\leq3', '\geq3'})
xlabel('amount of PFI');
ylabel('RESS log(SNR)')

colormap('viridis')
ylim(ylimsare)
shg
%             title([num2str(usehz) ' Hz (' num2str(hzis) 'f)'])
% axis tight
xlim([.5 4.5])


            counter=counter+1;
%%%%%%% ANOVA and LME stats if need be.             
            nppants=length(allppants);
            nconds = size(groupednPFISNR,2);
%             %RMANOVA?
            dataan= reshape(groupednPFISNR, [nppants*nconds,1]);
            conds = [ones(nppants,1); ones(nppants,1)*2;ones(nppants,1)*3;ones(nppants,1)*4];
                subs = repmat([1:length(allppants)]', [nconds,1]);
%             
%             [anret]=rmanova(dataan, conds, subs);
% %%
%             %LME to be sure:            
          tblA=table(dataan, conds, subs);
%%             %%
             UNterm1='dataan ~  conds + (1|subs)';
%              UNterm2='dataan ~  conds';
             UNterm2='dataan ~  1+ (1|subs)';
%              %creating models:
             lmeUN =fitlme(tblA, UNterm1);
             lmeUN2 =fitlme(tblA, UNterm2);
%              
    LMEout=compare(lmeUN2, lmeUN)
            %%
% %%
            else % plot time series instead.
%place Bar charts beneath SNR
            
                placeis = hzcounter+ 9;
            
            subplot(4,3,placeis);
                tbase = timeidDYN;
    lgsav=[];
% can also plot time-series:
for iPFI=1:4
    x= squeeze(groupednPFISNR_ts(:,iPFI,:));
    %adjust errorbars.
    
       
            mXppant =squeeze( nanmean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = nanmean(nanmean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
            %             %compute new stErr %which version?
            stE = nanstd(NEWdata)/sqrt(size(x,1));
            %
            hold on
            sh=shadedErrorBar(tbase, nanmean(x,1),stE,['-'],[1]);
                
                sh.mainLine.Color = col;
                sh.mainLine.LineWidth = iPFI*1.5;
                sh.patch.FaceColor = col;
                sh.edge(1).Color = col;
                sh.edge(2).Color = col;
                
   lgsav(iPFI)= sh.mainLine; 
end
%%
% legend([lgsav], [{'nPFI 0 \leq 1'}, {'nPFI 1 \leq 2'}, {'nPFI 2 \leq 3 or 4'}])
%%
set(gca,'fontsize', 15); set(gcf, 'color', 'w')
xlabel(['time from ' chtype])
ylabel('RESS log(SNR)')
axis tight
ylim([ylimsare(1) ylimsare(2)])

            end
%%
counter2=counter2+1;
    
            %%
           
       
        end 
        
        %stash bars for stacked plots:
        stackedB(iPFInCatch,:) = mbar;
        stackedB_sd(iPFInCatch,:) = stB;
        
      
         
      for i = 1:size(groupednPFISNR,2)
          d{i,iPFInCatch}= NEWdata(:,i);
      end
        end
        %%
           
        
           %%
           clf
           if hzis==1
           rb= [0,0,1; 1,0,0];
           else
               rb=[0,0,0; 1,0,0];
           end
           
           %plot rain clouds
        str=rm_raincloud(d, rb, 0,'ks',.18);
            
 set(gcf, 'Position',[-1320,219,537,564])
            set(gca, 'fontsize', 25)            
            xlabel('RESS log(SNR)')
            ylabel('amount of PFI or PMD');
            lg=legend([str.p{1,1}, str.p{2,2}] ,['f' num2str(hzis) ' during PFI'], ['f' num2str(hzis) ' during PMD']); 
            axis tight
            
            set(lg, 'location', 'SouthEast')
        %%
            
        hzcounter=hzcounter+1;
    end
end
    %%
    shg
    %%
    
    colormap('viridis')
    cd(basefol)
    cd('Figures')
    cd('GFX PFI trial by trial')
    shg
     print('-dpng', ['PFI graded SNR summary graded, ' ctype '.png'])
    
    
end

%% %%%%
if job.gradedchangesinERPimage_PFI==1
  %% now plot across ppants.
    % as above, with some tweaks.
    
    %%
    
    clf
    %can plot BARS next to trial-by-trial image, or snr-timecourse.
for     plotBarsorTimeseries=1:2 % plots bars beneath, then timeseries.

    
    getelocs
    rmvbase=0;
    
    counter2=1; % for plotting threeway PFI
    
    
    
    counter=1;
    hzcounter=1;
    
    %each Hz has a separate ylim range to optimize plots:
    ylimsareHz=[1.1,1.8;... %15 hz
             2.7,3.1;...%20 Hz
             .5, 1;... % 30 Hz
             2.8,3.3;... % 40 Hz
             0.4,.9;... % 45
             2.2,2.5;...
             .3, .6]; % 5 Hz
    
         hzlabel={'f1', 'f2', '2f1', '2f2', '3f1', '3f2*', 'f2-f1'};
         
         
    for hzis=2%:6%[1, 2, 7]
        cd(basefol)
        cd('EEG')
        cd('GFX_EEG-RESS')
        usehz= peakfreqsare(hzis);
        switch hzis
            case {1,3,5}
                                
                col='b';
                flickeris = 'Target';
            case {2,4}
%                 ylimsare=[2.5, 3.2];
                
                col=[.2 .2 .2];
                flickeris = 'Background';
            case 6
                col=[.2 1 .2];
                flickeris = 'Overlap';                
            case 7
%                 ylimsare=[2.5, 3.2];
                
                col='m';
                flickeris = 'Intermodulation';
        end
        
        ylimsare=ylimsareHz(hzis,:);
        
        load(['GFX_PFIperformance_withSNR_' num2str(usehz) '_min0_RESS'])
%         tgrm = timeidDYN;
        %%
        for itimezero=1%:2
            
            
            switch itimezero
                case 1
                    useBP=storeacrossPpant_onsetBP;
                    useSNR=storeacrossPpant_onsetSNR;
%                     useRTs=storeacrossPpant_onsetRTs;
%                     useTOPO=storeacrossPpant_onsetTOPO;
                    
                    ctype= 'target invisible';
                    chtype='button press';
                    
                case 2
                    
                    useBP=storeacrossPpant_offsetBP;
                    useSNR=storeacrossPpant_offsetSNR;
%                     useRTs=storeacrossPpant_offsetRTs;
%                     useTOPO=storeacrossPpant_offsetTOPO;
                    
                    ctype= 'target visible';
                    chtype='button release';
                    
                    
            end
            
            
            
            %                 for ip= 1:length(allppants)
            %                     trialind = [1:24] + 24*(ip-1);
            %                     acrBP(trialind,:,:) = useBP(ip,:,:);
            %                     acrSNR(trialind,:,:)=    useSNR(ip,:,:);
            %                     acrRTs(trialind)=    useRTs(ip,:);
            %
            %                 end
            
            %take mean across ppants for erp images..
            acrBP= squeeze(nanmean(useBP,1));
%             acrBP= squeeze(nanmean(useBP(21,:,:),1));
            
            acrSNR=squeeze(mean(useSNR,1));
%             acrSNR=squeeze(nanmean(useSNR(21,:,:),1));
            
            
            
            
            %sort across all
            %sort trial index by longest RT first.
            %%
            figure(1)
            hold on
            %%
            subplot(4,4,2:3)
            
            %         [sortedRTs, cid] = sort(acrRTs, 'descend');
            
            %
            imagesc([-3:1/60:3], 1:size(acrBP,1), acrBP);%(cid,:));
            %
            title({['PFI ' num2str(ctype)]}, 'fontsize', 15)
            c=colorbar;
            ylabel(c, 'buttons pressed')
            %                 ylabel('resampled catch trials')
            xlim([-1.75 1.75])
            set(gca, 'fontsize', 15, 'ytick', [])
            hold on
            plot([0 0] , ylim, ['w:'], 'linewidth', 3)
            plot([0 0] , ylim, ['k:'], 'linewidth', 2)
            %                 subplot(4,2,7)
            
            ppantMeanBP= squeeze(mean(useBP,2));
            %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            x = ppantMeanBP;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
            %compute new stErr %which version?
            %             stE = std(NEWdata)/sqrt(size(x,1));
            %                 shadedErrorBar([-3:1/60:3], mean(ppantMeanBP,1),stE)
            %                 set(gca, 'fontsize', 15)
            %                 hold on;
                            ylabel({['Normalized'];['trial count']})
            %                 plot([0 0] , ylim, ['k:'], 'linewidth', 1)
%             axis square
            %%
%             xlabel(['Time from ' ctype])
            
            %SNR
                placeis = hzcounter + 3;

            subplot(4,3, placeis);
            
            %reorder SNR
            acrSNRsort=acrSNR;
            
            %%
            %                 %smooth across trials
            sm_snrgrm20=zeros(size(acrSNRsort));
            for itime=1:size(acrSNRsort,2)
                sm_snrgrm20(:,itime)= smooth(acrSNRsort(:,itime),15);
            end
            
            
            %     imagesc(acrSNR);
            imagesc(timeidDYN, 1:size(acrSNRsort,1), acrSNRsort);
            c=colorbar;

%             ylabel(c, {['RESS log(SNR)']})

if hzis~=7    
            title([flickeris ' ' num2str(usehz) ' Hz (' hzlabel{hzis} ')'])
else
            title([flickeris ' ' num2str(usehz) ' Hz (f2-f1)'])
end
%             title(['BG ' num2str(hzis) 'F'])
            
            hold on
            plot([0 0] , ylim, ['k:'], 'linewidth', 3)
            plot([0 0] , ylim, ['w:'], 'linewidth', 2)
            
            
            xlim([-1.75 1.75])
            set(gca, 'fontsize', 15, 'ytick', [])
            %                 subplot(4,2,8)
            %
            %plot across ppant trace
            ppantMeanSNR= squeeze(mean(useSNR,2));
            
            %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            x = ppantMeanSNR;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
        
            xlabel(['Time from ' ctype])
            ylabel({'Normalized'; 'trial count'})
                caxis([ylimsare])
%                 axis square
                
            set(gca, 'fontsize', 15)
            cd(basefol)
            set(gcf, 'color', 'w')
            shg
            if hzis==7;
                ylabel(c, 'RESS log(SNR)')
            end
            
            % beside this, plot the bar!
            %grouped nPFI
            %%
            groupednPFISNR=nan(length(allppants), 3);  % for barchart
            groupednPFISNR_ts=nan(length(allppants), 3, size(acrSNR,2));  % for SNR time-seriestrace
            if itimezero==1 % interested in after BP                
            BPidx=181:361;
            SNRidx= 8:14;
            else
                BPidx=1:180;
                SNRidx= 1:7;
            end
                groupedBEHflpd = fliplr(squeeze(nanmean(acrBP(:,BPidx),2))'); %takes mean over all SNR timepoints.
            %%
            groupednPFISNR_thirds=nan(length(allppants),3);
            
            for ippant = 1:size(useBP,1); 
%                % per ppant, collect trial indices:
                ppant_BEH= squeeze(mean(useBP(ippant,:,BPidx),3)); %takes mean over all BP time points.
                ppant_SNR=squeeze(mean(useSNR(ippant,:,SNRidx),3)); %takes mean over only 3 second window.
%                 ppant_SNR=squeeze(mean(useSNR(ippant,:,:),3)); %takes mean over all SNR timepoints.
                ppant_SNR_ts=squeeze(useSNR(ippant,:,:)); %takes mean over all SNR timepoints.
               maxnPFI=ceil(max(ppant_BEH));
%                
               %flip BEH /SNR to correct orientation
               ppantBEHflpd=fliplr(ppant_BEH);
               ppantSNRflpd=fliplr(ppant_SNR);
               ppant_SNR_ts=flipud(ppant_SNR_ts);
               % only collect bin width.
%              
%                 nPFIidx = dsearchn(ppantBEHflpd', [0:3]');
               nPFIidx = dsearchn(groupedBEHflpd', [0:4]');
               
               nPFIidx=unique(nPFIidx);
               
               if nPFIidx(1) ~=1
                   nPFIidx(1)=1;
               end
               
               %take mean for these trials (actually rows in trialxtrial)
               for iPFI=1:length(nPFIidx)-1;
                   
                   try indices = nPFIidx(iPFI):nPFIidx(iPFI+1)-1;
                   catch %if max
                       indices = nPFIidx(iPFI):100;
                   end
                   
                   if iPFI==4;
                       shg
                   end
                   groupednPFISNR(ippant, iPFI)= squeeze(nanmean(ppantSNRflpd(indices)));
                   
                   
                   groupednPFISNR_ts(ippant, iPFI,:)= squeeze(nanmean(ppant_SNR_ts(indices,:)));
% %                    %
%                    figure(3)
%                    clf
%                    plot(ppantSNRflpd);
%                    hold on
%                    plot([indices(1) indices(1)], ylim, ['k:'])
%                    plot([indices(end) indices(end)], ylim, ['k:'])
%                    title([num2str(squeeze(nanmean(ppantSNRflpd(indices))))])
               end
               
% works well.  

groupednPFISNR_quarters(ippant,1)= squeeze(nanmean(ppant_SNR(1:24)));
groupednPFISNR_quarters(ippant,2)= squeeze(nanmean(ppant_SNR(25:50)));
groupednPFISNR_quarters(ippant,3)= squeeze(nanmean(ppant_SNR(51:75)));
groupednPFISNR_quarters(ippant,4)= squeeze(nanmean(ppant_SNR(76:100)));

%                
            end
            %% using raw thirds?
%             groupednPFISNR_thirds=fliplr(groupednPFISNR_thirds);
%                         groupednPFISNR=fliplr(groupednPFISNR_thirds);
%             groupednPFISNR
            
            %%
            figure(1); 
            if plotBarsorTimeseries==1
   %place Bar charts beneath SNR
            
                placeis = hzcounter + 6;
            
         mbar=(nanmean(groupednPFISNR,1));
         
         
         %adjust errorbars:
          %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            x = groupednPFISNR;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant)
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
        stB=nanstd(NEWdata,1)./sqrt(length(allppants));
        
%             stB=nanstd(groupednPFISNR,1)./sqrt(length(allppants));
        
subplot(4,3, placeis)

        b1=bar(mbar);        
         b1.FaceColor=col;
         hold on         
         errorbar(mbar, stB, 'color', 'k','LineStyle', 'none')
         
set(gca, 'fontsize', 15, 'xticklabels', {'0\leq1', '1\leq2', '2\leq3', '\geq3'})
xlabel('amount of PFI');
         ylabel('RESS log(SNR)')
         
         colormap('viridis')
ylim(ylimsare)            
shg
%             title([num2str(usehz) ' Hz (' num2str(hzis) 'f)'])
% axis tight
xlim([.5 4.5])
            counter=counter+1;
%%%%%%% ANOVA and LME stats if need be.             
            nppants=length(allppants);
            nconds = size(groupednPFISNR,2);
%             %RMANOVA?
            dataan= reshape(groupednPFISNR, [nppants*nconds,1]);
            conds = [ones(nppants,1); ones(nppants,1)*2;ones(nppants,1)*3;ones(nppants,1)*4];
                subs = repmat([1:length(allppants)]', [nconds,1]);
%             
%             [anret]=rmanova(dataan, conds, subs);
% %%
%             %LME to be sure:            
          tblA=table(dataan, conds, subs);
%%             %%
             UNterm1='dataan ~  conds + (1|subs)';
%              UNterm2='dataan ~  conds';
             UNterm2='dataan ~  1+ (1|subs)';
%              %creating models:
             lmeUN =fitlme(tblA, UNterm1);
             lmeUN2 =fitlme(tblA, UNterm2);
%              
    LMEout=compare(lmeUN2, lmeUN)
            %%
% %%
            else % plot time series instead.
%place Bar charts beneath SNR
            
                placeis = hzcounter+ 9;
            
            subplot(4,3,placeis);
                tbase = timeidDYN;
    lgsav=[];
% can also plot time-series:
for iPFI=1:4
    x= squeeze(groupednPFISNR_ts(:,iPFI,:));
    %adjust errorbars.
    
       
            mXppant =squeeze( nanmean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = nanmean(nanmean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
            %             %compute new stErr %which version?
            stE = nanstd(NEWdata)/sqrt(size(x,1));
            %
            hold on
            sh=shadedErrorBar(tbase, nanmean(x,1),stE,['-'],[1]);
                
                sh.mainLine.Color = col;
                sh.mainLine.LineWidth = iPFI*1.5;
                sh.patch.FaceColor = col;
                sh.edge(1).Color = col;
                sh.edge(2).Color = col;
                
   lgsav(iPFI)= sh.mainLine; 
end
%%
% legend([lgsav], [{'nPFI 0 \leq 1'}, {'nPFI 1 \leq 2'}, {'nPFI 2 \leq 3 or 4'}])
%%
set(gca,'fontsize', 15); set(gcf, 'color', 'w')
xlabel(['time from ' chtype])
ylabel('RESS log(SNR)')
axis tight
ylim([ylimsare(1) ylimsare(2)])

            end
%%
counter2=counter2+1;
    
            %%
           
       
        end 
        hzcounter=hzcounter+1;
    end
end
    %%
    shg
    %%
    
    colormap('viridis')
    cd(basefol)
    cd('Figures')
    cd('GFX PFI trial by trial')
    shg
     print('-dpng', ['PFI graded SNR summary graded, ' ctype '.png'])
    
    
end



if job.gradedchangesinERPimage_PFI_raincloudver==1
  %% now plot across ppants.
    % as above, with some tweaks.
    
    %%
    
    try
    % get nice colours from colorbrewer
    % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)
    [cb] = cbrewer('qual', 'Set3', 12, 'pchip');
catch
    % if you don't have colorbrewer, accept these far more boring colours
    cb = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4; 0.5 0.7 0.8; 1 0.8 0.5; 0.7 1 0.4; 1 0.7 1; 0.6 0.6 0.6; 0.7 0.5 0.7; 0.8 0.9 0.8; 1 1 0.4];
end

cl(1, :) = cb(4, :);
cl(2, :) = cb(1, :);
%%
    
    clf
    %can plot BARS next to trial-by-trial image, or snr-timecourse.
for     plotBarsorTimeseries=1:2 % plots bars beneath, then timeseries.

    
    getelocs
    rmvbase=0;
    
    counter2=1; % for plotting threeway PFI
    
    
    
    counter=1;
    hzcounter=1;
    
    %each Hz has a separate ylim range to optimize plots:
    ylimsareHz=[1.1,1.8;... %15 hz
             2.7,3.1;...%20 Hz
             .5, 1;... % 30 Hz
             2.8,3.3;... % 40 Hz
             0.4,.9;... % 45
             2.2,2.5;...
             .3, .6]; % 5 Hz
    
         hzlabel={'f1', 'f2', '2f1', '2f2', '3f1', '3f2*', 'f2-f1'};
         
         
    for hzis=1%:6%[1, 2, 7]
        cd(basefol)
        cd('EEG')
        cd('GFX_EEG-RESS')
        usehz= peakfreqsare(hzis);
        switch hzis
            case {1,3,5}
                                
                col='b';
                flickeris = 'Target';
            case {2,4}
%                 ylimsare=[2.5, 3.2];
                
                col=[.2 .2 .2];
                flickeris = 'Background';
            case 6
                col=[.2 1 .2];
                flickeris = 'Overlap';                
            case 7
%                 ylimsare=[2.5, 3.2];
                
                col='m';
                flickeris = 'Intermodulation';
        end
        
        ylimsare=ylimsareHz(hzis,:);
        
        load(['GFX_PFIperformance_withSNR_' num2str(usehz) '_min0_RESS'])
%         tgrm = timeidDYN;
        %%
        for itimezero=1%:2
            
            
            switch itimezero
                case 1
                    useBP=storeacrossPpant_onsetBP;
                    useSNR=storeacrossPpant_onsetSNR;
%                     useRTs=storeacrossPpant_onsetRTs;
%                     useTOPO=storeacrossPpant_onsetTOPO;
                    
                    ctype= 'target invisible';
                    chtype='button press';
                    
                case 2
                    
                    useBP=storeacrossPpant_offsetBP;
                    useSNR=storeacrossPpant_offsetSNR;
%                     useRTs=storeacrossPpant_offsetRTs;
%                     useTOPO=storeacrossPpant_offsetTOPO;
                    
                    ctype= 'target visible';
                    chtype='button release';
                    
                    
            end
            
            
            
            %                 for ip= 1:length(allppants)
            %                     trialind = [1:24] + 24*(ip-1);
            %                     acrBP(trialind,:,:) = useBP(ip,:,:);
            %                     acrSNR(trialind,:,:)=    useSNR(ip,:,:);
            %                     acrRTs(trialind)=    useRTs(ip,:);
            %
            %                 end
            
            %take mean across ppants for erp images..
            acrBP= squeeze(nanmean(useBP,1));
%             acrBP= squeeze(nanmean(useBP(21,:,:),1));
            
            acrSNR=squeeze(mean(useSNR,1));
%             acrSNR=squeeze(nanmean(useSNR(21,:,:),1));
            
            
            
            
            %sort across all
            %sort trial index by longest RT first.
            %%
            figure(1)
            hold on
            %%
            subplot(4,4,2:3)
            
            %         [sortedRTs, cid] = sort(acrRTs, 'descend');
            
            %
            imagesc([-3:1/60:3], 1:size(acrBP,1), acrBP);%(cid,:));
            %
            title({['PFI ' num2str(ctype)]}, 'fontsize', 15)
            c=colorbar;
            ylabel(c, 'buttons pressed')
            %                 ylabel('resampled catch trials')
            xlim([-1.75 1.75])
            set(gca, 'fontsize', 15, 'ytick', [])
            hold on
            plot([0 0] , ylim, ['w:'], 'linewidth', 3)
            plot([0 0] , ylim, ['k:'], 'linewidth', 2)
            %                 subplot(4,2,7)
            
            ppantMeanBP= squeeze(mean(useBP,2));
            %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            x = ppantMeanBP;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
            %compute new stErr %which version?
            %             stE = std(NEWdata)/sqrt(size(x,1));
            %                 shadedErrorBar([-3:1/60:3], mean(ppantMeanBP,1),stE)
            %                 set(gca, 'fontsize', 15)
            %                 hold on;
                            ylabel({['Normalized'];['trial count']})
            %                 plot([0 0] , ylim, ['k:'], 'linewidth', 1)
%             axis square
            %%
%             xlabel(['Time from ' ctype])
            
            %SNR
                placeis = hzcounter + 3;

            subplot(4,3, placeis);
            
            %reorder SNR
            acrSNRsort=acrSNR;
            
            %%
            %                 %smooth across trials
            sm_snrgrm20=zeros(size(acrSNRsort));
            for itime=1:size(acrSNRsort,2)
                sm_snrgrm20(:,itime)= smooth(acrSNRsort(:,itime),15);
            end
            
            
            %     imagesc(acrSNR);
            imagesc(timeidDYN, 1:size(acrSNRsort,1), acrSNRsort);
            c=colorbar;

%             ylabel(c, {['RESS log(SNR)']})

if hzis~=7    
            title([flickeris ' ' num2str(usehz) ' Hz (' hzlabel{hzis} ')'])
else
            title([flickeris ' ' num2str(usehz) ' Hz (f2-f1)'])
end
%             title(['BG ' num2str(hzis) 'F'])
            
            hold on
            plot([0 0] , ylim, ['k:'], 'linewidth', 3)
            plot([0 0] , ylim, ['w:'], 'linewidth', 2)
            
            
            xlim([-1.75 1.75])
            set(gca, 'fontsize', 15, 'ytick', [])
            %                 subplot(4,2,8)
            %
            %plot across ppant trace
            ppantMeanSNR= squeeze(mean(useSNR,2));
            
            %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            x = ppantMeanSNR;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
        
            xlabel(['Time from ' ctype])
            ylabel({'Normalized'; 'trial count'})
                caxis([ylimsare])
%                 axis square
                
            set(gca, 'fontsize', 15)
            cd(basefol)
            set(gcf, 'color', 'w')
            shg
            if hzis==7;
                ylabel(c, 'RESS log(SNR)')
            end
            
            % beside this, plot the bar!
            %grouped nPFI
            %%
            groupednPFISNR=nan(length(allppants), 3);  % for barchart
            groupednPFISNR_ts=nan(length(allppants), 3, size(acrSNR,2));  % for SNR time-seriestrace
            if itimezero==1 % interested in after BP                
            BPidx=181:361;
            SNRidx= 8:14;
            else
                BPidx=1:180;
                SNRidx= 1:7;
            end
                groupedBEHflpd = fliplr(squeeze(nanmean(acrBP(:,BPidx),2))'); %takes mean over all SNR timepoints.
            %%
            groupednPFISNR_thirds=nan(length(allppants),3);
            
            for ippant = 1:size(useBP,1); 
%                % per ppant, collect trial indices:
                ppant_BEH= squeeze(mean(useBP(ippant,:,BPidx),3)); %takes mean over all BP time points.
                ppant_SNR=squeeze(mean(useSNR(ippant,:,SNRidx),3)); %takes mean over only 3 second window.
%                 ppant_SNR=squeeze(mean(useSNR(ippant,:,:),3)); %takes mean over all SNR timepoints.
                ppant_SNR_ts=squeeze(useSNR(ippant,:,:)); %takes mean over all SNR timepoints.
               maxnPFI=ceil(max(ppant_BEH));
%                
               %flip BEH /SNR to correct orientation
               ppantBEHflpd=fliplr(ppant_BEH);
               ppantSNRflpd=fliplr(ppant_SNR);
               ppant_SNR_ts=flipud(ppant_SNR_ts);
               % only collect bin width.
%              
%                 nPFIidx = dsearchn(ppantBEHflpd', [0:3]');
               nPFIidx = dsearchn(groupedBEHflpd', [0:4]');
               
               nPFIidx=unique(nPFIidx);
               
               if nPFIidx(1) ~=1
                   nPFIidx(1)=1;
               end
               
               %take mean for these trials (actually rows in trialxtrial)
               for iPFI=1:length(nPFIidx)-1;
                   
                   try indices = nPFIidx(iPFI):nPFIidx(iPFI+1)-1;
                   catch %if max
                       indices = nPFIidx(iPFI):100;
                   end
                   
                   if iPFI==4;
                       shg
                   end
                   groupednPFISNR(ippant, iPFI)= squeeze(nanmean(ppantSNRflpd(indices)));
                   
                   
                   groupednPFISNR_ts(ippant, iPFI,:)= squeeze(nanmean(ppant_SNR_ts(indices,:)));
% %                    %
%                    figure(3)
%                    clf
%                    plot(ppantSNRflpd);
%                    hold on
%                    plot([indices(1) indices(1)], ylim, ['k:'])
%                    plot([indices(end) indices(end)], ylim, ['k:'])
%                    title([num2str(squeeze(nanmean(ppantSNRflpd(indices))))])
               end
               
% works well.  

groupednPFISNR_quarters(ippant,1)= squeeze(nanmean(ppant_SNR(1:24)));
groupednPFISNR_quarters(ippant,2)= squeeze(nanmean(ppant_SNR(25:50)));
groupednPFISNR_quarters(ippant,3)= squeeze(nanmean(ppant_SNR(51:75)));
groupednPFISNR_quarters(ippant,4)= squeeze(nanmean(ppant_SNR(76:100)));

%                
            end
            %% using raw thirds?
%             groupednPFISNR_thirds=fliplr(groupednPFISNR_thirds);
%                         groupednPFISNR=fliplr(groupednPFISNR_thirds);
%             groupednPFISNR
            
            %%
            figure(1); 
            if plotBarsorTimeseries==1
   %place Bar charts beneath SNR
            
                placeis = hzcounter + 6;
            
         %create structure for raincloudplots:
         d=[];
         for icat=1:size(groupednPFISNR,2);
             d{icat}=groupednPFISNR(:,icat);
         end
           %% 
         rm_raincloud(d, cb(1:4,:))
         %%
         
         
         
         
         %adjust errorbars:
          %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            x = groupednPFISNR;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant)
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
        stB=nanstd(NEWdata,1)./sqrt(length(allppants));
        
%             stB=nanstd(groupednPFISNR,1)./sqrt(length(allppants));
        
subplot(4,3, placeis)

        b1=bar(mbar);        
         b1.FaceColor=col;
         hold on         
         errorbar(mbar, stB, 'color', 'k','LineStyle', 'none')
         
set(gca, 'fontsize', 15, 'xticklabels', {'0\leq1', '1\leq2', '2\leq3', '\geq3'})
xlabel('amount of PFI');
         ylabel('RESS log(SNR)')
         
         colormap('viridis')
ylim(ylimsare)            
shg
%             title([num2str(usehz) ' Hz (' num2str(hzis) 'f)'])
% axis tight
xlim([.5 4.5])
            counter=counter+1;
%%%%%%% ANOVA and LME stats if need be.             
            nppants=length(allppants);
            nconds = size(groupednPFISNR,2);
%             %RMANOVA?
            dataan= reshape(groupednPFISNR, [nppants*nconds,1]);
            conds = [ones(nppants,1); ones(nppants,1)*2;ones(nppants,1)*3;ones(nppants,1)*4];
                subs = repmat([1:length(allppants)]', [nconds,1]);
%             
%             [anret]=rmanova(dataan, conds, subs);
% %%
%             %LME to be sure:            
          tblA=table(dataan, conds, subs);
%%             %%
             UNterm1='dataan ~  conds + (1|subs)';
%              UNterm2='dataan ~  conds';
             UNterm2='dataan ~  1+ (1|subs)';
%              %creating models:
             lmeUN =fitlme(tblA, UNterm1);
             lmeUN2 =fitlme(tblA, UNterm2);
%              
    LMEout=compare(lmeUN2, lmeUN)
            %%
% %%
            else % plot time series instead.
%place Bar charts beneath SNR
            
                placeis = hzcounter+ 9;
            
            subplot(4,3,placeis);
                tbase = timeidDYN;
    lgsav=[];
% can also plot time-series:
for iPFI=1:4
    x= squeeze(groupednPFISNR_ts(:,iPFI,:));
    %adjust errorbars.
    
       
            mXppant =squeeze( nanmean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = nanmean(nanmean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
            %             %compute new stErr %which version?
            stE = nanstd(NEWdata)/sqrt(size(x,1));
            %
            hold on
            sh=shadedErrorBar(tbase, nanmean(x,1),stE,['-'],[1]);
                
                sh.mainLine.Color = col;
                sh.mainLine.LineWidth = iPFI*1.5;
                sh.patch.FaceColor = col;
                sh.edge(1).Color = col;
                sh.edge(2).Color = col;
                
   lgsav(iPFI)= sh.mainLine; 
end
%%
% legend([lgsav], [{'nPFI 0 \leq 1'}, {'nPFI 1 \leq 2'}, {'nPFI 2 \leq 3 or 4'}])
%%
set(gca,'fontsize', 15); set(gcf, 'color', 'w')
xlabel(['time from ' chtype])
ylabel('RESS log(SNR)')
axis tight
ylim([ylimsare(1) ylimsare(2)])

            end
%%
counter2=counter2+1;
    
            %%
           
       
        end 
        hzcounter=hzcounter+1;
    end
end
    %%
    shg
    %%
    
    colormap('viridis')
    cd(basefol)
    cd('Figures')
    cd('GFX PFI trial by trial')
    shg
     print('-dpng', ['PFI graded SNR summary graded, ' ctype '.png'])
    
    
end