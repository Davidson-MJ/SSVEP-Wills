% s3_CE_computeRAWSNR_TargandCatch

%having epoched tgs and catch, pre RESS. apply FFT and SNR to windows when
%all targets are present/ away from Catch stimulus onset.
clear all
addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/EEG')
basefol=pwd;

%%


job.calcppantSTATICperfreq=1; %sorts by hz x location.
% do we want to reject 'bad' trials?
rejtrials = 0;


job.concatTOPOacrossppanst=1;


%%% note above needs to be rerun. 

job.plotTOPOsacrossppants=1;

job.plotHzsacrossppants=1;




getelocs;

cd(basefol)

dirs= dir([pwd filesep '*_*' 'EEG']);
%% load data to determine physical catch timing.

allppants=[1,2,4,6,9:16,18]; %
% window=[-2 4];
window=[-3 3];

srate=250;

epochdur = sum(abs(window));

timeid = [0:1/srate:epochdur];
timeid= timeid-3;

onsetc = ceil(epochdur)/2;
% peakfreqsare=[20,40]; %hz
%timing
tt = 0:1/srate:60;
ntrials=48;


if job.calcppantSTATICperfreq==1
    %Collect all DATA
    
    exclTransientBP=0; %don't need to sort by dur this analysis.
    
    %%
    windowsmall=[];
    windowsmall(1,:) = [-3 -.1] ;%  window targets present
    windowsmall(2,:) = [.1 3] ;%  window targets present after BP
    %%
    %     usechans=[51:64, 17:32]; %occipital
    usechans=[1:64]; %occipital
    
    param_spctrm.tapers = [1 1];
    param_spctrm.Fs= [250];
    param_spctrm.fpass= [0 50];
    param_spctrm.trialave=0;
    for ifol =allppants
        
        
        %%
        cd([basefol])

cd(dirs(ifol).name)        
        
        %% load the relevant PFI data.
        
        load('ppant_PFI_Epoched');
        load('ppant_Catch_Epoched')
        load('TrialIndicesbyCatchLocationandNum.mat')
        peakfreqsare=[15,20,30,40,5,35];
        %%
        
%         rawSNR_static_byBGharms=zeros(6,4,64,513);
        rawSNR_static_byHzxLoc=zeros(length(peakfreqsare),4,64,205);
        rawSNR_topo_byHzxLoc=zeros(length(peakfreqsare),4,64);
        rawSPEC_topo_byHzxLoc=zeros(length(peakfreqsare),4,64);
        
                
                % collect relevant trials for each type of spatial
                % configuration/filter construction.
                
                %%
                for id=1:12% Use all epochs so as not to bias condition comparisons.
                    
                    switch id
                        case 1
                            dataIN=ppant_SNREEG_PFI_0_1;
                            BPscheck= BPs0_1;
                            usewindow=1:2;
                        case 2
                            dataIN=ppant_SNREEG_PFI_1_0;
                            BPscheck= BPs1_0;
                            usewindow=1:2;
                        case 3
                            dataIN=ppant_SNREEG_PFI_1_2;
                            BPscheck= BPs1_2;
                            usewindow=1:2;
                        case 4
                            dataIN=ppant_SNREEG_PFI_2_1;
                            
                            BPscheck= BPs2_1;
                            usewindow=1:2;
                        case 5
                            dataIN=ppant_SNREEG_PFI_2_3;
                            
                            BPscheck= BPs2_3;
                            usewindow=1:2;
                        case 6
                            dataIN=ppant_SNREEG_PFI_3_2;
                            
                            BPscheck= BPs3_2;
                            usewindow=1:2;
                            
                            %%%% NEW
                        case 7
                             dataIN=ppant_SNREEG_PFI_3_4;
                            
                            BPscheck= BPs3_4;
                            usewindow=1:2;
                        case 8
                             dataIN=ppant_SNREEG_PFI_4_3;
                            
                            BPscheck= BPs4_3;
                            usewindow=1:2;
                        case 9
                            dataIN=ppant_SNREEG_disapBPwithincatch;
                            
                            usewindow=1;
                        case 10
                            dataIN=ppant_SNREEG_reapBPaftercatch;
                            
                            usewindow=2;
                        case 11
                            dataIN=ppant_SNREEG_catchramponset;
                            
                            usewindow=1;
                        case 12
                            dataIN=ppant_SNREEG_catchrampoffset;
                            usewindow=2;
                            
                            
                            
                    end
                    
                    % which trials can we look into?
                        
                        relevanttrials=1:ntrials;
                    
                        
                        if id<9 % check the duration of disap was not a transient
                            %onset or offset window?
                            switch mod(id,2)
                                case 0 %even numbers (in id) are reappearing targets
                                    % find duration of pre zero bp. (and
                                    % flip)
                                    BPsshort= fliplr(BPscheck(:,1:180));
                                    
                                case 1 %odd numbers are disappearing target types
                                    %
                                    BPsshort = BPscheck(:,181:end);                                    
                            end
                            %first difference in 2nd dimension.
                                    tmp=(diff(BPsshort,1,2));
                                    dursCheck=zeros(size(BPsshort,1),1);
                                    for itrial=1:size(BPsshort,1)
                                        try dursCheck(itrial,:)= find(squeeze(tmp(itrial,:)),1, 'first');
                                        catch % the button press didn't end!
                                            dursCheck(itrial,:) = 180; %ie at least epoch length.
                                        end
                                    end
                                    
                        %retain which trials?                                    
                                    trialtypesperTargPresent = find(dursCheck>exclTransientBP);

                        else
                            trialtypesperTargPresent= 1:48; %all trials for pre/post catch periods.
                        end
                            
                  
                        
                        
                    % now that we have a datatype, get the right epochs that share stimulus configuration
                    
                    % also correct window (pre or post indication of target
                    % presence)
                    for windch=1:length(usewindow)
                        wind=usewindow(windch);
                        windcheck= windowsmall(wind,:);
                        
                        tidx=dsearchn(timeid', [windcheck]');
                        
                        
                        %reduce size.
                        datast= squeeze(dataIN(trialtypesperTargPresent,:,tidx(1):tidx(2)));
                        
                        
                        %% check for bad trials (noisy)
                        %std per trial(average over electrodes)
                        tmp=[];
                        if ndims(datast)<3
                            tmp(1,:,:)= datast;
                            datast=tmp;
                        end
                        datastSD = nanstd(squeeze(nanmean(datast,2)),0,2);
                        
                        if rejtrials==1
                        %remove those with 2.5*std from mean.
                        trialSD=nanstd(datastSD);
                        mSD=nanmean(datastSD);
                        keeptrials=1:size(datast,1);
                        
                        %remove the trials with excessive noise.
                        badtrials=find(datastSD>mSD+2.5*trialSD)';
                        
                        % also skip the trials which were only transient button
                        % presses. (less than one second).
                        shorttrials = find(durscheck<60);
                        
                        
                        %also remove NANs:
                        nantrials = find(isnan(datast(:,1,1)));
                        
                        badtrials = [badtrials, shorttrials, nantrials'];
                        
                        
                        % remove these from consideration.
                        keeptrials(badtrials)=[];
                        datast=datast(keeptrials,:,:);
                        
                        else %just remove nans prior to analysis:
                            keeptrials=1:size(datast,1);
                        
                        %also remove NANs:
                        nantrials = find(isnan(datast(:,1,end)));
                        keeptrials(nantrials)=[];
                        datast=datast(keeptrials,:,:);                        
                        end
                            
                            
                        
                        % concatenate all these types to increase quality of cov matrix.
                        if id==1 && windch==1%start here
                            if ndims(datast)<3
                                dataOUT(1,:,:)=datast;
                            else
                                dataOUT=datast;
                            end
                            
                        else
                            % append
                            if ndims(datast)<3
                                tmp(1,:,:)=datast;
                                dataOUT=cat(1,dataOUT,tmp);
                                
                            else
                                dataOUT=cat(1,dataOUT,datast);
                            end
                            
                        end
                        
                    end
                    
                end
                
                %now we have ALL the data this freqxLoc.
                %%
                sNOW= zeros(64,205); %size of spectrum.
%                 sNOW=zeros(64,513);
                SpecNOW=sNOW; %save raw (preSNR) also
                for ichan= 1:64
                    tmp= squeeze(dataOUT(:,ichan,:));
                    
                    [s,f]=mtspectrumc(tmp', param_spctrm)    ;
                    
                    %% comp SNR;
                    kernelw = [-1/8 -1/8 -1/8 -1/8 0 0  1  0 0 -1/8 -1/8 -1/8 -1/8];
                    %compute SNR at freq of interest.                    
                    stmp=zeros(size(s));
                    for tmptr=1:size(s,2) 
                        
                    stmp(:,tmptr)=conv(log(s(:,tmptr))', kernelw, 'same');
                    end
                    
                    sNOW(ichan,:)= squeeze(mean(stmp,2)); %mean over trials.
                    
                    SpecNOW(ichan,:) = squeeze(mean(log(s),2));
                end
            
            %% sanity check
%             usehz=20;
% % %             
%             [~, hzid] = min(abs(f-usehz));
% % 
%             figure(1)
%             clf
% %             subplot(4,2,1:2)
%             topoplot(sNOW(:,hzid), elocs(1:64), 'emarker2', {62, 'o', 'w'});
%             c=colorbar;
%             caxis([0 2])
% %             ylabel(c, 'SNR')
% %             title({[num2str(usehz) 'Hz SNR '];['-3000:-100ms']})
% % 
% %             subplot(4,2,3:4)
%             plot(f, sNOW(62,:)); hold on; 
%             plot(f, SpecNOW(1,:))
% %             shg
% % %             %%
%             
%%            
            rawSNR_static = sNOW; %channs x hz
%           
            rawSPEC_static= SpecNOW; %collapsed by freq
        
        
        
        save('Raw_statSNR', 'rawSNR_static', 'rawSPEC_static', 'f', 'kernelw')
%         save('Raw_statSNR_no_ds', 'rawSNR_static_byBGharms')
        disp(['fin ' num2str(ifol)])
    end
    
end


if job.concatTOPOacrossppanst==1
    
    %%
    acrossRAWSPEC_freqs=zeros(length(allppants),64,205);
    acrossRAWSNR_freqs=zeros(length(allppants),64,205);



    icounter=1;
    for ippant=allppants
        cd([basefol])

        cd(dirs(ippant).name)        
       
        load('Raw_statSNR')
        acrossRAWSPEC_freqs(icounter,:,:)=rawSPEC_static;
        acrossRAWSNR_freqs(icounter,:,:)=rawSNR_static;
        icounter=icounter+1;
    end
    %
    cd(basefol)
    
    save('GFX_rawSNR_static', 'acrossRAWSPEC_freqs', 'acrossRAWSNR_freqs','f')
    
end
%%

if job.plotTOPOsacrossppants==1
    %%
    cd(basefol)    
    load('GFX_rawSNR_static.mat');
    clf
    meanTOPOs=squeeze(mean(acrossRAWSNR_freqs,1));
    figure(1);
    clf
%     titlesare={'15 Hz' , '13 Hz', '15 Hz', '18 Hz'};
    usefreqs=[15,20,30,40,5,35];
    for ipl=1:6
        subplot(2,3,ipl)
        
        [~,usef] = min(abs(f-usefreqs(ipl)));
        topoplot(meanTOPOs(:,usef), elocs(1:64));
%         topoplot(squeeze(nanmean(meanTOPOs(ipl,:,:),2)), elocs(1:64));
        c=colorbar;
        title([num2str(usefreqs(ipl)) ' Hz SNR'])
        ylabel(c, 'log(SNR)')
        caxis([0 2])
        set(gca, 'fontsize', 25)
    end
    set(gcf, 'color', 'w')
    %%
    cd(basefol)
    cd ../
    cd('Figures')
    cd('Topoplots')
    print('-dpng', ['Topos across all, SNR preRESS'])
    %%
    clf
    for ifreq=usefreqs
    [~,usef] = min(abs(f-ifreq));
    
    for ippant = 1:size(acrossRAWSNR_freqs,1)
        
        subplot(3,5,ippant)
        title(['# ' num2str(allppants(ippant))])
                topoplot(squeeze(acrossRAWSNR_freqs(ippant,:,usef)), elocs(1:64));
                c=colorbar;
                caxis([0 2]);
                set(gca, 'fontsize', 25)
       
    end
             set(gcf, 'color', 'w')
             cd(basefol)
             cd ../
             cd('Figures')
             cd('Topoplots')
             print('-dpng', ['Topos individual ' num2str(ifreq) ' Hz SNR preRESS'])
    end


end



if job.plotHzsacrossppants==1
    %%
    cd(basefol)
    
%     load('GFX_rawSNR_static.mat');
    load('GFX_rawSNR_static.mat');
    clf
    usechan=30;
    %mean over ppants for plot.
    meanHz=squeeze(nanmean(acrossRAWSNR_freqs(:,usechan,:),1));       
    %mean over locs.

    
    plot(f, meanHz, 'k', 'linewidth', 3)
    %
    hold on
    cols={'r', 'k', 'r', 'k', 'm','m', 'r', 'r', 'k','m'};
    counter=1;
    p=[];
        for ifreq=[15,20,30,40,5,35, 45]
            [~,usef]= min(abs(f-ifreq));
            yis= meanHz(usef);
            xis= f(usef);
           p(counter)=plot(xis, yis, 'o','markersize', 25, 'linew', 15, 'color', cols{counter}) ;
           
           if ifreq==60
               p(counter)=plot(xis, yis, 'o', 'markersize', 25,'linew', 2, 'color', 'k') ;
           end
            counter=counter+1;
            
        end
        legend([p(1), p(2), p(5)], {'TG Hz', 'BG Hz', 'im'})
    ylim([0 2])
    ylim([-1 2])
    ylabel('logSNR')
    xlabel('Frequency (Hz)')
    set(gca, 'fontsize', 25)
    cd(basefol)
    cd ../
    cd('Figures')
    cd('SNR spectrum')
    print('-dpng','short window SNR, chan POz')
    %% for testing sig.
    %mean over locs.
    dataT=squeeze(mean(acrossRAWSNR_freqs(:,1,:,usechan,:),3));
    
    hold on 
    %plot stimulus SSVEP markers   
   
    plcount=1;
    
    %place details:
    p=[];
    
    hztocheck=[8,13,15,18,20,40, 16,26,30,36];
    for ihz=hztocheck
        [~,hzid]=min(abs(f-ihz));
    
    
    if plcount<6
        col='r';
    else
        col=[ 0 .5 0];
    end
    
    
    %check SIG
    [~, p(plcount)] = ttest(dataT(:,hzid), 0);
    
    
%     ps=plot(f(hzid), plotme(hzid), ['o' col], 'linewidth', 10);
    
    if mod(ihz,20) ==0
        typel='--';
    else
        typel ='-';
    end
        
    st(plcount)=plot([f(hzid) f(hzid)], [-10 plotme(hzid)] , [typel col], 'linewidth', 2);
    st(plcount).Color = col;
    plcount=plcount+1;
    end
    
    
    % now plot after correcting for MCP.
    q= .05/length(hztocheck);
    pcount=1;
    for ihz=hztocheck
        [~,hzid]=min(abs(f-ihz));
        
        if p(pcount)<q
%            tt= text(f(hzid)-.55, plotme(hzid)+.1 , '*', 'Color', 'k', 'FontSize', 45);
           tsig= plot(f(hzid), plotme(hzid)+.4 , ['*k'], 'markers', 15,'linew', 2);
        end
        pcount = pcount+1;
    end
    %
    
    
        ylabel( 'log(SNR)')        
        title(['log(SNR) at channel ' num2str(elocs(usechan).labels)])
        set(gca, 'fontsize', 25)
        
        ylim([-1.5 5])
    xlim([1 45])
    xlabel('Hz')
    lg=legend([st(1),st(5), st(7), st(6), tsig],{'Targets 1st harmonic', 'Background 1st harm.', 'Targets 2nd harmonic', 'Background 2nd harm.', ' \itp\rm_F_D_R < .05'});
    set(gcf, 'color', 'w')
    %%
    cd('Figures')
    %
    print('-dpng', ['Raw static SNR Hz at ' num2str(elocs(usechan).labels) '.png'])
   
end
