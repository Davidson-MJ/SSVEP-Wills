% Follows the format of s3_CE... only now applying to RESS timeseries

%having epoched tgs and catch, POST RESS. apply FFT and SNR to windows when
%all targets are present/ away from Catch stimulus onset.


clear all
addpath('/Users/MattDavidson/Desktop/SSVEP-feedbackproject')
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/EEG')
basefol=pwd;
pdirs = dir([pwd filesep '* EEG']);
%%

job.calcppantSNRperfreq=0; %sorts by hz x location. %topos in preivous script (s3_Da)

job.concatacrossppanst=1;
job.plotHzsacrossppants=1; % this is averaging across all periods (excluding catch).

job.compareHzxLocations_bar=0;

getelocs;

cd(basefol)


peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; % don't change!   

%remaining participants after behavioral data analysis/exclusion
allppants=[1,2,4,6,9:16,18]; %
%% load data to determine physical catch timing.

%remaining participants after behavioral data analysis/exclusion

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



if job.calcppantSNRperfreq==1
    %Collect all DATA
    %%
    windowsmall=[];
    windowsmall(1,:) = [-3 -.1] ;%  window targets present
    windowsmall(2,:) = [.1 3] ;%  window targets present after BP
%     windowsmall(1,:) = [-3 -1] ;%  window targets present
%     windowsmall(2,:) = [1 3] ;%  window targets present after BP
    %%
    
    
    param_spctrm.tapers = [1 1];
    param_spctrm.Fs= [250];
    param_spctrm.fpass= [0 65];
    param_spctrm.trialave=0;
    param_spctrm.pad=2; %padding for detecting 5 Hz.
    
    for ifol =allppants
        
        
        %%
        cd(basefol)
        cd(pdirs(ifol).name)
        
        
        %% load the relevant PFI data.
        load('ppant_PFI_Epoched_RESS');
        load('ppant_Catch_Epoched_RESS')
        load('ppant_PFI_Epoched') %for info of types missing.
%         load('TrialIndicesbyLocationandHz.mat')
        
        %%
        
        RESS_staticSNR_byHz=zeros(length(peakfreqsare),1065);


        
        for ifreq=1:length(peakfreqsare)
            
            
         usehz=peakfreqsare(ifreq);
            
            
                
                
                
                % collect relevant trials for each type of spatial
                % configuration/filter construction.
                
                %%
                for id=1:14% Use all epochs so as not to bias condition comparisons.
                    
                    switch id
                        case 1
                            dataIN=ress_PFI_0_1_Hz;
                            
                            usewindow=1:2;
                        case 2
                            dataIN=ress_PFI_1_0_Hz;
                            
                            usewindow=1:2;
                        case 3
                            dataIN=ress_PFI_1_2_Hz;
                            
                            usewindow=1:2;
                        case 4
                            dataIN=ress_PFI_2_1_Hz;
                            
                            usewindow=1:2;
                        case 5
                            dataIN=ress_PFI_3_2_Hz;
                            
                            usewindow=1:2;
                        case 6
                            dataIN=ress_PFI_2_3_Hz;
                            
                            usewindow=1:2;                                                       
                            
                            
                        case 7
                            
                            dataIN=ress_PFI_3_4_Hz;
                            
                            usewindow=1:2;                                                       
                        case 8
                            
                            dataIN=ress_PFI_4_3_Hz;
                            
                            usewindow=1:2;                                                       
                            
                        case 9
                            dataIN=ress_BPcatchonsetTGs;
                            
                            usewindow=1;
                        case 10
                            dataIN=ress_BPcatchonsetBGs;
                            
                            usewindow=1;
                        case 11
                            dataIN=ress_catchonsetTGs;
                            
                            usewindow=1;
                        case 12
                            dataIN=ress_catchonsetBGs;
                            
                            usewindow=1;
                        case 13
                            dataIN=ress_catchoffsetTGs;
                            
                            usewindow=2;
                            
                        case 14
                            dataIN=ress_catchoffsetBGs;                            
                            usewindow=2;
                            
                        case 15
                            dataIN=ress_BPcatchoffsetTGs;                            
                            usewindow=2;
                        case 16
                            dataIN=ress_BPcatchoffsetBGs;                            
                            usewindow=2;
                            
                            
                    end
                    
              
                    % now that we have a datatype, get the right epochs that share stimulus configuration
                    
                    % also correct window (pre or post indication of target
                    % presence)
                    for windch=1:length(usewindow)
                        wind=usewindow(windch);
                        windcheck= windowsmall(wind,:);
                        
                        tidx=dsearchn(timeid', [windcheck]');
                        
                        
                        %reduce size. 
                        if id<9
                            datanow=squeeze(dataIN(ifreq).ressTS);
                            
                        datast= squeeze(datanow(:,tidx(1):tidx(2)));
                        else %catch trials
                            if ifreq<7 %TG or BG
                                
                                if mod(ifreq,2)==0 %BG hz
                                    usehztmp = ifreq/2;
                                elseif mod(ifreq,2)~=0% TG hz
                                    usehztmp= ceil(ifreq/2);
                                end
                            elseif ifreq>=7
                                usehztmp = ifreq-6;
                            end
                            datast= squeeze(dataIN(usehztmp,:,tidx(1):tidx(2)));
                            
                        end
                        
                        %% check for bad trials (noisy)
                        %std per trial(average over timepoints)
%                         
%                         datastSD = nanstd(datast,0,2);
%                         
%                         %remove those with 2.5*std from mean.
%                         trialSD=nanstd(datastSD);
%                         mSD=nanmean(datastSD);
keeptrials=1:size(datast,1);
%
%                         %remove the trials with excessive noise.
badtrials=[];
%                         badtrials=find(datastSD>mSD+2.5*trialSD)';
%
%                         % also skip the trials which were only transient button
%                         % presses. (less than one second).
shorttrials=[];
%
%                         %also remove NANs:
nantrials = find(isnan(datast(:,1)));
zerotrials= find(datast(:,1)==0);
badtrials = [badtrials, shorttrials, nantrials', zerotrials'];
%
%
%                         % remove these from consideration.
if sum(badtrials)>0
    keeptrials(badtrials)=[];
    datast=datast(keeptrials,:);
end
                        %%
%                         %%
                        
                        if id==1 && windch==1%start here
                            if ndims(datast)<2
                                dataOUT(1,:,:)=datast;
                            else
                                dataOUT=datast;
                            end
                            
                        else
                            % append
                            if ndims(datast)<2
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
                
                %
                    
                    [s,f]=mtspectrumc(dataOUT', param_spctrm)    ;
                    
                    %TAKE mean
                    sMEAN = squeeze(nanmean(log(s),2));
                    
                    %% comp SNR;
                    kernelw = [-1/24 -1/24 -1/24 -1/24 -1/24 -1/24 -1/24 -1/24 -1/24 -1/24 -1/24 -1/24   0 0 1  0 0 -1/24 -1/24 -1/24 -1/24 -1/24 -1/24 -1/24 -1/24 -1/24 -1/24 -1/24 -1/24 ];
%                     kernelw = [-1/6 -1/6 -1/6   0 0 1  0 0 -1/6  -1/6 -1/6];
%                     kernelw = [-1/4 -1/4 0 0 1  0 0 -1/4  -1/4 ];
                    %                 
%                     %take mean
                    sNOW= conv(sMEAN, kernelw, 'same');
                %
%                     figure(1); clf
%                     hold on;
%                     plot(f,sNOW);
%             
            %%
            
            
            RESS_staticSNR_byHz(ifreq,:) = sNOW; %channs x hz
            
            
        end
        
        
        save('RESS_statSNR', 'RESS_staticSNR_byHz', 'f', 'kernelw')
    end
    
end


if job.concatacrossppanst==1
    
    %%
    
    acrossRESS_SNR_freqs=zeros(length(allppants), length(peakfreqsare),1065);
     acrossRESS_SNR = zeros(length(allppants), length(peakfreqsare));
    
    icounter=1;
    
    for ippant=allppants
        cd(basefol)
    cd(pdirs(ippant).name)
        
        load('RESS_statSNR')
        
        acrossRESS_SNR_freqs(icounter,:,:)=RESS_staticSNR_byHz;
        
        
        for ifreq= 1:length(peakfreqsare)
            [~, id]= min(abs(f-peakfreqsare(ifreq)));
            
            acrossRESS_SNR(icounter,ifreq)=squeeze(RESS_staticSNR_byHz(ifreq,id));
            
        end
        icounter=icounter+1;
    end
    %
    cd(basefol)
    cd('GFX_EEG-RESS')
    save('GFX_ressSNR_static', 'acrossRESS_SNR_freqs','f', 'acrossRESS_SNR', 'peakfreqsare')
    
end



if job.plotHzsacrossppants==1
    %%
    cd(basefol)
    cd('GFX_EEG-RESS')
    
      load('GFX_ressSNR_static.mat'); 
      
    %%
    figure(1)
    clf
    
    cols(5).c='r';
    cols(6).c='r';
    
    
    ressd=[];

% peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; % don't change!   


    %%
    clf
    lc={};
    for ipl=1:3 % TGs, BGs, IMs
        switch ipl
            case 1
                usehz=1;%[1,3,5];
                col= 'b';
                pname = 'Target SSVEPs';
            case 2
                usehz=[2];
                
                col= [0 .5 0];
                pname = 'Background SSVEPs';
                
            case 3
                usehz=7;%:9;
                
                col= ['k'];
                
                pname = 'Inter-modulaion components';
%                 clf
        end
        pcounter=1;
    for ihz=usehz
        
%         subplot(2,2,pcounter)
        hold on
        % first plot the RESS result        
    meanHz=(squeeze(nanmean(acrossRESS_SNR_freqs(:,ihz,:),1)));
    
    
     ressd(ihz)=plot(f, meanHz, 'color', col, 'linewidth', 2);
    hold on
    
%     xlim([ 5 20])   
%     ylim([ -.25 .5])
%     
%     if  ihz>7
%         xlim([20 40])
%     elseif ipl==3
%         xlim([peakfreqsare(ihz)-5 peakfreqsare(ihz)+5])
%         ylim([ -1.5 5]) 
%     end
    ylabel('RESS log(SNR)')
    xlabel('Hz')
%     then plot the traditional
 
 
    pcounter=pcounter+1;

lc = [lc {pname}];

pvals = zeros(1,size(acrossRESS_SNR_freqs,3));
for it=1:size(acrossRESS_SNR_freqs,3)
    
    
[~,pvals(it)]= ttest(squeeze(acrossRESS_SNR_freqs(:,ihz,it)),0);
end

[qsig, qmask] = fdr(pvals,.05);

sigpoints=find(pvals<=qsig);

% for isig= 1:length(sigpoints)
% text(f(sigpoints(isig)), 4, '*')
% end

    end
    
    
       legend(lc)
       
        
       
    set(gcf, 'color', 'w')
    %
    ylim([-1 3])
    xlim([0 65])
    set(gca, 'fontsize', 25)
%       print('-dpng', ['RESS target SNR, ' pname])
    end
     %%
  
   
   %% construct hz x loc for Raw
% now BAR equivalent to test for sig:
clf
acrossRAWSNR =[];
for ihz = 1:6
    for iloc = 1:4

 hzid = peakfreqsare(ihz);
        [~,hzid]=min(abs(f-hzid));
        
        acrossRAWSNR(:,ihz,iloc)=squeeze(acrossRAWSNR_freqs(:,ihz,iloc,usechan,hzid));
    end
end

%%
% construct data:
mRAW = squeeze(nanmean(acrossRAWSNR,1));
mRESS= squeeze(nanmean(acrossRESS_SNR,1));
tmpB=[];
tmpB(1,:,:)= mRAW;
tmpB(2,:,:)= mRESS;
clf
subplot(211)
bar(mRAW);% hold on; bar(mRESS)
subplot(212)
bar(mRESS)
shg
end

if job.compareHzxLocations_bar==1
 cd(basefol)
    cd('newplots-MD')
    %load raw for comparison
%     load('GFX_rawSNR_static.mat');
      load('GFX_ressSNR_static.mat'); 
      
    %%
    usefreqs=[1:4, 7:10]; %all TGs
% usefreqs=[1:4];

   storeall= squeeze(acrossRESS_SNR(:,[usefreqs],:));% TGs, ignoring BG.
    
   figure(1)   ; clf

    mbar= squeeze(nanmean(storeall,1));
    bh=bar(mbar') %use transpose to display with LOC on x axis
%     for ib=1:4
%     bh(ib).LineStyle=':';
%     bh(ib).LineWidth=4;
%     end
    ylabel('SNR(dB/Hz)')
%     xlabel('Hz')
%     legend TL TR BL BR
xlabel('Target Location')
% legend({'8 Hz', '13 Hz', '15 Hz', '18 Hz'})
legend({'8 Hz', '13 Hz', '15 Hz', '18 Hz', '16 Hz', '26 Hz', '30 Hz', '36 Hz'})
% legend({'16 Hz', '26 Hz', '30 Hz', '36 Hz'})
    shg
    %%
% end
%     storeall=
    X=storeall;
    pX=squeeze(mean(nanmean(storeall,2),3));
    meanOverall = mean(pX);
    
    tmpPX= ones(size(storeall,1),size(storeall,2), size(storeall,3));
    for i=1:size(pX)
        tmpPX(i,:,:) = tmpPX(i,:,:).*pX(i);
    end
    
    tmpGX= ones(size(storeall,1),size(storeall,2), size(storeall,3));
    adjustedST = X - tmpPX + (tmpGX.*meanOverall);
    %
    barstore=squeeze(nanstd(adjustedST))./sqrt(size(adjustedST,1));
    %
    barstore=barstore;
    hold on
    %
    numgroups = 4;
    numbars = size(storeall,2);%length(peakfreqsare);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    % determine locations for errorbars:
    
    for i = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
%         % use this row for Hz on x axis,%         
%                 errorbar(x, mbar(:,i), barstore(:,i), 'k', 'linestyle', 'none');
% use this row for LOC on x axis
        errorbar(x, mbar(i,:), barstore(i,:), 'k', 'linestyle', 'none');
    end
    %%
    %
%     title({['Pre and post ' dtypep ' onset']})
%     set(gca,'xticklabel',{'8 Hz', '13 Hz', '15 Hz' ,'18 Hz'}, 'fontsize', 20)
    set(gca,'xticklabel',{'Top Left', 'Top Right', 'Bottom Left' ,'Bottom Right'}, 'fontsize', 25)
    ylabel('SNR(dB/Hz)')

%      set(lga, 'location', 'northeast', 'fontsize', 15);
    set(gcf, 'color', 'w')
    ylim([ 0 6])
     shg
     cd(basefol)
     cd('newplots-MD')
     cd('figures')
     %%
     print('-dpng', 'ress stat-TG SNR by Hz x Loc')
     
     %% export for JASP
       %% export for stats analysis
    %
    TableforExport=table();
    thisentry=1; %initialize counter.
    %
    peakfreqsare=[8,13,15,18,16,26,30,36];
%     st
    for iloc = 1:4
        for ihz=1:length(peakfreqsare)
            
            %loc data =
            ppantdata = squeeze(storeall(:, ihz, iloc));
            
            %append horizontally to table.
            TableforExport=[TableforExport, table(ppantdata)];
            TableforExport.Properties.VariableNames(thisentry)={['LOC' num2str(iloc) '_' num2str(peakfreqsare(ihz)) 'Hz']};
            thisentry=thisentry+1;
        end
    end
    
    try cd('/Users/matthewdavidson/Desktop')
    catch
        cd('/Users/mattdavidson/Desktop')
    end
    
    writetable(TableforExport, 'RESS-TG hz x loc.csv')
end
