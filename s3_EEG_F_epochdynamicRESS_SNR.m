% Follows the format of s3_CE... only now applying to RESS timeseries

%having epoched tgs and catch, POST RESS. apply FFT and SNR to windows when
%all targets are present/ away from Catch stimulus onset.

%%%%% NEW FOR WILLS DATA!!!
clear all
% close all


% addpath('/Users/MatthewDavidson/Desktop/SSVEP-Wills')
% cd('/Users/MatthewDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')
% 
% addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/mDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')



basefol=pwd;
clearvars -except basefol allppants
dbstop if error
%%


job.calcppantDYNSNRperfreq=0; %sorts by hz x location. %topos in preivous script (s3_Da)

appendtradSNRtoDYN=0; % runs through the above, but performs traditional SNR on RESStimecourse for comparison.



job.concatacrossppanst=0;

job.plot_dynSSEP_acrossppants=0;


%separate analysis at image level.
% for PFI
job.erpimagePpantlevel=0; %using dynRESS SSVEP or SNR
job.concaterpdataacrossppants=0; % concat above
job.erpimageacrossppants=0;
job.gradedchangesinERPimage_PFI=0;


% also tr x tr of catch responses.
% for catch
job.erpimage_catchPpantlevel=0; %stores responses as 100 ydimension trial x trial
job.concaterpdataforcatchacrossppants=0;
job.erpimage_forcatchacrossppants=0;
job.gradedchangesinERPimage_Catch=0;




job.BPandSSVEPtimecourseacrossppants_group=0; %same format as previous script (epoch PFI etc.)
job.BPandSSVEPtimecourseacrossppants_group_CATCH=0; %same format as previous script (epoch PFI etc.)

%%Paper plots: (see also plotSmallFigs.m)
job.BPandSSVEPtimecourseacrossppants_group_combinePFIandCATCH=1; %same format as previous script (epoch PFI etc.)

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


%%%%%%% the below are not up to date! (old scripts)

job.calcPpantDYNSSVEP_crosspoint=0;
job.plotPpant_crosspoints=0;
job.FirstSIGintimecourse=0; %no longer in use.
job.concaterpdataacrossppants_keepnumberSeparate=0; %trying to replicate plots, now separates by either side of median PFI
job.BPandSSVEPtimecourseacrossppants_numSeparate=0; 



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


if job.calcppantDYNSNRperfreq==1
    %Collect all DATA
    %%
    %%
    
    
    for ifol =allppants(6:end)
        
        
        %%
        cd(basefol)
        cd('EEG')
        cd(pdirs(ifol).name)
       
        %% load the relevant PFI data.
       
        load('ppant_PFI_Epoched_RESS');
        load('ppant_Catch_Epoched_RESS');
        

        %%
        
        RESS_dynamicSNR_byTYPExHz=nan(13,length(peakfreqsare),1501);
        RESS_traditionalSNR_byTYPExHz=nan(13,length(peakfreqsare),Nf);

        
        for ifreq=1:length(peakfreqsare)
            
            usehz=peakfreqsare(ifreq);
            
           
                
                
                
                % collect relevant trials for each type of spatial
                % configuration/filter construction.
                
                %%
                
                for id=1:13% Use all epochs so as not to bias condition comparisons.
                    
                    switch id
                        case 1
                            dataIN=ress_PFI_0_1_Hz;
                            
                            dirsi= 1; %increasing PFI
                            %means we subtrat -3 -2 as bsrem.
                        case 2
                            dataIN=ress_PFI_1_0_Hz;
                           
                            dirsi= 0; %decre PFI
                            %means we subtrat 2 3 as bsrem.
                        case 3
                            dataIN=ress_PFI_1_2_Hz;
                           
                            
                            dirsi=1;
                        case 4
                            dataIN=ress_PFI_2_1_Hz;

                            
                            dirsi=0;
                        case 5
                            dataIN=ress_PFI_2_3_Hz;
                            
                      
                            
                            dirsi=1;
                            
                        case 6
                            dataIN=ress_PFI_3_2_Hz;
                          
                            dirsi=0;
                            
                        case 7
                              dataIN=ress_PFI_3_4_Hz;
                          
                            dirsi=1; %if increasing PFI
                        case 8
                            dataIN=ress_PFI_4_3_Hz;
                          
                            dirsi=0; %if increasing PFI
                            
                        case 9
                            if ifreq==1 || ifreq==3 || ifreq==5
                            dataIN=ress_BPcatchonsetTGs;
                            elseif  ifreq==2 || ifreq==4 || ifreq==6
                                dataIN=ress_BPcatchonsetBGs;
                            elseif   ifreq==7 || ifreq==8 || ifreq==9
                            dataIN=ress_BPcatchonsetIMs;
                            end    
                            usewindow=1;
                            dirsi=1;
                        
                        case 10
                            if ifreq==1 || ifreq==3 || ifreq==5
                                dataIN=ress_catchonsetTGs;
                            elseif  ifreq==2 || ifreq==4 || ifreq==6
                                dataIN=ress_catchonsetBGs;
                            elseif   ifreq==7 || ifreq==8 || ifreq==9
                                dataIN=ress_catchonsetIMs;
                            end
                            usewindow=1;
                            dirsi=1;
                        
                        
                        case 11
                            
                            if ifreq==1 || ifreq==3 || ifreq==5
                                dataIN=ress_BPcatchoffsetTGs;
                            elseif  ifreq==2 || ifreq==4 || ifreq==6
                                dataIN=ress_BPcatchoffsetBGs;
                            elseif   ifreq==7 || ifreq==8 || ifreq==9
                                dataIN=ress_BPcatchoffsetIMs;
                            end
                            
                            dirsi=0;
                            usewindow=2;
                        
                        case 12
                             if ifreq==1 || ifreq==3 || ifreq==5
                                 dataIN=ress_catchoffsetTGs;
                             elseif  ifreq==2 || ifreq==4 || ifreq==6
                                 dataIN=ress_catchoffsetBGs;
                            elseif   ifreq==7 || ifreq==8 || ifreq==9
                                dataIN=ress_catchoffsetIMs;
                            end    
                            
                            dirsi=0;
                            usewindow=2;
                            
                            
                            %%%%% NEW CASE, the invisible catch onsets!
                        case 13
                            %%%%% NEW CASE, the invisible catch onsets!
                             if ifreq==1 || ifreq==3 || ifreq==5
                                 dataIN=ress_invisiblecatchonsetTGs;
                             elseif  ifreq==2 || ifreq==4 || ifreq==6
                                 dataIN=ress_invisiblecatchonsetBGs;
                            elseif   ifreq==7 || ifreq==8 || ifreq==9
                                dataIN=ress_invisiblecatchonsetIMs;
                            end    
                            
                            dirsi=0;
                            usewindow=2;
                            
                            
                    end
                    
                    
                    % now that we have a datatype, get the right epochs that share stimulus configuration
                    
                    % also correct window (pre or post indication of target
                    % presence)
                    %reduce size.
                    if id<9
                        datast=squeeze(dataIN(ifreq).ressTS);
                        
                        
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
                        datast= squeeze(dataIN(usehztmp,:,:));
                    end
                        
                    
                    if size(datast,1)>1
                        %% check for bad trials (noisy)
                        %std per trial(average over timepoints)
                        %
                        datastSD = nanstd(datast,0,2);
                        %
                        %                     %remove those with 2.5*std from mean.
                        trialSD=nanstd(datastSD);
                        mSD=nanmean(datastSD);
                        keeptrials=1:size(datast,1);
                        %%
                        %remove the trials with excessive noise.
                        %                     badtrials=[];
                        badtrials=find(datastSD>mSD+2.5*trialSD)';
                        
                        % also skip the trials which were only transient button
                        % presses. (less than one second).
                        shorttrials=[];
                        
                        %also remove NANs:
                        nantrials = find(isnan(datast(:,1)));
                        zerotrials= [];%find(datast(:,1)==0);
                        badtrials = [badtrials, shorttrials, nantrials', zerotrials'];
                        
                        
                        % remove these from consideration.
                        keeptrials(badtrials)=[];
                        datast=datast(keeptrials,:);
                        %%
                                                                     
                        %now we have ALL the data this freqxLoc.
                        %%
                        if appendtradSNRtoDYN~=1 % then continue with
                            % dyn SSVEP (as per Cohen &
                            % Gulbinaite, 2017). or straight SNR.
                            
                            
                            %dynamic SSVEP given by abs(hilbert transform!
                            %first filter:
                            peakwidt=2;
                            %filter, but inlcude side lobes for temporal dynamics
                            fdatAt = filterFGx(datast,srate,usehz,peakwidt,0);
                            
                            
                            % take absolute of hilbert transform, to show dynamic SSEP
                            try dynSSEP= abs(hilbert(fdatAt'))';
                                
                                
                                %                 %TAKE mean per trial, careful with baseline subtraction
                                %                 timeid = [0:1/srate:epochdur];
                                %                 timeid= timeid-3;
                                %                 if dirsi==1 %increasing PFI / Target absence
                                %                     if ifreq<5 %targets, so low point post event.
                                %                         windcheck= [2 3];
                                %                     else %BG hz, so
                                %                         windcheck= [-3 -1];
                                %                     end
                                %                 elseif dirsi==0 %oppoite case, decreasing PFI/ targets are appearing
                                %                     if ifreq<5 %targets, so low point pre event.
                                %                         windcheck= [-3 -2];
                                %                     else %BG hz, so
                                %                         windcheck= [2 3];
                                %                     end
                                %                 end
                                %remove one second baseline
                                windcheck=[-2.9 -1.5];
                                %                 windcheck=[-3 3]; %whole epoch.
                                %
                                tidx=dsearchn(timeidDYN', [windcheck]');
                                %
                                %                %remove base
                                for itrial=1:size(dynSSEP,1)
                                    tmp= dynSSEP(itrial,:);
                                    mtmp = nanmean(tmp(tidx(1):tidx(2)));
                                    %                    mtmp = mean(tmp);
                                    
                                    dynSSEP(itrial,:)= tmp - repmat(mtmp, 1, length(tmp));
                                    
                                end
                                
                                
                                
                            catch
                                dynSSEP=nan(1,1501);
                            end
                            
                            
                            % store this type
                            RESS_dynamicSNR_byTYPExHz(id,ifreq,:) = squeeze(nanmean(dynSSEP,1));
                            
                        else % spectrogram of this RESS ts for comparison
                            %copied params from previous (EPOCH scripts).
                            
                            
                            
                         
                            
%                             if ifreq==7
%                                 pause
%                             end
                            [sgrm ,tgrm, fgrm] = mtspecgramc(datast', movingwin, param_spcgrm);
                            %%
                            %conv SNR
                            snr_sgrm =zeros(size(sgrm));
                            
                            %adjust for HBW
                            k = ((param_spcgrm.tapers(1,1)));%
                            hbw = (k+1)./ (2.*[movingwin(1,1)]);
                            neighb = 1; %hz to compare noisewidth with.
                            
                            %space out kernel appropriately
                            distz = dsearchn(fgrm', [hbw, neighb, neighb*2+hbw]');
                   
tmps= squeeze(nanmean(sgrm,3));
                                
                                
%                                 % SNR by convolution.
%                                 snr_sgrm=[];
%                                 for itime= 1:size(tmps,1)
%                                     checkput = conv(log(tmps(itime,:)), kernelw,'same');
%                                     if ~isreal(checkput)
%                                         snr_sgrm(itime,:)= nan(1, size(tmps,2));
%                                     else
%                                         snr_sgrm(itime,:)= conv(log(tmps(itime,:)), kernelw,'same');
%                                     end
%                                 end



% loop over frequencies and compute SNR
numbins = distz(2); 
skipbins = distz(1);
snr_sgrm= zeros(size(tmps));
for itime=1:size(tmps,1)
    for hzi=numbins+1:length(fgrm)-numbins-1
        numer = tmps(itime,hzi);
        denom = nanmean( tmps(itime,[hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
        snr_sgrm(itime,hzi) = numer./denom;
        
    end
end
                                
                                %%
                                
                                %                             %store just stim freq.
                            [~, hzid]= min(abs(fgrm-usehz));
                            %
                            %reduce size.
                            mdynSSEP=squeeze(snr_sgrm(:,hzid,:));

                            RESS_traditionalSNR_byTYPExHz(id,ifreq,:) =mdynSSEP';
%                            

timeidGRM=tgrm-3;
%                            

                                %% sanitycheck
%                                 subplot(221)
%                                 imagesc(tgrm,fgrm, log(tmps)')
%                                 ylim([0 50])
%                                 subplot(222)
%                                 imagesc(tgrm, fgrm, snr_sgrm')
%                                 ylim([0 50])
%                                 c=colorbar;
%                                 caxis([0 2]);
%                                 subplot(2,2,3:4); plot(timeidGRM, log(tmps(:,hzid))')
%                                 hold on;  plot(timeidGRM, log(snr_sgrm(:,hzid))')
%                                 
                        end
                        
                        %% store mean per type per ppant.
                    else %else store nan if no PFI this type.
                        dynSSEP=datast;
                    end
        end
        end
        
        DIMSare = {'0_1', '1_0', '1_2', '2_1', '2_3', '3_2','3_4', '4_3',...
            'BPonsets', 'catch onset', 'BPoffsets', 'catchoffset', 'invisible catchonset'};
         
        
        
        %save whether traditional or dynamic amplitude, first or second
        %harmonic.
        
        
        if appendtradSNRtoDYN~=1
            
            
            save('RESS_dynamicSNR_byTYPExHz', 'RESS_dynamicSNR_byTYPExHz', 'timeidDYN', 'DIMSare')
        else
            try save('RESS_dynamicSNR_byTYPExHz', 'RESS_traditionalSNR_byTYPExHz', 'timeidGRM', 'DIMSare','-append')
            catch
                save('RESS_dynamicSNR_byTYPExHz', 'RESS_traditionalSNR_byTYPExHz', 'timeidGRM', 'DIMSare')
            end
                
        end
        
        disp(['fin ppant ' num2str(ifol)]);
        
    end
end




if job.concatacrossppanst==1
    
    %%
    
    acrossRESS_DYNSNR_freqs=zeros(length(allppants),13, length(peakfreqsare),1501);
    acrossRESS_TRADSNR_freqs=zeros(length(allppants),13,  length(peakfreqsare),Nf);
%     [nppants, ndims, nfreqs,nlocs,nsamps]=size(acrossRESS_DYNSNR_freqs);
    icounter=1;
    for ippant=allppants
        cd(basefol)
        cd('EEG')
        cd(pdirs(ippant).name)
        
        load('RESS_dynamicSNR_byTYPExHz')
        %save both types.
%         acrossRESS_DYNSNR_freqs(icounter,:,:,:)=RESS_dynamicSNR_byTYPExHz;        
        acrossRESS_TRADSNR_freqs(icounter,:,:,:)=RESS_traditionalSNR_byTYPExHz;
        
        icounter=icounter+1;
    end
    %
    cd(basefol)
    cd('EEG')
    cd('GFX_EEG-RESS')
    save('GFX_ressSNR_Dynamic', 'acrossRESS_DYNSNR_freqs','acrossRESS_TRADSNR_freqs', 'timeidDYN', 'timeidGRM', 'DIMSare')
    
end

%%

if job.plot_dynSSEP_acrossppants==1
    %%
    cd(basefol)
    cd('EEG')    
    
    %%
    cd('GFX_Pre-RESS')
    load('GFX_Catchperformance_withSNR_15,BPaligned.mat', 'storeacrossPpant_onsetBP', 'storeacrossPpant_offsetBP');
    BPoncatch=storeacrossPpant_onsetBP;
    BPoffcatch=storeacrossPpant_offsetBP;
    load('GFX_PFIperformance_withSNR_15_min0.mat', 'storeacrossPpant_onsetBP', 'storeacrossPpant_offsetBP');
    BPonPFI=storeacrossPpant_onsetBP;
    BPoffPFI=storeacrossPpant_offsetBP;
   %%   
    cd(basefol)
    cd('EEG')    
    cd('GFX_EEG-RESS')
    
    % ThiS SNR is average across subjects    
    load('GFX_ressSNR_Dynamic.mat');
    
    useTRADorDYNSNR = 2; % plot traditional or dynamic SNR after RESS.
    
    
    
    switch useTRADorDYNSNR
        case 1
            usedata=acrossRESS_TRADSNR_freqs;
            timeid=timeidGRM;
        case 2
            usedata=acrossRESS_DYNSNR_freqs;
            timeid=timeidDYN;
    end
    %
    figure(1)
%     clf
   
  

    PFIincrease = [1,3,5,7]; % %%%% change me for WILLS
    PFIdecrease = [2,4,6,8];
    CatchONSET= 10;
    CatchOFFSET= 12;
    BPCatchONSET= 9;
    BPCatchOFFSET= 11;
    invisibleCATCH=13;
    
    checkcluster=0;
    counter=1;
  
pl=[];
legendPRINT={};
       ttestdata=[];
       
       
       printOVERLAY=1; % set if want to print ontop (same columns).
       clf
         
       %move ttest data if comparing types.
       ttestdata=[];
       
       %blue color for target, red for BG.
       hzcols={'b', 'r', 'b', 'r', 'b', 'm', 'k','k','k'};
       
       for ihz=[1]%
% {1,3,5}
%         typesnr='Target SNR'; 1f, 2f, 3f
% {2,4,6}
%         typesnr= 'Background SNR'; 1f, 2f, 3f
% {7,8,9}
%         typesnr= 'Intermodulation SNR'; f2-f1, 2f2+f1, f1+f2

          counter=1;
          
          itypech=[1,2]%[1,2]; %used for indexing ttests
          
          for itype=itypech
                
    
            switch itype
                case 1
                    EVENTdata = squeeze(nanmean(usedata(:,PFIincrease,:,:,:),2));
%                     col = [.5 .7 .5];% dark green
                    
                    col='b';
locis='NorthEast';
                    chtype = {['PFI Increase'];['Target Disappearance']};
                    xlabelis = 'Time from target invisible';
                    
                    mrks='-';
                    
                    %load BP data.
                    lge='target disappear';
                    useBP=BPonPFI;
                case 2
                    EVENTdata = squeeze(nanmean(usedata(:,PFIdecrease,:,:,:),2));
                    col='r';
                    chtype = {['PFI Decrease'];['Target Reappearance']};
                    locis='SouthEast';
                    xlabelis = 'Time from target visible';                    
                    useBP=BPoffPFI;
                    mrks=':';
                    lge='target reappear';
                
                case 3
                    EVENTdata = squeeze(nanmean(usedata(:,CatchONSET,:,:,:),2));
                    col='g';
                    chtype = 'catch onset';
                    xlabelis = 'Time from catch onset';
                    
                case 4
                    EVENTdata = squeeze(nanmean(usedata(:,CatchOFFSET,:,:,:),2));
                    col='r';
                    chtype = 'catch offset';
                    xlabelis = 'Time from catch offset';
%                     
                case 5
                    EVENTdata = squeeze(nanmean(usedata(:,BPCatchONSET,:,:,:),2));
                    col = [.5 .7 .5];
                    chtype = {['Catch'];['Target Disappearance']};
                    xlabelis = 'Time from reporting catch onset';
                    mrks='-';
                    useBP=BPoncatch;
                                        lge='catch onset';

                case 6
                    EVENTdata = squeeze(nanmean(usedata(:,BPCatchOFFSET,:,:,:),2));
                    col='r';
                    mrks=':';
                    chtype = {['Catch'];['Target Reappearance']};
                    xlabelis = 'Time from reporting catch offset';
                    lge='catch offset';
                    
                    useBP=BPoffcatch;
                    
                    
                    
                     case 7
                    EVENTdata = squeeze(nanmean(usedata(:,invisibleCATCH,:,:,:),2));
                    col='k';
                    mrks=':';
                    chtype = {['invisible catch']};
                    xlabelis = 'Time from reporting catch offset';
                    lge='catch offset';
                    
%                     useBP=BPoffcatch;
            end
            
            
               
% restrict to HZ of interest
                                        
                    plotd= squeeze(EVENTdata(:,ihz,:));
%                     plotd=plotd-mean(plotd(:));
                    
switch ihz
    case {1,3,5}
        typesnr='Target SNR';
    case {2,4,6}
        typesnr= 'Background SNR';
    case {7,8,9}
        typesnr= 'Intermodulation SNR';
end


%to place on the same figure!
if printOVERLAY==1 && mod(counter,2)==0
    counter=counter-1;
end

             subplot(2,2,counter)
            %plot BP.
            hold on
            tvector=[-3:1/60:3];
            
             ppantMeanBP = squeeze(mean(useBP,2));
             stE = std(ppantMeanBP)/sqrt(size(ppantMeanBP,1));
             sh=shadedErrorBar(tvector, mean(ppantMeanBP,1),stE,['-'],[1]);
              
             if itype==1 && printOVERLAY==1
             sh.mainLine.LineWidth=6;
             else
                 sh.mainLine.LineWidth=3;
             end
                    
              if printOVERLAY==11 && (itype==5 || itype==6)
                  colBP=[.3 .3 .3];
              else
                  colBP='b';
              end
                    sh.mainLine.Color = colBP;
                    sh.patch.FaceColor = colBP;
                    sh.edge(1).Color = colBP;
                    sh.edge(2).Color = colBP;
                    sh.mainLine.LineStyle=mrks;
             ylabel(['Buttons pressed'])
             set(gca, 'fontsize', 20)    
             
             xlim([-3 3])
             ylim([0 4])
             if printOVERLAY~=1
             xlabel(xlabelis)
             else
             xlabel('Time from subjective report')
             end
%              
             subplot(2,2,counter+2)
             
             if ihz~=2 && ihz~=3
                 ylimsare=[0 .45];
             else
                 ylimsare=[1 2.2];
             end


col= hzcols{ihz};

            %SUBRTACT overall MEAN
%                       plotd=plotd-nanmean(plotd(:));
            
            %now plot %mean over ppants.
            mplot = squeeze(nanmean(plotd,1)); 
            
            
            %adjust errorbars.
            X=plotd;
            pX = nanmean(plotd,2);
            mX= nanmean(pX);
            
            newX = X - repmat(pX, [1, size(X,2)]) + repmat(mX, [size(X,1), size(X,2)]);
            
            stE = nanstd(newX,0,1)./sqrt(size(newX,1));
            
            sh=shadedErrorBar(timeid, mplot, stE, [ mrks], 1);
            
           
           
            if printOVERLAY==11 &&  (itype==5  ||itype==6 )     
           
            sh.mainLine.Color='k';
            sh.patch.FaceColor=col;
            sh.edge(1).Color='k';
            sh.edge(2).Color='k';
            
           else
               sh.mainLine.Color=col;
            sh.patch.FaceColor=col;
            sh.edge(1).Color=col;
            sh.edge(2).Color=col;
            sh.mainLine.LineWidth=3;
           end
            
            pl(itype)=sh.mainLine;
            hold on
%        
            set(gca, 'fontsize', 20)            
             
            xlim([-3 3])
            
            if printOVERLAY~=1
                xlabel(xlabelis)
            else
                xlabel('Time from subjective report')
            end
            
    
        ylabel({['RESS log(SNR)']})
    ttestdata(itype,:,:)=plotd;
        
%         if ihz==1
        
%         end
        
        counter=counter+1;
        
    
%         ylim([-.05 .5])
        xlim([-2.5 2.5])
        
        
        legendPRINT = [legendPRINT {lge}];
          end            
       
          
          
          
    %
    hold on; plot([0 0 ], ylim, ['k-'])
set(gcf, 'color', 'w')
% ylim([1.1 2.25])
    xlim([timeid(1) timeid(end)])
          
          
          
          
          
          
          
          
          try lg=legend([pl(itypech(1)) pl(itypech(2)) ], legendPRINT);    
         set(lg, 'Location','SouthEast', 'fontsize',18)
          catch
          end
        title([typesnr ])
%      check for significance from zero
%check for sig
%         pvals=zeros(1,size(ttestdata,3));
if checkcluster==1
        tvals=zeros(1,size(ttestdata,3));
        for itime = 1:size(ttestdata,3)
            
            try [h,pvals(itime),~,stat]=ttest(ttestdata(itypech(1),:,itime), ttestdata(itypech(2),:,itime));
                shuffType=1;
            catch
                [h,pvals(itime),~,stat]=ttest(ttestdata(itype,:,itime)); %compares to zero.
                shuffType=2; %whether or not to skip the non-parametric test for sig.
            end
            
            tvals(itime)= stat.tstat;
        end
        sigs=find(pvals<.05);

%         sigs=[];
%             %perform cluster based correction.
            if length(sigs)>0
                % find biggest cluster:
                %finds adjacent time points
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                d = diff(v1);
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                %grab largest
%                 ignore bad points.
                
                
                  % find biggest cluster:
                %finds adjacent time points
                sigs = find(pvals<.05);
                
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                d = diff(v1);
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
               
                %%
%                 sigheight=.45
                for icl=1:size(clusterSTandEND,1)

                    %start and end are now:
                    % change icl to maxClust if we only want the largest
                    % cluster.
                    STC=sigs(clusterSTandEND(icl,1));
                    ENDC=sigs(clusterSTandEND(icl,2)+1);
                    checktimes =STC:ENDC;
                    observedCV = sum(abs(tvals(checktimes)));
                    % now shuffle condition labels to see if this cluster is
                    % sig (compared to chance).
                    nshuffs=500;
                    sumTestStatsShuff = zeros(1,nshuffs);
                    
                    for irand = 1:nshuffs
                        %testing the null that it isn't mismatched - matched at time 2
                        % which creates a diff. so select from either!
                        shD= zeros(2,length(checktimes),size(ttestdata,2));
                        
                        %change shuffle parameters based on test of
                        %interest. (ie between conditions, or temporal
                        %null).
                        if shuffType==1 %null is that no condition differences.
                        for ipartition = 1:2
                            for ippant = 1:size(ttestdata,2)
                                for itime=1:length(checktimes)
                                    
                                    if mod(randi(100),2)==0 %if random even number
                                        pdata = ttestdata(itypech(1),randi(size(ttestdata,2)), checktimes(itime)); %select both chans
                                    else
                                        pdata = ttestdata(itypech(2),randi(size(ttestdata,2)), checktimes(itime));
                                    end
                                    
                                    shD(ipartition,itime,ippant) = pdata;
                                end
                            end
                        end
                        else %null is that there are no temporal coincident sig values.
                            for ipartition = 1:2
                            for ippant = 1:size(ttestdata,2)
                                for itime=1:length(checktimes)
                                    
                                    %take random timepoint.
                                        pdata = ttestdata(1,ippant, randi(size(ttestdata,3)));
                                    
                                    
                                    shD(ipartition,itime,ippant) = pdata;
                                end
                            end
                            end
                        end
                        %now compute difference between out hypothetical topoplots,
                        % and test for sig, checking the accumulated test statistic at our
                        % times of interest
                        tvalspertimepoint = zeros(1,length(checktimes));
                        
                        testdata = squeeze(shD(1,:,:)) - squeeze(shD(2,:,:));
                        
                        for itest = 1:length(checktimes) %test each time point
                            
                            [~, p, ~,stat]= ttest(testdata(itest,:));
                            
                            tvalspertimepoint(1,itest) = stat.tstat;
                        end
                        
                        sumTestStatsShuff(1,irand) = sum(abs(tvalspertimepoint));
                    end %repeat nshuff times
                    
                    
                    %is the observed greater than CV?
                    % plot histogram:
                    figure(2);
                    
                        clf
                    
%                     subplot(2,1, plcount-1)
                    H=histogram(abs(sort(sumTestStatsShuff)));
                    % fit CDF
                    cdf= cumsum(H.Data)/ sum(H.Data);
                    %the X values (actual CV) corresponding to .01
                    [~,cv05uncorr] = (min(abs(cdf-.95)));
                    [~,cv01uncorr] = (min(abs(cdf-.99)));
                    [~,cv001uncorr] = (min(abs(cdf-.999)));
                    hold on
                    pCV=plot([observedCV observedCV], ylim, ['r-']);
                    
                    p05=plot([H.Data(cv05uncorr) H.Data(cv05uncorr)], ylim, ['k:']);
                    plot([H.Data(cv01uncorr) H.Data(cv01uncorr)], ylim, ['k:']);
                    plot([H.Data(cv001uncorr) H.Data(cv001uncorr)], ylim, ['k:']);
                    legend([pCV p05], {['observed'] ['p01'] })
                     
                    %%
                    if observedCV>H.Data(cv05uncorr)
                          [~, c2] = min(abs(H.Data-observedCV)); %closest value. 
                        pvalis= 1-cdf(c2);
                    title(['sum tvals = ' num2str(observedCV), 'p=' num2str(pvalis)]);
%                     title(['sum tvals = ' num2str(observedCV)]);
                    %              title('Spatial Cluster  Significant!')
                    figure(1); hold on;
                    %space out the clusterfor plotting
                    if length(checktimes)>30
                    tryt= downsample(checktimes,30);
                    else
                        tryt=checktimes;
                    end
                    timeid(checktimes)
                    figure(1); hold on;
                        for itime=tryt
                            figure(useTRADorDYNSNR);
                            hold on
                            ystretch=get(gca, 'ylim');
                            sigheight= (ystretch(1) + ((ystretch(2)-ystretch(1))*.1)); 
%                             plot(timeid(itime), sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', sh.mainLine.Color)
                            plot(timeid(itime), [sigheight], ['*' ],'markersize', 15, 'linewidth', 3, 'color', col)
                        end
                    end
            
%            
            end
            elseif length(sigs)>2
                
%                     tryt= downsample(sigs,30);
%                     
%                     
%                         for itime=tryt
%                             figure(1);
%                             hold on
%                             plot(timeidDYN(itime), sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'k')
% %                             plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'm')
%                         end
            end
end
%    
    end
        
    %%
    hold on; plot([0 0 ], ylim, ['k-'])
set(gcf, 'color', 'w')
% ylim([1.1 2.25])
    xlim([timeid(1) timeid(end)])
cd(basefol)
cd('Figures')

% cd('"X" over time')

% print('-dpng', ['Dynamic SSVEP during report PFI TG 2nd harm'])
%     print('-dpng', ['Dynamic SSVEP during CATCH_BG.png'])

%     print('-dpng', ['Dynamic SSVEP during PFI_BG.png'])
end



if job.erpimagePpantlevel==1
    %%
    rmvbase=0;
    snrmethod=2; % 1 for division method (MXC), 2 for dot product convolution.
    useSNRorhilbert=1;
    for ifol = allppants
        
        for ihz=[3,4,5,6]; %TGs1 TGs 2, 20hz 40 hz.
            
            
            usehz= peakfreqsare(ihz);
            
            
            
            icounter=1;
            
             cd(basefol)
        cd('EEG')
        cd(pdirs(ifol).name)
            load(['ppant_PFI_Epoched_RESS'])
            
            for itimezero = 1:2
                if itimezero==1
                    
                    %append them all. TARG-> Disappearing (more buttons
                    %pressed).

                    datatouse = cat(1, ress_PFI_0_1_Hz(ihz).ressTS,ress_PFI_1_2_Hz(ihz).ressTS,ress_PFI_2_3_Hz(ihz).ressTS, ress_PFI_3_4_Hz(ihz).ressTS);
%                     RTstouse = [durs0_1'; durs1_2'; durs2_3'];
                    BPstouse = cat(1, ress_PFI_0_1_Hz(ihz).BPs,ress_PFI_1_2_Hz(ihz).BPs,ress_PFI_2_3_Hz(ihz).BPs,ress_PFI_3_4_Hz(ihz).BPs);
                    ctype = 'PFI increase';
                    bsrem = [-3 -1]; %seconds
                    
                else
                    
                    datatouse = cat(1, ress_PFI_1_0_Hz(ihz).ressTS,ress_PFI_2_1_Hz(ihz).ressTS,ress_PFI_3_2_Hz(ihz).ressTS,ress_PFI_4_3_Hz(ihz).ressTS);
%                     RTstouse = [durs1_0'; durs2_1'; durs3_2'];
                    BPstouse= cat(1, ress_PFI_1_0_Hz(ihz).BPs,ress_PFI_2_1_Hz(ihz).BPs,ress_PFI_3_2_Hz(ihz).BPs,ress_PFI_4_3_Hz(ihz).BPs);
                    ctype = 'PFI decrease';
                    
                    
                    bsrem = [1 3]; %seconds
                    
                    
                end
                
                
                
                %plot the topo for sanity check: (pre disap)
                windcheck= [-3 -0.1];
                tidx=dsearchn(timeidDYN', [windcheck]');
                
                
                
                %restrict to certain PFI duration?
                allt= 1:size(datatouse,1);
                %                 baddurs = find(RTstouse<30);
                %remove those trials.
                %                 allt(baddurs)=[];
                
                
                
               
                
                %%
                %%
                clf
                colormap('viridis')
                %plot BP for comparison
                subplot(3,2, [1 3])
                %sort trial index by longest RT first.
                %                 [sortedRTs, cid] = sort(RTstouse(allt), 'descend');
                
                %sort by accum PFI in window after onset.
                if itimezero==1
                    accumPeriod = sum(BPstouse(:,180:end),2);
                else
                    accumPeriod = sum(BPstouse(:,1:180),2);
                    
                end
                
                [checkBP, cid] = sort(accumPeriod, 'descend');
                
                
                % MOVE nan to bottom of order (not top)
                nanind= find(isnan(checkBP));
                %new start point
                cid=[cid(length(nanind)+1:end) ; cid(nanind)];
                %                     checkBP=[checkBP(length(nanind)+1:end); checkBP(nanind)]
                
                %rearrange.
                BPstouse=BPstouse(cid,:);
                datais=squeeze(datatouse(cid,:));
%                 RTstouse=RTstouse(cid);
                %%

                tBP=-3:1/60:3;
                imagesc(-3:1/60:3,1:length(cid),BPstouse)
                c=colorbar;
                ylabel(c, 'buttons pressed')
                %                 set(gca, 'ytick', 1:length(cid), 'yticklabel', round(sortedRTs./60,2))
                ylabel('PFI events')
                xlabel('Time [sec]')
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                title({['sorted BP data for ' ctype ','];['ppant' num2str(ifol)]})
                set(gca, 'fontsize', 15)
                %%
                % show mean over time.
                subplot(3,2,5)
                plot(tBP, nanmean(BPstouse,1),['k-'], 'linewidth', 3)
                %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                set(gca, 'fontsize', 15)
                hold on;
                ylabel('nanmean BP ')
                %
                xlabel('Time secs')
                axis tight
                plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                
                
                set(gca, 'fontsize', 15)
                %%
                
                if useSNRorhilbert~=1
                %% now plot dynamic RESS:
                
                peakwidt=2;
                 %filter, but inlcude side lobes for temporal dynamics
                fdatAt = filterFGx(datais,srate,usehz,peakwidt,0);                
                
                % take absolute of hilbert transform, to show dynamic SSEP
                dynSSEP= abs(hilbert(fdatAt'))';
                
                
%                 %TAKE mean per trial, careful with baseline subtraction
%                 timeid = [0:1/srate:epochdur];
%                 timeid= timeid-3;
%                 if dirsi==1 %increasing PFI / Target absence
%                     if ifreq<5 %targets, so low point post event.
%                         windcheck= [2 3];
%                     else %BG hz, so
%                         windcheck= [-3 -1];
%                     end
%                 elseif dirsi==0 %oppoite case, decreasing PFI/ targets are appearing
%                     if ifreq<5 %targets, so low point pre event.
%                         windcheck= [-3 -2];
%                     else %BG hz, so
%                         windcheck= [2 3];
%                     end
%                 end

%                 
                else % use SNR
                     
                            
                           
%                             if ifreq==7
%                                 pause
%                             end
                            [sgrm ,tgrm, fgrm] = mtspecgramc(datais', movingwin, param_spcgrm);
                            %% %% Two SNR methods:
                            
                            %%
                            snr_sgrm =zeros(size(sgrm));
                            
                            tmps=sgrm;
                            %adjust for HBW 
                                k = ((param_spcgrm.tapers(1,1)));%
                                hbw = (k+1)./ (2.*[movingwin(1,1)]);
                                neighb = 2; %hz
                                
                                %space out kernel appropriately
                                distz = dsearchn(fgrm', [hbw, neighb, neighb*2+hbw]');
                                
                                % we don't want odd numbers, messes with
                                % calculations, so round allup to even.
                                ods = find(mod(distz,2)~=0);
                                %adjust.
                                distz(ods)= distz(ods)+1;
                                
                                %make kernel
                                kernelneighb = -1*(repmat(1/(distz(2)*2), [1, distz(2)]));
                                kernelskip = zeros(1,distz(1));
                                
                                kernelw= [kernelneighb, kernelskip, 1, kernelskip, kernelneighb];
                                if sum(kernelw)~=0
                                    error('')
                                end
                            if snrmethod==1%method 1, using division.
                                
                               
                                
                                
                                
                                % loop over frequencies and compute SNR
                                numbins = distz(2);
                                skipbins = distz(1);
                                snr_sgrm= zeros(size(sgrm));
                                for itrial = 1:size(sgrm,3)
                                for itime=1:size(tmps,1)
                                    for hzi=numbins+1:length(fgrm)-numbins-1
                                        numer = tmps(itime,hzi, itrial);
                                        denom = nanmean( tmps(itime,[hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
                                        snr_sgrm(itime,hzi,itrial) = numer./denom;
                                        
                                    end
                                end
                                end
                                
                            else
                                %                             %conv SNR
                                %                             snr_sgrm =zeros(size(sgrm));
                                %
                                
                                %%
%                                 kernelw= [-1/8 -1/8 -1/8 -1/8 0 0 1 0 0 -1/8 -1/8 -1/8 -1/8];
                                
                                %snr on trials or
                                for itrial=1:size(sgrm,3)
                                    %compute SNR
                                    tmps= squeeze(sgrm(:,:,itrial));
                                    
                                    for itime= 1:size(tmps,1)
                                        checkput = conv(log(tmps(itime,:)), kernelw,'same');
                                        if ~isreal(checkput)
                                            snr_sgrm(itime,:,itrial)= nan(1, size(tmps,2));
                                        else
                                            snr_sgrm(itime,:,itrial)= conv(log(tmps(itime,:)), kernelw,'same');
                                        end
                                    end
                                end
                            
                             %% sanity check:
%                               figure(1);
%                               %plot specgrm and spctrm.
%                               subplot(221)
%                             imagesc(tgrm, fgrm, squeeze(mean(log(sgrm),3))'); colorbar; caxis([0 2])
%                             subplot(223)
%                             mtrials=squeeze(mean(log(sgrm),3));
%                             plot(fgrm, squeeze(mean(mtrials,1)));
%                             %plot specgrm and spctrm.
%                               subplot(222)
%                             imagesc(tgrm, fgrm, squeeze(mean((snr_sgrm),3))'); colorbar; caxis([0 2])
%                             subplot(224)
%                             mtrials=squeeze(mean((snr_sgrm),3));
%                             plot(fgrm, squeeze(mean(mtrials,1)));
%                             hold on; plot(kernelw)
                            end
                            %%
                            
                            %store just stim freq.
                            [~, hzid]= min(abs(fgrm-usehz));
                            %
% %                             %reduce size.
                            dynSSEP=squeeze(snr_sgrm(:,hzid,:))';
% %                             
                          timeidDYN=tgrm-3;

                            
                    
                    
                end
%                         
%                %remove base
if rmvbase==1
    tidx=dsearchn(timeidDYN', [bsrem]');
    for itrial=1:size(dynSSEP,1)
        tmp= dynSSEP(itrial,:);
        mtmp = mean(tmp(tidx(1):tidx(2)));
        %                    mtmp = mean(tmp);
        
        dynSSEP(itrial,:)= tmp - repmat(mtmp, 1, length(tmp));
        
    end
    
end
                
               %
               
                %%\
                snrgrm20=dynSSEP;
                % smooth across trials
                sm_snrgrm20=zeros(size(snrgrm20));
                for itime=1:size(snrgrm20,2)
                    sm_snrgrm20(:,itime)= smooth(snrgrm20(:,itime),5);
                end
                %%
                hold on
                subplot(3,2,[2 4]);
                hold on
                imagesc(timeidDYN,  1:size(dynSSEP,1), sm_snrgrm20);
                c=colorbar;
                if useSNRorhilbert==1
                ylabel(c, 'log(SNR)')
                else
                ylabel(c, 'abs(hilbert)')
                end
                
                xlabel('Time secs')
                caxis([-1*max(max(sm_snrgrm20)) max(max(sm_snrgrm20))])
                caxis([0  max(max(sm_snrgrm20))-1])%                 set(gca, 'ytick', 1:24, 'yticklabel', round(sortedRTs./60,2))
% caxis([ 1 5])
                axis tight
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                title({['sorted ' num2str(usehz) 'Hz RESS data (smoothed) for  ' ctype ','];['ppant' num2str(ifol)]})
                set(gca, 'fontsize', 15, 'yticklabel', [])
                
                %         ylim([0 20])
                %% show mean over time.
                subplot(3,2,6)
                plot(timeidDYN, nanmean(snrgrm20,1),['k-'], 'linewidth', 3)
                %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                set(gca, 'fontsize', 15)
%                 ylim([1 3])
                hold on;
                ylabel('RESS logSNR')
                %%
                axis tight
                plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                shg
                xlabel('Time secs')
                set(gcf,'color', 'w');
                
%                         ylim([-.1 10])
                
                %%
                shg
                cd([basefol filesep 'Figures' filesep 'Participant summaries'])
                %%
                %%
                print('-dpng', ['figure_PFI_' ctype '_summary_BPandSSVEP_' num2str(usehz) '_RESS_subj' num2str(ifol) ' snr method ' num2str(snrmethod) '_longerwindow.png']);
                
                
                switch itimezero
                    case 1 %store for across ppant plots:
                        Ppant_onsetBP=BPstouse;
                        Ppant_onsetSNR=snrgrm20; %sorted.
%                         Ppant_onsetRTs=RTstouse;
%                         Ppant_onsetTOPO=snr20;
                    case 2
                        
                        Ppant_offsetBP=BPstouse;
                        Ppant_offsetSNR=snrgrm20; %always sorted in descending order of PFI.
%                         Ppant_offsetRTs=RTstouse;
%                         Ppant_offsetTOPO=snr20;
                end
            end
            
            
            %%
            
                    savename=['PFIperformance_withSNR_' num2str(usehz) '_RESS'];
            
            cd(basefol)
            cd('EEG')
            cd(pdirs(ifol).name)
            %%
            save(savename,...
                'Ppant_onsetBP','Ppant_offsetBP',...
                'Ppant_onsetSNR', 'Ppant_offsetSNR', 'timeidDYN', 'param_spcgrm', 'kernelw')
            
            
        end
    end
end


%%

if job.concaterpdataacrossppants==1
    
    
    dursMINIMUM= 0; %1 second. % or change to  %this works.
%     dursMINIMUM=123; % selects top third!.
    
    ppantsmoothing=1; % average across participants, after smoothing, or no.
    
    for ihz=[3:6];%[1,2,7]
        usehz=peakfreqsare(ihz);
        
            loadname=['PFIperformance_withSNR_' num2str(usehz) '_RESS'];
        
        storeacrossPpant_onsetBP=[];
        storeacrossPpant_onsetSNR=[];
        storeacrossPpant_onsetRTs=[];
        storeacrossPpant_onsetTOPO=[];
        storeacrossPpant_offsetBP=[];
        storeacrossPpant_offsetSNR=[];
        storeacrossPpant_offsetRTs=[];
        storeacrossPpant_offsetTOPO=[];
        
        icounter=1;
        
        for ippant = allppants
             cd(basefol)
             cd('EEG')
             cd([pdirs(ippant).name])
                
                %onset types
                
                load(loadname)
%                 Ppant_onsetSNR=Ppant_onsetSNR';
%                 Ppant_offsetSNR=Ppant_offsetSNR';
                
                %restrict to only 'long' PFIs if needed.
%                 try
                    if dursMINIMUM>0 && dursMINIMUM~=123
                        
                        %restrict to those trials of interest!
                        all_onset = 1:size(Ppant_onsetBP,1);
                        all_offset = 1:size(Ppant_offsetBP,1);
                        
                        shrton= find(Ppant_onsetRTs<dursMINIMUM); %note switch
                        all_onset(shrton)=[];
                        keepON=all_onset;
                        
                        %same for offsets.
                        shrtoff=find(Ppant_offsetRTs<dursMINIMUM);
                        all_offset(shrtoff)=[];
                        keepOFF=all_offset;
                        
                        Ppant_onsetBP=Ppant_onsetBP(keepON,:);
                        Ppant_offsetBP=Ppant_offsetBP(keepOFF,:);
                        
                    elseif dursMINIMUM==123
                        
                        %select top third of disap durations.
                         %restrict to those trials of interest!
                        all_onset = 1:size(Ppant_onsetBP,1);
                        all_offset = 1:size(Ppant_offsetBP,1);
                        
                        
                        %find top third of distribution.
                        [ndurs,nid] = sort(Ppant_onsetRTs, 'descend');
                        %top third
                        keepON = nid(1:ceil(length(ndurs)/4));
                        
                        %same for offsets.
                        %find top third of distribution.
                        [ndurs,nid] = sort(Ppant_offsetRTs, 'descend');
                        %top third
                        keepOFF = nid(1:ceil(length(ndurs)/3));
                        
                        Ppant_onsetBP=Ppant_onsetBP(keepON,:);
                        Ppant_offsetBP=Ppant_offsetBP(keepOFF,:);
                        
                    end
                        
                        
                    
                    
                   
                    
                    % we need to resample the BP and SNR data, to equate across
                    % trial types    (ignoring nan)
                    
                    % currently at 100 fs. 
                    
                    % beware edge artefacts!
                   
                    
                    
                    % also apply participant level smoothing.
                    if ppantsmoothing==1
                        %smooth across trials, per ppant for SNR
                        
                        
                        % note, need to adjust for edge artefacts in smoothig
                        % process, repmat the end trials so smoothing doesnt
                        % contaminate.
                        
                        
                        for ionoff=1:4
                            
                            switch ionoff
                                case 1
                                    pData = Ppant_onsetSNR;
                                case 2
                                    pData = Ppant_offsetSNR;
                                case 3
                                    pData = Ppant_onsetBP;
                                case 4
                                    pData = Ppant_offsetBP;
                            end
                            
                            % first find out if there are nans, and remove.
                           if any(isnan(mean(pData,2)))
                                nantr= find(isnan(mean(pData,2)));
                                % for the trials with a nan. we need to
                                % remove.
                                pData(nantr,:)=[];
                            end
                               lc=size(pData,1);  
                                     
                            %discount the nans in data. (since these are
                            %sorted, should be at bottom).
                            %%
                            pData= pData(1:lc,:);
                            
                            %repmat edges to avoid edge artefacts!
                            repsize = 42; %n trials to repeat (top and bottom)
                            
                            %new data with expanded width.
                            pdataN = zeros(size(pData,1)+repsize*2, size(pData,2));
                            %top of image
                            % repeat the top trial, 
                            
%                             pdataN(1:repsize,:) = repmat(pData(1,:), [repsize,1]);
                            
                            % or reflect
                            pdataN(1:repsize,:) = flipud(pData(1:repsize,:));
                            
                            
                            %fill the remaining (central)trials.
                            pdataN((repsize+1):(lc+repsize),:) = pData;

                                %and bottom, repeat last.                               
%                             pdataN((lc+repsize):(lc+2*repsize-1),:) = repmat(pData(lc+repsize,:), [repsize,1]);
                        
                                %reflect works best.
                                pdataN((lc+repsize):(lc+2*repsize-1),:) = flipud(pData((lc-repsize+1:lc),:));



%                             subplot(311); imagesc(pData); % original
%                             subplot(312); imagesc(pdataN); % stretched
  
                            
                            % now resample this large set, but smooth from
                            % middle.
                            
                           pdataNre= resample(pdataN,160,size(pdataN,1));
                            
% subplot(313); imagesc(pdataNre)
%
                            %now safe to smooth
                            
                            pdataOUT= zeros(size(pdataNre));
                            for itime=1:size(pdataNre,2)
                                pdataOUT(:,itime)= smooth(pdataNre(:,itime),8);
                            end
 
                            % Now extract central component, to remove edge
                            % artegacts
                            pdataOUT = pdataOUT;
                            
%                              
% %                             %% san check:
                             figure(1); clf
                            subplot(221);imagesc(pData); title('original'); subplot(2,2,2); imagesc(pdataN); title('extended')                            
                            subplot(223); imagesc(pdataNre); title('resample ext.'); subplot(224); imagesc(pdataOUT); title('output all'); hold on;
                            plot([xlim], [20 20], ['k:'])
                            plot([xlim], [140 140], ['k:'])
%                             % keep this data, note these extensions mark
%                             % the flip end points.

                            pdataOUT = pdataOUT(20:140,:);
                            %%
%                             %%
                            switch ionoff
                                case 1
                                    Ppant_onsetSNR = pdataOUT;
                                case 2
                                    Ppant_offsetSNR = pdataOUT;
                                case 3
                                    Ppant_onsetBP = pdataOUT;
                                case 4
                                    Ppant_offsetBP = pdataOUT;
                                    
                            end
                            
                        end
                        
                    end
                    
                    
                    
                    
                    storeacrossPpant_onsetBP(icounter,:,:)=Ppant_onsetBP;
                    storeacrossPpant_onsetSNR(icounter,:,:)=Ppant_onsetSNR;
                    %                 storeacrossPpant_onsetRTs(icounter,:)=Ppant_onsetRTs;
%                     storeacrossPpant_onsetTOPO(icounter,:)=Ppant_onsetTOPO';
                    
                    
                    storeacrossPpant_offsetBP(icounter,:,:)=Ppant_offsetBP;
                    storeacrossPpant_offsetSNR(icounter,:,:)=Ppant_offsetSNR;
                    %                 storeacrossPpant_offsetRTs(icounter,:)=Ppant_offsetRTs;
%                     storeacrossPpant_offsetTOPO(icounter,:)=Ppant_offsetTOPO';
                    
                    
                    icounter=icounter+1;
%                 catch SETUP
%                     if ippant~=27
%                     rethrow(SETUP)
%                     end
%                 end
        end
        %%
%         for i=1:size(storeacrossPpant_offsetSNR,1);
%            subplot(4,4,i);
%            imagesc(squeeze(storeacrossPpant_offsetSNR(i,:,:)));
%             
%         end
        %%
        %save appropriately
        cd(basefol)
        cd('EEG')
        cd('GFX_EEG-RESS')
        
        
        
                savename=['GFX_PFIperformance_withSNR_' num2str(usehz) '_min' num2str(dursMINIMUM) '_RESS'];
            
                
        
        
        
        save(savename,...
            'storeacrossPpant_onsetBP','storeacrossPpant_offsetBP',...
            'storeacrossPpant_onsetSNR', 'storeacrossPpant_offsetSNR', 'timeidDYN')
        
        
        
    end
    
    
    
end
    





if job.erpimageacrossppants==1
    %% now plot across ppants.
    %reshape manually.
    
    %%
    getelocs
    rmvbase=0;
    
    
    counter=1;
    for hzis=[1,2,7]%1,2]
        cd(basefol)
        cd('EEG')
        
        usehz= peakfreqsare(hzis);
        loadname=['GFX_PFIperformance_withSNR_' num2str(usehz) '_min0_RESS'];
        cd('GFX_EEG-RESS')
        load(loadname);
                
               
        switch hzis
            case 1
                ylimsare = [1.1 1.6];
            case 2
                ylimsare = [2.5 2.9];
            case 7
                ylimsare = [0 .3];
                
        end
        % this is the scale of the catch image:
% switch hzis
%     case 1
%         ylimsare = [.8 2];
%     case 2
%         ylimsare = [2.6 3.4];
%     case 7
%         ylimsare = [0 .8];
% end
%         tgrm = timeidDYN;
        %%
        for itimezero=1:2%1:2
            
            
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
            
            acrSNR=squeeze(nanmean(useSNR,1));
%             acrSNR=squeeze(nanmean(useSNR(21,:,:),1));
%%            
% for i=1:size(useSNR,1);
%     subplot(4,4,i);
%     imagesc(squeeze(useSNR(i,:,:)));
%     
% end
  %%          
            
            %sort across all
            %sort trial index by longest RT first.
            %%
            figure(1)

            %%
            subplot(4,2,itimezero)
            colormap('jet')
            %         [sortedRTs, cid] = sort(acrRTs, 'descend');
            
            %
            imagesc(-3:1/60:3, 1:size(acrBP,1), acrBP);%(cid,:));
            %
            title({['Buttons Pressed'] }, 'fontsize', 20)
            c=colorbar;
            ylabel(c, 'Total')
            %                 ylabel('resampled catch trials')
            xlim([-2 2])
            set(gca, 'fontsize', 20, 'ytick', [])
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
            axis square
            %%
%             xlabel(['Time from ' ctype])
            
            %SNR
            if hzis~=7
            placeis= itimezero+2 + 2*(hzis-1);
            else
                placeis = itimezero+6;
            end
            subplot(4,2, placeis);
            
            %reorder SNR
            acrSNRsort=acrSNR;
            
            %%
%             %                 %smooth across trials
%             sm_snrgrm20=zeros(size(acrSNRsort));
%             for itime=1:size(acrSNRsort,2)
%                 sm_snrgrm20(:,itime)= smooth(acrSNRsort(:,itime),15);
%             end
            
            
            %     imagesc(acrSNR);
            imagesc(timeidDYN, 1:size(acrSNRsort,1), acrSNRsort);
            c=colorbar;

            ylabel(c, {['RESS log(SNR)']})

            if hzis~=7
            title([num2str(usehz) ' Hz (f' num2str(hzis) ')'])
            else
                title([num2str(usehz) ' Hz (f2-f1)'])
            end
%             title(['BG ' num2str(hzis) 'F'])
            
            hold on
            plot([0 0] , ylim, ['k:'], 'linewidth', 3)
            plot([0 0] , ylim, ['w:'], 'linewidth', 2)
            xlim([-2 2])
            set(gca, 'fontsize', 15, 'ytick', [])
            %                 subplot(4,2,8)
            %
%             %plot across ppant trace
%             ppantMeanSNR= squeeze(mean(useSNR,2));
%             
%             %adjust standard error as per COusineau(2005)
%             %confidence interval for within subj designs.
%             % y = x - mXsub + mXGroup,
%             x = ppantMeanSNR;
%             
%             mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
%             mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
%             
%             %for each observation, subjtract the subj average, add
%             %the group average.
%             NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
%             
%             %             %compute new stErr %which version?
%             %             stE = std(NEWdata)/sqrt(size(x,1));
%             %
%             %                 shadedErrorBar(tgrm-3, mean(ppantMeanSNR,1),stE,'k',[])
%             %                 set(gca, 'fontsize', 15)
%             %                 hold on;
%             %                 %%
            %
            %                 ylabel('mean SNR')
            %                 axis tight
            %                 plot([0 0] , ylim, ['k:'], 'linewidth', 1)
            if hzis==7
            xlabel(['Time from ' ctype])
            end
            ylabel({'Normalized'; 'trial count'})
                caxis([ylimsare])
                axis square
                
            set(gca, 'fontsize', 20)
            cd(basefol)            
            cd('Figures')
            
            set(gcf, 'color', 'w')
            shg
            
           colormap('viridis')
            shg
            counter=counter+1;
        end
        
    end
    %%
    shg
    %%
    cd([basefol filesep 'Figures' filesep 'GFX PFI trial by trial'])
     print('-dpng', ['PFI SNR summary all.png'])
end



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


if job.erpimage_catchPpantlevel==1
    %%
    rmvbase=0;
    snrmethod=2; % 1 for division method (MXC), 2 for dot product convolution.
    useSNRorhilbert=1;
    for ifol = allppants
        
        for ihz=[3:6]; %TGs1 TGs 2, 20hz 40 hz.
            
            
            usehz= peakfreqsare(ihz);
            
            
            
            icounter=1;
            
             cd(basefol)
        cd('EEG')
        cd(pdirs(ifol).name)
            load(['ppant_Catch_Epoched_RESS'])
            load('ppant_Catch_Epoched'); %for BP
            clearvars ppant_SNREEG*
            
            for itype=1:5
                switch itype
                    case 1 % pure catch onset
                    
                        switch ihz
                            case {1,3,5} % TGs
                                %adjust to 1,2,3
                                ind = ceil(ihz/2);
                                datatouse = ress_catchonsetTGs(ind,:,:);
                            case {2,4,6} % BGs
                                ind=ihz/2;
                                datatouse = ress_catchonsetBGs(ind,:,:);
                            case 7                                
                                datatouse = ress_catchonsetIMs(1,:,:);
                        end
                    BPstouse = catchonsetBPs;
                    ctype = 'catch onset visible';                    
                    case 2 %pure catch offset
                    
                        switch ihz
                            case {1,3,5} % TGs
                                %adjust to 1,2,3
                                ind = ceil(ihz/2);
                                datatouse = ress_catchoffsetTGs(ind,:,:);
                            case {2,4,6} % BGs
                                ind=ihz/2;
                                datatouse = ress_catchoffsetBGs(ind,:,:);
                            case 7
                                datatouse = ress_catchoffsetIMs(1,:,:);
                        end

                    BPstouse = catchoffsetBPs;
                    ctype = 'catch offset visible';
                    
                    case 3 % BP to catch
                        
                        switch ihz
                            case {1,3,5} % TGs
                                %adjust to 1,2,3
                                ind = ceil(ihz/2);
                                datatouse = ress_BPcatchonsetTGs(ind,:,:);
                            case {2,4,6} % BGs
                                ind=ihz/2;
                                datatouse = ress_BPcatchonsetBGs(ind,:,:);
                            case 7
                                datatouse = ress_BPcatchonsetIMs(1,:,:);
                        end

                    BPstouse = withincatchonsetBPs;
                    ctype = 'within catch onset BP';
                    
                    case 4
                        
                        switch ihz
                            case {1,3,5} % TGs
                                %adjust to 1,2,3
                                ind = ceil(ihz/2);
                                datatouse = ress_BPcatchoffsetTGs(ind,:,:);
                            case {2,4,6} % BGs
                                ind=ihz/2;
                                datatouse = ress_BPcatchoffsetBGs(ind,:,:);
                            case 7
                                datatouse = ress_BPcatchoffsetIMs(1,:,:);
                        end

                    BPstouse = postcatchoffsetBPs;
                    ctype = 'post offset BP';
                    
                    case 5 % invisible catch
                        switch ihz
                            case {1,3,5} % TGs
                                %adjust to 1,2,3
                                ind = ceil(ihz/2);
                                datatouse = ress_invisiblecatchonsetTGs(ind,:,:);
                            case {2,4,6} % BGs
                                ind=ihz/2;
                                datatouse = ress_invisiblecatchonsetBGs(ind,:,:);
                            case 7
                                datatouse = ress_invisiblecatchonsetIMs(1,:,:);
                        end

                    BPstouse = invisibleonsetBPs;
                    ctype = 'invisible catch';
                    
                        
                    
                    
                end
                datatouse=squeeze(datatouse);
                
                
                %plot the topo for sanity check: (pre disap)
                windcheck= [-3 -0.1];
                tidx=dsearchn(timeidDYN', [windcheck]');
                
                
                
                %restrict to certain PFI duration?
                allt= 1:size(datatouse,1);
                %                 baddurs = find(RTstouse<30);
                %remove those trials.
                %                 allt(baddurs)=[];
               
                
                %%
                %%
                clf
                colormap('viridis')
                %plot BP for comparison
                subplot(3,2, [1 3])
                %sort trial index by longest RT first.
                %                 [sortedRTs, cid] = sort(RTstouse(allt), 'descend');
                
                %sort by accum PFI in window after onset.
                switch itype
                    case {1,3,5}
                    accumPeriod = sum(BPstouse(:,180:end),2);
                    case {2,4}
                    accumPeriod = sum(BPstouse(:,1:180),2);
                    
                end
                
                [checkBP, cid] = sort(accumPeriod, 'descend');
                
                
                % MOVE nan to bottom of order (not top)
                nanind= find(isnan(checkBP));
                %new start point
                cid=[cid(length(nanind)+1:end) ];%  ; cid(nanind)];
                
                
                %rearrange.
                BPstouse=BPstouse(cid,:);
                datais=squeeze(datatouse(cid,:));
%                 RTstouse=RTstouse(cid);
                %%

                tBP=-3:1/60:3;
                imagesc(-3:1/60:3,1:length(cid),BPstouse)
                c=colorbar;
                ylabel(c, 'buttons pressed')
                %                 set(gca, 'ytick', 1:length(cid), 'yticklabel', round(sortedRTs./60,2))
                ylabel('PFI events')
                xlabel('Time [sec]')
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                title({['sorted BP data for ' ctype ','];['ppant' num2str(ifol)]})
                set(gca, 'fontsize', 15)
                %%
                % show mean over time.
                subplot(3,2,5)
                plot(tBP, nanmean(BPstouse,1),['k-'], 'linewidth', 3)
                %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                set(gca, 'fontsize', 15)
                hold on;
                ylabel('nanmean BP ')
                %
                xlabel('Time secs')
                axis tight
                plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                
                
                set(gca, 'fontsize', 15)
                %%
                
                if useSNRorhilbert~=1
                %% now plot dynamic RESS:
                
                peakwidt=2;
                 %filter, but inlcude side lobes for temporal dynamics
                fdatAt = filterFGx(datais,srate,usehz,peakwidt,0);                
                
                % take absolute of hilbert transform, to show dynamic SSEP
                dynSSEP= abs(hilbert(fdatAt'))';
                
                
%                 %TAKE mean per trial, careful with baseline subtraction
%                 timeid = [0:1/srate:epochdur];
%                 timeid= timeid-3;
%                 if dirsi==1 %increasing PFI / Target absence
%                     if ifreq<5 %targets, so low point post event.
%                         windcheck= [2 3];
%                     else %BG hz, so
%                         windcheck= [-3 -1];
%                     end
%                 elseif dirsi==0 %oppoite case, decreasing PFI/ targets are appearing
%                     if ifreq<5 %targets, so low point pre event.
%                         windcheck= [-3 -2];
%                     else %BG hz, so
%                         windcheck= [2 3];
%                     end
%                 end

%                 
                else % use SNR
                     
                            
                        
                            
%                             if ifreq==7
%                                 pause
%                             end
                            [sgrm ,tgrm, fgrm] = mtspecgramc(datais', movingwin, param_spcgrm);
                            %% %% Two SNR methods:
                            
                            %%
                            snr_sgrm =zeros(size(sgrm));
                            
                            tmps=sgrm;
                            %adjust for HBW 
                                k = ((param_spcgrm.tapers(1,1)));%
                                hbw = (k+1)./ (2.*[movingwin(1,1)]);
                                neighb = 2; %hz
                                
                                %space out kernel appropriately
                                distz = dsearchn(fgrm', [hbw, neighb, neighb*2+hbw]');
                                
                                % we don't want odd numbers, messes with
                                % calculations, so round allup to even.
                                ods = find(mod(distz,2)~=0);
                                %adjust.
                                distz(ods)= distz(ods)+1;
                                
                                %make kernel
                                kernelneighb = -1*(repmat(1/(distz(2)*2), [1, distz(2)]));
                                kernelskip = zeros(1,distz(1));
                                
                                kernelw= [kernelneighb, kernelskip, 1, kernelskip, kernelneighb];
                                if sum(kernelw)~=0
                                    error('')
                                end
                            if snrmethod==1%method 1, using division.
                                
                               
                                
                                
                                
                                % loop over frequencies and compute SNR
                                numbins = distz(2);
                                skipbins = distz(1);
                                snr_sgrm= zeros(size(sgrm));
                                for itrial = 1:size(sgrm,3)
                                for itime=1:size(tmps,1)
                                    for hzi=numbins+1:length(fgrm)-numbins-1
                                        numer = tmps(itime,hzi, itrial);
                                        denom = nanmean( tmps(itime,[hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
                                        snr_sgrm(itime,hzi,itrial) = numer./denom;
                                        
                                    end
                                end
                                end
                                
                            else
                                %                             %conv SNR
                                %                             snr_sgrm =zeros(size(sgrm));
                                %
                                
                                %%
%                                 kernelw= [-1/8 -1/8 -1/8 -1/8 0 0 1 0 0 -1/8 -1/8 -1/8 -1/8];
                                
                                %snr on trials or
                                for itrial=1:size(sgrm,3)
                                    %compute SNR
                                    tmps= squeeze(sgrm(:,:,itrial));
                                    
                                    for itime= 1:size(tmps,1)
                                        checkput = conv(log(tmps(itime,:)), kernelw,'same');
                                        if ~isreal(checkput)
                                            snr_sgrm(itime,:,itrial)= nan(1, size(tmps,2));
                                        else
                                            snr_sgrm(itime,:,itrial)= conv(log(tmps(itime,:)), kernelw,'same');
                                        end
                                    end
                                end
                            
%                              %% sanity check:
%                               figure(2);
%                               %plot specgrm and spctrm.
%                               subplot(221)
%                             imagesc(tgrm, fgrm, squeeze(mean(log(sgrm),3))'); colorbar; caxis([0 2])
%                             subplot(223)
%                             mtrials=squeeze(mean(log(sgrm),3));
%                             plot(fgrm, squeeze(mean(mtrials,1)));
%                             %plot specgrm and spctrm.
%                               subplot(222)
%                             imagesc(tgrm, fgrm, squeeze(mean((snr_sgrm),3))'); colorbar; caxis([0 2])
%                             subplot(224)
%                             mtrials=squeeze(mean((snr_sgrm),3));
%                             plot(fgrm, squeeze(mean(mtrials,1)));
%                             hold on; plot(kernelw)
                            end
                            %%
                            
                            %store just stim freq.
                            [~, hzid]= min(abs(fgrm-usehz));
                            %
% %                             %reduce size.
                            dynSSEP=squeeze(snr_sgrm(:,hzid,:))';
% %                             
                          timeidDYN=tgrm-3;

                            
                    
                    
                end
%                         
%                %remove base
if rmvbase==1
    tidx=dsearchn(timeidDYN', [bsrem]');
    for itrial=1:size(dynSSEP,1)
        tmp= dynSSEP(itrial,:);
        mtmp = mean(tmp(tidx(1):tidx(2)));
        %                    mtmp = mean(tmp);
        
        dynSSEP(itrial,:)= tmp - repmat(mtmp, 1, length(tmp));
        
    end
    
end
                
               %
               
                %%\
                snrgrm20=dynSSEP;
                % smooth across trials
                sm_snrgrm20=zeros(size(snrgrm20));
                for itime=1:size(snrgrm20,2)
                    sm_snrgrm20(:,itime)= smooth(snrgrm20(:,itime),5);
                end
                %%
                hold on
                subplot(3,2,[2 4]);
                hold on
                imagesc(timeidDYN,  1:size(dynSSEP,1), sm_snrgrm20);
                c=colorbar;
                if useSNRorhilbert==1
                ylabel(c, 'log(SNR)')
                else
                ylabel(c, 'abs(hilbert)')
                end
                
                xlabel('Time secs')
%                 caxis([-1*max(max(sm_snrgrm20)) max(max(sm_snrgrm20))])
                try caxis([0  max(max(sm_snrgrm20))-1])%                 
                catch
                    caxis([0  max(max(sm_snrgrm20))])%                 
                end
                    
% caxis([ 1 5])
                axis tight
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                title({['sorted ' num2str(usehz) 'Hz RESS data (smoothed) for  ' ctype ','];['ppant' num2str(ifol)]})
                set(gca, 'fontsize', 15, 'yticklabel', [])
                
                %         ylim([0 20])
                %% show mean over time.
                subplot(3,2,6)
                plot(timeidDYN, nanmean(snrgrm20,1),['k-'], 'linewidth', 3)
                %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                set(gca, 'fontsize', 15)
%                 ylim([1 3])
                hold on;
                ylabel('RESS logSNR')
                %%
                axis tight
                plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                shg
                xlabel('Time secs')
                set(gcf,'color', 'w');
                
%                         ylim([-.1 10])
                
                %%
                shg
                cd([basefol filesep 'Figures' filesep 'Participant Catch summaries'])
                %%
                %%
                print('-dpng', ['figure_' ctype '_summary_BPandSSVEP_' num2str(usehz) '_RESS_subj' num2str(ifol) ' snr method ' num2str(snrmethod) '.png']);
                
                
                switch itype
                    case 1 %catch onset
                        catch_onsetSNR_RESS_BPs=BPstouse;
                        catch_onsetSNR_RESS=snrgrm20; %sorted (un smoothed)

                    case 2 %catch offset
                        catch_offsetSNR_RESS_BPs=BPstouse;
                        catch_offsetSNR_RESS=snrgrm20; %sorted (un smoothed)
                    case 3
                        BPcatch_onsetSNR_RESS_BPs=BPstouse;
                        BPcatch_onsetSNR_RESS=snrgrm20; %sorted (un smoothed)
                    case 4
                        BPcatch_offsetSNR_RESS_BPs=BPstouse;
                        BPcatch_offsetSNR_RESS=snrgrm20; %sorted (un smoothed)
                    case 5
                        invisiblecatch_onsetSNR_RESS_BPs=BPstouse;
                        invisiblecatch_onsetSNR_RESS=snrgrm20; %sorted (un smoothed)
                end
            end
            
            
            %%
            
                    savename=['Catch_performance_withSNR_' num2str(usehz) '_RESS'];
            
            cd(basefol)
            cd('EEG')
            cd(pdirs(ifol).name)
            %%
            save(savename,... 
                'catch_onsetSNR_RESS_BPs','catch_onsetSNR_RESS',...
                'catch_offsetSNR_RESS_BPs','catch_offsetSNR_RESS',...
            'BPcatch_onsetSNR_RESS_BPs','BPcatch_onsetSNR_RESS',...    
            'BPcatch_offsetSNR_RESS_BPs','BPcatch_offsetSNR_RESS',...
            'invisiblecatch_onsetSNR_RESS_BPs','invisiblecatch_onsetSNR_RESS',...
            'timeidDYN', 'param_spcgrm', 'kernelw')
            
            
        end
    end
end
%%

if job.concaterpdataforcatchacrossppants==1
    %%
    
    dursMINIMUM= 0; %1 second. % or change to  %this works.
%     dursMINIMUM=123; % selects top third!.
    
    ppantsmoothing=1; % average across participants, after smoothing, or no.
    
    for ihz=[3:6];
        usehz=peakfreqsare(ihz);
        
            loadname=['Catch_performance_withSNR_' num2str(usehz) '_RESS'];
        
   storeacrossPpant_catchEVENTS_BPs = zeros(length(allppants), 5, 132, 361);
   storeacrossPpant_catchEVENTS_SNR = zeros(length(allppants), 5, 132, length(Nt));
        icounter=1;
        
        for ippant = allppants
             cd(basefol)
             cd('EEG')
             cd([pdirs(ippant).name])
                
                %onset types
                
                load(loadname)
%                 Ppant_onsetSNR=Ppant_onsetSNR';
%                 Ppant_offsetSNR=Ppant_offsetSNR';
                
                %restrict to only 'long' PFIs if needed.
%                 try
                    if dursMINIMUM>0 && dursMINIMUM~=123
                        
                        %restrict to those trials of interest!
                        all_onset = 1:size(Ppant_onsetBP,1);
                        all_offset = 1:size(Ppant_offsetBP,1);
                        
                        shrton= find(Ppant_onsetRTs<dursMINIMUM); %note switch
                        all_onset(shrton)=[];
                        keepON=all_onset;
                        
                        %same for offsets.
                        shrtoff=find(Ppant_offsetRTs<dursMINIMUM);
                        all_offset(shrtoff)=[];
                        keepOFF=all_offset;
                        
                        Ppant_onsetBP=Ppant_onsetBP(keepON,:);
                        Ppant_offsetBP=Ppant_offsetBP(keepOFF,:);
                        
                    elseif dursMINIMUM==123
                        
                        %select top third of disap durations.
                         %restrict to those trials of interest!
                        all_onset = 1:size(Ppant_onsetBP,1);
                        all_offset = 1:size(Ppant_offsetBP,1);
                        
                        
                        %find top third of distribution.
                        [ndurs,nid] = sort(Ppant_onsetRTs, 'descend');
                        %top third
                        keepON = nid(1:ceil(length(ndurs)/4));
                        
                        %same for offsets.
                        %find top third of distribution.
                        [ndurs,nid] = sort(Ppant_offsetRTs, 'descend');
                        %top third
                        keepOFF = nid(1:ceil(length(ndurs)/3));
                        
                        Ppant_onsetBP=Ppant_onsetBP(keepON,:);
                        Ppant_offsetBP=Ppant_offsetBP(keepOFF,:);
                        
                    end
                        
                        
                    
                    
                   
                    
                    % we need to resample the BP and SNR data, to equate across
                    % trial types    (ignoring nan)
                    
                    % currently at 100 fs. 
                    
                    % beware edge artefacts!
                   
                    
                    
                    % also apply participant level smoothing.
                    if ppantsmoothing==1
                        %smooth across trials, per ppant for SNR
                        
                        
                        % note, need to adjust for edge artefacts in smoothig
                        % process, repmat the end trials so smoothing doesnt
                        % contaminate.
                        
                        
                        for ionoff=1:10 %all SNR and BP combos.
                            
                            switch ionoff
                                case 1
                                    pData = catch_onsetSNR_RESS;
                                    storeid=1;
                                case 2
                                    pData = catch_onsetSNR_RESS_BPs;
                                    
                                case 3
                                    pData = catch_offsetSNR_RESS;
                                    storeid=2;
                                case 4
                                    pData = catch_offsetSNR_RESS_BPs;
                                case 5
                                    pData = BPcatch_onsetSNR_RESS;
                                    storeid=3;
                                case 6
                                    pData = BPcatch_onsetSNR_RESS_BPs;
                                case 7
                                    pData = BPcatch_offsetSNR_RESS;
                                    storeid=4;
                                case 8
                                    pData = BPcatch_offsetSNR_RESS_BPs;
                                case 9
                                    pData = invisiblecatch_onsetSNR_RESS;
                                    storeid=5;
                                case 10
                                    pData = invisiblecatch_onsetSNR_RESS_BPs;
                            end
                            
                            % first find out if there are nans, and remove.
                            if any(isnan(mean(pData,2)))
                                nantr= find(isnan(mean(pData,2)));
                                % for the trials with a nan. we need to
                                % remove.
                                pData(nantr,:)=[];
                            end
                               lc=size(pData,1);      
                            %discount the nans in data. (since these are
                            %sorted, should be at bottom).
                            %%
                         %%
                            pData= pData(1:lc,:);
                            
                            %repmat edges to avoid edge artefacts!
                            repsize = 4; %n trials to repeat (top and bottom)
                            
                            %new data with expanded width.
                            pdataN = zeros(size(pData,1)+repsize*2, size(pData,2));
                            %top of image
                            % repeat the top trial, 
                            
%                             pdataN(1:repsize,:) = repmat(pData(1,:), [repsize,1]);
                            
                            % or reflect
                            pdataN(1:repsize,:) = flipud(pData(1:repsize,:));
                            
                            
                            %fill the remaining (central)trials.
                            pdataN((repsize+1):(lc+repsize),:) = pData;

                                %and bottom, repeat last.                               
%                             pdataN((lc+repsize):(lc+2*repsize-1),:) = repmat(pData(lc+repsize,:), [repsize,1]);
                        
                                %reflect works best.
                                pdataN((lc+repsize):(lc+2*repsize-1),:) = flipud(pData((lc-repsize+1:lc),:));



%                             subplot(311); imagesc(pData); % original
%                             subplot(312); imagesc(pdataN); % stretched
  
                            
                            % now resample this large set, but smooth from
                            % middle.
                            
                           pdataNre= resample(pdataN,160,size(pdataN,1));
                            
% subplot(313); imagesc(pdataNre)
%
                            %now safe to smooth
                            
                            pdataOUT= zeros(size(pdataNre));
                            for itime=1:size(pdataNre,2)
                                pdataOUT(:,itime)= smooth(pdataNre(:,itime),8);
                            end
 
                            
%                              
% %                             %% san check:
%                              figure(1); clf
%                             subplot(221);imagesc(pData); title('original'); subplot(2,2,2); imagesc(pdataN); title('extended')                            
%                             subplot(223); imagesc(pdataNre); title('resample ext.'); subplot(224); imagesc(pdataOUT); title('output all'); hold on;
%                             plot([xlim], [14 14], ['k:'])
%                             plot([xlim], [145 145], ['k:'])
%                             % keep this data, note these extensions mark
%                             % the flip end points.
% Now extract central component, to remove edge
                            % artegacts
                            
                            pdataOUT = pdataOUT(14:145,:);
%% %                             %% san check:
%                              figure(1); clf
%                             subplot(211);imagesc(pData);
%                             subplot(212); imagesc(pdataOUT)
                            %%
%                             %%
%                             switch ionoff
%                                 case 1
%                                     catch_onsetSNR_RESS= pdataOUT;
%                                 case 2
%                                      catch_onsetSNR_RESS_BPs= pdataOUT;
%                                 case 3
%                                      catch_offsetSNR_RESS= pdataOUT;
%                                 case 4
%                                     catch_offsetSNR_RESS_BPs= pdataOUT;
%                                 case 5
%                                     BPcatch_onsetSNR_RESS= pdataOUT;
%                                 case 6
%                                     BPcatch_onsetSNR_RESS_BPs= pdataOUT;
%                                 case 7
%                                     BPcatch_offsetSNR_RESS= pdataOUT;
%                                 case 8
%                                     BPcatch_offsetSNR_RESS_BPs= pdataOUT;
%                                 case 9
%                                     invisiblecatch_onsetSNR_RESS= pdataOUT;
%                                 case 10
%                                     invisiblecatch_onsetSNR_RESS_BPs= pdataOUT;
%                                     
%                             end
                              if mod(ionoff,2)==0 %even numbers are BP
                    storeacrossPpant_catchEVENTS_BPs(icounter,storeid,:,:)=pdataOUT;
                              else
                    storeacrossPpant_catchEVENTS_SNR(icounter,storeid,:,:)=pdataOUT;
                              end
                        end
                        
                    end
                    
                    
                    icounter=icounter+1;
        end
        %%
%         for i=1:size(storeacrossPpant_catchEVENTS_SNR,1)
%            subplot(4,4,i);
%            imagesc(squeeze(storeacrossPpant_catchEVENTS_SNR(i,3,:,:)));
%             
%         end
        %%
        %save appropriately
        cd(basefol)
        cd('EEG')
        cd('GFX_EEG-RESS')
        
        DIMSare = {'catch onset', 'catch offset', 'BPcatchonset', 'BPcatchoffset', 'invisibleonset'};
        
                savename=['GFX_Catch_performance_withSNR_' num2str(usehz) '_min' num2str(dursMINIMUM) '_RESS'];
            
                
        
        
        
        save(savename,...
            'storeacrossPpant_catchEVENTS_SNR','storeacrossPpant_catchEVENTS_BPs',...
            'DIMSare', 'timeidDYN')
        
        
        
    end
    
end

if job.erpimage_forcatchacrossppants==1
    %% now plot across ppants.
    %reshape manually.
    
    %%
    getelocs
    rmvbase=0;
    
    clf
    counter=1;
     for itype=[1:5]
         
         clf
    for hzis=[1,2,7]%1,2]
        cd(basefol)
        cd('EEG')
        
        usehz= peakfreqsare(hzis);
        loadname=['GFX_Catch_performance_withSNR_' num2str(usehz) '_min0_RESS'];
        cd('GFX_EEG-RESS')
        load(loadname);
                
               
        switch hzis
            case 1
                ylimsare = [.8 2];
            case 2
                ylimsare = [2.6 3.4];
                case 7
                ylimsare = [0 .8];
        end
%         tgrm = timeidDYN;
        %%
       
                        
        chtype=DIMSare{itype};
        useBP= squeeze(storeacrossPpant_catchEVENTS_BPs(:,itype,:,:));
        useSNR= squeeze(storeacrossPpant_catchEVENTS_SNR(:,itype,:,:));
        
            
            
            
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
            
            acrSNR=squeeze(nanmean(useSNR,1));
%             acrSNR=squeeze(nanmean(useSNR(21,:,:),1));
%%            
% figure(10)
% for i=1:size(useSNR,1);
%     subplot(4,4,i);
%     imagesc(squeeze(useSNR(i,:,:)));
%     
% end
  %%          
            
            %sort across all
            %sort trial index by longest RT first.
            %%
            figure(1)
            

            %%
            subplot(4,1,1)
            colormap('viridis')
            %         [sortedRTs, cid] = sort(acrRTs, 'descend');
            
            %
            imagesc([-3:1/60:3], 1:size(acrBP,1), acrBP);%(cid,:));
            %
            title({['Buttons Pressed'] }, 'fontsize', 20)
            c=colorbar;
            ylabel(c, 'Total')
            %                 ylabel('resampled catch trials')
            xlim([-2 2])
            set(gca, 'fontsize', 20, 'ytick', [])
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
            axis square
            %%
%             xlabel(['Time from ' chtype])
            
           
            if hzis~=7
            subplot(4,1, 1+hzis);
            else
            subplot(4,1, 4);
            end
            
            %reorder SNR
            acrSNRsort=acrSNR;
            
            %%
%             %                 %smooth across trials
%             sm_snrgrm20=zeros(size(acrSNRsort));
%             for itime=1:size(acrSNRsort,2)
%                 sm_snrgrm20(:,itime)= smooth(acrSNRsort(:,itime),15);
%             end
            
            
            %     imagesc(acrSNR);
            imagesc(timeidDYN, 1:size(acrSNRsort,1), acrSNRsort);
            c=colorbar;

            ylabel(c, {['RESS log(SNR)']})
if hzis~=7
            title([num2str(usehz) ' Hz (f' num2str(hzis) ')'])
else
            title([num2str(usehz) ' Hz (f2-f1)'])
end
%             title(['BG ' num2str(hzis) 'F'])
            
            hold on
            plot([0 0] , ylim, ['k:'], 'linewidth', 3)
            plot([0 0] , ylim, ['w:'], 'linewidth', 2)
            xlim([-2 2])
            set(gca, 'fontsize', 15, 'ytick', [])
            %                 subplot(4,2,8)
            %
%      
if hzis==7
    xlabel(['Time from ' chtype])
end
            ylabel({'Normalized'; 'trial count'})
                caxis([ylimsare])
                axis square
                axis tight
            set(gca, 'fontsize', 20)
            cd(basefol)            
            cd('Figures')
            
            set(gcf, 'color', 'w')
            shg
            
           colormap('viridis')
            shg
            counter=counter+1;
        end
        
    
    %%
    shg
    cd([basefol filesep 'Figures' filesep 'GFX catch trial by trial'])
    %%    
     print('-dpng', ['Catch SNR summary for f1 and f2, ' chtype '.png'])
     end
end

if job.gradedchangesinERPimage_Catch==1
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
    for hzis=[1, 2, 7]
        cd(basefol)
        cd('EEG')
        cd('GFX_EEG-RESS')
        usehz= peakfreqsare(hzis);
        switch hzis
            case 1
                
                ylimsare=[.5, 2];
                col='b';
                flickis = 'Target';
            case 2
%                 ylimsare=[2.5, 3.2];
                ylimsare=[1.9, 3.4];
                col=[.2 .2 .2];
                flickis = 'Background';
                
                case 7
%                 ylimsare=[2.5, 3.2];
                ylimsare=[.1, .8];
                col='m';
                flickis = 'Intermodulation';
        end
        load(['GFX_Catch_performance_withSNR_' num2str(usehz) '_min0_RESS'])
%         tgrm = timeidDYN;
        %%
              
      
%         tgrm = timeidDYN;
        %%
       
      for itype = 4     % 1:5              
        
          DIMSare{1} = 'catch onset';
          DIMSare{2} = 'catch offset';
          DIMSare{3} = 'reporting catch onset';
          DIMSare{4} = 'reporting catch offset';
          DIMSare{5} = 'invisible catch onset';
          
          
          
          chtype=DIMSare{itype};
        useBP= squeeze(storeacrossPpant_catchEVENTS_BPs(:,itype,:,:));
        useSNR= squeeze(storeacrossPpant_catchEVENTS_SNR(:,itype,:,:));
        
            
            
            
            
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
            title({[ upper(chtype(1)) num2str(chtype(2:end))]}, 'fontsize', 15)
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
%             xlabel(['Time from ' chtype])
            
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
            title([flickis ' ' num2str(usehz) ' Hz (f' num2str(hzis) ')'])
else
            title([flickis ' ' num2str(usehz) ' Hz (f2-f1)'])
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
            
        
            xlabel(['Time from ' chtype])
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
            
            if mod(itype,2)~=0 % odd nubers are for disappearances (
            BPidx=181:361;
            else
                BPidx=1:180;
            end
                groupedBEHflpd = fliplr(squeeze(nanmean(acrBP(:,BPidx),2))'); %takes mean over all SNR timepoints.
            %%
            groupednPFISNR_thirds=nan(length(allppants),3);
            
            for ippant = 1:size(useBP,1); 
%                % per ppant, collect trial indices:
                ppant_BEH= squeeze(mean(useBP(ippant,:,BPidx),3)); %takes mean over all BP time points.
                ppant_SNR=squeeze(mean(useSNR(ippant,:,:),3)); %takes mean over all SNR timepoints.
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
xlabel('amount of Catch');
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
%             
            dataan= reshape(groupednPFISNR, [nppants*nconds,1]);
            conds = [ones(nppants,1); ones(nppants,1)*2;ones(nppants,1)*3; ones(nppants,1)*4];
                subs = repmat([1:length(allppants)]', [nconds,1]);
%             %RMANOVA?
%             [anret]=rmanova(dataan, conds, subs);
% %%
            %LME to be sure:            
          tblA=table(dataan, conds, subs);
            
             UNterm1='dataan ~  conds + (1|subs)';
%              UNterm2='dataan ~  conds';
             UNterm2='dataan ~  1+ (1|subs)';
             %creating models:
             lmeUN =fitlme(tblA, UNterm1);
             lmeUN2 =fitlme(tblA, UNterm2);
             
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
xlabel(['Time from ' chtype])
ylabel('RESS log(SNR)')
axis tight
% ylim([ylimsare(1) ylimsare(2)])

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
    cd('GFX Catch trial by trial')
    shg
     print('-dpng', ['Catch graded SNR summary graded, ' chtype '.png'])
    
    
end


if job.BPandSSVEPtimecourseacrossppants_group==1

    
    getelocs
%     clf
    
    rmvbase=0;
    checksigON=1;
    checkcluster=1;
    
%     clf
    plcount=1;
    legendprint=[];
    cd(basefol)
    cd('EEG')
    cd('GFX_EEG-RESS')
    %
    %     cd('newplots-MD')
    colsare={'b' , 'k',[],[],[],[],'m'}; % blue for tg, black for BG.
    hzare={'TG (f1)' , 'BG (f2)',[],[],[],[],'IM (f2-f1)'}; % blue for tg, black for BG.
    linetypes={'--', '-', [],[],'--'};
    legendis=[];
    lgc=1;
%     clf
figure(1); hold on;
clearvars yl;
    for hzis=2%[1,2]%[1,2]%
        usehz=peakfreqsare(hzis);
        
        load(['GFX_PFIperformance_withSNR_' num2str(usehz) '_min0_RESS'])
        
        
        
        tbase = timeidDYN;
        col=colsare{hzis};
        hzp= hzare{hzis};
        
        ttestdata=[];
        %         legendprint=[];
        for itimezero=1:2%1:2
            
            switch itimezero
                case 1
                    useBP=storeacrossPpant_onsetBP;
                    useSNR=storeacrossPpant_onsetSNR;
                    
                    
                    chtype='target invisible';
                    %                     chtype='button press';
                    
                    %                     col=[0 .5 0];
                    linet='--';
                case 2
                    
                    useBP=storeacrossPpant_offsetBP;
                    useSNR=storeacrossPpant_offsetSNR;
                    
                    
                    
                    
                    chtype='target visible';
%                     chtype='time from subjective report';
                    %                     chtype='button release';
                    
                    
                    %                     col='r';
                    bsrem=[-3 -1];
                    linet='-';
                    
                case 3
                    useSNR=storeacrossPpant_onsetSNR+storeacrossPpant_offsetSNR;
                    chtype = 'subjective report';
                    linet=':';
                    
            end
            
          figure(1);
            
            used=useSNR;
            
            hold on
            
            
            %plot across ppant trace
            ppantMeanSNR= squeeze(mean(used,2));
%             subplot(2,1,itimezero);
%             plot(ppantMeanSNR');
            
            
            %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            
            x = ppantMeanSNR;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
            
            stE = std(NEWdata)/sqrt(length(allppants));
            
            
            if rmvbase==1
                ppantMeanSNR=ppantMeanSNR-mean(ppantMeanSNR(:));
            end
            %plot
                sh=shadedErrorBar(timeidDYN, mean(ppantMeanSNR,1),stE,[],[1]);
                
            sh.mainLine.LineWidth=3;
            
            sh.mainLine.Color = col;
            sh.mainLine.LineStyle=linet;
            sh.patch.FaceColor = col;
            sh.edge(1).Color = col;
            sh.edge(2).Color = col;
                      %
            
            xlabel(['Time from ' chtype ' [s]'])%           
%             xlabel(['Time from button release [s]'])
% xlabel('Time from reporting catch')

            if itimezero==3
                
                ylabel({['\Delta RESS log(SNR)']})
                
            else
                
                if rmvbase~=1
                    ylabel({['RESS log(SNR), ' hzp]})
                else
                    ylabel({['mean subtracted'];['RESS log(SNR)']})
                end
                    
            end
            
            
            
            set(gca, 'fontsize', 25)
            
       
            
            
            
            set(gcf, 'color', 'w')
            
            %store output for legend/stats
            ttestdata(lgc,:,:) = ppantMeanSNR;
            legendprint(lgc)=sh.mainLine;
            %legend based on type
            legendis = [legendis {chtype}];
            
            plcount=plcount+1;
            
             
            lgc=lgc+1;
            
        end % by type
        %legend based on Hz
%             legendis = [legendis {chtype}];
      
    end
         
%             xlim([timeidDYN(1) timeidDYN(end)])
               axis tight
            if rmvbase==1
                ylim([-.1 .15])
                 hold on;
                plot(xlim, [0 0], ['k:'])
            end
            if itimezero==3                
                plot(xlim, [0 0], ['k:'])
            end
                
    %%
         
            
    % %%%%%%%% END OF PLOTTING, SIG TESTS BELOW
    if checksigON==1
        %check for sig
        pvals=zeros(1,size(ttestdata,3));
        tvals=zeros(1,size(ttestdata,3));
        for itime = 1:size(ttestdata,3)
            
            try [h,pvals(itime),~,stat]=ttest(ttestdata(1,:,itime), ttestdata(2,:,itime));
                shuffType=1;
            catch
                [h,pvals(itime),~,stat]=ttest(ttestdata(plcount-1,:,itime)); %compares to zero.
                shuffType=2; %whether or not to skip the non-parametric test for sig.
            end
            
            tvals(itime)= stat.tstat;
        end
        sigs=find(pvals<.05);
        %
        %             %perform cluster based correction.
        if length(sigs)>2 &&checkcluster==1
            % find biggest cluster:
            %finds adjacent time points
            vect1 = diff(sigs);
            v1 = (vect1(:)==1);
            d = diff(v1);
            clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
            %grab largest
            %                 ignore bad points.
            
            
            % find biggest cluster:
            %finds adjacent time points
            sigs = find(pvals<.05);
            
            vect1 = diff(sigs);
            v1 = (vect1(:)==1);
            d = diff(v1);
            clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
            [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
            
            %
            for icl=1:size(clusterSTandEND,1)
                
                %start and end are now:
                % change icl to maxClust if we only want the largest
                % cluster.
                STC=sigs(clusterSTandEND(icl,1));
                ENDC=sigs(clusterSTandEND(icl,2)+1);
                checktimes =STC:ENDC;
                observedCV = sum(abs(tvals(checktimes)));
                % now shuffle condition labels to see if this cluster is
                % sig (compared to chance).
                
                
                
                %all trials
                alltr=reshape(ttestdata, [length(allppants)*2, length(Nt)]);
                alltrialsvec = 1:size(alltr,1);
                
                nshuff=2000;
                
                sumTestStatsShuff = zeros(1,nshuff);
                for irand = 1:nshuff
                    
                    if shuffType==1 %null is that no condition differences.
                        
                        
                        shD=zeros(size(ttestdata));
                        %since this is a within subjects design, we permute
                        %the subjet specific averages within each subject
                        %(as per Maris & Oostenveld (2007).
                        
                        % for each subject, randomly permute the averages.
                        %(Dsub1cond1,Datasub1cond2)
                        for ippant = 1:size(ttestdata,2)
                            
                            if mod(randi(100),2)==0 %if random even number
                                shD(1,ippant,:) = ttestdata(1,ippant,:); % all time points.
                                shD(2,ippant,:) = ttestdata(2,ippant,:); % all time points.
                            else
                                shD(1,ippant,:) = ttestdata(2,ippant,:); % all time points.
                                shD(2,ippant,:) = ttestdata(1,ippant,:); % all time points.
                            end
                            
                            %                             shD(ipartition,ippant,:) = pdata;
                        end
                        
                    else %null is that there are no temporal coincident sig values.
                        for ipartition = 1:2
                            for ippant = 1:size(ttestdata,2)
                                for itime=1:length(checktimes)
                                    
                                    %take random timepoint.
                                    pdata = ttestdata(1,ippant, randi(size(ttestdata,3)));
                                    
                                    
                                    shD(ipartition,itime,ippant) = pdata;
                                end
                            end
                        end
                    end
                    
                    %                     tvalspertimepoint = zeros(1,length(checktimes));
                    %%
                    % figure(3); clf; plot(squeeze(mean(shD(1,:,:),2))); hold on
                    %                     plot(squeeze(mean(shD(2,:,:),2))); ylim([2.4 3])
                    %%
                    testdata = squeeze(shD(1,:,:)) - squeeze(shD(2,:,:));
                    p=[];
                    for itest = 1:length(checktimes)
                        
                        [~, p(itest), ~,stat]= ttest(testdata(:,checktimes(itest)));
                        
                        tvalspertimepoint(1,itest) = stat.tstat;
                    end
                    
                    % the null hypothesis is that these prob distributions
                    % are exchangeable, so retain this permutation cluster-
                    % level stat.
                    sumTestStatsShuff(1,irand) = sum((tvalspertimepoint));
                    
                    
                end %repeat nshuff times
                
                
                %is the observed greater than CV?
                % plot histogram:
                %%
                figure(2);
                
                clf
                
%                    H=histogram(abs(sort(sumTestStatsShuff)));
                H=histogram((sort(sumTestStatsShuff)));
                % fit CDF
                cdf= cumsum(abs(H.Data))/ sum(abs(H.Data));
                %the X values (actual CV) corresponding to .01
                [~,cv05uncorr] = (min(abs(cdf-.95)));
                %                 [~,cv01uncorr] = (min(abs(cdf-.99)));
                %                 [~,cv001uncorr] = (min(abs(cdf-.999)));
                hold on
                pCV=plot([observedCV observedCV], ylim, ['r-']);
                
                p05=plot([H.Data(cv05uncorr) H.Data(cv05uncorr)], ylim, ['k:']);
                %                 plot([H.Data(cv01uncorr) H.Data(cv01uncorr)], ylim, ['k:']);
                %                 plot([H.Data(cv001uncorr) H.Data(cv001uncorr)], ylim, ['k:']);
                legend([pCV p05], {['observed'] ['95%'] })
                
                %%
               
                
                if observedCV>H.Data(cv05uncorr)
                    title(['sum tvals = ' num2str(observedCV)]);
                    %              title('Spatial Cluster  Significant!')
                    %                     timeidDYN=tgrm-3;
                  
                    
                    for itime=checktimes
                        figure(1);
                        if ~exist('yl', 'var')
                            yl=get(gca, 'ylim');
                        end
                        hold on
                        
                        
                        plot(timeidDYN(5), .6, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'k')
%                         plot(timeidDYN(itime), yl(2)+.3*(diff(yl)), ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'k')
%                         plot(timeidDYN(itime), 3.25, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'k')
                        %                             plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'm')
                    end
                end
                
                %
            end
        end
     
            %% adjust ylims.
    yl=get(gca, 'ylim');
    %extend by 10% (keeps sig points in relative space).
    dyl=diff(yl)*.1;
    ylim([yl(1)-dyl yl(2)+dyl])

    end
    
    
    

    
    %%
    cd([basefol filesep 'Figures' filesep 'GFX PFI SNR time course'])
    %
%     print('-dpng', ['PFI trace SSVEP summary, for msubtracted PFIreap.png'])
%     print('-dpng', ['PFI trace SSVEP summary, for PFIdisap-reap.png'])
end

if job.BPandSSVEPtimecourseacrossppants_group_CATCH==1


   getelocs
    
    rmvbase=0;
    checksigON=1;
    checkcluster=1;
    
    
    
    %perform subtraction = catch onset ioffset
    catchdiff=0;
    
%     clf
    plcount=1;
    legendprint=[];
    cd(basefol)
    cd('EEG')
    cd('GFX_EEG-RESS')
    %
%     clf
figure(1); hold on;
%     cd('newplots-MD')
    colsare={'b' , 'k',[],[],[],[],'m'}; % blue for tg, black for BG.
    hzare={'TG(f1)' , 'BG(f2)',[],[],[],[],'IM(f2-f1)'}; % blue for tg, black for BG.
    linetypes={'--', '-', '--','-','--'};
    legendis=[];
    lgc=1;
    for hzis=7%[1,2]%[1,2]%
        usehz=peakfreqsare(hzis);
        
        load(['GFX_Catch_performance_withSNR_' num2str(usehz) '_min0_RESS'])
        
        DIMSare{5} = 'invisible catch onset';
        %1=catch onset 2= offset, 3= BPonset 4=BP offset
        tbase = timeidDYN;
        col=colsare{hzis};
        hzp= hzare{hzis};
        ttestdata=zeros(2, length(allppants), length(Nt));
        %         legendprint=[];
        for itype=[3,4]%,4];%[1,5];%,3,4,5]%:5 6 
            
            linet=linetypes{itype};
            useBP=squeeze(storeacrossPpant_catchEVENTS_BPs(:,itype,:,:));
            useSNR=squeeze(storeacrossPpant_catchEVENTS_SNR(:,itype,:,:));
            
            
            %also take mean across trials (x100).
            %                     useBP=squeeze(mean(useBP,2));
            %                     useSNR=squeeze(mean(useSNR,2));
            
            chtype=DIMSare{itype};
            %                     chtype='button press';
            
            if catchdiff==1 %subtract offset from onset.
            
                            useSNR1 =squeeze(storeacrossPpant_catchEVENTS_SNR(:,itype,:,:));
                            useSNR2 =squeeze(storeacrossPpant_catchEVENTS_SNR(:,itype+1,:,:));
                            useSNR=useSNR1-useSNR2;
                            

                linet=':';
                
            end
                
            %%
            %
            %SNR
            figure(1);
            
            used=useSNR;
            
            hold on
            
            
            %plot across ppant trace
            ppantMeanSNR= squeeze(mean(used,2));
            
            %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            
            x = ppantMeanSNR;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
            
            stE = std(NEWdata)/sqrt(length(allppants));
            
            
            
            if rmvbase==1
                % normalise around zero.
                ppantMeanSNR = ppantMeanSNR - mean(ppantMeanSNR(:));
            end
            
            %plot
                sh=shadedErrorBar(timeidDYN, mean(ppantMeanSNR,1),stE,[],[1]);
                if rmvbase~=1
                ylabel({['RESS log(SNR), ' hzp]})
                else
                    ylabel({['mean subtracted'];['RESS log(SNR)']})
                    hold on
                    plot([xlim ], [0 0], ['k:'])
                end
            
                if catchdiff==1;
                       ylabel({['\Delta RESS log(SNR)']})
                end
                
                
            sh.mainLine.LineWidth=3;
            if itype==5
                %change color of shading, invisible catches.
                col=['r'];
            end
            sh.mainLine.Color = col;
            sh.mainLine.LineStyle=linet;
            sh.patch.FaceColor = col;
            sh.edge(1).Color = col;
            sh.edge(2).Color = col;
                      %
            %             title({[num2str(usehz) ' Hz SSVEP']})
            xlabel(['Time from ' chtype ' [s]'])
            xlabel(['Time from reporting catch [s]'])
%             xlabel(['Time from start/end catch [s]'])
%             xlabel(['Time from visible/invisible catch onset [s]'])
            %             xlabel('Time from perceptual report')
            set(gca, 'fontsize', 25)
            
            
         
            %                 ylim( [ylimsare])
            axis tight
            
            
            
            
            set(gcf, 'color', 'w')
            
            %store output for legend/stats
            ttestdata(lgc,:,:) = ppantMeanSNR;
            legendprint(lgc)=sh.mainLine;
            %legend based on type
            legendis = [legendis {chtype}];
            
            plcount=plcount+1;
            
             
            lgc=lgc+1;
            
        end % by type
        %legend based on Hz
%             legendis = [legendis {chtype}];
        
        
        axis tight
       
           xlim([timeidDYN(1) timeidDYN(end)])
    end
    % 
    if rmvbase==1
        ylim([-.6 .5])
    end
    %
%     legend(legendprint, legendis)
    %%
    % %%%%%%%% END OF PLOTTING, SIG TESTS BELOW
    if checksigON==1
        %check for sig
        pvals=zeros(1,size(ttestdata,3));
        tvals=zeros(1,size(ttestdata,3));
        for itime = 1:size(ttestdata,3)
            
            try [h,pvals(itime),~,stat]=ttest(ttestdata(1,:,itime), ttestdata(2,:,itime));
                shuffType=1;
            catch
                [h,pvals(itime),~,stat]=ttest(ttestdata(plcount-1,:,itime)); %compares to zero.
                shuffType=2; %whether or not to skip the non-parametric test for sig.
            end
            
            tvals(itime)= stat.tstat;
        end
        sigs=find(pvals<.05);
        %
        %             %perform cluster based correction.
        if length(sigs)>2 &&checkcluster==1
            % find biggest cluster:
            %finds adjacent time points
            vect1 = diff(sigs);
            v1 = (vect1(:)==1);
            d = diff(v1);
            clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
            %grab largest
            %                 ignore bad points.
            
            
            % find biggest cluster:
            %finds adjacent time points
            sigs = find(pvals<.05);
            
            vect1 = diff(sigs);
            v1 = (vect1(:)==1);
            d = diff(v1);
            clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
            [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
            
            %
            for icl=1:size(clusterSTandEND,1)
                
                %start and end are now:
                % change icl to maxClust if we only want the largest
                % cluster.
                STC=sigs(clusterSTandEND(icl,1));
                ENDC=sigs(clusterSTandEND(icl,2)+1);
                checktimes =STC:ENDC;
                observedCV = sum(abs(tvals(checktimes)));
                % now shuffle condition labels to see if this cluster is
                % sig (compared to chance).
                
                
                
                %all trials
                alltr=reshape(ttestdata, [length(allppants)*2, length(Nt)]);                
                alltrialsvec = 1:size(alltr,1);
                
                nshuff=2000;
                
                sumTestStatsShuff = zeros(1,nshuff);
                for irand = 1:nshuff
                    
                    if shuffType==1 %null is that no condition differences.
                        
                        
                        shD=zeros(size(ttestdata));
                        %since this is a within subjects design, we permute
                        %the subjet specific averages within each subject
                        %(as per Maris & Oostenveld (2007).
                        
                        % for each subject, randomly permute the averages.
                        %(Dsub1cond1,Datasub1cond2)
                        for ippant = 1:size(ttestdata,2)
                            
                            if mod(randi(100),2)==0 %if random even number
                                shD(1,ippant,:) = ttestdata(1,ippant,:); % all time points.
                                shD(2,ippant,:) = ttestdata(2,ippant,:); % all time points.
                            else
                                shD(1,ippant,:) = ttestdata(2,ippant,:); % all time points.
                                shD(2,ippant,:) = ttestdata(1,ippant,:); % all time points.
                            end
                            
%                             shD(ipartition,ippant,:) = pdata;
                        end
                        
                    else %null is that there are no temporal coincident sig values.
                        for ipartition = 1:2
                            for ippant = 1:size(ttestdata,2)
                                for itime=1:length(checktimes)
                                    
                                    %take random timepoint.
                                    pdata = ttestdata(1,ippant, randi(size(ttestdata,3)));
                                    
                                    
                                    shD(ipartition,itime,ippant) = pdata;
                                end
                            end
                        end
                    end
                    
%                     tvalspertimepoint = zeros(1,length(checktimes));
%%                    
% figure(3); clf; plot(squeeze(mean(shD(1,:,:),2))); hold on
%                     plot(squeeze(mean(shD(2,:,:),2))); ylim([2.4 3])
  %%                  
                    testdata = squeeze(shD(1,:,:)) - squeeze(shD(2,:,:));
                    p=[];
                    for itest = 1:length(checktimes)
                        
                        [~, p(itest), ~,stat]= ttest(testdata(:,checktimes(itest)));
                        
                        tvalspertimepoint(1,itest) = stat.tstat;
                    end
                    
                    % the null hypothesis is that these prob distributions
                    % are exchangeable, so retain this permutation cluster-
                    % level stat.
                    sumTestStatsShuff(1,irand) = sum((tvalspertimepoint));
                 
                        
                end %repeat nshuff times
                
                
                %is the observed greater than CV?
                % plot histogram:
                %%
                figure(2);
                
                clf
                
                
%                 H=histogram(abs(sort(sumTestStatsShuff)));
                H=histogram((sort(sumTestStatsShuff)));
                % fit CDF
                cdf= cumsum(abs(H.Data))/ sum(abs(H.Data));
                %the X values (actual CV) corresponding to .01
                [~,cv05uncorr] = (min(abs(cdf-.95)));
%                 [~,cv01uncorr] = (min(abs(cdf-.99)));
%                 [~,cv001uncorr] = (min(abs(cdf-.999)));
                hold on
                pCV=plot([observedCV observedCV], ylim, ['r-']);
                
                p05=plot([H.Data(cv05uncorr) H.Data(cv05uncorr)], ylim, ['k:']);
%                 plot([H.Data(cv01uncorr) H.Data(cv01uncorr)], ylim, ['k:']);
%                 plot([H.Data(cv001uncorr) H.Data(cv001uncorr)], ylim, ['k:']);
                legend([pCV p05], {['observed'] ['95%'] })
                
                %%
                if observedCV>H.Data(cv05uncorr)
                    title(['sum tvals = ' num2str(observedCV)]);
                    %              title('Spatial Cluster  Significant!')
%                     timeidDYN=tgrm-3;
     figure(1);
     if icl==1 % otherwise we get stepped sigheight.
yl=get(gca, 'ylim');
     end
     timeidDYN(checktimes);
                    for itime=checktimes
                   
                        hold on
                        
                        
                        plot(timeidDYN(itime), yl(1)-.1*(diff(yl)), ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'r')
                        %                             plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'm')
                    end
                end
                
                %
            end
        end
        
    end
    %% adjust ylims.
    yl=get(gca, 'ylim');
        %extend by 10% (keeps sig points in relative space).
        dyl=diff(yl)*.1;
        ylim([yl(1)-dyl yl(2)])
        
    %%
    cd(basefol)
    cd('Figures')
    cd('GFX catch SNR timecourse')
    %%
%     print('-dpng', ['Catchtrace Bground SSVEP summary, for ' hzp 'reporting catch onset.png'])
% print('-dpng', ['Catchtrace Bground SSVEP summary, for msubtracted reporting catch offset.png'])
print('-dpng', ['Catchtrace Bground SSVEP summary, for BPcatchonset - BPcatchoffset.png'])
    
    %%
    
    % print('-dpng', ['PFI trace Bground SSVEP summary, during both, 40Hz.png'])
end
  


if job.BPandSSVEPtimecourseacrossppants_group_combinePFIandCATCH==1
    getelocs

    clf
    rmvbase=0; % normalize
    checksigON=1; % check sigs between directions
    checkcluster=1; % clusters only

    collectBarwindow = [0, 1]; % seconds to collect SNR difference between disappearance/reappearance.
    
    % we have 6 combinations for Disap - Reap : ( 5/15/20 hz x PFI/Catch)
    
    collectBardata = zeros(length(peakfreqsare)*2, length(allppants));
    
    
    % plot specs:
    plcount=1;
    legendprint=[];
    cd(basefol)
    cd('EEG')
    cd('GFX_EEG-RESS')
    colsare={'b' , 'k','b','k','b','g','m'}; % blue for tg, black for BG.
    hzare={'TG (f1)' , 'BG (f2)','TG (2f1)','BG (2f2)','TG (3f1)','TGBG','IM (f2-f1)'}; % blue for tg, black for BG.
    hzare={'Target (f1)' , 'Backgroud (f2)','TG (2f1)','BG (2f2)','TG (3f1)','TGBG','IM (f2-f1)'}; % blue for tg, black for BG.
    linetypes={'--', '-', [],[],'--'};
    legendis=[];
    lgc=1;

figure(1); hold on;
clearvars yl;

lsh=[];
% per frequency, plot both the PFI data AND then CATCH data, and produce
% barchart output.
plotBARoutput=0; % plot the output or no, not need all hzis = 1:7 for BAR output.
barcounter=1;
%%
% grabmyPFIcols; % returns my cols, first dim = hz, second dim = pfi;pmd


for hzis=7%1:2%:7% need all 7 for bar data.

    %plot colour:
    col=colsare{hzis};
    
    for idtype=1:2
        
        figure(idtype); clf
        
        switch idtype
            case 1
                useD= 'PFI';
                ltype = useD; % for legend
                falpha = .25;  % alpha value to darken plot.
                
                
            case 2
                useD= 'Catch_';
                ltype = 'PMD'; % for legend
                sigcolor='r';  % color of sig points on figure,
                falpha = .25;  % alpha value to lighten plot.
              
        end
        cd(basefol)
        cd('EEG')
        cd('GFX_EEG-RESS')
        
        usehz=peakfreqsare(hzis);
        
        load(['GFX_' useD 'performance_withSNR_' num2str(usehz) '_min0_RESS'])
        
        sigcolor=col;
        
        tbase = timeidDYN;
        
        hzp= hzare{hzis};
        
        ttestdata=[];
        %         legendprint=[];
        
        
        
        
%         col =[ mycols(hzis, idtype).c];
%         col='m';
        diffTEMP = []; % we want the difference between disap/reap for each iteration.
        for itimezero=1:2%1%:2%1:2
            
            switch itimezero
                case 1
                    % different design for PFI vs Catch
                    if idtype==1
                    useSNR=storeacrossPpant_onsetSNR;
                    else
                        % BP catch onset = 
                        useSNR = squeeze(storeacrossPpant_catchEVENTS_SNR(:,3,:,:));
                    end
                    
                    chtype='button press/release';                    
                    linet=':';
                    
                    % new marker types:
                    if idtype==1
                    markert ='o';
                    else
                        markert='square';
                    end
                case 2
                    
                    if idtype==1
                    useSNR=storeacrossPpant_offsetSNR;
                    else
                        %BPcatch offset = dim 4
                    useSNR = squeeze(storeacrossPpant_catchEVENTS_SNR(:,4,:,:));
                    end
                    chtype='button press/release';
                                       
                    linet='-';
                    markert='none';
                    
                case 3
                    useSNR=storeacrossPpant_onsetSNR+storeacrossPpant_offsetSNR;
                    chtype = 'subjective report';
                    linet='.';
                    
            end
            
            figure(1);
            subplot(1,2,idtype)
            used=useSNR;
            
            hold on
            
            
            %plot across ppant trace
            ppantMeanSNR= squeeze(mean(used,2));
            %             subplot(2,1,itimezero);
            %             plot(ppantMeanSNR');
            
            
            %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            
            x = ppantMeanSNR;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
            
            stE = std(NEWdata)/sqrt(length(allppants));
            
            
            if rmvbase==1
                ppantMeanSNR=ppantMeanSNR-mean(ppantMeanSNR(:));
            end
            %plot
            sh=shadedErrorBar(timeidDYN, mean(ppantMeanSNR,1),stE,[],[1]);
            
            sh.mainLine.LineWidth=3;
            
            sh.mainLine.Color = col;
            sh.mainLine.LineStyle=linet;
            sh.patch.FaceColor = col;
            sh.patch.FaceAlpha= falpha;
            % change PMD shading.
            if idtype==2
%                 sh.mainLine.Color = 'r';
            sh.patch.FaceColor = 'r';            
            end
            
            sh.edge(1).Color = col;
            sh.edge(2).Color = col;
            
            
            sh.mainLine.Marker = markert;
            sh.mainLine.MarkerSize=25;
            
            lsh(itimezero) = sh.mainLine;
            
            xlabel(['Time from ' chtype])%
            
            
            if itimezero==3
                
                ylabel({['\Delta RESS log(SNR)']})
                
            else
                
                if rmvbase~=1
                    ylabel({['RESS log(SNR), ' hzp]})
                    ylabel({['RESS log(SNR)']})
                else
                    ylabel({['mean subtracted'];['RESS log(SNR)']})
                end
                
            end
            
       
            %store output for legend/stats
            ttestdata(itimezero,:,:) = ppantMeanSNR;
            legendprint(lgc)=sh.mainLine;
            %legend based on type
            legendis = [legendis {chtype}];
            
            
            % also collect bar data for plottting
            % over which temporal window?
            twind= dsearchn(timeidDYN', collectBarwindow');
          
            
             
            %average within ppant:
            storedata = squeeze(mean(ppantMeanSNR(:, twind), 2));
            
            diffTEMP(itimezero,:) = storedata; % we want the difference between disap/reap for each iteration.
            
            if itimezero==2; % then perform subtraction and store:                
            % we have 6 combinations (PFIx Catch, 5,15,20 hz)
            collectBardata(barcounter, :) =diffTEMP(1,:) - diffTEMP(2,:);
            barcounter=barcounter+1;
            end
            
            
            plcount=plcount+1;
            
            
            lgc=lgc+1;
            
            
        end % by type
       
        axis tight
        if rmvbase==1
            ylim([-.1 .15])
            hold on;
            plot(xlim, [0 0], ['k:'])
        end
        if itimezero==3
            plot(xlim, [0 0], ['k:'])
        end
        
      
        % %%%%%%%% END OF PLOTTING, SIG TESTS BELOW
        if checksigON==1
          temporalclustercheck
        end
        
        
         
            set(gca, 'fontsize', 25)
            
            set(gcf, 'color', 'w')
            
            yl=get(gca, 'ylim');  
            
            %Percent of ydim.
            PTy=(yl(2)-yl(1))*.35;
            ylim([yl(1)-PTy/2 yl(2)+PTy]) % f1            
            
 plot([0 0 ], ylim, ['k-']);   
 lg=legend([lsh(1), lsh(2)], [{['Disappearance (' ltype ')']}, {'Reappearance'}], 'autoupdate', 'off');
    
    set(lg, 'Location', 'Northeast'); shg
    
    end % PFI vs Catch
    
    hold on 
    

    %%
  %% adjust ylims.
       
           
    %% set same dimensions as the raincloud plots (plotted in plotSmallFigs.m)    
    set(gcf, 'units', 'normalized', 'position', [.25 .41 .25 .4], 'color', 'w');
    set(gcf, 'units', 'normalized', 'position', [.25 .41 .55 .4], 'color', 'w');
      cd(basefol)
    cd('Figures')
    cd('GFX Dynamic RESS')
    print('-dpng', ['final result at f ' num2str(hzis)]);
    
end
    
%%
if plotBARoutput ==1
    plotBARoutput_compareHz
    
end
end