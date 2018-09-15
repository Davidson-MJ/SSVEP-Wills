
% WILLS DATA
addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')


basefol=pwd;
clearvars -except basefol allppants
dbstop if error

cd('EEG');
pdirs = dir([pwd filesep '*EEG']);


%%



%%% some new analyses (at top)
% calc and store Whole trial RESS SNR, for correlation with behaviour.
job.wholetrialRESS=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try a cross validation procedure, within trials, epoch any change
%first we predict the SNR for PFI, using subset of PFI trials.
%training using PFI.

getelocs;
allppants=[1,2,4,6,9:16,18]; %


if job.wholetrialRESS==1

    
    
    
    
   
    param_spcgrm.tapers = [1 1];
    param_spcgrm.Fs= [250];
    param_spcgrm.fpass= [0 50];
    param_spcgrm.trialave=0;
    movingwin=[1,.15];
   %%
    cd(basefol)
    cd('Behaviour')
    %% load data to determine physical catch timing.
    load('MD_AllBP_perppant')
    load('PFI_data')
    %%
    
    srate=250;
    
    
    
    peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ];
    %remaining participants after behavioral data analysis/exclusion
  
    
    icounter=1;
    for ifol = allppants%(12)
        cd(basefol)
        cd('EEG')
        cd(pdirs(ifol).name)
        %load new rereferenced data.
        load(['P' num2str(ifol) '_autopreprocd'])
        
      
        
        snrOUT_TS=[];
        
        snr_wholetrRESS=zeros(length(peakfreqsare),48,389); % hz, trials, samps)
        
        
         load('RESSfilterdata.mat');
        %apply RESS to whole trial:
        
        for ifreq=1:length(peakfreqsare)
           usefreq=peakfreqsare(ifreq);
                         
            
            datast= EEG_preprocd;
          
            for itrial = 1:size(datast,3)
            
                
                datapre= squeeze(datast(:,:,itrial));

            
            % need to account for TG hz location.
            evecs = squeeze(ressEVEC_byHz(ifreq,:));
            
            %reduce data single time series.
            
            %apply filter 
            
                ress_ts1 = evecs*(datapre);
            
            
        % 
          
          
          
        [sgrm, tgrm_wholetrial, fgrm_wholetrial]= mtspecgramc(ress_ts1, movingwin, param_spcgrm);
        
        %now take SNR>
        kernelw = [-.25 -.25 0 0 1 0 0 -.25 -.25];


         snr_sgrm =zeros(size(sgrm));
                
                
         %compute SNR
         
             tmps=sgrm;
             
             for itime= 1:size(tmps,1)
                 checkput = conv(log(tmps(itime,:)), kernelw,'same');
                 if ~isreal(checkput)
                     snr_sgrm(itime,:)= nan(1, size(tmps,2));
                 else
                     
                     %% Note the +1
                     
                     snr_sgrm(itime,:)= conv(log(tmps(itime,:)), kernelw,'same');
%                      clf;
%                      plot(fgrm_wholetrial, snr_sgrm(itime,:));
%                      shg
                     
                 end
             end
         
                
%                 
%                 
%         %reduce to two freqs.
%         
%         
        [~, idF]= min(abs(fgrm_wholetrial-usefreq));
%         [~, id40]= min(abs(fgrm-40));
%         
        snrOUT_TS= squeeze(snr_sgrm(:,[idF]))';
        
          

snr_wholetrRESS(ifreq, itrial,:) = snrOUT_TS;
            end
        end
        
        
        
        
        save('P_wholetrialRESS', 'snr_wholetrRESS','fgrm_wholetrial', 'tgrm_wholetrial')
        
        
        disp(['fin ppant' num2str(ifol)])
        icounter=icounter+1;
    end
    
end
%should have all targets now.


%plot whole trial spectrum?



