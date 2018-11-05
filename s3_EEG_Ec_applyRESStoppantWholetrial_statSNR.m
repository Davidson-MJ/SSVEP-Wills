% s3_B_rereferenceEEG
%new rereference script.
clear all
addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/EEG')
basefol=pwd;
%%
pdirs = dir([pwd filesep '*_*' 'EEG']);

allppants=[1,2,4,6,9:16,18]; %

job.crunchacrossppants=1; %saves also within participant folders.




if job.crunchacrossppants==1;
params.Fs=250;
params.tapers=[1 1];
% params.fpass=[0 100];
% params.fpass=[0 55];

acrossPPSNR=zeros(length(allppants), 9, 8193);
% acrossPPSNR=zeros(length(allppants), 64, 3605);
counter=1;
for ippant = 1:length(allppants)
    cd(basefol)
    cd(pdirs(allppants(ippant)).name)
    
    load(['P' num2str(allppants(ippant)) '_autopreprocd.mat']);
    load('RESSfilterdata')
    
    %for each frequency, apply ress to collapse across channels.
    datast=EEG_preprocd;
    
    
    
    for ifreq=1:length(peakfreqsare)

        
        evecs = squeeze(ressEVEC_byHz(ifreq,:)); %
        
        ress_ts1=zeros(size(datast,3), size(datast,2));
        %apply RESS per rtial,
        for ti=1:size(datast,3)
            ress_ts1(ti,:) = evecs*(squeeze(datast(:,:,ti)));
        
            %then take SNR
            chdata=ress_ts1(ti,:);            
        
            [s,f]= mtspectrumc(chdata, params);
        
            kernelw = [-1/8 -1/8 -1/8 -1/8 0 0 0 0 0 1  0 0 0 0 0 -1/8 -1/8 -1/8 -1/8];
        
            snr=conv(log(s), kernelw, 'same');
%             plot(f,snr)
            
            stmp(:,ti) = snr;
        end
        
        %now store the SNR
        
        %take average
        pSNR(ifreq,:) = nanmean(stmp,2);
    end
    
    save('WholetrialmSNR_RESSxfreq', 'pSNR', 'f')
    
    acrossPPSNR(counter,:,:)=pSNR;
    counter=counter+1;
end
cd(basefol)
cd('GFX_EEG-RESS')
save('GFX_ressSNR_static_wholetrial', 'acrossPPSNR', 'f')

end

        
        
        