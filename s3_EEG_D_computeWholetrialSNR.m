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
job.plotacrossppants=1;




if job.crunchacrossppants==1;
params.Fs=250;
params.tapers=[1 1];
% params.fpass=[0 100];
% params.fpass=[0 55];

acrossPPSNR=zeros(length(allppants), 64, 8193);
% acrossPPSNR=zeros(length(allppants), 64, 3605);
counter=1;
for ippant = 1:length(allppants)
    cd(basefol)
    cd(pdirs(allppants(ippant)).name)
    
    load(['P' num2str(allppants(ippant)) '_autopreprocd.mat']);
    
    pSNR= zeros(64,8193);
%     pSNR= zeros(64,3605);
    for ichan = 1:64
        chdata = squeeze(EEG_preprocd(ichan,:,:));
        
        [s,f]= mtspectrumc(chdata, params);
        
        kernelw = [-1/8 -1/8 -1/8 -1/8 0 0 0 0 0 1  0 0 0 0 0 -1/8 -1/8 -1/8 -1/8];
        
        stmp=zeros(size(s));
        
        for itrial = 1:size(s,2)
%             subplot(311)
%             plot(f,s(:,itrial));
%             subplot(312)
%             plot(f,log(s(:,itrial)));
%             subplot(313)
            snr=conv(log(s(:,itrial))', kernelw, 'same');
%             plot(f,snr)
            
            stmp(:,itrial) = snr;
        end
        
        %take average
        pSNR(ichan,:) = nanmean(stmp,2);
    end
    
    save('WholetrialmSNR_chanXfreq', 'pSNR', 'f')
    
    acrossPPSNR(counter,:,:)=pSNR;
end
cd(basefol)
save('GFX_rawSNR_static_wholetrial', 'acrossPPSNR', 'f')

end
if job.plotacrossppants==1
    %%
    cd(basefol)
    load('GFX_rawSNR_static_wholetrial');
    cd ../
    cd('Figures')
    cd('SNR spectrum')
    %%
    figure(1);
    clf
    mSNRchan=squeeze(mean(acrossPPSNR(:,62,:),1));
    plot(f, mSNRchan)
    hold on
    cols={'r', 'k', 'r', 'k', 'm','m', 'r', 'r', 'k','m'};
    counter=1;
    p=[];
        for ifreq=[15,20,30,40,5,35, 45,60,80, 25]
            [~,usef]= min(abs(f-ifreq));
            yis= mSNRchan(usef);
            xis= f(usef);
           p(counter)=plot(xis, yis, 'o', 'linew', 2, 'color', cols{counter}) ;
           
           if ifreq==60
               p(counter)=plot(xis, yis, 'o', 'markersize', 15,'linew', 2, 'color', 'k') ;
           end
            counter=counter+1;
            
        end
        legend([p(1), p(2), p(5)], {'TG Hz', 'BG Hz', 'im'})
        axis tight
        xlabel('Frequency (Hz)')
        ylabel('logSNR')
        set(gcf, 'color', 'w')
        set(gca, 'fontsize', 25)
        ylim([0 .4])
        %%
        print('-dpng', 'Whole trial SNR, chan POz')
        %%
        cd(basefol)
        cd ../
        cd('Figures')
    cd('Topoplots')
    %%
    %topos first.
    getelocs;
    counter=1;
    mTOPOs = squeeze(mean(acrossPPSNR,1));
    clf
    for ifreq=[15,30,20,40,5,25,35,60]
        
        [~,usef]= min(abs(f-ifreq));
        subplot(2,4,counter)
        
        topoplot(squeeze(mTOPOs(:,usef)), elocs(1:64))
        c=colorbar;
        caxis([0 .2])
        title([num2str(ifreq) ' Hz'])
        set(gca, 'fontsize', 25)
        counter=counter+1;
    end
    %%
    set(gcf, 'color', 'w')
    print('-dpng', 'Topos across all, SNR preRESS, whole trial')
    %%
    
end
        
        
        