% s3_B_rereferenceEEG
%new rereference script.
clear all
addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/EEG')
basefol=pwd;
%%
pdirs = dir([pwd filesep '*_*' 'EEG']);

% allppants=[1,2,4,6,9:16,18]; %
 allppants=[1,2,4,6,7,9:19]; %

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
    counter=counter+1;
end
cd(basefol)
save('GFX_rawSNR_static_wholetrial', 'acrossPPSNR', 'f')

end
if job.plotacrossppants==1
    %%
    cd(basefol)
    load('GFX_rawSNR_static_wholetrial');
    %%
    cd ../
    cd('Figures')
    cd('SNR spectrum')
    %%
    figure(1);
    clf
    mSNRchan=squeeze(mean(acrossPPSNR(:,62,:),1));
    plot(f, mSNRchan, 'k')
    hold on
    cols={'b', 'k', 'b', 'k', 'm','m', 'b', 'b', 'k','m','b'};
    counter=1;
    p=[];
        for ifreq=[15,20,30,40,5,35, 45,60,80, 25, 75]
            [~,usef]= min(abs(f-ifreq));
            yis= mSNRchan(usef);
            xis= f(usef);
           p(counter)=plot(xis, yis, 'o', 'markersize', 15, 'linew', 5, 'color', cols{counter}) ;
%            p(counter)=text(xis, yis, 'o', 'fontsize', 25, 'color', cols{counter}) ;
           
           if ifreq==60
               p(counter)=plot(xis, yis, 'o', 'markersize', 30,'linew', 5, 'color', 'k') ;
           end
            counter=counter+1;
            
        end
        %%
       lg= legend([p(1), p(2), p(5)], {'Target flicker', 'Background flicker', 'Intermodulation'});
       set(lg, 'fontsize', 35)
       
        axis tight
        xlabel('Frequency (Hz)')
        ylabel('log(SNR)')
        set(gcf, 'color', 'w')
        set(gca, 'fontsize', 35)
%         ylim([0 .4])
xlim([ 0 85])
ylim([-1 6])
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
    %%
    counter=1;
    colormap('viridis')
    mTOPOs = squeeze(mean(acrossPPSNR,1));
    clf
     titles={'(TG) f1', '(TG) 2f1', '(TG) 3f1',...
        '(BG) f2', '(BG) 2f2', '(BG) 3f2',...
        '(IM) f2-f1', '(IM) 2f2-f1', '(IM) f1 +f2'};
    for ifreq=[15,30,45,20,40,60,5,25,35]
        
        [~,usef]= min(abs(f-ifreq));
        subplot(3,3,counter)
        
        plotme=squeeze(mTOPOs(:,usef));
        topoplot(squeeze(mTOPOs(:,usef)), elocs(1:64))
        c=colorbar;
        cm=round(max(plotme));
        if cm<1
            cm=1;
        end
%         caxis([0 cm])
        caxis([0 3])
%         title([num2str(ifreq) ' Hz'])
title(titles{counter})
        set(gca, 'fontsize', 25)
        counter=counter+1;
       ylabel(c,'log(SNR)')
    end
    colormap('viridis')
    %
    shg
    %%
    set(gcf, 'color', 'w')
    print('-dpng', 'Topos across all, SNR preRESS, whole trial')
    %%
    % also alternate plot
    counter=1
     for ifreq=[5,15,20,25,30,35,40,45,60]
        
        [~,usef]= min(abs(f-ifreq));
        subplot(1,9,counter)
        
        plotme=squeeze(mTOPOs(:,usef));
        topoplot(squeeze(mTOPOs(:,usef)), elocs(1:64))
      
        
        cm=round(max(plotme));
        if cm<1
            cm=1;
        end

        caxis([0 2])
%         title([num2str(ifreq) ' Hz'])
% title(titles{counter})
title([num2str(ifreq)])
        set(gca, 'fontsize', 25)
        counter=counter+1;
%           if ifreq==60
%         c=colorbar;
        
%        ylabel(c,'log(SNR)')
%           end
    end
    colormap('viridis')
    %
    shg
    %%
    set(gcf, 'color', 'w')
    %%
    print('-dpng', 'Topos across all, SNR preRESS, whole trial, singlefile caxis 2')
end
        
        
        