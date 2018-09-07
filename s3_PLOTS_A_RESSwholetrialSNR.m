
clear all
addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/EEG')
basefol=pwd;
    %%
    cd(basefol)    
    cd('GFX_EEG-RESS')
    %
    load('GFX_ressSNR_static_wholetrial');
    cd ../../
    cd('Figures')
    %%
%     cd('SNR spectrum')
    %%
    figure(1);
    clf
        counter=1;

    for ifreq=1:size(acrossPPSNR,2);
    mSNRchan=squeeze(mean(acrossPPSNR(:,ifreq,:),1));
    
    plot(f, mSNRchan,'k')
    p=[];
    hold on
%     peakfreqsare= 15,20,30,40,45,60,5,25,35]
    cols={'r', 'k', 'r', 'k', 'r','k', 'm', 'm', 'm'};
            
    [~,usef]= min(abs(f-peakfreqsare(ifreq)));
            
    yis= mSNRchan(usef);
    
    xis= f(usef);
    
    p(counter)=plot(xis, yis, 'o', 'linew', 2, 'color', cols{counter}) ;
    
           if peakfreqsare(ifreq)==60
               p(counter)=plot(xis, yis, 'o', 'markersize', 15,'linew', 2, 'color', 'k') ;
           end
            counter=counter+1;
            
        
    end
    %%
    ylim([-1 6]); xlim([0 85])
    %%
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