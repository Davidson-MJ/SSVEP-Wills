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

job.crunchacrossppants=0; %saves also within participant folders.
job.plotacrossppants=0;


job.crunchacrossppants_compareWindow=1;

if job.crunchacrossppants==1
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
        



        
if job.crunchacrossppants_compareWindow==1;
    %% also plot change in SNR with time window considered.
%space time kernels.
lt = [1,2,4,6,10,12,15,20,30,60]; % for averaging.        
params.Fs=250;
params.tapers=[1 1];
params.pad=2;
params.fpass=[0 55];
% params.fpass=[0 100];
% params.fpass=[0 55];

hztest = [5,15,20, 30, 40];
acrossPPSNR=zeros(length(allppants), length(lt), length(hztest));
acrossPPSNR_spectrum= [];
acrossPPsignalpower= acrossPPSNR;
acrossPPnoisepower = acrossPPSNR;
counter=1;

for ippant = 1:length(allppants)
    cd(basefol)
    cd(pdirs(allppants(ippant)).name)
    
    load(['P' num2str(allppants(ippant)) '_autopreprocd.mat']);

    ichan = 62; % POz
    
    
    chdata = squeeze(EEG_preprocd(ichan,:,:));
    
    %for each step size, use the specgram function.
    counter=1;
   figure(1); clf
   
        for ilt = 1:length(lt)
       %%
       
            movingwin= [lt(ilt), lt(ilt)];
        
        [spg,tgrm,fgrm,]= mtspecgramc(chdata, movingwin, params);
        
        % take mean over segments.
        if ~ismatrix(spg)
            mspg= squeeze(mean(spg,1));
        else
            mspg=spg;
        end
            
        %
        %adjust kernel shape appropriately.
         k = 1; %tapers
         hbw = (k+1)./ (2.*[movingwin(1,1)]);
         neighb = 1; %hz
         
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
%         kernelw = [-1/8 -1/8 -1/8 -1/8 0 0 0 0 0 1  0 0 0 0 0 -1/8 -1/8 -1/8 -1/8];
        % take SNR>
        stmp=zeros(size(mspg));
        for itrial = 1:size(mspg,2)
            snr=conv(log(mspg(:,itrial))', kernelw, 'same');

            %sanity check
%             plot(fgrm,log(mspg(:,itrial))); hold on
%             plot(fgrm,snr)
            
            stmp(:,itrial) = snr;
        end
        
        
        %take average at correct freqs, this timestep.
        
        hzid= dsearchn(fgrm', [hztest]');
    
    
%     save('WholetrialmSNR_chanXfreq', 'pSNR', 'f')
    
    acrossPPSNR(ippant,ilt,:)=squeeze(mean(stmp([hzid],:),2));
    
    % also retain power at signal
    acrossPPsignalpower(ippant,ilt,:)=squeeze(mean(mspg([hzid],:),2));
    
    %and noise. note that for noise calculations, the location (index) is
    %computed above.
    for itmp= 1:length(hzid)
    centref= (hzid(itmp));
    nbelowvec = (centref- (distz(1)+distz(2))) : (centref- distz(1));
    nabovevec = (centref+ distz(1) : (centref+ distz(1)+distz(2)));
    
    noiseatfreq = squeeze(mean(mspg([nbelowvec,nabovevec], :) ,2));
    acrossPPnoisepower(ippant,ilt,itmp)= squeeze(mean(noiseatfreq));
    end
    
    acrossPPSNR_spectrum(ilt).data(ippant,:) = squeeze(mean(stmp,2));
    
    
%     if ippant == length(allppants) % plot across all snr.
    figure(1);
    subplot(5,2, counter)
    tmp=acrossPPSNR_spectrum(ilt).data;
    plot(fgrm, squeeze(mean(tmp,1)), 'k');
    title([num2str(lt(ilt)) 's, hbw = ' sprintf('%.2f',hbw) ' Hz']);
ylim([-2 5]);     hold on;
    % center kernel at 20 Hz
%     hlfkern= (length(kernelw)-1)/2;
%     fvector = (hzid(3)- hlfkern):(hzid(3)+hlfkern);
%     plot(fgrm(fvector), kernelw*2, 'r')
    xlim([1 40])
    counter=counter+1;
    xlabel('Frequency (Hz)')
    ylabel('log(SNR)')
    set(gca, 'fontsize', 15)
%     end
    end

%%
set(gcf, 'color', 'w')
cd(basefol)
cd ../
cd(['Figures'  filesep 'SNR by different temporal windows'])
%%
print('-dpng', [' SNR spectrum by window length, ppant ' num2str(ippant)]);

%%
figure(2); clf
for  iplotype=1:3 % signal, noise, snr.
subplot(3,1,iplotype)
switch iplotype
    case 1
        plotd= log(acrossPPsignalpower(1,:,:));
        yis= 'log signal power';
    case 2
        plotd= log(acrossPPnoisepower(1,:,:));
        yis= 'log noise power';
    case 3
        plotd= acrossPPSNR(1,:,:);
        yis='log(SNR)';
end

macrossPPSNR = squeeze(mean(plotd,1));
plot(lt,macrossPPSNR);
ylabel(yis); xlabel('window length [s]')
legend('IM', 'TG f1', 'BG f2', 'TG 2f1', 'BG 2f2')
set(gca, 'xtick', lt, 'fontsize', 15); set(gcf, 'color', 'w')
xlim([1 10]);
%% pvals

for ihz=1%:5;
    pvals=[];
for it = 1:size(plotd,2)
    
[~,pvals(it)]= ttest(squeeze(plotd(:,it,ihz)));
end
end
%%

if iplotype==3
%     hold on;
%     plot(xlim, [0 0 ], ['k:'])
%     
%     % also plot at 2 and 4, 
%     plot([2 2], [-1 mean(macrossPPSNR(2,1))], ['b-'])
%     plot([4 4], [-1 mean(macrossPPSNR(3,1))], ['b-'])
end
end
%%

suptitle(['Participant ' num2str(ippant)])
print('-dpng', [' window length summary, ppant ' num2str(ippant)]);

% save('GFX_rawSNR_static_wholetrial', 'acrossPPSNR', 'f')
end
end      
