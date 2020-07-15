
%CALLED from s3_EEG_alphala_evoked


job.performClusterstat_infieldtrip=0;

% once the above is done, progress:
%% >>>>>>
job.plotHILB_final =1; % paper figure, comparing PFI to PMD>
job.plotAlpha_f1corr =0;


usePFIorPMD=1; % for the above.


%% run field trip cluster correction, to identify significant electrodes.

if job.performClusterstat_infieldtrip==1
    try load('Ft-clusterresults.mat')
    catch
        clusterAlpha_statsinFT;
    end
end
 

 
 %%
 
 % LARGE posterior electrode subset.
if  job.plotHILB_final ==1; % paper figure, comparing PFI to PMD>
    
    %laod both PFI and PMD>
     %first, trim all, 
     %%
     cd(basefol)
     cd('EEG')
     for iPFIvsPMD=1:2
     
         if iPFIvsPMD==1
             load('GFX_LeftvsRight_TFdecomp.mat')
             load('Ft-clusterresults.mat', 'clustCHANS');
             %load channels for plots.
             clustCHANS_PFI = clustCHANS;
         else
             load('GFX_LeftvsRight_TFdecomp_PMD.mat')
             load('Ft-clusterresultsPMD.mat', 'clustCHANS');
             clustCHANS_PMD = clustCHANS;
         end
      tstamps = [1:1:size(all_leftON_hilb,3)]/250 - 3;

     %same var names
     trim = dsearchn(tstamps', [-2.5 2.5]');
     timeax = tstamps(trim(1):trim(2));
   GFX_leftPFI = all_leftON_hilb(:,:,trim(1):trim(2));
   GFX_rightPFI = all_rightON_hilb(:,:,trim(1):trim(2));
   GFX_leftReap = all_leftOFF_hilb(:,:,trim(1):trim(2));
   GFX_rightReap = all_rightOFF_hilb(:,:,trim(1):trim(2));
   
   
     %data should be chansxframesxn (for comparison.
     if iPFIvsPMD==1
         
     allPFI = cat(4, GFX_leftPFI, GFX_rightPFI);
     allPFI= squeeze(mean(allPFI,4));
     allReap = cat(4, GFX_leftReap, GFX_rightReap);
     allReap= squeeze(mean(allReap,4));
     else
         allPFI_PMD = cat(4, GFX_leftPFI, GFX_rightPFI);
        allPFI_PMD= squeeze(mean(allPFI_PMD,4));
        allReap_PMD = cat(4, GFX_leftReap, GFX_rightReap);
     allReap_PMD= squeeze(mean(allReap_PMD,4));
     end
     end
     %%
     % now that we have the data, plot both on the same axes:
     clf
     plot_PFIvsPMD_Alpha;
         
         
    
end
%%
if job.plotAlpha_f1corr ==1;

    %% for both, lets compare the same window. Cluster start end in seconds, 
    % is :
    
    %load the RESS f1 data
    
    clf
    for iPFIvsPMD=1%1:2
        cd(basefol)
        cd('EEG')
        cd('GFX_EEG-RESS')        
        %%
        if iPFIvsPMD==1
        load(['GFX_PFIperformance_withSNR_15_min0_RESS'])
         
        cd(basefol)
        cd('EEG')        
            load('Ft-clusterresults.mat')
            load('GFX_LeftvsRight_TFdecomp.mat')
            coln='b';
        else
            
            %load PMD f1 and convert to same size.
            load('GFX_Catch_performance_withSNR_15_min0_RESS.mat')
            %third dimension is the BP locked catch onset data.
               
        cd(basefol)
        cd('EEG')   
            storeacrossPpant_onsetSNR = squeeze(storeacrossPpant_catchEVENTS_SNR(:,3,:,:));
            load('Ft-clusterresultsPMD.mat')            
            load('GFX_LeftvsRight_TFdecomp_PMD.mat')
            coln='r';
        end
        
        twin= [];
        twin(1) = staterp.time(find(plength,1, 'first'));
        twin(2) = staterp.time(find(plength,1, 'last'));
        
%         twin = [-1  1.5]
        %extract the SNR over the cluster window:
        avgg= dsearchn(timeidDYN', twin');
        RESSsnr = mean(squeeze(nanmean(storeacrossPpant_onsetSNR(:, :, avgg(1):avgg(2)),3)),2);
        RESSsnr_t = squeeze(nanmean(storeacrossPpant_onsetSNR,2));
        figure(3);
        hold on
        plot(squeeze(mean(RESSsnr_t)))
        %now do the same for alpha, per subject.
        %%
        % load alpha
        cd(basefol)
        cd('EEG')
        
        
        tstamps = [1:1:size(all_leftON_hilb,3)]/250 - 3;
        
        %same var names for alpha in both cases.
        trim = dsearchn(tstamps', [-2.5 2.5]');
        timeax = tstamps(trim(1):trim(2));
        GFX_leftPFI = all_leftON_hilb(:,:,trim(1):trim(2));
        GFX_rightPFI = all_rightON_hilb(:,:,trim(1):trim(2));
        
        
        allPFI = cat(4, GFX_leftPFI, GFX_rightPFI);
        allPFI= squeeze(mean(allPFI,4));
        
        %convert to log scale to avoid non-normal distribution.
        
        
        
        %perform baseline normalization:
        tmp = zeros(size(allPFI));
        %divide by baseline per subject.
        baselineis = dsearchn(timeax', [-2.5 -2]');
        
        for isub = 1:size(allPFI,1)
            for ichan = 1:64
                tmpchan = squeeze(allPFI(isub, ichan, :));
                tmpbase = nanmean(tmpchan(baselineis(1):baselineis(2)));
                tmp(isub,ichan,:)= (tmpchan./ tmpbase) ;
            end
        end
        allPFI=tmp; % this is the alpha data.
        
        
        
        clustCHANSn = abs(clustCHANS-1);
        clustwindow = find(plength);

        Alphacorr= mean(squeeze(nanmean(allPFI(:,find(clustCHANSn),clustwindow),2)),2);
        
        figure(1); 
%         subplot(1,2, iPFIvsPMD);  
        hold on;
        set(gcf, 'units', 'normalized', 'position', [0 1 .29 .29], 'color', 'w')
        
        %% before scatter, check for outliers:
        rma = isoutlier(Alphacorr);
        rmp = isoutlier(RESSsnr);
        rmall = [find(rma), find(rmp)];
        
        outlierplots = [RESSsnr(rmall), Alphacorr(rmall)];
        if ~isempty(rmall)
            Alphacorr(rmall) = [];
            RESSsnr(rmall) = [];
        end
        
%         %log transform?
%         RESSsnr = log(RESSsnr);
%         Alphacorr= log(Alphacorr);
        %% plot remaining data points.
        H=scatter( (RESSsnr),(Alphacorr), 100, coln); shg
        H.Marker = '*';
        H.LineWidth = 3; 
        hold on 
        %now plot outliers.
         H=scatter(outlierplots(:,1), outlierplots(:,2), 100, coln); shg
        H.Marker = 'o';
        H.LineWidth = 3; 
        %%
        
        [r,p] = corrcoef([Alphacorr, RESSsnr]);
%         [rho, pval] = corr(Alphacorr, RESSsnr, 'type', 'Spearman');
%               title(['rho = ' sprintf('%.2f', rho) ', p=' sprintf('%.2f', pval)])
        P= polyfit(RESSsnr, Alphacorr,1);
        
        Y_fit= polyval(P, RESSsnr);
        
        hold on
        lg= plot(RESSsnr, Y_fit, 'linew', 4, 'color', coln);
        xlabel({['f1 RESS log(SNR)']})
        ylabel({['Evoked Alpha during PFI']});
        set(gca, 'fontsize', 20)

        axis square
if iPFIvsPMD ==1         
    ylim([ 0.75 1.05])
else
    ylim([ 0.55 1.25])
end
        
        legh=legend(lg, ['r = ' sprintf('%.2f', r(1,2)) ', p=' sprintf('%.3f', p(1,2))]);
legh=legend(lg, ['\itr\rm = .58, \itp\rm = .019']);
        legh.FontSize = 20;
        legh.Location = 'SouthWest';
        shg
    end
end
     
   