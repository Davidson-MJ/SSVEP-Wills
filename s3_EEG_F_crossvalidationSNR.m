
% WILLS DATA
addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')


basefol=pwd;
clearvars -except basefol allppants
dbstop if error

cd('EEG');
pdirs = dir([pwd filesep '*EEG']);



%
job.processTrainingData_singleparticipant=0;
%predicted based on PFI.
job.processPredictedData_singleparticipant=0; 
%can also predict the catch from  remaining trials (untrained), to show
%difference.
%predicted Catch based on PFI
job.processPredictedandObservedCATCHData_singleparticipant=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%also flip this process, and use the catch data (SNR) as training, 
%to attempt to predict PFI. Useful for plots to show dissociation between
%subjective and physical SNR time-course. 
%  *not completed*

job.processTrainingData_singleparticipant_basedonCatch=0;
job.processPredictedData_singleparticipant_basedonCatch=0; 
job.processPredObsPFIData_singleparticipant_basedonCatch=0;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%concat results of above.
job.concatPredandObsacrossppants=0;
    %%%% plotting
%plot results of above.
job.plotPredvsObserved_acrossparticipants=0;
job.plotPredvsObserved_acrossparticipants_ver2=1;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Second tier of jobs. this time to check classification accuracy of
%background harmonic SNR.

job.processharmonicClassificationacc_perppant=0;
job.processharmonicClassificationacc_Predicted_perppant=0;
job.concatClassificationacc=0;
job.plotClassificationaccuracyacrossppants=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 




peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; 

if job.processTrainingData_singleparticipant==1
    
    cd(basefol)
    %first load BP data to extract disap and reapp.
    cd('Behaviour')
    
    load('MD_AllBP_perppant')
    load('PFI_data')
    
    
    
    for  ifol=allppants %this is ppant of interest.
        
        for ihz=1:length(peakfreqsare)
        
            useHz=peakfreqsare(ihz);
            
                    
        cd(basefol)
        cd('EEG')
        cd(pdirs(ifol).name)
        
        load('P_wholetrialRESS'); %load snr data
        
        
        
        window=[-3 3];
        
        srate=1/min(unique(diff(tgrm_wholetrial)));
        epochdur = sum(abs(window))*srate;
        onsetc = ceil(epochdur)/2;
        
        Setlist=[];
        TrainingDatatmp=[];
                
        for iSet=1:10

        
        
        %Random selection of 2/3 cross validation, each set list.
        randp= randi(48,[ 1 36]);
        while unique(randp)<48
            randp=randi(48,[1 36]);
        end
        
        Setlist(iSet).training_tr=randp;
        %this script captures the real SNR, so perform for both the
            %training and remaining data.
            for iobs=1:2
                if iobs==1
                    %first time, collect ERP for training data
                    trials = Setlist(iSet).training_tr;
                else
                    %second time, collect ERPfor remaining trials
                    tmp=1:48;
                    tmp(Setlist(iSet).training_tr)=[];
                    trials=tmp;
                end
                
                
                
              %separate by number involved
        %disappearances:
        ppant_SNREEG_PFI_0_1 = [];
        counter0_1=1;
        ppant_SNREEG_PFI_1_2 = [];
        counter1_2=1;
        ppant_SNREEG_PFI_2_3 = [];
        counter2_3=1;
        %Reappearances
        ppant_SNREEG_PFI_3_2 = [];
        counter3_2=1;
        ppant_SNREEG_PFI_2_1 = [];
        counter2_1=1;
        ppant_SNREEG_PFI_1_0 = [];
        counter1_0=1;
        
        ppant_SNREEG_PFI_3_4 = [];
        ppant_SNREEG_PFI_4_3 = [];
        counter3_4=1;
        counter4_3=1;
        
        durs0_1=0;
        durs1_0=0;
        durs1_2=0;
        durs2_1=0;
        durs2_3=0;
        durs3_2=0;
        durs3_4=0;
        durs4_3=0;
        
        BPs0_1=[];
        BPs1_0=[];
        BPs1_2=[];
        BPs2_1=[];
        BPs2_3=[];
        BPs3_2=[];
        BPs3_4=[];
        BPs4_3=[];
        
                
                Locationwas=[]; %same for location of PFI events BPs.
                
              ehz=ihz;
                
                for itrial=trials
                    
                    
                    SNRdata = squeeze(snr_wholetrRESS(ehz,itrial,:))';
                    
                    trialdata= PFI_only(ifol).Trial(itrial);
                    if trialdata.Goodtrial==1
                        %first establish the direction of changes (disap /reappear
                        %targets).
                        accumBP = nansum(trialdata.allBPs,1);
                        directionBP = diff(accumBP);
                        directionBPtimes = find(directionBP~=0);
                        directionBPtimes=directionBPtimes+1; %account for diff function
                        %now remove zeros, leaving time of all disap/reap
                        directionBP(directionBP==0)=[];
                        
                        
                        
                        
                        %we want to only keep the BPs that survived our previous,
                        %exclusion criteria.
                        %so create vector of all the relevant BP frames already stored.
                        survivingPFI = [trialdata.PFI_disap_0target_framestart,...
                            trialdata.PFI_disap_1target_framestart,...
                            trialdata.PFI_disap_2target_framestart,...
                            trialdata.PFI_disap_4target_framestart];
                        
                        
                        
                        %so now only continue by assessing those BPs (ie remove
                        %transients)
                        keepdirectionBPtimes = find(ismember(directionBPtimes, survivingPFI));
                        directionBPtimes= directionBPtimes(keepdirectionBPtimes);
                        
                        
                        
                        for iPFI = 1:length(directionBPtimes)
                            %Collect event info.
                            timeis = directionBPtimes(iPFI);
                            perceptis = accumBP(directionBPtimes(iPFI));
                            perceptwas = accumBP(directionBPtimes(iPFI)-1);
                            
                            
                            outgoingPFI=[];
                            %disap or reap.
                            if perceptwas<perceptis
                                PFIdir=1; %disappearing, an extra button is being pressed
                            else
                                PFIdir=-1; %Reappearing
                            end
                            
                            
                            %gather EEG time index.
                            [~,timePFI]= min(abs(tgrm_wholetrial-timeis/60));
                            
                            
                            
                            
                            % ignore events at start/end of block.
                            
                            if  timePFI-onsetc>1 && timePFI+onsetc<size(SNRdata,2)
                                %collect data
                                outgoingPFI= SNRdata(:,timePFI-onsetc:timePFI+onsetc);
                                
                                
                                
                                %now store this PFI epoch according to type.
                                switch perceptis %num pressed
                                    case 0 %currently NO buttons pressed.
                                        if PFIdir==1
                                            error('checkcode') %cannot go from less than zero button press.
                                        else
                                            if perceptwas==1 % gone from 1 to 0 BP (reappearing)
                                                try
                                                    ppant_SNREEG_PFI_1_0(counter1_0,:,:) = outgoingPFI;
                                                    
                                                    %how long was this event?
                                                    %find the correct framestamp and duration.
                                                    [~,indx]=min(abs(trialdata.PFI_disap_1target_framestart - timeis));
                                                    durs=trialdata.PFI_disap_1target_durs;
                                                    
                                                    %save duration of this event for averaging
                                                    %'best' later.
                                                    %                                         durs1_0(counter1_0) = durs(indx);
                                                    
                                                    
                                                    
                                                    %Store which frequency was involved:
                                                    % timestap across BPs.
                                                    tmpBPs=trialdata.allBPs(:,timeis-1); % prev BP since reappearing.
                                                    
                                                    
                                                    if ~isnan(sum(tmpBPs))
                                                        
                                                        
                                                        Locationwas.dir1_0(counter1_0).d= (find(tmpBPs));
                                                    else %it was during catch period, so don't count.
                                                        
                                                        
                                                        Locationwas.dir1_0(counter1_0).d= nan;
                                                    end
                                                    
                                                    % store BP also
                                                    
                                                    BPs1_0(counter1_0,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                                    
                                                    
                                                    
                                                    counter1_0=counter1_0+1;
                                                catch
                                                    disp('skipped a 1_0')
                                                end
                                            elseif perceptwas==2
                                                try
                                                    ppant_SNREEG_PFI_2_1(counter2_1,:,:) = outgoingPFI;
                                                    
                                                    [~,indx]=min(abs(trialdata.PFI_disap_2target_framestart - timeis));
                                                    %                                         durs=trialdata.PFI_disap_2target_durs;
                                                    
                                                    %save duration of this event for averaging
                                                    %'best' later.
                                                    %                                         durs2_1(counter2_1) = durs(indx);
                                                    
                                                    
                                                    
                                                    %Store which frequency was involved:
                                                    % timestap across BPs.
                                                    tmpBPs1=trialdata.allBPs(:,timeis-1);
                                                    tmpBPs2=trialdata.allBPs(:,timeis);
                                                    dtmp= tmpBPs1-tmpBPs2;
                                                    
                                                    
                                                    
                                                    if ~isnan(sum(dtmp))
                                                        
                                                        Locationwas.dir2_1(counter2_1).d= (find(tmpBPs));
                                                    else %it was during catch period, so don't count.
                                                        
                                                        Locationwas.dir2_1(counter2_1).d= (find(tmpBPs));
                                                    end
                                                    
                                                    
                                                    
                                                    %also store BP trace.
                                                    BPs2_1(counter2_1,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                                    
                                                    
                                                    counter2_1=counter2_1+1;
                                                catch
                                                    disp('skipped a 2_1')
                                                end
                                            end
                                        end
                                        
                                    case 1 %single target
                                        if PFIdir==1 % PFI numbers increasing
                                            try
                                                ppant_SNREEG_PFI_0_1(counter0_1,:,:) = outgoingPFI;
                                                
                                                %how long was this event?
                                                %find the correct framestamp and duration.
                                                [~,indx]=min(abs(trialdata.PFI_disap_1target_framestart - timeis));
                                                %                                     durs=trialdata.PFI_disap_1target_durs;
                                                
                                                %save duration of this event for averaging
                                                %'best' later.
                                                %                                     durs0_1(counter0_1) = durs(indx);
                                                
                                                
                                                
                                                
                                                
                                                %also store BP trace.
                                                BPs0_1(counter0_1,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                                
                                                
                                                counter0_1=counter0_1+1;
                                            catch
                                                disp('skipped a 0_1')
                                            end
                                        else %decreasing numbers being pressed
                                            if perceptwas==2
                                                try
                                                    ppant_SNREEG_PFI_2_1(counter2_1,:,:) = outgoingPFI;
                                                    
                                                    [~,indx]=min(abs(trialdata.PFI_disap_2target_framestart - timeis));
                                                    %                                         durs=trialdata.PFI_disap_2target_durs;
                                                    %
                                                    %                                         %save duration of this event for averaging
                                                    %                                         %'best' later.
                                                    %                                         durs2_1(counter2_1) = durs(indx);
                                                    %
                                                    
                                                    %also store BP trace.
                                                    BPs2_1(counter2_1,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                                    
                                                    counter2_1=counter2_1+1;
                                                catch
                                                    disp('skipped a 2_1')
                                                end
                                            elseif perceptwas==3
                                                try
                                                    ppant_SNREEG_PFI_3_2(counter3_2,:,:) = outgoingPFI;
                                                    
                                                    
                                                    [~,indx]=min(abs(trialdata.PFI_disap_3target_framestart - timeis));
                                                    %                                         durs=trialdata.PFI_disap_3target_durs;
                                                    %
                                                    %                                         %save duration of this event for averaging
                                                    %                                         %'best' later.
                                                    %                                         durs3_2(counter3_2) = durs(indx);
                                                    
                                                    %Store which frequency was involved:
                                                    % timestap across BPs.
                                                    tmpBPs1=trialdata.allBPs(:,timeis-1);
                                                    tmpBPs2=trialdata.allBPs(:,timeis);
                                                    dtmp=tmpBPs1-tmpBPs2;
                                                    
                                                    
                                                    
                                                    
                                                    if ~isnan(sum(dtmp))
                                                        
                                                        Locationwas.dir3_2(counter3_2).d= (find(tmpBPs));
                                                    else %it was during catch period, so don't count.
                                                        
                                                        Locationwas.dir3_2(counter3_2).d= nan;
                                                    end
                                                    
                                                    
                                                    
                                                    %also store BP trace.
                                                    BPs3_2(counter3_2,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                                    
                                                    
                                                    
                                                    counter3_2=counter3_2+1;
                                                catch
                                                    disp('skipped a 3_2')
                                                end
                                            end
                                            
                                        end
                                    case 2 %double target
                                        if PFIdir==1
                                            try
                                                ppant_SNREEG_PFI_1_2(counter1_2,:,:) = outgoingPFI;
                                                [~,indx]=min(abs(trialdata.PFI_disap_2target_framestart - timeis));
                                                %                                     durs=trialdata.PFI_disap_2target_durs;
                                                %
                                                %                                     %save duration of this event for averaging
                                                %                                     %'best' later.
                                                %                                     durs1_2(counter1_2) = durs(indx);
                                                %
                                                
                                                %Store which frequency was involved:
                                                % timestap across BPs.
                                                tmpBPs1=trialdata.allBPs(:,timeis-1);
                                                tmpBPs2=trialdata.allBPs(:,timeis);
                                                dtmp=tmpBPs2-tmpBPs1;
                                                
                                                if ~isnan(sum(dtmp))
                                                    
                                                    Locationwas.dir1_2(counter1_2).d= (find(dtmp));
                                                else %it was during catch period, so don't count.
                                                    
                                                    Locationwas.dir1_2(counter1_2).d= nan;
                                                end
                                                
                                                
                                                
                                                %also store BP trace.
                                                BPs1_2(counter1_2,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                                
                                                
                                                counter1_2=counter1_2+1;
                                            catch
                                                disp('skipped a 1_2')
                                            end
                                        else
                                            try
                                                ppant_SNREEG_PFI_3_2(counter3_2,:,:) = outgoingPFI;
                                                
                                                
                                                [~,indx]=min(abs(trialdata.PFI_disap_3target_framestart - timeis));
                                                %                                     durs=trialdata.PFI_disap_3target_durs;
                                                
                                                %save duration of this event for averaging
                                                %'best' later.
                                                %                                     durs3_2(counter3_2) = durs(indx);
                                                
                                                %Store which frequency was involved:
                                                % timestamp across BPs.
                                                tmpBPs1=trialdata.allBPs(:,timeis-1);
                                                tmpBPs2=trialdata.allBPs(:,timeis);
                                                dBPs= tmpBPs1-tmpBPs2;
                                                
                                                
                                                if ~isnan(sum(dBPs))
                                                    
                                                    Locationwas.dir3_2(counter3_2).d= (find(dBPs));
                                                else %it was during catch period, so don't count.
                                                    
                                                    Locationwas.dir3_2(counter3_2).d= nan;
                                                end
                                                
                                                
                                                %also store BP trace.
                                                BPs3_2(counter3_2,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                                
                                                counter3_2=counter3_2+1;
                                            catch
                                                disp('skipped a 3_2')
                                            end
                                        end
                                        
                                    case 3 %triple target
                                        if PFIdir==1
                                            try
                                                ppant_SNREEG_PFI_2_3(counter2_3,:,:) = outgoingPFI;
                                                
                                                [~,indx]=min(abs(trialdata.PFI_disap_3target_framestart - timeis));
                                                %                                     durs=trialdata.PFI_disap_3target_durs;
                                                %
                                                %                                     %save duration of this event for averaging
                                                %                                     %'best' later.
                                                %                                     durs2_3(counter2_3) = durs(indx);
                                                %
                                                
                                                
                                                %Store which frequency was involved:
                                                % timestap across BPs.
                                                tmpBPs1=trialdata.allBPs(:,timeis-1);
                                                tmpBPs2=trialdata.allBPs(:,timeis);
                                                dtmp = tmpBPs2-tmpBPs1;
                                                
                                                if ~isnan(sum(dtmp))
                                                    
                                                    Locationwas.dir2_3(counter2_3).d= (find(dtmp));
                                                else %it was during catch period, so don't count.
                                                    Locationwas.dir2_3(counter2_3).d= nan;
                                                end
                                                
                                                %also store BP trace.
                                                BPs2_3(counter2_3,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                                
                                                Freqwas.dir2_3(counter2_3).Trialind= itrial;
                                                
                                                counter2_3=counter2_3+1;
                                            catch
                                                disp('skipped a 2_3')
                                            end
                                        else %decreasing.
                                            try
                                                ppant_SNREEG_PFI_4_3(counter4_3,:,:) = outgoingPFI;
                                                
                                                [~,indx]=min(abs(trialdata.PFI_disap_4target_framestart - timeis));
                                                %                                     durs=trialdata.PFI_disap_4target_durs;
                                                %
                                                %                                     %save duration of this event for averaging
                                                %                                     %'best' later.
                                                %                                     durs4_3(counter4_3) = durs(indx);
                                                %
                                                
                                                
                                                %Store which frequency was involved:
                                                % timestap across BPs.
                                                tmpBPs1=trialdata.allBPs(:,timeis-1);
                                                tmpBPs2=trialdata.allBPs(:,timeis);
                                                dtmp = tmpBPs2-tmpBPs1;
                                                
                                                if ~isnan(sum(dtmp))
                                                    
                                                    Locationwas.dir4_3(counter4_3).d= (find(dtmp));
                                                else %it was during catch period, so don't count.
                                                    
                                                    Locationwas.dir4_3(counter4_3).d= nan;
                                                end
                                                
                                                %also store BP trace.
                                                BPs4_3(counter4_3,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                                
                                                
                                                
                                                counter4_3=counter4_3+1;
                                            catch
                                                disp('skipped 4_3')
                                            end
                                        end
                                    case 4
                                        if PFIdir==1
                                            try
                                                ppant_SNREEG_PFI_3_4(counter3_4,:,:) = outgoingPFI;
                                                
                                                [~,indx]=min(abs(trialdata.PFI_disap_4target_framestart - timeis));
                                                %                                     durs=trialdata.PFI_disap_4target_durs;
                                                %
                                                %                                     %save duration of this event for averaging
                                                %                                     %'best' later.
                                                %                                     durs3_4(counter3_4) = durs(indx);
                                                %
                                                
                                                
                                                %Store which frequency was involved:
                                                % timestap across BPs.
                                                tmpBPs1=trialdata.allBPs(:,timeis-1);
                                                tmpBPs2=trialdata.allBPs(:,timeis);
                                                dtmp = tmpBPs2-tmpBPs1;
                                                
                                                if ~isnan(sum(dtmp))
                                                    
                                                    Locationwas.dir3_4(counter3_4).d= (find(dtmp));
                                                else %it was during catch period, so don't count.
                                                    Locationwas.dir3_4(counter3_4).d= nan;
                                                end
                                                
                                                %also store BP trace.
                                                BPs3_4(counter3_4,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                                
                                                counter3_4=counter3_4+1;
                                            catch
                                                disp('skipped a 3_4')
                                            end
                                        else %decreasing.
                                            try
                                                ppant_SNREEG_PFI_4_3(counter4_3,:,:) = outgoingPFI;
                                                
                                                [~,indx]=min(abs(trialdata.PFI_disap_4target_framestart - timeis));
                                                %                                     durs=trialdata.PFI_disap_4target_durs;
                                                %
                                                %                                     %save duration of this event for averaging
                                                %                                     %'best' later.
                                                %                                     durs4_3(counter4_3) = durs(indx);
                                                %
                                                
                                                
                                                %Store which frequency was involved:
                                                % timestap across BPs.
                                                tmpBPs1=trialdata.allBPs(:,timeis-1);
                                                tmpBPs2=trialdata.allBPs(:,timeis);
                                                dtmp = tmpBPs2-tmpBPs1;
                                                
                                                if ~isnan(sum(dtmp))
                                                    
                                                    Locationwas.dir4_3(counter4_3).d= (find(dtmp));
                                                else %it was during catch period, so don't count.
                                                    
                                                    Locationwas.dir4_3(counter4_3).d= nan;
                                                end
                                                
                                                %also store BP trace.
                                                BPs4_3(counter4_3,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                                
                                                
                                                
                                                counter4_3=counter4_3+1;
                                                %
                                            catch
                                                disp('skipped a 4_3')
                                            end
                                        end
                                        
                                        
                                end
                            end
                        end
                        
                        
                        
                        
                        
                    end
                    
                    
                    
                end
                
                
%                         %%
%                         keeptrials=[];
%                         for item = 1:size(tmp,2) %this is needed, as some disaps had simultaneous BP, messes with indexing.
%                             hzc=tmp(item).d;
%                             if any(hzc==useHz)
%                         keeptrials=[keeptrials, (item)];
%                             end
%                         end
                        %%
%                       
                      
                
                ppantOnsetSNR = squeeze(cat(1, ppant_SNREEG_PFI_0_1,ppant_SNREEG_PFI_1_2,ppant_SNREEG_PFI_2_3,ppant_SNREEG_PFI_3_4));
                
                
                ppantOffsetSNR= squeeze(cat(1, ppant_SNREEG_PFI_1_0,ppant_SNREEG_PFI_2_1,ppant_SNREEG_PFI_3_2,ppant_SNREEG_PFI_4_3));
                
                
                
                
                
                %store this average.
                switch iobs
                    case 1
                        TrainingDatatmp(iSet).Training_ppantOnsetSNR = ppantOnsetSNR;
                        TrainingDatatmp(iSet).Training_ppantOffsetSNR = ppantOffsetSNR;
                    case 2
                        %store the observed data for the remaining trials, to
                        %compare next with 'predicted'.
                        TrainingDatatmp(iSet).Observed_ppantOnsetSNR = ppantOnsetSNR;
                        TrainingDatatmp(iSet).Observed_ppantOffsetSNR = ppantOffsetSNR;
                end
            end
            xt_plot= linspace(window(1), window(2), size(ppantOnsetSNR,2));
            %
        end
        
        %save results.
        savename= ['PFI_SNR_trainingData_allHz'];
        

% peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; 
        
        switch ihz % order as per peakfreqsare
            case 1
                TrainingData_15Hz = TrainingDatatmp;
            case 2
                TrainingData_20Hz = TrainingDatatmp;
            case 3
                TrainingData_30Hz = TrainingDatatmp;
            case 4
                TrainingData_40Hz = TrainingDatatmp;
            case 5
                TrainingData_45Hz = TrainingDatatmp;            
            case 6                
                TrainingData_60Hz = TrainingDatatmp;
            case 7
                TrainingData_5Hz = TrainingDatatmp;
            case 8
                TrainingData_25Hz = TrainingDatatmp;
            case 9
                TrainingData_35Hz = TrainingDatatmp;
            
                
        end
        
        end
    
    save(savename,  'TrainingData_15Hz',...
'TrainingData_20Hz',...
'TrainingData_30Hz',...
'TrainingData_40Hz',...
'TrainingData_45Hz',...
'TrainingData_60Hz',...
'TrainingData_5Hz',...
'TrainingData_25Hz',...
'TrainingData_35Hz',...
'Setlist', 'xt_plot');

    end
end

%now create predicted (from remaining trials)
if job.processPredictedData_singleparticipant
    %Now we artificially model each disap and reapp.
    cd(basefol)
    %first load BP data to extract disap and reapp.
    cd('Behaviour')
    load('MD_AllBP_perppant')
    load('PFI_data')
    %%
    %the next is modelled on job 1 in s3_C_EpochPFIEEG:
    
    for  ifol=allppants %this is ppant of interest.
        cd(basefol)
        cd('EEG')        
cd(pdirs(ifol).name)
load('PFI_SNR_trainingData_allHz.mat')
        
        load('P_wholetrialRESS.mat', 'tgrm_wholetrial');
        
        window=[-3 3];
        
        srate=1/min(unique(diff(tgrm_wholetrial)));
        epochdur = sum(abs(window))*srate;
        onsetc = ceil(epochdur)/2;
        
        for ihz=1:9
           switch ihz % order as per peakfreqsare
            case 1
                TrainingDatatmp=  TrainingData_15Hz ;
            case 2
                TrainingDatatmp=  TrainingData_20Hz ;
            case 3
                TrainingDatatmp=  TrainingData_30Hz ;
            case 4
                TrainingDatatmp=  TrainingData_40Hz ;
            case 5
                TrainingDatatmp=  TrainingData_45Hz ;            
            case 6                
                TrainingDatatmp=  TrainingData_60Hz ;
            case 7
                TrainingDatatmp=  TrainingData_5Hz ;
            case 8
              TrainingDatatmp=    TrainingData_25Hz ;
            case 9
              TrainingDatatmp=  TrainingData_35Hz;
            
           end
           
                    
        
        predictedSSVEP=[]; %for storage.
        
        for iSet=1:size(Setlist,2)
            
            %load the training onset and offset ERPs
            onsetERP=squeeze(nanmean(TrainingDatatmp(iSet).Training_ppantOnsetSNR,1));
            offsetERP=squeeze(nanmean(TrainingDatatmp(iSet).Training_ppantOffsetSNR,1));
            
            %need ERP mid point to align with BP data.
            midpERP= dsearchn(xt_plot',[0]);
            
            %create vector for our model
            mvector = zeros(size(tgrm_wholetrial));
            
            %adjust to baseline SNR level of observed data
            mvector= mvector + [mean(onsetERP(:))+mean(offsetERP(:))]/2;
            
            
            
            %find which trials were not used in this list:
            allt=[1:48];
            %rem trained trials.
            allt(Setlist(iSet).training_tr)=[];
            itrialcount=1;
            
            
            
              outgoingPFIons=[];
                    outgoingPFIoffs=[];
                    inc=1;
                    outc=1;
                    
            for itrial = allt
                %load trial data.
                trialdata= PFI_only(ifol).Trial(itrial);
                if trialdata.Goodtrial==1
                    %first establish the direction of changes (disap /reappear
                    %targets).
                    accumBP = nansum(trialdata.allBPs,1);
                    directionBP = diff(accumBP);
                    directionBPtimes = find(directionBP~=0);
                    directionBPtimes=directionBPtimes+1; %account for diff function
                    %now remove zeros, leaving time of all disap/reap
                    directionBP(directionBP==0)=[];
                    
                    
                    
                    
                    %we want to only keep the BPs that survived our previous,
                    %exclusion criteria.
                    %so create vector of all the relevant BP frames already stored.
                    survivingPFI = [trialdata.PFI_disap_0target_framestart,...
                        trialdata.PFI_disap_1target_framestart,...
                        trialdata.PFI_disap_2target_framestart,...
                        trialdata.PFI_disap_3target_framestart,...
                        trialdata.PFI_disap_4target_framestart];
                    
                    
                    
                    
                    %so now only continue by assessing those BPs (ie remove
                    %transients)
                    keepdirectionBPtimes = find(ismember(directionBPtimes, survivingPFI));
                    directionBPtimes= directionBPtimes(keepdirectionBPtimes);
                    
                    %for each pfi. input a
                    for iPFI = 1:length(directionBPtimes)
                        %Collect event info.
                        timeis = directionBPtimes(iPFI);
                        perceptis = accumBP(directionBPtimes(iPFI));
                        perceptwas = accumBP(directionBPtimes(iPFI)-1);
                        
                        
                        outgoingPFI=[];
                        %disap or reap.
                        if perceptwas<perceptis
                            PFIdir=1; %disappearing, an extra button is being pressed
                            useERP=onsetERP;
                        else
                            PFIdir=-1; %Reappearing
                            useERP=offsetERP;
                        end
                        
                        %zero the ERP
                        useERP=useERP-mean(useERP(:));
                        
                        %gather EEG time index.
                        
                        [~,timePFI]= min(abs(tgrm_wholetrial-timeis/60));
                        
                        % this is the time of change (center)
                        
                        
                        % ignore events at start/end of block.
                        
                        if  timePFI-onsetc>1 && timePFI+onsetc<length(tgrm_wholetrial)
                            
                            %input onset/offset ERP at each point.
                            %first rescale
                            
                            
                            
                            try
                                mvector(timePFI-midpERP:timePFI+(length(onsetERP)-midpERP)-1) = (mvector(timePFI-midpERP:timePFI+(length(onsetERP)-midpERP)-1) + useERP);
                                %adjust remaining timeline. 
                                lp = find(mvector,1,'last');
                                mvector(lp:end)=mvector(lp-1);
                            catch
                            end
                        end
                        
                        
                        
                        
                        
                    end
                    
                    %at the end of this trial, store the predicted SSVEP.
                    predictedSSVEP(iSet,itrialcount,:) = mvector;
                    
                    %now what are the average ERPs for this data?
                    
                  
                     for iPFI = 1:length(directionBPtimes)
                        %Collect event info.
                        timeis = directionBPtimes(iPFI);
                        perceptis = accumBP(directionBPtimes(iPFI));
                        perceptwas = accumBP(directionBPtimes(iPFI)-1);
                        
                        %gather EEG time index.                        
                        [~,timePFI]= min(abs(tgrm_wholetrial-timeis/60));
                        
                        
                        if  timePFI-onsetc>1 && timePFI+onsetc<size(mvector,2)
                            
                            %disap or reap.
                            if perceptwas<perceptis
                                PFIdir=1; %disappearing, an extra button is being pressed
                                outgoingPFIons(inc,:)= mvector(:,timePFI-onsetc:timePFI+onsetc);
                                inc=inc+1;
                            else
                                PFIdir=-1; %Reappearing
                                outgoingPFIoffs(outc,:)= mvector(:,timePFI-onsetc:timePFI+onsetc);
                                outc=outc+1;
                            end
                        end
                        
                        
                     end
                    
                  
                    itrialcount=itrialcount+1;
                end
            end
            
                     
            %now that we have new 'predicted' SSVEP time-course, what are
            %the predicted averages ERPs?
          
            
            %add to training Data.
            TrainingDatatmp(iSet).predicted_ppantonsetSNR= outgoingPFIons;
            TrainingDatatmp(iSet).predicted_ppantoffsetSNR= outgoingPFIoffs;
            
        end
        
        %rename before saving
        
        switch ihz % order as per peakfreqsare
            case 1
                TrainingData_15Hz = TrainingDatatmp;
            case 2
                TrainingData_20Hz = TrainingDatatmp;
            case 3
                TrainingData_30Hz = TrainingDatatmp;
            case 4
                TrainingData_40Hz = TrainingDatatmp;
            case 5
                TrainingData_45Hz = TrainingDatatmp;            
            case 6                
                TrainingData_60Hz = TrainingDatatmp;
            case 7
                TrainingData_5Hz = TrainingDatatmp;
            case 8
                TrainingData_25Hz = TrainingDatatmp;
            case 9
                TrainingData_35Hz = TrainingDatatmp;
            
                
        end
        
        end
        
           savename= ['PFI_SNR_trainingData_allHz'];
         save(savename,  'TrainingData_15Hz',...
'TrainingData_20Hz',...
'TrainingData_30Hz',...
'TrainingData_40Hz',...
'TrainingData_45Hz',...
'TrainingData_60Hz',...
'TrainingData_5Hz',...
'TrainingData_25Hz',...
'TrainingData_35Hz',...
'Setlist', 'xt_plot');
    
    end
end

if job.processPredictedandObservedCATCHData_singleparticipant
    %Now we artificially model each disap and reapp.
    cd(basefol)
    cd('Behaviour')
    %first load BP data to extract disap and reapp.
    
    load('MD_AllBP_perppant')
    load('PFI_data')
    %%
    

    
    for  ifol=allppants %this is ppant of interest.
    cd(basefol)
    cd('EEG')    
        cd(pdirs(ifol).name)
        load('PFI_SNR_trainingData_allHz.mat') %to concatenate date.
        load('P_wholetrialRESS.mat'); %snr for epoching.
        load('TrialIndicesbyCatchLocationandNum.mat'); %for catch type.
        %         load('Trial
        window=[-3 3];
        
        srate=1/min(unique(diff(tgrm_wholetrial)));
        epochdur = sum(abs(window))*srate;
        onsetc = ceil(epochdur)/2;
        
        for ihz=1:9
       
           switch ihz % order as per peakfreqsare
            case 1
                TrainingDatatmp=  TrainingData_15Hz ;
            case 2
                TrainingDatatmp=  TrainingData_20Hz ;
            case 3
                TrainingDatatmp=  TrainingData_30Hz ;
            case 4
                TrainingDatatmp=  TrainingData_40Hz ;
            case 5
                TrainingDatatmp=  TrainingData_45Hz ;            
            case 6                
                TrainingDatatmp=  TrainingData_60Hz ;
            case 7
                TrainingDatatmp=  TrainingData_5Hz ;
            case 8
              TrainingDatatmp=    TrainingData_25Hz ;
            case 9
              TrainingDatatmp=  TrainingData_35Hz;
            
           end
            
            %
   
        
        predictedSSVEP=[]; %for storage.
        
        for iSet=1:size(Setlist,2)
            
            %since no catch events were used in generating predicted data,
            %we can use ALL catch events, per freq, in this run.
            %find which trials were not used in this list:
            
            %per freq, we need to collect the trials that had a catch
            %disappearance per freq.
               allt=[1:48];
            %rem trained trials.
            allt(Setlist(iSet).training_tr)=[];
            
            nonSettrials=allt;
            
            
%             if mod(ihz,2)~=0
%                 %then not BG
% %                usehz=
%                 
%                 %find the approriate trials.
%                 allt=squeeze([TopLeftCatchindexbyHz(usehz,:),TopRightCatchindexbyHz(usehz,:),...
%                     BottomLeftCatchindexbyHz(usehz,:),BottomRightCatchindexbyHz(usehz,:)]);
%                 allt(allt<1)=[];
%             else
%                 allt=1:24; %use all trials
%             end
%             
allt=1:48;            
            
            %which trials can we use in this set list, given the Hz
            %restriction?
            
            usetrials = ismember(allt, nonSettrials);
            allt= allt(usetrials);
            
            
            itrialcount=1;
            
            
            
              outgoingPFIons=[];
                    outgoingPFIoffs=[];
                    outgoing_realons=[];
                    outgoing_realoffs=[];
                    
                    inc=1;
                    outc=1;
                    
                    %ok, so for non-trained trials, calculate predicted and
                    %observed catch SNR (
            for itrial = allt 
                
                realSNR = squeeze(snr_wholetrRESS(ihz, itrial,:))';
                
                trialdata= ppantTrialDatawithDetails(ifol).TrialDetails(itrial);
                trialGB = PFI_only(ifol).Trial(itrial); %need for 'good trial' sub field.
                
                
                
                if trialGB.Goodtrial==1 && sum(trialdata.CatchBPs_exact(:))>0
                    % % % prepare predicted SNR.
                    %load the training onset and offset ERPs
                    onsetERP=squeeze(nanmean(TrainingDatatmp(iSet).Training_ppantOnsetSNR,1));
                    offsetERP=squeeze(nanmean(TrainingDatatmp(iSet).Training_ppantOffsetSNR,1));
                    
                    %need ERP mid point to align with BP data.
                    midpERP= dsearchn(xt_plot',[0]);
                    
                    %create vector for our model
                    mvector = zeros(size(tgrm_wholetrial));
                    
                    %adjust to baseline SNR level of observed data
                    mvector= mvector + [mean(onsetERP(:))+mean(offsetERP(:))]/2;
                    
                       
                    % % % determine BP timing of catch responses 
                        
                        catchrampstart = trialdata.Catchstart_frames/60;
                        catchrampend =  trialdata.Catchend_frames/60; %perhaps remove the ramp?
            
                        CatchBPs = nansum(trialdata.CatchBPs_exact,1);
                        firstBPwithincatch= min(find(diff(CatchBPs)>0));
                        BP_disapwithincatch_secs = firstBPwithincatch/60 + catchrampstart;

                        %also collect and epoch around first BP after catch
                        % return.
                    
                    %note that the 'long epoch' is 2s before catch start and
                    %8seconds after catch START
                    
                    catchdur=trialdata.totalCatchdur;
                    
                    longcatchEpoch=nansum(trialdata.CatchBPs_longepoch,1);
                    
                    %catch end frames, within this epoch = diff +2s for the
                    %onset
                    postOffset=longcatchEpoch(1,(catchdur+120):end);
                    
                    %first change in BP, (could also use return to zero?)
                    
                    firstBPaftercatch= min(find(diff(postOffset)<0));
                    catchreturn_secs = trialdata.Catchend_frames/60 + firstBPaftercatch/60;
                    
                    %now collect physical onset/offet, as well as BP
                    %onset/offset (per trial).
                    
                    for iPFI = 3:4
                        %Collect event info.
                        
                        switch iPFI
%                             case 1
%                                 centrep=catchrampstart;
%                             case 2
%                                 centrep=catchrampend;
                            case 3
                                 %if BPdata during catch removal
                                 centrep=BP_disapwithincatch_secs;
                            case 4
                                centrep=catchreturn_secs;
                        end
                        
                        %we want an onset or offset ERP
                        if mod(iPFI,2)==0 %even numbers (target return)
                            useERP=offsetERP;
                        else
                            useERP=onsetERP;
                        end
                        
                        %zero the ERP
                        useERP=useERP-mean(useERP(:));
                        
                        %gather EEG time index.
                        try [~,timePFI]= min(abs(tgrm_wholetrial-centrep));
                        catch
                            continue
                        end
                        % this is the time of change (center)
                        
                        
                        % ignore events at start/end of block.
                        
                        if  timePFI-onsetc>1 && timePFI+onsetc<length(tgrm_wholetrial)
                            
                            %input onset/offset ERP at each point.
                            %first rescale
                            
                            
                            
                            try
                                mvector(timePFI-midpERP:timePFI+(length(onsetERP)-midpERP)-1) = (mvector(timePFI-midpERP:timePFI+(length(onsetERP)-midpERP)-1) + useERP);
                                %adjust remaining timeline. 
                                lp = find(mvector,1,'last');
                                mvector(lp:end)=mvector(lp-1);
                            catch
                            end
                        end
                        
                            
                    end
                    %%
                    
                    
                    %at the end of this trial, store the predicted SSVEP.
%                     predictedSSVEP(iSet,itrialcount,:) = mvector;
                    
                    %now what are the average ERPs for this data?                  
                     for iPFI = 3:4
                        switch iPFI
%                             case 1
%                                 centrep=catchrampstart;
%                             case 2
%                                 centrep=catchrampend;
                            case 3
                                 %if BPdata during catch removal
                                 centrep=BP_disapwithincatch_secs;
                            case 4
                                centrep=catchreturn_secs;
                        end
                         
                        %gather EEG time index.  
                        try
                        [~,timePFI]= min(abs(tgrm_wholetrial-centrep));
                        catch %if no BP to epoch around, skip this trial.
%                             
                                
                            timePFI=0;
                        end
                        
                        if  timePFI-onsetc>1 && timePFI+onsetc<size(mvector,2)
                            
                            %disap or reap.
                            switch iPFI
                                case 3
                                    
                                    outgoingPFIons(inc,:)= mvector(:,timePFI-onsetc:timePFI+onsetc);
                                    outgoing_realons(inc,:)=realSNR(:,timePFI-onsetc:timePFI+onsetc);
                                    inc=inc+1;
                                    
                                case 4
                                    outgoingPFIoffs(outc,:)= mvector(:,timePFI-onsetc:timePFI+onsetc);
                                    outgoing_realoffs(outc,:)=realSNR(:,timePFI-onsetc:timePFI+onsetc);
                                    outc=outc+1;
                            end
                        
                          
                        else %store a nan.
                            if iPFI==3
                                outgoingPFIons(inc,:)= nan(1,2*onsetc+1);
                                outgoing_realons(inc,:)=nan(1,2*onsetc +1);
                                inc=inc+1;
                            elseif iPFI==4
                                
                                outgoingPFIoffs(outc,:)= nan(1,2*onsetc+1);
                                outgoing_realoffs(outc,:)=nan(1,2*onsetc +1);
                                outc=outc+1;
                            end
                            
                        
                        
                        end
                     end
                  
                    itrialcount=itrialcount+1;
                end
            end
            %its possible all trials were skipped (if none available):
            if isempty(outgoingPFIons)
                outgoingPFIons=nan(1,length(onsetERP));
            end
            if isempty(outgoingPFIoffs)
                outgoingPFIoffs=nan(1,length(onsetERP));
            end
            if isempty(outgoing_realons)
                outgoing_realons=nan(1,length(onsetERP));
            end
            if isempty(outgoing_realoffs)
                outgoing_realoffs=nan(1,length(onsetERP));
            end
                     
            %now that we have new 'predicted' SSVEP time-course, what are
            %the predicted averages ERPs?
          
            if isempty(outgoingPFIons)
                break
            end
            %add to training Data.
            TrainingDatatmp(iSet).predicted_BPcatchonsetSNR= outgoingPFIons;
            TrainingDatatmp(iSet).predicted_BPcatchoffsetSNR= outgoingPFIoffs;
            TrainingDatatmp(iSet).Observed_BPcatchonsetSNR= outgoing_realons;
            TrainingDatatmp(iSet).Observed_BPcatchoffsetSNR= outgoing_realoffs;
            
        end
        
        %rename before saving
          switch ihz % order as per peakfreqsare
            case 1
                TrainingData_15Hz = TrainingDatatmp;
            case 2
                TrainingData_20Hz = TrainingDatatmp;
            case 3
                TrainingData_30Hz = TrainingDatatmp;
            case 4
                TrainingData_40Hz = TrainingDatatmp;
            case 5
                TrainingData_45Hz = TrainingDatatmp;            
            case 6                
                TrainingData_60Hz = TrainingDatatmp;
            case 7
                TrainingData_5Hz = TrainingDatatmp;
            case 8
                TrainingData_25Hz = TrainingDatatmp;
            case 9
                TrainingData_35Hz = TrainingDatatmp;
            
                
        end
        
        end
        
           savename= ['PFI_SNR_trainingData_allHz'];
         save(savename,  'TrainingData_15Hz',...
'TrainingData_20Hz',...
'TrainingData_30Hz',...
'TrainingData_40Hz',...
'TrainingData_45Hz',...
'TrainingData_60Hz',...
'TrainingData_5Hz',...
'TrainingData_25Hz',...
'TrainingData_35Hz',...
'Setlist', 'xt_plot');
    
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% now for catch based analyses.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









if job.concatPredandObsacrossppants==1
    
  %%
  
            acrossallPP_TrainingData=[];
            for ihz=1:9
                
                ifolc=1; %initiate counter.
                for ifol=allppants
                    
                    cd(basefol)
                    cd('EEG')
                    cd(pdirs(ifol).name)
                    load('PFI_SNR_trainingData_allHz.mat');
                    
                    
                    switch ihz % order as per peakfreqsare
                        case 1
                            TrainingDatatmp=  TrainingData_15Hz ;
                        case 2
                            TrainingDatatmp=  TrainingData_20Hz ;
                        case 3
                            TrainingDatatmp=  TrainingData_30Hz ;
                        case 4
                            TrainingDatatmp=  TrainingData_40Hz ;
                        case 5
                            TrainingDatatmp=  TrainingData_45Hz ;
                        case 6
                            TrainingDatatmp=  TrainingData_60Hz ;
                        case 7
                            TrainingDatatmp=  TrainingData_5Hz ;
                        case 8
                            TrainingDatatmp=    TrainingData_25Hz ;
                        case 9
                            TrainingDatatmp=  TrainingData_35Hz;
                            
                    end
                    
                    
                    for iset=1:size(Setlist,2)
                        %take ppant mean and store.
                        acrossallPP_TrainingData.predicted_onsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).predicted_ppantonsetSNR,1));
                        acrossallPP_TrainingData.predicted_offsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).predicted_ppantoffsetSNR,1));
                        acrossallPP_TrainingData.observed_onsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).Observed_ppantOnsetSNR,1));
                        acrossallPP_TrainingData.observed_offsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).Observed_ppantOffsetSNR,1));
                        %also add catch data
                        acrossallPP_TrainingData.predicted_BPcatchonsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).predicted_BPcatchonsetSNR,1));
                        acrossallPP_TrainingData.predicted_BPcatchoffsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).predicted_BPcatchoffsetSNR,1));
                        acrossallPP_TrainingData.observed_BPcatchonsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).Observed_BPcatchonsetSNR,1));
                        acrossallPP_TrainingData.observed_BPcatchoffsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).Observed_BPcatchoffsetSNR,1));
                    end
                    ifolc=ifolc+1; %folder counter.
                end
                    
                


% peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; 
                    switch ihz
                        case 1
                            allppTrainingData_15Hz=acrossallPP_TrainingData;
                        case 2
                            allppTrainingData_20Hz=acrossallPP_TrainingData;
                        case 3
                            allppTrainingData_30Hz=acrossallPP_TrainingData;
                        case 4
                            allppTrainingData_40Hz=acrossallPP_TrainingData;
                        case 5
                            allppTrainingData_45Hz=acrossallPP_TrainingData;
                        case 6

                            allppTrainingData_60Hz=acrossallPP_TrainingData;
                        case 7
                            allppTrainingData_5Hz=acrossallPP_TrainingData;
                        case 8
                            allppTrainingData_25Hz=acrossallPP_TrainingData;
                        case 9
                            allppTrainingData_35Hz=acrossallPP_TrainingData;
                        
                    end
                    
                
            end
    cd(basefol)
    cd('EEG')
    cd('GFX_EEG-RESS')
    save('GFX_predictedvsobserved', 'allppTrainingData_15Hz',...
        'allppTrainingData_20Hz',...
'allppTrainingData_30Hz',...
'allppTrainingData_40Hz',...
'allppTrainingData_45Hz',...
'allppTrainingData_60Hz',...
'allppTrainingData_5Hz',...
'allppTrainingData_25Hz',...
'allppTrainingData_35Hz',...
'Setlist', 'xt_plot')
    
end

if job.plotPredvsObserved_acrossparticipants==1
    %%
    cd(basefol)
    cd('EEG')
    cd('GFX_EEG-RESS')
    
    %%
  load('GFX_predictedvsobserved')
 
  %% which HZ to look at?



% peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; 

  pcount=1;
%   clf
  usePFIorCatch =2;
  checkcluster=1;
  
  figure(1)
  clf
for ihz= [1,2]
    



  if usePFIorCatch==1
            eventtype='PFI';
            
            %also load BP
  load('GFX_PFIperformance_withSNR_20_min0_RESS.mat');
        else
            eventtype='Catch';
            %load BP data
            load('GFX_Catchperformance_withSNR_20,BPaligned.mat')
  end
  %PLOT correct hz data.
  acrossallPP_TrainingData=[];
  switch ihz
      case 1
          acrossallPP_TrainingData = allppTrainingData_15Hz;
          
      case 2
          acrossallPP_TrainingData = allppTrainingData_20Hz;
      case 3
          acrossallPP_TrainingData = allppTrainingData_30Hz;
      case 4
          acrossallPP_TrainingData = allppTrainingData_40Hz;
      case 5
          acrossallPP_TrainingData = allppTrainingData_45Hz;                    
      case 6
          acrossallPP_TrainingData = allppTrainingData_60Hz;
          
      case 7
          acrossallPP_TrainingData = allppTrainingData_5Hz;
      case 8
          acrossallPP_TrainingData = allppTrainingData_25Hz;
      case 9
          acrossallPP_TrainingData = allppTrainingData_35Hz;
      
      
  end
  
%   if ihz==11 || ihz==12 %average of TG1f or TG2f
%       stacktypes = [1,2,3,4; 5,6,7,8];
%       
%                 ylimsare=[-.05 .5];
% 
%       for ist=1:4
%           itype=stacktypes(ihz-10,ist);
%           
%           switch itype
%               case 1
%                   useme = allppTrainingData_8Hz;
%                   TGcase='TG 1f';
%               case 2
%                   useme = allppTrainingData_13Hz;
%               case 3
%                   useme = allppTrainingData_15Hz;
%               case 4
%                   useme = allppTrainingData_18Hz;
%                   
%                   %or second harms.
%               case 5
%                   useme = allppTrainingData_16Hz;
%                   TGcase='TG 2f';
%               case 6
%                   useme = allppTrainingData_26Hz;
%               case 7
%                   useme = allppTrainingData_30Hz;
%               case 8
%                   useme = allppTrainingData_36Hz;
%           end
%           
%           if usePFIorCatch==1
%               acrossallPP_TrainingData.predicted_onsetSNR(itype,:,:,:)= useme.predicted_onsetSNR;
%               acrossallPP_TrainingData.predicted_offsetSNR(itype,:,:,:)= useme.predicted_offsetSNR;
%               acrossallPP_TrainingData.observed_onsetSNR(itype,:,:,:)= useme.observed_onsetSNR;
%               acrossallPP_TrainingData.observed_offsetSNR(itype,:,:,:)= useme.observed_offsetSNR;
%               
%           else
%               acrossallPP_TrainingData.predicted_BPcatchonsetSNR(itype,:,:,:)= useme.predicted_BPcatchonsetSNR;
%               acrossallPP_TrainingData.observed_BPcatchonsetSNR(itype,:,:,:)=useme.observed_BPcatchonsetSNR;
%               acrossallPP_TrainingData.predicted_BPcatchoffsetSNR(itype,:,:,:)= useme.predicted_BPcatchoffsetSNR;
%               acrossallPP_TrainingData.observed_BPcatchoffsetSNR(itype,:,:,:)=useme.observed_BPcatchoffsetSNR;
%               
%               
%           end
%       end
%       
%       
%       
% %   end
%          
%   
% 
%   %take average if necessary
%   if ihz==11 || ihz ==12
%       if usePFIorCatch==1
%       acrossallPP_TrainingData.predicted_onsetSNR= squeeze(mean(acrossallPP_TrainingData.predicted_onsetSNR,1));
%         acrossallPP_TrainingData.predicted_offsetSNR= squeeze(mean(acrossallPP_TrainingData.predicted_offsetSNR,1));
%         acrossallPP_TrainingData.observed_onsetSNR= squeeze(mean(acrossallPP_TrainingData.observed_onsetSNR,1));
%         acrossallPP_TrainingData.observed_offsetSNR= squeeze(mean(acrossallPP_TrainingData.observed_offsetSNR,1));
%       else
%          acrossallPP_TrainingData.predicted_BPcatchonsetSNR= squeeze(mean(acrossallPP_TrainingData.predicted_BPcatchonsetSNR,1));
%               acrossallPP_TrainingData.observed_BPcatchonsetSNR= squeeze(mean(acrossallPP_TrainingData.observed_BPcatchonsetSNR,1));
%               acrossallPP_TrainingData.predicted_BPcatchoffsetSNR= squeeze(mean(acrossallPP_TrainingData.predicted_BPcatchoffsetSNR,1));
%               acrossallPP_TrainingData.observed_BPcatchoffsetSNR= squeeze(mean(acrossallPP_TrainingData.observed_BPcatchoffsetSNR,1));
%       end 
%           
%           
%   end
%   
  

  
        
%         clf
for ions_off=1:2
    switch ions_off
        case 1 %plot onsets:
            if usePFIorCatch==1
                pred_plotme = acrossallPP_TrainingData.predicted_onsetSNR;
                obs_plotme = acrossallPP_TrainingData.observed_onsetSNR;
            else
                pred_plotme = acrossallPP_TrainingData.predicted_BPcatchonsetSNR;
                obs_plotme = acrossallPP_TrainingData.observed_BPcatchonsetSNR;
            end
            %                       xttimeis= 'button press';
            if usePFIorCatch==1
                xttimeis= 'target invisible';
            else
                xttimeis= 'reporting catch onset';
            end
            targetaction = 'target disappears';
            
            
            
            useBP=storeacrossPpant_onsetBP;
            
            
        case 2 %plot offsets
            if usePFIorCatch==1
                    pred_plotme = acrossallPP_TrainingData.predicted_offsetSNR;
                    obs_plotme = acrossallPP_TrainingData.observed_offsetSNR;
                    else
                                        pred_plotme = acrossallPP_TrainingData.predicted_BPcatchoffsetSNR;
                    obs_plotme = acrossallPP_TrainingData.observed_BPcatchoffsetSNR;
                    end

                    if usePFIorCatch==1
                        xttimeis='target visible';
                    else
                        xttimeis='catch reporting offset';
                    end
    
                    targetaction = 'target reappears';
                    
                    
                      useBP=storeacrossPpant_offsetBP;
            end
            
            %take mean over cross-validations (within ppant).
  pred_plotme = squeeze(nanmean(pred_plotme,2));
  obs_plotme = squeeze(nanmean(obs_plotme,2));
            pleg=[];
            
            
          % %% %%%%%% %%%%% interp bad points if necessary.  
          
       figure(1);       
        if pcount<3
            %plot BP
            subplot(3,2,ions_off)
            
            newDATA=squeeze(mean(useBP,2));
            pMean=squeeze(mean(newDATA,1));
            stE = std(newDATA)/sqrt(size(newDATA,1));
            
        sh=shadedErrorBar(-3:1/60:3, pMean, stE, [], 1);
        colp='r';       
        sh.mainLine.Color=colp;
        sh.mainLine.LineWidth=2;
        sh.patch.FaceColor=colp;
        sh.edge(1).Color=colp;
        sh.edge(2).Color=colp;
        xlabel({['Time from ' xttimeis ]})
        set(gca, 'fontsize', 25)
        axis tight
        ylim([0 2.2])
        ylabel('Buttons Pressed')
        end
            
            anovatestdata=[];
            for plotme=1:2
                switch plotme
                    case 1
                        datais=pred_plotme;
                        colp= 'k';
                        linet=':';
                      
                    case 2
                        datais=obs_plotme;
                        if ihz<6 || ihz==11
                        colp='b';
                        else
                            colp=['m'];
                        end
                        
                        linet='-';
                end
            
                
                
                %adjust for w/in subj errorbars
                newx=datais;
                pmean = nanmean(newx,2); %mean across time points
                gmean = nanmean(pmean); %overall mean
                
                newx = newx - repmat(pmean,[1 size(newx,2)]) + repmat(gmean, size(newx));
                
                
                stE= nanstd(newx)/sqrt(size(newx,1));
                pM= squeeze(nanmean(newx,1));
                
                
        figure(1);       
        subplot(3,2,ions_off+2+2*(pcount-1))
        sh=shadedErrorBar(xt_plot, pM, stE, [], 1);
        sh.mainLine.Color=colp;
        sh.mainLine.LineWidth=2;
        sh.mainLine.LineStyle=linet;
        sh.patch.FaceColor=colp;
        sh.edge(1).Color=colp;
        sh.edge(2).Color=colp;
        hold on
        pleg(plotme)=sh.mainLine;
        
        anovatestdata(plotme,:,:)=datais;
            end

            
            %some stats:
            
               %store ANOVA results
                pvals=nan(1,size(anovatestdata,3));
                Fvals=pvals;
               
                for itime=2:size(anovatestdata,3) %at each timeponit
        
                 [~, pvals(itime),~,stat]=ttest(squeeze(anovatestdata(1,:,itime)),squeeze(anovatestdata(2,:,itime)), .05);
        Fvals(itime)=stat.tstat;
                end

            
               sigs=find(pvals<.05);
%         
%             %perform cluster based correction.
            if length(sigs)>1 && checkcluster==1 
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
                for icl=1:size(clusterSTandEND,1)

                    %start and end are now:
                    STC=sigs(clusterSTandEND(icl,1));
                    ENDC=sigs(clusterSTandEND(icl,2)+1);
                    checktimes =STC:ENDC;
                    observedCV = sum(abs(Fvals(checktimes)));
                    % now shuffle condition labels to see if this cluster is
                    % sig (compared to chance).
                    sumTestStatsShuff = zeros(1,2000);
                    
                    for irand = 1:2000
                        %testing the null that it isn't mismatched - matched at time 2
                        % which creates a diff. so select from either!
                        shD= zeros(size(anovatestdata,1),length(checktimes),size(anovatestdata,2));
                        
                        
                        for ipartition = 1:size(anovatestdata,1)
                            for ippant = 1:size(anovatestdata,2)
                                for itime=1:length(checktimes)
                                    
                                    if mod(randi(100),2)==0 %if random even number
                                        pdata = anovatestdata(1,randi(size(anovatestdata,2)), checktimes(itime)); %select both chans
                                    else
                                        pdata = anovatestdata(2,randi(size(anovatestdata,2)), checktimes(itime));
                                    end
                                    
%                                         pdata = anovatestdata(randi(3),randi(size(anovatestdata,2)), checktimes(itime)); %select both chans
                                    
                                    shD(ipartition,itime,ippant) = pdata;
                                end
                            end
                        end
                       
                        %now compute difference between out hypothetical
                        %data
                        % and test for sig, checking the accumulated test statistic at our
                        % times of interest
                        tvalspertimepoint = zeros(1,length(checktimes));
                        
                        testdata = squeeze(shD(1,:,:)) - squeeze(shD(2,:,:));
                        pvalsNew=[];
                        FvalsNew=[];
                        for itest = 1:length(checktimes) %test each time point
                            
%                             
%                             datat= [squeeze(shD(1,:,itest))'; squeeze(shD(2,:,itest))';...
%                         squeeze(shD(3,:,itest))'];
%                     
%                     anvret=rmanova(datat,factorarray,subjects);
% %                     pvalsNew(itime)=anvret.p(1);
%                     FvalsNew(itime)=anvret.table{2,6};
                    
%                             
%                             tvalspertimepoint(1,itest) = anvret.table{2,6};
                            
                            
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
                    
                    %observed pvalue in distribution=
                    [~, c2] = min(abs(H.Data-observedCV)); %closest value.
                    pvalis= 1-cdf(c2);
                    title(['sum tvals = ' num2str(observedCV), 'p=' num2str(pvalis)]);
                    
                    %              title('Spatial Cluster  Significant!')
                        for itime=checktimes
                            figure(1);
                            hold on
                            plot(xt_plot(itime), sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', colp)
                        end
                    end
            
%            
            end
            else
                try %just plotem
                for itime=sigs
                            figure(1);
                            hold on
%                             plot(timeidDYN(itime), sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', col)
                end
                catch
                end
                
            end
            
            
            
            
            
if ions_off==2
       legend([pleg(1) ,pleg(2) ], {['Predicted'],['Observed']})
end
       set(gca, 'fontsize', 15)
       if ihz<=10
        ylabel({['RESS log(SNR)'];[ num2str(peakfreqsare(ihz)) ' Hz']})
       else
           ylabel(['RESS log(SNR) ' TGcase  ])
       end
%         xlabel({['Time from ' xttimeis ' (' eventtype ' ' targetaction ')']})
        xlabel({['Time from ' xttimeis ]})
        axis tight
ytix=get(gca, 'ytick');

% ylim([.8 1.3])
axis tight
if usePFIorCatch==2
%     ylimsare(2)= ylimsare(2)+.2;
end
ylim([ylimsare])

hold on; plot([0 0 ], ylim, ['k:'], 'linew', 2)
        set(gca, 'fontsize', 25)
        
        
        end
  %%
        
      
        %%
        set(gcf, 'color', 'w')
        cd(basefol)
        cd('newplots-MD')
        cd('Figures')
        cd('Predicted vs Observed')
%         if ihz<=10
%         print('-dpng', [eventtype ' results ' num2str(peakfreqsare(ihz)) ' Hz excl<30, interp'])
%         else
%             print('-dpng', [eventtype ' results ' TGcase ' Hz excl<30'])
%         end
print('-dpng', 'Catch Pred vs Observed summary')

        pcount=pcount+1;
        end
end
        %%
if job.plotPredvsObserved_acrossparticipants_ver2==1
    
    
  
  cd(basefol)
    cd('EEG')
    cd('GFX_EEG-RESS')
    
    %%
  load('GFX_predictedvsobserved')
 
  %% which HZ to look at?



% peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; 

  

  %% which HZ to look at?
  clf
  rsqacross=zeros(4,22);
  rsqcount=1;
  pcount=1;
  
for usePFIorCatch =1:2
  
  %select one of th ebelow
  useScatter=1;
  useCorrovertime=0;
  
  figure(1)
%   clf
  
  overlay=1; %same plots
  barcount=1;

    
  


%%
cd(basefol)  
cd('EEG')
cd('GFX_Pre-RESS')
%%
if usePFIorCatch==1
    eventtype='PFI';
    
    %also load BP
    load('GFX_PFIperformance_withSNR_15_min0.mat');
else
    eventtype='Catch';
    %load BP data
    load('GFX_Catchperformance_withSNR_15,BPaligned.mat')
end
  %PLOT correct hz data.
  acrossallPP_TrainingData=[];

  for ihz= [1,3]
  switch ihz
      case 1
          acrossallPP_TrainingData = allppTrainingData_15Hz;
%           ylimsare=[.75 1.2];
      case 2
          acrossallPP_TrainingData = allppTrainingData_20Hz;
      case 3
          acrossallPP_TrainingData = allppTrainingData_30Hz;
      case 4
          acrossallPP_TrainingData = allppTrainingData_40Hz;
      case 5
          acrossallPP_TrainingData = allppTrainingData_45Hz;
          
      case 6
          acrossallPP_TrainingData = allppTrainingData_60Hz;
          
      case 7
          acrossallPP_TrainingData = allppTrainingData_5Hz;
      case 8
          acrossallPP_TrainingData = allppTrainingData_25Hz;
      case 9
          acrossallPP_TrainingData = allppTrainingData_35Hz;
      
  end
%   
%   if ihz==11 || ihz==12 %average of TG1f or TG2f
%       stacktypes = [1,2,3,4; 5,6,7,8];
%       
%                 ylimsare=[-.05 .5];
% 
%       for ist=1:4
%           itype=stacktypes(ihz-10,ist);
%           
%           switch itype
%               case 1
%                   useme = allppTrainingData_8Hz;
%                   TGcase='TG 1f';
%               case 2
%                   useme = allppTrainingData_13Hz;
%               case 3
%                   useme = allppTrainingData_15Hz;
%               case 4
%                   useme = allppTrainingData_18Hz;
%                   
%                   %or second harms.
%               case 5
%                   useme = allppTrainingData_16Hz;
%                   TGcase='TG 2f';
%               case 6
%                   useme = allppTrainingData_26Hz;
%               case 7
%                   useme = allppTrainingData_30Hz;
%               case 8
%                   useme = allppTrainingData_36Hz;
%           end
%           
%           if usePFIorCatch==1
%               acrossallPP_TrainingData.predicted_onsetSNR(itype,:,:,:)= useme.predicted_onsetSNR;
%               acrossallPP_TrainingData.predicted_offsetSNR(itype,:,:,:)= useme.predicted_offsetSNR;
%               acrossallPP_TrainingData.observed_onsetSNR(itype,:,:,:)= useme.observed_onsetSNR;
%               acrossallPP_TrainingData.observed_offsetSNR(itype,:,:,:)= useme.observed_offsetSNR;
%               
%           else
%               acrossallPP_TrainingData.predicted_BPcatchonsetSNR(itype,:,:,:)= useme.predicted_BPcatchonsetSNR;
%               acrossallPP_TrainingData.observed_BPcatchonsetSNR(itype,:,:,:)=useme.observed_BPcatchonsetSNR;
%               acrossallPP_TrainingData.predicted_BPcatchoffsetSNR(itype,:,:,:)= useme.predicted_BPcatchoffsetSNR;
%               acrossallPP_TrainingData.observed_BPcatchoffsetSNR(itype,:,:,:)=useme.observed_BPcatchoffsetSNR;
%               
%               
%           end
%       end
%       
%       
%       
%   end
%          
%   
% 
%   %take average if necessary
%   if ihz==11 || ihz ==12
%       if usePFIorCatch==1
%       acrossallPP_TrainingData.predicted_onsetSNR= squeeze(mean(acrossallPP_TrainingData.predicted_onsetSNR,1));
%         acrossallPP_TrainingData.predicted_offsetSNR= squeeze(mean(acrossallPP_TrainingData.predicted_offsetSNR,1));
%         acrossallPP_TrainingData.observed_onsetSNR= squeeze(mean(acrossallPP_TrainingData.observed_onsetSNR,1));
%         acrossallPP_TrainingData.observed_offsetSNR= squeeze(mean(acrossallPP_TrainingData.observed_offsetSNR,1));
%       else
%          acrossallPP_TrainingData.predicted_BPcatchonsetSNR= squeeze(mean(acrossallPP_TrainingData.predicted_BPcatchonsetSNR,1));
%               acrossallPP_TrainingData.observed_BPcatchonsetSNR= squeeze(mean(acrossallPP_TrainingData.observed_BPcatchonsetSNR,1));
%               acrossallPP_TrainingData.predicted_BPcatchoffsetSNR= squeeze(mean(acrossallPP_TrainingData.predicted_BPcatchoffsetSNR,1));
%               acrossallPP_TrainingData.observed_BPcatchoffsetSNR= squeeze(mean(acrossallPP_TrainingData.observed_BPcatchoffsetSNR,1));
%       end 
%           
%           
%   end
%   
%   

  
        
%         clf

for ions_off=1%:2
    switch ions_off
        case 1 %plot onsets:
            if usePFIorCatch==1
                pred_plotme = acrossallPP_TrainingData.predicted_onsetSNR;
                obs_plotme = acrossallPP_TrainingData.observed_onsetSNR;
            else
                pred_plotme = acrossallPP_TrainingData.predicted_BPcatchonsetSNR;
                obs_plotme = acrossallPP_TrainingData.observed_BPcatchonsetSNR;
            end
            %                       xttimeis= 'button press';
            if usePFIorCatch==1
                xttimeis= 'target invisible';
            else
                xttimeis= 'reporting catch onset';
            end
            targetaction = 'target disappears';
            
            
            
            useBP=storeacrossPpant_onsetBP;
            
            
        case 2 %plot offsets
            if usePFIorCatch==1
                    pred_plotme = acrossallPP_TrainingData.predicted_offsetSNR;
                    obs_plotme = acrossallPP_TrainingData.observed_offsetSNR;
                    else
                                        pred_plotme = acrossallPP_TrainingData.predicted_BPcatchoffsetSNR;
                    obs_plotme = acrossallPP_TrainingData.observed_BPcatchoffsetSNR;
                    end

                    if usePFIorCatch==1
                        xttimeis='target visible';
                    else
                        xttimeis='catch reporting offset';
                    end
    
                    targetaction = 'target reappears';
                    
                    
                      useBP=storeacrossPpant_offsetBP;
            end
            
            %take mean over cross-validations (within ppant).
  pred_plotme = squeeze(nanmean(pred_plotme,2));
  obs_plotme = squeeze(nanmean(obs_plotme,2));
            pleg=[];
            
%             
%           % %% %%%%%% %%%%% interp bad points if necessary.  
%           if interpCatch==1 && (ihz ==10) && ions_off==1 && usePFIorCatch ==2
%               %may need to interpolate at trial level. see what this looks like.
%               %time vector
%               tvector=xt_plot;
%               points= 1:length(tvector);
%               %remove bad section from trace.
%               %we used a one second sliding window, so:
%               
%                   badsec = dsearchn(tvector', [-1.45 -.36]');
%               
%               
%               tmp=points;
%               tmp(badsec(1):badsec(2))=[];
%               %input x vector with points missing.
%               xIN=tmp;
%               for iptrial=1:size(obs_plotme,1)
%                   reald= obs_plotme(iptrial,:);
%                   
%                   %remove 'bad' points
%                   realrem = [reald(1, 1:(badsec(1)-1)), reald(1,badsec(2)+1:end)];
%                   
%                   %interp
%                   querypoints = badsec(1):badsec(2);
%                   y=interp1(xIN, realrem, querypoints, 'spline');
%                   
%                   %now replace the data with the interpolated values
%                   
%                   reald2=reald;
%                   reald2(1,badsec(1):badsec(2))=y;
%                   
%                   %plot for comparison.
%                   %       clf;
%                   %       plot(reald); hold on; plot(reald2)
%                   
%                   
%                   %      replace for plotting
%                   
%                   obs_plotme(iptrial,:)=reald2;
%                   
%               end
%               
%           end
%           
       figure(1);       
        if pcount<3
            %plot BP
%             subplot(4,2,ions_off)
            subplot(3,2, usePFIorCatch)
            
            newDATA=squeeze(mean(useBP,2));
            pMean=squeeze(mean(newDATA,1));
            stE = std(newDATA)/sqrt(size(newDATA,1));
            
        sh=shadedErrorBar(-3:1/60:3, pMean, stE, [], 1);
        colp='r';       
        sh.mainLine.Color=colp;
        sh.mainLine.LineWidth=2;
        sh.patch.FaceColor=colp;
        sh.edge(1).Color=colp;
        sh.edge(2).Color=colp;
        xlabel({['Time from ' xttimeis ]})
        set(gca, 'fontsize', 20)
        axis tight
        ylim([0 3.5])
        ylabel('Buttons Pressed')
        end
            
            anovatestdata=[];
            
            for plotme=1:2
                switch plotme
                    case 1
                        datais=pred_plotme;
                        colp= 'k';
                        linet=':';
                      
                    case 2
                        datais=obs_plotme;
                        if ihz<7 && ihz>2
                        colp='m';
                        elseif ihz>=7
                            colp=[.5 .5 .5];
                        elseif ihz<=2
                        colp='b';
                        end
                        linet='-';
                end
            
                
                
                %adjust for w/in subj errorbars
                newx=datais;
                pmean = nanmean(newx,2); %mean across time points
                gmean = nanmean(pmean); %overall mean
                
                newx = newx - repmat(pmean,[1 size(newx,2)]) + repmat(gmean, size(newx));
                
                
                stE= nanstd(newx)/sqrt(size(newx,1));
                pM= squeeze(nanmean(newx,1));
                
                
        figure(1);       
        if overlay~=1
        subplot(3,2,ions_off+2+2*(pcount-1))
        else
%             subplot(4,2,[ions_off+2,ions_off+4])
%             subplot(4,2,[ions_off+2])
            subplot(3,2, usePFIorCatch+2)
        end
            
        sh=shadedErrorBar(xt_plot, pM, stE, [], 1);
        sh.mainLine.Color=colp;
        sh.mainLine.LineWidth=2;
        sh.mainLine.LineStyle=linet;
        sh.patch.FaceColor=colp;
        sh.edge(1).Color=colp;
        sh.edge(2).Color=colp;
        hold on
        pleg(plotme)=sh.mainLine;
        
%         anovatestdata(plotme,:,:)=datais;
            if ihz==5 && plotme==2% 2f
                firstharm=pleg(2);
            end
            
            
            % add plot details
            
            
            
            if ions_off==2 && ihz==10 && plotme==2
%                 legend([pleg(1) , firstharm, pleg(2) ], {['Predicted'],['1f'], ['2f']})
            end
            set(gca, 'fontsize', 15)
            if ihz<=10
%                 ylabel({['RESS log(SNR)'];[ num2str(peakfreqsare(ihz)) ' Hz']})
                ylabel({['RESS log(SNR)']})
            else
                ylabel(['RESS log(SNR) ' TGcase  ])
            end
            %         xlabel({['Time from ' xttimeis ' (' eventtype ' ' targetaction ')']})
            xlabel({['Time from ' xttimeis ]})
            axis tight
            ytix=get(gca, 'ytick');
            
            % ylim([.8 1.3])
            axis tight
            if usePFIorCatch==2
                %     ylimsare(2)= ylimsare(2)+.2;
            end
%             ylim([ylimsare])
            
            hold on; plot([0 0 ], ylim, ['k:'], 'linew', 2)
            set(gca, 'fontsize', 20)
            
            
           
            %% now add the correlation values.
            if plotme==2 && useScatter==1
%             subplot(4,2, ions_off+4)
              subplot(3,2, 5:6)
              
              
            hold on;
            
            rsqsacrossppants= zeros(22,1);
            for ippant =  1:size(pred_plotme)
                
            scY = pred_plotme(ippant,:);
            scX = obs_plotme(ippant,:);
            
            
            % calculate R^2
            %SSTotal = SSR + SSE
            
            SSR = sum((scY-mean(scX)).^2);
            %SSErr = error sum of squares, quant. how much the data vary around regression line.
            
            SSerr = sum((scX-scY).^2);
            %SStot = titak sum of squares, data points about mean.
            SStot = SSR+SSerr;
            
            
            rsq = 1-SSerr/SStot;
            
            rsqsacrossppants(ippant)=rsq;
            end
            
            rsq=rsqsacrossppants;
%             %%
%             
%             sch= scatter(scX, scY);
%             sch.MarkerFaceColor = colp;
%             sch.MarkerEdgeColor=colp;
%             
%             [p, S] = polyfit(scX, scY, 1); %linear fit
%              f1=polyval(p,scX);
%              hold on;
% %              plot(scX, f1, 'color',colp, 'linew', 6);
%             
%             xlabel('Observed RESS logSNR')
%             ylabel('Predicted');
%               set(gca, 'fontsize', 20)
%             
%                subplot(4,2, ions_off+6)
%                % plot results from across ppant poly fit.
%                rsq=zeros(22,1);
%                for ippant = 1:size(pred_plotme)
%                    
%                    cofs =corrcoef(pred_plotme(ippant,:), obs_plotme(ippant,:));
%                    
%                    rsq(ippant)= cofs(2,1)^2;                   
%                end
%                
%                 hold on
                barcount = ihz/5 + 2*(usePFIorCatch-1);
               bh=bar(barcount, mean(rsq,1));
               bh.FaceColor=colp;
               errorbar(barcount, mean(rsq,1), std(rsq), 'k');
               barcount=barcount+1;
               
%%               

try set(gca, 'xtick', [1,2,3,4], 'xticklabel', {'1f', '2f', '1f', '2f'})
catch
end
               xlabel('Prediction accuracy')
               
               ylabel('R^2')
               ylim([-.05 1])
                 set(gca, 'fontsize', 20)
                 
                 
                 %store the rsq for sig tests.
                 rsqacross(rsqcount,:) = rsq;
                 rsqcount=rsqcount+1;
            elseif useCorrovertime==1 && plotme==2
                 
%                 subplot(4,2, ions_off+4)                
                    subplot(4,2, usePFIorCatch+4)
                    
                %plot mean correlation over time
                mCorrR = zeros(1, size(pred_plotme,2));
                %%
                
                %subtract mean first.
                pred_plotme = pred_plotme- repmat(mean(pred_plotme,2), 1,41);
                obs_plotme = obs_plotme- repmat(mean(obs_plotme,2), 1,41);
                %
%                 figure();
%                 plot(mean(pred_plotme,1)); hold on; plot(mean(obs_plotme,1));
                %%
                
                for it=1:size(pred_plotme,2)
                   tmpX=pred_plotme(:,it);
                   tmpY=obs_plotme(:,it);
                   
                    
                   [r,p ]=corrcoef(tmpX,tmpY);
                    
                   mCorrR(it)= r(1,2);
                    
                end
                
               hold on; 
               %%%%%%%% plot CORR
               plot(xt_plot, mCorrR, 'color', colp, 'linew', 5); 
               %%
           %no construct shuffled time course from bootstrap.
           shufftimecourse=zeros(200,size(pred_plotme,2));    
           for ishuff=1:200
                  shindx = randi(22, [1,22]); 
                  %plot the new
                  pred_plotmetmp= pred_plotme(shindx,:);
                  obs_plotmetmp = obs_plotme(shindx,:);
                  
                   mCorr = zeros(1, size(pred_plotme,2));
                
                   for it=1:size(pred_plotme,2)
                       tmpX=pred_plotmetmp(:,it);
                       tmpY=obs_plotmetmp(:,it);
                       
                       [r,p ]=corrcoef(tmpX,tmpY);
                       
                       shufftimecourse(ishuff,it)= r(1,2);
                       
                                                                    
                   end
                   
                   
           end
           
           % now plot this as the CI.
           %take 95% CI
           Shuffupperbound=zeros(1,size(shufftimecourse,2));
           Shufflowerbound= Shuffupperbound;
            for ip = 1:size(shufftimecourse,2)
                     
            %%
            shd= sort(shufftimecourse(:,ip)');
            
            %note not normally distributed!
            %try logit transform
%             lp = log(shd)-log(1-shd);
                
            %now rearrange so that we can place correct scores.
%               [lpn,id]=sort(lp);%convert to zscore

              %sort first.
            tz= zscore((shd));
                        
            %find closest to zscores to 99% CI (2.57), 95% CI= 1.96;
%             AB=dsearchn(tz, [-2.57 2.57]');
            AB=dsearchn(tz', [-1.96 1.96]');              
            %this is upper and lower bound.
            CI= [shd(AB(1)), shd(AB(2))];
            
            %now convert back from logit.
%CI upper bound (95%)
% Shuffupperbound(1,ip) = exp(CI(2)) ./ (exp(CI(2))+1);
% Shufflowerbound(1,ip) = exp(CI(1)) ./ (exp(CI(1))+1);
Shuffupperbound(1,ip) = CI(2);
Shufflowerbound(1,ip) = CI(1);
          

            end
            
            x=xt_plot;
            %now plot using patch!
           %Calculate the error bars
            uE=Shuffupperbound;
            lE=Shufflowerbound;
            
            %Make the patch
            yP=[lE,fliplr(uE)];
            xP=[x,fliplr(x)];
            
            %remove nans otherwise patch won't work
            xP(isnan(yP))=[];
            yP(isnan(yP))=[];
            
            
            H.patch=patch(xP,yP,1,'facecolor',colp,...
                'edgecolor','none',...
                'facealpha',.15);
            
            
            %Make pretty edges around the patch.
            H.edge(1)=plot(x,lE,'-','color',colp);
            H.edge(2)=plot(x,uE,'-','color',colp);
               
            ylim([-1 1])
            xlim([-3 3])
            ylabel('corrcoef')
            set(gca, 'fontsize', 20);
            end
                
                
                
                
                
                
                
            end

            
            

        
        
        
        
        
        
        end
  %%
        
      
        %%
      
end
end

% now test for rsq significance.
%% set up rmanova
subjects = repmat([1:22], [size(rsqacross,1),1]);
%one column vector with factors
freqfactors = repmat([1,2,1,2]', [1, 22]);
modelfactors = repmat([1,1,2,2]',[1,22]);
%reshape all into single column.
datacol= reshape(rsqacross, [size(subjects,1)*size(subjects,2),1]);
subjectscol= reshape(subjects, [size(subjects,1)*size(subjects,2),1]);
withinfactor= reshape(freqfactors, [size(subjects,1)*size(subjects,2),1]);
bwssfactor = reshape(modelfactors, [size(subjects,1)*size(subjects,2),1]);
factorarray = [withinfactor,bwssfactor];
%perform rmANOVA

ret=rmanova(datacol, factorarray, subjectscol,{'Hz', 'Model'},2);
%%
set(gcf, 'color', 'w')
        cd(basefol)
        cd('newplots-MD')
        cd('Figures')
        cd('Predicted vs Observed')
if  ions_off==1
print('-dpng', 'PFI Pred vs Observed summary, for onsets')
else
    print('-dpng', 'Catch Pred vs Observed summary, for offsets')
end

        pcount=pcount+1;
end
        %%
          
        
        
        %%%%% %%%%%% modelled similarly to the previous pred vs observed.
        % now we use classification accuracy as our data. 
        
        %first, we just need the median SNR per freq, per whole trial, from training data
        %(75% subset of trials).
        if job.processharmonicClassificationacc_perppant==1
            
            cd(basefol)
            %first load BP data to extract disap and reapp.
            cd('newplots-MD')
            load('MD_AllBP_perppant')
            load('PFI_data')
            peakfreqsare=[8,13,15,18,20,40,16,26,30,36];
            %the next is modelled on job 1 in s3_C_EpochPFIEEG:
            
            for  ifol=allppants %this is ppant of interest.
                
                for ihz=1:10%length(peakfreqsare)
                    
                    useHz=peakfreqsare(ihz);
                    if ihz>6
                        useHz=peakfreqsare(ihz)/2; %ie second harmonic response, but first harmonic target flicker.
                    end
                    
                    cd(basefol)
                    cd(num2str(ifol))
                    load('P_wholetrialRESS'); %load snr data
                    
                    
                    
                    window=[-3 3];
                    
                    srate=1/min(unique(diff(tgrm_wholetrial)));
                    epochdur = sum(abs(window))*srate;
                    onsetc = ceil(epochdur)/2;
                    
                    Setlist=[];
                    TrainingDatatmp=[];
                    
                    for iSet=1:10
                        
                        
                        
                        %Random selection of 2/3 cross validation, each set list.
                        randp= randi(24,[ 1 18]);
                        while unique(randp)<24
                            randp=randi(24,[1 18]);
                        end
                        
                        Setlist(iSet).training_tr=randp;
                        %this script captures the real SNR, so perform for both the
                        %training and remaining data.
                        for iobs=1:2
                            if iobs==1
                                %first time, collect ERP for training data
                                trials = Setlist(iSet).training_tr;
                            else
                                %second time, collect ERPfor remaining trials
                                tmp=1:24;
                                tmp(Setlist(iSet).training_tr)=[];
                                trials=tmp;
                            end
                            
                            
                            %make sure we index the correct SNR per harmonic.
                            if ihz>6
                                ehz=ihz-6;
                            else
                                ehz=ihz;
                            end
                            tmp=zeros(1,length(trials));
                            tcount=1;
                            for itrial=trials
                                
                                
                                SNRdata = squeeze(snr_wholetrRESS(ehz,itrial,:))';
                                
                                trialdata= PFI_only(ifol).Trial(itrial);
                                if trialdata.Goodtrial==1
                                    %this is all we need. we will use the cross validation of 'median' values,
                                    % to assess classification accuracy.
                                    medianSNR=nanmedian(SNRdata);
                                    
                                    
                                    tmp(tcount)=medianSNR;
                                    tcount=tcount+1;
                                end
                            end
                            
                            %store this average.
                            switch iobs
                                case 1
                                    TrainingDatatmp(iSet).Training_medianSNR=squeeze(nanmean(tmp));
                                    
                                case 2
                                    %store the observed data for the remaining trials, to
                                    %compare next with 'predicted'.
                                    TrainingDatatmp(iSet).Observed_medianSNR=squeeze(nanmean(tmp));
                                    
                            end
                        end
                        
                        %
                    end
                    
                    %save results.
                    savename= ['PFI_SNR_ClassificationAccuracy_allHz'];
                    switch ihz
                        case 1
                            TrainingData_8Hz_median = TrainingDatatmp;
                        case 2
                            TrainingData_13Hz_median = TrainingDatatmp;
                        case 3
                            TrainingData_15Hz_median = TrainingDatatmp;
                        case 4
                            TrainingData_18Hz_median = TrainingDatatmp;
                        case 5
                            TrainingData_20Hz_median = TrainingDatatmp;
                        case 6
                            TrainingData_40Hz_median = TrainingDatatmp;
                        case 7
                            TrainingData_16Hz_median = TrainingDatatmp;
                        case 8
                            TrainingData_26Hz_median = TrainingDatatmp;
                        case 9
                            TrainingData_30Hz_median = TrainingDatatmp;
                        case 10
                            TrainingData_36Hz_median = TrainingDatatmp;
                            
                    end
                    
                end
                
                save(savename,  'TrainingData_8Hz_median',...
                    'TrainingData_13Hz_median',...
                    'TrainingData_15Hz_median',...
                    'TrainingData_18Hz_median',...
                    'TrainingData_20Hz_median',...
                    'TrainingData_40Hz_median',...
                    'TrainingData_16Hz_median',...
                    'TrainingData_26Hz_median',...
                    'TrainingData_30Hz_median',...
                    'TrainingData_36Hz_median',...
                    'Setlist');
                
            end
        end
        
        
        
        % now use these subsets of median SNR, to predict classification
        % accuracy in remaining trials:

if job.processharmonicClassificationacc_Predicted_perppant==1
            
          cd(basefol)
            %first load BP data to extract disap and reapp.
            cd('newplots-MD')
            load('MD_AllBP_perppant')
            load('PFI_data')
            peakfreqsare=[8,13,15,18,20,40,16,26,30,36];
            %the next is modelled on job 1 in s3_C_EpochPFIEEG:
            
            for  ifol=allppants %this is ppant of interest.
                cd(basefol)
              cd(num2str(ifol))
        load('PFI_SNR_ClassificationAccuracy_allHz.mat')
        
        load('P_wholetrialRESS.mat', 'snr_wholetrRESS','tgrm_wholetrial');
        
        window=[-3 3];
        
        srate=1/min(unique(diff(tgrm_wholetrial)));
        epochdur = sum(abs(window))*srate;
        onsetc = ceil(epochdur)/2;
        
        for ihz=1:10
            switch ihz
                case 1
                    TrainingDatatmp= TrainingData_8Hz_median;
                case 2
                    
                    TrainingDatatmp= TrainingData_13Hz_median;
                case 3
                    TrainingDatatmp= TrainingData_15Hz_median;
                case 4
                    
                    TrainingDatatmp= TrainingData_18Hz_median;
                case 5
                    TrainingDatatmp= TrainingData_20Hz_median;
                case 6
                    TrainingDatatmp= TrainingData_40Hz_median;
                case 7
                    TrainingDatatmp= TrainingData_16Hz_median;
                case 8
                    
                    TrainingDatatmp= TrainingData_26Hz_median;
                case 9
                    TrainingDatatmp= TrainingData_30Hz_median;
                case 10
                    
                    TrainingDatatmp= TrainingData_36Hz_median;
            
            end
            useHz=peakfreqsare(ihz);
            if ihz>6
                useHz=peakfreqsare(ihz)/2; %ie second harmonic response, but first harmonic target flicker.
            end
            
        
        cd(basefol)
        cd(num2str(ifol))
        %load the observed training data and used Set list.
        
        
        predictedSSVEP=[]; %for storage.
        
        for iSet=1:size(Setlist,2)
            
            
            %find which trials were not used in this list:
            allt=[1:24];
            %rem trained trials.
            allt(Setlist(iSet).training_tr)=[];
            itrialcount=1;
                            
            
            
            %what is the median SNR from this subset of trials we are
            %testing?
            medianSNR=TrainingDatatmp(iSet).Training_medianSNR;
            
            
            %get ready for data collecton (across trials)
            %separate by number involved
            %disappearances:
            ppant_SNREEG_PFI_0_1 = [];
            counter0_1=1;
            ppant_SNREEG_PFI_1_2 = [];
            counter1_2=1;
            
            ppant_SNREEG_PFI_2_3 = [];
            counter2_3=1;
            %Reappearances
            ppant_SNREEG_PFI_3_2 = [];
            counter3_2=1;
            ppant_SNREEG_PFI_2_1 = [];
            counter2_1=1;
            ppant_SNREEG_PFI_1_0 = [];
            counter1_0=1;
            
            durs0_1=0;
            durs1_0=0;
            durs1_2=0;
            durs2_1=0;
            durs2_3=0;
            durs3_2=0;
            
            BPs0_1=[];
            BPs1_0=[];
            BPs1_2=[];
            BPs2_1=[];
            BPs2_3=[];
            BPs3_2=[];
            
            Freqwas=[]; % make structure for later extraction of relevant freqs in analysis.
            Locationwas=[]; %same for location of PFI events BPs.
            
            %make sure we index the correct SNR per harmonic.
            if ihz>6
                ehz=ihz-6;
            else
                ehz=ihz;
            end
            
            for itrial=allt
                
                
                
                %we look at the real SNR data, in untrained trials, to see
                %if median SNR from trained trials predicts Button press.
                                SNRdata = squeeze(snr_wholetrRESS(ehz,itrial,:))';
                
                trialdata= PFI_only(ifol).Trial(itrial);
                if trialdata.Goodtrial==1
                    % start with some trial/BP
                    % classification.
                    
                    accumBP = nansum(trialdata.allBPs,1);
                    
                    BPtimes= [0:length(accumBP)]/60;
                    %                                     %disap is all >0
                    %                                     disapBP= find(accumBP>0);
                    %                                     reapBP= find(accumBP==0);
                    %
                    %                                     %find these events in seconds:
                    %                                     disapBP_secs = BPtimes(disapBP);
                    %                                     reapBP_secs = BPtimes(reapBP);
                    %
                    %now we want to check whether </>
                    %median SNR predicts BP.
                    
                    
                    above_medSNR=find(SNRdata>medianSNR);
                    below_medSNR=find(SNRdata<=medianSNR);
                    
                    %find these events in seconds:
                    above_medSNR_secs= tgrm_wholetrial(above_medSNR);
                    below_medSNR_secs= tgrm_wholetrial(below_medSNR);
                    
                    
                    % now for both cases, assign class=0 if
                    % congruent with BP.
                    %first for disap.
                    classAcc=zeros(1,length(SNRdata));
                    
                    for idtype=1%:2
                        switch idtype
                            case 1
                                dSNR=above_medSNR_secs;
                            case 2
                                dSNR=below_medSNR_secs;
                        end
                        
                        for i=1:length(above_medSNR_secs) %shorter vector
                            % find nearest time point and check
                            % classification.
                            timetmp=dSNR(i);
                            
                            nearpointBEH = dsearchn([BPtimes'], timetmp');
                            nearpointEEG= dsearchn([tgrm_wholetrial'], timetmp');
                            
                            %check!
                            %want above vs below
                            switch idtype
                                case 1
                                    if accumBP(nearpointBEH) >=0
                                        classAcc(nearpointEEG)=1;
                                    end
                                case 2
                                    if accumBP(nearpointBEH) ==0
                                        classAcc(nearpointEEG)=1;
                                    end
                            end
                        end
                        
                        
                    end
                    
                    %now we have our vector based on
                    %classification of BP, using median
                    %split of SNR.
                    
                    % use this as 'predicted SNR'
                    
                    
                    %first establish the direction of changes (disap /reappear
                    %targets).
                    
                    directionBP = diff(accumBP);
                    directionBPtimes = find(directionBP~=0);
                    directionBPtimes=directionBPtimes+1; %account for diff function
                    %now remove zeros, leaving time of all disap/reap
                    directionBP(directionBP==0)=[];
                    
                    
                    
                    %Find Hz locations.
                    Hzare = [trialdata.TL_Freq, trialdata.TR_Freq, trialdata.BL_Freq, trialdata.BR_Freq];
                    
                    
                    %we want to only keep the BPs that survived our previous,
                    %exclusion criteria.
                    %so create vector of all the relevant BP frames already stored.
                    survivingPFI = [trialdata.PFI_disap_0target_framestart,...
                        trialdata.PFI_disap_1target_framestart,...
                        trialdata.PFI_disap_2target_framestart,...
                        trialdata.PFI_disap_3ormoretarget_framestart];
                    
                    
                    
                    %so now only continue by assessing those BPs (ie remove
                    %transients)
                    keepdirectionBPtimes = find(ismember(directionBPtimes, survivingPFI));
                    directionBPtimes= directionBPtimes(keepdirectionBPtimes);
                    
                    
                    
                    for iPFI = 1:length(directionBPtimes)
                        %Collect event info.
                        timeis = directionBPtimes(iPFI);
                        perceptis = accumBP(directionBPtimes(iPFI));
                        perceptwas = accumBP(directionBPtimes(iPFI)-1);
                        
                        
                        outgoingPFI=[];
                        %disap or reap.
                        if perceptwas<perceptis
                            PFIdir=1; %disappearing, an extra button is being pressed
                        else
                            PFIdir=-1; %Reappearing
                        end
                        
                        %gather EEG time index.
                        
                        [~,timePFI]= min(abs(tgrm_wholetrial-timeis/60));
                        
                        
                        
                        
                        %SNRdata to use now is the class
                        %accuracy
                        SNRdata = classAcc;
                        
                        
                        % ignore events at start/end of block.
                        
                        if  timePFI-onsetc>1 && timePFI+onsetc<size(SNRdata,2)
                            %collect data
                            outgoingPFI= SNRdata(:,timePFI-onsetc:timePFI+onsetc);
                            
                            
                            
                            %now store this PFI epoch according to type.
                            switch perceptis %num pressed
                                case 0 %currently NO buttons pressed.
                                    if PFIdir==1
                                        error('checkcode') %cannot go from less than zero button press.
                                    else
                                        if perceptwas==1 % gone from 1 to 0 BP (reappearing)
                                            ppant_SNREEG_PFI_1_0(counter1_0,:,:) = outgoingPFI;
                                            
                                            %how long was this event?
                                            %find the correct framestamp and duration.
                                            
                                            %Store which frequency was involved:
                                            % timestap across BPs.
                                            tmpBPs=trialdata.allBPs(:,timeis-1); % prev BP since reappearing.
                                            
                                            
                                            if ~isnan(sum(tmpBPs))
                                                
                                                Freqwas.dir1_0(counter1_0).d= Hzare(find(tmpBPs));
                                                Locationwas.dir1_0(counter1_0).d= (find(tmpBPs));
                                            else %it was during catch period, so don't count.
                                                
                                                Freqwas.dir1_0(counter1_0).d= nan;
                                                Locationwas.dir1_0(counter1_0).d= nan;
                                            end
                                            Freqwas.dir1_0(counter1_0).Trialind= itrial;
                                            
                                            
                                            
                                            
                                            
                                            
                                            % store BP also
                                            
                                            BPs1_0(counter1_0,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                            
                                            
                                            
                                            counter1_0=counter1_0+1;
                                        elseif perceptwas==2
                                            ppant_SNREEG_PFI_2_1(counter2_1,:,:) = outgoingPFI;
                                            
                                            
                                            
                                            %Store which frequency was involved:
                                            % timestap across BPs.
                                            tmpBPs1=trialdata.allBPs(:,timeis-1);
                                            tmpBPs2=trialdata.allBPs(:,timeis);
                                            dtmp= tmpBPs1-tmpBPs2;
                                            
                                            
                                            
                                            if ~isnan(sum(dtmp))
                                                
                                                Freqwas.dir2_1(counter2_1).d= Hzare(find(dtmp));
                                                Locationwas.dir2_1(counter2_1).d= (find(dtmp));
                                            else %it was during catch period, so don't count.
                                                
                                                Freqwas.dir2_1(counter2_1).d= nan;
                                                Locationwas.dir2_1(counter2_1).d= (find(dtmp));
                                            end
                                            Freqwas.dir2_1(counter2_1).Trialind= itrial;
                                            
                                            
                                            
                                            %also store BP trace.
                                            BPs2_1(counter2_1,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                            
                                            
                                            counter2_1=counter2_1+1;
                                        end
                                    end
                                    
                                case 1 %single target
                                    if PFIdir==1 % PFI numbers increasing
                                        ppant_SNREEG_PFI_0_1(counter0_1,:,:) = outgoingPFI;
                                        
                                        
                                        
                                        
                                        %Store which frequency was involved:
                                        % timestap across BPs.
                                        tmpBPs=trialdata.allBPs(:,timeis);
                                        if ~isnan(sum(tmpBPs))
                                            
                                            Freqwas.dir0_1(counter0_1).d= Hzare(find(tmpBPs));
                                            
                                        else %it was during catch period, so don't count.
                                            
                                            Freqwas.dir0_1(counter0_1).d= nan;
                                        end
                                        Freqwas.dir0_1(counter0_1).Trialind= itrial;
                                        
                                        
                                        %also store BP trace.
                                        BPs0_1(counter0_1,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                        
                                        
                                        counter0_1=counter0_1+1;
                                    else %decreasing numbers being pressed
                                        if perceptwas==2
                                            
                                            ppant_SNREEG_PFI_2_1(counter2_1,:,:) = outgoingPFI;
                                            
                                            
                                            
                                            
                                            %Store which frequency was involved:
                                            % timestap across BPs.
                                            tmpBPs=trialdata.allBPs(:,timeis);
                                            if ~isnan(sum(tmpBPs))
                                                
                                                Freqwas.dir2_1(counter2_1).d= Hzare(find(tmpBPs));
                                                Locationwas.dir2_1(counter2_1).d= (find(tmpBPs));
                                            else
                                                Freqwas.dir2_1(counter2_1).d= nan;
                                                Locationwas.dir2_1(counter2_1).d= nan;
                                            end
                                            
                                            Freqwas.dir2_1(counter2_1).Trialind= itrial;
                                            
                                            
                                            %also store BP trace.
                                            BPs2_1(counter2_1,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                            
                                            counter2_1=counter2_1+1;
                                        elseif perceptwas==3
                                            ppant_SNREEG_PFI_3_2(counter3_2,:,:) = outgoingPFI;
                                            
                                            
                                            
                                            
                                            %Store which frequency was involved:
                                            % timestap across BPs.
                                            tmpBPs1=trialdata.allBPs(:,timeis-1);
                                            tmpBPs2=trialdata.allBPs(:,timeis);
                                            dtmp=tmpBPs1-tmpBPs2;
                                            
                                            
                                            if ~isnan(sum(dtmp))
                                                
                                                Freqwas.dir3_2(counter3_2).d= Hzare(find(dtmp));
                                                Locationwas.dir3_2(counter3_2).d= (find(dtmp));
                                            else %it was during catch period, so don't count.
                                                
                                                Freqwas.dir3_2(counter3_2).d= nan;
                                                Locationwas.dir3_2(counter3_2).d= nan;
                                            end
                                            
                                            Freqwas.dir3_2(counter3_2).Trialind= nan;
                                            
                                            
                                            %also store BP trace.
                                            BPs3_2(counter3_2,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                            
                                            
                                            
                                            counter3_2=counter3_2+1;
                                        end
                                        
                                    end
                                case 2 %double target
                                    if PFIdir==1
                                        ppant_SNREEG_PFI_1_2(counter1_2,:,:) = outgoingPFI;
                                        
                                        
                                        %Store which frequency was involved:
                                        % timestap across BPs.
                                        tmpBPs1=trialdata.allBPs(:,timeis-1);
                                        tmpBPs2=trialdata.allBPs(:,timeis);
                                        dtmp=tmpBPs2-tmpBPs1;
                                        
                                        if ~isnan(sum(dtmp))
                                            
                                            Freqwas.dir1_2(counter1_2).d= Hzare(find(dtmp));
                                            Locationwas.dir1_2(counter1_2).d= (find(dtmp));
                                        else %it was during catch period, so don't count.
                                            
                                            Freqwas.dir1_2(counter1_2).d= nan;
                                            Locationwas.dir1_2(counter1_2).d= nan;
                                        end
                                        
                                        Freqwas.dir1_2(counter1_2).Trialind= itrial;
                                        
                                        
                                        %also store BP trace.
                                        BPs1_2(counter1_2,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                        
                                        
                                        counter1_2=counter1_2+1;
                                    else
                                        ppant_SNREEG_PFI_3_2(counter3_2,:,:) = outgoingPFI;
                                        
                                        
                                        %Store which frequency was involved:
                                        % timestamp across BPs.
                                        tmpBPs1=trialdata.allBPs(:,timeis-1);
                                        tmpBPs2=trialdata.allBPs(:,timeis);
                                        dBPs= tmpBPs1-tmpBPs2;
                                        
                                        
                                        if ~isnan(sum(dBPs))
                                            
                                            Freqwas.dir3_2(counter3_2).d= Hzare(find(dBPs));
                                            Locationwas.dir3_2(counter3_2).d= (find(dBPs));
                                        else %it was during catch period, so don't count.
                                            
                                            Freqwas.dir3_2(counter3_2).d= nan;
                                            Locationwas.dir3_2(counter3_2).d= nan;
                                        end
                                        Freqwas.dir3_2(counter3_2).Trialind= itrial;
                                        
                                        
                                        %also store BP trace.
                                        BPs3_2(counter3_2,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                        
                                        counter3_2=counter3_2+1;
                                        
                                    end
                                    
                                case 3 %triple target
                                    if PFIdir==1
                                        ppant_SNREEG_PFI_2_3(counter2_3,:,:) = outgoingPFI;
                                        
                                        
                                        
                                        
                                        %Store which frequency was involved:
                                        % timestap across BPs.
                                        tmpBPs1=trialdata.allBPs(:,timeis-1);
                                        tmpBPs2=trialdata.allBPs(:,timeis);
                                        dtmp = tmpBPs2-tmpBPs1;
                                        
                                        if ~isnan(sum(dtmp))
                                            
                                            Freqwas.dir2_3(counter2_3).d= Hzare(find(dtmp));
                                            Locationwas.dir2_3(counter2_3).d= (find(dtmp));
                                        else %it was during catch period, so don't count.
                                            
                                            Freqwas.dir2_3(counter2_3).d= nan;
                                            Locationwas.dir2_3(counter2_3).d= nan;
                                        end
                                        
                                        %also store BP trace.
                                        BPs2_3(counter2_3,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                        
                                        Freqwas.dir2_3(counter2_3).Trialind= itrial;
                                        
                                        counter2_3=counter2_3+1;
                                        
                                    else
                                        disp(['Found a 4 ppant ' num2str(ifol) ' trial ' num2str(itrial)])
                                        %
                                    end
                            end
                        end
                    end
                    
                    
                    
                    
                    
                end
                
                
                
            end
            
            
            %critically for TG Hz, we only want to retain the
            %disappearances (ie. BP that involved a particular HZ)
            if ihz~=5 && ihz ~=6 %skip 20 and 40 hz.
                
                for icycle = [1:6] %6 types of disap/reap to check
                    tmp=[];
                    snrd=[];
                    listHz=[];
                    switch icycle
                        case 1
                            % which hz were involved in this each case?
                            try tmp=Freqwas.dir0_1;
                                snrd=ppant_SNREEG_PFI_0_1;
                            catch
                                continue
                            end
                        case 2
                            
                            try tmp=Freqwas.dir1_0;
                                snrd=ppant_SNREEG_PFI_1_0;
                            catch
                                continue
                            end
                        case 3
                            try tmp=Freqwas.dir1_2;
                                snrd=ppant_SNREEG_PFI_1_2;
                            catch
                                continue
                            end
                            
                        case 4
                            try tmp=Freqwas.dir2_1;
                                snrd=ppant_SNREEG_PFI_2_1;
                            catch
                                continue
                            end
                            
                        case 5
                            try tmp=Freqwas.dir2_3;
                                snrd=ppant_SNREEG_PFI_2_3;
                            catch
                                continue
                            end
                            
                        case 6
                            try tmp=Freqwas.dir3_2;
                                snrd=ppant_SNREEG_PFI_3_2;
                            catch
                                continue
                            end
                    end
                    %%
                    keeptrials=[];
                    for item = 1:size(tmp,2) %this is needed, as some disaps had simultaneous BP, messes with indexing.
                        hzc=tmp(item).d;
                        if any(hzc==useHz)
                            keeptrials=[keeptrials, (item)];
                        end
                    end
                    %%
                    %
                    switch icycle
                        case 1
                            ppant_SNREEG_PFI_0_1=snrd(keeptrials,:,:);
                        case 2
                            ppant_SNREEG_PFI_1_0=snrd(keeptrials,:,:);
                        case 3
                            ppant_SNREEG_PFI_1_2=snrd(keeptrials,:,:);
                        case 4
                            ppant_SNREEG_PFI_2_1=snrd(keeptrials,:,:);
                        case 5
                            ppant_SNREEG_PFI_2_3=snrd(keeptrials,:,:);
                        case 6
                            ppant_SNREEG_PFI_3_2=snrd(keeptrials,:,:);
                            
                    end
                    
                    
                end
            end
            
            
            ppantOnsetSNR = squeeze(cat(1, ppant_SNREEG_PFI_0_1,ppant_SNREEG_PFI_1_2,ppant_SNREEG_PFI_2_3));
            
            
            ppantOffsetSNR= squeeze(cat(1, ppant_SNREEG_PFI_1_0,ppant_SNREEG_PFI_2_1,ppant_SNREEG_PFI_3_2));
            
            
            TrainingDatatmp(iSet).PredictedOnsetClassification = squeeze(mean(ppantOnsetSNR,1));
            TrainingDatatmp(iSet).PredictedOffsetClassification = squeeze(mean(ppantOffsetSNR,1));
            
            
          
        xt_plot= linspace(window(1), window(2), size(ppantOnsetSNR,2));
        %
        end
        
        %save results.
        savename= ['PFI_SNR_ClassificationAccuracy_allHz'];
        switch ihz
            case 1
                TrainingData_8Hz_median = TrainingDatatmp;
            case 2
                TrainingData_13Hz_median = TrainingDatatmp;
            case 3
                TrainingData_15Hz_median = TrainingDatatmp;
            case 4
                TrainingData_18Hz_median = TrainingDatatmp;
            case 5
                TrainingData_20Hz_median = TrainingDatatmp;
            case 6
                TrainingData_40Hz_median = TrainingDatatmp;
            case 7
                TrainingData_16Hz_median = TrainingDatatmp;
            case 8
                TrainingData_26Hz_median = TrainingDatatmp;
            case 9
                TrainingData_30Hz_median = TrainingDatatmp;
            case 10
                TrainingData_36Hz_median = TrainingDatatmp;
                
        end
        end


save(savename,  'TrainingData_8Hz_median',...
    'TrainingData_13Hz_median',...
    'TrainingData_15Hz_median',...
    'TrainingData_18Hz_median',...
    'TrainingData_20Hz_median',...
    'TrainingData_40Hz_median',...
    'TrainingData_16Hz_median',...
    'TrainingData_26Hz_median',...
    'TrainingData_30Hz_median',...
    'TrainingData_36Hz_median',...
    'Setlist', 'xt_plot', '-append');

            end
end

if job.concatClassificationacc==1

            acrossallPP_TrainingData=[];
            for ihz=1:10
                
                ifolc=1; %initiate counter.
                for ifol=allppants
                    
                    cd(basefol)
                    cd(num2str(ifol))
                    load('PFI_SNR_ClassificationAccuracy_allHz.mat');
                    
                    
                    switch ihz
                        case 1
                            TrainingDatatmp= TrainingData_8Hz_median;
                        case 2
                            
                            TrainingDatatmp= TrainingData_13Hz_median;
                        case 3
                            TrainingDatatmp= TrainingData_15Hz_median;
                        case 4
                            
                            TrainingDatatmp= TrainingData_18Hz_median;
                        case 5
                            TrainingDatatmp= TrainingData_20Hz_median;
                        case 6
%                               TrainingDatatmp= TrainingData_40Hz;
                            TrainingDatatmp= TrainingData_16Hz_median;
                        case 7
                            
                            TrainingDatatmp= TrainingData_26Hz_median;
                        case 8
                            TrainingDatatmp= TrainingData_30Hz_median;
                        case 9                                                        
                            TrainingDatatmp= TrainingData_36Hz_median;
                        case 10
                            TrainingDatatmp= TrainingData_40Hz_median;
                    end
                    
                    
                    
                    for iset=1:size(Setlist,2)
                        %take ppant mean and store.
                        acrossallPP_TrainingData.predicted_onsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).PredictedOnsetClassification,1));
                        acrossallPP_TrainingData.predicted_offsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).PredictedOffsetClassification,1));
%                         acrossallPP_TrainingData.observed_onsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).Observed_ppantOnsetSNR,1));
%                         acrossallPP_TrainingData.observed_offsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).Observed_ppantOffsetSNR,1));
                        %also add catch data
%                         acrossallPP_TrainingData.predicted_BPcatchonsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).predicted_BPcatchonsetSNR,1));
%                         acrossallPP_TrainingData.predicted_BPcatchoffsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).predicted_BPcatchoffsetSNR,1));
%                         acrossallPP_TrainingData.observed_BPcatchonsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).Observed_BPcatchonsetSNR,1));
%                         acrossallPP_TrainingData.observed_BPcatchoffsetSNR(ifolc,iset,:)= squeeze(nanmean(TrainingDatatmp(1,iset).Observed_BPcatchoffsetSNR,1));
                    end
                    ifolc=ifolc+1; %folder counter.
                end
                    
                
                    switch ihz
                        case 1
                            allppTrainingData_8Hz_median=acrossallPP_TrainingData;
                        case 2
                            allppTrainingData_13Hz_median=acrossallPP_TrainingData;
                        case 3
                            allppTrainingData_15Hz_median=acrossallPP_TrainingData;
                        case 4
                            allppTrainingData_18Hz_median=acrossallPP_TrainingData;
                        case 5
                            allppTrainingData_20Hz_median=acrossallPP_TrainingData;
                        case 6
%                             allppTrainingData_40Hz=acrossallPP_TrainingData;
                            allppTrainingData_16Hz_median=acrossallPP_TrainingData;
                        case 7
                            allppTrainingData_26Hz_median=acrossallPP_TrainingData;
                        case 8
                            allppTrainingData_30Hz_median=acrossallPP_TrainingData;
                        case 9
                            allppTrainingData_36Hz_median=acrossallPP_TrainingData;
                        case 10
                            allppTrainingData_40Hz_median=acrossallPP_TrainingData;
                    end
                    
                
            end
    cd(basefol)
    cd('newplots-md')
    save('GFX_observedClassificationAccuracy', 'allppTrainingData_8Hz_median',...
        'allppTrainingData_13Hz_median',...
'allppTrainingData_15Hz_median',...
'allppTrainingData_18Hz_median',...
'allppTrainingData_20Hz_median',...
'allppTrainingData_16Hz_median',...
'allppTrainingData_26Hz_median',...
'allppTrainingData_30Hz_median',...
'allppTrainingData_36Hz_median',...
'allppTrainingData_40Hz_median','Setlist', 'xt_plot')
    
end

if job.plotClassificationaccuracyacrossppants==1;
    %%
    cd(basefol)
    cd('newplots-md')    
    %%
  load('GFX_observedClassificationAccuracy')
 
  
  interpCatch=0; %can also interpolate bad points in catch onset.
  %% which HZ to look at?
%   peakfreqsare=[8,13,15,18,20,40,16,26,30,36]; %11 = average over targets (1st harm) 12 = 2nd harm
peakfreqsare=[8,13,15,18,20,16,26,30,36,40]; %11 = average over targets (1st harm) 12 = 2nd harm
  clf


pleg=[];
    hzcounter=1; %for legend.        
for ihz= [5,10]%11:12
  
    
    
usePFIorCatch = 1;



  if usePFIorCatch==1
            eventtype='PFI';
        else
            eventtype='Catch';
  end
  %PLOT correct hz data.
  acrossallPP_TrainingData=[];
  switch ihz
      case 1
          acrossallPP_TrainingData = allppTrainingData_8Hz_median;
          
      case 2
          acrossallPP_TrainingData = allppTrainingData_13Hz_median;
      case 3
          acrossallPP_TrainingData = allppTrainingData_15Hz_median;
      case 4
          acrossallPP_TrainingData = allppTrainingData_18Hz_median;
      case 5
          acrossallPP_TrainingData = allppTrainingData_20Hz_median;
          colp='r';
      case 6
          acrossallPP_TrainingData = allppTrainingData_16Hz_median;
      case 7
          acrossallPP_TrainingData = allppTrainingData_26Hz_median;
      case 8
          acrossallPP_TrainingData = allppTrainingData_30Hz_median;
      case 9
          acrossallPP_TrainingData = allppTrainingData_36Hz_median;
      case 10
          acrossallPP_TrainingData = allppTrainingData_40Hz_median;
          colp=[0 .5 0];
      
  end
  
  if ihz==11 || ihz==12 %average of TG1f or TG2f
      stacktypes = [1,2,3,4; 5,6,7,8];
      
                ylimsare=[-.05 .5];

      for ist=1:4
          itype=stacktypes(ihz-10,ist);
          
          switch itype
              case 1
                  useme = allppTrainingData_8Hz;
                  TGcase='TG 1f';
              case 2
                  useme = allppTrainingData_13Hz;
              case 3
                  useme = allppTrainingData_15Hz;
              case 4
                  useme = allppTrainingData_18Hz;
                  
                  %or second harms.
              case 5
                  useme = allppTrainingData_16Hz;
                  TGcase='TG 2f';
              case 6
                  useme = allppTrainingData_26Hz;
              case 7
                  useme = allppTrainingData_30Hz;
              case 8
                  useme = allppTrainingData_36Hz;
          end
          
          if usePFIorCatch==1
              acrossallPP_TrainingData.predicted_onsetSNR(itype,:,:,:)= useme.predicted_onsetSNR;
              acrossallPP_TrainingData.predicted_offsetSNR(itype,:,:,:)= useme.predicted_offsetSNR;
%               acrossallPP_TrainingData.observed_onsetSNR(itype,:,:,:)= useme.observed_onsetSNR;
%               acrossallPP_TrainingData.observed_offsetSNR(itype,:,:,:)= useme.observed_offsetSNR;
              
          else
              acrossallPP_TrainingData.predicted_BPcatchonsetSNR(itype,:,:,:)= useme.predicted_BPcatchonsetSNR;
%               acrossallPP_TrainingData.observed_BPcatchonsetSNR(itype,:,:,:)=useme.observed_BPcatchonsetSNR;
              acrossallPP_TrainingData.predicted_BPcatchoffsetSNR(itype,:,:,:)= useme.predicted_BPcatchoffsetSNR;
%               acrossallPP_TrainingData.observed_BPcatchoffsetSNR(itype,:,:,:)=useme.observed_BPcatchoffsetSNR;
              
              
          end
      end
      
      
      
  end
         
  

  %take average if necessary
  if ihz==11 || ihz ==12
      if usePFIorCatch==1
      acrossallPP_TrainingData.predicted_onsetSNR= squeeze(mean(acrossallPP_TrainingData.predicted_onsetSNR,1));
        acrossallPP_TrainingData.predicted_offsetSNR= squeeze(mean(acrossallPP_TrainingData.predicted_offsetSNR,1));
        acrossallPP_TrainingData.observed_onsetSNR= squeeze(mean(acrossallPP_TrainingData.observed_onsetSNR,1));
        acrossallPP_TrainingData.observed_offsetSNR= squeeze(mean(acrossallPP_TrainingData.observed_offsetSNR,1));
      else
         acrossallPP_TrainingData.predicted_BPcatchonsetSNR= squeeze(mean(acrossallPP_TrainingData.predicted_BPcatchonsetSNR,1));
              acrossallPP_TrainingData.observed_BPcatchonsetSNR= squeeze(mean(acrossallPP_TrainingData.observed_BPcatchonsetSNR,1));
              acrossallPP_TrainingData.predicted_BPcatchoffsetSNR= squeeze(mean(acrossallPP_TrainingData.predicted_BPcatchoffsetSNR,1));
              acrossallPP_TrainingData.observed_BPcatchoffsetSNR= squeeze(mean(acrossallPP_TrainingData.observed_BPcatchoffsetSNR,1));
      end 
          
          
  end
  
  

  
        
        
        for ions_off=1:2
            switch ions_off
                case 1 %plot onsets:
                    if usePFIorCatch==1
                    pred_plotme = acrossallPP_TrainingData.predicted_onsetSNR;
%                     obs_plotme = acrossallPP_TrainingData.observed_onsetSNR;
                    else
                    pred_plotme = acrossallPP_TrainingData.predicted_BPcatchonsetSNR;
%                     obs_plotme = acrossallPP_TrainingData.observed_BPcatchonsetSNR;
                    end
                      xttimeis= 'button press';
                      targetaction = 'target disappears';
                      linep=4;
                case 2 %plot offsets
                    if usePFIorCatch==1
                    pred_plotme = acrossallPP_TrainingData.predicted_offsetSNR;
%                     obs_plotme = acrossallPP_TrainingData.observed_offsetSNR;
                    else
                                        pred_plotme = acrossallPP_TrainingData.predicted_BPcatchoffsetSNR;
%                     obs_plotme = acrossallPP_TrainingData.observed_BPcatchoffsetSNR;
                    end
                    xttimeis= 'button release';
                    targetaction = 'target reappears';
                    linep=2;
            end
            
            %take mean over cross-validations (within ppant).
  pred_plotme = squeeze(nanmean(pred_plotme,2));
%   obs_plotme = squeeze(nanmean(obs_plotme,2));
            
            
          % %% %%%%%% %%%%% interp bad points if necessary.  
          if interpCatch==1 && (ihz ==10) && ions_off==1
              %may need to interpolate at trial level. see what this looks like.
              %time vector
              tvector=xt_plot;
              points= 1:length(tvector);
              %remove bad section from trace.
              %we used a one second sliding window, so:
              
                  badsec = dsearchn(tvector', [-1.45 -.36]');
              
              
              tmp=points;
              tmp(badsec(1):badsec(2))=[];
              %input x vector with points missing.
              xIN=tmp;
              for iptrial=1:size(obs_plotme,1)
                  reald= obs_plotme(iptrial,:);
                  
                  %remove 'bad' points
                  realrem = [reald(1, 1:(badsec(1)-1)), reald(1,badsec(2)+1:end)];
                  
                  %interp
                  querypoints = badsec(1):badsec(2);
                  y=interp1(xIN, realrem, querypoints, 'spline');
                  
                  %now replace the data with the interpolated values
                  
                  reald2=reald;
                  reald2(1,badsec(1):badsec(2))=y;
                  
                  %plot for comparison.
                  %       clf;
                  %       plot(reald); hold on; plot(reald2)
                  
                  
                  %      replace for plotting
                  
                  obs_plotme(iptrial,:)=reald2;
                  
              end
              
          end
          
          
            
            
            
            for plotme=1%:2
                switch plotme
                    case 1
                        datais=pred_plotme;
                        
                      
                    case 2
                        datais=obs_plotme;
                        if ihz<6 || ihz==11
                        colp='r';
                        else
                            colp=[0 .5 0];
                        end
                        
                end
            
                
                
                %adjust for w/in subj errorbars
                newx=datais;
                pmean = nanmean(newx,2); %mean across time points
                gmean = nanmean(pmean); %overall mean
                
                newx = newx - repmat(pmean,[1 size(newx,2)]) + repmat(gmean, size(newx));
                
                
                stE= nanstd(newx)/sqrt(size(newx,1));
                pM= squeeze(nanmean(newx,1));
                
                
        figure(1);
% can plot across ppants:
% for ippant=1:size(newx,1)
    
% subplot(3,8,ippant)
% hold on;
%         sh=shadedErrorBar(xt_plot, squeeze(newx(ippant,:)), ones(1,length(xt_plot))./1000, [], 1);

%or acrossall
        sh=shadedErrorBar(xt_plot, pM, stE, [], 1);
        sh.mainLine.Color=colp;
        sh.mainLine.LineWidth=linep;
        sh.patch.FaceColor=colp;
        sh.edge(1).Color=colp;
        sh.edge(2).Color=colp;
        hold on
        pleg(hzcounter)=sh.mainLine;
        axis tight
        ylim([.4 1])
        plot(xlim,[.5 .5], ['k:'])
%         title(['ppant' num2str(allppants(ippant))], 'fontsize', 15)
% end
            end
%        p1= plot(pred_plotme); hold on;
%        p2= plot(obs_plotme);
%


% how about some stats? compare to class accuracy= 0.5
% p=zeros(1,size(newx,2));
% h=p;
% for it=1:size(datais,2)
%     
%     [h(it),p(it)] = ttest(datais(:,it), 0.5);
%     
%     
% end

% %plot results (basic)
% sigs = find(p<.05);
% 
% for ip=1:length(sigs)
%     placement = xt_plot(sigs(ip));
%     hold on
%     plot(placement, .55-.05*hzcounter, '*', 'linew', 15, 'color', colp)
%     
% end






      ylabel({['Classification accuracy'];['SNR>median=Button press']})
        
xlabel(['Time from button press or release'])
%         axis tight
% ytix=get(gca, 'ytick');

% ylim([.8 1.3])
% axis tight
ylim([ .45 .8])
hold on; plot([0 0 ], ylim, ['k:'], 'linew', 2)
        set(gca, 'fontsize', 25)
        end
  %%
        
      
        %%
        set(gcf, 'color', 'w')
        cd(basefol)
        cd('newplots-MD')
        cd('Figures')
        cd('Predicted Classification Accuracy')
%         if ihz<=10
%         print('-dpng', [eventtype ' results ' num2str(peakfreqsare(ihz)) ' Hz excl<30, interp'])
%         else
%             print('-dpng', [eventtype ' results ' TGcase ' Hz excl<30'])
%         end


hzcounter=hzcounter+1;
end
         legend([pleg(1) ,pleg(2) ], {['20 Hz'],['40 Hz']})
         %%
         xlim([-2 2])
print('-dpng', 'classification accuracy')
end
