%% File one: making the raw EEG file
%based on IG's Raw_EEG but with corrections for catch timing.

cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/EEG Will SSVEP Data 01')
basefol=pwd;
pdirs=dir([pwd filesep '*_*' 'EEG']);

%%
for ippant = 1:19
    %%
    
    cd(basefol)
    % find ppant folder
    
    %change into that directory
    cd(pdirs(ippant).name);
    
% allcatchEEG = [];
    allEEG = dir([pwd filesep '*.vhdr']);
for ieeg = 1:length(allEEG)
    
    
    clear trialtrack catchinfo trialstart adjustedcatch EpochedEEG EpochEEG
    loadfile = allEEG(ieeg).name;
    
    EEG= pop_loadbv(pwd, num2str(loadfile));
    
    eegchans = EEG.data(1:64,:);
    
    Timing = (1:length(eegchans))/EEG.srate; %in seconds.
    %%
    trialend = [];
    trialstart = [];
    for ievent = 1:length(EEG.event)
        eventype = EEG.event(ievent).type;
        if strcmp(eventype, 'S 20'); %catchmarker
            trialstart=[trialstart EEG.event(ievent).latency];
        end
        if strcmp(eventype, 'S 88') %trial end
            trialend =  [trialend EEG.event(ievent).latency];
        end
    end
%     trialstart = trialend - (60*EEG.srate); %60 second trials
%     adjustedcatch = catchinfo-trialstart;
    
   
    %% Now epoch
    EpochedEEG=zeros(length(trialstart), 64, 60*EEG.srate+1); %60 seconds trials.
    
    for itrial = 1:length(trialstart)
        
        %beware negatives
        if trialstart(itrial)<0
            %supplement with NaN in duration effected
            trial= EEG.data(:, 1:trialend(itrial));
            nantrain = nan(1,abs(trialstart(itrial))+1);
            newd= zeros(64,60*EEG.srate+1);
            for ichan=1:64
                newd(ichan,1:length(nantrain)) = nantrain;
                newd(ichan,length(nantrain)+1:end) = squeeze(trial(ichan,:));
            end
            
            EpochedEEG(itrial,:,:) = newd;
        else
            trial= EEG.data(:, trialstart(itrial):trialend(itrial));
            EpochedEEG(itrial,:,:) = trial;
        end
    end
    
    
        if ieeg==1
            EpochedEEGdata =EpochedEEG;
        else
            
            %concat.
            EpochedEEGdata = cat(1, EpochedEEGdata, EpochedEEG);
            
        end
    
end
    
    %EpochedEEGdata = cat(1, data1, data2);
   cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/EEG')
   try cd(pdirs(ippant).name)
   catch
       mkdir(pdirs(ippant).name)
       cd(pdirs(ippant).name)
   end
   savename =  ['P' num2str(ippant) 'RawEEG'];
    save(savename, 'EpochedEEGdata')
    clearvars -except basefol pdirs
%     save(['P' num2str(ippant) 'Trigger_Catch_Timing'], 'allcatchEEG')
end