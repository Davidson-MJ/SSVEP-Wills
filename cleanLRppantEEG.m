% function cleanLRppantEEG(basefol,dirs, allppants)

cd(basefol)        
    cd('EEG')
    [ALLEEG,EEG,CURRENTSET]=eeglab; % open eeglab in background.
    getelocs;
    for ifol=allppants(9):allppants(end);
        cd(basefol)
        cd('EEG')
        cd(dirs(ifol).name)
        
        load('ppant_PFI_Epoched_LeftRight');
        
        %for each participant, we will actually combine all data together
        %(more efficient, better for ICA).
        
        datatmp= cat(1,  ppant_EEG_PFI_Lefton,ppant_EEG_PFI_Leftoff,ppant_EEG_PFI_Righton,ppant_EEG_PFI_Rightoff);
        %% to reshape later, record trial counts.
        trialindices = [ones(1,size(ppant_EEG_PFI_Lefton,1)), ...
            ones(1,size(ppant_EEG_PFI_Leftoff,1))+1,...
            ones(1,size(ppant_EEG_PFI_Righton,1))+2,...
            ones(1,size(ppant_EEG_PFI_Rightoff,1))+3];
        % now can find to extract trial #s.
        %e.g. left on = find(trialindices==1)
        
        
        p_data = permute(datatmp,[2,3,1]);
        
        % before converting to EEGLab (for ICA), need to demean:
        
        for ichan=1:64
            for itrial=1:size(p_data,3)
                
                %demean
                trdata= squeeze(p_data(ichan,:,itrial));
                tmp = trdata- mean(trdata);
                
                
                
                %and filter from broadband
                fdata = eegfilt(tmp, 250, 1,40, 1501,100,0,'fir1');
                
                p_data(ichan,:,itrial)=fdata;
            end
        end
        
        
        %% remove bad channels
        
        %inputs:
        rejtrialcounts=zeros(64,1);
        SDchancounts=zeros(64,1);
        % reject trials with large voltage spikes, and or
        % SD greater than 5* channel average.
        
        for ichan= 1:64
            rejcount=0;
            sdpertrial=[];
            for itrial = 1:size(p_data,3)
                
                datain = squeeze(p_data(ichan,:,itrial));
                
                if range(datain)>1000
                    rejcount=rejcount+1;
                end
                sdpertrial=[sdpertrial, std(datain)];
            end
            %how many trials rejected?
            rejtrialcounts(ichan)=rejcount;
            SDchancounts(ichan)=mean(sdpertrial);
        end
        %%
        
        %identify (iteratively) outliers based on SD:
        TF= isoutlier(SDchancounts', 'gesd', 'MaxNumOutliers', 10);
        rejchan1= find(TF);
        %else identify channels with extreme SD.
        MallchanSD= std(SDchancounts);
        rejchan2 = find(SDchancounts>(3*MallchanSD));
        
        %identify all bad channels.
        allchans_torej = union(rejchan1,rejchan2');
        
        %%
        %convert matlab now to EEG.
        
        EEG= pop_importdata('setname', 'tmp','data', p_data, 'dataformat', 'matlab', 'srate', 250,'chanlocs', elocs(1:64),'nbchan', 64);
        %save data set here.
        EEG=pop_saveset(EEG,'filename', ['participant ' num2str(ifol) 'EEGset'], ...
            'filepath', num2str(pwd));
        %store also (in workspace)
        [ALLEEG, EEG]=eeg_store(ALLEEG,EEG,CURRENTSET);
        %load to eeglab:
        pop_loadset('filename',['participant ' num2str(ifol) 'EEGset.set'], 'filepath',num2str(pwd));
        eeglab redraw
     %%
        %remove bad channels (interpolate).
        EEG=pop_interp(EEG,allchans_torej, 'spherical');
        
        %% run ica
        EEG=pop_runica(EEG, 'icatype', 'runica');
        %% now  we have ICA components, import to EEG lab, reject bad components and save
       
%         [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
        eeglab redraw
%
        %select components for rejection (manually), drawing top 15 at a time
        EEG=pop_selectcomps(EEG,1:15);
        
        %breaks into GUI.
        %%
        yn=input(['ready to continue ? ']);
        
        %automate rejection:       
        rejcomps = find(EEG.reject.gcompreject);    
        
        EEG= pop_subcomp(EEG,rejcomps);
        
        %update
        eeglab redraw;
        
        %% rerun ICA (free from initial bad datapoints).
        disp('running ICA2...')
        EEG=pop_runica(EEG, 'icatype', 'runica');
        %throw new to eeglab.
        eeglab redraw
        EEG=pop_selectcomps(EEG,1:15);
         %breaks into GUI.
        %%
        yn=input(['ready to continue ? ']);
        
        eeglab redraw;
        
        %automate rejection:       
        rejcomps = find(EEG.reject.gcompreject)  
        %automate rejection:
            
        EEG= pop_subcomp(EEG,rejcomps);
        
        %save final set
        EEG=pop_saveset(EEG,'filename', ['participant ' num2str(ifol) 'EEGset'], ...
            'filepath', num2str(pwd));
        
        %Now reshape and save matlab vars for working (next script).
        alltrials = size(EEG.data,3);
        
        
        for id=1:4
            
            dataind = find(trialindices==id);
            
            dataE=EEG.data(:,:,dataind);
            
            switch id
                case 1
                    ppant_EEG_PFI_Lefton=dataE;
                case 2
                    ppant_EEG_PFI_Leftoff=dataE;
                case 3
                    ppant_EEG_PFI_Righton=dataE;
                case 4
                    ppant_EEG_PFI_Rightoff=dataE;
            end
        
        end
        
        %save new clean data:
        save('ppant_PFI_Epoched_LeftRight_cleaned', 'ppant_EEG_PFI_Lefton', ...
            'ppant_EEG_PFI_Leftoff', 'ppant_EEG_PFI_Righton',...
            'ppant_EEG_PFI_Rightoff')
    end