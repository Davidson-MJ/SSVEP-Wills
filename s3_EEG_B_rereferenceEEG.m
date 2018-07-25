% s3_B_rereferenceEEG
%new rereference script.
addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/EEG')
basefol=pwd;

%%
getelocs
dbstop if error

allppants=[1,2,4,6,9:16,18]; %
dirs = dir([pwd filesep '*_*' 'EEG']);

%%
for ifol = allppants
    cd(basefol)
    cd(dirs(ifol).name)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         f=dir();
        
%         f={f.name};
        %continue if the output file doesn't exist
%         if isempty(find(strcmp(f,['P' num2str(ifol) '_autopreprocd.mat'])))
            
            %don't delete this file
%             n=find(strcmp(f,['P' num2str(ifol) 'RawEEG.mat']));
            
            
%             f{n}=[];
            
            %clear unnecessary vars/ scripts. we are starting again
            
%                     for k=1:numel(f);
%                         delete([ f{k}])
%                     end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %load Epoched EEG;
            load(['P' num2str(ifol) 'RawEEG'])
            
            %first trial is a dud.
            %
            
%             outgoing=zeros(48,64,60001);
            outgoing=zeros(48,64,15001);
            
            for itrial=1:48 %since first was a practice.
                if size(EpochedEEGdata,1)==49
                realtrial=itrial+1;
                elseif  size(EpochedEEGdata,1)==48
                    realtrial=itrial;
                end
                tmp=squeeze(EpochedEEGdata(realtrial,1:64,:));
%                 
%                 %reref to average of chans.
                re_tmp=reref(tmp);
%                 
%                 %demean, detrend, and down sample to 250Hz
%                 for ichan = 1:64
%                     chandata = squeeze(re_tmp(ichan,:));
%                     
%                     %demean
%                     chandata= chandata-mean(chandata);
%                     
%                     %detrend
%                     tmp= detrend(chandata,'linear');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
%                     %downsample
                    re_tmp = downsample(re_tmp',4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     outgoing(itrial,ichan,:) = dstmp;
%                 end
                        outgoing(itrial,:,:) = re_tmp';
            end
            
            %perform laplacian
            outgoing = permute(outgoing, [2 3 1]);
            EEG_preprocd=laplacian_perrinX(outgoing, [elocs(1:64).X], [elocs(1:64).Y], [elocs(1:64).Z]);
            % EEG_preprocd=outgoing;
            
            
            save(['P' num2str(ifol) '_autopreprocd'], 'EEG_preprocd')
            %accidentally saved all as P1 on first run through.
            %save(['P1RereferencedEEG_new'], 'RerefEEG')
            disp(['Fin ppant ' num2str(ifol)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %         end
        
end
