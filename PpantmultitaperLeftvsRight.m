function PpantmultitaperLeftvsRight(basefol,dirs,allppants)
 %take multitaper spectrogram, zero-pad.
    params.tapers = [1,1];
    params.Fs=250;
    params.pad=2;
    params.fpass=[0 45];
    movingwin=[.5,.2];
    params.trialave=1;
    
    
    
    % prep for data outgoing:
    % nfreqs:
    N=1501;
    Nwin=round(params.Fs*movingwin(1)); % number of samples in window
    Nstep=round(movingwin(2)*params.Fs); % number of samples to step through
    nfft=max(2^(nextpow2(Nwin)+params.pad),Nwin);
    f=getfgrid(params.Fs,nfft,params.fpass); Nf=length(f);
    
    % and time steps
    winstart=1:Nstep:N-Nwin+1;
    nw=length(winstart);
    winmid=winstart+round(Nwin/2);
    t=winmid/params.Fs;
    Nt=length(t);
    
        for ifol=allppants
           cd(basefol)
           cd('EEG')
        cd(dirs(ifol).name)
        
        load('ppant_PFI_Epoched_LeftRight');
        
        
        for idatatype=1:4 
            switch idatatype
                case 1
                    datatmp=ppant_EEG_PFI_Lefton;
                case 2
                    datatmp=ppant_EEG_PFI_Leftoff;
                case 3                    
                    datatmp=ppant_EEG_PFI_Righton;
                case 4
                    
                    datatmp=ppant_EEG_PFI_Rightoff;
            end
        
            
            %for each type, prep the outgoin:
            TF_outgoing = zeros(64, Nf,Nt); %chans, freqs, timep.
        
        for ichan=1:64
                        
            chan_data= squeeze(datatmp(:,ichan,:));
            
            %%
            
            [s,t,f]=mtspecgramc(chan_data', movingwin, params);
            %
            s=log(s);
            
            tstamp=t-3;
            %sanity check
            % for each trial, remove the baseline first 1000 ms before 
            [baseend]=dsearchn(tstamp',-2');
           
            brem_specgram=zeros(size(s));
            
            for ifreq=1:size(s,2)
               
              
               lp= s(:,ifreq)';
               brem=mean(lp(1:baseend));
               brep= repmat(brem,[1, length(lp)]);
               
               brem_specgram(:,ifreq)= lp-brep;
            end
           TF_outgoing(ichan,:,:)=brem_specgram';
           
           
%             % check baseline subtraction:
%             subplot(2,1,1);
%             imagesc(t,f,(s)');
%             axis xy; colorbar;
%             shg; ylim([0 20]);
%             subplot(212)
%             imagesc(t,f,(brem_specgram)');
%             axis xy; colorbar;
%             shg; ylim([0 20]);
%         
        end % trial
        
        %% display TF across chans:
%         figure(10);
%         set(gcf, 'units', 'normalized', 'position', [0 0 1 1 ]);
%         for ichan=1:64
%             subplot(8,8,ichan);
%             imagesc(tstamp,f,squeeze(TF_outgoing(ichan,:,:)));
%             title(elocs(ichan).labels);
%             axis xy; ylim([0 20])
%             caxis([-1 1])
%         end
        
        %%
        
        switch idatatype
            case 1
                ppant_TF_EEG_Lefton=TF_outgoing;
            case 2
                ppant_TF_EEG_Leftoff=TF_outgoing;
            case 3
                ppant_TF_EEG_Righton=TF_outgoing;
            case 4
                ppant_TF_EEG_Rightoff=TF_outgoing;
        end
        
        end % idatatype
    
        
        %save in ppant fol:
        save('ppant_PFI_Epoched_LeftRight',...
            'ppant_TF_EEG_Rightoff','ppant_TF_EEG_Righton',...
            'ppant_TF_EEG_Leftoff','ppant_TF_EEG_Lefton', 'tstamp','f','-append');
        end