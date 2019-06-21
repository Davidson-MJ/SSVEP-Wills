% s3_EEG_alphalateralization
clear all, close all; clc
%% UPDATED for Wills data:
try cd('/Users/MDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')
catch
end
basefol=pwd;
cd('EEG')
dirs = dir([pwd filesep '*_*' 'EEG']);

%Epoch by number of targets involved, and direction (disapp and reapp).
% also stores duration per event.
% also stores cumularive BP (single array) per event, for next job plots.

dbstop if error


% epoch per side (left/right disappearance).

job.epochperppant_PFI=0; % epochs as either left or right target disappearances.

job.filterandrejecttrialEEG=0; % perform basic cleaning.
    


job.participant_TFdecomp=0; % average change in T-F for Left and Right disap, per ppant.
job.participant_FilterandHilbert=1; % average change in T-F for Left and Right disap, per ppant.

job.plotppant_leveltopos=1;

job.combineAcrossppants=0; % combine TF per channel, across participants.

job.plotacrossppants=0;

%load common variables:
window=[-3 3];
srate=250;
allppants=[1,2,4,6,7,9:19]; % new
getelocs;


%%%%%%% function calls:
                                % preprocessing
% epoch EEG left vs Right
if job.epochperppant_PFI==1
   epochperppantLeftvsRight(basefol,dirs, allppants)
end
if job.filterandrejecttrialEEG==1
    %filter and run ICA component rejection per ppant:
    
    %NOTE THAT EEG REDRAW DOESNT WORK INSIDE FUNCTIONS!!
    cleanLRppantEEG;   %(basefol,dirs,allppants);
end
%%                              % analysis
%crunch Time-frequency decomposition using multitaper method.
if job.participant_TFdecomp==1    
    PpantmultitaperLeftvsRight(basefol,dirs,allppants)      
end

if job.participant_FilterandHilbert==1 %filter at alpha range     
    Ppant_hilbert_LeftvsRight(basefol,dirs,allppants)          
end
 

if job.combineAcrossppants==1       
       combineLRdataacrossppants(basefol,dirs,allppants)  
end
%%                              %plotting
if job.plotppant_leveltopos==1
   

    for ifol=allppants(1:9)
    %%
        cd(basefol)
cd('EEG')

cd(dirs(ifol).name)

load('ppant_PFI_Epoched_LeftRight_cleaned')
% first thing to plot is the topoplot, 
%which time window toplot alpha power difference?
 
 
 tstamps = [1:1:size(ppant_HILB_EEG_Lefton,2)]/250 - 3;
 %
 topowindow = dsearchn(tstamps',[-.5 1]');

 %topoplot
%%
%plot topo of difference:
clf; figure(1);  set(gcf, 'color', 'w', 'units', 'normalized', 'position',[.5 .5 .5 .5]);
titlesare={'PFI on Left', 'PFI on right', 'Left minus right', ...
    'reapp on left', 'reapp on right', 'left minus right',...
    'Left(Disap-reap)', 'Right(Disap-Reap)'};
for itype=1:8
    switch itype
        case 1
       % plot PFI on left
       %mean over time points.
            d1a= squeeze(mean(ppant_HILB_EEG_Lefton(:, topowindow(1):topowindow(2)),2));
            
            % mean over ppants.
            d=d1a;
        case 2
            % plot PFI on right
            d2a= squeeze(mean(ppant_HILB_EEG_Righton(:, topowindow(1):topowindow(2)),2));
            
            d=d2a;
        case 3
            % difference between them.
            d=d1a-d2a;
            
        case 4
            d1b= squeeze(mean(ppant_HILB_EEG_Leftoff(:, topowindow(1):topowindow(2)),2));
            d=d1b;
        case 5
            d2b= squeeze(mean(ppant_HILB_EEG_Rightoff(:, topowindow(1):topowindow(2)),2));
            d=d2b;
        case 6
            d=d1b-d2b;
            
        case 7 %left disap - reap
            d=d1a-d1b;
           
        case 8
            d=d2a-d2b; % right PFI - right reapp.
           
    end
            
subplot(3,3,itype);
topoplot(d, elocs(1:64)); c=colorbar;
ylabel(c,'alpha power (uV)', 'fontsize', 15)
title([titlesare{itype} ', ' num2str(tstamps(topowindow(1))) ',' num2str((tstamps(topowindow(2))))]) 
% caxis([-10 10])
end
set(gcf, 'color', 'w')
cd(basefol)
cd('Figures')
cd('Participant L vs R Target PFI');
print('-dpng', ['Participant ' num2str(ifol) 'cleaned']);
    end
end



if job.plotacrossppants==1
cd(basefol)
cd('EEG')
load('GFX_LeftvsRight_TFdecomp.mat')


% first thing to plot is the topoplot, 
%which time window toplot alpha power difference?
 
 
 tstamps = [1:1:size(datatmp,3)]/250 - 3;
 %
 topowindow = dsearchn(tstamps',[-.5 1]');

 %topoplot

%plot topo of difference:
clf; figure(1);  set(gcf, 'color', 'w', 'units', 'normalized', 'position',[.5 .5 .5 .5]);
titlesare={'PFI on Left', 'PFI on right', 'Left minus right', ...
    'reapp on left', 'reapp on right', 'left minus right',...
    'Left(Disap-reap)', 'Right(Disap-Reap)'};
for itype=1:8
    switch itype
        case 1
       % plot PFI on left
       %mean over time points.
            d1a= squeeze(mean(all_leftON_hilb(:,:, topowindow(1):topowindow(2)),3));
            
            % mean over ppants.
            d=mean(d1a,1);
        case 2
            % plot PFI on right
            d2a= squeeze(mean(all_rightON_hilb(:,:, topowindow(1):topowindow(2)),3));
            
            d=mean(d2a,1);
        case 3
            % difference between them.
            d=d1a-d2a;
            d=mean(d,1);
        case 4
            d1b= squeeze(mean(all_leftOFF_hilb(:,:, topowindow(1):topowindow(2)),3));
            d=mean(d1b,1);
        case 5
            d2b= squeeze(mean(all_rightOFF_hilb(:,:, topowindow(1):topowindow(2)),3));
            d=mean(d2b,1);
        case 6
            d=d1b-d2b;
            d=mean(d,1);
        case 7 %left disap - reap
            d=d1a-d1b;
            d=mean(d,1);
        case 8
            d=d2a-d2b; % right PFI - right reapp.
            d=mean(d,1);
    end
            
subplot(3,3,itype);
topoplot(d, elocs(1:64)); c=colorbar;
ylabel(c,'alpha power (uV)', 'fontsize', 15)
title([titlesare{itype} ', ' num2str(tstamps(topowindow(1))) ',' num2str((tstamps(topowindow(2))))]) 
% caxis([-10 10])
end
set(gcf, 'color', 'w')
% suptitle('1')
     
% for ippant=1:size(Diffdata,1)
%     figure(ippant); 
%     %subtract one from the other
%     subplot(2,2,1);        
% topoData = mean(squeeze(mean(Diffdata(ippant,:,[topowindow(1):topowindow(2)]),3)),1);
% topoplot(topoData, elocs(1:64)); shg
% end
% end
%% % %% plot trace?
figure(1); clf
%  
 
        set(gcf, 'units', 'normalized', 'position', [0 0 1 1 ]);
% TF_outgoing = squeeze(mean(all_leftON,1));
% %%        
clf
for ichan=1:64
            subplot(8,8,ichan);
            colsare={'r', 'b', 'k', 'm'};
for id=[1,2]%[1,3]
    switch id
        case 1
            datahilb=all_leftON_hilb;
        case 2
            datahilb=all_leftOFF_hilb;
        case 3
            datahilb=all_rightON_hilb;
        case 4
            datahilb=all_rightOFF_hilb;
    end
    hold on;
    pdata=smooth(squeeze(mean(datahilb(:,ichan,:),1)),100)';
    plot(tstamps,pdata, colsare{id});
    title([elocs(ichan).labels]);
    ylim([-2 2])
    xlim([-.200 .500]);
    
end
end
%%
%now show topo of diff.
 
tstamps = [1:1:size(datahilb,3)]/250 - 3;
 %
topowindow = dsearchn(tstamps',[0 1]');


%plot topo of difference:
figure(2); clf
Leftgone= squeeze(mean(all_leftON_hilb,1));
subplot(221);
d1= squeeze(mean(Leftgone(:, topowindow(1):topowindow(2)),2));

topoplot(d1, elocs(1:64)); colorbar;
title(['Left TG disap, ' num2str(tstamps(topowindow(1))) ',' num2str((tstamps(topowindow(2))))]) 
subplot(224)
d2= squeeze(mean(ppant_HILB_EEG_Righton(:, topowindow(1):topowindow(2)),2));
title(['Right TG disap, ' num2str(tstamps(topowindow(1))) ',' num2str((tstamps(topowindow(2))))]) 
topoplot(d2, elocs(1:64)); colorbar;

%diff between
subplot(222)
Dd=d1-d2;
topoplot(Dd, elocs(1:64)); colorbar;




end
%%