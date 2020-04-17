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
%%

% epoch per side (left/right disappearance).

job.epochperppant_PFI=0; % epochs as either left or right target disappearances.
job.epochperppant_PMD=0; % epochs as either left or right target disappearances.
   

job.filterandrejecttrialEEG=0; % perform basic cleaning.
    

job.participant_TFdecomp=0; % average change in T-F for Left and Right disap, per ppant.
job.participant_FilterandHilbert=0; % average change in T-F for Left and Right disap, per ppant.

job.plotppant_leveltopos=0;

job.combineAcrossppants=0; % combine TF per channel, across participants.



%% >>>>>> actual plot of hilbert results:
job.plotHILBresults =1;
%%
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

if job.epochperppant_PMD==1
   epochperppantLeftvsRight_PMD;%(basefol,dirs, allppants)
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
   plotPPANT_LvsRtopos
 
end

if job.plotHILBresults ==1
    plotALPHAHILBresults;
end


%%