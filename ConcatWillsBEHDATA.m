clear all
%gather folder directories.
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/Behavioural WILL SSVEP Data')
allppantfols = dir([pwd filesep '*_*']);
%%
%preallocate variables (using same var names as later scripts (ithink)
AllTrialsAllPPTopLeft=[];
AllTrialsAllPPTopRight=AllTrialsAllPPTopLeft;
AllTrialsAllPPBottomLeft=AllTrialsAllPPTopLeft;
AllTrialsAllPPBottomRight=AllTrialsAllPPTopLeft;

%cycle through and append data.
basefol=pwd;
for ifol = 1:length(allppantfols)
    cd(basefol)
    try
    cd(allppantfols(ifol).name)
    
    
    
    %first we need to rename so that the order is corrected!!!
    for filen=1:9
        try
        movefile(['AllTarget15Hz_Trial' num2str(filen) '_data.mat'], ['AllTarget15Hz_Trial0' num2str(filen) '_data.mat'])
        catch
        end
    end
    %how many trials?
    trdata = dir([pwd filesep 'AllTarget*']);
    
    
    
    
    for itrial = 1:length(trdata)
        
        %first clear previous variables, so we overwrite.
        clearvars data*
        
        load(trdata(itrial).name)
        %now append each.
        %please note that this is a less than ideal way to append and store
        %data, but i'm using it to ease the transition to using Irene's
        %analysis pipeline (which uses these variable formats).
        AllTrialsAllPPTopLeft=[AllTrialsAllPPTopLeft;dataTopLeft];
        AllTrialsAllPPTopRight=[AllTrialsAllPPTopRight; dataTopRight];
        AllTrialsAllPPBottomLeft=[AllTrialsAllPPBottomLeft; dataBottomLeft];
        AllTrialsAllPPBottomRight=[AllTrialsAllPPBottomRight; dataBottomRight];
    end
    catch
        disp(['skipped fol: ' num2str(allppantfols(ifol).name)])
    end
end
%%
cd(basefol)
cd ../
cd('AA_ANALYZED DATA')
cd('Behaviour')
%%
save('Raw_Behavioral_AllPP', 'AllTrialsAllPPBottomLeft',...
    'AllTrialsAllPPBottomRight',...
    'AllTrialsAllPPTopLeft',...
    'AllTrialsAllPPTopRight');