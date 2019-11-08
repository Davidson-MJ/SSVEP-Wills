function combineLRdataacrossppants(basefol,dirs,allppants)

[all_leftON,all_leftOFF,all_rightON,all_rightOFF]=deal(zeros(length(allppants), 64, size(ppant_TF_EEG_Lefton,2),size(ppant_TF_EEG_Lefton,3)));
[all_leftON_hilb,all_leftOFF_hilb,all_rightON_hilb,all_rightOFF_hilb]=deal(zeros(length(allppants), 64,1501));

counter=1;
for ifol=1:allppants
    
    cd(basefol)
    cd('EEG')
    cd(dirs(ifol).name)
    
    load('ppant_PFI_Epoched_LeftRight');
    all_leftON(counter,:,:,:)=ppant_TF_EEG_Lefton;
    all_leftOFF(counter,:,:,:)=ppant_TF_EEG_Leftoff;
    all_rightON(counter,:,:,:)=ppant_TF_EEG_Righton;
    all_rightOFF(counter,:,:,:)=ppant_TF_EEG_Rightoff;
    
    
    all_leftON_hilb(counter,:,:,:)=ppant_HILB_EEG_Lefton;
    all_leftOFF_hilb(counter,:,:,:)=ppant_HILB_EEG_Leftoff;
    all_rightON_hilb(counter,:,:,:)=ppant_HILB_EEG_Righton;
    all_rightOFF_hilb(counter,:,:,:)=ppant_HILB_EEG_Rightoff;
    
    counter=counter+1;
end

cd(basefol)
cd('EEG');
save('GFX_LeftvsRight_TFdecomp', 'all_leftOFF', 'all_leftON',...
    'all_rightOFF','all_rightON','tstamp', 'f',...;
    'all_leftON_hilb', 'all_leftOFF_hilb',...
    'all_rightON_hilb','all_rightOFF_hilb');