%called from s3_EEG_alphalateralization  
for ifol=allppants
    %%
        cd(basefol)
cd('EEG')

cd(dirs(ifol).name)

load('ppant_PFI_Epoched_LeftRight')
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
% print('-dpng', ['Participant ' num2str(ifol) 'cleaned']);
    end