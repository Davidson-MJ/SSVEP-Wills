

cd(basefol)
cd('EEG')
load('GFX_LeftvsRight_TFdecomp.mat')


% first thing to plot is the topoplot, 
%which time window toplot alpha power difference?
 
 
 tstamps = [1:1:size(all_leftON_hilb,3)]/250 - 3;
 %
 topowindow = dsearchn(tstamps',[-3 3]');

 %topoplot

%plot topo of difference:
clf; figure(1);  set(gcf, 'color', 'w', 'units', 'normalized', 'position',[.5 .5 .5 .5]);
titlesare={'PFI on Left', 'PFI on right', 'Left minus right', ...
    'reapp on left', 'reapp on right', 'left minus right',...
    'Left(Disap-reap)', 'Right(Disap-Reap)'};

for itype=1%:8
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


