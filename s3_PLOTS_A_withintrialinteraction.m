
%% first check each ppants 48 trials, if necessary. Most/all probably have clear interaction

clearvars -except basefol allppants
addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')


basefol=pwd;
cd('Behaviour')
%%
load('MD_AllBP_perppant.mat')


dbstop if error
cd ../
cd('EEG');
pdirs = dir([pwd filesep '*EEG']);


peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; 
%
% % check all trials
xt=[1:3600]/60+1;

%%

catchcol= 'r';
fontsize=25;
%
%ippant=3, itrial = 2 in the paper
%also load the example SNR.
useshuffled=0;
for ippant=1%allppants(15)
    
    cd(basefol)
    cd('EEG')    
    %load whole trial RESS
    cd(pdirs(ippant).name)
    
    load('P_wholetrialRESS'); %
    cd('/Users/MattDavidson/Desktop')
    cd('ppant trial-SNR correlation')
    %%
    for    itrial=10%:24%
        
        figure(1); clf
        subplot(311)
        
        hold on;
        colormap('jet')
        cmap = [ .9 .9 .9; 0 .9 0 ];
%         cmap = [ .9 .9 .9; 0 .9 0 ; 0 .6 0; 0 .3 0; 0 .1 0];
        colormap(cmap)
        plotall = zeros(5,3600);
        if useshuffled==1
            plotall = squeeze(allRandomAllPP(ippant,itrial,:,:));
        else
            plottmp=ppantTrialDatawithDetails(ippant).AllBPData;
            plot1=squeeze(plottmp(itrial,:,:));
            plott=sum(plot1,1);
plotall(1:4,:) = plot1;
% plotall(5,:) = plott;
%             plotall=plott;
        end
        
        imagesc(xt, 1:1:size(plotall,1), plotall);
     
        set(gca, 'yticklabel', [])
        % set(gca, 'yticklabel', yt)
        hold on
      
            % title(['Participant ' num2str(ippant) ', Trial ' num2str(itrial) ', Button Press Response'], 'fontsize', fontsize)
            catchstart = ppantTrialDatawithDetails(ippant).TrialDetails(itrial).Catchstart_frames/60;
            catchend = ppantTrialDatawithDetails(ippant).TrialDetails(itrial).Catchend_frames/60;
            xvec = catchstart:.1:catchend;
       
            hzlocsthistrial=[ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TL_Freq,...
                ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TR_Freq,...
                ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BL_Freq,...
                ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BR_Freq];
            
            %catchlocs this trial:
            catchlocsthistrial=[ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TL_Catch,...
                ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TR_Catch,...
                ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BL_Catch,...
                ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BR_Catch];
            
        % ylabel('Button Press location');
        xlabel('Time(s)')
        catchexist=0; % unless plotted below.
        try
        if useshuffled==0
            % box around catch
            if ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TL_Catch==1;
                %drawbox around catch times
                plot([catchstart catchstart], [.5 1.5], [catchcol], 'linewidth',5);%box end
                plot([catchend catchend], [0.5 1.5], [catchcol], 'linewidth',5);
                plot([catchstart catchend], [.5 .5],  [catchcol], 'linewidth',5);
                l=plot([catchstart catchend], [1.5 1.5],  [catchcol], 'linewidth',5);
                catchexist=1;
            end
            %
            if ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TR_Catch==1;
                l=plot([catchstart catchstart], [1.5 2.5], [catchcol], 'linewidth',5);
                l=plot([catchend catchend], [1.5 2.5], [catchcol], 'linewidth',5);
                plot([catchstart catchend], [1.5 1.5],  [catchcol], 'linewidth',5);
                l=plot([catchstart catchend], [2.5 2.5],  [catchcol], 'linewidth',5);
                catchexist=1;
            end
            %
            if ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BL_Catch==1;
                plot([catchstart catchstart], [2.5 3.5], [catchcol], 'linewidth',5);
                plot([catchend catchend], [2.5 3.5], [catchcol], 'linewidth',5);
                plot([catchstart catchend], [2.5 2.5],  [catchcol], 'linewidth',5);
                l=plot([catchstart catchend], [3.5 3.5],  [catchcol], 'linewidth',5);
                catchexist=1;
            end
            
            if ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BR_Catch==1;
                plot([catchstart catchstart], [3.5 4.5], [catchcol], 'linewidth',5);
                plot([catchstart catchend], [3.5 3.5],  [catchcol], 'linewidth',5);
                %box end
                plot([catchend catchend], [3.5 4.5], [catchcol], 'linewidth',5);
                l=plot([catchstart catchend], [4.5 4.5],  [catchcol], 'linewidth',5);
                catchexist=1;
            end
%             lg=legend(l, 'Catch Period');
%             set(lg, 'fontsize', fontsize)
        end
        %
        catch
        end
        set(gca, 'ydir', 'reverse')
        axis tight
        xlim([tgrm_wholetrial(1) tgrm_wholetrial(end)]) %adjust behavioural window to match SNR.
        plot(xlim, [4.5 4.5], 'b')
        
        %%
        set(gca, 'xtick', 0:5:60, 'fontsize', fontsize)
        
%         set(gca, 'ytick', 1:5, 'yticklabel', [{['TL(' num2str(hzlocsthistrial(1)) ' Hz)']},...
%             {['TR(' num2str(hzlocsthistrial(2)) ' Hz)']}, {['BL(' num2str(hzlocsthistrial(3)) ' Hz)']},...
%             {['BR(' num2str(hzlocsthistrial(4)) ' Hz)']}, 'Total'])
        set(gca, 'ytick', 1:5, 'yticklabel', [{['Top Left']},...
            {['Top Right']}, {['Bottom Left']},...
            {['Bottom Right']}, 'Total'])
        %%
        % end
        set(gcf, 'color', 'w')
                title(['Participant ' num2str(ippant) ', trial ' num2str(itrial)])

        hold on;
%         ylim([0 8])
        %
        
        %% load this trial SNR and compare
        figure(1);  hold on
%         clf
        pt=[];
        ptlabel={};
        freqsare= peakfreqsare;
        %%
        peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; 
        cols = {'b', 'm', 'b', 'm', 'm', 'g', 'k', 'k', 'k'};
        linews = [1,1,2,2,3,3,1,1,1]*2;
%         clf
        for ihz = [1,2,7]
            
            subplot(3,1,2:3)
            xlim([1 60])
            
            hold on % plot catch start and end
            
            % plot both harms?
            
            plotme=squeeze(snr_wholetrRESS(ihz,itrial,:))';
            
            plw=3;
            
            
            pt(ihz)=plot(tgrm_wholetrial, smooth(plotme), 'linew', linews(ihz), 'color', cols{ihz});
            ptlabel{ihz}= [num2str(peakfreqsare(ihz)) 'Hz'];
            
            
            
        end
        axis tight
            %%
% now plot catch locations.
         if catchexist==1
             if ihz>=5
                 hold on
                 pch=patch([catchstart catchstart catchend catchend], [-1 7 7 -1], catchcol);
                 pch.FaceAlpha=.1;
             else % place on target plot
                 tindx = find(hzlocsthistrial==usehz);
                 
                 if catchlocsthistrial(tindx)==1
                     pch=patch([catchstart catchstart catchend catchend], [mean(plotme(:))-4 mean(plotme(:))+4 mean(plotme(:))+4 mean(plotme(:))-4], catchcol);
                     pch.FaceAlpha=.1;
                 end
                 
                 
             end
             %and add legend / yticklabels.
             if ihz==4
                 
                 set(gca, 'ytick', [-40:10:-10], 'yticklabels',  {ptlabel{4}, ptlabel{3},ptlabel{2},ptlabel{1}}) %reverse order.
                 xlabel('Time(s)')
                 ylabel('Target SSVEPs')
             elseif ihz==6
                 %              lg=legend([pt(5), pt(6)], {ptlabel{5}, ptlabel{6}});
                 %              set(lg, 'fontsize', fontsize)
                 % set(gca, 'ytick', [-60:10:-50], 'yticklabels',  {ptlabel{6}, ptlabel{5}})
                 ylabel({['Background SSVEPs,'];['log(SNR)']})
                 % legend ('1F', '2F')
                 lg=legend([pt(5), pt(6)], {'1F', '2F'});
                 set(lg, 'location', 'SouthEast')
             end
         end

            
        set(gca, 'fontsize', fontsize, 'xtick', 0:5:60)

%         xlim([tgrm_wholetrial(1) tgrm_wholetrial(end)]) %adjust behavioural window to match SNR.
        
        
%         print('-dpng', ['ppant ' num2str(ippant) ', trial ' num2str(itrial)  ])% ', smoothed 1s'])
    end
end
%%
figure(2);
  imagesc(xt, 1:4, squeeze(plotall(1:4,:))); 
  
 