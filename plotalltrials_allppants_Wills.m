dbstop if error
%plot all trials all ppants for easy quick reference and sanity checks.
% try cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/Behaviour')
% catch cd('/Users/MatthewDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/Behaviour')
% end

load('MD_AllBP_perppant.mat')
load('Catchperformance.mat', 'catchStruct')

%
cd ../
cd('Figures')
cd('All Participant Trials')
%%
fontsize=10;

linewidth=5;
dbstop if error
%which colour for the figure?
catchcol= 'r';

%
useshuffled=1;
 
%
xt=0:1/60:60;
dbstop if error
for ippant = 15%:29
    % check all trials
    
    figure(1); clf
    
    for itrial=31%:24
%         subplot(8,3, itrial)
        
%         colormap('Gray(4)')
%         colormap('Viridis')
if useshuffled~=1
        colormap([0.95, 0.95, 0.95; 0, 0, 1])
else
        colormap([0.95, 0.95, 0.95; .5, .5, .5])
end
    rowis = (itrial+1) + 24*(ippant-1);
    
        if useshuffled==1
            plotall = squeeze(allRandomAllPP(2,itrial,1:4,:));
        else
        plot1=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BPData;
        plotall = [plot1];%; sum(plot1,1)];
        end
        imagesc(xt, 1:1:4, plotall);
        imagesc(1:length(plotall)/60, 1:1:4, plotall);
        hold on
%         plot([1 length(plot1)], [4.5 4.5], ['r'], 'linewidth',1);
        catchstart = ppantTrialDatawithDetails(ippant).TrialDetails(itrial).Catchstart_frames/60;
        catchend = ppantTrialDatawithDetails(ippant).TrialDetails(itrial).Catchend_frames/60;
        % catch locations?
        
%         set(gca, 'fontsize', 40)
        %
        xvec = catchstart:.1:catchend;
%         c=colorbar;
%         ylabel(c, 'Total Buttons Pressed')
%         set(c,'Ytick', [ 0, 1, 2, 3])
%         caxis([0 3])
        if str2num(catchStruct{rowis,6})==1
            title('Reject(?)', 'Color', 'r');
        elseif isnan(str2num(catchStruct{rowis,6}))
            title('Ignore(?)', 'Color', 'b');
        else
%         title(['Participant ' num2str(ippant) ', Trial ' num2str(itrial) ], 'fontsize', fontsize)
        end
        
        % box around catch
        if useshuffled~=1
        if ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TL_Catch==1
            %drawbox around catch times
            plot([catchstart catchstart], [0.5 1.5], [catchcol], 'linewidth',linewidth);%box end
            plot([catchend catchend], [0.5 1.5], [catchcol], 'linewidth',linewidth);
            plot([catchstart catchend], [0.5 0.5],  [catchcol], 'linewidth',linewidth);
            l=plot([catchstart catchend], [1.5 1.5],  [catchcol], 'linewidth',linewidth);
        end
        %
        if ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TR_Catch==1
            l=plot([catchstart catchstart], [1.5 2.5], [catchcol], 'linewidth',linewidth);
            l=plot([catchend catchend], [1.5 2.5], [catchcol], 'linewidth',linewidth);
            plot([catchstart catchend], [1.5 1.5], [catchcol], 'linewidth',linewidth);
            l=plot([catchstart catchend], [2.5 2.5],  [catchcol], 'linewidth',linewidth);
        end
        %
        if ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BL_Catch==1
            plot([catchstart catchstart], [2.5 3.5], [catchcol], 'linewidth',linewidth);
            plot([catchend catchend], [2.5 3.5], [catchcol], 'linewidth',linewidth);
            plot([catchstart catchend], [2.5 2.5],  [catchcol], 'linewidth',linewidth);
            l=plot([catchstart catchend], [3.5 3.5],  [catchcol], 'linewidth',linewidth);
        end
        
        if ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BR_Catch==1
            plot([catchstart catchstart], [3.5 4.5], [catchcol], 'linewidth',linewidth);
            plot([catchstart catchend], [3.5 3.5],  [catchcol], 'linewidth',linewidth);
            %box end
            plot([catchend catchend], [3.5 4.5], [catchcol], 'linewidth',linewidth);
            l=plot([catchstart catchend], [4.5 4.5],  [catchcol], 'linewidth',linewidth);                        
        end
%         lg=legend(l, 'catch period');
%         set(lg, 'fontsize', 35, 'Location', 'northeastOutside')
    
    
                set(gca, 'ytick', 1:4, 'yticklabel', {'Top Left', 'Top Right', 'Bottom Left', 'Bottom Right'})
        else
            
            set(gca, 'ytick', 1:4, 'yticklabel', {'Top Left \bf(Trial 1)' , 'Top Right \bf(Trial 18)', 'Bottom Left \bf(Trial 7)', 'Bottom Right \bf(Trial 21)'})
        end
        
    end
    set(gcf, 'color', 'w');
    
%       print('-dpng', ['All Trials Participant ' num2str(ippant)])
end
%
    set(gca, 'fontsize', 25)
    xlabel(['seconds'])
  %% collect number of trials rejected. make sure 'calcmissed catches x2' has been run in s1
  prej = [];
for ippant = 1:19
    %
    itrials = [2:49] + (ippant-1)*48;
    trialsare=[];
    for i=1:48
    trialsare =[trialsare, str2num(catchStruct{itrials(i),6})];
    end    
    %
    prej(ippant) = length(find(trialsare>0));
    
end
   
%
retainedtrials = repmat(48, [1 length(prej)]);

 allppants=[1,2,4,6,7,9:19]; % new
 adjustedrej= retainedtrials -prej;
   catchfailpercent = 1- adjustedrej./retainedtrials;
   
   plotgoodppants = catchfailpercent(allppants);
   
   histogram(plotgoodppants,6);
   xlabel('proportion of trials marked for rejection')
ylabel('# of participants')
set(gca, 'fontsize', 15*1.5)
set(gcf, 'color', 'w')
% ylabel('# of subjects')
xlim([0 .1])
shg
    %%
    percrej=nansum(prej)/( 48);
    %%
    iqr(prej)
    
    
    