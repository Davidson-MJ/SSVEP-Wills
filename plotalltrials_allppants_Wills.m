dbstop if error
%plot all trials all ppants for easy quick reference and sanity checks.
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/Behaviour')

load('MD_AllBP_perppant.mat')
load('Catchperformance.mat', 'catchStruct')
fontsize=10;
linewidth=1;
%
cd ../
cd('Figures')
cd('All Participant Trials')
%%
dbstop if error
for ippant = 1:20
    % check all trials
    
    figure(1); clf
    % each ppant has 48 trials.
    for ifig=1:2
        %different indexing based on data type:
        thesetrialsCATCHstruct = 1+(1:24)+ (24*(ifig-1)) + 48*(ippant-1);
        
        thesetrialsPPANTstruct= (1:24)+ (24*(ifig-1));
    for itrial=1:24
        subplot(8,3, itrial)
        
        colormap('Gray(4)')
    
    
        usetrial = thesetrialsPPANTstruct(itrial);
        plot1=ppantTrialDatawithDetails(ippant).TrialDetails(usetrial).BPData;
        plotall = [plot1];%; sum(plot1,1)];
        imagesc(1:length(plot1)/60, 1:1:4, plotall);
        hold on
%         plot([1 length(plot1)], [4.5 4.5], ['r'], 'linewidth',1);
        catchstart = ppantTrialDatawithDetails(ippant).TrialDetails(usetrial).Catchstart_frames/60;
        catchend = ppantTrialDatawithDetails(ippant).TrialDetails(usetrial).Catchend_frames/60;
        % catch locations?
        
%         set(gca, 'fontsize', 40)
        %
        xvec = (catchstart:.1:catchend);
%         c=colorbar;
%         ylabel(c, 'Total Buttons Pressed')
%         set(c,'Ytick', [ 0, 1, 2, 3])
        caxis([0 3])
        rowis=thesetrialsCATCHstruct(itrial);
        if str2num(catchStruct{rowis,6})==1
            title('Reject(?)', 'Color', 'r');
        elseif isnan(str2num(catchStruct{rowis,6}))
            title('Ignore(?)', 'Color', 'b');
        else
        title(['Participant ' num2str(ippant) ', Trial ' num2str(usetrial) ], 'fontsize', fontsize)
        end
        
        % box around catch
        if ppantTrialDatawithDetails(ippant).TrialDetails(usetrial).TL_Catch==1
            %drawbox around catch times
            plot([catchstart catchstart], [0.5 1.5], ['g'], 'linewidth',linewidth);%box end
            plot([catchend catchend], [0.5 1.5], ['g'], 'linewidth',linewidth);
            plot([catchstart catchend], [0.5 0.5],  ['g'], 'linewidth',linewidth);
            l=plot([catchstart catchend], [1.5 1.5],  ['g'], 'linewidth',linewidth);
        end
        %
        if ppantTrialDatawithDetails(ippant).TrialDetails(usetrial).TR_Catch==1
            l=plot([catchstart catchstart], [1.5 2.5], ['g'], 'linewidth',linewidth);
            l=plot([catchend catchend], [1.5 2.5], ['g'], 'linewidth',linewidth);
            plot([catchstart catchend], [1.5 1.5],  ['g'], 'linewidth',linewidth);
            l=plot([catchstart catchend], [2.5 2.5],  ['g'], 'linewidth',linewidth);
        end
        %
        if ppantTrialDatawithDetails(ippant).TrialDetails(usetrial).BL_Catch==1
            plot([catchstart catchstart], [2.5 3.5], ['g'], 'linewidth',linewidth);
            plot([catchend catchend], [2.5 3.5], ['g'], 'linewidth',linewidth);
            plot([catchstart catchend], [2.5 2.5],  ['g'], 'linewidth',linewidth);
            l=plot([catchstart catchend], [3.5 3.5],  ['g'], 'linewidth',linewidth);
        end
        
        if ppantTrialDatawithDetails(ippant).TrialDetails(usetrial).BR_Catch==1
            plot([catchstart catchstart], [3.5 4.5], ['g'], 'linewidth',linewidth);
            plot([catchstart catchend], [3.5 3.5],  ['g'], 'linewidth',linewidth);
            %box end
            plot([catchend catchend], [3.5 4.5], ['g'], 'linewidth',linewidth);
            l=plot([catchstart catchend], [4.5 4.5],  ['g'], 'linewidth',linewidth);                        
        end
    
        
        set(gca, 'ytick', 1:4,'yticklabel', {'TL', 'TR', 'BL', 'BR'}) 
        xlabel('sec')
    end
    set(gcf, 'color', 'w');
    print('-dpng', ['All Trials Participant ' num2str(ippant) '-' num2str(ifig) ])
    end
%       
end
    
  %% collect number of trials rejected. make sure 'calcmissed catches x2' has been run in s1
  prej = [];
for ippant = 1:20
    %
    itrials = [2:49] + (ippant-1)*48;
    trialsare=[];
    for i=1:48
    trialsare =[trialsare, str2num(catchStruct{itrials(i),6})];
    end    
    %
    prej(ippant) = length(find(trialsare>0));
    
end
   
%%
% badppants = [8, 15,28, 4 , 5 6,10, 7,11]; %8, 15 ,28 based on catch. 7, based on trials, 
% badppants=[];
% prej = prej(goodppants);

% prej(badppants)=nan;
   %%
   clf
   hist(prej);
   xlabel('Trials marked for rejection')
ylabel('Subject count')
set(gca, 'fontsize', 25)
set(gcf, 'color', 'w')
% ylabel('# of subjects')
xlim([0 7])
shg

    %%
    percrej=nansum(prej)/( 24*26);
    %%
    iqr(prej)
    
    
    